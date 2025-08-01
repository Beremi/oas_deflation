#include "indirect_control.h"
#include "model.h"

using namespace std;

//////////////////////////////////////////////////////////
IndirectControl :: IndirectControl() {
    name = "indirect Controler";
    funcnum = -1;
    func = nullptr;
    nummaxunit = 0;
    requiref = false;
};

//////////////////////////////////////////////////////////
void IndirectControl :: readFromStream(unsigned num, ifstream &inputfile) {
    nummaxunit++;

    nodes.resize(nummaxunit);
    dirs.resize(nummaxunit);
    r_weights.resize(nummaxunit);
    f_weights.resize(nummaxunit);
    DoFs.resize(nummaxunit);
    xcoords.resize(nummaxunit);
    ycoords.resize(nummaxunit);
    zcoords.resize(nummaxunit);
    coords_active.resize(nummaxunit);
    coords_active [ nummaxunit - 1 ] = false;
    nodes_active.resize(nummaxunit);
    nodes_active [ nummaxunit - 1 ] = false;

    nodes [ nummaxunit - 1 ].resize(num, 0);
    dirs [ nummaxunit - 1 ].resize(num, 0);
    r_weights [ nummaxunit - 1 ].resize(num, 0);
    f_weights [ nummaxunit - 1 ].resize(num, 0);
    xcoords [ nummaxunit - 1 ].resize(num, 0);
    ycoords [ nummaxunit - 1 ].resize(num, 0);
    zcoords [ nummaxunit - 1 ].resize(num, 0);

    string param, line;

    streampos oldpos = inputfile.tellg();  // stores the position
    while ( getline(inputfile >> std :: ws, line) ) {
        if ( line.empty() || ( line.at(0) == '#' ) ) {
            continue;
        }
        istringstream iss(line);
        iss >> param;
        if ( param.compare("ic_nodes") == 0 || param.compare("idc_nodes") == 0 ) {
            nodes_active [ nummaxunit - 1 ] = true;
            for ( unsigned j = 0; j < num; j++ ) {
                iss >> nodes [ nummaxunit - 1 ] [ j ];
                //cout << "IDC node " << j << c_nodes [ nummaxunit - 1 ] [ j ] << endl;
            }
        } else if ( param.compare("ic_xcoords") == 0 || param.compare("idc_xcoords") == 0 ) {
            coords_active [ nummaxunit - 1 ] = true;
            //cout << "IDC xcoords: ";
            for ( unsigned j = 0; j < num; j++ ) {
                iss >> xcoords [ nummaxunit - 1 ] [ j ];
                //cout << "IDC xcoords " << xcoords [ nummaxunit - 1 ] [ j ] << endl;
            }
            //cout << endl;
        } else if ( param.compare("ic_ycoords") == 0 || param.compare("idc_ycoords") == 0 ) {
            //cout << "IDC ycoords: ";
            for ( unsigned j = 0; j < num; j++ ) {
                iss >> ycoords [ nummaxunit - 1 ] [ j ];
                //cout << "IDC ycoords " << xcoords [ nummaxunit - 1 ] [ j ] << endl;
            }
            //cout << endl;
        } else if ( param.compare("ic_zcoords") == 0 || param.compare("idc_zcoords") == 0 ) {
            for ( unsigned j = 0; j < num; j++ ) {
                iss >> zcoords [ nummaxunit - 1 ] [ j ];
            }
        } else if ( param.compare("ic_directions") == 0 || param.compare("idc_directions") == 0) {
            for ( unsigned j = 0; j < num; j++ ) {
                iss >> dirs [ nummaxunit - 1 ] [ j ];
            }
        } else if ( param.compare("ic_displ_weights") == 0 || param.compare("idc_weights") == 0 ) {
            for ( unsigned j = 0; j < num; j++ ) {
                iss >> r_weights [ nummaxunit - 1 ] [ j ];
            }
        } else if ( param.compare("ic_force_weights") == 0 ) {
            for ( unsigned j = 0; j < num; j++ ) {
                iss >> f_weights [ nummaxunit - 1 ] [ j ];
            }
        } else if ( param.compare("ic_function") == 0 || param.compare("idc_function") == 0 ) {
            iss >> funcnum;
        } else {
            inputfile.seekg(oldpos);    // get back to the position
            return;
        }
        oldpos = inputfile.tellg();  // stores the position
    }
}



//////////////////////////////////////////////////////////
void IndirectControl :: init(NodeContainer *mnodes, FunctionContainer *funcs, bool initial) {
    unsigned clength;
    Node *n;
    double dist;
    if ( funcnum >= 0 ) {
        func = funcs->giveFunction(funcnum);
    }
    for ( unsigned c = 0; c < nummaxunit; c++ ) {
        clength = r_weights [ c ].size(); //todo: warning C4267: '=': conversion from 'size_t' to 'unsigned int', possible loss of data
        if ( clength == 0 ) {
            cerr << "Error: Indirect displacement Control weights were not set" << endl;
            exit(1);
        }
        clength = f_weights [ c ].size(); //todo: warning C4267: '=': conversion from 'size_t' to 'unsigned int', possible loss of data
        if ( clength == 0 ) {
            cerr << "Error: Indirect force Control weights were not set" << endl;
            exit(1);
        }
        for ( unsigned i = 0; i < clength; i++ ) {
            if (f_weights [ c ] [ i ] != 0) requiref = true;
        }
        if ( DoFs [ c ].size() < clength || !initial ) {
            // JK: in adaptivity, number of control DoFs remain, but DoFs from updated geometry are used
            if ( initial ) {
                DoFs [ c ].resize(clength);
            }
            if ( nodes_active [ c ] ) {
                for ( unsigned i = 0; i < clength; i++ ) {
                    DoFs [ c ] [ i ] = mnodes->giveNode(nodes [ c ] [ i ])->giveStartingDoF() + dirs [ c ] [ i ];
                }
            } else if ( coords_active [ c ] ) {
                for ( unsigned i = 0; i < clength; i++ ) {
                    n = mnodes->findClosestMechanicalNode(Point(xcoords [ c ] [ i ], ycoords [ c ] [ i ], zcoords [ c ] [ i ]), & dist);
                    DoFs [ c ] [ i ] = n->giveStartingDoF() + dirs [ c ] [ i ];
                }
            } else {
                cerr << "Error: Indirect displacement Control was not correctly set" << endl;
                exit(1);
            }
        }
    }
}

//////////////////////////////////////////////////////////
double IndirectControl :: giveMultiplierCorrection(Vector &prev_displ, Vector &prev_force, Vector &diff_displ, Vector &diff_force, double time) {
    double df, dd;
    double pdispl = givePrescribedValue(time);
    double lambda = INFINITY;
    double lambda_temp;
    for ( unsigned c = 0; c < nummaxunit; c++ ) {
        dd = 0;
        df = 0;
        for ( unsigned i = 0; i < r_weights [ c ].size(); i++ ) {
            dd += prev_displ [ DoFs [ c ] [ i ] ] * r_weights [ c ] [ i ] + prev_force [ DoFs [ c ] [ i ] ] * f_weights [ c ] [ i ];
            df += diff_displ [ DoFs [ c ] [ i ] ] * r_weights [ c ] [ i ] + diff_force [ DoFs [ c ] [ i ] ] * f_weights [ c ] [ i ];
        }
        lambda_temp = ( pdispl - dd ) / df;
        if ( lambda_temp < lambda ) {
            lambda = lambda_temp;
        }
    }
    return lambda;
}

//////////////////////////////////////////////////////////
double IndirectControl :: givePrescribedValue(double time) {
    if ( func ) {
        return func->giveY(time);
    }
    return time;
}

//////////////////////////////////////////////////////////
bool IndirectControl :: requireForces() const {
    return requiref;
}

//////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////
IndirectControlSumOfSquares :: IndirectControlSumOfSquares() {
    name = "indirect Controler using sum of squares";
};

//////////////////////////////////////////////////////////
double IndirectControlSumOfSquares :: giveMultiplierCorrection(Vector &prev_displ, Vector &prev_force, Vector &diff_displ, Vector &diff_force, double time) {
    double df, dd;
    double pdispl = givePrescribedValue(time);
    double lambda = INFINITY;
    double lambda_temp;
    for ( unsigned c = 0; c < nummaxunit; c++ ) {
        dd = 0;
        df = 0;
        for ( unsigned i = 0; i < r_weights [ c ].size(); i++ ) {
            dd += pow(prev_displ [ DoFs [ c ] [ i ] ], 2) * r_weights [ c ] [ i ] + pow(prev_force [ DoFs [ c ] [ i ] ], 2) * f_weights [ c ] [ i ];
            df += pow(diff_displ [ DoFs [ c ] [ i ] ], 2) * r_weights [ c ] [ i ] + pow(diff_force [ DoFs [ c ] [ i ] ], 2) * f_weights [ c ] [ i ];
        }        
        lambda_temp = abs( ( pdispl - sqrt(dd) ) / sqrt(df) );
        if (pdispl > sqrt(dd)) lambda_temp *= -1.;
        if ( lambda_temp < lambda ) {
            lambda = lambda_temp;
        }
    }
    return lambda;
}

