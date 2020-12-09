#include "indirect_displ_control.h"


//////////////////////////////////////////////////////////
IndirectDC :: IndirectDC() {
    name = "indirect displacement controller";
    funcnum = -1;
    func = nullptr;
    nummaxunit = 0;
};

//////////////////////////////////////////////////////////
void IndirectDC :: readFromStream(unsigned num, ifstream &inputfile) {
    nummaxunit++;

    c_nodes.resize(nummaxunit);
    c_dirs.resize(nummaxunit);
    c_weights.resize(nummaxunit);
    c_DoFs.resize(nummaxunit);
    xcoords.resize(nummaxunit);
    ycoords.resize(nummaxunit);
    zcoords.resize(nummaxunit);
    coords_active.resize(nummaxunit);
    coords_active [ nummaxunit - 1 ] = false;
    nodes_active.resize(nummaxunit);
    nodes_active [ nummaxunit - 1 ] = false;

    c_nodes [ nummaxunit - 1 ].resize(num, 0);
    c_dirs [ nummaxunit - 1 ].resize(num, 0);
    c_weights [ nummaxunit - 1 ].resize(num, 0);
    xcoords [ nummaxunit - 1 ].resize(num, 0);
    ycoords [ nummaxunit - 1 ].resize(num, 0);
    zcoords [ nummaxunit - 1 ].resize(num, 0);

    string param, line;

    streampos oldpos = inputfile.tellg();  // stores the position
    while ( getline(inputfile >> std :: ws, line) ) {
        if ( line.empty() ) {
            continue;
        }
        if ( line.at(0) == '#' ) {
            continue;
        }
        istringstream iss(line);
        iss >> param;
        if ( param.compare("idc_nodes") == 0 ) {
            coords_active [ nummaxunit - 1 ] = true;
            for ( unsigned j = 0; j < num; j++ ) {
                iss >> c_nodes [ nummaxunit - 1 ] [ j ];
            }
        } else if ( param.compare("idc_xcoords") == 0 ) {
            coords_active [ nummaxunit - 1 ] = true;
            for ( unsigned j = 0; j < num; j++ ) {
                iss >> xcoords [ nummaxunit - 1 ] [ j ];
            }
        } else if ( param.compare("idc_ycoords") == 0 ) {
            for ( unsigned j = 0; j < num; j++ ) {
                iss >> ycoords [ nummaxunit - 1 ] [ j ];
            }
        } else if ( param.compare("idc_zcoords") == 0 ) {
            for ( unsigned j = 0; j < num; j++ ) {
                iss >> zcoords [ nummaxunit - 1 ] [ j ];
            }
        } else if ( param.compare("idc_directions") == 0 ) {
            for ( unsigned j = 0; j < num; j++ ) {
                iss >> c_dirs [ nummaxunit - 1 ] [ j ];
            }
        } else if ( param.compare("idc_weights") == 0 ) {
            for ( unsigned j = 0; j < num; j++ ) {
                iss >> c_weights [ nummaxunit - 1 ] [ j ];
            }
        } else if ( param.compare("idc_function") == 0 ) {
            iss >> funcnum;
        } else {
            inputfile.seekg(oldpos);    // get back to the position
            return;
        }
        oldpos = inputfile.tellg();  // stores the position
    }
}



//////////////////////////////////////////////////////////
void IndirectDC :: init(NodeContainer *nodes, FunctionContainer *funcs) {
    unsigned clength;
    Node *n;
    if ( funcnum >= 0 ) {
        func = funcs->giveFunction(funcnum);
    }
    for ( unsigned c = 0; c < nummaxunit; c++ ) {
        clength = c_weights [ c ].size(); //todo: warning C4267: '=': conversion from 'size_t' to 'unsigned int', possible loss of data
        if ( clength == 0 ) {
            cerr << "Error: Indirect displacement controll weights were not set" << endl;
            exit(1);
        }
        if ( c_DoFs [ c ].size() < clength ) {
            c_DoFs [ c ].resize(clength);
            if ( nodes_active [ c ] ) {
                for ( unsigned i = 0; i < clength; i++ ) {
                    c_DoFs [ c ] [ i ] = nodes->giveNode(c_nodes [ c ] [ i ])->giveStartingDoF() + c_dirs [ c ] [ i ];
                }
            } else if ( coords_active [ c ] ) {
                for ( unsigned i = 0; i < clength; i++ ) {
                    n = nodes->findClosestMechanicalNode(Point(xcoords [ c ] [ i ], ycoords [ c ] [ i ], zcoords [ c ] [ i ]) );
                    c_DoFs [ c ] [ i ] = n->giveStartingDoF() + c_dirs [ c ] [ i ];
                }
            } else {
                cerr << "Error: Indirect displacement controll was not correctly set" << endl;
                exit(1);
            }
        }
    }
}

//////////////////////////////////////////////////////////
double IndirectDC :: giveMultiplierCorrection(Vector &prev_displ, Vector &displ_d, Vector &displ_f, double time) {
    double dd = -INFINITY;
    double df = 0;
    double m;
    for ( unsigned c = 0; c < nummaxunit; c++ ) {
        m = 0;
        for ( unsigned i = 0; i < c_weights [ c ].size(); i++ ) {
            m += ( prev_displ [ c_DoFs [ c ] [ i ] ] ) * c_weights [ c ] [ i ];
        }
        if ( m > dd ) {
            dd = m;
            df = 0;
            for ( unsigned i = 0; i < c_weights [ c ].size(); i++ ) {
                df += displ_f [ c_DoFs [ c ] [ i ] ] * c_weights [ c ] [ i ];
            }
        }
    }
    return ( givePrescribedDisplacement(time) - dd ) / df;
}

//////////////////////////////////////////////////////////
double IndirectDC :: giveControlValue(Vector &displ) {
    double dd = -INFINITY;
    double m;
    for ( unsigned c = 0; c < nummaxunit; c++ ) {
        m = 0;
        for ( unsigned i = 0; i < c_weights [ c ].size(); i++ ) {
            m += displ [ c_DoFs [ c ] [ i ] ] * c_weights [ c ] [ i ];
        }
        if ( m > dd ) {
            dd = m;
        }
    }
    return dd;
}

//////////////////////////////////////////////////////////
double IndirectDC :: givePrescribedDisplacement(double time) {
    if ( func ) {
        return func->giveY(time);
    }
    return time;
}
