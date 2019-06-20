#include "boundary_condition.h"
#include "node_container.h"

//////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////
// PIECE-WISE LINEAR FUNCTION
void PieceWiseLinearFunction :: readFromLine(istringstream &iss) {
    unsigned num;
    iss >> num;
    x.resize(num);
    y.resize(num);
    for ( unsigned i = 0; i < num; i++ ) {
        iss >> x [ i ];
    }
    for ( unsigned i = 0; i < num; i++ ) {
        iss >> y [ i ];
    }
}

//////////////////////////////////////////////////////////
double PieceWiseLinearFunction :: giveY(double t) const {
    unsigned i = 0;
    while ( x [ i ] < t && i < x.size() ) {
        i++;
    }
    if ( i == 0 ) {
        return 0.;
    } else if ( i == x.size() )   {
        return y [ x.size() - 1 ];
    } else                                                       {
        return y [ i - 1 ] + ( y [ i ] - y [ i - 1 ] ) / ( x [ i ] - x [ i - 1 ] ) * ( t - x [ i - 1 ] );
    }
}

//////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////
// CONTAINER FOR FUNCTIONS
FunctionContainer :: ~FunctionContainer() {
    for ( vector< Function * > :: iterator f = functions.begin(); f != functions.end(); ++f ) {
        delete * f;
    }
}

//////////////////////////////////////////////////////////
void FunctionContainer :: readFromFile(const string filename) {
    unsigned origsize = functions.size();
    string line, ftype;
    ifstream inputfile(filename.c_str() );
    if ( inputfile.is_open() ) {
        while ( getline(inputfile, line) ) {
            if ( line.at(0) == '#' ) {
                continue;
            }
            istringstream iss(line);
            iss >> ftype;
            if ( !ftype.rfind("#", 0) == 0 ) {
                if ( ftype.compare("PWLFunction") == 0 ) {
                    PieceWiseLinearFunction *newf = new PieceWiseLinearFunction();
                    newf->readFromLine(iss);
                    functions.push_back(newf);
                } else  {
                    cerr << "Error: function '" <<  ftype <<  "' is not implemented yet." << endl;
                    exit(0);
                }
            }
        }
        inputfile.close();
        cout << "Input file '" <<  filename << "' succesfully loaded; " << functions.size() - origsize << " functions found" << endl;
    } else {
        cerr << "Error: unable to open input file '" <<  filename <<  "'" << endl;
        exit(0);
    }
}

//////////////////////////////////////////////////////////
double FunctionContainer :: giveY(unsigned f, double t) const {
    return functions [ f ]->giveY(t);
}

//////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////
// DIRRICHLET AND NEUMANN BOUNDARY CONDITION
void BoundaryCondition :: init() {
    blockedDoFNum = 0;
    loadedDoFNum = 0;
    for ( vector< int > :: const_iterator i = dirrichBC.begin(); i != dirrichBC.end(); ++i ) {
        if ( * i >= 0 ) {
            blockedDoFNum++;
        }
    }
    for ( vector< int > :: const_iterator i = neumannBC.begin(); i != neumannBC.end(); ++i ) {
        if ( * i >= 0 ) {
            loadedDoFNum++;
        }
    }
}

//////////////////////////////////////////////////////////
vector< unsigned >BoundaryCondition :: giveBlockedDoFs() const {
    vector< unsigned >blocked;
    blocked.resize(blockedDoFNum);
    unsigned k = 0;
    unsigned s = 0;
    for ( vector< int > :: const_iterator i = dirrichBC.begin(); i != dirrichBC.end(); ++i, k++ ) {
        if ( * i >= 0 ) {
            blocked [ s ] = node->giveStartingDoF() + k;
            s++;
        }
    }
    ;
    return blocked;
}

//////////////////////////////////////////////////////////
vector< unsigned >BoundaryCondition :: giveLoadedDoFs() const {
    vector< unsigned >loaded;
    loaded.resize(loadedDoFNum);
    unsigned k = 0;
    unsigned s = 0;
    for ( vector< int > :: const_iterator i = neumannBC.begin(); i != neumannBC.end(); ++i, k++ ) {
        if ( * i >= 0 ) {
            loaded [ s ] = node->giveStartingDoF() + k;
            s++;
        }
    }
    ;
    return loaded;
}

//////////////////////////////////////////////////////////
vector< unsigned >BoundaryCondition :: giveBlockedFunctions() const {
    vector< unsigned >blocked;
    blocked.resize(blockedDoFNum);
    unsigned k = 0;
    unsigned s = 0;
    for ( vector< int > :: const_iterator i = dirrichBC.begin(); i != dirrichBC.end(); ++i, k++ ) {
        if ( * i >= 0 ) {
            blocked [ s ] = dirrichBC [ k ];
            s++;
        }
    }
    ;
    return blocked;
}

//////////////////////////////////////////////////////////
vector< unsigned >BoundaryCondition :: giveLoadedFunctions() const {
    vector< unsigned >loaded;
    loaded.resize(loadedDoFNum);
    unsigned k = 0;
    unsigned s = 0;
    for ( vector< int > :: const_iterator i = neumannBC.begin(); i != neumannBC.end(); ++i, k++ ) {
        if ( * i >= 0 ) {
            loaded [ s ] = neumannBC [ k ];
            s++;
        }
    }
    ;
    return loaded;
}

//////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////
// CONTAINER FOR BOUNDARY CONDITIONS
BCContainer :: ~BCContainer() {
    for ( vector< BoundaryCondition * > :: iterator bc = BC.begin(); bc != BC.end(); ++bc ) {
        delete * bc;
    }
}

//////////////////////////////////////////////////////////
void BCContainer :: readFromFile(const string filename, NodeContainer *nodes) {
    unsigned origsize = BC.size();
    string line, aux;
    unsigned intnum, nDoFs;
    vector< int >dirrichBC, neumannBC;
    Node *node;
    ifstream inputfile(filename.c_str() );

    if ( inputfile.is_open() ) {
        while ( getline(inputfile, line) ) {
            if ( line.at(0) == '#' ) {
                continue;
            }
            istringstream iss(line);
            iss >> intnum;
            node = nodes->giveNode(intnum);
            nDoFs = node->giveNumberOfDoFs();
            dirrichBC.resize(nDoFs);
            neumannBC.resize(nDoFs);
            for ( unsigned i = 0; i < nDoFs; i++ ) {
                iss >> dirrichBC [ i ];
            }
            for ( unsigned i = 0; i < nDoFs; i++ ) {
                iss >> neumannBC [ i ];
            }
            for ( unsigned i = 0; i < nDoFs; i++ ) {
                if ( neumannBC [ i ] >= 0 && dirrichBC [ i ] >= 0 ) {
                    cerr << "Error: Dirrichlet and Neumann boundary conditions assigned simulatneuosly" << endl;
                    cerr << line << endl;
                    exit(0);
                }
            }
            ;
            BoundaryCondition *newBC = new BoundaryCondition(node, dirrichBC, neumannBC);
            BC.push_back(newBC);
        }
        inputfile.close();
        for ( unsigned i = 0; i < BC.size() - origsize; i++ ) {
            BC [ i + origsize ]->giveNode()->setBC(BC [ i + origsize ]);
        }
        cout << "Input file '" <<  filename << "' succesfully loaded; " << BC.size() - origsize << " boundary conditions found" << endl;
    } else {
        cerr << "Error: unable to open input file '" <<  filename <<  "'" << endl;
        exit(0);
    }
}

//////////////////////////////////////////////////////////
void BCContainer :: init() {
    dirrichDoFs.resize(0);
    neumannDoFs.resize(0);
    for ( vector< BoundaryCondition * > :: iterator bc = BC.begin(); bc != BC.end(); ++bc ) {
        ( * bc )->init();
    }
}

//////////////////////////////////////////////////////////
void BCContainer :: calculateDoFfields() {
    vector< unsigned >help;
    for ( vector< BoundaryCondition * > :: iterator bc = BC.begin(); bc != BC.end(); ++bc ) {
        help = ( * bc )->giveBlockedDoFs();
        dirrichDoFs.insert(dirrichDoFs.end(), help.begin(), help.end() );
        help = ( * bc )->giveLoadedDoFs();
        neumannDoFs.insert(neumannDoFs.end(), help.begin(), help.end() );
        help = ( * bc )->giveBlockedFunctions();
        dirrichF.insert(dirrichF.end(), help.begin(), help.end() );
        help = ( * bc )->giveLoadedFunctions();
        neumannF.insert(neumannF.end(), help.begin(), help.end() );
    }
}

//////////////////////////////////////////////////////////
vector< double >BCContainer :: giveBlockedDoFValues(double t) const {
    vector< double >blocked(dirrichDoFs.size() );
    for ( unsigned h = 0; h < dirrichDoFs.size(); h++ ) {
        blocked [ h ] = functions->giveY(dirrichF [ h ], t);
    }
    return blocked;
}

//////////////////////////////////////////////////////////
vector< double >BCContainer :: giveLoadedDoFValues(double t) const {
    vector< double >loaded(neumannDoFs.size() );
    for ( unsigned h = 0; h < neumannDoFs.size(); h++ ) {
        loaded [ h ] = functions->giveY(neumannF [ h ], t);
    }
    return loaded;
}
