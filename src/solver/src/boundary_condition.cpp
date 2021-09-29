#include "boundary_condition.h"
#include "node_container.h"
#include "element_container.h"
#include "linear_algebra.h"


//////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////
// DIRICHLET AND NEUMANN BOUNDARY CONDITION
void BoundaryCondition :: init(FunctionContainer *funcs) {
    blockedDoFNum = 0;
    loadedDoFNum = 0;

    dirichF.resize(dirichBC.size() );
    for ( unsigned i = 0; i < dirichBC.size(); i++ ) {
        if ( dirichBC [ i ] >= 0 ) {
            blockedDoFNum++;
            dirichF [ i ] = funcs->giveFunction(dirichBC [ i ]);
        } else {
            dirichF [ i ] = nullptr;
        }
    }

    neumannF.resize(neumannBC.size() );
    for ( unsigned i = 0; i < neumannBC.size(); i++ ) {
        if ( neumannBC [ i ] >= 0 ) {
            loadedDoFNum++;
            neumannF [ i ] = funcs->giveFunction(neumannBC [ i ]);
        } else {
            neumannF [ i ] = nullptr;
        }
    }
}

//////////////////////////////////////////////////////////
vector< unsigned >BoundaryCondition :: giveBlockedDoFs() const {
    vector< unsigned >blocked;
    blocked.resize(blockedDoFNum);
    unsigned k = 0;
    unsigned s = 0;
    for ( vector< int > :: const_iterator i = dirichBC.begin(); i != dirichBC.end(); ++i, k++ ) {
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
vector< double >BoundaryCondition :: giveBlockedDoFValues(double t) const {
    vector< double >blocked;
    blocked.resize(blockedDoFNum);
    unsigned s = 0;
    unsigned i = 0;
    for ( vector< Function * > :: const_iterator f = dirichF.begin(); f != dirichF.end(); ++f, i++ ) {
        if ( * f ) {
            blocked [ s ] = ( * f )->giveY(t) * multipliers [ i ];
            s++;
        }
    }
    return blocked;
}

//////////////////////////////////////////////////////////
vector< double >BoundaryCondition :: giveLoadedDoFValues(double t) const {
    vector< double >loaded;
    loaded.resize(loadedDoFNum);
    unsigned s = 0;
    unsigned i = 0;
    for ( vector< Function * > :: const_iterator f = neumannF.begin(); f != neumannF.end(); ++f, i++ ) {
        if ( * f ) {
            loaded [ s ] = ( * f )->giveY(t) * multipliers [ i ];
            s++;
        }
    }
    return loaded;
}

//////////////////////////////////////////////////////////
void BoundaryCondition :: readFromLine(istringstream &iss, NodeContainer *nodes) {
    unsigned intnum, nDoFs;
    iss >> intnum;
    node = nodes->giveNode(intnum);
    nDoFs = node->giveNumberOfDoFs();
    dirichBC.resize(nDoFs);
    neumannBC.resize(nDoFs);
    for ( unsigned i = 0; i < nDoFs; i++ ) {
        iss >> dirichBC [ i ];
    }
    for ( unsigned i = 0; i < nDoFs; i++ ) {
        iss >> neumannBC [ i ];
    }
    for ( unsigned i = 0; i < nDoFs; i++ ) {
        if ( neumannBC [ i ] >= 0 && dirichBC [ i ] >= 0 ) {
            cerr << "Error: Dirichlet and Neumann boundary conditions assigned simulatneuosly" << endl;
            exit(EXIT_FAILURE);
        }
    }
    multipliers.resize(nDoFs, 1.);
}

//////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////
// BODY LOAD
void BodyLoad :: readFromLine(istringstream &iss, ElementContainer *elems) {
    unsigned numelems, k;
    string code;

    while ( !iss.eof() ) {
        iss >> code;
        if ( code.compare("elems") == 0 ) {
            iss >> code;
            if ( code.compare("all") == 0 ) {
                numelems = elems->giveSize();
                els.resize(numelems);
                for ( unsigned i = 0; i < numelems; i++ ) {
                    els [ i ] = elems->giveElement(i);
                }
            } else {
                numelems = stoi(code);
                els.resize(numelems);
                for ( unsigned i = 0; i < numelems; i++ ) {
                    iss >> k;
                    els [ i ] = elems->giveElement(k);
                }
            }
        } else if ( code.compare("spatialFunction") == 0 ) {
            iss >> spatialFunctionNum;
        } else if ( code.compare("timeFunction") == 0 ) {
            iss >> timeFunctionNum;
        } else if ( code.compare("direction") == 0 ) {
            iss >> dir;
        }
    }
}

//////////////////////////////////////////////////////////
void BodyLoad :: init(FunctionContainer *funcs) {
    timeFunction = funcs->giveFunction(timeFunctionNum);
    spatialFunction = funcs->giveFunction(spatialFunctionNum);
}

//////////////////////////////////////////////////////////
double BodyLoad :: giveValue(const Point *xyz, double t) {
    return timeFunction->giveY(t) * spatialFunction->giveY(xyz);
}

//////////////////////////////////////////////////////////
vector< double >BodyLoad :: giveBodyForceDoFValues(double t) {
    vector< double >load, elemLoad;
    for ( auto &e: els ) {
        elemLoad = e->integrateLoad(this, t);
        load.insert(load.end(), elemLoad.begin(), elemLoad.end() );
    }
    return load;
}

//////////////////////////////////////////////////////////
vector< unsigned >BodyLoad :: giveArrayOfBodyForceDoFs() const {
    vector< unsigned >DoFs, elemDoFs;
    for ( auto &e: els ) {
        elemDoFs = e->giveDoFsInDirection(dir);
        DoFs.insert(DoFs.end(), elemDoFs.begin(), elemDoFs.end() );
    }
    return DoFs;
}

//////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////
// CONTAINER FOR BOUNDARY CONDITIONS
BCContainer :: ~BCContainer() {
    for ( auto &bc: BC ) {
        if(bc!=nullptr) delete bc;
    }
    for ( auto &vl: loads ) {
        if(vl!=nullptr) delete vl;
    }
}

void BCContainer :: clear() {
    for ( auto &bc: BC ) {
        if(bc!=nullptr) delete bc;
    }
    for ( auto &vl: loads ) {
        if(vl!=nullptr) delete vl;
    }
}

//////////////////////////////////////////////////////////
void BCContainer :: readFromFile(const string filename, NodeContainer *nodes, ElementContainer *elems) {
    size_t origBCsize = BC.size();
    size_t origLoadsize = loads.size();
    string line, aux;
    ifstream inputfile(filename.c_str() );
    if ( inputfile.is_open() ) {
        while ( getline(inputfile >> std :: ws, line) ) {
            if ( line.empty() ) {
                continue;
            }
            if ( line.at(0) == '#' ) {
                continue;
            }
            istringstream iss(line);
            iss >> aux;
            if ( aux.compare("NodalBC") == 0 ) {
                BoundaryCondition *newbc = new BoundaryCondition();
                newbc->readFromLine(iss, nodes);
                BC.push_back(newbc);
            } else if ( aux.compare("BodyLoad") == 0 ) {
                BodyLoad *newBodyLoad = new BodyLoad();
                newBodyLoad->readFromLine(iss, elems);
                loads.push_back(newBodyLoad);
            } else {
                cerr << "Error: boundary condition '" <<  aux <<  "' in not implemented" << endl;
                cerr << "Did you forget keyword 'NodalBC'?" << endl;
                exit(EXIT_FAILURE);
            }
        }
        inputfile.close();
        cout << "Input file '" <<  filename << "' succesfully loaded; " << BC.size() - origBCsize << " boundary conditions and " << loads.size() - origLoadsize << " volume loads found" << endl;
    } else {
        cerr << "Error: unable to open input file '" <<  filename <<  "'" << endl;
        exit(EXIT_FAILURE);
    }
}

//////////////////////////////////////////////////////////
void BCContainer :: init() {
    dirichDoFs.resize(0);
    neumannDoFs.resize(0);

    for ( auto &bc: BC ) {
        bc->init(functions);
    }

    for ( auto &l: loads ) {
        l->init(functions);
    }
}

//////////////////////////////////////////////////////////
void BCContainer :: calculateDoFfields() {
    vector< unsigned >help;
    vector< double >help2;
    for ( vector< BoundaryCondition * > :: iterator bc = BC.begin(); bc != BC.end(); ++bc ) {
        help = ( * bc )->giveBlockedDoFs();
        dirichDoFs.insert(dirichDoFs.end(), help.begin(), help.end() );
        help = ( * bc )->giveLoadedDoFs();
        neumannDoFs.insert(neumannDoFs.end(), help.begin(), help.end() );
    }

    // NOTE know which fns are actually used //JE WHY? What is this information for? Nobody cares. // JK to prevent restricting time step in extreme points of unused fns, it is probably not necessary and if anyone does not comment (or remove fn from fn file) fn that is not used it is his problem, can be removed then
    // for (auto const &f_id : dirichF ){
    //   if ( !functions->isActive(f_id) ){
    //     functions->setActive(f_id);
    //   }
    // }
    // for (auto const &f_id : neumannF ){
    //   if ( !functions->isActive(f_id) ){
    //     functions->setActive(f_id);
    //   }
    // }
}

//////////////////////////////////////////////////////////
vector< double >BCContainer :: giveBlockedDoFValues(double t) const {
    vector< double >blocked(dirichDoFs.size() );
    unsigned i, s = 0;
    vector< double >b;
    for ( auto &bc: BC ) {
        b = bc->giveBlockedDoFValues(t);
        for ( i = 0; i < b.size(); i++, s++ ) {
            blocked [ s ] = b [ i ];
        }
    }
    return blocked;
}

//////////////////////////////////////////////////////////
vector< double >BCContainer :: giveLoadedDoFValues(double t) const {
    vector< double >loaded(neumannDoFs.size() );
    unsigned i, s = 0;
    vector< double >b;
    for ( auto &bc: BC ) {
        b = bc->giveLoadedDoFValues(t);
        for ( i = 0; i < b.size(); i++, s++ ) {
            loaded [ s ] = b [ i ];
        }
    }
    return loaded;
}

//////////////////////////////////////////////////////////
vector< unsigned >BCContainer :: giveArrayOfBodyForceDoFs() const {
    vector< unsigned >DoFs, elemDoFs;
    for ( auto &l: loads ) {
        elemDoFs = l->giveArrayOfBodyForceDoFs();
        DoFs.insert(DoFs.end(), elemDoFs.begin(), elemDoFs.end() );
    }
    return DoFs;
}

//////////////////////////////////////////////////////////
vector< double >BCContainer :: giveBodyForceDoFValues(double t) {
    vector< double >structLoads, elemLoads;
    for ( auto &l: loads ) {
        elemLoads = l->giveBodyForceDoFValues(t);
        structLoads.insert(structLoads.end(), elemLoads.begin(), elemLoads.end() );
    }
    return structLoads;
}


//////////////////////////////////////////////////////////
void BCContainer :: removeBoundaryCondition(unsigned i) {
    if ( i > BC.size() - 1 ) {
        cerr << "BCContainer Error: requester BC number " << i << " out of " << BC.size() << endl;
        exit(1);
    }
    delete BC [ i ];
    BC [ i ] = nullptr;
    BC.erase(BC.begin() + i);
}
