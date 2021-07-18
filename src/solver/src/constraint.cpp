#include "constraint.h"
#include "node_container.h"
#include "element_container.h"
#include "element_discrete.h"
#include "boundary_condition.h"
#include <fstream>

//////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////
// Joint Degree of Freedom
JointDoF :: JointDoF(Node *s, const unsigned &dir, const std :: vector< Node * > &m, const std :: vector< unsigned > &dirs, const std :: vector< double > &mult, const std :: vector< Function * > &fns, const std :: vector< double > &time_mult) {
    slaveNode = s;
    direction = dir;
    masters = m;
    directions = dirs;
    multipliers = mult;
    if ( fns.empty() ) {
        this->time_fns = std :: vector< Function * >(multipliers.size(), nullptr);
    } else {
        this->time_fns = fns;
        this->additional_term = time_mult;
    }
}

//////////////////////////////////////////////////////////
void JointDoF :: readFromLine(istringstream &iss, NodeContainer *nodes) {
    // # type      node    direction   numMasters master1 multiplier1 master2...
    unsigned intn, intmas;
    double mult;
    iss >> intn;
    slaveNode = nodes->giveNode(intn);
    iss >> intn;
    direction = intn;
    iss >> intn;
    for ( unsigned i = 0; i < intn; i++ ) {
        iss >> intmas >> mult;
        masters.push_back( nodes->giveNode(intmas) );
        directions.push_back(0);
        multipliers.push_back(mult);
    }
}

//////////////////////////////////////////////////////////
unsigned JointDoF :: giveSlaveDoF() const {
    return slaveNode->giveStartingDoF() + direction;
}

//////////////////////////////////////////////////////////
unsigned JointDoF :: giveMasterDoF(unsigned k) const {
    return masters [ k ]->giveStartingDoF() + directions [ k ];
}
//////////////////////////////////////////////////////////
void JointDoF :: print() {
    std :: cout << "slave DoF = " << giveSlaveDoF() << '\n';
    for ( unsigned i = 0; i < masters.size(); i++ ) {
        std :: cout << "master DoF[ " << i << " ] = " << masters [ i ]->giveStartingDoF() + directions [ i ] << " with multiplier = " << multipliers [ i ];
        if ( this->time_fns [ i ] != nullptr ) {
            std :: cout << " and time_dep_mult = " << additional_term [ i ];
        }
        std :: cout << '\n';
    }
}

//////////////////////////////////////////////////////////
void JointDoF :: init() {
    //consolidate (do not allow repetition of masters)
    for ( unsigned k = masters.size() - 1; k > 0; k-- ) {
        for ( unsigned l = 0; l < k; l++ ) {
            if ( masters [ l ] == masters [ k ] && directions [ l ] == directions [ k ] ) {
                multipliers [ l ] += multipliers [ k ];
                masters.erase(masters.begin() + k);
                directions.erase(directions.begin() + k);
                multipliers.erase(multipliers.begin() + k);
                break;
            }
        }
    }
}

//////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////
// Volumetric Average
VolumetricAverage :: VolumetricAverage(vector< Node * > &n, std :: vector< unsigned > &d, Node *mn, unsigned md, ElementContainer *ec, ConstraintContainer *cc) {
    nodes = n;
    dirs = d;
    masternode = mn;
    masterdir = md;
    elems = ec;
    constraints = cc;


    //collect all slaves from other joint DoFs and also all involved DOFs
    vector< Node * >excludedNodes;
    vector< unsigned >excudedDirs;
    vector< Node * >jdm;
    vector< unsigned >jdd;
    JointDoF *jd;
    VolumetricAverage *va;
    for ( unsigned j = 0; j < constraints->giveSize(); j++ ) {
        jd = constraints->giveConstraint(j);
        va = dynamic_cast< VolumetricAverage * >( jd );
        if ( !va ) {
            excludedNodes.push_back( jd->giveSlaveNode() );
            jdm = jd->giveMasterNodes();
            excludedNodes.insert( excludedNodes.end(), jdm.begin(), jdm.end() );
            excudedDirs.push_back( jd->giveSlaveDir() );
            jdd = jd->giveMasterDirs();
            excudedDirs.insert( excudedDirs.end(), jdd.begin(), jdd.end() );
        }
    }


    // find DoF to become a slave (cannot be slave or master of any other constraint)
    unsigned k;
    bool found;
    for ( k = 0; k < nodes.size(); k++ ) {
        found = false;
        for ( unsigned l = 0; l < excludedNodes.size(); l++ ) {
            if ( excludedNodes [ l ] == nodes [ k ] && excudedDirs [ l ] == dirs [ k ] ) {
                found = true;
                break;
            }
        }
        if ( !found ) {
            break;
        }
    }
    if ( found ) { //Todo potentially uninitialized variable 'found' used!
        cerr << "volumetric average didn't find any suitable slave DoF" << endl;
        exit(1);
    }
    slaveNode = nodes [ k ];
    direction = dirs [ k ];
}

//////////////////////////////////////////////////////////
void VolumetricAverage :: init() {
    //collect all slaves from other joint DoFs
    vector< unsigned >otherslaves;
    otherslaves.resize( constraints->giveSize() );
    JointDoF *jd;
    for ( unsigned j = 0; j < constraints->giveSize(); j++ ) {
        jd = constraints->giveConstraint(j);
        otherslaves [ j ] = jd->giveSlaveDoF();
    }

    //find slave position in array
    unsigned slaveid = find(nodes.begin(), nodes.end(), slaveNode) - nodes.begin();
    if ( slaveid == nodes.size() ) {
        cerr << "Volumetric average Error: slave node not included in node list" << endl;
        exit(1);
    }

    // calculate volumes associated with nodes
    unsigned r;
    unsigned s = nodes.size();
    vector< double >m;
    m.resize( nodes.size() );
    Transp1D *et;
    RigidBodyContact *em;
    for ( unsigned e = 0; e < elems->giveSize(); e++ ) {
        et = dynamic_cast< Transp1D * >( elems->giveElement(e) );
        em = dynamic_cast< RigidBodyContact * >( elems->giveElement(e) );
        if ( et ) {
            for ( unsigned p = 0; p < 2; p++ ) {
                r = ( std :: find( nodes.begin(), nodes.end(), et->giveNode(p) ) - nodes.begin() );
                if ( r < s ) {
                    m [ r ] += et->giveVolume(p);
                }
            }
        } else if ( em ) {
            for ( unsigned p = 0; p < 2; p++ ) {
                r = ( std :: find( nodes.begin(), nodes.end(), em->giveNode(p) ) - nodes.begin() );
                if ( r < s ) {
                    m [ r ] += em->giveVolume(p);
                }
            }
        }
    }

    //rearange everything and update weights
    double factor = m [ slaveid ];
    cout << slaveid << " " << factor << m.size() << endl;
    double fullVolume = 0;
    for ( auto &ss: m ) {
        fullVolume += ss;
        ss /= -factor;
    }
    cout << "volumetric Average: check of volume " << fullVolume << endl;
    m [ slaveid ] = fullVolume / factor;
    nodes [ slaveid ] = masternode;
    dirs [ slaveid ] = masterdir;

    vector< Node * >jdmasters;
    vector< unsigned >jddirs;
    vector< double >jdmults;
    for ( unsigned i = 0; i < nodes.size(); i++ ) {
        r = ( std :: find(otherslaves.begin(), otherslaves.end(), nodes [ i ]->giveStartingDoF() + dirs [ i ]) - otherslaves.begin() );
        if ( r == otherslaves.size() ) { //true masters
            masters.push_back(nodes [ i ]);
            directions.push_back(dirs [ i ]);
            multipliers.push_back(m [ i ]);
        } else {
            jd = constraints->giveConstraint(r);
            jdmasters = jd->giveMasterNodes();
            jddirs = jd->giveMasterDirs();
            jdmults = jd->giveMasterMultipliers();
            for ( auto &k:jdmults ) {
                k *= m [ i ];
            }
            masters.insert( masters.end(), jdmasters.begin(), jdmasters.end() );
            directions.insert( directions.end(), jddirs.begin(), jddirs.end() );
            multipliers.insert( multipliers.end(), jdmults.begin(), jdmults.end() );
        }
    }
    JointDoF :: init();
}

//////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////
// Container
//////////////////////////////////////////////////////////
void ConstraintContainer :: readFromFile(const string filename, const unsigned ndim, NodeContainer *nodecont) {
    ( void ) ndim;
    unsigned origsize = constraints.size();
    string line, ConstrType;
    ifstream inputfile( filename.c_str() );
    if ( inputfile.is_open() ) {
        while ( getline(inputfile >> std :: ws, line) ) {
            if ( line.empty() ) {
                continue;
            }
            if ( line.at(0) == '#' ) {
                continue;
            }
            istringstream iss(line);
            iss >> std :: ws >> ConstrType;
            if ( !ConstrType.rfind("#", 0) == 0 ) {
                if ( ConstrType.compare("jointDoF") == 0 ) {
                    JointDoF *newJD = new JointDoF();
                    newJD->readFromLine(iss, nodecont);
                    constraints.push_back(newJD);
                } else {
                    cerr << "Error: constraint '" <<  ConstrType <<  "' is not implemented yet." << endl;
                    exit(EXIT_FAILURE);
                }
                inputfile.close();
                cout << "Input file '" <<  filename << "' succesfully loaded; " << constraints.size() - origsize << " dependent DoFs found" << endl;
            } else {
                cerr << "Error: unable to open input file '" <<  filename <<  "'" << endl;
                exit(1);
            }
        }
    }
}

// void ConstraintContainer :: calculateSlaveDoFfield(NodeContainer *nodes){
//     for (auto &jD : constraints){
//
//     }
// }

//////////////////////////////////////////////////////////
ConstraintContainer :: ~ConstraintContainer() {
    for ( auto &c: constraints ) delete c;
}

//////////////////////////////////////////////////////////
void ConstraintContainer :: init(NodeContainer *nodecont, BCContainer *bccont) {
    //initiate volumetric averages

    unsigned numFreeDoFs = nodecont->giveTotalNumDoFs() - bccont->giveNumBlockedDoFs();

    nodes = nodecont;
    bconds = bccont;

    for ( auto const &jD : constraints ) {
        jD->init();
    }
    VolumetricAverage *va;
    for ( int k = constraints.size() - 1; k >= 0; k-- ) {
        va = dynamic_cast< VolumetricAverage * >( constraints [ k ] );
        if ( va ) {
            constraints.push_back(va);
            constraints.erase(constraints.begin() + k);
        }
    }

    map< pair< size_t, size_t >, double >indeces11;
    // map<pair<size_t, size_t>, double> indeces12;
    // // this should remain empty, unless you apply BC on Constrained DoF (Do not do it!)
    // // map<pair<size_t, size_t>, double> indeces21;
    // map<pair<size_t, size_t>, double> indeces22;
    if ( !this->isActive() ) {
        return;
    }

    ///////////////////////////////////////////////////
    // fill the matrix with corresponding multipliers for slaveDoFs
    unsigned i; // row index =  slave DoF
    unsigned j, numM; // column index = master DoF
    for ( auto const &jD : constraints ) {
        // jD->print();
        i = nodes->giveDoFid( jD->giveSlaveDoF() );
        //std::cout << jD->giveSlaveDoF() << " " << nodes->giveTotalNumDoFs() << " i = " << i << ", numFreeDoFs = " << numFreeDoFs << '\n';
        // auto res = std::find(nodes->begin(), nodes->end(), jD->giveSlaveNode());
        // std::cout << "node ID = " << std::distance(nodes->begin(), res) << '\n';
        // std::cout << "DoF num = " << jD->giveSlaveDoF() << '\n';

        if ( i < numFreeDoFs - constraints.size() ) {
            std :: cerr << "CONSTRAINT error: should never come here, constraint application unsuccesfull (hint: you are applying bondary conditions on constrained DoF) " << endl;
            cout << i << " " << jD->giveSlaveDoF() << " " <<  numFreeDoFs << " " << numFreeDoFs - constraints.size() << endl;
            exit(1);
        } else if ( i >= numFreeDoFs ) {
            // Point A = jD->giveSlaveNode()->givePoint();
            // std::cout << "node name = " << jD->giveSlaveNode()->giveName() << '\n';
            // std::cout << "Point(" << A.getX() << ", " << A.getY() << ", " << A.getZ() << ")" << '\n';
            std :: cerr << "constraint applied simultaneously with boundary conditions" << '\n';
            exit(1);
        }
        numM = jD->giveNumOfMasters();
        for ( unsigned ind = 0; ind < numM; ind++ ) {
            j = nodes->giveDoFid( jD->giveMasterDoF(ind) );
            if ( j < numFreeDoFs - constraints.size() ) {
                // master DoF is free
                indeces11.insert( pair< pair< size_t, size_t >, double >
                                      ( pair< size_t, size_t >(i, j),
                                      jD->giveMasterMultiplier(ind) ) );
            }
        }
    }

    // here fill in value 1 for all other DoFs
    ///////////////////////////////////////////////////
    for ( i = 0; i < numFreeDoFs - constraints.size(); i++ ) {
        // fill the matrix with 1 for each unrestrained DoF (diagonal)
        indeces11.insert( pair< pair< size_t, size_t >, double >(pair< size_t, size_t >(i, i), 1) );
    }
    X = CoordinateIndexedSparseMatrix( indeces11, numFreeDoFs, numFreeDoFs - constraints.size() );

    // // JK - left for testing:
    // for ( auto const &cn : this->constraints ){
    //     cn->print();
    // }
}


//////////////////////////////////////////////////////////
void ConstraintContainer :: transformToConstraintSpace(CoordinateIndexedSparseMatrix &K, const double time_now) {
    // if ( this->isTimeDependent() ) {
    //   // std::cout << "transformToConstraintSpace, time = " << time_now << '\n';
    //   this->init(this->nodes, this->bconds);
    // }
    ( void ) time_now;
    if ( X.ColumnCount > 0 ) {
        CoordinateIndexedSparseMatrix Knew;
        Knew = X.transpose() * K * X;
        K = Knew;
    } else {
        map< pair< size_t, size_t >, double >indices11;
        K = CoordinateIndexedSparseMatrix(indices11, 0, 0);
    }
}


//////////////////////////////////////////////////////////
void ConstraintContainer :: calculateDependentDoFs(Vector &fullDoFs, const double time_now, const bool all) {
    // bool all explanation: if false (by default), only multipliers are taken into account, if true,  also timeFunction-dependent parts
    if ( this->isActive() ) {
        // std::cout << "calculateDependentDoFs, time = " << time_now << '\n';
        for ( auto const &jD : constraints ) {
            fullDoFs [ jD->giveSlaveDoF() ] = 0; // to be sure that there is zero value
            for ( unsigned i = 0; i < jD->giveNumOfMasters(); i++ ) {
                fullDoFs [ jD->giveSlaveDoF() ] += fullDoFs [ jD->giveMasterDoF(i) ] * jD->giveMasterMultiplier(i) + ( all ? jD->giveFnDepPart(i, time_now) : 0.0 );
            }
        }
    }
}

//////////////////////////////////////////////////////////
void ConstraintContainer :: calculateMasterForces(Vector &fullForces) {
    if ( this->isActive() ) {
        // std::cout << "calculateMasterForces, time = " << time_now << '\n';
        for ( auto const &jD : constraints ) {
            for ( unsigned i = 0; i < jD->giveNumOfMasters(); i++ ) {
                fullForces [ jD->giveMasterDoF(i) ] += fullForces [ jD->giveSlaveDoF() ] * ( jD->giveMasterMultiplier(i) );
                // JK TODO clear fullForces[ jD->giveSlaveDoF() ] * jD->giveMultipliers()[ i ];
                // JK do not!! because they are needed for export, the same values are added to external forces, so the equilibrium is kept
                // fullForces [ jD->giveSlaveDoF() ] = 0.0;
            }
        }
    }
}


//////////////////////////////////////////////////////////
void ConstraintContainer :: removeConstraint(unsigned i) {
    if ( i > constraints.size() - 1 ) {
        cerr << "ConstraintContainer Error: requester constraint number " << i << " out of " << constraints.size() << endl;
        exit(1);
    }
    delete constraints[i];
    constraints[i] = nullptr;
    constraints.erase(constraints.begin() + i);
}



