#include "constraint.h"
#include "node_container.h"
#include "element_container.h"
#include "element_discrete.h"
#include "boundary_condition.h"
#include "solver.h"
#include "solver_implicit.h"
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
void JointDoF :: readFromLine(std :: istringstream &iss, NodeContainer *nodes) {
    // # type      node    direction   numMasters master1 multiplier1 master2...
    unsigned intn, intmas, dir;
    double mult;
    iss >> intn;
    slaveNode = nodes->giveNode(intn);
    iss >> direction;
    iss >> intn;
    for ( unsigned i = 0; i < intn; i++ ) {
        iss >> intmas >> dir >> mult;
        masters.push_back(nodes->giveNode(intmas) );
        directions.push_back(dir);
        time_fns.push_back(nullptr);
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
void JointDoF :: init(Solver *solver) {
    ( void ) solver;
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
VolumetricAverage :: VolumetricAverage(std :: vector< Node * > &n, std :: vector< unsigned > &d, Node *mn, unsigned md, ElementContainer *ec, ConstraintContainer *cc) {
    nodes = n;
    dirs = d;
    masternode = mn;
    masterdir = md;
    elems = ec;
    constraints = cc;


    //collect all slaves from other joint DoFs and also all involved DOFs
    std :: vector< Node * >excludedNodes;
    std :: vector< unsigned >excudedDirs;
    std :: vector< Node * >jdm;
    std :: vector< unsigned >jdd;
    JointDoF *jd;
    VolumetricAverage *va;
    for ( unsigned j = 0; j < constraints->giveSize(); j++ ) {
        jd = constraints->giveConstraint(j);
        va = dynamic_cast< VolumetricAverage * >( jd );
        if ( !va ) {
            excludedNodes.push_back(jd->giveSlaveNode() );
            jdm = jd->giveMasterNodes();
            excludedNodes.insert(excludedNodes.end(), jdm.begin(), jdm.end() );
            excudedDirs.push_back(jd->giveSlaveDir() );
            jdd = jd->giveMasterDirs();
            excudedDirs.insert(excudedDirs.end(), jdd.begin(), jdd.end() );
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
        std :: cerr << "volumetric average didn't find any suitable slave DoF" << std :: endl;
        exit(1);
    }
    slaveNode = nodes [ k ];
    direction = dirs [ k ];
}

//////////////////////////////////////////////////////////
void VolumetricAverage :: init(Solver *solver) {
    //collect all slaves from other joint DoFs
    std :: vector< unsigned >otherslaves;
    otherslaves.resize(constraints->giveSize() );
    JointDoF *jd;
    for ( unsigned j = 0; j < constraints->giveSize(); j++ ) {
        jd = constraints->giveConstraint(j);
        otherslaves [ j ] = jd->giveSlaveDoF();
    }

    //find slave position in array
    unsigned slaveid = find(nodes.begin(), nodes.end(), slaveNode) - nodes.begin();
    if ( slaveid == nodes.size() ) {
        std :: cerr << "Volumetric average Error: slave node not included in node list" << std :: endl;
        exit(1);
    }

    // calculate volumes associated with nodes
    unsigned r;
    unsigned s = nodes.size();
    std :: vector< double >m;
    m.resize(nodes.size() );
    DiscreteTrsprtElem *et;
    RigidBodyContact *em;
    for ( unsigned e = 0; e < elems->giveSize(); e++ ) {
        et = dynamic_cast< DiscreteTrsprtElem * >( elems->giveElement(e) );
        em = dynamic_cast< RigidBodyContact * >( elems->giveElement(e) );
        if ( et ) {
            for ( unsigned p = 0; p < 2; p++ ) {
                r = ( std :: find(nodes.begin(), nodes.end(), et->giveNode(p) ) - nodes.begin() );
                if ( r < s ) {
                    m [ r ] += et->giveVolumeAssociatedWithNode(p);
                }
            }
        } else if ( em ) {
            for ( unsigned p = 0; p < 2; p++ ) {
                r = ( std :: find(nodes.begin(), nodes.end(), em->giveNode(p) ) - nodes.begin() );
                if ( r < s ) {
                    m [ r ] += em->giveVolumeAssociatedWithNode(p);
                }
            }
        }
    }

    //rearange everything and update weights
    double factor = m [ slaveid ];
    std :: cout << slaveid << " " << factor << m.size() << std :: endl;
    double fullVolume = 0;
    for ( auto &ss: m ) {
        fullVolume += ss;
        ss /= -factor;
    }
    std :: cout << "volumetric Average: check of volume " << fullVolume << std :: endl;
    m [ slaveid ] = fullVolume / factor;
    nodes [ slaveid ] = masternode;
    dirs [ slaveid ] = masterdir;

    std :: vector< Node * >jdmasters;
    std :: vector< unsigned >jddirs;
    std :: vector< double >jdmults;
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
            masters.insert(masters.end(), jdmasters.begin(), jdmasters.end() );
            directions.insert(directions.end(), jddirs.begin(), jddirs.end() );
            multipliers.insert(multipliers.end(), jdmults.begin(), jdmults.end() );
        }
    }
    JointDoF :: init(solver);
}


//////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////
// BC DEPENDENT ON OTHER DOFS OR CONJUGATE VARIABLES
DoFDependentOnConjugates :: DoFDependentOnConjugates(Node *s, const unsigned &dir, const std :: vector< Node * > &m, const std :: vector< unsigned > &dirs, const std :: vector< double > &mult, const std :: vector< Function * > &fns, const std :: vector< double > &time_mult) : JointDoF(s, dir, m, dirs, mult, fns, time_mult) {}


//////////////////////////////////////////////////////////
void DoFDependentOnConjugates :: init(Solver *solver) {
    JointDoF :: init(solver);
    TransientLinearTransportSolver *ts = dynamic_cast< TransientLinearTransportSolver * >( solver );
    if ( ts ) {
        std :: cerr << "*******************************" << std :: endl;
        std :: cerr << "Warning: using DoFDependentOnConjugates is strongly discouraged as the time integration might not take place at the end of time step so the update of primary variables at the end will be incorrect" << std :: endl;
        std :: cerr << "*******************************" << std :: endl;
    }
}

//////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////
// Container
//////////////////////////////////////////////////////////
void ConstraintContainer :: readFromFile(const std :: string filename, const unsigned ndim, NodeContainer *nodecont) {
    ( void ) ndim;
    unsigned origsize = constraints.size();
    std :: string line, ConstrType;
    std :: ifstream inputfile(filename.c_str() );
    if ( inputfile.is_open() ) {
        while ( getline(inputfile >> std :: ws, line) ) {
            if ( line.empty() || ( line.at(0) == '#' ) ) {
                continue;
            }
            std :: istringstream iss(line);
            iss >> std :: ws >> ConstrType;
            if ( !( ConstrType.rfind("#", 0) == 0 ) ) {
                if ( ConstrType.compare("jointDoF") == 0 ) {
                    JointDoF *newJD = new JointDoF();
                    newJD->readFromLine(iss, nodecont);
                    constraints.push_back(newJD);
                } else {
                    std :: cerr << "Error: constraint '" <<  ConstrType <<  "' is not implemented yet." << std :: endl;
                    exit(EXIT_FAILURE);
                }
            }
        }
        std :: cout << "Input file '" <<  filename << "' succesfully loaded; " << constraints.size() - origsize << " dependent DoFs found" << std :: endl;
        inputfile.close();
    } else {
        std :: cerr << "Error: unable to open input file '" <<  filename <<  "'" << std :: endl;
        exit(1);
    }
}

// void ConstraintContainer :: calculateSlaveDoFfield(NodeContainer *nodes){
//     for (auto &jD : constraints){
//
//     }
// }

//////////////////////////////////////////////////////////
ConstraintContainer :: ~ConstraintContainer() {
    for ( auto &c: constraints ) {
        if ( c != nullptr ) {
            delete c;
        }
    }
}


void ConstraintContainer :: clear() {
    for ( auto &c: constraints ) {
        if ( c != nullptr ) {
            delete c;
        }
    }
}

//////////////////////////////////////////////////////////
void ConstraintContainer :: init(NodeContainer *nodecont, BCContainer *bccont, Solver *solver) {
    //initiate volumetric averages

    unsigned numFreeDoFs = nodecont->giveTotalNumDoFs() - bccont->giveNumBlockedDoFs();

    nodes = nodecont;
    bconds = bccont;

    for ( auto const &jD : constraints ) {
        jD->init(solver);
    }
    VolumetricAverage *va;
    for ( int k = constraints.size() - 1; k >= 0; k-- ) {
        va = dynamic_cast< VolumetricAverage * >( constraints [ k ] );
        if ( va ) {
            constraints.push_back(va);
            constraints.erase(constraints.begin() + k);
        }
    }

    std :: vector< Ttripletd >tripletList;
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
        i = nodes->giveDoFid(jD->giveSlaveDoF() );
        //std::cout << jD->giveSlaveDoF() << " " << nodes->giveTotalNumDoFs() << " i = " << i << ", numFreeDoFs = " << numFreeDoFs << '\n';
        // auto res = std::find(nodes->begin(), nodes->end(), jD->giveSlaveNode());
        // std::cout << "node ID = " << std::distance(nodes->begin(), res) << '\n';
        // std::cout << "DoF num = " << jD->giveSlaveDoF() << '\n';

        if ( i < numFreeDoFs - constraints.size() ) {
            std :: cerr << "CONSTRAINT error: should never come here, constraint application unsuccesfull (hint: you are applying bondary conditions on constrained DoF) " << std :: endl;
            std :: cout << i << " " << jD->giveSlaveDoF() << " " <<  numFreeDoFs << " " << numFreeDoFs - constraints.size() << std :: endl;
            exit(1);
        } else if ( i >= numFreeDoFs ) {
            // Point A = jD->giveSlaveNode()->givePoint();
            // std::cout << "node name = " << jD->giveSlaveNode()->giveName() << '\n';
            // std::cout << "Point(" << A.getX() << ", " << A.getY() << ", " << A.getZ() << ")" << '\n';
            std :: cerr << "constraint applied simultaneously with boundary conditions" << '\n';
            exit(1);
        }
        numM = jD->giveNumOfDoFMasters();
        for ( unsigned ind = 0; ind < numM; ind++ ) {
            j = nodes->giveDoFid(jD->giveMasterDoF(ind) );
            if ( j < numFreeDoFs - constraints.size() ) {
                // master DoF is free
                tripletList.push_back(Ttripletd(i, j, jD->giveMasterMultiplier(ind) ) );
            }
        }
    }

    // here fill in value 1 for all other DoFs
    ///////////////////////////////////////////////////
    for ( i = 0; i < numFreeDoFs - constraints.size(); i++ ) {
        // fill the matrix with 1 for each unrestrained DoF (diagonal)
        tripletList.push_back(Ttripletd(i, i, 1) );
    }

    X.resize(numFreeDoFs, numFreeDoFs - constraints.size() );
    X.setFromTriplets(tripletList.begin(), tripletList.end() );
    X.makeCompressed();
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
    // TODO: its possible that it will work
    CoordinateIndexedSparseMatrix Knew;
    Knew = X.transpose() * K * X;
    K = Knew;
    /*
     * if ( X.cols() > 0 ) {
     *  CoordinateIndexedSparseMatrix Knew;
     *  Knew = X.transpose() * K * X;
     *  K = Knew;
     * } else {
     *  //map< pair< size_t, size_t >, double >indices11;
     *  //K = CoordinateIndexedSparseMatrix(indices11, 0, 0);
     *  //Ttripletd tripletList;
     *  //K.setFromTriplets(tripletList.begin(), tripletList.end() );
     *  //K.makeCompressed();
     * }*/
}

//////////////////////////////////////////////////////////
void ConstraintContainer :: calculateDependentDoFs(Vector &fullDoFs, const double time_now, const bool all) const {
    // bool all explanation: if false (by default), only multipliers are taken into account, if true,  also timeFunction-dependent parts
    if ( this->isActive() ) {
        for ( auto const &jD : constraints ) {
            fullDoFs [ jD->giveSlaveDoF() ] = 0; // to be sure that there is zero value
            //standard constraint, imposed relations between dofs
            for ( unsigned i = 0; i < jD->giveNumOfDoFMasters(); i++ ) {
                fullDoFs [ jD->giveSlaveDoF() ] += fullDoFs [ jD->giveMasterDoF(i) ] * jD->giveMasterMultiplier(i) + ( all ? jD->giveFnDepPart(i, time_now) : 0.0 );
            }
        }
    }
}

//////////////////////////////////////////////////////////
void ConstraintContainer :: calculateDoFsDependentOnConjugates(Vector &full_ddr, const Vector &trial_r, const Vector &fullFExt) const {
    // bool all explanation: if false (by default), only multipliers are taken into account, if true,  also timeFunction-dependent parts
    if ( this->isActive() ) {
        for ( auto const &jD : constraints ) {
            //full_ddr [ jD->giveSlaveDoF() ] = 0; // to be sure that there is zero value
            //constraint relating DoFs to conjugate variables
            for ( unsigned i = 0; i < jD->giveNumOfConjugateMasters(); i++ ) {
                full_ddr [ jD->giveSlaveDoF() ] += fullFExt [ jD->giveMasterDoF(i) ] * jD->giveMasterMultiplier(i);
            }
            if ( jD->giveNumOfConjugateMasters() > 0 ) {
                full_ddr [ jD->giveSlaveDoF() ] -= trial_r [ jD->giveSlaveDoF() ];   //this is needed because ddr will be later summed with trial_r, but the calculated pressures are already absolute (after the summation)
            }
        }
    }
}

//////////////////////////////////////////////////////////
void ConstraintContainer :: calculateMasterForces(Vector &fullForces) {
    if ( this->isActive() ) {
        // std::cout << "calculateMasterForces, time = " << time_now << '\n';
        for ( auto const &jD : constraints ) {
            for ( unsigned i = 0; i < jD->giveNumOfDoFMasters(); i++ ) {
                fullForces [ jD->giveMasterDoF(i) ] += fullForces [ jD->giveSlaveDoF() ] * ( jD->giveMasterMultiplier(i) );
            }
            fullForces [ jD->giveSlaveDoF() ] = 0;
        }
    }
}


//////////////////////////////////////////////////////////
void ConstraintContainer :: removeConstraint(unsigned i) {
    if ( i > constraints.size() - 1 ) {
        std :: cerr << "ConstraintContainer Error: requester constraint number " << i << " out of " << constraints.size() << std :: endl;
        exit(1);
    }
    delete constraints [ i ];
    constraints [ i ] = nullptr;
    constraints.erase(constraints.begin() + i);
}
