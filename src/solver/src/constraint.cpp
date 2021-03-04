#include "constraint.h"
#include "node_container.h"
#include "element_container.h"
#include "element_container.h"
#include "boundary_condition.h"
#include <fstream>

bool containsChar(const std :: string &str, char c)
{
    // std::cout << "str = " << str << ", char = " << c << '\n';
    return str.find(c) != std :: string :: npos;
}

//////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////
// Joint Degree of Freedom
JointDoF :: JointDoF(Node *s, unsigned &dir, std :: vector< Node * > &m, std :: vector< unsigned > &dirs, std :: vector< double > &mult) {
    slaveNode = s;
    direction = dir;
    masters = m;
    directions = dirs;
    multipliers = mult;
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
        masters.push_back(nodes->giveNode(intmas) );
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
        std :: cout << "master DoF[ " << i << " ] = " << masters [ i ]->giveStartingDoF() + directions [ i ] << " with multiplier = " << multipliers [ i ] <<  '\n';
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
    otherslaves.resize(constraints->giveSize() );
    JointDoF *jd;
    for ( unsigned j = 0; j < constraints->giveSize(); j++ ) {
        jd = constraints->giveConstraint(j);
        otherslaves [ j ] = jd->giveSlaveDoF();
    }

    //find slave position in array
    unsigned slaveid = find(nodes.begin(), nodes.end(), slaveNode) - nodes.begin();
    if(slaveid==nodes.size()){
        cerr << "Volumetric average Error: slave node not included in node list" << endl;
        exit(1);
    } 

    // calculate volumes associated with nodes
    unsigned r;
    unsigned s=nodes.size();
    vector< double >m;
    m.resize(nodes.size() );
    Transp1D *et;
    RigidBodyContact *em;
    for ( unsigned e = 0; e < elems->giveSize(); e++ ) {
        et = dynamic_cast< Transp1D * >( elems->giveElement(e) );
        em = dynamic_cast< RigidBodyContact * >( elems->giveElement(e) );
        if ( et ) {
            for ( unsigned p = 0; p < 2; p++ ) {
                r = ( std :: find(nodes.begin(), nodes.end(), et->giveNode(p) ) - nodes.begin() );    
                if (r<s) m [ r ] += et->giveVolume(p);
            }
        } else if ( em )    {
            for ( unsigned p = 0; p < 2; p++ ) {
                r = ( std :: find(nodes.begin(), nodes.end(), em->giveNode(p) ) - nodes.begin() );
                if (r<s) m [ r ] += em->giveVolume(p);
            }
        }
    }

    //rearange everything and update weights
    double factor = m [ slaveid ];
    cout << slaveid << " " << factor << m.size() << endl;
    double fullVolume = 0;
    for ( auto &s: m ) {
        fullVolume += s;
        s /= -factor;
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
            masters.insert(masters.end(), jdmasters.begin(), jdmasters.end() );
            directions.insert(directions.end(), jddirs.begin(), jddirs.end() );
            multipliers.insert(multipliers.end(), jdmults.begin(), jdmults.end() );
        }
    }
    JointDoF :: init();

}

//////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////
// Container
void ConstraintContainer :: connectSlaveMaster(Node *slave, Node *master, unsigned const &ndim, const string &which, const bool &trsp) {
    unsigned nDoFsPerNode;
    if ( trsp ) {
        // if master is transport, slave must be transport
        if ( !dynamic_cast< TrsNode * >( slave ) ) {
            return;
        }
        nDoFsPerNode = 1;
    } else {
        if ( !dynamic_cast< Particle * >( slave ) ) {
            return;                                         // NOTE could be MechNode, but so far, nDoFs corresponds to Particles
        }
        nDoFsPerNode = 3 * ( ndim - 1 );
    }


    if ( slave->giveNumberOfDoFs() != master->giveNumberOfDoFs() ) {
        std :: cerr << "slave and master must have the same number of DoFs, slave numDoFs = " << slave->giveNumberOfDoFs() << ", master numDoFs = " << master->giveNumberOfDoFs() << '\n';
        exit(1);
    }

    // calculate multipliers and construct jointDof for every slaveDof
    vector< vector< double > >tableOfMultipliers;
    vector< Node * >masterNodes;
    vector< double >multipliers;
    vector< unsigned >directions;

    tableOfMultipliers.resize(nDoFsPerNode); // first index (j) is for master dof

    for ( auto &mul : tableOfMultipliers ) {
        mul.resize(nDoFsPerNode); // second index (i) is for slave dof
    }

    for ( unsigned i = 0; i < nDoFsPerNode; i++ ) {
        for ( unsigned j = 0; j < nDoFsPerNode; j++ ) {
            if ( i == j ) {
                tableOfMultipliers [ j ] [ i ] = 1;
            } else if ( i == 0 && j == nDoFsPerNode - 1 ) {
                tableOfMultipliers [ j ] [ i ] = -( slave->givePoint().getY() - master->givePoint().getY() );
            } else if ( i == 1 && j == nDoFsPerNode - 1 ) {
                tableOfMultipliers [ j ] [ i ] = slave->givePoint().getX() - master->givePoint().getX();
            } else if ( i == 0 && j == 4 ) {
                tableOfMultipliers [ j ] [ i ] = slave->givePoint().getZ() - master->givePoint().getZ();
            } else if ( i == 1 && j == 3 ) {
                tableOfMultipliers [ j ] [ i ] = -( slave->givePoint().getZ() - master->givePoint().getZ() );
            } else if ( i == 2 && j == 3 ) {
                tableOfMultipliers [ j ] [ i ] = slave->givePoint().getY() - master->givePoint().getY();
            } else if ( i == 2 && j == 4 ) {
                tableOfMultipliers [ j ] [ i ] =  -( slave->givePoint().getX() - master->givePoint().getX() );
            } else {
                tableOfMultipliers [ j ] [ i ] =  0;
            }
        }
    }
    // read the table and put the information to jointDoFs
    // NOTE this could be done in the previous loop
    for ( unsigned i = 0; i < nDoFsPerNode; i++ ) {
        // for transport nodes, only ones are in tableOfMultipliers
        if ( !trsp ) {
            if ( containsChar(which, 'x') && ( i == 0 || i == 4 || i == nDoFsPerNode - 1 ) ) {
                // std::cout << "fixed in x dir" << '\n';
            } else if ( containsChar(which, 'y') && ( i == 1 || i == 3 || i == nDoFsPerNode - 1 ) ) {
                // std::cout << "fixed in y dir" << '\n';
            } else if ( containsChar(which, 'z') && ( i == 2 || i == 3 || i == 4 ) && ndim > 2 ) {
                // std::cout << "fixed in z dir" << '\n';
            } else {
                continue;
            }
        }
        for ( unsigned j = 0; j < nDoFsPerNode; j++ ) {
            if ( tableOfMultipliers [ j ] [ i ] != 0 ) {
                masterNodes.push_back(master);
                multipliers.push_back(tableOfMultipliers [ j ] [ i ]);
                directions.push_back(j);
            }
        }
        JointDoF *newJD = new JointDoF(slave, i, masterNodes, directions, multipliers);
        constraints.push_back(newJD);
        masterNodes.clear();
        multipliers.clear();
        directions.clear();
    }
}

//////////////////////////////////////////////////////////
void ConstraintContainer :: readFromFile(const string filename, const unsigned ndim, NodeContainer *nodes) {
    unsigned origsize = constraints.size();
    string line, ConstrType;
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
            iss >> std :: ws >> ConstrType;
            if ( !ConstrType.rfind("#", 0) == 0 ) {
                if ( ConstrType.compare("jointDoF") == 0 ) {
                    JointDoF *newJD = new JointDoF();
                    newJD->readFromLine(iss, nodes);
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
void ConstraintContainer :: init(NodeContainer *nodes, BCContainer *bconds) {
    //initiate volumetric averages

    unsigned numFreeDoFs = nodes->giveTotalNumDoFs() - bconds->giveNumBlockedDoFs();


    for ( auto const &jD : constraints ) {
        jD->init();
    }

    //reorganize, move volumetric averages at the end
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
        i = nodes->giveDoFid(jD->giveSlaveDoF() );
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
            j = nodes->giveDoFid(jD->giveMasterDoF(ind) );
            if ( j < numFreeDoFs - constraints.size() ) {
                // master DoF is free
                indeces11.insert(pair< pair< size_t, size_t >, double >(pair< size_t, size_t >(i, j), jD->giveMasterMultiplier(ind) ) );
            }
        }
    }

    // here fill in value 1 for all other DoFs
    ///////////////////////////////////////////////////
    for ( i = 0; i < numFreeDoFs - constraints.size(); i++ ) {
        // fill the matrix with 1 for each unrestrained DoF (diagonal)
        indeces11.insert(pair< pair< size_t, size_t >, double >(pair< size_t, size_t >(i, i), 1) );
    }
    X = CoordinateIndexedSparseMatrix(indeces11, numFreeDoFs, numFreeDoFs - constraints.size() );
}


//////////////////////////////////////////////////////////
void ConstraintContainer :: transformToConstraintSpace(CoordinateIndexedSparseMatrix &K) {
    CoordinateIndexedSparseMatrix Kold = K;
    K = X.transpose() * Kold * X;
}


//////////////////////////////////////////////////////////
void ConstraintContainer :: calculateDependentDoFs(Vector &fullDoFs) {
    for ( auto const &jD : constraints ) {
        fullDoFs [ jD->giveSlaveDoF() ] = 0; // to be sure that there is zero value
        for ( unsigned i = 0; i < jD->giveNumOfMasters(); i++ ) {
            fullDoFs [ jD->giveSlaveDoF() ] += fullDoFs [ jD->giveMasterDoF(i) ] * jD->giveMasterMultiplier(i);
        }
    }
}

//////////////////////////////////////////////////////////
void ConstraintContainer :: calculateMasterForces(Vector &fullForces) {
    for ( auto const &jD : constraints ) {
        for ( unsigned i = 0; i < jD->giveNumOfMasters(); i++ ) {
            fullForces [ jD->giveMasterDoF(i) ] += fullForces [ jD->giveSlaveDoF() ] * jD->giveMasterMultiplier(i);
            // JK TODO clear fullForces[ jD->giveSlaveDoF() ] * jD->giveMultipliers()[ i ];
            // JK do not!! because they are needed for export, the same values are added to external forces, so the equilibrium is kept
        }
    }
}
