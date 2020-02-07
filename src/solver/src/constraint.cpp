#include "constraint.h"
#include "node_container.h"
#include <fstream>

JointDoF :: JointDoF(Node *s, unsigned &dir, std :: vector< Node * > &m, std :: vector< unsigned > &dirs, std :: vector< double > &mult) {
    slaveNode = s;
    direction = dir;
    masters = m;
    directions = dirs;
    multipliers = mult;
}

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


unsigned JointDoF :: giveSlaveDoF() const {
    return slaveNode->giveStartingDoF() + direction;
}

void JointDoF :: print() {
    std :: cout << "slave DoF = " << giveSlaveDoF() << '\n';
    for ( unsigned i = 0; i < masters.size(); i++ ) {
        std :: cout << "master DoF[ " << i << " ] = " << masters [ i ]->giveStartingDoF() + directions [ i ] << " with multiplier = " << multipliers [ i ] <<  '\n';
    }
}


bool containsChar(const std :: string &str, char c)
{
    // std::cout << "str = " << str << ", char = " << c << '\n';
    return str.find(c) != std :: string :: npos;
}

void ConstraintContainer :: connectSlaveMaster(Node *slave, Node *master, unsigned const &ndim, const string &which) {
    unsigned nDoFsPerNode = 3 * ( ndim - 1 );

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
        if ( containsChar(which, 'x') && ( i == 0 || i == 4 || i == nDoFsPerNode - 1 ) ) {
            // std::cout << "fixed in x dir" << '\n';
        } else if ( containsChar(which, 'y') && ( i == 1 || i == 3 || i == nDoFsPerNode - 1 ) ) {
            // std::cout << "fixed in y dir" << '\n';
        } else if ( containsChar(which, 'z') && ( i == 2 || i == 3 || i == 4 ) && ndim > 2 ) {
            // std::cout << "fixed in z dir" << '\n';
        } else {
            continue;
        }
        for ( unsigned j = 0; j < nDoFsPerNode; j++ ) {
            if ( tableOfMultipliers [ j ] [ i ] != 0 ) {
                // std::cout << "master->giveStartingDoF() = " << master->giveStartingDoF() << '\n';
                // THE PROBLEM IS IN SETTING THE STARTING DOF FOR MASTER
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


void ConstraintContainer :: readRigidPlate(istringstream &iss, const unsigned ndim, NodeContainer *nodes) {
    // jointDoF jD;
    unsigned nodeid, nslaves;
    Node *master;
    Node *slave;
    //////////////////////////////////////////////////////////
    // read the line "masterId numSlaves slaveId1, slaveId2...."
    iss >> nodeid >> nslaves;
    master = nodes->giveNode(nodeid);
    // check if it is master node
    MasterNode *n = dynamic_cast< MasterNode * >( master );
    if ( !n ) {
        cerr << "Error in " << __func__ << ": node must be MasterDoF, " << master->giveName() << " provided" << endl;
        exit(1);
    }
    if ( n->giveNumberOfDoFs() != ( 3 * ( ndim - 1 ) ) ) {
        cerr << "Error in " << __func__ << ": MasterDoF for RigidPlate must have " << ( 3 * ( ndim - 1 ) ) << " DoFs, " << n->giveNumberOfDoFs() << " provided" << endl;
        exit(1);
    }

    for ( unsigned i = 0; i < nslaves; i++ ) {
        iss >> nodeid;
        slave = nodes->giveNode(nodeid);
        connectSlaveMaster(slave, master, ndim, "xyz");
    }
}


bool isInBlock(const Point &P, const Point &leftBottom, const Point &rightTop) {
    double tol = 1e-9; // NOTE this should be raltive to lmin or something
    if ( ( leftBottom.getX() - P.getX() ) < tol && ( P.getX() - rightTop.getX() ) < tol ) {
        if ( ( leftBottom.getY() - P.getY() ) < tol && ( P.getY() - rightTop.getY() ) < tol ) {
            if ( ( leftBottom.getZ() - P.getZ() ) < tol && ( P.getZ() - rightTop.getZ() ) < tol ) {
                return true;
            }
        }
    }
    return false;
}


void ConstraintContainer :: readCoordRigidPlate(istringstream &iss, const unsigned ndim, NodeContainer *nodes) {
    // jointDoF jD;
    unsigned nodeid;
    double x0, x1, y0, y1, z0, z1;
    Node *master;
    Point leftBottom, rightTop;
    //////////////////////////////////////////////////////////
    // read the line "masterId numSlaves slaveId1, slaveId2...."
    if ( ndim == 2 ) {
        iss >> nodeid >> x0 >> x1 >> y0 >> y1;
        leftBottom = Point(x0, y0, 0);
        rightTop = Point(x1, y1, 0);
    } else if ( ndim == 3 ) {
        iss >> nodeid >> x0 >> x1 >> y0 >> y1 >> z0 >> z1;
        leftBottom = Point(x0, y0, z0);
        rightTop = Point(x1, y1, z1);
    } else {
        std :: cerr << "dimension " << ndim << " not implemented yet" << '\n';
        exit(1);
    }
    master = nodes->giveNode(nodeid);
    // check if it is master node
    MasterNode *n = dynamic_cast< MasterNode * >( master );
    if ( !n ) {
        cerr << "Error in " << __func__ << ": node must be MasterDoF, " << master->giveName() << " provided" << endl;
        exit(1);
    }
    if ( n->giveNumberOfDoFs() != ( 3 * ( ndim - 1 ) ) ) {
        cerr << "Error in " << __func__ << ": MasterDoF for RigidPlate must have " << ( 3 * ( ndim - 1 ) ) << " DoFs, " << n->giveNumberOfDoFs() << " provided" << endl;
        exit(1);
    }

    for ( auto const &nod : * nodes ) {
        if ( isInBlock(nod->givePoint(), leftBottom, rightTop) ) {
            // NOTE this is quite unefficient, could be done checking num of DoFs (...?)
            Particle *nn = dynamic_cast< Particle * >( nod );
            if ( nn ) {
                connectSlaveMaster(nod, master, ndim, "xyz");
            }
        }
    }
}

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
                } else if ( ConstrType.compare("RigidPlate") == 0 ) {
                    readRigidPlate(iss, ndim, nodes);
                } else if ( ConstrType.compare("CoordRigidPlate") == 0 ) {
                    readCoordRigidPlate(iss, ndim, nodes);
                } else {
                    cerr << "Error: constraint '" <<  ConstrType <<  "' is not implemented yet." << endl;
                    exit(0);
                }
            }
        }
        inputfile.close();
        cout << "Input file '" <<  filename << "' succesfully loaded; " << constraints.size() - origsize << " dependent DoFs found" << endl;
    } else {
        cerr << "Error: unable to open input file '" <<  filename <<  "'" << endl;
        exit(1);
    }
}

// void ConstraintContainer :: calculateSlaveDoFfield(NodeContainer *nodes){
//     for (auto &jD : constraints){
//
//     }
// }


void ConstraintContainer :: init(NodeContainer *nodes) {
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
        //std::cout << jD->giveSlaveDoF() << " " << nodes->giveTotalNumDoFs() << " i = " << i << ", giveNumFreeDoFs() = " << nodes->giveNumFreeDoFs() << '\n';
        // auto res = std::find(nodes->begin(), nodes->end(), jD->giveSlaveNode());
        // std::cout << "node ID = " << std::distance(nodes->begin(), res) << '\n';
        // std::cout << "DoF num = " << jD->giveSlaveDoF() << '\n';
        if ( i < nodes->giveNumFreeDoFs() - constraints.size() ) {
            std :: cerr << "should never come here, constraint application unsuccesfull (hint: you are applying bondary conditions on constrained DoF) " << '\n';
            exit(1);
        } else if ( i >= nodes->giveNumFreeDoFs() ) {
            // Point A = jD->giveSlaveNode()->givePoint();
            // std::cout << "node name = " << jD->giveSlaveNode()->giveName() << '\n';
            // std::cout << "Point(" << A.getX() << ", " << A.getY() << ", " << A.getZ() << ")" << '\n';
            std :: cerr << "constraint applied simultaneously with boundary conditions" << '\n';
            exit(1);
        }
        numM = jD->giveMasterDoFs().size();
        for ( unsigned ind = 0; ind < numM; ind++ ) {
            j = nodes->giveDoFid( ( jD->giveMasterDoFs() [ ind ] )->giveStartingDoF() + jD->giveDirs() [ ind ] );
            if ( j < nodes->giveNumFreeDoFs() - constraints.size() ) {
                // master DoF is free
                indeces11.insert(pair< pair< size_t, size_t >, double >(pair< size_t, size_t >(i, j), jD->giveMultipliers() [ ind ]) );
            }
        }
    }

    // here fill in value 1 for all other DoFs
    ///////////////////////////////////////////////////
    for ( i = 0; i < nodes->giveNumFreeDoFs() - constraints.size(); i++ ) {
        // fill the matrix with 1 for each unrestrained DoF (diagonal)
        indeces11.insert(pair< pair< size_t, size_t >, double >(pair< size_t, size_t >(i, i), 1) );
    }
    X = CoordinateIndexedSparseMatrix(indeces11, nodes->giveNumFreeDoFs(), nodes->giveNumFreeDoFs() - constraints.size() );
    // X.print();
}


void ConstraintContainer :: transformToConstraintSpace(CoordinateIndexedSparseMatrix &K) {
    CoordinateIndexedSparseMatrix Kold = K;
    K = X.transpose() * Kold * X;
}


void ConstraintContainer :: calculateDependentDoFs(Vector &fullDoFs) {
    for ( auto const &jD : constraints ) {
        fullDoFs [ jD->giveSlaveDoF() ] = 0; // to be sure that there is zero value
        for ( unsigned i = 0; i < jD->giveMasterDoFs().size(); i++ ) {
            fullDoFs [ jD->giveSlaveDoF() ] += fullDoFs [ jD->giveMasterDoFs() [ i ]->giveStartingDoF() + jD->giveDirs() [ i ] ] * jD->giveMultipliers() [ i ];
        }
    }
}

void ConstraintContainer :: calculateMasterForces(Vector &fullForces) {
    for ( auto const &jD : constraints ) {
        for ( unsigned i = 0; i < jD->giveMasterDoFs().size(); i++ ) {
            fullForces [ jD->giveMasterDoFs() [ i ]->giveStartingDoF() + jD->giveDirs() [ i ] ] += fullForces [ jD->giveSlaveDoF() ] * jD->giveMultipliers() [ i ];
            // JK TODO clear fullForces[ jD->giveSlaveDoF() ] * jD->giveMultipliers()[ i ];
            // JK do not!! because they are needed for export, the same values are added to external forces, so the equilibrium is kept
        }
    }
}
