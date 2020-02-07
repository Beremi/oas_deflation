#include "preprocessing_block.h"

//////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////
// Basic Periodic BC
void BasicPeriodicBC :: apply(NodeContainer *nodes, ElementContainer *e, BCContainer *bcs, ConstraintContainer *constrs, FunctionContainer *funcs, ExporterContainer *ex) {
    ( void ) e;
    volume = 1;
    for ( auto const a : PUCsize ) {
        volume *= a;
    }

    unsigned const_num = constrs->giveSize();
    unsigned funcs_num = funcs->giveSize();
    unsigned bcs_num = bcs->giveSize();
    unsigned ex_num = ex->giveSize();

    //create new degrees of freedom representing strains ex, ey, gammaxy=2exy or ex, ey, ez, gammyz, gammaxz, gammaxy,
    MasterDoF *mn;
    unsigned intialNodeNum = nodes->giveSize();
    for ( unsigned i = 0; i < 3 * ( dim - 1 ); i++ ) {
        mn = new MasterDoF(Point(0, 0, 0), 1);
        nodes->addNode(mn);
    }

    //apply contraints, connect periodic images
    JointDoF *jd;
    vector< Node * >vm;
    vector< unsigned >dirs;
    vector< double >mults;
    Node *s, *m;
    Point diff;
    for ( unsigned i = 0; i < masters.size(); i++ ) {
        m = nodes->giveNode(masters [ i ]);
        s = nodes->giveNode(slaves [ i ]);
        //connect rotations
        if ( dynamic_cast< Particle * >( s ) && dynamic_cast< Particle * >( m ) ) {
            dirs.resize(1);
            mults.resize(1);
            vm.resize(1);
            mults [ 0 ] = 1;
            vm [ 0 ] = m;
            for ( unsigned k = 0; k < 2 * ( dim - 1 ) - 1; k++ ) {
                dirs [ 0 ] = dim + k;
                jd = new JointDoF(s, dirs [ 0 ], vm, dirs, mults);
                constrs->addConstraint(jd);
            }
        }

        //connect translations
        diff = s->givePoint() - m->givePoint();

        if ( use_half_gammas ) {
            //direction X  (all gammaxy and gammaxy realized here)
            if ( dim == 3 ) {
                vm.resize(4);
                mults.resize(4);
                dirs.resize(4, 0);
                dirs [ 2 ] = 0;
                dirs [ 3 ] = 0;
                vm [ 2 ] = nodes->giveNode(intialNodeNum + 5); //gamma xy
                vm [ 3 ] = nodes->giveNode(intialNodeNum + 4); //gamma xz
                mults [ 3 ] = diff.z / 2;
            } else if ( dim == 2 ) {
                vm.resize(3);
                mults.resize(3);
                dirs.resize(3, 0);
                dirs [ 2 ] = 0;
                vm [ 2 ] = nodes->giveNode(intialNodeNum + 2); //gamma xy
            }
            dirs [ 0 ] = 0;
            dirs [ 1 ] = 0;
            vm [ 0 ] = m; //master
            vm [ 1 ] = nodes->giveNode(intialNodeNum); //eps x
            mults [ 0 ] = 1;
            mults [ 1 ] = diff.x / 2;
            mults [ 2 ] = diff.y / 2;
            jd = new JointDoF(s, dirs [ 0 ], vm, dirs, mults);
            constrs->addConstraint(jd);

            //direction Y  (gammaxz realized here)
            if ( dim == 3 ) {
                vm.resize(4);
                mults.resize(4);
                dirs.resize(4, 0);
                dirs [ 2 ] = 0;
                dirs [ 3 ] = 0;
                vm [ 2 ] = nodes->giveNode(intialNodeNum + 5); //gamma xy
                vm [ 3 ] = nodes->giveNode(intialNodeNum + 3); //gamma yz
                mults [ 3 ] = diff.z / 2;
            } else if ( dim == 2 ) {
                vm.resize(3);
                mults.resize(3);
                dirs.resize(3, 0);
                dirs [ 2 ] = 0;
                vm [ 2 ] = nodes->giveNode(intialNodeNum + 2); //gamma xy
            }
            dirs [ 0 ] = 1;
            dirs [ 1 ] = 0;
            vm [ 0 ] = m; //master
            vm [ 1 ] = nodes->giveNode(intialNodeNum + 1); //eps y
            mults [ 0 ] = 1;
            mults [ 1 ] = diff.y;
            mults [ 2 ] = diff.x / 2;
            jd = new JointDoF(s, dirs [ 0 ], vm, dirs, mults);
            constrs->addConstraint(jd);

            //direction Z  (gammaxz realized here)
            if ( dim == 3 ) {
                vm.resize(4);
                mults.resize(4);
                dirs.resize(4, 0);
                dirs [ 2 ] = 0;
                dirs [ 3 ] = 0;
                vm [ 2 ] = nodes->giveNode(intialNodeNum + 4); //gamma xy
                vm [ 3 ] = nodes->giveNode(intialNodeNum + 3); //gamma yz
                mults [ 3 ] = diff.y / 2;

                dirs [ 0 ] = 2;
                dirs [ 1 ] = 0;
                vm [ 0 ] = m; //master
                vm [ 1 ] = nodes->giveNode(intialNodeNum + 2); //eps y
                mults [ 0 ] = 1;
                mults [ 1 ] = diff.z;
                mults [ 2 ] = diff.x / 2;
                jd = new JointDoF(s, dirs [ 0 ], vm, dirs, mults);
                constrs->addConstraint(jd);
            }
        } else {
            //direction X  (all gammaxy and gammaxy realized here)
            if ( dim == 3 ) {
                vm.resize(4);
                mults.resize(4);
                dirs.resize(4, 0);
                dirs [ 2 ] = 0;
                dirs [ 3 ] = 0;
                vm [ 2 ] = nodes->giveNode(intialNodeNum + 5); //gamma xy
                vm [ 3 ] = nodes->giveNode(intialNodeNum + 4); //gamma xz
                mults [ 3 ] = diff.z;
            } else if ( dim == 2 ) {
                vm.resize(3);
                mults.resize(3);
                dirs.resize(3, 0);
                dirs [ 2 ] = 0;
                vm [ 2 ] = nodes->giveNode(intialNodeNum + 2); //gamma xy
            }
            dirs [ 0 ] = 0;
            dirs [ 1 ] = 0;
            vm [ 0 ] = m; //master
            vm [ 1 ] = nodes->giveNode(intialNodeNum); //eps x
            mults [ 0 ] = 1;
            mults [ 1 ] = diff.x;
            mults [ 2 ] = diff.y;
            jd = new JointDoF(s, dirs [ 0 ], vm, dirs, mults);
            constrs->addConstraint(jd);

            //direction Y  (gammaxz realized here)
            if ( dim == 3 ) {
                vm.resize(3);
                mults.resize(3);
                dirs.resize(3, 0);
                vm [ 2 ] = nodes->giveNode(intialNodeNum + 3); //gamma yz
                mults [ 2 ] = diff.z;
            } else if ( dim == 2 ) {
                vm.resize(2);
                mults.resize(2);
                dirs.resize(2, 0);
            }
            dirs [ 0 ] = 1;
            vm [ 0 ] = m; //master
            vm [ 1 ] = nodes->giveNode(intialNodeNum + 1); //eps y
            mults [ 0 ] = 1;
            mults [ 1 ] = diff.y;
            dirs [ 0 ] = 1;  // JK: here should be dirs[1] = 0; shouldn't it...? probably stays form dir X, so it does not matter, but dirs[0] is aleready  defined 5 lines upper
            jd = new JointDoF(s, dirs [ 0 ], vm, dirs, mults);
            constrs->addConstraint(jd);

            //direction Z  (gammaxz realized here)
            if ( dim == 3 ) {
                vm.resize(2);
                mults.resize(2);
                dirs.resize(2, 0);
                dirs [ 0 ] = 2;
                vm [ 0 ] = m; //master
                vm [ 1 ] = nodes->giveNode(intialNodeNum + 2); //eps z
                mults [ 0 ] = 1;
                mults [ 1 ] = diff.z;
                dirs [ 0 ] = 2;
                jd = new JointDoF(s, dirs [ 0 ], vm, dirs, mults);
                constrs->addConstraint(jd);
            }
        }
    }


    //boundary conditions
    //last master node cannot move
    BoundaryCondition *bc;
    vector< int >dBC, nBC;
    dBC.resize(m->giveNumberOfDoFs(), -1);
    nBC.resize(m->giveNumberOfDoFs(), -1);
    for ( unsigned i = 0; i < dim; i++ ) {
        dBC [ i ] = funcs->giveSize();                        //only translation
    }
    bc = new BoundaryCondition(m, dBC, nBC);
    bcs->addBoundaryCondition(bc);

    //set prescribed strain and stress
    vector< double >bcmults;
    dBC.resize(1, -1);
    nBC.resize(1, -1);
    bcmults.resize(1, 1);
    for ( unsigned i = 0; i < strainFunc.size(); i++ ) {
        if ( strainFunc [ i ] >= 0 ) {
            nBC [ 0 ] = -1;
            dBC [ 0 ] = strainFunc [ i ];
            bc = new BoundaryCondition(nodes->giveNode(intialNodeNum + i), dBC, nBC, bcmults);
            bcs->addBoundaryCondition(bc);
        }
        bcmults [ 0 ] = volume;
        if ( stressFunc [ i ] >= 0 ) {
            nBC [ 0 ] = stressFunc [ i ];
            dBC [ 0 ] = -1;
            bc = new BoundaryCondition(nodes->giveNode(intialNodeNum + i), dBC, nBC, bcmults);
            bcs->addBoundaryCondition(bc);
        }
        if ( strainFunc [ i ] >= 0 && stressFunc [ i ] >= 0 ) {
            cerr << "Error in Periodic boundary condition: cannot prescribe both stress and strain for the same direction" << endl;
        }
    }


    //add constant function
    vector< double >x, y;
    x.resize(1, 0);
    y.resize(1, 0);
    PieceWiseLinearFunction *newf = new PieceWiseLinearFunction(x, y);
    funcs->addFunction(newf);

    //export data
    string name = "PUCstrain_stress";
    vector< string >gname;
    vector< unsigned >n;
    vector< string >codes;
    n.resize(1);
    codes.resize(1);
    ForceGauge *fg;
    codes [ 0 ] = "fx";
    if ( dim == 2 ) {
        gname.resize(3);
        gname [ 0 ] = "sigma_x";
        gname [ 1 ] = "sigma_y";
        gname [ 2 ] = "tau_xy";
    } else if ( dim == 3 ) {
        gname.resize(3);
        gname [ 0 ] = "sigma_x";
        gname [ 1 ] = "sigma_y";
        gname [ 2 ] = "sigma_z";
        gname [ 3 ] = "tau_yz";
        gname [ 4 ] = "tau_xz";
        gname [ 5 ] = "tau_xy";
    }
    for ( unsigned i = 0; i < gname.size(); i++ ) {
        n [ 0 ] = intialNodeNum + i;
        fg = new ForceGauge(name, gname [ i ], codes, n, nodes, 1. / volume, dim);
        ex->addExporter(fg);
    }
    DoFGauge *dg;
    if ( dim == 2 ) {
        gname [ 0 ] = "eps_x";
        gname [ 1 ] = "eps_y";
        gname [ 2 ] = "gamma_xy";
    } else if ( dim == 3 ) {
        gname [ 0 ] = "eps_x";
        gname [ 1 ] = "eps_y";
        gname [ 2 ] = "eps_z";
        gname [ 3 ] = "gamma_yz";
        gname [ 4 ] = "gamma_xz";
        gname [ 5 ] = "gamma_xy";
    }
    for ( unsigned i = 0; i < gname.size(); i++ ) {
        dg = new DoFGauge(name, gname [ i ], intialNodeNum + i, 0, nodes, 1., dim);
        ex->addExporter(dg);
    }



    cout << "Applied periodic boundary conditions: " << nodes->giveSize() - intialNodeNum << " new DoFs (nodes " << intialNodeNum << " - " <<  nodes->giveSize() - 1 << "); " << constrs->giveSize() - const_num << " new constraints; " << bcs->giveSize() - bcs_num << " new boundary conditions; " << funcs->giveSize() - funcs_num << " new function; " << ex->giveSize() - ex_num << " new exporters; " << "created" << endl;
}

//////////////////////////////////////////////////////////
void BasicPeriodicBC :: readFromLine(istringstream &iss, unsigned d) {
    use_half_gammas = false;
    dim = d;
    string param;
    unsigned num, hnum;

    bool sizeB, loadB, pairsB;
    sizeB = loadB = pairsB = false;

    strainFunc.resize(3 * ( dim - 1 ), -1);
    stressFunc.resize(3 * ( dim - 1 ), -1);

    while ( !iss.eof() ) {
        iss >> param;
        if ( param.compare("pairs") == 0 ) {
            pairsB = true;
            iss >> num;
            masters.resize(num);
            slaves.resize(num);
            for ( unsigned i = 0; i < num; i++ ) {
                iss >> slaves [ i ] >> masters [ i ];
            }
        } else if ( param.compare("use_half_gammas") == 0 ) {
            use_half_gammas = true;
        } else if ( param.compare("load") == 0 ) {
            loadB = true;
            iss >> num;
            for ( unsigned i = 0; i < num; i++ ) {
                iss >> param >> hnum;
                if ( dim == 2 ) {
                    std :: size_t found = param.find("z");
                    if ( found != std :: string :: npos ) {
                        cout << "Error in PeriodicBoundaryCondition: cannot load by " << param << " in two dimensional setup" << '\n';
                        exit(1);
                    }
                    if      ( param.compare("ex") == 0 ) {
                        strainFunc [ 0 ] = hnum;
                    } else if ( param.compare("ey") == 0 ) {
                        strainFunc [ 1 ] = hnum;
                    } else if ( param.compare("gxy") == 0 ) {
                        strainFunc [ 2 ] = hnum;
                    } else if ( param.compare("sx") == 0 ) {
                        stressFunc [ 0 ] = hnum;
                    } else if ( param.compare("sy") == 0 ) {
                        stressFunc [ 1 ] = hnum;
                    } else if ( param.compare("txy") == 0 ) {
                        stressFunc [ 2 ] = hnum;
                    } else {
                        cout << "Error in PeriodicBoundaryCondition: loaded by " << param << " not implemented yet" << '\n';
                        exit(1);
                    }
                } else if ( dim == 3 ) {
                    if      ( param.compare("ex") == 0 ) {
                        strainFunc [ 0 ] = hnum;
                    } else if ( param.compare("ey") == 0 ) {
                        strainFunc [ 1 ] = hnum;
                    } else if ( param.compare("ez") == 0 ) {
                        strainFunc [ 2 ] = hnum;
                    } else if ( param.compare("gyz") == 0 ) {
                        strainFunc [ 3 ] = hnum;
                    } else if ( param.compare("gxz") == 0 ) {
                        strainFunc [ 4 ] = hnum;
                    } else if ( param.compare("gxy") == 0 ) {
                        strainFunc [ 5 ] = hnum;
                    } else if ( param.compare("sx") == 0 ) {
                        stressFunc [ 0 ] = hnum;
                    } else if ( param.compare("sy") == 0 ) {
                        stressFunc [ 1 ] = hnum;
                    } else if ( param.compare("sz") == 0 ) {
                        stressFunc [ 2 ] = hnum;
                    } else if ( param.compare("tyz") == 0 ) {
                        stressFunc [ 3 ] = hnum;
                    } else if ( param.compare("txz") == 0 ) {
                        stressFunc [ 4 ] = hnum;
                    } else if ( param.compare("txy") == 0 ) {
                        stressFunc [ 5 ] = hnum;
                    } else {
                        cout << "Error in PeriodicBoundaryCondition: loaded by " << param << " not implemented yet" << '\n';
                        exit(1);
                    }
                }
            }
        } else if ( param.compare("size") == 0 ) {
            sizeB = true;
            iss >> num;
            PUCsize.resize(num);
            for ( unsigned i = 0; i < num; i++ ) {
                iss >> PUCsize [ i ];
            }
        }
    }

    if ( !sizeB ) {
        cout << "Error BasicPeriodicBC: size was not specified" << endl;
        exit(1);
    }
    if ( !pairsB ) {
        cout << "Error BasicPeriodicBC: pairs were not specified" << endl;
        exit(1);
    }
    if ( !loadB ) {
        cout << "Error BasicPeriodicBC: load was not specified" << endl;
        exit(1);
    }
}


void RigidPlate :: readFromLine(istringstream &iss, unsigned d) {
    // jointDoF jD;
    ndim = d;
    unsigned nslaves, nodeid;
    //////////////////////////////////////////////////////////
    // read the line "masterId numSlaves slaveId1, slaveId2...."
    iss >> master_id >> nslaves;

    for ( unsigned i = 0; i < nslaves; i++ ) {
        iss >> nodeid;
        slave_ids.push_back(nodeid);
    }

    bool bw = false;
    string param;
    while ( !iss.eof() ) {
        iss >> param;
        if ( param.compare("which") == 0 ) {
            iss >> which;
            // std::cout << "which" << '\n';
            bw = true;
            std :: cout << "using RigidPlate rigid in " << which << " direction, use proper BC for masterDoF (fix unused DoFs to zero)" << '\n';
        }
    }
    if ( !bw ) {
        which = "xyz";
    }
}

void RigidPlate :: apply(NodeContainer *nodes, ElementContainer *e, BCContainer *bcs, ConstraintContainer *constrs, FunctionContainer *funcs, ExporterContainer *ex) {
    ( void ) e;
    ( void ) bcs;
    ( void ) funcs;
    ( void ) ex;
    Node *master, *slave;
    //////////////////////////////////////////////////////////
    // read the line "masterId numSlaves slaveId1, slaveId2...."
    master = nodes->giveNode(master_id);
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

    for ( auto const &sl_id : slave_ids ) {
        slave = nodes->giveNode(sl_id);
        constrs->connectSlaveMaster(slave, master, ndim, which);
    }
}

void CoordRigidPlate :: readFromLine(istringstream &iss, unsigned d) {
    ndim = d;
    if ( d == 2 ) {
        double x0, x1, y0, y1;
        iss >> master_id >> x0 >> x1 >> y0 >> y1;
        leftBottom = Point(x0, y0, 0);
        rightTop = Point(x1, y1, 0);
    } else if ( d == 3 ) {
        double x0, x1, y0, y1, z0, z1;
        iss >> master_id >> x0 >> x1 >> y0 >> y1 >> z0 >> z1;
        leftBottom = Point(x0, y0, z0);
        rightTop = Point(x1, y1, z1);
    } else {
        std :: cerr << "dimension " << d << " not implemented yet" << '\n';
        exit(1);
    }
    bool bw = false;
    string param;
    while ( !iss.eof() ) {
        iss >> param;
        if ( param.compare("which") == 0 ) {
            iss >> which;
            // std::cout << "which" << '\n';
            bw = true;
            std :: cout << "using RigidPlate rigid in " << which << " direction, use proper BC for masterDoF (fix unused DoFs to zero)" << '\n';
        }
    }
    if ( !bw ) {
        which = "xyz";
    }
}

void CoordRigidPlate :: apply(NodeContainer *nodes, ElementContainer *e, BCContainer *bcs, ConstraintContainer *constrs, FunctionContainer *funcs, ExporterContainer *ex) {
    ( void ) e;
    ( void ) bcs;
    ( void ) funcs;
    ( void ) ex;
    // jointDoF jD;
    Node *master;

    master = nodes->giveNode(master_id);
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
                constrs->connectSlaveMaster(nod, master, ndim, which);
            }
        }
    }
}



//////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////
// CONTAINER FOR PBLOCKS
PBlockContainer :: ~PBlockContainer() {
    for ( vector< PBlock * > :: iterator f = blocks.begin(); f != blocks.end(); ++f ) {
        delete * f;
    }
}

//////////////////////////////////////////////////////////
void PBlockContainer :: setContainers(NodeContainer *n, ElementContainer *e, BCContainer *b, ConstraintContainer *c, FunctionContainer *f, ExporterContainer *ex) {
    nodes = n;
    elems = e;
    bcs = b;
    constrs = c;
    funcs = f;
    exporters = ex;
}

//////////////////////////////////////////////////////////
void PBlockContainer :: apply() {
    for ( auto b: blocks ) {
        b->apply(nodes, elems, bcs, constrs, funcs, exporters);
    }
}

//////////////////////////////////////////////////////////
void PBlockContainer :: readFromFile(const string filename, unsigned dim) {
    unsigned origsize = blocks.size();
    string line, ftype;
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
            iss >> ftype;
            if ( !ftype.rfind("#", 0) == 0 ) {
                if ( ftype.compare("BasicPeriodicBC") == 0 ) {
                    BasicPeriodicBC *newblock = new BasicPeriodicBC();
                    newblock->readFromLine(iss, dim);
                    blocks.push_back(newblock);
                } else if ( ftype.compare("RigidPlate") == 0 ) {
                    RigidPlate *newblock = new RigidPlate();
                    newblock->readFromLine(iss, dim);
                    blocks.push_back(newblock);
                } else if ( ftype.compare("CoordRigidPlate") == 0 ) {
                    CoordRigidPlate *newblock = new CoordRigidPlate();
                    newblock->readFromLine(iss, dim);
                    blocks.push_back(newblock);
                } else {
                    cerr << "Error: preprocessor block '" <<  ftype <<  "' is not implemented yet." << endl;
                    exit(EXIT_FAILURE);
                }
            }
        }
        inputfile.close();
        cout << "Input file '" <<  filename << "' succesfully loaded; " << blocks.size() - origsize << " preprocessing blocks found" << endl;
    } else {
        cerr << "Error: unable to open input file '" <<  filename <<  "'" << endl;
        exit(EXIT_FAILURE);
    }
}
