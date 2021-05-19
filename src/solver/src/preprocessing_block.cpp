#include "preprocessing_block.h"
#include "geometry.h"
// #include "misc.h"


std :: vector< double >PointToStdVector(const Point &p, unsigned dim = 3) {
    std :: vector< double >vect;
    vect.push_back(p.getX() );
    vect.push_back(p.getY() );
    if ( dim == 3 ) {
        vect.push_back(p.getZ() );
    }
    return vect;
}

bool containsChar(const std :: string &str, char c)
{
    // std::cout << "str = " << str << ", char = " << c << '\n';
    return str.find(c) != std :: string :: npos;
}


//////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////
// Mechanical Periodic BC
void MechanicalPeriodicBC :: generateNewDoFs(NodeContainer *nodes) {
    //create new degrees of freedom representing strains ex, ey, gammaxy=2exy or ex, ey, ez, gammyz, gammaxz, gammaxy,
    MechDoF *mn;
    initalNodeNum = nodes->giveSize();
    mn = new MechDoF(3 * ( dim - 1 ) );
    nodes->addNode(mn);
}

//////////////////////////////////////////////////////////
void MechanicalPeriodicBC :: generateConstraints(NodeContainer *nodes, ConstraintContainer *constrs) {
    //apply contraints, connect periodic images
    JointDoF *jd;
    vector< Node * >vm;
    vector< unsigned >dirs;
    vector< double >mults;
    Node *s = nullptr;
    Node *m = nullptr;
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
                dirs [ 2 ] = 5; //gamma xy
                dirs [ 3 ] = 4; //gamma xz
                vm [ 2 ] = nodes->giveNode(initalNodeNum);
                vm [ 3 ] = nodes->giveNode(initalNodeNum); //gamma xz
                mults [ 3 ] = diff.z / 2;
            } else if ( dim == 2 ) {
                vm.resize(3);
                mults.resize(3);
                dirs.resize(3, 0);
                dirs [ 2 ] = 2; //gamma xy
                vm [ 2 ] = nodes->giveNode(initalNodeNum);
            }
            dirs [ 0 ] = 0;
            dirs [ 1 ] = 0; //eps x
            vm [ 0 ] = m; //master
            vm [ 1 ] = nodes->giveNode(initalNodeNum);
            mults [ 0 ] = 1;
            mults [ 1 ] = diff.x;
            mults [ 2 ] = diff.y / 2;
            jd = new JointDoF(s, dirs [ 0 ], vm, dirs, mults);
            constrs->addConstraint(jd);

            //direction Y  (gammaxz realized here)
            if ( dim == 3 ) {
                dirs [ 2 ] = 5; //gamma xy
                dirs [ 3 ] = 3; //gamma yz
                mults [ 3 ] = diff.z / 2;
            }
            dirs [ 0 ] = 1;
            dirs [ 1 ] = 1; //eps y
            mults [ 1 ] = diff.y;
            mults [ 2 ] = diff.x / 2;
            jd = new JointDoF(s, dirs [ 0 ], vm, dirs, mults);
            constrs->addConstraint(jd);

            //direction Z  (gammaxz realized here)
            if ( dim == 3 ) {
                dirs [ 2 ] = 4; //gamma xz
                dirs [ 3 ] = 3; //gamma yz
                mults [ 3 ] = diff.y / 2;

                dirs [ 0 ] = 2;
                dirs [ 1 ] = 2; //eps z
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
                dirs [ 2 ] = 5;  //gamma xy
                dirs [ 3 ] = 4;  //gamma xz
                vm [ 2 ] = nodes->giveNode(initalNodeNum);
                vm [ 3 ] = nodes->giveNode(initalNodeNum);
                mults [ 3 ] = diff.z;
            } else if ( dim == 2 ) {
                vm.resize(3);
                mults.resize(3);
                dirs.resize(3, 0);
                dirs [ 2 ] = 2; //gamma xy
                vm [ 2 ] = nodes->giveNode(initalNodeNum);
            }
            dirs [ 0 ] = 0;
            dirs [ 1 ] = 0; //eps x
            vm [ 0 ] = m; //master
            vm [ 1 ] = nodes->giveNode(initalNodeNum);
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
                dirs [ 2 ] = 3; //gamma yz
                vm [ 2 ] = nodes->giveNode(initalNodeNum);
                mults [ 2 ] = diff.z;
            } else if ( dim == 2 ) {
                vm.resize(2);
                mults.resize(2);
                dirs.resize(2, 0);
            }
            dirs [ 0 ] = 1;
            dirs [ 1 ] = 1; //eps y
            vm [ 0 ] = m; //master
            vm [ 1 ] = nodes->giveNode(initalNodeNum);
            mults [ 0 ] = 1;
            mults [ 1 ] = diff.y;
            jd = new JointDoF(s, dirs [ 0 ], vm, dirs, mults);
            constrs->addConstraint(jd);

            //direction Z  (gammaxz realized here)
            if ( dim == 3 ) {
                vm.resize(2);
                mults.resize(2);
                dirs.resize(2, 0);
                dirs [ 0 ] = 2;
                dirs [ 1 ] = 2;  //eps z
                vm [ 0 ] = m; //master
                vm [ 1 ] = nodes->giveNode(initalNodeNum);
                mults [ 0 ] = 1;
                mults [ 1 ] = diff.z;
                jd = new JointDoF(s, dirs [ 0 ], vm, dirs, mults);
                constrs->addConstraint(jd);
            }
        }
    }
}

//////////////////////////////////////////////////////////
void MechanicalPeriodicBC :: generateExporters(NodeContainer *nodes, ExporterContainer *ex) {
    //export data
    string export_name = "PUCstrain_stress";
    vector< unsigned >n(1, initalNodeNum);
    vector< string >gname(3 * dim - 3);
    vector< string >codes(3 * dim - 3);
    ForceGauge *fg;
    gname [ 0 ] = "sigma_x";
    gname [ 1 ] = "sigma_y";
    gname [ 2 ] = "tau_xy";
    codes [ 0 ] = "0";
    codes [ 1 ] = "1";
    codes [ 2 ] = "2";
    if ( dim == 3 ) {
        gname [ 2 ] = "sigma_z";
        gname [ 3 ] = "tau_yz";
        gname [ 4 ] = "tau_xz";
        gname [ 5 ] = "tau_xy";
        codes [ 2 ] = "2";
        codes [ 3 ] = "3";
        codes [ 4 ] = "4";
        codes [ 5 ] = "5";
    }
    for ( unsigned i = 0; i < gname.size(); i++ ) {
        fg = new ForceGauge(export_name, gname [ i ], codes [ i ], n, nodes, 1. / volume, dim);
        ex->addExporter(fg);
    }

    DoFGauge *dg;
    gname [ 0 ] = "eps_x";
    gname [ 1 ] = "eps_y";
    gname [ 2 ] = "gamma_xy";
    codes [ 0 ] = "ux";
    codes [ 1 ] = "uy";
    codes [ 2 ] = "rz";
    if ( dim == 3 ) {
        gname [ 2 ] = "eps_z";
        gname [ 3 ] = "gamma_yz";
        gname [ 4 ] = "gamma_xz";
        gname [ 5 ] = "gamma_xy";
        codes [ 2 ] = "uz";
        codes [ 3 ] = "rx";
        codes [ 4 ] = "ry";
        codes [ 5 ] = "rz";
    }
    for ( unsigned i = 0; i < gname.size(); i++ ) {
        dg = new DoFGauge(export_name, gname [ i ], codes [ i ], n, nodes, 1., dim);
        ex->addExporter(dg);
    }
}

//////////////////////////////////////////////////////////
void MechanicalPeriodicBC :: generateRigidBodyBC(NodeContainer *nodes, ElementContainer *elems, BCContainer *bcs, ConstraintContainer *constrs, FunctionContainer *funcs) {
    if ( volumetricAverageRigidBC < 0 ) {  //last master node cannot move
        Node *m = constrs->giveConstraint(constrs->giveSize() - 1)->giveMasterNode(0);// warning C4267: 'argument': conversion from 'size_t' to 'const unsigned int', possible loss of data
        BoundaryCondition *bc;
        vector< int >dBC, nBC;
        dBC.resize(m->giveNumberOfDoFs(), funcs->giveSize() ); //todo: warning C4267: 'argument': conversion from 'size_t' to 'const _Ty', possible loss of data
        nBC.resize(m->giveNumberOfDoFs(), -1);
        bc = new BoundaryCondition(m, dBC, nBC);
        bcs->addBoundaryCondition(bc);

        //add constant function
        vector< double >x, y;
        x.resize(1, 0);
        y.resize(1, 0);
        PieceWiseLinearFunction *newf = new PieceWiseLinearFunction(x, y);
        funcs->addFunction(newf);
    } else {  //volumetric average
        VolumetricAverage *va;
        vector< Node * >vm;

        for ( unsigned n = 0; n < nodes->giveSize(); n++ ) {
            if ( nodes->giveNode(n)->doesMechanics() && ( dynamic_cast< MechDoF * >( nodes->giveNode(n) ) == nullptr ) ) {
                vm.push_back(nodes->giveNode(n) );
            }
        }
        if ( vm.size() > 0 ) {
            unsigned nDoFs = 3;
            if ( dim == 3 ) {
                nDoFs = 6;
            }
            MechDoF *pn = new MechDoF(nDoFs);
            nodes->addNode(pn);

            vector< unsigned >dirs(vm.size() );

            for ( unsigned vi = 0; vi < nDoFs; vi++ ) {
                fill(dirs.begin(), dirs.end(), vi);
                va = new VolumetricAverage(vm, dirs, pn, vi, elems, constrs);
                constrs->addConstraint(va);
            }

            BoundaryCondition *bc;
            vector< int >dBC, nBC;
            dBC.resize(nDoFs, volumetricAverageRigidBC);
            nBC.resize(nDoFs, -1);
            bc = new BoundaryCondition(pn, dBC, nBC);
            bcs->addBoundaryCondition(bc);
        }
    }
}

//////////////////////////////////////////////////////////
void MechanicalPeriodicBC :: apply(NodeContainer *nodes, ElementContainer *e, BCContainer *bcs, ConstraintContainer *constrs, FunctionContainer *funcs, ExporterContainer *ex) {
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
    generateNewDoFs(nodes);

    //apply contraints, connect periodic images
    generateConstraints(nodes, constrs);

    //boundary conditions
    generateRigidBodyBC(nodes, e, bcs, constrs, funcs);

    //set prescribed strain and stress
    vector< double >bcmults;
    BoundaryCondition *bc;
    vector< int >dBC, nBC;

    unsigned n = strainFunc.size();
    dBC.resize(n, -1);
    nBC.resize(n, -1);
    bcmults.resize(n, 1);
    for ( unsigned i = 0; i < strainFunc.size(); i++ ) {
        if ( strainFunc [ i ] >= 0 ) {
            dBC [ i ] = strainFunc [ i ];
        }
        if ( stressFunc [ i ] >= 0 ) {
            bcmults [ i ] = volume;
            nBC [ i ] = stressFunc [ i ];
        }
        if ( strainFunc [ i ] >= 0 && stressFunc [ i ] >= 0 ) {
            cerr << "Error in Periodic boundary condition: cannot prescribe both stress and strain for the same direction" << endl;
        }
    }
    bc = new BoundaryCondition(nodes->giveNode(initalNodeNum), dBC, nBC, bcmults);
    bcs->addBoundaryCondition(bc);

    //export data
    generateExporters(nodes, ex);

    cout << "Applied periodic boundary conditions: " << nodes->giveSize() - initalNodeNum << " new DoFs (nodes " << initalNodeNum << " - " <<  nodes->giveSize() - 1 << "); " << constrs->giveSize() - const_num << " new constraints; " << bcs->giveSize() - bcs_num << " new boundary conditions; " << funcs->giveSize() - funcs_num << " new function; " << ex->giveSize() - ex_num << " new exporters; " << "created" << endl;
}

//////////////////////////////////////////////////////////
void MechanicalPeriodicBC :: readLoading(istringstream &iss) {
    unsigned num, hnum;
    string param;
    volumetricAverageRigidBC = -1;
    iss >> num;

    strainFunc.resize(3 * ( dim - 1 ), -1);
    stressFunc.resize(3 * ( dim - 1 ), -1);

    for ( unsigned i = 0; i < num; i++ ) {
        iss >> param >> hnum;
        if ( dim == 2 ) {
            std :: size_t found = param.find("z");
            if ( found != std :: string :: npos ) {
                cout << "Error in " << name << " : cannot load by " << param << " in two dimensional setup" << '\n';
                exit(1);
            }
            if ( param.compare("ex") == 0 ) {
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
            } else if ( param.compare("volumetricAverage") == 0 ) {
                volumetricAverageRigidBC = hnum;
            } else {
                cout << "Error in " << name << " : loading by " << param << " not implemented yet" << '\n';
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
            } else if ( param.compare("volumetricAverage") == 0 ) {
                volumetricAverageRigidBC = hnum;
            } else {
                cout << "Error in " << name << " : loading by " << param << " not implemented yet" << '\n';
                exit(1);
            }
        }
    }
}

//////////////////////////////////////////////////////////
void MechanicalPeriodicBC :: readFromLine(istringstream &iss, unsigned d) {
    use_half_gammas = false;
    dim = d;
    string param;
    unsigned num;

    bool sizeB, loadB, pairsB;
    sizeB = loadB = pairsB = false;

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
            readLoading(iss);
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
        cout << "Error " << name << " : size was not specified" << endl;
        exit(1);
    }
    if ( !pairsB ) {
        cout << "Error " << name << " : pairs were not specified" << endl;
        exit(1);
    }
    if ( !loadB ) {
        cout << "Error " << name << " : load was not specified" << endl;
        exit(1);
    }
}


//////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////
// Transport Periodic BC
//////////////////////////////////////////////////////////
void TransportPeriodicBC :: generateNewDoFs(NodeContainer *nodes) {
    //create new degrees of freedom representing strains ex, ey, gammaxy=2exy or ex, ey, ez, gammyz, gammaxz, gammaxy,
    TrsDoF *mn;
    initalNodeNum = nodes->giveSize(); //todo warning C4267: '=': conversion from 'size_t' to 'unsigned int', possible loss of data
    mn = new TrsDoF(dim);
    nodes->addNode(mn);
}

//////////////////////////////////////////////////////////
void TransportPeriodicBC :: generateConstraints(NodeContainer *nodes, ConstraintContainer *constrs) {
    //apply contraints, connect periodic images
    JointDoF *jd;
    vector< Node * >vm;
    vector< unsigned >dirs;
    vector< double >mults;
    Node *s = nullptr;
    Node *m = nullptr;
    Point diff;
    for ( unsigned i = 0; i < masters.size(); i++ ) {
        m = nodes->giveNode(masters [ i ]);
        s = nodes->giveNode(slaves [ i ]);

        //connect pressure
        diff = s->givePoint() - m->givePoint();

        if ( dim == 3 ) {
            vm.resize(4);
            mults.resize(4);
            dirs.resize(4, 0);
            dirs [ 3 ] = 2;
            mults [ 3 ] = diff.z;
            vm [ 3 ] = nodes->giveNode(initalNodeNum);
        } else if ( dim == 2 ) {
            vm.resize(3);
            mults.resize(3);
            dirs.resize(3, 0);
        }
        dirs [ 0 ] = 0;
        dirs [ 1 ] = 0;
        dirs [ 2 ] = 1;
        vm [ 0 ] = m; //master
        vm [ 1 ] = nodes->giveNode(initalNodeNum);
        vm [ 2 ] = nodes->giveNode(initalNodeNum);
        mults [ 0 ] = 1;
        mults [ 1 ] = diff.x;
        mults [ 2 ] = diff.y;
        jd = new JointDoF(s, dirs [ 0 ], vm, dirs, mults);
        constrs->addConstraint(jd);
    }
}

//////////////////////////////////////////////////////////
void TransportPeriodicBC :: generateExporters(NodeContainer *nodes, ExporterContainer *ex) {
    //export data
    string export_name = "PUCgrad_flux";
    vector< string >gname;
    vector< unsigned >n(1, initalNodeNum);
    vector< string >codes;
    codes.resize(dim);
    gname.resize(dim);

    ForceGauge *fg;
    gname [ 0 ] = "flux_x";
    gname [ 1 ] = "flux_y";
    codes [ 0 ] = "0";
    codes [ 1 ] = "1";
    if ( dim == 3 ) {
        gname [ 2 ] = "flux_z";
        codes [ 2 ] = "2";
    }
    for ( unsigned i = 0; i < gname.size(); i++ ) {
        fg = new ForceGauge(export_name, gname [ i ], codes [ i ], n, nodes, 1. / volume, dim);
        ex->addExporter(fg);
    }

    DoFGauge *dg;
    gname [ 0 ] = "grad_x";
    gname [ 1 ] = "grad_y";
    codes [ 0 ] = "ux";
    codes [ 1 ] = "uy";
    if ( dim == 3 ) {
        gname [ 2 ] = "grad_z";
        codes [ 2 ] = "uz";
    }
    for ( unsigned i = 0; i < gname.size(); i++ ) {
        dg = new DoFGauge(export_name, gname [ i ], codes [ i ], n, nodes, 1., dim);
        ex->addExporter(dg);
    }
}

//////////////////////////////////////////////////////////
void TransportPeriodicBC :: readLoading(istringstream &iss) {
    unsigned num, hnum;
    string param;
    strainFunc.resize(dim, -1);
    stressFunc.resize(dim, -1);
    volumetricAverageRigidBC = -1;
    microscaleSources.resize(dim, -1);

    iss >> num;
    for ( unsigned i = 0; i < num; i++ ) {
        iss >> param >> hnum;
        if ( dim == 2 ) {
            std :: size_t found = param.find("z");
            if ( found != std :: string :: npos ) {
                cout << "Error in " << name << " : cannot load by " << param << " in two dimensional setup" << '\n';
                exit(1);
            }
            if ( param.compare("gx") == 0 ) {
                strainFunc [ 0 ] = hnum;
            } else if ( param.compare("gy") == 0 ) {
                strainFunc [ 1 ] = hnum;
            } else if ( param.compare("jx") == 0 ) {
                stressFunc [ 0 ] = hnum;
            } else if ( param.compare("jy") == 0 ) {
                stressFunc [ 1 ] = hnum;
            } else if ( param.compare("volumetricAverage") == 0 ) {
                volumetricAverageRigidBC = hnum;
            } else if ( param.compare("microSources") == 0 ) {
                microscaleSources [ 0 ] = hnum;
                iss >> microscaleSources [ 1 ];
            } else {
                cout << "Error in " << name << " : loading by " << param << " not implemented yet" << '\n';
                exit(1);
            }
        } else if ( dim == 3 ) {
            if      ( param.compare("gx") == 0 ) {
                strainFunc [ 0 ] = hnum;
            } else if ( param.compare("gy") == 0 ) {
                strainFunc [ 1 ] = hnum;
            } else if ( param.compare("gz") == 0 ) {
                strainFunc [ 2 ] = hnum;
            } else if ( param.compare("jx") == 0 ) {
                stressFunc [ 0 ] = hnum;
            } else if ( param.compare("jy") == 0 ) {
                stressFunc [ 1 ] = hnum;
            } else if ( param.compare("jz") == 0 ) {
                stressFunc [ 2 ] = hnum;
            } else if ( param.compare("volumetricAverage") == 0 ) {
                volumetricAverageRigidBC = hnum;
            } else if ( param.compare("microSources") == 0 ) {
                microscaleSources [ 0 ] = hnum;
                iss >> microscaleSources [ 1 ];
                iss >> microscaleSources [ 2 ];
            } else {
                cout << "Error in " << name << " : loading by " << param << " not implemented yet" << '\n';
                exit(1);
            }
        }
    }
}

//////////////////////////////////////////////////////////
void TransportPeriodicBC :: generateRigidBodyBC(NodeContainer *nodes, ElementContainer *elems, BCContainer *bcs, ConstraintContainer *constrs, FunctionContainer *funcs) {
    if ( volumetricAverageRigidBC < 0 ) {  //last master node cannot move
        Node *m = constrs->giveConstraint(constrs->giveSize() - 1)->giveMasterNode(0); //todo:  warning C4267: 'argument': conversion from 'size_t' to 'const unsigned int', possible loss of data
        BoundaryCondition *bc;
        vector< int >dBC, nBC;
        dBC.resize(m->giveNumberOfDoFs(), funcs->giveSize() ); //todo: conversion from 'size_t' to 'const _Ty', possible loss of data
        nBC.resize(m->giveNumberOfDoFs(), -1);
        bc = new BoundaryCondition(m, dBC, nBC);
        bcs->addBoundaryCondition(bc);

        //add constant function
        vector< double >x, y;
        x.resize(1, 0);
        y.resize(1, 0);
        PieceWiseLinearFunction *newf = new PieceWiseLinearFunction(x, y);
        funcs->addFunction(newf);
    } else {  //volumetric average
        VolumetricAverage *va;
        vector< Node * >vm;
        for ( unsigned n = 0; n < nodes->giveSize(); n++ ) {
            if ( nodes->giveNode(n)->doesTransport() && ( dynamic_cast< TrsDoF * >( nodes->giveNode(n) ) == nullptr ) ) {
                vm.push_back(nodes->giveNode(n) );
            }
        }
        if ( vm.size() > 0 ) {
            TrsDoF *tn = new TrsDoF(1);
            nodes->addNode(tn);

            vector< unsigned >dirs(vm.size() );
            va = new VolumetricAverage(vm, dirs, tn, 0, elems, constrs);
            constrs->addConstraint(va);

            BoundaryCondition *bc;
            vector< int >dBC, nBC;
            dBC.resize(1, volumetricAverageRigidBC);
            nBC.resize(1, -1);
            bc = new BoundaryCondition(tn, dBC, nBC);
            bcs->addBoundaryCondition(bc);
        }
    }
}

//////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////

void RigidPlate :: setDirectionToFix(istringstream &iss) {
    bool bw = false;
    string param;
    while ( !iss.eof() ) {
        iss >> param;
        if ( param.compare("which") == 0 ) {
            iss >> which;
            // std::cout << "which" << '\n';
            bw = true;
            std :: cout << "using RigidPlate rigid in " << which << " direction, use proper BC for Particle (fix unused DoFs to zero)" << '\n';
        }
    }
    if ( !bw ) {
        this->which = "xyz"; //abc for rotations
    }
}


void RigidPlate :: readFromLine(istringstream &iss, unsigned d) {
    // jointDoF jD;
    this->dim = d;
    unsigned nslaves, nodeid;
    //////////////////////////////////////////////////////////
    // read the line "masterId numSlaves slaveId1, slaveId2...."
    iss >> master_id >> nslaves;

    for ( unsigned i = 0; i < nslaves; i++ ) {
        iss >> nodeid;
        slave_ids.push_back(nodeid);
    }

    this->setDirectionToFix(iss);
}

void RigidPlate :: checkMechTransport(Node *master) {
    // in case of rigid plate, master is a virtual virtual node and not a physical particle or
    master->setName(master->giveName().append("-virtual") );

    if ( dynamic_cast< MechNode * >( master ) ) {
        if ( master->giveNumberOfDoFs() != ( 3 * ( this->dim - 1 ) ) ) {
            cerr << "Error in " << __func__ << ": Master for RigidPlate in mechnics must have " << ( 3 * ( this->dim - 1 ) ) << " DoFs, " << master->giveNumberOfDoFs() << " provided" << '\n';
            exit(EXIT_FAILURE);
        }
    } else if ( dynamic_cast< TrsNode * >( master ) ) {
        this->transport = true;
        if ( master->giveNumberOfDoFs() != 1 ) {
            cerr << "Error in " << __func__ << ": Master for RigidPlate in transport must have " << 1 << " DoFs, " << master->giveNumberOfDoFs() << " provided" << '\n';
            exit(EXIT_FAILURE);
        }
    } else {
        cerr << "Error in " << __func__ << ": Master for RigidPlate is niether mechanical nor transport " << '\n';
        exit(EXIT_FAILURE);
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
    this->checkMechTransport(master);

    for ( auto const &sl_id : slave_ids ) {
        slave = nodes->giveNode(sl_id);
        connectSlaveMasterRigid(constrs, slave, master, this->dim, which, this->transport);
    }
}

void CoordRigidPlate :: readFromLine(istringstream &iss, unsigned d) {
    this->dim = d;
    if ( d == 2 ) {
        double x0, x1, y0, y1;
        iss >> master_id >> x0 >> x1 >> y0 >> y1;
        leftBottom = Point(x0, y0, -100);
        rightTop = Point(x1, y1, 100);
    } else if ( d == 3 ) {
        double x0, x1, y0, y1, z0, z1;
        iss >> master_id >> x0 >> x1 >> y0 >> y1 >> z0 >> z1;
        leftBottom = Point(x0, y0, z0);
        rightTop = Point(x1, y1, z1);
    } else {
        std :: cerr << "dimension " << d << " not implemented yet" << '\n';
        exit(EXIT_FAILURE);
    }
    RigidPlate :: setDirectionToFix(iss);
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
    RigidPlate :: checkMechTransport(master);

    for ( auto const &nod : * nodes ) {
        if ( isInBlock(nod->givePoint(), leftBottom, rightTop) ) {
            if ( nod == master ) {
                continue;
            }
            connectSlaveMasterRigid(constrs, nod, master, this->dim, which, this->transport);
        }
    }
}


void RingRigidPlate :: readFromLine(istringstream &iss, unsigned d) {
    this->dim = d;
    double x, y, z, rI, rO;
    string dir;
    iss >> master_id >> x >> y;
    if ( d == 3 ) {
        iss >> z;
    } else {
        z = 0.0;
    }
    this->center = Point(x, y, z);
    iss >> rI >> rO;
    this->r_inner = rI;
    this->r_outer = rO;
    if ( dim == 3 ) {
        iss >> dir;
        if ( dir.compare("x") ) {
            direction = 0;
        } else if ( dir.compare("y") ) {
            direction = 1;
        } else if ( dir.compare("z") ) {
            direction = 2;
        }
        // iss >> x >> y >> z;
        // axis = Point(x, y, z);
    } else {
        // for 2D case, normal
        // axis = Point(0, 0, 1);
        direction = 2;
    }

    RigidPlate :: setDirectionToFix(iss);
}

void RingRigidPlate :: apply(NodeContainer *nodes, ElementContainer *e, BCContainer *bcs, ConstraintContainer *constrs, FunctionContainer *funcs, ExporterContainer *ex) {
    ( void ) e;
    ( void ) bcs;
    ( void ) funcs;
    ( void ) ex;
    // jointDoF jD;
    Node *master;

    master = nodes->giveNode(master_id);
    // check if it is master node
    RigidPlate :: checkMechTransport(master);

    Point node_point;
    int xm, ym, zm;
    xm = 1;
    ym = 1;
    zm = 1;
    if ( direction == 0 ) {
        xm = 0;
    } else if ( direction == 1 ) {
        ym = 0;
    } else if ( direction == 2 ) {
        zm = 0;
    }
    this->center = Point(this->center.getX() * xm, this->center.getY() * ym, this->center.getZ() * zm);
    for ( auto const &nod : * nodes ) {
        node_point = Point(nod->givePoint().getX() * xm, nod->givePoint().getY() * ym, nod->givePoint().getZ() * zm);
        if ( isInCircle(node_point, this->center, this->r_outer) ) {
            if ( !isInCircle(node_point, this->center, this->r_inner) ) {
                if ( nod == master ) {
                    continue;
                }
                connectSlaveMasterRigid(constrs, nod, master, this->dim, which, this->transport);
            }
        }
    }
}



void ExpansionRing :: readFromLine(istringstream &iss, unsigned d) {
    RingRigidPlate :: readFromLine(iss, d);
    if ( typeid(this) != typeid(ExpansionRing) ) {
        // do not perform following for derived classes
        return;
    }
    iss.clear(); // clear string stream
    iss.seekg(0, iss.beg); //reset position in string stream
    bool bf = false;
    string param;
    while ( !iss.eof() ) {
        iss >> param;
        if ( param.compare("volExpFn") == 0 ) {
            iss >> this->fn_id;
            bf = true;
        }
    }
    if ( !bf ) {
        std :: cerr << "Error: no function governing volumetric expansion specified" << '\n';
        exit(EXIT_FAILURE);
    }
}

void ExpansionRing :: apply(NodeContainer *nodes, ElementContainer *e, BCContainer *bcs, ConstraintContainer *constrs, FunctionContainer *funcs, ExporterContainer *ex) {
    ( void ) e;
    ( void ) bcs;
    ( void ) funcs;
    ( void ) ex;
    // jointDoF jD;
    Node *master;

    master = nodes->giveNode(master_id);
    // check if it is master node
    RigidPlate :: checkMechTransport(master);

    Point node_point;
    int xm, ym, zm;
    xm = 1;
    ym = 1;
    zm = 1;
    if ( direction == 0 ) {
        xm = 0;
    } else if ( direction == 1 ) {
        ym = 0;
    } else if ( direction == 2 ) {
        zm = 0;
    }
    this->center = Point(this->center.getX() * xm, this->center.getY() * ym, this->center.getZ() * zm);
    for ( auto const &nod : * nodes ) {
        node_point = Point(nod->givePoint().getX() * xm, nod->givePoint().getY() * ym, nod->givePoint().getZ() * zm);
        if ( isInCircle(node_point, this->center, this->r_outer) ) {
            if ( !isInCircle(node_point, this->center, this->r_inner) ) {
                if ( nod == master ) {
                    continue;
                }
                connectSlaveMasterExpansion(constrs, nod, master, this->dim, this->transport, funcs->giveFunction(this->fn_id) );
            }
        }
    }
}


void ExpansionRingDoFLoad :: readFromLine(istringstream &iss, unsigned d) {
    RingRigidPlate :: readFromLine(iss, d);
    if ( typeid(this) != typeid(ExpansionRing) ) {
        // do not perform following for derived classes
        return;
    }
    iss.clear(); // clear string stream
    iss.seekg(0, iss.beg); //reset position in string stream
    bool bf = false;
    string param;
    while ( !iss.eof() ) {
        iss >> param;
        if ( param.compare("expansionMaster") == 0 ) {
            iss >> this->expansion_master_id;
            bf = true;
        }
    }
    if ( !bf ) {
        std :: cerr << "Error: no master DoF governing volumetric expansion specified" << '\n';
        exit(EXIT_FAILURE);
    }
}

void ExpansionRingDoFLoad :: apply(NodeContainer *nodes, ElementContainer *e, BCContainer *bcs, ConstraintContainer *constrs, FunctionContainer *funcs, ExporterContainer *ex) {
    ( void ) e;
    ( void ) bcs;
    ( void ) funcs;
    ( void ) ex;
    // jointDoF jD;
    Node *master;
    Node *expMaster;

    master = nodes->giveNode(this->master_id);
    RigidPlate :: checkMechTransport(master);

    expMaster = nodes->giveNode(this->expansion_master_id);


    Point node_point;
    int xm, ym, zm;
    xm = 1;
    ym = 1;
    zm = 1;
    if ( direction == 0 ) {
        xm = 0;
    } else if ( direction == 1 ) {
        ym = 0;
    } else if ( direction == 2 ) {
        zm = 0;
    }
    this->center = Point(this->center.getX() * xm, this->center.getY() * ym, this->center.getZ() * zm);
    for ( auto const &nod : * nodes ) {
        node_point = Point(nod->givePoint().getX() * xm, nod->givePoint().getY() * ym, nod->givePoint().getZ() * zm);
        if ( isInCircle(node_point, this->center, this->r_outer) ) {
            if ( !isInCircle(node_point, this->center, this->r_inner) ) {
                if ( nod == master ) {
                    continue;
                }
                connectSlaveMasterExpansionFLoad(constrs, nod, master, expMaster, this->dim);
            }
        }
    }
}


void ExpansionRingSingleDoFLoad :: apply(NodeContainer *nodes, ElementContainer *e, BCContainer *bcs, ConstraintContainer *constrs, FunctionContainer *funcs, ExporterContainer *ex) {
    ( void ) e;
    ( void ) bcs;
    ( void ) funcs;
    ( void ) ex;
    // jointDoF jD;
    Node *slave;
    Node *expMaster;

    expMaster = nodes->giveNode(this->master_id);

    Point node_point;
    int xm, ym, zm;
    xm = 1;
    ym = 1;
    zm = 1;
    if ( direction == 0 ) {
        xm = 0;
    } else if ( direction == 1 ) {
        ym = 0;
    } else if ( direction == 2 ) {
        zm = 0;
    }

    // // NOTE JK this is intended for preferred choice of slave node (dof) according to max/min position in dir of user choice
    // unsigned x_max_id, y_max_id, z_max_id, x_min_id, y_min_id, z_min_id;
    // double x_max, y_max, z_max, x_min, y_min, z_min;

    unsigned slave_dir;
    unsigned num_nodes = 0;
    bool slave_used = false;

    double l_i;
    Point n_i;

    vector< Node * >masterNodes;
    vector< double >multipliers;
    vector< unsigned >directions;

    std :: vector< double >n_vect;
    double slave_dir_vect_value;

    // n_vect_slave uložit dopředu
    // dělit všechno n_vect_slave dir
    //
    // zjistit na lineárním výpočtu to jak má být contraint

    this->center = Point(this->center.getX() * xm, this->center.getY() * ym, this->center.getZ() * zm);
    for ( auto const &nod : * nodes ) {
        if ( nod->giveName().compare("Particle") != 0){
            continue;
        }
        node_point = Point(nod->givePoint().getX() * xm, nod->givePoint().getY() * ym, nod->givePoint().getZ() * zm);
        if ( isInCircle(node_point, this->center, this->r_outer) ) {
            if ( !isInCircle(node_point, this->center, this->r_inner) ) {
                // if ( nod == expMaster ) {
                //     continue;
                // }
                num_nodes++;
                // std::cout << "nod name = " << nod->giveName() << '\n';
                n_i = ( nod->givePoint() - this->center );
                l_i = n_i.norm();
                n_i.normalize();
                // std::cout << "dim = " << this->dim << '\n';
                n_vect = PointToStdVector(n_i, this->dim);
                if ( slave_used ) {
                    for ( unsigned i = 0; i < this->dim; i++ ){
                        // std::cout << "i = " << i << ", num_nodes = " << num_nodes << '\n';
                        masterNodes.push_back(nod);
                        directions.push_back(i);
                        multipliers.push_back(- n_vect[ i ]);
                    }
                } else {
                    // first node is taken as a slave
                    slave = nod;
                    if ( * std :: max_element(n_vect.begin(), n_vect.end()) >  abs( * std :: min_element(n_vect.begin(), n_vect.end() )) ) {
                        slave_dir = std :: distance(n_vect.begin(), std :: max_element(n_vect.begin(), n_vect.end() ) );
                    } else {
                        slave_dir = std :: distance(n_vect.begin(), std :: min_element(n_vect.begin(), n_vect.end() ) );
                    }
                    slave_dir_vect_value = n_vect[slave_dir];
                    for ( unsigned i = 0; i < this->dim; i++ ) {
                        if ( i != slave_dir ) {
                            // multiplying self other DoFs
                            multipliers.push_back(-n_vect [ i ]);
                            masterNodes.push_back(slave);
                            directions.push_back(i);
                        }
                    }
                    slave_used = true;
                }
            }
        }
    }

    // adding master DoF governing the expansion
    multipliers.push_back(double(num_nodes) * this->r_outer);
    masterNodes.push_back(expMaster);
    directions.push_back(0);

    if ( !masterNodes.empty() ) {
        for ( unsigned j = 0; j < multipliers.size(); j++ ) {
            multipliers[ j ] /= slave_dir_vect_value;
        }
        JointDoF *newJD = new JointDoF(slave, slave_dir, masterNodes, directions, multipliers);
        constrs->addConstraint(newJD);
    }

    masterNodes.clear();
    multipliers.clear();
    directions.clear();
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
    unsigned origsize = blocks.size(); //todo: warning C4267: 'initializing': conversion from 'size_t' to 'unsigned int', possible loss of dat
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
                if ( ftype.compare("MechanicalPeriodicBC") == 0 ) {
                    MechanicalPeriodicBC *newblock = new MechanicalPeriodicBC();
                    newblock->readFromLine(iss, dim);
                    blocks.push_back(newblock);
                } else if ( ftype.compare("TransportPeriodicBC") == 0 ) {
                    MechanicalPeriodicBC *newblock = new TransportPeriodicBC();
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
                } else if ( ftype.compare("RingRigidPlate") == 0 ) {
                    RingRigidPlate *newblock = new RingRigidPlate();
                    newblock->readFromLine(iss, dim);
                    blocks.push_back(newblock);
                } else if ( ftype.compare("ExpansionRing") == 0 ) {
                    ExpansionRing *newblock = new ExpansionRing();
                    newblock->readFromLine(iss, dim);
                    blocks.push_back(newblock);
                } else if ( ftype.compare("ExpansionRingDoFLoad") == 0 ) {
                    ExpansionRingDoFLoad *newblock = new ExpansionRingDoFLoad();
                    newblock->readFromLine(iss, dim);
                    blocks.push_back(newblock);
                } else if ( ftype.compare("ExpansionRingSingleDoFLoad") == 0 ) {
                    ExpansionRingSingleDoFLoad *newblock = new ExpansionRingSingleDoFLoad();
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



void connectSlaveMasterRigid(ConstraintContainer *constrs, Node *slave, Node *master, unsigned const &ndim, const string &which, const bool trsp) {
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
            if ( containsChar(which,'x') && ( i == 0 || i == 4 || i == nDoFsPerNode - 1 ) ) {
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
        constrs->addConstraint(newJD);
        masterNodes.clear();
        multipliers.clear();
        directions.clear();
    }
}



void connectSlaveMasterExpansionFLoad(ConstraintContainer *constrs, Node *slave, Node *master, Node *expMaster, unsigned const &ndim) {
    if ( !dynamic_cast< Particle * >( slave ) ) {
        return;                                         // NOTE could be MechNode, but so far, nDoFs corresponds to Particles
    }


    if ( slave->giveNumberOfDoFs() != master->giveNumberOfDoFs() ) {
        std :: cerr << "slave and master must have the same number of DoFs, slave numDoFs = " << slave->giveNumberOfDoFs() << ", master numDoFs = " << master->giveNumberOfDoFs() << '\n';
        exit(1);
    }

    vector< Node * >masterNodes;
    vector< double >multipliers;
    vector< double >time_multipliers;
    vector< Function * >time_fns;
    vector< unsigned >directions;


    Point n = ( slave->givePoint() - master->givePoint() );
    double l = n.norm();
    n.normalize();


    std :: vector< double >n_vect = PointToStdVector(n, ndim);

    unsigned slave_dir;
    if ( * std :: max_element(n_vect.begin(), n_vect.end()) >  abs( * std :: min_element(n_vect.begin(), n_vect.end() )) ) {
        slave_dir = std :: distance(n_vect.begin(), std :: max_element(n_vect.begin(), n_vect.end() ) );
    } else {
        slave_dir = std :: distance(n_vect.begin(), std :: min_element(n_vect.begin(), n_vect.end() ) );
    }

    for ( unsigned i = 0; i < ndim; i++ ) {
        if ( i == slave_dir ) {
            multipliers.push_back(1.0);
            masterNodes.push_back(master);
            directions.push_back(i);
            time_multipliers.push_back(0.0);
            time_fns.push_back(nullptr);
        } else {
            multipliers.push_back(n_vect [ i ] / n_vect [ slave_dir ]);
            masterNodes.push_back(master);
            directions.push_back(i);
            time_multipliers.push_back(0.0);
            time_fns.push_back(nullptr);

            // multiplying self other DoFs
            multipliers.push_back(-n_vect [ i ] / n_vect [ slave_dir ]);
            masterNodes.push_back(slave);
            directions.push_back(i);
            time_multipliers.push_back(0.0);
            time_fns.push_back(nullptr);
        }
    }

    // adding master DoF governing the expansion
    multipliers.push_back(l / n_vect [ slave_dir ]);
    masterNodes.push_back(expMaster);
    directions.push_back(0);
    time_multipliers.push_back(0.0);
    time_fns.push_back(nullptr);

    if ( !masterNodes.empty() ) {
        JointDoF *newJD = new JointDoF(slave, slave_dir, masterNodes, directions, multipliers, time_fns, time_multipliers);
        constrs->addConstraint(newJD);
    }

    masterNodes.clear();
    multipliers.clear();
    directions.clear();
    time_multipliers.clear();
    time_fns.clear();
}




void connectSlaveMasterExpansion(ConstraintContainer *constrs, Node *slave, Node *master, unsigned const &ndim, const bool trsp, Function *fn) {
    ( void ) trsp;

    if ( !dynamic_cast< Particle * >( slave ) ) {
        return;                                         // NOTE could be MechNode, but so far, nDoFs corresponds to Particles
    }


    if ( slave->giveNumberOfDoFs() != master->giveNumberOfDoFs() ) {
        std :: cerr << "slave and master must have the same number of DoFs, slave numDoFs = " << slave->giveNumberOfDoFs() << ", master numDoFs = " << master->giveNumberOfDoFs() << '\n';
        exit(1);
    }

    vector< Node * >masterNodes;
    vector< double >multipliers;
    vector< double >time_multipliers;
    vector< Function * >time_fns;
    vector< unsigned >directions;


    Point n = ( slave->givePoint() - master->givePoint() );
    double l = n.norm();
    n.normalize();


    std :: vector< double >n_vect = PointToStdVector(n, ndim);

    unsigned slave_dir;
    if ( * std :: max_element(n_vect.begin(), n_vect.end()) >  abs( * std :: min_element(n_vect.begin(), n_vect.end() )) ) {
        slave_dir = std :: distance(n_vect.begin(), std :: max_element(n_vect.begin(), n_vect.end() ) );
    } else {
        slave_dir = std :: distance(n_vect.begin(), std :: min_element(n_vect.begin(), n_vect.end() ) );
    }

    for ( unsigned i = 0; i < ndim; i++ ) {
        if ( i == slave_dir ) {
            multipliers.push_back(1.0);
            masterNodes.push_back(master);
            directions.push_back(i);
            time_multipliers.push_back(l / n_vect [ i ]);
            time_fns.push_back(fn);
        } else {
            multipliers.push_back(n_vect [ i ] / n_vect [ slave_dir ]);
            masterNodes.push_back(master);
            directions.push_back(i);
            time_multipliers.push_back(0.0);
            time_fns.push_back(nullptr);

            // multiplying self other DoFs
            multipliers.push_back(-n_vect [ i ] / n_vect [ slave_dir ]);
            masterNodes.push_back(slave);
            directions.push_back(i);
            time_multipliers.push_back(0.0);
            time_fns.push_back(nullptr);
        }
    }

    if ( !masterNodes.empty() ) {
        JointDoF *newJD = new JointDoF(slave, slave_dir, masterNodes, directions, multipliers, time_fns, time_multipliers);
        constrs->addConstraint(newJD);
    }

    masterNodes.clear();
    multipliers.clear();
    directions.clear();
    time_multipliers.clear();
    time_fns.clear();
}

// void connectSlaveMasterExpansionOLD(ConstraintContainer * constrs, Node *slave, Node *master, unsigned const &ndim, const string &which, const bool trsp, Function * fn) {
//     unsigned nDoFsPerNode;
//     if ( trsp ) {
//         // if master is transport, slave must be transport
//         if ( !dynamic_cast< TrsNode * >( slave ) ) {
//             return;
//         }
//         nDoFsPerNode = 1;
//     } else {
//         if ( !dynamic_cast< Particle * >( slave ) ) {
//             return;                                         // NOTE could be MechNode, but so far, nDoFs corresponds to Particles
//         }
//         nDoFsPerNode = 3 * ( ndim - 1 );
//     }
//
//
//     if ( slave->giveNumberOfDoFs() != master->giveNumberOfDoFs() ) {
//         std :: cerr << "slave and master must have the same number of DoFs, slave numDoFs = " << slave->giveNumberOfDoFs() << ", master numDoFs = " << master->giveNumberOfDoFs() << '\n';
//         exit(1);
//     }
//
//     // calculate multipliers and construct jointDof for every slaveDof
//     vector< vector< double > >tableOfMultipliers;
//     vector< vector< double > >tableOfTimeDepMultipliers;
//     vector< Node * >masterNodes;
//     vector< double >multipliers;
//     vector< double >time_multipliers;
//     vector< Function *> time_fns;
//     vector< unsigned >directions;
//
//     tableOfMultipliers.resize(nDoFsPerNode); // first index (j) is for master dof
//     tableOfTimeDepMultipliers.resize(nDoFsPerNode); // first index (j) is for master dof
//
//     for ( auto &mul : tableOfMultipliers ) {
//         mul.resize(nDoFsPerNode); // second index (i) is for slave dof
//     }
//     for ( auto &mul : tableOfTimeDepMultipliers ) {
//         mul = std :: vector< double > (nDoFsPerNode, 0.0); // second index (i) is for slave dof, initially zero
//     }
//     for ( unsigned i = 0; i < nDoFsPerNode; i++ ) {
//         for ( unsigned j = 0; j < nDoFsPerNode; j++ ) {
//             if ( i == j && i < ndim) {
//                 tableOfMultipliers [ j ] [ i ] = 1;
//                 if ( j == 0 ){
//                   // std::cout << "push_back_X" << '\n';
//                   tableOfTimeDepMultipliers [ j ] [ i ] = slave->givePoint().getX() - master->givePoint().getX();
//                 } else if ( j == 1 ){
//                   // std::cout << "push_back_Y" << '\n';
//                   tableOfTimeDepMultipliers [ j ] [ i ] = slave->givePoint().getY() - master->givePoint().getY();
//                 } else if ( j == 2 && j != nDoFsPerNode - 1){
//                   // std::cout << "push_back_Z" << '\n';
//                   tableOfTimeDepMultipliers [ j ] [ i ] = slave->givePoint().getZ() - master->givePoint().getZ();
//                 }
//             }
//             // else if ( i == 0 && j == nDoFsPerNode - 1 ) {
//             //     tableOfMultipliers [ j ] [ i ] = -( slave->givePoint().getY() - master->givePoint().getY() );
//             // } else if ( i == 1 && j == nDoFsPerNode - 1 ) {
//             //     tableOfMultipliers [ j ] [ i ] = slave->givePoint().getX() - master->givePoint().getX();
//             // } else if ( i == 0 && j == 4 ) {
//             //     tableOfMultipliers [ j ] [ i ] = slave->givePoint().getZ() - master->givePoint().getZ();
//             // } else if ( i == 1 && j == 3 ) {
//             //     tableOfMultipliers [ j ] [ i ] = -( slave->givePoint().getZ() - master->givePoint().getZ() );
//             // } else if ( i == 2 && j == 3 ) {
//             //     tableOfMultipliers [ j ] [ i ] = slave->givePoint().getY() - master->givePoint().getY();
//             // } else if ( i == 2 && j == 4 ) {
//             //     tableOfMultipliers [ j ] [ i ] =  -( slave->givePoint().getX() - master->givePoint().getX() );
//             // }
//             else {
//                 tableOfMultipliers [ j ] [ i ] =  0;
//             }
//         }
//     }
//     // read the table and put the information to jointDoFs
//     // NOTE this could be done in the previous loop
//     for ( unsigned i = 0; i < nDoFsPerNode; i++ ) {
//         // for transport nodes, only ones are in tableOfMultipliers
//         if ( !trsp ) {
//             if ( containsChar(which, 'x') && ( i == 0 || i == 4 || i == nDoFsPerNode - 1 ) ) {
//                 // std::cout << "fixed in x dir" << '\n';
//             } else if ( containsChar(which, 'y') && ( i == 1 || i == 3 || i == nDoFsPerNode - 1 ) ) {
//                 // std::cout << "fixed in y dir" << '\n';
//             } else if ( containsChar(which, 'z') && ( i == 2 || i == 3 || i == 4 ) && ndim > 2 ) {
//                 // std::cout << "fixed in z dir" << '\n';
//             } else {
//                 continue;
//             }
//         }
//         for ( unsigned j = 0; j < nDoFsPerNode; j++ ) {
//             // std::cout << "j = " << j << '\n';
//             if ( tableOfMultipliers [ j ] [ i ] != 0 || tableOfTimeDepMultipliers [ j ] [ i ] != 0 ) {
//                 masterNodes.push_back(master);
//                 multipliers.push_back(tableOfMultipliers [ j ] [ i ]);
//                 directions.push_back(j);
//                 if ( tableOfTimeDepMultipliers [ j ] [ i ] != 0 ) {
//                   // std::cout << "push_back " << '\n';
//                   time_multipliers.push_back( tableOfTimeDepMultipliers [ j ] [ i ] );
//                   time_fns.push_back( fn );
//                 } else {
//                   time_multipliers.push_back( 0.0 );
//                   time_fns.push_back( nullptr );
//                 }
//             }
//         }
//         if ( !masterNodes.empty() ) {
//             JointDoF *newJD = new JointDoF(slave, i, masterNodes, directions, multipliers, time_fns, time_multipliers);
//             constrs->addConstraint(newJD);
//         }
//         masterNodes.clear();
//         multipliers.clear();
//         directions.clear();
//         time_multipliers.clear();
//         time_fns.clear();
//     }
// }
