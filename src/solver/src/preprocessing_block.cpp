#include "preprocessing_block.h"
#include "geometry.h"
#include "model.h"
// #include "misc.h"  // TODO JK: this include causes linking error


std :: vector< double >PointToStdVector(const Point &p, unsigned dim = 3) {
    std :: vector< double >vect;
    vect.push_back( p.getX() );
    vect.push_back( p.getY() );
    if ( dim == 3 ) {
        vect.push_back( p.getZ() );
    }
    return vect;
}

bool containsChar(const std :: string &str, char c)
{
    // std::cout << "str = " << str << ", char = " << c << '\n';
    return str.find(c) != std :: string :: npos;
}

/*
 * Case Sensitive Implementation of endsWith()
 * It checks if the string 'mainStr' ends with given string 'toMatch'
 */
bool endsWith(const std :: string &mainStr, const std :: string &toMatch)
{
    if ( mainStr.size() >= toMatch.size() &&
         mainStr.compare(mainStr.size() - toMatch.size(), toMatch.size(), toMatch) == 0 ) {
        return true;
    } else {
        return false;
    }
}

//////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////
// MaterialRegion
void MaterialRegion :: readFromLine(istringstream &iss, unsigned d) {
    // enter line in pblockFile.inp:
    // MaterialRegion mech/trsp material_id id region 'block' x1 y1 z1 x2 y2 z2
    // coordinates refer to minimum and maximum point of block bounding box
    // TODO JK distinguish betwwen 2D and 3D and enable for use of cilinder region in 3D
    // dim = d;
    ( void ) d;
    std :: string param, region_name;
    while ( !iss.eof() ) {
        iss >> param;
        if ( param.compare("trsp") == 0 ) {
            this->transport = true;
        } else if ( param.compare("mech") == 0 ) {
            this->transport = false;
        } else if ( param.compare("material_id") == 0 || param.compare("materialID") == 0 ) {
            iss >> material_id;
        } else if ( param.compare("region") == 0 ) {
            iss >> region_name;
            if ( region_name.compare("block") == 0 ){
                double x, y, z, x2, y2, z2;
                iss >> x >> y >> z >> x2 >> y2 >> z2;
                block = Block(Point(x, y, z), Point(x2, y2, z2));
            // } else if ( region_name.compare("shpere") == 0 ){
            //     double x, y, z, r;
            //     reg = Sphere(Point(x, y, z), r);
            } else {
                std::cerr << "region named '" << region_name << "' not implemented yet" << '\n';
                exit(EXIT_FAILURE);
            }
        }
    }
}

void MaterialRegion :: apply(NodeContainer *nodes, ElementContainer *e, BCContainer *bcs, ConstraintContainer *constrs, FunctionContainer *funcs, ExporterContainer *ex, MaterialContainer *mats, Solver *solver) {
    ( void ) solver;
    ( void ) ex;
    ( void ) funcs;
    ( void ) constrs;
    ( void ) bcs;
    ( void ) nodes;

    // TODO make sure material is sutable for mech / transport

    for (auto & el : *e) {
        for (auto const & nod : el->giveNodes()){
            if (( !this->transport && nod->doesMechanics() )
                || ( this->transport && nod->doesTransport() )) {
                // TODO JK: are there any elems that are mechanical and connect transport nodes and same for transport elems?
                if ( this->block.isInside(nod->givePoint() ) ){
                    el->changeMaterial(mats->giveMaterial(this->material_id));
                    break;
                }
            }
        }
    }
}


//////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////
// Mechanical Periodic BC
void MechanicalPeriodicBC :: generateNewDoFs(NodeContainer *nodes) {
    //create new degrees of freedom representing strains ex, ey, gammaxy=2exy or ex, ey, ez, gammyz, gammaxz, gammaxy,
    MechDoF *mn;
    initalNodeNum = nodes->giveSize();
    mn = new MechDoF( dim, 3 * ( dim - 1 ) );
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

        if ( !nonsymmetric_shear ) {
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
        dBC.resize( m->giveNumberOfDoFs(), funcs->giveSize() );  //todo: warning C4267: 'argument': conversion from 'size_t' to 'const _Ty', possible loss of data
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
                vm.push_back( nodes->giveNode(n) );
            }
        }
        if ( vm.size() > 0 ) {
            unsigned nDoFs = 3;
            if ( dim == 3 ) {
                nDoFs = 6;
            }
            MechDoF *pn = new MechDoF(dim, nDoFs);
            nodes->addNode(pn);

            vector< unsigned >dirs( vm.size() );

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
void MechanicalPeriodicBC :: apply(NodeContainer *nodes, ElementContainer *e, BCContainer *bcs, ConstraintContainer *constrs, FunctionContainer *funcs, ExporterContainer *ex, MaterialContainer *mats, Solver *solver) {
    ( void ) e;
    ( void ) mats;
    ( void ) solver;

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
        } else if ( param.compare("nonsymmetric_shear") == 0 ) {
            nonsymmetric_shear = true;
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
double MechanicalPeriodicBC :: giveVolume() const {
    return volume;
}

//////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////
// Mechanical Periodic BC with Voigt's constraint
//////////////////////////////////////////////////////////
void MechanicalPeriodicBCwithVoigtConstraint :: generateConstraints(NodeContainer *nodes, ConstraintContainer *constrs) {
    //apply contraints, connect periodic images
    JointDoF *jd;
    vector< Node * >vm;
    vector< unsigned >dirs;
    vector< double >mults;
    Node *s = nullptr;
    Point diff;
    for ( unsigned n = 0; n < nodes->giveSize(); n++ ) {
        s = nodes->giveNode(n);
        if ( s->doesMechanics() && ( dynamic_cast< MechDoF * >( s ) == nullptr ) ) {
            //connect translations
            diff = s->givePoint();

            if ( !nonsymmetric_shear ) {
                //direction X  (all gammaxy and gammaxy realized here)
                if ( dim == 3 ) {
                    vm.resize(3);
                    mults.resize(3);
                    dirs.resize(3, 0);
                    dirs [ 1 ] = 5; //gamma xy
                    dirs [ 2 ] = 4; //gamma xz
                    vm [ 1 ] = nodes->giveNode(initalNodeNum);
                    vm [ 2 ] = nodes->giveNode(initalNodeNum); //gamma xz
                    mults [ 2 ] = diff.z / 2;
                } else if ( dim == 2 ) {
                    vm.resize(2);
                    mults.resize(2);
                    dirs.resize(2, 0);
                    dirs [ 1 ] = 2; //gamma xy
                    vm [ 1 ] = nodes->giveNode(initalNodeNum);
                }
                dirs [ 0 ] = 0; //eps x
                vm [ 0 ] = nodes->giveNode(initalNodeNum);
                mults [ 0 ] = diff.x;
                mults [ 1 ] = diff.y / 2;
                jd = new JointDoF(s, dirs [ 0 ], vm, dirs, mults);
                constrs->addConstraint(jd);

                //direction Y  (gammaxz realized here)
                if ( dim == 3 ) {
                    dirs [ 1 ] = 5; //gamma xy
                    dirs [ 2 ] = 3; //gamma yz
                    mults [ 2 ] = diff.z / 2;
                }
                dirs [ 0 ] = 1; //eps y
                mults [ 0 ] = diff.y;
                mults [ 1 ] = diff.x / 2;
                jd = new JointDoF(s, dirs [ 0 ], vm, dirs, mults);
                constrs->addConstraint(jd);

                //direction Z  (gammaxz realized here)
                if ( dim == 3 ) {
                    dirs [ 1 ] = 4; //gamma xz
                    dirs [ 2 ] = 3; //gamma yz
                    mults [ 2 ] = diff.y / 2;

                    dirs [ 0 ] = 2; //eps z
                    mults [ 0 ] = diff.z;
                    mults [ 1 ] = diff.x / 2;
                    jd = new JointDoF(s, dirs [ 0 ], vm, dirs, mults);
                    constrs->addConstraint(jd);
                }
            } else {
                //direction X  (all gammaxy and gammaxy realized here)
                if ( dim == 3 ) {
                    vm.resize(3);
                    mults.resize(3);
                    dirs.resize(3, 0);
                    dirs [ 1 ] = 5;  //gamma xy
                    dirs [ 2 ] = 4;  //gamma xz
                    vm [ 1 ] = nodes->giveNode(initalNodeNum);
                    vm [ 2 ] = nodes->giveNode(initalNodeNum);
                    mults [ 2 ] = diff.z;
                } else if ( dim == 2 ) {
                    vm.resize(2);
                    mults.resize(2);
                    dirs.resize(2, 0);
                    dirs [ 1 ] = 2; //gamma xy
                    vm [ 1 ] = nodes->giveNode(initalNodeNum);
                }
                dirs [ 0 ] = 0; //eps x
                vm [ 0 ] = nodes->giveNode(initalNodeNum);
                mults [ 0 ] = diff.x;
                mults [ 1 ] = diff.y;
                jd = new JointDoF(s, dirs [ 0 ], vm, dirs, mults);
                constrs->addConstraint(jd);

                //direction Y  (gammaxz realized here)
                if ( dim == 3 ) {
                    vm.resize(2);
                    mults.resize(2);
                    dirs.resize(2, 0);
                    dirs [ 1 ] = 3; //gamma yz
                    vm [ 1 ] = nodes->giveNode(initalNodeNum);
                    mults [ 1 ] = diff.z;
                } else if ( dim == 2 ) {
                    vm.resize(1);
                    mults.resize(1);
                    dirs.resize(1, 0);
                }
                dirs [ 0 ] = 1; //eps y
                vm [ 0 ] = nodes->giveNode(initalNodeNum);
                mults [ 0 ] = diff.y;
                jd = new JointDoF(s, dirs [ 0 ], vm, dirs, mults);
                constrs->addConstraint(jd);

                //direction Z  (gammaxz realized here)
                if ( dim == 3 ) {
                    vm.resize(1);
                    mults.resize(1);
                    dirs.resize(1, 0);
                    dirs [ 0 ] = 2;  //eps z
                    vm [ 0 ] = nodes->giveNode(initalNodeNum);
                    mults [ 0 ] = diff.z;
                    jd = new JointDoF(s, dirs [ 0 ], vm, dirs, mults);
                    constrs->addConstraint(jd);
                }
            }
        }
    }
}

//////////////////////////////////////////////////////////
void MechanicalPeriodicBCwithVoigtConstraint :: generateRigidBodyBC(NodeContainer *nodes, ElementContainer *elems, BCContainer *bcs, ConstraintContainer *constrs, FunctionContainer *funcs) {
    ( void ) elems;
    ( void ) constrs;

    if ( !nonsymmetric_shear ) {
        //all rotations are zero

        //add constant function
        vector< double >x, y;
        x.resize(1, 0);
        y.resize(1, 0);
        PieceWiseLinearFunction *newf = new PieceWiseLinearFunction(x, y);
        funcs->addFunction(newf);

        Node *node;
        vector< int >dBC, nBC;
        if ( dim == 2 ) {
            dBC.resize(3);
            nBC.resize(3);
            for ( unsigned i = 0; i < 3; i++ ) {
                dBC [ i ] = -1;
                nBC [ i ] = -1;
            }
            dBC [ 2 ] = funcs->giveSize() - 1; //rotation
        } else if ( dim == 3 ) {
            dBC.resize(6);
            nBC.resize(6);
            for ( unsigned i = 0; i < 6; i++ ) {
                if ( i < 3 ) {
                    dBC [ i ] = -1;
                } else {
                    dBC [ i ] = funcs->giveSize() - 1; //rotations
                }
                nBC [ i ] = -1;
            }
        }
        BoundaryCondition *bc;
        for ( unsigned n = 0; n < nodes->giveSize(); n++ ) {
            node = nodes->giveNode(n);
            if ( dynamic_cast< Particle * >( node ) != nullptr ) {
                bc = new BoundaryCondition(node, dBC, nBC);
                bcs->addBoundaryCondition(bc);
            }
        }
    } else {
        cerr << "Warning: MechanicalPeriodicBCwithVoigtConstraint should use avoid parameter 'nonsymmetric_shear' as the resulting rotations of the bodies might not be set to correct nonzero value and one can experience some differences" << endl;
        //exit(1);

        //nonzero rotations compensating shear strain applied assymetrically - Cosseart continuum
        //create new degrees of freedom connected to rotations
        MechDoF *mn;
        unsigned rotDoF = nodes->giveSize();
        unsigned rotnum =  2 * ( dim - 1 ) - 1;
        mn = new MechDoF( dim, rotnum);
        nodes->addNode(mn);
        //it is better to leave rotations free, setting them is not clear in case of shear stress loading
        /*
        BoundaryCondition *bc;
        vector< int >dBC, nBC;
        vector< double > bcmults;
        dBC.resize(rotnum);
        nBC.resize(rotnum);
        bcmults.resize(rotnum);
        for ( unsigned i = 0; i < rotnum; i++ ) {
            dBC[i] = strainFunc[dim+i];
            nBC[i] = stressFunc[dim+i];
            if(nBC[i]<0) bcmults[i] = 1;
            else bcmults[i] = volume;
        }
        bc = new BoundaryCondition(nodes->giveNode(rotDoF), dBC, nBC, bcmults);
        //bcs->addBoundaryCondition(bc);
        */

        JointDoF *jd;
        vector< Node * >vm(1);
        vm [ 0 ] = nodes->giveNode(rotDoF);
        vector< unsigned >dirs(1);
        vector< double >mults(1);
        Node *s = nullptr;
        for ( unsigned n = 0; n < nodes->giveSize(); n++ ) {
            s = nodes->giveNode(n);
            if ( dynamic_cast< Particle * >( s ) != nullptr ) {
                for ( unsigned i = 0; i < rotnum; i++ ) {
                    dirs [ 0 ] = i;
                    mults [ 0 ] = 0.5*pow(-1.,i+1);
                    jd = new JointDoF(s, dim + i, vm, dirs, mults);
                    constrs->addConstraint(jd);
                }
            }
        }
    }
}

//////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////
// Mechanical Periodic BC with Elastic constraint
//////////////////////////////////////////////////////////
void MechanicalPeriodicBCwithElasticConstraint :: apply(NodeContainer *nodes, ElementContainer *e, BCContainer *bcs, ConstraintContainer *constrs, FunctionContainer *funcs, ExporterContainer *ex, MaterialContainer *mats, Solver *solver) {
    ( void ) e;
    ( void ) mats;
    ( void ) solver;

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

    //new functions
    unsigned cfunc = funcs->giveSize() - 1;
    unsigned lfunc = funcs->giveSize();
    vector< double >x, y;
    x.resize(2);
    x [ 0 ] = 0;
    x [ 1 ] = 1;
    y.resize(2);
    y [ 0 ] = 0;
    y [ 1 ] = 1e-6;
    PieceWiseLinearFunction *newf = new PieceWiseLinearFunction(x, y);
    funcs->addFunction(newf);

    //set prescribed strain and stress
    vector< double >bcmults;
    BoundaryCondition *bc;
    vector< int >dBC, nBC;
    unsigned n = nodes->giveNode(initalNodeNum)->giveNumberOfDoFs();
    dBC.resize(n, cfunc);
    nBC.resize(n, -1);
    bcmults.resize(n, 1);
    bc = new BoundaryCondition(nodes->giveNode(initalNodeNum), dBC, nBC, bcmults);
    bcs->addBoundaryCondition(bc);

    //compute elastic solutions
    //Solver* oldSolver = masterModel->giveSolver();
    cout << "*** computing elastic solution on the periodic model" << endl;
    double dt = 1.;
    SteadyStateLinearSolver *linS = new SteadyStateLinearSolver();
    linS->setContainers( masterModel->giveElements(), masterModel->giveNodes(), masterModel->giveFunctions() );
    linS->setTimeStep(dt);
    linS->setInitialTimeStep(dt);
    masterModel->setSolver(linS);
    vector< Vector >elastSol(n);
    for ( unsigned i = 0; i < n; i++ ) {
        dBC [ i ] = lfunc;
        bc->replaceDirichBC(dBC);
        masterModel->init();
        linS->runBeforeEachStep();
        linS->solve();
        elastSol [ i ] = linS->giveTrialDoFValues()/(linS->giveTime()*1e-6);
        dBC [ i ] = cfunc;
    }
    delete linS; linS = nullptr;

    //remove added BC
    for ( int p = int( bcs->giveSize() ) - 1; p >= int( bcs_num ); p-- ) {
        bcs->removeBoundaryCondition(p);
    }
    //remove added constraints (this removes the master-slave constraint)
    for ( int p = int( constrs->giveSize() ) - 1; p >= int( const_num ); p-- ) {
        constrs->removeConstraint(p);
    }
    //remove added functions
    for ( int p = int( funcs->giveSize() ) - 1; p >= int( funcs_num ); p-- ) {
        funcs->removeFunction(p);
    }

    //set true BC
    for ( unsigned i = 0; i < n; i++ ) {
        if ( strainFunc [ i ] >= 0 ) {
            dBC [ i ] = strainFunc [ i ];
        } else {
            dBC [ i ] = -1;
        }
        if ( stressFunc [ i ] >= 0 ) {
            bcmults [ i ] = volume;
            nBC [ i ] = stressFunc [ i ];
        } else {
            nBC [ i ] = -1;
        }
        if ( strainFunc [ i ] >= 0 && stressFunc [ i ] >= 0 ) {
            cerr << "Error in Periodic boundary condition: cannot prescribe both stress and strain for the same direction" << endl;
        }
    }
    bc = new BoundaryCondition(nodes->giveNode(initalNodeNum), dBC, nBC, bcmults);
    bcs->addBoundaryCondition(bc);

    //create new constraints
    JointDoF *jd;
    vector< Node * >vm(n);
    vector< unsigned >dirs(n);
    vector< double >mults(n);
    for ( unsigned dir = 0; dir < n; dir++ ) {
        vm  [ dir ] = nodes->giveNode(initalNodeNum);
        dirs [ dir ] = dir;
    }

    Node *s = nullptr;
    unsigned DoFnum = 0;
    unsigned nodeDoFs;
    for ( unsigned nn = 0; nn < nodes->giveSize(); nn++ ) {
        s = nodes->giveNode(nn);
        nodeDoFs = s->giveNumberOfDoFs();
        if ( s->doesMechanics() && ( dynamic_cast< MechDoF * >( s ) == nullptr ) ) {
            for ( unsigned dir = 0; dir < nodeDoFs; dir++ ) {
                for ( unsigned k = 0; k < n; k++ ) {
                    mults [ k ] = elastSol [ k ] [ DoFnum + dir ];
                }
                jd = new JointDoF(s, dir, vm, dirs, mults);
                constrs->addConstraint(jd);
            }
        }
        DoFnum += nodeDoFs;
    }

    //masterModel->setSolver(oldSolver);
    masterModel->setSolver(solver);
    cout << "*** reseting solver and leaving preprocessing block" << endl;

    //export data
    generateExporters(nodes, ex);

    cout << "Applied periodic boundary conditions: " << nodes->giveSize() - initalNodeNum << " new DoFs (nodes " << initalNodeNum << " - " <<  nodes->giveSize() - 1 << "); " << constrs->giveSize() - const_num << " new constraints; " << bcs->giveSize() - bcs_num << " new boundary conditions; " << funcs->giveSize() - funcs_num << " new function; " << ex->giveSize() - ex_num << " new exporters; " << "created" << endl;
}


//////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////
// Transport Periodic BC
//////////////////////////////////////////////////////////
void TransportPeriodicBC :: generateNewDoFs(NodeContainer *nodes) {
    //create new degrees of freedom representing strains ex, ey, gammaxy=2exy or ex, ey, ez, gammyz, gammaxz, gammaxy,
    TrsDoF *mn;
    initalNodeNum = nodes->giveSize(); //todo warning C4267: '=': conversion from 'size_t' to 'unsigned int', possible loss of data
    mn = new TrsDoF(dim, dim);
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
        dBC.resize( m->giveNumberOfDoFs(), funcs->giveSize() );  //todo: conversion from 'size_t' to 'const _Ty', possible loss of data
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
                vm.push_back( nodes->giveNode(n) );
            }
        }
        if ( vm.size() > 0 ) {
            TrsDoF *tn = new TrsDoF(dim, 1);
            nodes->addNode(tn);

            vector< unsigned >dirs( vm.size() );
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
// Voigt's constraint
//////////////////////////////////////////////////////////
VoigtConstraint :: VoigtConstraint() {
    volume = 1;
}

//////////////////////////////////////////////////////////
VoigtConstraint :: ~VoigtConstraint() {}

//////////////////////////////////////////////////////////
void VoigtConstraint :: apply(NodeContainer *nodes, ElementContainer *elems, BCContainer *bcs, ConstraintContainer *constrs, FunctionContainer *funcs, ExporterContainer *ex, MaterialContainer *mats, Solver *solver) {
    ( void ) elems;
    ( void ) mats;
    ( void ) solver;

    MechDoF *master;
    unsigned masterNodeNum = nodes->giveSize();
    master = new MechDoF( dim, 3 * ( dim - 1 ) );
    nodes->addNode(master);

    //export data
    string export_name = "strain_stress";
    vector< unsigned >nn(1, masterNodeNum);
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
        fg = new ForceGauge(export_name, gname [ i ], codes [ i ], nn, nodes, 1. / volume, dim);
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
        dg = new DoFGauge(export_name, gname [ i ], codes [ i ], nn, nodes, 1., dim);
        ex->addExporter(dg);
    }

    //create new zero function for rotations
    unsigned funcnum = funcs->giveSize();
    vector< double >x, y;
    x.resize(1, 0);
    y.resize(1, 0);
    PieceWiseLinearFunction *newf = new PieceWiseLinearFunction(x, y);
    funcs->addFunction(newf);

    //set prescribed strain and stress
    vector< double >bcmults;
    BoundaryCondition *bc;
    vector< int >dBC, nBC;

    unsigned rotbcnum = ( 3 * ( dim - 1 ) );
    dBC.resize(rotbcnum, -1);
    nBC.resize(rotbcnum, -1);
    bcmults.resize(rotbcnum, 1);
    for ( unsigned i = dim; i < rotbcnum; i++ ) {
        dBC [ i ] = funcnum;
    }

    //apply contraints
    JointDoF *jd;
    vector< Node * >vm;
    vector< unsigned >dirs;
    vector< double >mults;
    Node *s = nullptr;
    Point diff;
    for ( unsigned n = 0; n < nodes->giveSize(); n++ ) {
        s = nodes->giveNode(n);
        if ( s->doesMechanics() && ( dynamic_cast< MechDoF * >( s ) == nullptr ) ) {
            //connect rotations
            if ( dynamic_cast< Particle * >( s ) ) {
                bc = new BoundaryCondition(s, dBC, nBC, bcmults);
                bcs->addBoundaryCondition(bc);
            }

            //connect translations
            diff = s->givePoint();

            //direction X  (all gammaxy and gammaxy realized here)
            if ( dim == 3 ) {
                vm.resize(3);
                mults.resize(3);
                dirs.resize(3, 0);
                dirs [ 1 ] = 5; //gamma xy
                dirs [ 2 ] = 4; //gamma xz
                vm [ 1 ] = nodes->giveNode(masterNodeNum);
                vm [ 2 ] = nodes->giveNode(masterNodeNum); //gamma xz
                mults [ 2 ] = diff.z / 2;
            } else if ( dim == 2 ) {
                vm.resize(2);
                mults.resize(2);
                dirs.resize(2, 0);
                dirs [ 1 ] = 2; //gamma xy
                vm [ 1 ] = nodes->giveNode(masterNodeNum);
            }
            dirs [ 0 ] = 0; //eps x
            vm [ 0 ] = nodes->giveNode(masterNodeNum);
            mults [ 0 ] = diff.x;
            mults [ 1 ] = diff.y / 2;
            jd = new JointDoF(s, dirs [ 0 ], vm, dirs, mults);
            constrs->addConstraint(jd);

            //direction Y  (gammaxz realized here)
            if ( dim == 3 ) {
                dirs [ 1 ] = 5; //gamma xy
                dirs [ 2 ] = 3; //gamma yz
                mults [ 2 ] = diff.z / 2;
            }
            dirs [ 0 ] = 1; //eps y
            mults [ 0 ] = diff.y;
            mults [ 1 ] = diff.x / 2;
            jd = new JointDoF(s, dirs [ 0 ], vm, dirs, mults);
            constrs->addConstraint(jd);

            //direction Z  (gammaxz realized here)
            if ( dim == 3 ) {
                dirs [ 1 ] = 4; //gamma xz
                dirs [ 2 ] = 3; //gamma yz
                mults [ 2 ] = diff.y / 2;

                dirs [ 0 ] = 2; //eps z
                mults [ 0 ] = diff.z;
                mults [ 1 ] = diff.x / 2;
                jd = new JointDoF(s, dirs [ 0 ], vm, dirs, mults);
                constrs->addConstraint(jd);
            }
        }
    }
    //set prescribed strain and stress
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
            cout << strainFunc [ i ] << " " << stressFunc [ i ] << endl;
            cerr << "Error in Voigt's constraint: cannot prescribe both stress and strain for the same direction" << endl;
        }
    }
    bc = new BoundaryCondition(nodes->giveNode(masterNodeNum), dBC, nBC, bcmults);
    bcs->addBoundaryCondition(bc);
}

//////////////////////////////////////////////////////////
void VoigtConstraint :: readFromLine(istringstream &iss, unsigned d) {
    dim = d;
    string param;
    unsigned num, hnum;

    bool volumeB, loadB;
    volumeB = loadB = false;

    while ( !iss.eof() ) {
        iss >> param;
        if ( param.compare("load") == 0 ) {
            loadB = true;
            iss >> num;

            strainFunc.resize(3 * ( dim - 1 ), -1);
            stressFunc.resize(3 * ( dim - 1 ), -1);

            for ( unsigned i = 0; i < num; i++ ) {
                iss >> param >> hnum;
                if ( dim == 2 ) {
                    std :: size_t found = param.find("z");
                    if ( found != std :: string :: npos ) {
                        cout << "Error in VoigtConstraint: cannot load by " << param << " in two dimensional setup" << '\n';
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
                    } else {
                        cout << "Error in VoigtConstraint: cannot load by " << param << " not implemented yet" << '\n';
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
                        cout << "Error in VoigtConstraint: cannot load by " << param << " not implemented yet" << '\n';
                        exit(1);
                    }
                }
            }
        } else if ( param.compare("volume") == 0 ) {
            volumeB = true;
            iss >> volume;
        }
    }

    if ( !volumeB ) {
        cout << "Error in VoigtConstraint: volume was not specified" << endl;
        exit(1);
    }
    if ( !loadB ) {
        cout << "Error in VoigtConstraint:  load was not specified" << endl;
        exit(1);
    }
}

//////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////
PressureFromMechanicalLoad :: PressureFromMechanicalLoad() {
    multiplier = 0.;
    master = 0;
    direction = masterdirection = 0;
    materialnum = 0;
};

//////////////////////////////////////////////////////////
PressureFromMechanicalLoad :: ~PressureFromMechanicalLoad() {};

//////////////////////////////////////////////////////////
void PressureFromMechanicalLoad :: apply(NodeContainer *n, ElementContainer *e, BCContainer *b, ConstraintContainer *c, FunctionContainer *funcs, ExporterContainer *ex, MaterialContainer *mats, Solver *solver) {
    ( void ) e;
    ( void ) b;
    ( void ) funcs;
    ( void ) ex;
    ( void ) solver;

    DiscreteTrsprtCoupledMaterial *dtcm = dynamic_cast< DiscreteTrsprtCoupledMaterial * >( mats->giveMaterial(materialnum) );
    if ( dtcm == nullptr ) {
        cerr << "Error in PressureFromMechanicalLoad: material is not of DiscreteTrsprtCoupledMaterial type" << endl;
        exit(1);
    }
    double biot = dtcm->giveBiotCoeff();
    if ( biot < 1e-5 ) {
        cerr << "Error in PressureFromMechanicalLoad: Biot coefficient " <<  biot << " is too small" << endl;
        exit(1);
    }


    vector< unsigned >directions(1);
    directions [ 0 ] = masterdirection;
    vector< Node * >masters(1);
    masters [ 0 ] = n->giveNode(master);
    vector< double >mults(1);
    mults [ 0 ] = multiplier / biot;

    DoFDependentOnConjugates *ddc;
    Node *trs;
    for ( auto &p: trsprtnodes ) {
        trs = n->giveNode(p);
        if ( !trs->isDoFTransport(direction) ) {
            cerr << "Error in PressureFromMechanicalLoad: node " << p << " direction " << direction <<  " is not pressure DoF" << endl;
            exit(1);
        }
        ddc = new DoFDependentOnConjugates(trs, direction, masters, directions, mults);
        c->addConstraint(ddc);
    }
}

//////////////////////////////////////////////////////////
void PressureFromMechanicalLoad :: readFromLine(istringstream &iss, unsigned d) {
    ( void ) d;
    iss >> master >> masterdirection >> materialnum >> multiplier >> direction;
    unsigned s;
    iss >> s;
    trsprtnodes.resize(s);
    for ( unsigned i = 0; i < s; i++ ) {
        iss >> trsprtnodes [ i ];
    }
}

//////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////

void RigidPlate :: setDirectionToFix(istringstream &iss) {

    string param;
    unsigned num;
    bool b;
    while ( !iss.eof() ) {
        iss >> param;
        if ( param.compare("dirs") == 0 ) {
            iss >> num;
            activeDirs.resize(num);
            for (unsigned i=0; i<num; i++) {
                iss >> b;
                activeDirs[i] = b;
            }
        }else if ( param.compare("which") == 0 ) {
            cerr << "RigidPlate error: you are using old unsupported way of specifying active directions" << endl;
            exit(1);
        }
    }
}

//////////////////////////////////////////////////////////
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

//////////////////////////////////////////////////////////
void RigidPlate :: checkMechTransport(Node *master) {
    // in case of rigid plate, master is a virtual virtual node and not a physical particle or
    if ( !endsWith(master->giveName(), "virtual") ) {
        master->setName( master->giveName().append("-virtual") );
    }

    if ( dynamic_cast< MechNode * >( master ) ) {
        if ( master->giveNumberOfDoFs() != ( 3 * ( this->dim - 1 ) ) ) {
            cerr << "Error in " << __func__ << ": Master for RigidPlate in mechnics must have " << ( 3 * ( this->dim - 1 ) ) << " DoFs, " << master->giveNumberOfDoFs() << " provided" << '\n';
            exit(EXIT_FAILURE);
        }
        if (activeDirs.size()==0) {
            activeDirs.resize( 3 * ( this->dim - 1 ), true );
        }
    } else if ( dynamic_cast< TrsNode * >( master ) ) {
        this->transport = true;
        if ( master->giveNumberOfDoFs() != 1 ) {
            cerr << "Error in " << __func__ << ": Master for RigidPlate in transport must have " << 1 << " DoFs, " << master->giveNumberOfDoFs() << " provided" << '\n';
            exit(EXIT_FAILURE);
        }
        if (activeDirs.size()==0) {
            activeDirs.resize( 1, true );
        }
    } else {
        cerr << "Error in " << __func__ << ": Master for RigidPlate is niether mechanical nor transport " << '\n';
        exit(EXIT_FAILURE);
    }
}

//////////////////////////////////////////////////////////
void RigidPlate :: apply(NodeContainer *nodes, ElementContainer *e, BCContainer *bcs, ConstraintContainer *constrs, FunctionContainer *funcs, ExporterContainer *ex, MaterialContainer *mats, Solver *solver) {
    ( void ) e;
    ( void ) bcs;
    ( void ) funcs;
    ( void ) ex;
    ( void ) mats;
    ( void ) solver;

    Node *master, *slave;
    //////////////////////////////////////////////////////////
    // read the line "masterId numSlaves slaveId1, slaveId2...."
    master = nodes->giveNode(master_id);
    // check if it is master node
    this->checkMechTransport(master);

    for ( auto const &sl_id : slave_ids ) {
        slave = nodes->giveNode(sl_id);
        connectSlaveMasterRigid(constrs, slave, master, this->dim, activeDirs, this->transport);
    }
}

//////////////////////////////////////////////////////////
void CoordRigidPlate :: readFromLine(istringstream &iss, unsigned d) {
    dim = d;
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

//////////////////////////////////////////////////////////
void CoordRigidPlate :: apply(NodeContainer *nodes, ElementContainer *e, BCContainer *bcs, ConstraintContainer *constrs, FunctionContainer *funcs, ExporterContainer *ex, MaterialContainer *mats, Solver *solver) {
    ( void ) e;
    ( void ) bcs;
    ( void ) funcs;
    ( void ) ex;
    ( void ) mats;
    ( void ) solver;

    // jointDoF jD;
    Node *master;

    master = nodes->giveNode(master_id);
    // check if it is master node
    RigidPlate :: checkMechTransport(master);

    for ( auto const &nod : * nodes ) {
        if ( isInBlock(nod->givePoint(), leftBottom, rightTop) ) {
            if ( nod == master || endsWith(nod->giveName(), "virtual")) {
                continue;
            }
            if (transport){
                if (!dynamic_cast< TrsNode* >(nod) || dynamic_cast< TrsDoF* >(nod)) continue;
            }else if (!dynamic_cast< MechNode* >(nod) || dynamic_cast< MechDoF* >(nod)) continue;
            connectSlaveMasterRigid(constrs, nod, master, this->dim, activeDirs, this->transport);
        }
    }
}


//////////////////////////////////////////////////////////
void RingRigidPlate :: readFromLine(istringstream &iss, unsigned d) {
    this->direction = 2;
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
        if ( dir.compare("dir") == 0 ) {
            iss >> this->direction;
        }
        // iss >> x >> y >> z;
        // axis = Point(x, y, z);
    }
    // else {
    //     // for 2D case, normal
    //     // axis = Point(0, 0, 1);
    // }

    RigidPlate :: setDirectionToFix(iss);
}

//////////////////////////////////////////////////////////
void RingRigidPlate :: apply(NodeContainer *nodes, ElementContainer *e, BCContainer *bcs, ConstraintContainer *constrs, FunctionContainer *funcs, ExporterContainer *ex, MaterialContainer *mats, Solver *solver) {
    ( void ) e;
    ( void ) bcs;
    ( void ) funcs;
    ( void ) ex;
    ( void ) mats;
    ( void ) solver;

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
        if ( isInCircle(node_point, this->center, this->r_outer, this->direction) ) {
            if ( !isInCircle(node_point, this->center, this->r_inner, this->direction) ) {
                if ( nod == master ) {
                    continue;
                }
                connectSlaveMasterRigid(constrs, nod, master, this->dim, activeDirs, this->transport);
            }
        }
    }
}


//////////////////////////////////////////////////////////
void ExpansionRing :: readFromLine(istringstream &iss, unsigned d) {
    RingRigidPlate :: readFromLine(iss, d);
    if ( typeid( this ) != typeid( ExpansionRing ) ) {
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

//////////////////////////////////////////////////////////
void ExpansionRing :: apply(NodeContainer *nodes, ElementContainer *e, BCContainer *bcs, ConstraintContainer *constrs, FunctionContainer *funcs, ExporterContainer *ex, MaterialContainer *mats, Solver *solver) {
    ( void ) e;
    ( void ) bcs;
    ( void ) funcs;
    ( void ) ex;
    ( void ) mats;
    ( void ) solver;

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
        if ( isInCircle(node_point, this->center, this->r_outer, this->direction) ) {
            if ( !isInCircle(node_point, this->center, this->r_inner, this->direction) ) {
                if ( nod == master || !dynamic_cast< MechNode * >( nod ) || dynamic_cast< MechDoF * >( nod )) {
                    continue;
                }
                connectSlaveMasterExpansion( constrs, nod, master, this->dim, this->transport, funcs->giveFunction(this->fn_id) );
            }
        }
    }
}


void ExpansionRingDoFLoad :: readFromLine(istringstream &iss, unsigned d) {
    RingRigidPlate :: readFromLine(iss, d);
    if ( typeid( this ) != typeid( ExpansionRing ) ) {
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

void ExpansionRingDoFLoad :: apply(NodeContainer *nodes, ElementContainer *e, BCContainer *bcs, ConstraintContainer *constrs, FunctionContainer *funcs, ExporterContainer *ex, MaterialContainer *mats, Solver *solver) {
    ( void ) e;
    ( void ) bcs;
    ( void ) funcs;
    ( void ) ex;
    ( void ) mats;
    ( void ) solver;

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
        if ( isInCircle(node_point, this->center, this->r_outer, this->direction) ) {
            if ( !isInCircle(node_point, this->center, this->r_inner, this->direction) ) {
                if ( nod == master ) {
                    continue;
                }
                connectSlaveMasterExpansionFLoad(constrs, nod, master, expMaster, this->dim);
            }
        }
    }
}


void ExpansionRingSingleDoFLoad :: apply(NodeContainer *nodes, ElementContainer *e, BCContainer *bcs, ConstraintContainer *constrs, FunctionContainer *funcs, ExporterContainer *ex, MaterialContainer *mats, Solver *solver) {
    ( void ) e;
    ( void ) bcs;
    ( void ) funcs;
    ( void ) ex;
    ( void ) mats;
    ( void ) solver;

    // jointDoF jD;
    Node *slave;
    Node *expMaster;

    expMaster = nodes->giveNode(this->master_id);

    // std::cout << expMaster->giveNumberOfDoFs() << " master point " << expMaster->giveName() << '\n';
    // expMaster->givePoint().print();

    Point node_point;
    int xm, ym, zm;
    xm = 1.;
    ym = 1.;
    zm = 1.;
    if ( this->direction == 0 ) {
        xm = 0;
    } else if ( this->direction == 1 ) {
        ym = 0;
    } else if ( this->direction == 2 ) {
        zm = 0;
    }

    // // NOTE JK this is intended for preferred choice of slave node (dof) according to max/min position in dir of user choice
    // unsigned x_max_id, y_max_id, z_max_id, x_min_id, y_min_id, z_min_id;
    // double x_max, y_max, z_max, x_min, y_min, z_min;

    unsigned slave_dir;
    unsigned num_nodes = 0;
    bool slave_used = false;

    double l_i;
    ( void ) l_i;           //It is used later but compiler still says it is not. Void is used to prevent warning.
    Point n_i;

    vector< Node * >masterNodes;
    vector< double >multipliers;
    vector< unsigned >directions;

    std :: vector< double >n_vect;
    double slave_dir_vect_value = 0;

    // n_vect_slave uložit dopředu
    // dělit všechno n_vect_slave dir
    //
    // zjistit na lineárním výpočtu to jak má být contraint

    MechNode * mn;
    MechDoF * md;
    this->center = Point(this->center.getX() * xm, this->center.getY() * ym, this->center.getZ() * zm);

    for ( auto const &nod : * nodes ) {
        mn = dynamic_cast < MechNode * >(nod);
        if ( mn==nullptr ) {
            continue;
        }
    	md = dynamic_cast < MechDoF * >(nod);
    	if ( md != nullptr ) continue;
    	if ( endsWith(nod->giveName(), "virtual") ) continue;

        node_point = Point(nod->givePoint().getX() * xm, nod->givePoint().getY() * ym, nod->givePoint().getZ() * zm);
        if ( isInCircle(node_point, this->center, this->r_outer, this->direction) ) {
            if ( !isInCircle(node_point, this->center, this->r_inner, this->direction) ) {
                // if ( nod == expMaster ) {
                //     continue;
                // }
                num_nodes++;
                // std::cout << "nod name = " << nod->giveName() << '\n';

                // n_i = ( Point((this->direction == 0) ? 0 : nod->givePoint().getX(),
                //               (this->direction == 1) ? 0 : nod->givePoint().getY(),
                //               (this->direction == 2) ? 0 : nod->givePoint().getZ()) -
                //       Point((this->direction == 0) ? 0 : this->center.getX(),
                //                     (this->direction == 1) ? 0 : this->center.getY(),
                //                     (this->direction == 2) ? 0 : this->center.getZ())
                //         );
                n_i = (node_point - this->center);
                l_i = n_i.norm();
                n_i.normalize();
                // std::cout << "dim = " << this->dim << '\n';
                n_vect = PointToStdVector(n_i, this->dim);
                if ( slave_used ) {
                    for ( unsigned i = 0; i < this->dim; i++ ) {
                        // std::cout << "i = " << i << ", num_nodes = " << num_nodes << '\n';
                        masterNodes.push_back(nod);
                        directions.push_back(i);
                        multipliers.push_back(-n_vect [ i ]);
                    }
                    //cout << "MASTER " << nod << " " << nod->giveName() << " " << n_vect [ 0 ] << " " << n_vect [ 1 ] << " " << n_vect [ 2 ] << endl;
                } else {
                    // first node is taken as a slave
                    slave = nod;
                    if ( * std :: max_element( n_vect.begin(), n_vect.end() ) >  abs(* std :: min_element( n_vect.begin(), n_vect.end() ) ) ) {
                        slave_dir = std :: distance( n_vect.begin(), std :: max_element( n_vect.begin(), n_vect.end() ) );
                    } else {
                        slave_dir = std :: distance( n_vect.begin(), std :: min_element( n_vect.begin(), n_vect.end() ) );
                    }
                    slave_dir_vect_value = n_vect [ slave_dir ];
                    for ( unsigned i = 0; i < this->dim; i++ ) {
                        if ( i != slave_dir ) {
                            // multiplying self other DoFs
                            multipliers.push_back(-n_vect [ i ]);
                            masterNodes.push_back(slave);
                            directions.push_back(i);
                        }
                    }
                    slave_used = true;
                    //cout << "SLAVE " << slave->giveName() << " " << n_vect [ 0 ] << " " << n_vect [ 1 ] << " " << n_vect [ 2 ] << " " <<  slave_dir << endl;
                }
            }
        }
    }

    // adding master DoF governing the expansion
    multipliers.push_back(double( num_nodes ) * (r_outer+r_inner) / 4.);
    masterNodes.push_back(expMaster);
    directions.push_back(0);

    if ( !masterNodes.empty() ) {
        for ( unsigned j = 0; j < multipliers.size(); j++ ) {
            multipliers [ j ] /= slave_dir_vect_value;
        }
        JointDoF *newJD = new JointDoF(slave, slave_dir, masterNodes, directions, multipliers);
            for(unsigned u=0; u<multipliers.size(); u++) if(multipliers[u]!=multipliers[u]) {cout << "multipliersX ERROR" << endl; exit(1);}
        constrs->addConstraint(newJD);
    }

    masterNodes.clear();
    multipliers.clear();
    directions.clear();
    std::cout << "Expansion volumetric load applied: center(" << this->center.getX() << ", " << this->center.getY() << ", " << this->center.getZ() << ") rI = " << this->r_inner << ", rO = " << this->r_outer << ", direction: " << this->direction <<  '\n';
    // exit(0);
}

//////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////
// CONTAINER FOR PBLOCKS
PBlockContainer :: ~PBlockContainer() {
    for ( vector< PBlock * > :: iterator f = blocks.begin(); f != blocks.end(); ++f ) {
        if ( * f != nullptr ) {
            delete * f;
        }
    }
}

//////////////////////////////////////////////////////////
void PBlockContainer :: clear() {
    for ( vector< PBlock * > :: iterator f = blocks.begin(); f != blocks.end(); ++f ) {
        if ( * f != nullptr ) {
            delete * f;
        }
    }
}

//////////////////////////////////////////////////////////
void PBlockContainer :: setContainers(NodeContainer *n, ElementContainer *e, BCContainer *b, ConstraintContainer *c, FunctionContainer *f, ExporterContainer *ex, MaterialContainer *mats, Solver *solv) {
    nodes = n;
    elems = e;
    bcs = b;
    constrs = c;
    funcs = f;
    exporters = ex;
    materials = mats;
    solver = solv;
}

//////////////////////////////////////////////////////////
void PBlockContainer :: apply() {
    for ( auto b: blocks ) {
        b->apply(nodes, elems, bcs, constrs, funcs, exporters, materials, solver);
    }
}

//////////////////////////////////////////////////////////
void PBlockContainer :: readFromFile(const string filename, unsigned dim) {
    unsigned origsize = blocks.size(); //todo: warning C4267: 'initializing': conversion from 'size_t' to 'unsigned int', possible loss of dat
    string line, ftype;
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
            iss >> ftype;
            if ( !(ftype.rfind("#", 0) == 0) ) {
                if ( ftype.compare("MechanicalPeriodicBC") == 0 ) {
                    MechanicalPeriodicBC *newblock = new MechanicalPeriodicBC();
                    newblock->readFromLine(iss, dim);
                    blocks.push_back(newblock);
                } else if ( ftype.compare("MechanicalPeriodicBCwithVoigtConstraint") == 0 ) {
                    MechanicalPeriodicBCwithVoigtConstraint *newblock = new MechanicalPeriodicBCwithVoigtConstraint();
                    newblock->readFromLine(iss, dim);
                    blocks.push_back(newblock);
                } else if ( ftype.compare("MechanicalPeriodicBCwithElasticConstraint") == 0 ) {
                    MechanicalPeriodicBCwithElasticConstraint *newblock = new MechanicalPeriodicBCwithElasticConstraint();
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
                } else if ( ftype.compare("VoigtConstraint") == 0 ) {
                    VoigtConstraint *newblock = new VoigtConstraint();
                    newblock->readFromLine(iss, dim);
                    blocks.push_back(newblock);
                } else if ( ftype.compare("PressureFromMechanicalLoad") == 0 ) {
                    PressureFromMechanicalLoad *newblock = new PressureFromMechanicalLoad();
                    newblock->readFromLine(iss, dim);
                    blocks.push_back(newblock);
                } else if ( ftype.compare("MaterialRegion") == 0 ) {
                    MaterialRegion *newblock = new MaterialRegion();
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



void connectSlaveMasterRigid(ConstraintContainer *constrs, Node *slave, Node *master, unsigned const &ndim, const vector < bool > &activeDirs, const bool transport) {
    unsigned nDoFsPerNode;
    if ( transport ) {
        nDoFsPerNode = 1;
    } else {
        nDoFsPerNode = 3 * ( ndim - 1 );
    }


    if ( nDoFsPerNode != master->giveNumberOfDoFs()) {
        std :: cerr << "RigidPlate Error: master node with " << master->giveNumberOfDoFs() << " DoFs should have " << nDoFsPerNode << "DoFs instead" << endl;
        exit(1);
    }
    if ( nDoFsPerNode != activeDirs.size() ) {
        std :: cerr << "RigidPlate Error: active directions with " << nDoFsPerNode << " entries should have " << nDoFsPerNode << " entries instead" << endl;
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
        if ( !activeDirs[i] ) continue;
        //if (i==1 || i==2) continue;
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
    if ( * std :: max_element( n_vect.begin(), n_vect.end() ) >  abs(* std :: min_element( n_vect.begin(), n_vect.end() ) ) ) {
        slave_dir = std :: distance( n_vect.begin(), std :: max_element( n_vect.begin(), n_vect.end() ) );
    } else {
        slave_dir = std :: distance( n_vect.begin(), std :: min_element( n_vect.begin(), n_vect.end() ) );
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
    if ( * std :: max_element( n_vect.begin(), n_vect.end() ) >  abs(* std :: min_element( n_vect.begin(), n_vect.end() ) ) ) {
        slave_dir = std :: distance( n_vect.begin(), std :: max_element( n_vect.begin(), n_vect.end() ) );
    } else {
        slave_dir = std :: distance( n_vect.begin(), std :: min_element( n_vect.begin(), n_vect.end() ) );
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
