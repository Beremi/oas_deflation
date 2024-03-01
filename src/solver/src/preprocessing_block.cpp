#include "preprocessing_block.h"
#include "model.h"
#include "material_vectorial.h"
#include "element_ldpm.h"
#include "solver_implicit.h"

// #include "misc.h"  // TODO JK: this include causes linking error
using namespace std;


std :: vector< double >PointToStdVector(const Point &p, unsigned dim = 3) {
    std :: vector< double >vect;
    vect.push_back(p.x() );
    vect.push_back(p.y() );
    if ( dim == 3 ) {
        vect.push_back(p.z() );
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
    // MaterialRegion mech/trsp material_id id region 'sphere' x1 y1 z1 r
    // MaterialRegion mech/trsp material_id id region 'cylinder' x1 y1 r
    // important NOTE: in 2D model, use only x1 y1 x2 y2


    std :: string param, region_name;
    while (  iss >> param ) {
        if ( param.compare("trsp") == 0 ) {
            this->transport = true;
        } else if ( param.compare("mech") == 0 ) {
            this->transport = false;
        } else if ( param.compare("material_id") == 0 || param.compare("materialID") == 0 ) {
            iss >> material_id;
        } else if ( param.compare("all_nodes") == 0 ) {
            this->all_nodes = true;
        } else if ( param.compare("region") == 0 ) {
            iss >> region_name;
            if ( region_name.compare("block") == 0 ) {
                if ( d == 2 ) {
                    double x, y, x2, y2;
                    iss >> x >> y >> x2 >> y2;
                    reg = new Block(Point(x, y, -1e10), Point(x2, y2, 1e10), 2);
                } else if ( d == 3 ) {
                    double x, y, z, x2, y2, z2;
                    iss >> x >> y >> z >> x2 >> y2 >> z2;
                    reg = new Block(Point(x, y, z), Point(x2, y2, z2), 3);
                }
            } else if ( region_name.compare("shpere") == 0 ) {
                if ( d == 2 ) {
                    double x, y, r;
                    iss >> x >> y >> r;
                    reg = new Sphere(Point(x, y, 0), r);
                } else if ( d == 3 ) {
                    double x, y, z, r;
                    iss >> x >> y >> z >> r;
                    reg = new Sphere(Point(x, y, z), r);
                }
            } else if ( region_name.compare("cylinder") == 0 || region_name.compare("circle") == 0 ) {
                // aplicable in 2D or in 3D as cylinder along z axis
                double x, y, r;
                iss >> x >> y >> r;
                reg = new Circle(Point(x, y, 0), r);
            } else {
                std :: cerr << "region named '" << region_name << "' not implemented yet" << '\n';
                exit(EXIT_FAILURE);
            }
        }
    }
}

void MaterialRegion :: apply(NodeContainer *nodes, ElementContainer *elems, BCContainer *bcs, ConstraintContainer *constrs, FunctionContainer *funcs, ExporterContainer *ex, MaterialContainer *mats, RegionContainer *regions, Solver *solver) {
    ( void ) solver;
    ( void ) regions;
    ( void ) ex;
    ( void ) funcs;
    ( void ) constrs;
    ( void ) bcs;
    ( void ) nodes;

    // TODO make sure material is sutable for mech / transport
    unsigned num_nodes_inside;
    for ( auto &el : * elems ) {
        num_nodes_inside = 0;
        for ( auto const &nod : el->giveNodes() ) {
            if ( ( !this->transport && nod->doesMechanics() ) ||
                 ( this->transport && nod->doesTransport() ) ) {
                // TODO JK: are there any elems that are mechanical and connect transport nodes and same for transport elems?
                if ( this->reg->isInside( nod->givePoint() ) ) {
                    if ( this->all_nodes ) {
                        num_nodes_inside++;
                    } else {
                        el->changeMaterial( mats->giveMaterial(this->material_id) );
                        break;
                    }
                }
            }
        }
        if ( this->all_nodes && num_nodes_inside == el->giveNodes().size() ) {
            el->changeMaterial( mats->giveMaterial(this->material_id) );
        }
    }
}


//////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////
// Mechanical Periodic BC
//////////////////////////////////////////////////////////
void MechanicalPeriodicBC :: generateNewDoFs(NodeContainer *nodes) {
    //create new degrees of freedom representing strains ex, ey, gammaxy=2exy or ex, ey, ez, gammyz, gammaxz, gammaxy,
    MechDoF *mn;
    initalNodeNum = nodes->giveSize();
    mn = new MechDoF(dim, 3 * ( dim - 1 ) );
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
                mults [ 3 ] = diff.z() / 2;
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
            mults [ 1 ] = diff.x();
            mults [ 2 ] = diff.y() / 2;
            jd = new JointDoF(s, dirs [ 0 ], vm, dirs, mults);
            constrs->addConstraint(jd);

            //direction Y  (gammaxz realized here)
            if ( dim == 3 ) {
                dirs [ 2 ] = 5; //gamma xy
                dirs [ 3 ] = 3; //gamma yz
                mults [ 3 ] = diff.z() / 2;
            }
            dirs [ 0 ] = 1;
            dirs [ 1 ] = 1; //eps y
            mults [ 1 ] = diff.y();
            mults [ 2 ] = diff.x() / 2;
            jd = new JointDoF(s, dirs [ 0 ], vm, dirs, mults);
            constrs->addConstraint(jd);

            //direction Z  (gammaxz realized here)
            if ( dim == 3 ) {
                dirs [ 2 ] = 4; //gamma xz
                dirs [ 3 ] = 3; //gamma yz
                mults [ 3 ] = diff.y() / 2;

                dirs [ 0 ] = 2;
                dirs [ 1 ] = 2; //eps z
                mults [ 1 ] = diff.z();
                mults [ 2 ] = diff.x() / 2;
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
                mults [ 3 ] = diff.z();
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
            mults [ 1 ] = diff.x();
            mults [ 2 ] = diff.y();
            jd = new JointDoF(s, dirs [ 0 ], vm, dirs, mults);
            constrs->addConstraint(jd);

            //direction Y  (gammaxz realized here)
            if ( dim == 3 ) {
                vm.resize(3);
                mults.resize(3);
                dirs.resize(3, 0);
                dirs [ 2 ] = 3; //gamma yz
                vm [ 2 ] = nodes->giveNode(initalNodeNum);
                mults [ 2 ] = diff.z();
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
            mults [ 1 ] = diff.y();
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
                mults [ 1 ] = diff.z();
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
        dBC.resize(m->giveNumberOfDoFs(), funcs->giveSize() );       //todo: warning C4267: 'argument': conversion from 'size_t' to 'const _Ty', possible loss of data
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
            MechDoF *pn = new MechDoF(dim, nDoFs);
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
void MechanicalPeriodicBC :: apply(NodeContainer *nodes, ElementContainer *elems, BCContainer *bcs, ConstraintContainer *constrs, FunctionContainer *funcs, ExporterContainer *ex, MaterialContainer *mats, RegionContainer *regions, Solver *solver) {
    ( void ) elems;
    ( void ) mats;
    ( void ) solver;
    ( void ) regions;

    calculateVolume(elems);

    unsigned const_num = constrs->giveSize();
    unsigned funcs_num = funcs->giveSize();
    unsigned bcs_num = bcs->giveSize();
    unsigned ex_num = ex->giveSize();

    //create new degrees of freedom representing strains ex, ey, gammaxy=2exy or ex, ey, ez, gammyz, gammaxz, gammaxy,
    generateNewDoFs(nodes);

    //apply contraints, connect periodic images
    generateConstraints(nodes, constrs);

    //boundary conditions
    generateRigidBodyBC(nodes, elems, bcs, constrs, funcs);

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

    while (  iss >> param ) {
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
void MechanicalPeriodicBC :: calculateVolume(ElementContainer *elems) {
    /*volume = 1;
     * for ( auto const a : PUCsize ) {
     *  volume *= a;
     * }
     */
    //THIS IS INCORRECT - WILL NOT WORK FOR COUPLED PROBLEMS WHERE ELEMENTS MIGHT OVERLAP
    volume = 0;
    for ( unsigned i = 0; i < elems->giveSize(); i++ ) {
        volume += elems->giveElement(i)->giveVolume();
    }
}


//////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////
// Mechanical Periodic BC on SPHERE
//////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////
void MechanicalSphericalPeriodicBC :: generateConstraints(NodeContainer *nodes, ConstraintContainer *constrs) {
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
            mults [ 0 ] = -1;
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
                mults [ 3 ] = diff.z() / 2;
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
            mults [ 1 ] = diff.x();
            mults [ 2 ] = diff.y() / 2;
            jd = new JointDoF(s, dirs [ 0 ], vm, dirs, mults);
            constrs->addConstraint(jd);

            //direction Y  (gammaxz realized here)
            if ( dim == 3 ) {
                dirs [ 2 ] = 5; //gamma xy
                dirs [ 3 ] = 3; //gamma yz
                mults [ 3 ] = diff.z() / 2;
            }
            dirs [ 0 ] = 1;
            dirs [ 1 ] = 1; //eps y
            mults [ 1 ] = diff.y();
            mults [ 2 ] = diff.x() / 2;
            jd = new JointDoF(s, dirs [ 0 ], vm, dirs, mults);
            constrs->addConstraint(jd);

            //direction Z  (gammaxz realized here)
            if ( dim == 3 ) {
                dirs [ 2 ] = 4; //gamma xz
                dirs [ 3 ] = 3; //gamma yz
                mults [ 3 ] = diff.y() / 2;

                dirs [ 0 ] = 2;
                dirs [ 1 ] = 2; //eps z
                mults [ 1 ] = diff.z();
                mults [ 2 ] = diff.x() / 2;
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
                mults [ 3 ] = diff.z();
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
            mults [ 1 ] = diff.x();
            mults [ 2 ] = diff.y();
            jd = new JointDoF(s, dirs [ 0 ], vm, dirs, mults);
            constrs->addConstraint(jd);

            //direction Y  (gammaxz realized here)
            if ( dim == 3 ) {
                vm.resize(3);
                mults.resize(3);
                dirs.resize(3, 0);
                dirs [ 2 ] = 3; //gamma yz
                vm [ 2 ] = nodes->giveNode(initalNodeNum);
                mults [ 2 ] = diff.z();
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
            mults [ 1 ] = diff.y();
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
                mults [ 1 ] = diff.z();
                jd = new JointDoF(s, dirs [ 0 ], vm, dirs, mults);
                constrs->addConstraint(jd);
            }
        }
    }
}

//////////////////////////////////////////////////////////
void MechanicalSphericalPeriodicBC :: calculateVolume(ElementContainer *elems) {
    MechanicalPeriodicBC :: calculateVolume(elems);
}

//////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////
// Mechanical Periodic BC on SPHERE experimental
// uM.n = uS.n & uM.t = uS.(-t)
// rotations not finished
//////////////////////////////////////////////////////////
void MechanicalSphericalPeriodicBCExperimental :: generateConstraints(NodeContainer *nodes, ConstraintContainer *constrs) {
    //apply contraints, connect periodic images
    Node *s = nullptr;
    Node *m = nullptr;

    int constrained_rots = 0;   // counter so that constrainRotations is applied to one pair only

    for ( unsigned i = 0; i < masters.size(); i++ ) {
        // get coords
        m = nodes->giveNode(masters [ i ]);
        s = nodes->giveNode(slaves [ i ]);
        Point m_coords = m->givePoint();
        Point s_coords = s->givePoint();
        // get direction vectors n & t
        Point n = m_coords - s_coords;
        n /= n.norm();
        Point t;
        t [ 0 ] = -n [ 1 ];
        t [ 1 ] = n [ 0 ];
        if ( 0.1 < abs(n [ 0 ]) && abs(n [ 0 ]) < 0.4  && constrained_rots == 0 ) {         // applied to only one pair!
            constrained_rots = 1;
            constrainRotation(nodes, constrs, m, s, n, t);
        } else {
            constrainRegular(nodes, constrs, m, s, n, t);
        }
    }
}

//////////////////////////////////////////////////////////
void MechanicalSphericalPeriodicBCExperimental :: generateRigidBodyBC(NodeContainer *nodes, ElementContainer *elems, BCContainer *bcs, ConstraintContainer *constrs, FunctionContainer *funcs) {
    ( void ) nodes;
    ( void ) elems;
    ( void ) bcs;
    ( void ) constrs;
    ( void ) funcs;
    // substituted by constraints on one master-slave pair in generateConstraints which gets restricted in tangentional direction (constrainRotations)
}

//////////////////////////////////////////////////////////
void MechanicalSphericalPeriodicBCExperimental :: constrainRegular(NodeContainer *nodes, ConstraintContainer *constrs, Node *m, Node *s, Point n, Point t) {
    Point m_coords = m->givePoint();
    Point s_coords = s->givePoint();
    JointDoF *jd;
    vector< Node * >vm;
    vector< unsigned >dirs;
    vector< double >mults;
    // connect rotations
    if ( dynamic_cast< Particle * >( s ) && dynamic_cast< Particle * >( m ) ) {
        dirs.resize(1);
        mults.resize(1);
        vm.resize(1);
        mults [ 0 ] = -1;
        vm [ 0 ] = m;
        for ( unsigned k = 0; k < 2 * ( dim - 1 ) - 1; k++ ) {  // general, for 2D not necessary
            dirs [ 0 ] = dim + k;
            jd = new JointDoF(s, dirs [ 0 ], vm, dirs, mults);
            constrs->addConstraint(jd);
        }
    }

    double diam = sqrt(pow(m_coords [ 0 ] - s_coords [ 0 ], 2) + pow(m_coords [ 1 ] - s_coords [ 1 ], 2) );

    dirs.resize(4);
    mults.resize(4);
    vm.resize(4);

    // calculate weights
    // note: A = - D and B = C (perpendicular)
    double A = ( n [ 0 ] * t [ 1 ] + n [ 1 ] * t [ 0 ] ) / ( n [ 0 ] * t [ 1 ] - n [ 1 ] * t [ 0 ] );
    double B = ( 2 * n [ 1 ] * t [ 1 ] ) / ( n [ 0 ] * t [ 1 ] - n [ 1 ] * t [ 0 ] );
    double C = ( 2 * n [ 0 ] * t [ 0 ] ) / ( n [ 1 ] * t [ 0 ] - n [ 0 ] * t [ 1 ] );
    double D = ( n [ 1 ] * t [ 0 ] + n [ 0 ] * t [ 1 ] ) / ( n [ 1 ] * t [ 0 ] - n [ 0 ] * t [ 1 ] );

    // 2D ONLY
    // M refers to macroscopic strain vector w/ symmetrical shear (epsx, epsy, gammayxy), gotten via nodes->giveNode(initalNodeNum)
    // m refers to displacement vector of the same node (ux, uy, ..., rotx,... ) 6 members?
    vm [ 0 ] = m;
    vm [ 1 ] = m;
    vm [ 2 ] = nodes->giveNode(initalNodeNum);  // M
    vm [ 3 ] = nodes->giveNode(initalNodeNum);  // M
    dirs [ 0 ] = 0;     // position where it takes ux value from m
    dirs [ 1 ] = 1;     // position where it takes uy value from m

    // direction X
    // ux_S = A * ux_M + B * uy_M + (epsx*nx*2*r + gammaxy*ny*2*r)
    //      = mults [ 0 ] * m ( 0 ) + mults [ 1 ] * m ( 1 ) - (M ( 0 ) * mults [ 2 ] + M ( 2 ) * mults [ 3 ])
    dirs [ 2 ] = 0;     // position where it takes epsx value from M
    dirs [ 3 ] = 2;     // position where it takes gammaxy value from M
    mults [ 0 ] = A;
    mults [ 1 ] = B;
    mults [ 2 ] = -n [ 0 ] * diam;
    mults [ 3 ] = -n [ 1 ] * diam;
    jd = new JointDoF(s, 0, vm, dirs, mults);
    constrs->addConstraint(jd);

    // direction Y
    // uy_S = C * ux_M + D * uy_M + (epsy*ny*2*r + gammaxy*nx*2*r)
    //      = mults [ 0 ] * m ( 0 ) + mults [ 1 ] * m ( 1 ) - (M ( 0 ) * mults [ 2 ] + M ( 2 ) * mults [ 3 ])
    dirs [ 2 ] = 1;     // position where it takes epsx value from M
    dirs [ 3 ] = 2;     // position where it takes gammaxy value from M
    mults [ 0 ] = C;
    mults [ 1 ] = D;
    mults [ 2 ] = -n [ 1 ] * diam;
    mults [ 3 ] = -n [ 0 ] * diam;
    jd = new JointDoF(s, 1, vm, dirs, mults);
    constrs->addConstraint(jd);
}

//////////////////////////////////////////////////////////
void MechanicalSphericalPeriodicBCExperimental :: constrainRotation(NodeContainer *nodes, ConstraintContainer *constrs, Node *m, Node *s, Point n, Point t) {
    JointDoF *jd;
    Point m_coords = m->givePoint();
    Point s_coords = s->givePoint();
    vector< Node * >vm;
    vector< unsigned >dirs;
    vector< double >mults;
    double diam = sqrt(pow(m_coords [ 0 ] - s_coords [ 0 ], 2) + pow(m_coords [ 1 ] - s_coords [ 1 ], 2) );
    // constrain movement of master node in the direction of n
    dirs.resize(1);
    mults.resize(1);
    vm.resize(1);
    dirs [ 0 ] = 1;
    mults [ 0 ] = n [ 0 ] / n [ 1 ];    // |u|*nx = ux, |u|*ny = uy -> ux = (nx/ny) * uy
    vm [ 0 ] = m;
    jd = new JointDoF(m, 0, vm, dirs, mults);
    constrs->addConstraint(jd);

    // constrain rotations (tangentional movement of the master-slave pair)
    dirs.resize(3);
    mults.resize(3);
    vm.resize(3);
    // calculate weights
    double AB = ( n [ 0 ] / n [ 1 ] ) * ( ( n [ 0 ] * t [ 1 ] + n [ 1 ] * t [ 0 ] ) / ( n [ 0 ] * t [ 1 ] - n [ 1 ] * t [ 0 ] ) ) +  ( 2 * n [ 1 ] * t [ 1 ] ) / ( n [ 0 ] * t [ 1 ] - n [ 1 ] * t [ 0 ] );
    double CD = ( n [ 0 ] / n [ 1 ] ) * ( 2 * n [ 0 ] * t [ 0 ] ) / ( n [ 1 ] * t [ 0 ] - n [ 0 ] * t [ 1 ] ) + ( n [ 1 ] * t [ 0 ] + n [ 0 ] * t [ 1 ] ) / ( n [ 1 ] * t [ 0 ] - n [ 0 ] * t [ 1 ] );
    // 2D ONLY
    // M refers to macroscopic strain vector w/ symmetrical shear (epsx, epsy, gammayxy), gotten via nodes->giveNode(initalNodeNum)
    // m refers to displacement vector of the same node (ux, uy, ..., rotx,... )
    vm [ 0 ] = m;
    vm [ 1 ] = nodes->giveNode(initalNodeNum);  // M
    vm [ 2 ] = nodes->giveNode(initalNodeNum);
    dirs [ 0 ] = 1;     // takes only uy value from m (ux value is constrained by weights above)
    dirs [ 2 ] = 2;     // position where it takes gammaxy value from M

    // direction X
    // ux_S = ( A * (nx/ny) + B ) * uy_M + (epsx*nx*2*r + gammaxy*ny*2*r)
    dirs [ 1 ] = 0;     // position where it takes epsx value from M
    mults [ 0 ] = AB;
    mults [ 1 ] = -n [ 0 ] * diam;
    mults [ 2 ] = -n [ 1 ] * diam;
    jd = new JointDoF(s, 0, vm, dirs, mults);
    constrs->addConstraint(jd);

    // direction Y
    // uy_S = ( C * (nx/ny) + D ) * uy_M + (epsy*ny*2*r + gammaxy*nx*2*r)
    dirs [ 1 ] = 1;     // position where it takes epsy value from M
    mults [ 0 ] = CD;
    mults [ 1 ] = -n [ 1 ] * diam;
    mults [ 2 ] = -n [ 0 ] * diam;
    jd = new JointDoF(s, 1, vm, dirs, mults);
    constrs->addConstraint(jd);
}

//////////////////////////////////////////////////////////
void MechanicalSphericalPeriodicBCExperimental :: calculateVolume(ElementContainer *elems) {
    MechanicalPeriodicBC :: calculateVolume(elems);
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
                    mults [ 2 ] = diff.z() / 2;
                } else if ( dim == 2 ) {
                    vm.resize(2);
                    mults.resize(2);
                    dirs.resize(2, 0);
                    dirs [ 1 ] = 2; //gamma xy
                    vm [ 1 ] = nodes->giveNode(initalNodeNum);
                }
                dirs [ 0 ] = 0; //eps x
                vm [ 0 ] = nodes->giveNode(initalNodeNum);
                mults [ 0 ] = diff.x();
                mults [ 1 ] = diff.y() / 2;
                jd = new JointDoF(s, dirs [ 0 ], vm, dirs, mults);
                constrs->addConstraint(jd);

                //direction Y  (gammaxz realized here)
                if ( dim == 3 ) {
                    dirs [ 1 ] = 5; //gamma xy
                    dirs [ 2 ] = 3; //gamma yz
                    mults [ 2 ] = diff.z() / 2;
                }
                dirs [ 0 ] = 1; //eps y
                mults [ 0 ] = diff.y();
                mults [ 1 ] = diff.x() / 2;
                jd = new JointDoF(s, dirs [ 0 ], vm, dirs, mults);
                constrs->addConstraint(jd);

                //direction Z  (gammaxz realized here)
                if ( dim == 3 ) {
                    dirs [ 1 ] = 4; //gamma xz
                    dirs [ 2 ] = 3; //gamma yz
                    mults [ 2 ] = diff.y() / 2;

                    dirs [ 0 ] = 2; //eps z
                    mults [ 0 ] = diff.z();
                    mults [ 1 ] = diff.x() / 2;
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
                    mults [ 2 ] = diff.z();
                } else if ( dim == 2 ) {
                    vm.resize(2);
                    mults.resize(2);
                    dirs.resize(2, 0);
                    dirs [ 1 ] = 2; //gamma xy
                    vm [ 1 ] = nodes->giveNode(initalNodeNum);
                }
                dirs [ 0 ] = 0; //eps x
                vm [ 0 ] = nodes->giveNode(initalNodeNum);
                mults [ 0 ] = diff.x();
                mults [ 1 ] = diff.y();
                jd = new JointDoF(s, dirs [ 0 ], vm, dirs, mults);
                constrs->addConstraint(jd);

                //direction Y  (gammaxz realized here)
                if ( dim == 3 ) {
                    vm.resize(2);
                    mults.resize(2);
                    dirs.resize(2, 0);
                    dirs [ 1 ] = 3; //gamma yz
                    vm [ 1 ] = nodes->giveNode(initalNodeNum);
                    mults [ 1 ] = diff.z();
                } else if ( dim == 2 ) {
                    vm.resize(1);
                    mults.resize(1);
                    dirs.resize(1, 0);
                }
                dirs [ 0 ] = 1; //eps y
                vm [ 0 ] = nodes->giveNode(initalNodeNum);
                mults [ 0 ] = diff.y();
                jd = new JointDoF(s, dirs [ 0 ], vm, dirs, mults);
                constrs->addConstraint(jd);

                //direction Z  (gammaxz realized here)
                if ( dim == 3 ) {
                    vm.resize(1);
                    mults.resize(1);
                    dirs.resize(1, 0);
                    dirs [ 0 ] = 2;  //eps z
                    vm [ 0 ] = nodes->giveNode(initalNodeNum);
                    mults [ 0 ] = diff.z();
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
        mn = new MechDoF(dim, rotnum);
        nodes->addNode(mn);
        //it is better to leave rotations free, setting them is not clear in case of shear stress loading
        /*
         * BoundaryCondition *bc;
         * vector< int >dBC, nBC;
         * vector< double > bcmults;
         * dBC.resize(rotnum);
         * nBC.resize(rotnum);
         * bcmults.resize(rotnum);
         * for ( unsigned i = 0; i < rotnum; i++ ) {
         *  dBC[i] = strainFunc[dim+i];
         *  nBC[i] = stressFunc[dim+i];
         *  if(nBC[i]<0) bcmults[i] = 1;
         *  else bcmults[i] = volume;
         * }
         * bc = new BoundaryCondition(nodes->giveNode(rotDoF), dBC, nBC, bcmults);
         * //bcs->addBoundaryCondition(bc);
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
                    mults [ 0 ] = 0.5 * pow(-1., i + 1);
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
void MechanicalPeriodicBCwithElasticConstraint :: apply(NodeContainer *nodes, ElementContainer *elems, BCContainer *bcs, ConstraintContainer *constrs, FunctionContainer *funcs, ExporterContainer *ex, MaterialContainer *mats, RegionContainer *regions, Solver *solver) {
    ( void ) elems;
    ( void ) mats;
    ( void ) solver;
    ( void ) regions;

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
    generateRigidBodyBC(nodes, elems, bcs, constrs, funcs);

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
    linS->setContainers(masterModel->giveElements(), masterModel->giveNodes(), masterModel->giveFunctions(),  masterModel->giveBC() );
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
        elastSol [ i ] = linS->giveTrialDoFValues() / ( linS->giveTime() * 1e-6 );
        dBC [ i ] = cfunc;
    }
    delete linS;
    linS = nullptr;

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
            mults [ 3 ] = diff.z();
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
        mults [ 1 ] = diff.x();
        mults [ 2 ] = diff.y();
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
        dBC.resize(m->giveNumberOfDoFs(), funcs->giveSize() );       //todo: conversion from 'size_t' to 'const _Ty', possible loss of data
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
            TrsDoF *tn = new TrsDoF(dim, 1);
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
// Voigt's constraint
//////////////////////////////////////////////////////////
VoigtConstraint :: VoigtConstraint() {
    volume = 1;
}

//////////////////////////////////////////////////////////
VoigtConstraint :: ~VoigtConstraint() {}

//////////////////////////////////////////////////////////
void VoigtConstraint :: apply(NodeContainer *nodes, ElementContainer *elems, BCContainer *bcs, ConstraintContainer *constrs, FunctionContainer *funcs, ExporterContainer *ex, MaterialContainer *mats, RegionContainer *regions, Solver *solver) {
    ( void ) elems;
    ( void ) mats;
    ( void ) solver;
    ( void ) regions;

    MechDoF *master;
    unsigned masterNodeNum = nodes->giveSize();
    master = new MechDoF(dim, 3 * ( dim - 1 ) );
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
                mults [ 2 ] = diff.z() / 2;
            } else if ( dim == 2 ) {
                vm.resize(2);
                mults.resize(2);
                dirs.resize(2, 0);
                dirs [ 1 ] = 2; //gamma xy
                vm [ 1 ] = nodes->giveNode(masterNodeNum);
            }
            dirs [ 0 ] = 0; //eps x
            vm [ 0 ] = nodes->giveNode(masterNodeNum);
            mults [ 0 ] = diff.x();
            mults [ 1 ] = diff.y() / 2;
            jd = new JointDoF(s, dirs [ 0 ], vm, dirs, mults);
            constrs->addConstraint(jd);

            //direction Y  (gammaxz realized here)
            if ( dim == 3 ) {
                dirs [ 1 ] = 5; //gamma xy
                dirs [ 2 ] = 3; //gamma yz
                mults [ 2 ] = diff.z() / 2;
            }
            dirs [ 0 ] = 1; //eps y
            mults [ 0 ] = diff.y();
            mults [ 1 ] = diff.x() / 2;
            jd = new JointDoF(s, dirs [ 0 ], vm, dirs, mults);
            constrs->addConstraint(jd);

            //direction Z  (gammaxz realized here)
            if ( dim == 3 ) {
                dirs [ 1 ] = 4; //gamma xz
                dirs [ 2 ] = 3; //gamma yz
                mults [ 2 ] = diff.y() / 2;

                dirs [ 0 ] = 2; //eps z
                mults [ 0 ] = diff.z();
                mults [ 1 ] = diff.x() / 2;
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

    while (  iss >> param ) {
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
void PressureFromMechanicalLoad :: apply(NodeContainer *nodes, ElementContainer *elems, BCContainer *bcs, ConstraintContainer *constrs, FunctionContainer *funcs, ExporterContainer *ex, MaterialContainer *mats, RegionContainer *regions, Solver *solver) {
    ( void ) elems;
    ( void ) bcs;
    ( void ) funcs;
    ( void ) ex;
    ( void ) solver;
    ( void ) regions;

    VectTrsprtCoupledMaterial *dtcm = dynamic_cast< VectTrsprtCoupledMaterial * >( mats->giveMaterial(materialnum) );
    if ( dtcm == nullptr ) {
        cerr << "Error in PressureFromMechanicalLoad: material is not of VectTrsprtCoupledMaterial type" << endl;
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
    masters [ 0 ] = nodes->giveNode(master);
    vector< double >mults(1);
    mults [ 0 ] = multiplier / biot;

    DoFDependentOnConjugates *ddc;
    Node *trs;
    for ( auto &p: trsprtnodes ) {
        trs = nodes->giveNode(p);
        if ( trs->giveRelativeDoFPhysicalFieldNum(direction) != 1 ) { //physical field 1 is transport
            cerr << "Error in PressureFromMechanicalLoad: node " << p << " direction " << direction <<  " is not pressure DoF" << endl;
            exit(1);
        }
        ddc = new DoFDependentOnConjugates(trs, direction, masters, directions, mults);
        constrs->addConstraint(ddc);
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
void RigidPlate :: readFromLine(istringstream &iss, unsigned d) {
    this->dim = d;
    unsigned nslaves, nodeid, num;
    bool b;
    string param;
    iss >> master_id;

    while ( iss >> param ) {
        if ( param.compare("dependent_nodes") == 0 ) {
            iss >> nslaves;
            for ( unsigned i = 0; i < nslaves; i++ ) {
                iss >> nodeid;
                slave_ids.push_back(nodeid);
            }
        } else if ( param.compare("inside") == 0 ) {
            iss >> num;
            insideRegions.resize(num);
            for ( unsigned i = 0; i < num; i++ ) {
                iss >> insideRegions [ i ];
            }
        } else if ( param.compare("outside") == 0 ) {
            iss >> num;
            outsideRegions.resize(num);
            for ( unsigned i = 0; i < num; i++ ) {
                iss >> outsideRegions [ i ];
            }
        } else if ( param.compare("dirs") == 0 ) {
            iss >> num;
            activeDirs.resize(num);
            for ( unsigned i = 0; i < num; i++ ) {
                iss >> b;
                activeDirs [ i ] = b;
            }
        }
    }
}

//////////////////////////////////////////////////////////
void RigidPlate :: checkPhysicalField(Node *master) {
    if ( activeDirs.size() == 0 ) {
        activeDirs.resize(master->giveNumberOfDoFs(), true);
    }
}

//////////////////////////////////////////////////////////
void RigidPlate :: apply(NodeContainer *nodes, ElementContainer *elems, BCContainer *bcs, ConstraintContainer *constrs, FunctionContainer *funcs, ExporterContainer *ex, MaterialContainer *mats, RegionContainer *regions, Solver *solver) {
    ( void ) elems;
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
    this->checkPhysicalField(master);

    for ( auto const &sl_id : slave_ids ) {
        slave = nodes->giveNode(sl_id);
        connectSlaveMasterRigid(constrs, slave, master, this->dim, activeDirs);
    }

    if ( insideRegions.size() > 0 ) {
        for ( auto const &nod : * nodes ) {
            //test for regions
            if ( !regions->isLocationValid(nod->givePoint(), insideRegions, outsideRegions) ) {
                continue;
            }
            if ( nod == master || endsWith(nod->giveName(), "virtual") ) {
                continue;
            }
            if ( nod->giveNumberOfDoFs() == 0 ) {
                continue;
            }
            if ( dynamic_cast< FreeDoF * >( nod ) ) {
                continue;
            }
            connectSlaveMasterRigid(constrs, nod, master, this->dim, activeDirs);
        }
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

    bool b;
    unsigned num;
    string param;
    while (  iss >> param ) {
        if ( param.compare("dirs") == 0 ) {
            iss >> num;
            activeDirs.resize(num);
            for ( unsigned i = 0; i < num; i++ ) {
                iss >> b;
                activeDirs [ i ] = b;
            }
        }
    }
}

//////////////////////////////////////////////////////////
void CoordRigidPlate :: apply(NodeContainer *nodes, ElementContainer *elems, BCContainer *bcs, ConstraintContainer *constrs, FunctionContainer *funcs, ExporterContainer *ex, MaterialContainer *mats, RegionContainer *regions, Solver *solver) {
    ( void ) elems;
    ( void ) bcs;
    ( void ) funcs;
    ( void ) ex;
    ( void ) mats;
    ( void ) regions;
    ( void ) solver;

    // jointDoF jD;
    Node *master;

    master = nodes->giveNode(master_id);
    // check if it is master node
    RigidPlate :: checkPhysicalField(master);

    for ( auto const &nod : * nodes ) {
        if ( isInBlock(nod->givePoint(), leftBottom, rightTop) ) {
            if ( nod == master || endsWith(nod->giveName(), "virtual") ) {
                continue;
            }
            if ( nod->giveNumberOfDoFs() == 0 ) {
                continue;
            }
            if ( dynamic_cast< FreeDoF * >( nod ) ) {
                continue;
            }
            connectSlaveMasterRigid(constrs, nod, master, this->dim, activeDirs);
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
    iss >> rI >> rO >> this->w0 >> this->w1;
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

    bool b;
    unsigned num;
    string param;
    while (  iss >> param ) {
        if ( param.compare("dirs") == 0 ) {
            iss >> num;
            activeDirs.resize(num);
            for ( unsigned i = 0; i < num; i++ ) {
                iss >> b;
                activeDirs [ i ] = b;
            }
        }
    }
}

//////////////////////////////////////////////////////////
void RingRigidPlate :: apply(NodeContainer *nodes, ElementContainer *elems, BCContainer *bcs, ConstraintContainer *constrs, FunctionContainer *funcs, ExporterContainer *ex, MaterialContainer *mats, RegionContainer *regions, Solver *solver) {
    ( void ) elems;
    ( void ) bcs;
    ( void ) funcs;
    ( void ) ex;
    ( void ) mats;
    ( void ) regions;
    ( void ) solver;

    // jointDoF jD;
    Node *master;

    master = nodes->giveNode(master_id);
    // check if it is master node
    RigidPlate :: checkPhysicalField(master);

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
    this->center = Point(this->center.x() * xm, this->center.y() * ym, this->center.z() * zm);
    for ( auto const &nod : * nodes ) {
        node_point = Point(nod->givePoint().x() * xm, nod->givePoint().y() * ym, nod->givePoint().z() * zm);
        if ( isInCircle(node_point, this->center, this->r_outer, this->direction) ) {
            if ( !isInCircle(node_point, this->center, this->r_inner, this->direction) ) {
                if ( node_point(direction) < w0 || node_point(direction) > w1 ) {
                    continue;
                }
                if ( nod == master ) {
                    continue;
                }
                connectSlaveMasterRigid(constrs, nod, master, this->dim, activeDirs);
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
    while (  iss >> param ) {
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
void ExpansionRing :: apply(NodeContainer *nodes, ElementContainer *elems, BCContainer *bcs, ConstraintContainer *constrs, FunctionContainer *funcs, ExporterContainer *ex, MaterialContainer *mats, RegionContainer *regions, Solver *solver) {
    ( void ) elems;
    ( void ) bcs;
    ( void ) funcs;
    ( void ) ex;
    ( void ) mats;
    ( void ) regions;
    ( void ) solver;

    // jointDoF jD;
    Node *master;

    master = nodes->giveNode(master_id);
    // check if it is master node
    RigidPlate :: checkPhysicalField(master);

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
    this->center = Point(this->center.x() * xm, this->center.y() * ym, this->center.z() * zm);
    for ( auto const &nod : * nodes ) {
        node_point = Point(nod->givePoint().x() * xm, nod->givePoint().y() * ym, nod->givePoint().z() * zm);
        if ( isInCircle(node_point, this->center, this->r_outer, this->direction) ) {
            if ( !isInCircle(node_point, this->center, this->r_inner, this->direction) ) {
                if ( nod == master || !dynamic_cast< MechNode * >( nod ) || dynamic_cast< MechDoF * >( nod ) ) {
                    continue;
                }
                connectSlaveMasterExpansion(constrs, nod, master, this->dim, this->transport, funcs->giveFunction(this->fn_id) );
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
    while (  iss >> param ) {
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

void ExpansionRingDoFLoad :: apply(NodeContainer *nodes, ElementContainer *elems, BCContainer *bcs, ConstraintContainer *constrs, FunctionContainer *funcs, ExporterContainer *ex, MaterialContainer *mats, RegionContainer *regions, Solver *solver) {
    ( void ) elems;
    ( void ) bcs;
    ( void ) funcs;
    ( void ) ex;
    ( void ) mats;
    ( void ) regions;
    ( void ) solver;

    // jointDoF jD;
    Node *master;
    Node *expMaster;

    master = nodes->giveNode(this->master_id);
    RigidPlate :: checkPhysicalField(master);

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
    this->center = Point(this->center.x() * xm, this->center.y() * ym, this->center.z() * zm);
    for ( auto const &nod : * nodes ) {
        node_point = Point(nod->givePoint().x() * xm, nod->givePoint().y() * ym, nod->givePoint().z() * zm);
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


void ExpansionRingSingleDoFLoad :: apply(NodeContainer *nodes, ElementContainer *elems, BCContainer *bcs, ConstraintContainer *constrs, FunctionContainer *funcs, ExporterContainer *ex, MaterialContainer *mats, RegionContainer *regions, Solver *solver) {
    ( void ) elems;
    ( void ) bcs;
    ( void ) funcs;
    ( void ) ex;
    ( void ) mats;
    ( void ) regions;
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

    MechNode *mn;
    MechDoF *md;
    this->center = Point(this->center.x() * xm, this->center.y() * ym, this->center.z() * zm);

    for ( auto const &nod : * nodes ) {
        mn = dynamic_cast< MechNode * >( nod );
        if ( mn == nullptr ) {
            continue;
        }
        md = dynamic_cast< MechDoF * >( nod );
        if ( md != nullptr ) {
            continue;
        }
        if ( endsWith(nod->giveName(), "virtual") ) {
            continue;
        }

        node_point = Point(nod->givePoint().x() * xm, nod->givePoint().y() * ym, nod->givePoint().z() * zm);
        if ( isInCircle(node_point, this->center, this->r_outer, this->direction) ) {
            if ( !isInCircle(node_point, this->center, this->r_inner, this->direction) ) {
                // if ( nod == expMaster ) {
                //     continue;
                // }
                num_nodes++;
                // std::cout << "nod name = " << nod->giveName() << '\n';

                // n_i = ( Point((this->direction == 0) ? 0 : nod->givePoint().x(),
                //               (this->direction == 1) ? 0 : nod->givePoint().y(),
                //               (this->direction == 2) ? 0 : nod->givePoint().z()) -
                //       Point((this->direction == 0) ? 0 : this->center.x(),
                //                     (this->direction == 1) ? 0 : this->center.y(),
                //                     (this->direction == 2) ? 0 : this->center.z())
                //         );
                n_i = ( node_point - this->center );
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
                    if ( * std :: max_element(n_vect.begin(), n_vect.end() ) >  abs( * std :: min_element(n_vect.begin(), n_vect.end() ) ) ) {
                        slave_dir = std :: distance(n_vect.begin(), std :: max_element(n_vect.begin(), n_vect.end() ) );
                    } else {
                        slave_dir = std :: distance(n_vect.begin(), std :: min_element(n_vect.begin(), n_vect.end() ) );
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
    multipliers.push_back(double( num_nodes ) * ( r_outer + r_inner ) / 4.);
    masterNodes.push_back(expMaster);
    directions.push_back(0);

    if ( !masterNodes.empty() ) {
        for ( unsigned j = 0; j < multipliers.size(); j++ ) {
            multipliers [ j ] /= slave_dir_vect_value;
        }
        JointDoF *newJD = new JointDoF(slave, slave_dir, masterNodes, directions, multipliers);
        for ( unsigned u = 0; u < multipliers.size(); u++ ) {
            if ( multipliers [ u ] != multipliers [ u ] ) {
                cout << "multipliers ERROR" << endl;
                exit(1);
            }
        }
        constrs->addConstraint(newJD);
    }

    masterNodes.clear();
    multipliers.clear();
    directions.clear();
    std :: cout << "Expansion volumetric load applied: center(" << this->center.x() << ", " << this->center.y() << ", " << this->center.z() << ") rI = " << this->r_inner << ", rO = " << this->r_outer << ", direction: " << this->direction <<  '\n';
    // exit(0);
}

/////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////
// NormalSurfaceLoad
//////////////////////////////////////////////////////////
void NormalSurfaceLoad :: readFromLine(istringstream &iss, unsigned d) {
    dim = d;
    bool bf = false;
    string param;
    unsigned num;
    while (  iss >> param ) {
        if ( param.compare("load_function") == 0 ) {
            iss >> this->fnID;
            bf = true;
        } else if ( param.compare("inside") == 0 ) {
            iss >> num;
            insideRegions.resize(num);
            for ( unsigned i = 0; i < num; i++ ) {
                iss >> insideRegions [ i ];
            }
        } else if ( param.compare("outside") == 0 ) {
            iss >> num;
            outsideRegions.resize(num);
            for ( unsigned i = 0; i < num; i++ ) {
                iss >> outsideRegions [ i ];
            }
        }
    }



    if ( !bf ) {
        std :: cerr << "Error: no function governing surface load specified" << '\n';
        exit(EXIT_FAILURE);
    }
}

//////////////////////////////////////////////////////////
void NormalSurfaceLoad :: apply(NodeContainer *n, ElementContainer *e, BCContainer *b, ConstraintContainer *c, FunctionContainer *funcs, ExporterContainer *ex, MaterialContainer *mats, RegionContainer *regions, Solver *solver) {
    //add boundary conditions
    ( void ) n;
    ( void ) c;
    ( void ) funcs;
    ( void ) ex;
    ( void ) mats;
    ( void ) solver;
    RigidBodyBoundary *rbb;
    LDPMTetra *tet;
    BoundaryCondition *bc;
    vector< int >dBC, nBC;
    Point diff;
    dBC.resize( ( dim - 1 ) * 3, -1 );
    nBC.resize( ( dim - 1 ) * 3, fnID );
    vector< double >mult;
    mult.resize( ( dim - 1 ) * 3, 0. );
    Vector normal;
    double area;
    for ( auto &el : * e ) {
        rbb = dynamic_cast< RigidBodyBoundary * >( el );
        if ( rbb ) {
            //test for regions
            if ( !regions->isLocationValid(rbb->giveIPLoc(0), insideRegions, outsideRegions) ) {
                continue;
            }
            area = rbb->giveArea();
            normal = rbb->giveNormal();
            diff = rbb->giveIPLoc(0) - rbb->giveNode(0)->givePoint();
            //forces
            for ( unsigned i = 0; i < dim; i++ ) {
                mult [ i ] = normal [ i ] * area;
            }
            //moments
            if ( dim == 2 ) {
                mult [ 2 ] = mult [ 1 ] * diff [ 0 ] - mult [ 0 ] * diff [ 1 ];
            } else if ( dim == 3 ) {
                mult [ 3 ] = mult [ 2 ] * diff [ 1 ] - mult [ 1 ] * diff [ 2 ];
                mult [ 4 ] = mult [ 0 ] * diff [ 2 ] - mult [ 2 ] * diff [ 0 ];
                mult [ 5 ] = mult [ 1 ] * diff [ 0 ] - mult [ 0 ] * diff [ 1 ];
            }

            bc = new BoundaryCondition(rbb->giveNode(0), dBC, nBC, mult);
            b->addBoundaryCondition(bc);
        }
    }
    for ( unsigned i = 0; i < mult.size(); i++ ) {
        mult [ i ] = 0;
    }
    for ( auto &el : * e ) {
        tet = dynamic_cast< LDPMTetra * >( el );
        if ( tet ) {
            vector< Node * >innodes;
            vector< Node * >tetnodes = tet->giveNodes();
            for ( auto &nn:tetnodes ) {
                if ( regions->isLocationValid(nn->givePoint(), insideRegions, outsideRegions) ) {
                    innodes.push_back(nn);
                }
            }
            if ( innodes.size() != 3 ) {
                continue;
            }

            area = triArea3D( innodes [ 0 ]->givePointPointer(), innodes [ 1 ]->givePointPointer(), innodes [ 2 ]->givePointPointer() );
            normal = ( innodes [ 1 ]->givePoint() - innodes [ 0 ]->givePoint() ).cross( innodes [ 2 ]->givePoint() - innodes [ 0 ]->givePoint() );
            normal /= normal.norm();
            if ( ( innodes [ 0 ]->givePoint() - tet->giveVertex(0)->givePoint() ).dot(normal) < 0. ) {
                normal *= -1; //must be outward normal
            }
            //forces
            for ( unsigned i = 0; i < dim; i++ ) {
                mult [ i ] = normal [ i ] * area / 3.;
            }

            for ( unsigned k = 0; k < 3; k++ ) {
                bc = new BoundaryCondition(innodes [ k ], dBC, nBC, mult);
                b->addBoundaryCondition(bc);
            }
        }
    }
}

/////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////
// Hanging Node
//////////////////////////////////////////////////////////
void MechHangingNode :: readFromLine(istringstream &iss, unsigned d) {
    ( void ) d;
    if ( iss >> nodeid ) {} else {
        std :: cerr << "Error in MechHangingNode: nodeid not found" << '\n';
        exit(EXIT_FAILURE);
    }
    if ( iss >> elemid ) {} else {
        std :: cerr << "Error in MechHangingNode: elemid not found" << '\n';
        exit(EXIT_FAILURE);
    }
}

//////////////////////////////////////////////////////////
void MechHangingNode :: apply(NodeContainer *n, ElementContainer *e, BCContainer *b, ConstraintContainer *c, FunctionContainer *funcs, ExporterContainer *ex, MaterialContainer *mats, RegionContainer *regions, Solver *solver) {
    //add boundary conditions
    ( void ) funcs;
    ( void ) ex;
    ( void ) mats;
    ( void ) solver;
    ( void ) b;
    ( void ) regions;

    Element *ee = e->giveElement(elemid);
    unsigned ndim = ee->giveDimension();
    Point natcoords = ee->findNaturalCoords( n->giveNode(nodeid)->givePointPointer() );
    Vector weights = ee->giveShapeFunctions(& natcoords);
    vector< Node * >masters = ee->giveNodes();
    vector< unsigned >dirs( masters.size() );
    vector< double >mults( masters.size() );

    for ( unsigned k = 0; k < masters.size(); k++ ) {
        mults [ k ] = weights [ k ];
    }
    for ( unsigned i = 0; i < ndim; i++ ) {
        for ( unsigned k = 0; k < masters.size(); k++ ) {
            dirs [ k ] = i;
        }
        JointDoF *newJD = new JointDoF(n->giveNode(nodeid), i, masters, dirs, mults);
        c->addConstraint(newJD);
    }
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
void PBlockContainer :: setContainers(NodeContainer *n, ElementContainer *e, BCContainer *b, ConstraintContainer *c, FunctionContainer *f, ExporterContainer *ex, MaterialContainer *mats, RegionContainer *r, Solver *solv) {
    nodes = n;
    elems = e;
    bcs = b;
    constrs = c;
    funcs = f;
    exporters = ex;
    materials = mats;
    regions = r;
    solver = solv;
}

//////////////////////////////////////////////////////////
void PBlockContainer :: init() {
    for ( auto b: blocks ) {
        b->apply(nodes, elems, bcs, constrs, funcs, exporters, materials, regions, solver);
    }
}

//////////////////////////////////////////////////////////
void PBlockContainer :: readFromFile(const string filename, unsigned dim) {
    unsigned origsize = blocks.size(); //todo: warning C4267: 'initializing': conversion from 'size_t' to 'unsigned int', possible loss of dat
    string line, ftype;
    ifstream inputfile(filename.c_str() );
    if ( inputfile.is_open() ) {
        while ( getline(inputfile >> std :: ws, line) ) {
            if ( line.empty() || ( line.at(0) == '#' ) ) {
                continue;
            }
            istringstream iss(line);
            iss >> ftype;
            if ( !( ftype.rfind("#", 0) == 0 ) ) {
                if ( ftype.compare("MechanicalPeriodicBC") == 0 ) {
                    MechanicalPeriodicBC *newblock = new MechanicalPeriodicBC();
                    newblock->readFromLine(iss, dim);
                    blocks.push_back(newblock);
                } else if ( ftype.compare("MechanicalSphericalPeriodicBCExperimental") == 0 ) {
                    MechanicalSphericalPeriodicBCExperimental *newblock = new MechanicalSphericalPeriodicBCExperimental();
                    newblock->readFromLine(iss, dim);
                    blocks.push_back(newblock);
                } else if ( ftype.compare("MechanicalSphericalPeriodicBC") == 0 ) {
                    MechanicalSphericalPeriodicBC *newblock = new MechanicalSphericalPeriodicBC();
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
                } else if ( ftype.compare("NormalSurfaceLoad") == 0 ) {
                    NormalSurfaceLoad *newblock = new NormalSurfaceLoad();
                    newblock->readFromLine(iss, dim);
                    blocks.push_back(newblock);
                } else if ( ftype.compare("MechHangingNode") == 0 ) {
                    MechHangingNode *newblock = new MechHangingNode();
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



void connectSlaveMasterRigid(ConstraintContainer *constrs, Node *slave, Node *master, unsigned const &ndim, const vector< bool > &activeDirs, bool includeRigidBodyRotation) {
    vector< unsigned >slaveDoFNuM = slave->givePhysicalFieldsDoFNum();
    vector< unsigned >masterDoFNuM = master->givePhysicalFieldsDoFNum();
    vector< unsigned >physInd;  // Mechanics, Transport, Thermal, Humidity
    physInd.resize(4, 0);

    // calculate multipliers and construct jointDof for every slaveDof
    vector< Node * >masterNodes;
    vector< double >multipliers;
    vector< unsigned >directions;


    unsigned kmaster = 0;
    unsigned kslave = 0;
    unsigned kslavelast = 0;
    unsigned stop = 0;
    for ( unsigned i = 0; i < masterDoFNuM.size(); i++ ) {
        kslave = kslavelast;
        kmaster = stop;
        stop += masterDoFNuM [ i ];
        kslavelast += slaveDoFNuM [ i ];
        for ( ; kmaster < stop && kslave < kslavelast; kmaster++, kslave++ ) {
            if ( activeDirs [ kmaster ] ) {
                masterNodes.push_back(master);
                directions.push_back(kmaster);
                multipliers.push_back(1);
                //rigid body rotation
                if ( i == 0 && includeRigidBodyRotation && masterDoFNuM [ 0 ] == 3 * ( ndim - 1 ) && kslave < ndim ) {
                    if ( ndim == 2 ) {
                        masterNodes.push_back(master);
                        directions.push_back(2);
                        if ( kslave == 0 ) {
                            multipliers.push_back(-( slave->givePoint().y() - master->givePoint().y() ) );
                        } else {
                            multipliers.push_back( ( slave->givePoint().x() - master->givePoint().x() ) );
                        }
                    } else if ( ndim == 3 ) {
                        masterNodes.push_back(master);
                        masterNodes.push_back(master);
                        if ( kslave == 0 ) {
                            directions.push_back(4);
                            multipliers.push_back( ( slave->givePoint().z() - master->givePoint().z() ) );
                            directions.push_back(5);
                            multipliers.push_back(-( slave->givePoint().y() - master->givePoint().y() ) );
                        } else if ( kslave == 1 ) {
                            directions.push_back(3);
                            multipliers.push_back(-( slave->givePoint().z() - master->givePoint().z() ) );
                            directions.push_back(5);
                            multipliers.push_back( ( slave->givePoint().x() - master->givePoint().x() ) );
                        } else {
                            directions.push_back(3);
                            multipliers.push_back( ( slave->givePoint().y() - master->givePoint().y() ) );
                            directions.push_back(4);
                            multipliers.push_back(-( slave->givePoint().x() - master->givePoint().x() ) );
                        }
                    }
                }
                JointDoF *newJD = new JointDoF(slave, kslave, masterNodes, directions, multipliers);
                constrs->addConstraint(newJD);
                masterNodes.clear();
                multipliers.clear();
                directions.clear();
            }
        }
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
    if ( * std :: max_element(n_vect.begin(), n_vect.end() ) >  abs( * std :: min_element(n_vect.begin(), n_vect.end() ) ) ) {
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
    if ( * std :: max_element(n_vect.begin(), n_vect.end() ) >  abs( * std :: min_element(n_vect.begin(), n_vect.end() ) ) ) {
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
