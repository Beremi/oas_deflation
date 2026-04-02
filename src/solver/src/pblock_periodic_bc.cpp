#include "pblock_periodic_bc.h"
#include "model.h"
#include "material_vectorial.h"
#include "element_ldpm.h"
#include "solver_implicit.h"

using namespace std;


//////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////
// Mechanical Periodic BC
//////////////////////////////////////////////////////////

MechanicalPeriodicBC::~MechanicalPeriodicBC(){
    //if ( crackOpeningDoF != nullptr ) {
    //        delete crackOpeningDoF;
    //}
}

//////////////////////////////////////////////////////////
void MechanicalPeriodicBC :: generateNewDoFs(NodeContainer *nodes) {
    //create new degrees of freedom representing strains ex, ey, gammaxy=2exy or ex, ey, ez, gammyz, gammaxz, gammaxy,
    MechDoF *mn;
    initalNodeNum = nodes->giveSize();
    mn = new MechDoF( dim, 3 * ( dim - 1 ) );
    nodes->addNode(mn);
}

//////////////////////////////////////////////////////////
void MechanicalPeriodicBC :: setEigenStrain(Vector &eigstr) {
    if ( strainFunc.size() != ( unsigned ) eigstr.size() ) {
        cout << giveName() << " uses " << strainFunc.size() << " strain values but " << eigstr.size() << " were submitted" << endl;
        exit(1);
    }
    strainBC->setInitialDoFFields(eigstr);
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
        Node *m = constrs->giveConstraint(constrs->giveConstraintsSize() - 1)->giveMasterNode(0);// warning C4267: 'argument': conversion from 'size_t' to 'const unsigned int', possible loss of data
        BoundaryCondition *bc;
        vector< int >dBC, nBC;
        dBC.resize(m->giveNumberOfDoFs(), -1);         //todo: warning C4267: 'argument': conversion from 'size_t' to 'const _Ty', possible loss of data
        for ( unsigned k = 0; k < dim; k++ ) {
            dBC [ k ] = funcs->giveSize();
        }
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
            if ( nodes->giveNode(n)->doesMechanics() &&  dynamic_cast< MechDoF * >( nodes->giveNode(n) ) ) {
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
void MechanicalPeriodicBC :: apply(Model *model) {
    NodeContainer *nodes = model->giveNodes();
    ElementContainer *elems = model->giveElements();
    BCContainer *bcs = model->giveBoundaryConditions();
    ConstraintContainer *constrs = model->giveConstraints();
    FunctionContainer *funcs = model->giveFunctions();
    ExporterContainer *ex = model->giveExporters();

    calculateVolume(elems);

    unsigned const_num = constrs->giveConstraintsSize();
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
            cerr << "Error in Periodic boundary condition: cannot prescribe both stress and strain for the same direction, direction " << i << endl;
        }
    }
    strainBC = new BoundaryCondition(nodes->giveNode(initalNodeNum), dBC, nBC, bcmults);
    bcs->addBoundaryCondition(strainBC);

    //export data
    generateExporters(nodes, ex);

    cout << "Applied periodic boundary conditions: " << nodes->giveSize() - initalNodeNum << " new DoFs (nodes " << initalNodeNum << " - " <<  nodes->giveSize() - 1 << "); " << constrs->giveConstraintsSize() - const_num << " new constraints; " << bcs->giveSize() - bcs_num << " new boundary conditions; " << funcs->giveSize() - funcs_num << " new function; " << ex->giveSize() - ex_num << " new exporters; " << "created" << endl;
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
    unsigned num, init_num;
    double gamma = volume/PUCsize[0]; //PUCsize is diameter    
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
                constrs->addConstraint(jd); //rotations
            }
        }

        //connect translations
        diff = s->givePoint() - m->givePoint();
        if ( !nonsymmetric_shear ) {
            //direction X  (all gammaxy and gammaxy realized here)
            if ( dim == 3 ) {
                init_num = 4; 
                num = init_num;
                if (crack_active) num +=3;
                vm.resize(num);
                mults.resize(num);
                dirs.resize(num, 0);
                dirs [ 2 ] = 5; //gamma xy
                dirs [ 3 ] = 4; //gamma xz
                vm [ 2 ] = nodes->giveNode(initalNodeNum);
                vm [ 3 ] = nodes->giveNode(initalNodeNum); //gamma xz
                mults [ 3 ] = diff.z() / 2;
            } else if ( dim == 2 ) {
                init_num = 3; 
                num = init_num;
                if (crack_active) num +=2;         
                vm.resize(num);
                mults.resize(num);
                dirs.resize(num, 0);
                dirs [ 2 ] = 2; //gamma xy
                vm [ 2 ] = nodes->giveNode(initalNodeNum);
            }
            if (crack_active){
                for (unsigned d=0; d<dim; d++){
                    vm[init_num+d] = crackOpeningDoF;
                    dirs[init_num+d] = d;                    
                }
            }                
            
            dirs [ 0 ] = 0;
            dirs [ 1 ] = 0; //eps x
            vm [ 0 ] = m; //master
            vm [ 1 ] = nodes->giveNode(initalNodeNum);
            mults [ 0 ] = 1;
            mults [ 1 ] = diff.x();
            mults [ 2 ] = diff.y() / 2;
            if (crack_active){
                mults[3] = - (diff.x()*crack_normal.x()+diff.y()*crack_normal.y()/4.)/gamma + copysign(1.,diff.x()*crack_normal.x()+diff.y()*crack_normal.y()); 
                mults[4] = - (diff.y()*crack_normal.x()/4.)/gamma;                
            }            
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
            if (crack_active){
                mults[3] = - (diff.x()*crack_normal.y()/4.)/gamma;             
                mults[4] = - (diff.x()*crack_normal.x()/4.+diff.y()*crack_normal.y())/gamma + copysign(1.,diff.x()*crack_normal.x()+diff.y()*crack_normal.y());     
            }                    
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
// MechanicalSphericalPeriodicBCWithCrack
//////////////////////////////////////////////////////////

void MechanicalSphericalPeriodicBCWithCrack:: readFromLine(std :: istringstream &iss, unsigned d){
    MechanicalSphericalPeriodicBC :: readFromLine(iss, d);

    iss.clear(); // clear string stream
    iss.seekg(0, iss.beg); //reset position in string stream
    
    bool bcrack_normal = false;
    bool bcrack_DoF = false;
    
    unsigned num;
    string param;
    while (  iss >> param ) {
        if ( param.compare("crack_normal") == 0 ) {
            iss >> num;
            if (num!=dim){
                cerr << "MechanicalSphericalPeriodicBCWithCrack Error: crack_normal size is " << num << " but dimension of the problem is " << dim << endl;
                exit(1);
            }
            crack_normal.resize(num);
            for ( unsigned i = 0; i < num; i++ ) {
                iss >> crack_normal [ i ];
            }
            bcrack_normal = true;
        } else if ( param.compare("crack_DoF_node") == 0 ) {
            iss >> crackDoFNodeNumber;
            bcrack_DoF = true;
        }         
    }
    if (!bcrack_normal || !bcrack_DoF){
        cout << "Crack nor activated, missing either crack_normal or crack_DoF_node" << endl;   
    } else {
        crack_active = true;
    }
}    
    
//////////////////////////////////////////////////////////
void MechanicalSphericalPeriodicBCWithCrack :: generateConstraints(NodeContainer *nodes, ConstraintContainer *constrs) {
    MechanicalSphericalPeriodicBC :: generateConstraints(nodes, constrs);

    
    
    /*
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
    */
}

//////////////////////////////////////////////////////////
void MechanicalSphericalPeriodicBCWithCrack :: apply(Model *model) {

    
    crackOpeningDoF = dynamic_cast<MechDoF*>(model->giveNodes()->giveNode(crackDoFNodeNumber));
    if (!crackOpeningDoF){
      cerr << "MechanicalSphericalPeriodicBCWithCrack Error: expected MechDoF node, received " << model->giveNodes()->giveNode(crackDoFNodeNumber)->giveName() << endl;
      exit(1);
    }
    if (crackOpeningDoF->giveNumberOfDoFs()!=dim){
      cerr << "MechanicalSphericalPeriodicBCWithCrack Error: crack MechDoF does not have " << dim << " degrees of freedom but " <<  crackOpeningDoF->giveNumberOfDoFs() << endl;
      exit(1);
    }    
            
    MechanicalSphericalPeriodicBC :: apply(model);
    
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
void MechanicalPeriodicBCwithElasticConstraint :: apply(Model *model) {
    NodeContainer *nodes = model->giveNodes();
    ElementContainer *elems = model->giveElements();
    BCContainer *bcs = model->giveBoundaryConditions();
    ConstraintContainer *constrs = model->giveConstraints();
    FunctionContainer *funcs = model->giveFunctions();
    Solver *solver = model->giveSolver();
    ExporterContainer *ex = model->giveExporters();

    volume = 1;
    for ( auto const a : PUCsize ) {
        volume *= a;
    }

    unsigned const_num = constrs->giveConstraintsSize();
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
    linS->setContainers( masterModel->giveElements(), masterModel->giveNodes(), masterModel->giveFunctions(),  masterModel->giveBoundaryConditions(), masterModel->giveExporters() );
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
    for ( int p = int( constrs->giveConstraintsSize() ) - 1; p >= int( const_num ); p-- ) {
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
            cerr << "Error in Periodic boundary condition: cannot prescribe both stress and strain for the same direction, direction " << i << endl;
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

    cout << "Applied periodic boundary conditions: " << nodes->giveSize() - initalNodeNum << " new DoFs (nodes " << initalNodeNum << " - " <<  nodes->giveSize() - 1 << "); " << constrs->giveConstraintsSize() - const_num << " new constraints; " << bcs->giveSize() - bcs_num << " new boundary conditions; " << funcs->giveSize() - funcs_num << " new function; " << ex->giveSize() - ex_num << " new exporters; " << "created" << endl;
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
        unsigned p = 0;
        while ( constrs->isDependent( constrs->giveConstraint(p)->giveMasterNode(0) ) ) {
            p++;
        }
        Node *m = constrs->giveConstraint(p)->giveMasterNode(0); //todo:  warning C4267: 'argument': conversion from 'size_t' to 'const unsigned int', possible loss of data
        BoundaryCondition *bc;
        vector< int >dBC, nBC;
        dBC.resize( m->giveNumberOfDoFs(), funcs->giveSize() );          //todo: conversion from 'size_t' to 'const _Ty', possible loss of data
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
            if ( nodes->giveNode(n)->doesTransport() && ( dynamic_cast< TrsDoF * >( nodes->giveNode(n) ) != nullptr ) ) {
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
// Cosserat Mechanical Periodic BC
//////////////////////////////////////////////////////////
void CosseratMechanicalPeriodicBC :: generateNewDoFs(NodeContainer *nodes) {
    //create new degrees of freedom representing strains
    MechDoF *mn;
    initalNodeNum = nodes->giveSize();
    mn = new MechDoF(dim, ( dim == 3 )? 18: 6);
    nodes->addNode(mn);
}

//////////////////////////////////////////////////////////
void CosseratMechanicalPeriodicBC :: generateConstraints(NodeContainer *nodes, ConstraintContainer *constrs) {
    //apply contraints, connect periodic images
    JointDoF *jd;
    vector< Node * >vm;
    vector< unsigned >dirs;
    vector< double >mults;
    Node *s = nullptr;
    Node *m = nullptr;
    Point diff, cdiff2, cdiff1;
    Point centroid = Point(PUCsize [ 0 ], PUCsize [ 1 ], ( dim == 3 ) ? PUCsize [ 2 ] : 0) / 2.;
    for ( unsigned i = 0; i < masters.size(); i++ ) {
        m = nodes->giveNode(masters [ i ]);
        s = nodes->giveNode(slaves [ i ]);

        if ( ( not dynamic_cast< Particle * >( s ) ) || ( not dynamic_cast< Particle * >( m ) ) ) {
            cerr << "CosseratMechanicalPeriodicBC should use only Particle type of nodes, it is " << s->giveName() << " and " << m->giveName() << endl;
        }


        diff = s->givePoint() - m->givePoint();
        cdiff1 = m->givePoint() - centroid;
        cdiff2 = s->givePoint() - centroid;
        //connect translations
        if ( dim == 3 ) {
            vm.resize( 10, nodes->giveNode(initalNodeNum) );
            mults.resize(10);
            dirs.resize(10, 0);
            vm [ 0 ] = m;
            mults [ 0 ] = 1.;
            mults [ 1 ] = diff.x();
            mults [ 2 ] = diff.y();
            mults [ 3 ] = diff.z();
        } else {
            vm.resize(5);
            mults.resize(5);
            dirs.resize(5, 0);
            vm [ 0 ] = m;
            vm [ 1 ] = nodes->giveNode(initalNodeNum);
            vm [ 2 ] = nodes->giveNode(initalNodeNum);
            vm [ 3 ] = nodes->giveNode(initalNodeNum);
            vm [ 4 ] = nodes->giveNode(initalNodeNum);
            mults [ 0 ] = 1;
            mults [ 1 ] = diff.x();
            mults [ 2 ] = diff.y();
        }
        //X DIR
        dirs [ 0 ] = 0;
        if ( dim == 3 ) {
            dirs [ 1 ] = 0;
            dirs [ 2 ] = 7;
            dirs [ 3 ] = 5;
            dirs [ 4 ] = 15;
            dirs [ 5 ] = 13;
            dirs [ 6 ] = 11;
            dirs [ 7 ] = 17;
            dirs [ 8 ] = 10;
            dirs [ 9 ] = 12;
            mults [ 4 ] =  ( cdiff1.x() * cdiff1.y() - cdiff2.x() * cdiff2.y() );
            mults [ 5 ] =  ( cdiff1.y() * cdiff1.y() - cdiff2.y() * cdiff2.y() );
            mults [ 6 ] =  ( cdiff1.z() * cdiff1.y() - cdiff2.z() * cdiff2.y() );
            mults [ 7 ] = -( cdiff1.x() * cdiff1.z() - cdiff2.x() * cdiff2.z() );
            mults [ 8 ] = -( cdiff1.z() * cdiff1.y() - cdiff2.z() * cdiff2.y() );
            mults [ 9 ] = -( cdiff1.z() * cdiff1.z() - cdiff2.z() * cdiff2.z() );
        } else {
            dirs [ 1 ] = 0;
            dirs [ 2 ] = 2;
            dirs [ 3 ] = 4;
            dirs [ 4 ] = 5;
            mults [ 3 ] = ( cdiff1.x() * cdiff1.y() - cdiff2.x() * cdiff2.y() );
            mults [ 4 ] = ( cdiff1.y() * cdiff1.y() - cdiff2.y() * cdiff2.y() );
        }
        jd = new JointDoF(s, 0, vm, dirs, mults);
        constrs->addConstraint(jd);
        //Y DIR
        dirs [ 0 ] = 1;
        if ( dim == 3 ) {
            dirs [ 1 ] = 8;
            dirs [ 2 ] = 1;
            dirs [ 3 ] = 3;
            dirs [ 4 ] = 9;
            dirs [ 5 ] = 16;
            dirs [ 6 ] = 14;
            dirs [ 7 ] = 15;
            dirs [ 8 ] = 13;
            dirs [ 9 ] = 11;
            mults [ 4 ] =  ( cdiff1.x() * cdiff1.z() - cdiff2.x() * cdiff2.z() );
            mults [ 5 ] =  ( cdiff1.y() * cdiff1.z() - cdiff2.y() * cdiff2.z() );
            mults [ 6 ] =  ( cdiff1.z() * cdiff1.z() - cdiff2.z() * cdiff2.z() );
            mults [ 7 ] = -( cdiff1.x() * cdiff1.x() - cdiff2.x() * cdiff2.x() );
            mults [ 8 ] = -( cdiff1.x() * cdiff1.y() - cdiff2.x() * cdiff2.y() );
            mults [ 9 ] = -( cdiff1.x() * cdiff1.z() - cdiff2.x() * cdiff2.z() );
        } else {
            dirs [ 1 ] = 3;
            dirs [ 2 ] = 1;
            dirs [ 3 ] = 4;
            dirs [ 4 ] = 5;
            mults [ 3 ] = -( cdiff1.x() * cdiff1.x() - cdiff2.x() * cdiff2.x() );
            mults [ 4 ] = -( cdiff1.x() * cdiff1.y() - cdiff2.x() * cdiff2.y() );
        }
        jd = new JointDoF(s, 1, vm, dirs, mults);
        constrs->addConstraint(jd);
        //Z DIR
        dirs [ 0 ] = 2;
        if ( dim == 3 ) {
            dirs [ 1 ] = 6;
            dirs [ 2 ] = 4;
            dirs [ 3 ] = 2;
            dirs [ 4 ] = 17;
            dirs [ 5 ] = 10;
            dirs [ 6 ] = 12;
            dirs [ 7 ] = 9;
            dirs [ 8 ] = 16;
            dirs [ 9 ] = 14;
            mults [ 4 ] =  ( cdiff1.x() * cdiff1.x() - cdiff2.x() * cdiff2.x() );
            mults [ 5 ] =  ( cdiff1.y() * cdiff1.x() - cdiff2.y() * cdiff2.x() );
            mults [ 6 ] =  ( cdiff1.z() * cdiff1.x() - cdiff2.z() * cdiff2.x() );
            mults [ 7 ] = -( cdiff1.x() * cdiff1.y() - cdiff2.x() * cdiff2.y() );
            mults [ 8 ] = -( cdiff1.y() * cdiff1.y() - cdiff2.y() * cdiff2.y() );
            mults [ 9 ] = -( cdiff1.y() * cdiff1.z() - cdiff2.y() * cdiff2.z() );
            jd = new JointDoF(s, 2, vm, dirs, mults);
            constrs->addConstraint(jd);
        }
        //connect rotations
        if ( dim == 3 ) {
            vm.resize(4);
            mults.resize(4);
            dirs.resize(4, 0);
            vm [ 0 ] = m;
            vm [ 1 ] = nodes->giveNode(initalNodeNum);
            vm [ 2 ] = nodes->giveNode(initalNodeNum);
            vm [ 3 ] = nodes->giveNode(initalNodeNum);
            mults [ 0 ] = 1.;
            mults [ 1 ] = diff.x();
            mults [ 2 ] = diff.y();
            mults [ 3 ] = diff.z();
        } else {
            vm.resize(3);
            mults.resize(3);
            dirs.resize(3, 0);
            vm [ 0 ] = m;
            vm [ 1 ] = nodes->giveNode(initalNodeNum);
            vm [ 2 ] = nodes->giveNode(initalNodeNum);
            mults [ 0 ] = 1;
            mults [ 1 ] = diff.x();
            mults [ 2 ] = diff.y();
        }
        //X DIR
        if ( dim == 3 ) {
            dirs [ 0 ] = 3;
            dirs [ 1 ] = 9;
            dirs [ 2 ] = 16;
            dirs [ 3 ] = 14;
        } else {
            dirs [ 0 ] = 2;
            dirs [ 1 ] = 4;
            dirs [ 2 ] = 5;
        }
        jd = new JointDoF(s, dirs [ 0 ], vm, dirs, mults);
        constrs->addConstraint(jd);
        //Y DIR
        if ( dim == 3 ) {
            dirs [ 0 ] = 4;
            dirs [ 1 ] = 17;
            dirs [ 2 ] = 10;
            dirs [ 3 ] = 12;
            jd = new JointDoF(s, 4, vm, dirs, mults);
            constrs->addConstraint(jd);
        }
        //Z DIR
        if ( dim == 3 ) {
            dirs [ 0 ] = 5;
            dirs [ 1 ] = 15;
            dirs [ 2 ] = 13;
            dirs [ 3 ] = 11;
            jd = new JointDoF(s, 5, vm, dirs, mults);
            constrs->addConstraint(jd);
        }
    }
}

//////////////////////////////////////////////////////////
void CosseratMechanicalPeriodicBC :: generateExporters(NodeContainer *nodes, ExporterContainer *ex) {
    //export data
    string export_name = "PUCstrain_stress";
    vector< unsigned >n(1, initalNodeNum);
    vector< string >gname( ( dim == 3 )? 18: 6 );
    vector< string >codes( ( dim == 3 )? 18: 6 );
    ForceGauge *fg;
    if ( dim == 2 ) {
        gname [ 0 ] = "sigma_xx";
        gname [ 1 ] = "sigma_yy";
        gname [ 2 ] = "sigma_yx";
        gname [ 3 ] = "sigma_xy";
        gname [ 4 ] = "m_xz";
        gname [ 5 ] = "m_yz";
        codes [ 0 ] = "0";
        codes [ 1 ] = "1";
        codes [ 2 ] = "2";
        codes [ 3 ] = "3";
        codes [ 4 ] = "4";
        codes [ 5 ] = "5";
    } else {
        gname [ 0 ] = "sigma_xx";
        gname [ 1 ] = "sigma_yy";
        gname [ 2 ] = "sigma_zz";
        gname [ 3 ] = "sigma_zy";
        gname [ 4 ] = "sigma_yz";
        gname [ 5 ] = "sigma_zx";
        gname [ 6 ] = "sigma_xz";
        gname [ 7 ] = "sigma_yx";
        gname [ 8 ] = "sigma_xy";
        gname [ 9 ] = "m_xx";
        gname [ 10 ] = "m_yy";
        gname [ 11 ] = "m_zz";
        gname [ 12 ] = "m_zy";
        gname [ 13 ] = "m_yz";
        gname [ 14 ] = "m_zx";
        gname [ 15 ] = "m_xz";
        gname [ 16 ] = "m_yx";
        gname [ 17 ] = "m_xy";
        codes [ 0 ] = "0";
        codes [ 1 ] = "1";
        codes [ 2 ] = "2";
        codes [ 3 ] = "3";
        codes [ 4 ] = "4";
        codes [ 5 ] = "5";
        codes [ 6 ] = "6";
        codes [ 7 ] = "7";
        codes [ 8 ] = "8";
        codes [ 9 ] = "9";
        codes [ 10 ] = "10";
        codes [ 11 ] = "11";
        codes [ 12 ] = "12";
        codes [ 13 ] = "13";
        codes [ 14 ] = "14";
        codes [ 15 ] = "15";
        codes [ 16 ] = "16";
        codes [ 17 ] = "17";
    }
    for ( unsigned i = 0; i < gname.size(); i++ ) {
        fg = new ForceGauge(export_name, gname [ i ], codes [ i ], n, nodes, 1. / volume, dim);
        ex->addExporter(fg);
    }

    DoFGauge *dg;
    if ( dim == 2 ) {
        gname [ 0 ] = "gamma_xx";
        gname [ 1 ] = "gamma_yy";
        gname [ 2 ] = "gamma_yx";
        gname [ 3 ] = "gamma_xy";
        gname [ 4 ] = "kappa_xz";
        gname [ 5 ] = "kappa_yz";
        codes [ 0 ] = "0";
        codes [ 1 ] = "1";
        codes [ 2 ] = "2";
        codes [ 3 ] = "3";
        codes [ 4 ] = "4";
        codes [ 5 ] = "5";
    } else {
        gname [ 0 ] = "gamma_xx";
        gname [ 1 ] = "gamma_yy";
        gname [ 2 ] = "gamma_zz";
        gname [ 3 ] = "gamma_zy";
        gname [ 4 ] = "gamma_yz";
        gname [ 5 ] = "gamma_zx";
        gname [ 6 ] = "gamma_xz";
        gname [ 7 ] = "gamma_yx";
        gname [ 8 ] = "gamma_xy";
        gname [ 9 ] = "kappa_xx";
        gname [ 10 ] = "kappa_yy";
        gname [ 11 ] = "kappa_zz";
        gname [ 12 ] = "kappa_zy";
        gname [ 13 ] = "kappa_yz";
        gname [ 14 ] = "kappa_zx";
        gname [ 15 ] = "kappa_xz";
        gname [ 16 ] = "kappa_yx";
        gname [ 17 ] = "kappa_xy";
        codes [ 0 ] = "0";
        codes [ 1 ] = "1";
        codes [ 2 ] = "2";
        codes [ 3 ] = "3";
        codes [ 4 ] = "4";
        codes [ 5 ] = "5";
        codes [ 6 ] = "6";
        codes [ 7 ] = "7";
        codes [ 8 ] = "8";
        codes [ 9 ] = "9";
        codes [ 10 ] = "10";
        codes [ 11 ] = "11";
        codes [ 12 ] = "12";
        codes [ 13 ] = "13";
        codes [ 14 ] = "14";
        codes [ 15 ] = "15";
        codes [ 16 ] = "16";
        codes [ 17 ] = "17";
    }
    for ( unsigned i = 0; i < gname.size(); i++ ) {
        dg = new DoFGauge(export_name, gname [ i ], codes [ i ], n, nodes, 1., dim);
        ex->addExporter(dg);
    }
}

//////////////////////////////////////////////////////////
void CosseratMechanicalPeriodicBC :: readLoading(istringstream &iss) {
    unsigned num, hnum;
    string param;
    volumetricAverageRigidBC = -1;
    iss >> num;

    strainFunc.resize( ( dim == 3 )? 18: 6, -1 );
    stressFunc.resize( ( dim == 3 )? 18: 6, -1 );

    for ( unsigned i = 0; i < num; i++ ) {
        iss >> param >> hnum;
        if ( dim == 2 ) {
            if ( param.compare("gamma_xx") == 0 ) {
                strainFunc [ 0 ] = hnum;
            } else if ( param.compare("gamma_yy") == 0 ) {
                strainFunc [ 1 ] = hnum;
            } else if ( param.compare("gamma_yx") == 0 ) {
                strainFunc [ 2 ] = hnum;
            } else if ( param.compare("gamma_xy") == 0 ) {
                strainFunc [ 3 ] = hnum;
            } else if ( param.compare("kappa_xz") == 0 ) {
                strainFunc [ 4 ] = hnum;
            } else if ( param.compare("kappa_yz") == 0 ) {
                strainFunc [ 5 ] = hnum;
            } else if ( param.compare("sigma_xx") == 0 ) {
                stressFunc [ 0 ] = hnum;
            } else if ( param.compare("sigma_yy") == 0 ) {
                stressFunc [ 1 ] = hnum;
            } else if ( param.compare("sigma_yx") == 0 ) {
                stressFunc [ 2 ] = hnum;
            } else if ( param.compare("sigma_xy") == 0 ) {
                stressFunc [ 3 ] = hnum;
            } else if ( param.compare("m_xz") == 0 ) {
                stressFunc [ 4 ] = hnum;
            } else if ( param.compare("m_yz") == 0 ) {
                stressFunc [ 5 ] = hnum;
            } else if ( param.compare("volumetricAverage") == 0 ) {
                volumetricAverageRigidBC = hnum;
            } else {
                cout << "Error in " << name << " : loading by " << param << " not implemented yet" << '\n';
                exit(1);
            }
        } else if ( dim == 3 ) {
            if ( param.compare("gamma_xx") == 0 ) {
                strainFunc [ 0 ] = hnum;
            } else if ( param.compare("gamma_yy") == 0 ) {
                strainFunc [ 1 ] = hnum;
            } else if ( param.compare("gamma_zz") == 0 ) {
                strainFunc [ 2 ] = hnum;
            } else if ( param.compare("gamma_zy") == 0 ) {
                strainFunc [ 3 ] = hnum;
            } else if ( param.compare("gamma_yz") == 0 ) {
                strainFunc [ 4 ] = hnum;
            } else if ( param.compare("gamma_zx") == 0 ) {
                strainFunc [ 5 ] = hnum;
            } else if ( param.compare("gamma_xz") == 0 ) {
                strainFunc [ 6 ] = hnum;
            } else if ( param.compare("gamma_yx") == 0 ) {
                strainFunc [ 7 ] = hnum;
            } else if ( param.compare("gamma_xy") == 0 ) {
                strainFunc [ 8 ] = hnum;
            } else if ( param.compare("kappa_xx") == 0 ) {
                strainFunc [ 9 ] = hnum;
            } else if ( param.compare("kappa_yy") == 0 ) {
                strainFunc [ 10 ] = hnum;
            } else if ( param.compare("kappa_zz") == 0 ) {
                strainFunc [ 11 ] = hnum;
            } else if ( param.compare("kappa_zy") == 0 ) {
                strainFunc [ 12 ] = hnum;
            } else if ( param.compare("kappa_yz") == 0 ) {
                strainFunc [ 13 ] = hnum;
            } else if ( param.compare("kappa_zx") == 0 ) {
                strainFunc [ 14 ] = hnum;
            } else if ( param.compare("kappa_xz") == 0 ) {
                strainFunc [ 15 ] = hnum;
            } else if ( param.compare("kappa_yx") == 0 ) {
                strainFunc [ 16 ] = hnum;
            } else if ( param.compare("kappa_xy") == 0 ) {
                strainFunc [ 17 ] = hnum;
            } else if ( param.compare("sigma_xx") == 0 ) {
                stressFunc [ 0 ] = hnum;
            } else if ( param.compare("sigma_yy") == 0 ) {
                stressFunc [ 1 ] = hnum;
            } else if ( param.compare("sigma_zz") == 0 ) {
                stressFunc [ 2 ] = hnum;
            } else if ( param.compare("sigma_zx") == 0 ) {
                stressFunc [ 3 ] = hnum;
            } else if ( param.compare("sigma_yz") == 0 ) {
                stressFunc [ 4 ] = hnum;
            } else if ( param.compare("sigma_zx") == 0 ) {
                stressFunc [ 5 ] = hnum;
            } else if ( param.compare("sigma_xz") == 0 ) {
                stressFunc [ 6 ] = hnum;
            } else if ( param.compare("sigma_yx") == 0 ) {
                stressFunc [ 7 ] = hnum;
            } else if ( param.compare("sigma_xy") == 0 ) {
                stressFunc [ 8 ] = hnum;
            } else if ( param.compare("m_xx") == 0 ) {
                stressFunc [ 9 ] = hnum;
            } else if ( param.compare("m_yy") == 0 ) {
                stressFunc [ 10 ] = hnum;
            } else if ( param.compare("m_zz") == 0 ) {
                stressFunc [ 11 ] = hnum;
            } else if ( param.compare("m_zx") == 0 ) {
                stressFunc [ 12 ] = hnum;
            } else if ( param.compare("m_yz") == 0 ) {
                stressFunc [ 13 ] = hnum;
            } else if ( param.compare("m_zx") == 0 ) {
                stressFunc [ 14 ] = hnum;
            } else if ( param.compare("m_xz") == 0 ) {
                stressFunc [ 15 ] = hnum;
            } else if ( param.compare("m_yx") == 0 ) {
                stressFunc [ 16 ] = hnum;
            } else if ( param.compare("m_xy") == 0 ) {
                stressFunc [ 17 ] = hnum;
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
void CosseratMechanicalPeriodicBC :: generateRigidBodyBC(NodeContainer *nodes, ElementContainer *elems, BCContainer *bcs, ConstraintContainer *constrs, FunctionContainer *funcs) {
    if ( 1 ) {  //last master node cannot move
        Node *m = constrs->giveConstraint(constrs->giveConstraintsSize() - 1)->giveMasterNode(0);// warning C4267: 'argument': conversion from 'size_t' to 'const unsigned int', possible loss of data
        BoundaryCondition *bc;
        vector< int >dBC, nBC;
        dBC.resize(m->giveNumberOfDoFs(), -1);         //todo: warning C4267: 'argument': conversion from 'size_t' to 'const _Ty', possible loss of data
        for ( unsigned k = 0; k < dim; k++ ) {
            dBC [ k ] = funcs->giveSize();
        }
        nBC.resize(m->giveNumberOfDoFs(), -1);
        bc = new BoundaryCondition(m, dBC, nBC);
        bcs->addBoundaryCondition(bc);

        //add constant function
        vector< double >x, y;
        x.resize(1, 0);
        y.resize(1, 0);
        PieceWiseLinearFunction *newf = new PieceWiseLinearFunction(x, y);
        funcs->addFunction(newf);
    }
    if ( 1 ) {  //volumetric average for rotations
        VolumetricAverage *va;
        vector< Node * >vm;

        for ( unsigned n = 0; n < nodes->giveSize(); n++ ) {
            if ( nodes->giveNode(n)->doesMechanics() && dynamic_cast< Particle * >( nodes->giveNode(n) ) ) {
                vm.push_back( nodes->giveNode(n) );
            }
        }
        if ( vm.size() > 0 ) {
            unsigned skip = dim;
            unsigned nDoFs = 3;
            if ( dim == 3 ) {
                nDoFs = 6;
            }
            MechDoF *pn = new MechDoF(dim, nDoFs - skip);
            nodes->addNode(pn);

            vector< unsigned >dirs( vm.size() );

            for ( unsigned vi = 0; vi < nDoFs - skip; vi++ ) {//only rotations
                fill(dirs.begin(), dirs.end(), vi + skip);
                va = new VolumetricAverage(vm, dirs, pn, vi, elems, constrs);
                constrs->addConstraint(va);
            }

            BoundaryCondition *bc;
            vector< int >dBC, nBC;
            dBC.resize( nDoFs - skip, funcs->giveSize() );
            nBC.resize(nDoFs - skip, -1);
            bc = new BoundaryCondition(pn, dBC, nBC);
            bcs->addBoundaryCondition(bc);

            //add constant function
            vector< double >x, y;
            x.resize(1, 0);
            y.resize(1, 0);
            PieceWiseLinearFunction *newf = new PieceWiseLinearFunction(x, y);
            funcs->addFunction(newf);
        }
    }
}
