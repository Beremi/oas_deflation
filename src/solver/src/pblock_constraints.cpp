#include "pblock_constraints.h"
#include "model.h"
#include "material_vectorial.h"
#include "element_ldpm.h"
#include "solver_implicit.h"
#include "pblock_periodic_bc.h"

// #include "misc.h"  // TODO JK: this include causes linking error
using namespace std;


std :: vector< double >PointToStdVector(const Point &p, unsigned dim = 3) {
    std :: vector< double >vect;
    vect.push_back( p.x() );
    vect.push_back( p.y() );
    if ( dim == 3 ) {
        vect.push_back( p.z() );
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

void MaterialRegion :: apply(Model *model) {
    MaterialContainer *mats = model->giveMaterials();
    ElementContainer *elems = model->giveElements();

    // TODO make sure material is sutable for mech / transport
    unsigned num_nodes_inside;
    for ( auto &el : * elems ) {
        num_nodes_inside = 0;
        for ( auto const &nod : el->giveNodes() ) {
            if ( ( !this->transport && nod->doesMechanics() ) ||
                 ( this->transport && nod->doesTransport() ) ) {
                // TODO JK: are there any elems that are mechanical and connect transport nodes and same for transport elems?
                if ( this->reg->isInside(nod->givePoint() ) ) {
                    if ( this->all_nodes ) {
                        num_nodes_inside++;
                    } else {
                        el->changeMaterial(mats->giveMaterial(this->material_id) );
                        break;
                    }
                }
            }
        }
        if ( this->all_nodes && num_nodes_inside == el->giveNodes().size() ) {
            el->changeMaterial(mats->giveMaterial(this->material_id) );
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
void VoigtConstraint :: apply(Model *model) {
    BCContainer *bcs = model->giveBoundaryConditions();
    NodeContainer *nodes = model->giveNodes();
    ExporterContainer *ex = model->giveExporters();
    ConstraintContainer *constrs = model->giveConstraints();
    FunctionContainer *funcs = model->giveFunctions();

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
void PressureFromMechanicalLoad :: apply(Model *model) {
    ConstraintContainer *constrs = model->giveConstraints();
    MaterialContainer *mats = model->giveMaterials();
    NodeContainer *nodes = model->giveNodes();

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
void RigidPlate :: apply(Model *model) {
    NodeContainer *nodes = model->giveNodes();
    ConstraintContainer *constrs = model->giveConstraints();
    RegionContainer *regions = model->giveRegions();

    Node *master, *slave;
    // read the line "masterId numSlaves slaveId1, slaveId2...."
    master = nodes->giveNode(master_id);
    // check if it is master node
    this->checkPhysicalField(master);

    for ( auto const &sl_id : slave_ids ) {
        slave = nodes->giveNode(sl_id);
        connectSlaveMasterRigid(constrs, slave, master, this->dim, activeDirs);
    }


    if ( insideRegions.size() > 0 || outsideRegions.size() > 0 ) {
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
void CoordRigidPlate :: apply(Model *model) {
    NodeContainer *nodes = model->giveNodes();
    ConstraintContainer *constrs = model->giveConstraints();

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
void RingRigidPlate :: apply(Model *model) {
    NodeContainer *nodes = model->giveNodes();
    ConstraintContainer *constrs = model->giveConstraints();

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
void ExpansionRing :: apply(Model *model) {
    NodeContainer *nodes = model->giveNodes();
    ConstraintContainer *constrs = model->giveConstraints();
    FunctionContainer *funcs = model->giveFunctions();

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

void ExpansionRingDoFLoad :: apply(Model *model) {
    NodeContainer *nodes = model->giveNodes();
    ConstraintContainer *constrs = model->giveConstraints();

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


void ExpansionRingSingleDoFLoad :: apply(Model *model) {
    NodeContainer *nodes = model->giveNodes();
    ConstraintContainer *constrs = model->giveConstraints();

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
void NormalSurfaceLoad :: apply(Model *model) {
    BCContainer *b = model->giveBoundaryConditions();
    RegionContainer *regions = model->giveRegions();
    ElementContainer *e = model->giveElements();

    RigidBodyBoundary *rbb;
    LDPMTetra *tet;
    BoundaryCondition *bc;
    vector< int >dBC, nBC;
    Point diff;
    dBC.resize( ( dim - 1 ) * 3, -1);
    nBC.resize( ( dim - 1 ) * 3, fnID);
    vector< double >mult;
    mult.resize( ( dim - 1 ) * 3, 0.);
    Point normal;
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

            area = triArea3D(innodes [ 0 ]->givePointPointer(), innodes [ 1 ]->givePointPointer(), innodes [ 2 ]->givePointPointer() );
            normal = ( innodes [ 1 ]->givePoint() - innodes [ 0 ]->givePoint() ).cross(innodes [ 2 ]->givePoint() - innodes [ 0 ]->givePoint() );
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
void MechHangingNode :: apply(Model *model) {
    Element *ee = model->giveElements()->giveElement(elemid);
    Node *nn = model->giveNodes()->giveNode(nodeid);
    model->giveConstraints()->addHangingNodeConstraint(nn, ee);

    /*
     *  ElementContainer * e = model->giveElements();
     *  NodeContainer * n = model->giveNodes();
     *  ConstraintContainer * c = model->giveConstraints();
     *
     *  Element *ee = e->giveElement(elemid);
     *  unsigned ndim = ee->giveDimension();
     *  Point natcoords = ee->findNaturalCoords(n->giveNode(nodeid)->givePointPointer() );
     *  Vector weights = ee->giveShapeFunctions(& natcoords);
     *  vector< Node * >masters = ee->giveNodes();
     *  vector< unsigned >dirs(masters.size() );
     *  vector< double >mults(masters.size() );
     *
     *  //displacements
     *  for ( unsigned k = 0; k < masters.size(); k++ ) {
     *      mults [ k ] = weights [ k ];
     *  }
     *  unsigned i = 0;
     *  for ( ; i < ndim; i++ ) {
     *      for ( unsigned k = 0; k < masters.size(); k++ ) {
     *          dirs [ k ] = i;
     *      }
     *      JointDoF *newJD = new JointDoF(n->giveNode(nodeid), i, masters, dirs, mults);
     *      c->addConstraint(newJD);
     *  }
     */
}

//////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////


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
                            multipliers.push_back( -( slave->givePoint().y() - master->givePoint().y() ) );
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
                            multipliers.push_back( -( slave->givePoint().y() - master->givePoint().y() ) );
                        } else if ( kslave == 1 ) {
                            directions.push_back(3);
                            multipliers.push_back( -( slave->givePoint().z() - master->givePoint().z() ) );
                            directions.push_back(5);
                            multipliers.push_back( ( slave->givePoint().x() - master->givePoint().x() ) );
                        } else {
                            directions.push_back(3);
                            multipliers.push_back( ( slave->givePoint().y() - master->givePoint().y() ) );
                            directions.push_back(4);
                            multipliers.push_back( -( slave->givePoint().x() - master->givePoint().x() ) );
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
