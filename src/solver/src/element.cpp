#include "element.h"
#include "element_container.h"

//////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////
// BASIC ELEMENT - MASTER CLASS
Element :: ~Element() {
    for ( vector< MaterialStatus * > :: iterator e = stats.begin(); e != stats.end(); ++e ) {
        delete * e;
    }
}

//////////////////////////////////////////////////////////
void Element :: init() {
    unsigned totalDoFs = 0;
    for ( vector< Node * > :: const_iterator n = nodes.begin(); n != nodes.end(); ++n ) {
        totalDoFs += ( * n )->giveNumberOfDoFs();
    }
    DoFids.resize(totalDoFs);
    unsigned i = 0;
    unsigned k;
    for ( vector< Node * > :: const_iterator n = nodes.begin(); n != nodes.end(); ++n ) {
        k = ( * n )->giveStartingDoF();
        for ( unsigned s = 0; s < ( * n )->giveNumberOfDoFs(); s++, i++ ) {
            DoFids [ i ] = k + s;
        }
    }
    outDoFs = totalDoFs; //basic elems will alway have input = output
}

//////////////////////////////////////////////////////////
void Element :: initMaterialStatuses() {
    for ( vector< MaterialStatus * > :: iterator m = stats.begin(); m != stats.end(); ++m ) {
        ( * m )->init();
    }
}

//////////////////////////////////////////////////////////
void Element :: updateMaterialStatuses() {
    for ( vector< MaterialStatus * > :: iterator m = stats.begin(); m != stats.end(); ++m ) {
        ( * m )->update();
    }
}

//////////////////////////////////////////////////////////
double Element :: giveValue(string code) const {
    ( void ) code;
    return 0;
};

//////////////////////////////////////////////////////////
double Element :: giveIPValue(string code, unsigned ipnum) const {
    if ( code.compare("x") == 0 ) {
        return ip_locs [ ipnum ].x;
    } else if ( code.compare("y") == 0 ) {
        return ip_locs [ ipnum ].y;
    } else if ( code.compare("z") == 0 ) {
        return ip_locs [ ipnum ].z;
    } else {
        return stats [ ipnum ]->giveValue(code);
    }
};

//////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////
// RBSN ELEMENT
RigidBodyContact :: RigidBodyContact(const unsigned dim) {
    ndim = dim;
    nodes.resize(2);
    // tangs.resize(dim-1);
    ip_locs.resize(1);
    stats.resize(1);
    name = "RigidBodyContact";
}

//////////////////////////////////////////////////////////
double RigidBodyContact :: giveValue(string code) const {
    if ( code.compare("damage") == 0 ) {
        DisMechMaterialStatus *dmstats = static_cast< DisMechMaterialStatus * >( stats [ 0 ] );
        return dmstats->giveValue(code);
    } else {
        return MechanicalElement :: giveValue(code);
    }
}

//////////////////////////////////////////////////////////
double RigidBodyContact :: giveIPValue(string code, unsigned ipnum) const {
    if ( code.compare("normal_x") == 0 ) {
        return normal.getX();
    } else if ( code.compare("normal_y") == 0 ) {
        return normal.getY();
    } else if ( code.compare("normal_z") == 0 ) {
        return normal.getZ();
    } else if ( code.compare("t1_x") == 0 ) {
        return R [ 1 ] [ 0 ];
    } else if ( code.compare("t1_y") == 0 ) {
        return R [ 1 ] [ 1 ];
    } else if ( code.compare("t1_z") == 0 ) {
        return R [ 1 ] [ 2 ];
    } else if ( code.compare("t2_x") == 0 ) {
        return R [ 2 ] [ 0 ];
    } else if ( code.compare("t2_y") == 0 ) {
        return R [ 2 ] [ 1 ];
    } else if ( code.compare("t2_z") == 0 ) {
        return R [ 2 ] [ 2 ];
    } else if ( code.compare("volume") == 0 ) {
        return area * length / ndim;
    } else if ( code.compare("area") == 0 ) {
        return area;
    } else if ( code.compare("length") == 0 ) {
        return length;
    } else if ( code.compare("energy_total_elem") == 0 ) {
        // TODO all following repair in the way that anything ending by _elem would be multiplied by area mat->giveValue(code without '_elem') * ( area * length / ndim )
        MaterialStatus *fmstat = static_cast< MaterialStatus * >( stats [ ipnum ] );
        return fmstat->giveValue("energy_total")  * ( area * length / ndim );
    } else if ( code.compare("energy_totalT_elem") == 0 ) {
        MaterialStatus *fmstat = static_cast< MaterialStatus * >( stats [ ipnum ] );
        return fmstat->giveValue("energy_totalT")  * ( area * length / ndim );
    } else if ( code.compare("energy_totalN_elem") == 0 ) {
        MaterialStatus *fmstat = static_cast< MaterialStatus * >( stats [ ipnum ] );
        return fmstat->giveValue("energy_totalN")  * ( area * length / ndim );
    } else if ( code.compare("energy_PLN_elem") == 0 ) {
        MaterialStatus *fmstat = static_cast< MaterialStatus * >( stats [ ipnum ] );
        return fmstat->giveValue("energy_PLN")  * ( area * length / ndim );
    } else if ( code.compare("energy_DN_elem") == 0 ) {
        MaterialStatus *fmstat = static_cast< MaterialStatus * >( stats [ ipnum ] );
        return fmstat->giveValue("energy_DN")  * ( area * length / ndim );
    } else if ( code.compare("energy_KinN_elem") == 0 ) {
        MaterialStatus *fmstat = static_cast< MaterialStatus * >( stats [ ipnum ] );
        return fmstat->giveValue("energy_KinN")  * ( area * length / ndim );
    } else if ( code.compare("energy_IsoN_elem") == 0 ) {
        MaterialStatus *fmstat = static_cast< MaterialStatus * >( stats [ ipnum ] );
        return fmstat->giveValue("energy_IsoN")  * ( area * length / ndim );
    } else if ( code.compare("energy_PLT_elem") == 0 ) {
        MaterialStatus *fmstat = static_cast< MaterialStatus * >( stats [ ipnum ] );
        return fmstat->giveValue("energy_PLT")  * ( area * length / ndim );
    } else if ( code.compare("energy_DT_elem") == 0 ) {
        MaterialStatus *fmstat = static_cast< MaterialStatus * >( stats [ ipnum ] );
        return fmstat->giveValue("energy_DT")  * ( area * length / ndim );
    } else if ( code.compare("energy_KinT_elem") == 0 ) {
        MaterialStatus *fmstat = static_cast< MaterialStatus * >( stats [ ipnum ] );
        return fmstat->giveValue("energy_KinT")  * ( area * length / ndim );
    } else if ( code.compare("energy_IsoT_elem") == 0 ) {
        MaterialStatus *fmstat = static_cast< MaterialStatus * >( stats [ ipnum ] );
        return fmstat->giveValue("energy_IsoT")  * ( area * length / ndim );
    } else if ( code.compare("work_dissip_elem") == 0 ) {
        MaterialStatus *fmstat = static_cast< MaterialStatus * >( stats [ ipnum ] );
        return fmstat->giveValue("work_dissip")  * ( area * length / ndim );
    } else if ( code.compare("work_dissipN_elem") == 0 ) {
        MaterialStatus *fmstat = static_cast< MaterialStatus * >( stats [ ipnum ] );
        return fmstat->giveValue("work_dissipN")  * ( area * length / ndim );
    } else if ( code.compare("work_dissipT_elem") == 0 ) {
        MaterialStatus *fmstat = static_cast< MaterialStatus * >( stats [ ipnum ] );
        return fmstat->giveValue("work_dissipT")  * ( area * length / ndim );
    } else if ( code.compare("work_total_elem") == 0 ) {
        MaterialStatus *fmstat = static_cast< MaterialStatus * >( stats [ ipnum ] );
        return fmstat->giveValue("work_total")  * ( area * length / ndim );
    } else if ( code.compare("work_ela_elem") == 0 ) {
        MaterialStatus *fmstat = static_cast< MaterialStatus * >( stats [ ipnum ] );
        return fmstat->giveValue("work_ela")  * ( area * length / ndim );
    } else if ( code.compare("work_totalT_elem") == 0 ) {
        MaterialStatus *fmstat = static_cast< MaterialStatus * >( stats [ ipnum ] );
        return fmstat->giveValue("work_totalT")  * ( area * length / ndim );
    } else if ( code.compare("work_elaT_elem") == 0 ) {
        MaterialStatus *fmstat = static_cast< MaterialStatus * >( stats [ ipnum ] );
        return fmstat->giveValue("work_elaT")  * ( area * length / ndim );
    } else if ( code.compare("work_totalN_elem") == 0 ) {
        MaterialStatus *fmstat = static_cast< MaterialStatus * >( stats [ ipnum ] );
        return fmstat->giveValue("work_totalN")  * ( area * length / ndim );
    } else if ( code.compare("work_elaN_elem") == 0 ) {
        MaterialStatus *fmstat = static_cast< MaterialStatus * >( stats [ ipnum ] );
        return fmstat->giveValue("work_elaN")  * ( area * length / ndim );
    } else {
        return MechanicalElement :: giveIPValue(code, ipnum);
    }
}

//////////////////////////////////////////////////////////
void RigidBodyContact :: readFromLine(istringstream &iss, NodeContainer *fullnodes, MaterialContainer *fullmatrs) {
    unsigned num, num2;
    iss >> num;
    nodes [ 0 ] = fullnodes->giveNode(num);
    iss >> num;
    nodes [ 1 ] = fullnodes->giveNode(num);
    iss >> num;
    vert.resize(num);
    for ( unsigned i = 0; i < num; i++ ) {
        iss >> num2;
        vert [ i ] = fullnodes->giveNode(num2);
    }
    iss >> num;
    mat = fullmatrs->giveMaterial(num);
}

//////////////////////////////////////////////////////////

void RigidBodyContact :: checkNodeType() const {
    //check that nodes are particles
    for ( unsigned i = 0; i < 2; i++ ) {
        Particle *p = dynamic_cast< Particle * >( nodes [ i ] );
        if ( !p ) {
            cerr << "Error in " << name << ": nodes must be inherited from Particle, " << nodes [ i ]->giveName() << " provided" << endl;
            exit(1);
        }
    }
}

//////////////////////////////////////////////////////////

Matrix RigidBodyContact :: giveBMatrix() const {
    //Matrix B
    Matrix B = Matrix(ndim, 6 * ( ndim - 1 ) );
    Matrix Aa = giveAMatrix(nodes [ 0 ]->givePoint(), ip_locs [ 0 ]) * ( -1. );
    Matrix Ab = giveAMatrix(nodes [ 1 ]->givePoint(), ip_locs [ 0 ]);
    for ( unsigned i = 0; i < ndim; i++ ) {
        for ( unsigned j = 0; j < 3 * ( ndim - 1 ); j++ ) {
            B [ i ] [ j ] = Aa [ i ] [ j ];
            B [ i ] [ j + 3 * ( ndim - 1 ) ] = Ab [ i ] [ j ];
        }
    }
    return B / length;
}

//////////////////////////////////////////////////////////
void RigidBodyContact :: init() {
    Element :: init(); //calling base class method;

    checkNodeType();

    //check that material is DisMechMat
    DisMechMaterial *p = dynamic_cast< DisMechMaterial * >( mat );
    if ( !p ) {
        cerr << "Error in " << name << ": material must be inherited from DisMechMaterial, " << mat->giveName() << " provided" << endl;
        exit(1);
    }

    Point t;
    if ( ndim == 2 ) {
        if ( !( vert.size() == 2 ) ) {
            cerr << "Error: exactly 2 vertices must be involved, " << vert.size() << " provided" << endl;
            exit(1);
        }

        ip_locs [ 0 ] = ( vert [ 0 ]->givePoint() + vert [ 1 ]->givePoint() ) / 2.;
        t = vert [ 1 ]->givePoint() - vert [ 0 ]->givePoint();
        area = t.norm();
        t = t / area;
    } else {
        //JM: Coplanarity check for vertices on the face
        //JM: checking coplanarity of every consecutive 4 nodes
        double maxErr = 0.0;
        double currErr = 0.0;
        //
        for ( unsigned int i = 0; i < vert.size() - 3; i++ ) {
            currErr = checkCoplanarity(vert [ i ]->givePoint(), vert [ i + 1 ]->givePoint(), vert [ i + 2 ]->givePoint(), vert [ i + 3 ]->givePoint() );
            if ( abs(currErr) > maxErr ) {
                maxErr = abs(currErr);
            }
        }
        //JM: also checking if the beam midpoint is coplanar with the face
        Point midPoint = ( nodes [ 1 ]->givePoint() + nodes [ 0 ]->givePoint() ) / 2.;
        currErr = checkCoplanarity(vert [ 0 ]->givePoint(), vert [ 1 ]->givePoint(), vert [ 2 ]->givePoint(), midPoint);
        if ( abs(currErr) > maxErr ) {
            maxErr = abs(currErr);
        }
        //
        if ( maxErr > 1e-10 ) {
            cerr << "Vertices are not coplanar!!! Coplanarity error: " << maxErr << endl;
            exit(1);
        }

        //JM: face normal vector made from first 3 vertices
        //JM: coordinate swap for tangential vector according to https://orbit.dtu.dk/files/126824972/onb_frisvad_jgt2012_v2.pdf
        Point n = cross(vert [ 1 ]->givePoint() - vert [ 0 ]->givePoint(), vert [ 2 ]->givePoint() - vert [ 0 ]->givePoint() );
        n /= n.norm();
        Point t2;
        if ( fabs(n.x) > fabs(n.z) ) {
            t2 = Point(-n.y, n.x, 0.0f);
        } else {
            t2 = Point(0.0f, -n.z, n.y);
        }
        t = cross(t2, n);
        t /= t.norm();

        //JM: Perpendicularity check of the beam and face directions
        //JM: normal of the face surface taken from first 3 vertices is (B - A) x (C - A)
        //JM: perpendicularity check: cross (beam, face)=>0
        Point prp = ( nodes [ 1 ]->givePoint() - nodes [ 0 ]->givePoint() ) * t;
        if ( prp.norm() > 1e-10 ) {
            cerr << "Face surface is not perpendicular to beam direction!!! Error: " << prp.norm() << endl;
            //  exit(1);
        }

        //JM: finding position of the SINGLE integration point -> center of gravity of the face polygon
        //JM: average point of the polygon for triangulation
        Point avgPoint = Point(0.0, 0.0, 0.0);
        for ( unsigned int i = 0; i < vert.size(); i++ ) {
            avgPoint += vert [ i ]->givePoint();
        }
        avgPoint /= vert.size();

        //JM: integration point coordinates as an average of CGs of face triangles weighted by areas
        ip_locs [ 0 ] = Point(0.0, 0.0, 0.0);
        area = 0.0;
        double ai = 0.0;
        unsigned int j = 0;
        for ( unsigned int i = 0; i < vert.size(); i++ ) {
            j = i + 1;
            if ( i == vert.size() - 1 ) {
                j = 0;
            }
            //triangle area computed as a_i = norm(cross(AB, AC)) / 2
            ai = ( cross(vert [ i ]->givePoint() - avgPoint,   vert [ j ]->givePoint() - avgPoint) ).norm();
            area += ai;
            //triangle cg_i is an average of simplex vertices, adding to CG coordinates multiplied by a_i weight
            ip_locs [ 0 ] += ( avgPoint + vert [ i ]->givePoint() + vert [ j ]->givePoint() ) / 3.0 * ai;
        }
        ip_locs [ 0 ] /= area;

        //JM: Check if integration point is coplanar with face
        currErr = checkCoplanarity(vert [ 0 ]->givePoint(), vert [ 1 ]->givePoint(), vert [ 2 ]->givePoint(), ip_locs [ 0 ]);
        if ( abs(currErr) > 1e-6 ) {
            cerr << "Integration point is not coplanar with the face!!! Coplanarity error: " << currErr << endl;
            exit(1);
        }
    }

    stats [ 0 ] = mat->giveNewMaterialStatus(this);
    normal = nodes [ 1 ]->givePoint() - nodes [ 0 ]->givePoint();
    length = normal.norm();
    normal = normal / length;
    if ( abs(normal * t) > 1e-8 ) {
        cout << vert [ 0 ]->givePoint().x << " " <<  vert [ 0 ]->givePoint().y <<  " X " << vert [ 1 ]->givePoint().x << " " <<  vert [ 1 ]->givePoint().y << endl;
        cout << nodes [ 0 ]->givePoint().x << " " <<  nodes [ 0 ]->givePoint().y <<  " X " << nodes [ 1 ]->givePoint().x << " " <<  nodes [ 1 ]->givePoint().y << endl;
        cerr << "Error: normal and contact vector are not parallel, error " << normal * t << " normal v." << normal.x << " " << normal.y << " contact v. " << t.x << " " << t.y << endl;
        // exit(1);
    }

    //Matrices according to habilitation of Jan Elias (2017, page 42): https://www.vutbr.cz/www_base/vutdisk.php?i=103116a130

    Matrix B = giveBMatrix();

    // Matrix R;
    if ( ndim == 2 ) {
        Point t1 = Point(-normal.y, normal.x);
        R = Matrix(2, 2);
        R [ 0 ] [ 0 ] = normal.x;
        R [ 0 ] [ 1 ] = normal.y;
        R [ 1 ] [ 0 ] = t1.x;
        R [ 1 ] [ 1 ] = t1.y;
    } else if ( ndim == 3 ) {
        Point t1, t2;
        Point arbit(sqrt(2.), -sqrt(3.), M_PI);
        if ( ( normal - arbit ).norm() < 1e-3 ) {
            t1 = cross(arbit, normal);
            t1.normalize();
            t2 = cross(normal, t1);
            t2.normalize();
        } else {
            // the following results in zeros in stiffness matrix in case of normal in direstion of any of global base axes
            if ( abs(normal.x) > 1e-3 ) {
                t1 = Point(-normal.y / normal.x, 1, 0);
            } else if ( abs(normal.y) > 1e-3 ) {
                t1 = Point(0, -normal.z / normal.y, 1);
            } else {
                t1 = Point(1, 0, -normal.x / normal.z);
            }
        }
        t1 = t1 / t1.norm();
        t2 = cross(normal, t1);
        R = Matrix(3, 3);
        R [ 0 ] [ 0 ] = normal.x;
        R [ 0 ] [ 1 ] = normal.y;
        R [ 0 ] [ 2 ] = normal.z;
        R [ 1 ] [ 0 ] = t1.x;
        R [ 1 ] [ 1 ] = t1.y;
        R [ 1 ] [ 2 ] = t1.z;
        R [ 2 ] [ 0 ] = t2.x;
        R [ 2 ] [ 1 ] = t2.y;
        R [ 2 ] [ 2 ] = t2.z;
    } else {
        cerr << "Error - RigidBodyContact: dimension " << ndim << "not implemented" << endl;
        exit(EXIT_FAILURE);
    }
    GeomM = R * B;
}

//////////////////////////////////////////////////////////
Matrix RigidBodyContact :: giveAMatrix(Point a, Point x) const {
    Matrix A(ndim, 3 * ( ndim - 1 ) );
    if ( ndim == 3 ) {
        A [ 0 ] [ 0 ] = A [ 1 ] [ 1 ] = A [ 2 ] [ 2 ] = 1;
        A [ 1 ] [ 3 ] = a.z - x.z;
        A [ 0 ] [ 4 ] = -A [ 1 ] [ 3 ];
        A [ 2 ] [ 3 ] = x.y - a.y;
        A [ 0 ] [ 5 ] = -A [ 2 ] [ 3 ];
        A [ 2 ] [ 4 ] = a.x - x.x;
        A [ 1 ] [ 5 ] = -A [ 2 ] [ 4 ];
    } else if ( ndim == 2 ) {
        A [ 0 ] [ 0 ] = A [ 1 ] [ 1 ] = 1;
        A [ 0 ] [ 2 ] = a.y - x.y;
        A [ 1 ] [ 2 ] = x.x - a.x;
    } else {
        cerr << "Error - RigidBodyContact: dimension " << ndim << "not implemented" << endl;
        exit(EXIT_FAILURE);
    }
    return A;
}

//////////////////////////////////////////////////////////
Matrix RigidBodyContact :: giveStiffnessMatrix(string matrixType) const {
    DisMechMaterialStatus *dmstats = static_cast< DisMechMaterialStatus * >( stats [ 0 ] );
    Vector matstiffness;
    matstiffness = dmstats->giveNormalShearStiffness(matrixType);
    Matrix Alpha(ndim, ndim);
    Alpha [ 0 ] [ 0 ] = matstiffness [ 0 ]; //E0
    for ( unsigned i = 1; i < ndim; i++ ) {
        Alpha [ i ] [ i ] = matstiffness [ 1 ];                  //E0*alpha
    }
    Matrix K = GeomM.transpose() * Alpha * GeomM;
    return K * ( length * area );
}

//////////////////////////////////////////////////////////
Matrix RigidBodyContact :: giveMassMatrix() const {
    Matrix M(12, 12);
    DisMechMaterialStatus *tstats = static_cast< DisMechMaterialStatus * >( stats [ 0 ] );
    ( void ) tstats;
    //TO BE DONE
    return M;
}

//////////////////////////////////////////////////////////
Vector RigidBodyContact :: giveInternalForces(const Vector &DoFs, bool frozen) const {
    if ( frozen ) { //frozen internal variables
        return giveStiffnessMatrix("secant") * DoFs;
    } else  {    //evalution of internal variables
        Vector strainNT = giveContactStrainNT(DoFs);
        DisMechMaterialStatus *dmstats = static_cast< DisMechMaterialStatus * >( stats [ 0 ] );
        Vector stressNT = dmstats->giveStress(strainNT);
        return GeomM.transpose() * ( stressNT * ( length * area ) );
    }
};

//////////////////////////////////////////////////////////
Vector RigidBodyContact :: giveContactStrainNT(const Vector &DoFs) const {
    return GeomM * DoFs;
};

//////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////
// TRUSS ELEMENT
Matrix Truss :: giveAMatrix(Point a, Point x) const {
    ( void ) a;
    ( void ) x;
    Matrix A(ndim, ndim);
    if ( ndim == 3 ) {
        A [ 0 ] [ 0 ] = A [ 1 ] [ 1 ] = A [ 2 ] [ 2 ] = 1;
    } else if ( ndim == 2 ) {
        A [ 0 ] [ 0 ] = A [ 1 ] [ 1 ] = 1;
    } else {
        cerr << "Error - Truss: dimension " << ndim << "not implemented" << endl;
        exit(EXIT_FAILURE);
    }
    return A;
}

//////////////////////////////////////////////////////////
Matrix Truss :: giveStiffnessMatrix(string matrixType) const {
    DisMechMaterialStatus *dmstats = static_cast< DisMechMaterialStatus * >( stats [ 0 ] );
    Vector matstiffness;
    matstiffness = dmstats->giveNormalShearStiffness(matrixType);
    Matrix Alpha(ndim, ndim);
    Alpha [ 0 ] [ 0 ] = matstiffness [ 0 ]; //E0
    for ( unsigned i = 1; i < ndim; i++ ) {
        Alpha [ i ] [ i ] = 0;                  //no shear stiffness
    }
    Matrix K = GeomM.transpose() * Alpha * GeomM;
    return K * ( length * area );
}

//////////////////////////////////////////////////////////
Matrix Truss :: giveMassMatrix() const {
    Matrix S(6, 6);
    //TrsprtMaterialStatus *tstats = static_cast<TrsprtMaterialStatus *>(stats[0]);
    //TO BE DONE
    return S;
}

//////////////////////////////////////////////////////////
void Truss :: checkNodeType() const {
    //check that nodes are mechanical nodes
    for ( unsigned i = 0; i < 2; i++ ) {
        MechNode *p = dynamic_cast< MechNode * >( nodes [ i ] );
        if ( !p ) {
            cerr << "Error in " << name << ": nodes must be inherited from MechNode, " << nodes [ i ]->giveName() << " provided" << endl;
            exit(1);
        }
    }
}

//////////////////////////////////////////////////////////
Matrix Truss :: giveBMatrix() const {
    //Matrix B
    Matrix B = Matrix(ndim, 2 * ndim);
    Matrix Aa = giveAMatrix(nodes [ 0 ]->givePoint(), ip_locs [ 0 ]) * ( -1. );
    Matrix Ab = giveAMatrix(nodes [ 1 ]->givePoint(), ip_locs [ 0 ]);
    for ( unsigned i = 0; i < ndim; i++ ) {
        for ( unsigned j = 0; j < ndim; j++ ) {
            B [ i ] [ j ] = Aa [ i ] [ j ];
            B [ i ] [ j + ndim ] = Ab [ i ] [ j ];
        }
    }
    return B / length;
}

//////////////////////////////////////////////////////////
Vector Truss :: giveContactStrainNT(const Vector &DoFs) const {
    Vector strain = GeomM * DoFs;
    for ( size_t k = 1; k < strain.size(); k++ ) {
        strain [ k ] = 0;                                  //only normal strain active in truss
    }
    return strain;
};

//////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////
// 1D TRANSPORT ELEMENT
Transp1D :: Transp1D(const unsigned dim) {
    ndim = dim;
    nodes.resize(2);
    ip_locs.resize(1);
    stats.resize(1);
    bound = false;
    name = "Transp1D";
    reducedCapacityMatrix = false;
}

//////////////////////////////////////////////////////////
void Transp1D :: readFromLine(istringstream &iss, NodeContainer *fullnodes, MaterialContainer *fullmatrs) {
    unsigned num, num2;

    iss >> num;
    nodes [ 0 ] = fullnodes->giveNode(num);

    iss >> num;
    nodes [ 1 ] = fullnodes->giveNode(num);

    iss >> num;
    vert.resize(num);
    for ( unsigned i = 0; i < num; i++ ) {
        iss >> num2;
        vert [ i ] = fullnodes->giveNode(num2);
    }
    iss >> num;
    mat = fullmatrs->giveMaterial(num);

    string code;
    while ( iss >> code ) {
        if ( code.compare("reducedCapacityMatrix") == 0 ) {
            reducedCapacityMatrix = true;
        }
    }

    //  cout<< "Loaded 1D trsprt: "<<nodes.size()<<" nodes, "<<vert.size()<<" vertices"<<endl;
}

//////////////////////////////////////////////////////////
void Transp1D :: init() {
    Element :: init(); //calling base class method;
    //check that nodes are TrsNodes
    for ( unsigned i = 0; i < 2; i++ ) {
        TrsNode *p = dynamic_cast< TrsNode * >( nodes [ i ] );
        if ( !p ) {
            cerr << "Error in " << name << ": nodes must be inherited from TrsNode, " << nodes [ i ]->giveName() << " provided" << endl;
            exit(1);
        }
    }
    //check that material is DisMechMat
    TrsprtMaterial *p = dynamic_cast< TrsprtMaterial * >( mat );
    if ( !p ) {
        cerr << "Error in " << name << ": material must be inherited from TrsprtMaterial, " << mat->giveName() << " provided" << endl;
        exit(1);
    }

    Point t;
    if ( ndim == 2 ) {
        if ( !( vert.size() == 2 ) ) {
            cerr << "TRSPRT: exactly 2 vertices must be involved, " << vert.size() << " provided" << endl;
            exit(1);
        }

        ip_locs [ 0 ] = ( vert [ 0 ]->givePoint() + vert [ 1 ]->givePoint() ) / 2.;
        t = vert [ 1 ]->givePoint() - vert [ 0 ]->givePoint();
        area = t.norm();
        t = t / area;
    } else {
        //JM: Coplanarity check for vertices on the face
        //JM: checking coplanarity of every consecutive 4 nodes
        double maxErr = 0.0;
        double currErr = 0.0;
        //
        if ( vert.size() > 3 ) {
            for ( unsigned int i = 0; i < vert.size() - 3; i++ ) {
                // JM Zakomentoval cout << i <<  " " << endl;
                currErr = checkCoplanarity(vert [ i ]->givePoint(), vert [ i + 1 ]->givePoint(), vert [ i + 2 ]->givePoint(), vert [ i + 3 ]->givePoint() );
                if ( abs(currErr) > maxErr ) {
                    maxErr = abs(currErr);
                }
            }
            // JM Zakomentoval  cout << maxErr << endl;
        }
        if ( maxErr > 1e-10 ) {
            cerr << "TRSPRT: Vertices are not coplanar!!! Coplanarity error: " << maxErr << endl;
            exit(1);
        }

        //JM: face normal vector made from first 3 vertices
        //JM: coordinate swap for tangential vector according to https://orbit.dtu.dk/files/126824972/onb_frisvad_jgt2012_v2.pdf
        Point n = cross(vert [ 1 ]->givePoint() - vert [ 0 ]->givePoint(), vert [ 2 ]->givePoint() - vert [ 0 ]->givePoint() );
        n /= n.norm();
        Point t2;
        if ( fabs(n.x) > fabs(n.z) ) {
            t2 = Point(-n.y, n.x, 0.0f);
        } else {
            t2 = Point(0.0f, -n.z, n.y);
        }
        t = cross(t2, n);
        t /= t.norm();

        //JM: Perpendicularity check of the beam and face directions
        //JM: normal of the face surface taken from first 3 vertices is (B - A) x (C - A)
        //JM: perpendicularity check: cross (beam, face)=>0
        Point prp = ( nodes [ 1 ]->givePoint() - nodes [ 0 ]->givePoint() ) * t;
        if ( prp.norm() > 1e-10 ) {
            cerr << "TRSPRT: Face surface is not perpendicular to beam direction!!! Error: " << prp.norm() << endl;
            //  exit(1);
        }

        //JM: finding position of the SINGLE integration point -> center of gravity of the face polygon
        //JM: average point of the polygon for triangulation
        Point avgPoint = Point(0.0, 0.0, 0.0);
        for ( unsigned int i = 0; i < vert.size(); i++ ) {
            avgPoint += vert [ i ]->givePoint();
        }
        avgPoint /= vert.size();

        //JM: integration point coordinates as an average of CGs of face triangles weighted by areas
        ip_locs [ 0 ] = Point(0.0, 0.0, 0.0);
        area = 0.0;
        double ai = 0.0;
        unsigned int j = 0;
        for ( unsigned int i = 0; i < vert.size(); i++ ) {
            j = i + 1;
            if ( i == vert.size() - 1 ) {
                j = 0;
            }
            //triangle area computed as a_i = norm(cross(AB, AC)) / 2
            ai = ( cross(vert [ i ]->givePoint() - avgPoint,   vert [ j ]->givePoint() - avgPoint) ).norm();
            area += ai;
            //triangle cg_i is an average of simplex vertices, adding to CG coordinates multiplied by a_i weight
            ip_locs [ 0 ] += ( avgPoint + vert [ i ]->givePoint() + vert [ j ]->givePoint() ) / 3.0 * ai;
        }
        ip_locs [ 0 ] /= area;

        //JM: Check if integration point is coplanar with face
        currErr = checkCoplanarity(vert [ 0 ]->givePoint(), vert [ 1 ]->givePoint(), vert [ 2 ]->givePoint(), ip_locs [ 0 ]);
        if ( abs(currErr) > 1e-10 ) {
            cerr << "TRSPRT: Integration point is not coplanar with the face!!! Coplanarity error: " << currErr << endl;
            exit(1);
        }
    }

    stats [ 0 ] = mat->giveNewMaterialStatus(this);
    normal = nodes [ 1 ]->givePoint() - nodes [ 0 ]->givePoint();
    length = normal.norm();

    normal = normal / length;
    if ( abs(normal * t) > 1e-5 ) {
        cout << vert [ 0 ]->givePoint().x << " " <<  vert [ 0 ]->givePoint().y <<  " X " << vert [ 1 ]->givePoint().x << " " <<  vert [ 1 ]->givePoint().y << endl;
        cout << nodes [ 0 ]->givePoint().x << " " <<  nodes [ 0 ]->givePoint().y <<  " X " << nodes [ 1 ]->givePoint().x << " " <<  nodes [ 1 ]->givePoint().y << endl;
        cerr << "TRSPRT: normal and contact vector are not parallel, error " << normal * t << " normal v." << normal.x << " " << normal.y << " contact v. " << t.x << " " << t.y << endl;
        exit(1);
    }
}

//////////////////////////////////////////////////////////
Matrix Transp1D :: giveConductivityMatrix(string matrixType) const {
    ( void ) matrixType;
    TrsprtMaterialStatus *tstat = static_cast< TrsprtMaterialStatus * >( stats [ 0 ] );
    double c = area * tstat->giveEffectiveConductivity(matrixType) / length;
    Matrix C(2, 2);
    C [ 0 ] [ 0 ] = C [ 1 ] [ 1 ] = c;
    C [ 1 ] [ 0 ] = C [ 0 ] [ 1 ] = -c;
    return C;
}

//////////////////////////////////////////////////////////
Matrix Transp1D :: giveCapacityMatrix() const {
    Matrix S(2, 2);
    TrsprtMaterial *tmat = static_cast< TrsprtMaterial * >( mat );
    double s = area * tmat->giveCapacity() * tmat->giveDensity() * length / ( 12. );

    if ( reducedCapacityMatrix ) {    //derived based on assumption of constant pressure
        S [ 0 ] [ 0 ] = S [ 1 ] [ 1 ] = 3 * s;
    } else {    //derived based on linear assumption of pressure
        S [ 0 ] [ 0 ] = S [ 1 ] [ 1 ] = 2 * s;
        S [ 1 ] [ 0 ] = S [ 0 ] [ 1 ] = s;
    }
    return S;
}

//////////////////////////////////////////////////////////
Vector Transp1D :: giveInternalForces(const Vector &DoFs, bool frozen) const {
    ( void ) frozen;
    Vector pressureGrad(1);
    pressureGrad [ 0 ] = ( DoFs [ 1 ] - DoFs [ 0 ] ) / length;
    Vector flux = stats [ 0 ]->giveStress(pressureGrad);
    Vector intf(2);
    intf [ 0 ] = flux [ 0 ] * area;
    intf [ 1 ] = -intf [ 0 ];
    return intf;
};

//////////////////////////////////////////////////////////
double Transp1D :: giveVolume() const {
    return area * length / ndim;
};

//////////////////////////////////////////////////////////
double Transp1D :: giveVolume(unsigned nodenum) const {
    if ( nodenum == 0 ) {
        return dot(vert [ 0 ]->givePoint() - nodes [ 0 ]->givePoint(), normal) * area / ndim;
    } else if ( nodenum == 1 )        {
        return -dot(vert [ 0 ]->givePoint() - nodes [ 1 ]->givePoint(), normal) * area / ndim;
    } else                                                                                                                       {
        cerr << "Error in " << name << ": attempting to reach node number different form 0 or 1." << endl;
        exit(1);
    }
};

//////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////
// 1D TRANSPORT ELEMENT COUPLED WITH MECHANICS
//////////////////////////////////////////////////////////
void Transp1DCoupled :: init() {
    Transp1D :: init(); //calling base class method;
}

//////////////////////////////////////////////////////////
void Transp1DCoupled :: findFriends2D(ElementContainer *elemcont) {
    for ( unsigned k = 0; k < elemcont->giveSize(); k++ ) {
        RigidBodyContact *rbc = dynamic_cast< RigidBodyContact * >( elemcont->giveElement(k) );
        if ( rbc ) {
            if ( ( rbc->giveNode(0) == vert [ 0 ] && rbc->giveNode(1) == vert [ 1 ] ) || ( rbc->giveNode(0) == vert [ 1 ] && rbc->giveNode(1) == vert [ 0 ] ) ) {
                friends.push_back(rbc);
                friendsweight.push_back(rbc->giveArea() );
                break; //only one friend avalilable in 2D
            }
        }
    }
}

//////////////////////////////////////////////////////////
void Transp1DCoupled :: findFriends3D(ElementContainer *elemcont) {
    for ( vector< Node * > :: const_iterator b = vert.begin(); b != vert.end(); ++b ) {}
    //for(unsigned k)
}

//////////////////////////////////////////////////////////
void Transp1DCoupled :: findElementFriends(ElementContainer *elemcont) {
    if ( ndim == 2 ) {
        findFriends2D(elemcont);
    } else if ( ndim == 3 )   {
        findFriends3D(elemcont);
    }

    //collect nodes from friend elements
}

//////////////////////////////////////////////////////////
Vector Transp1DCoupled :: giveInternalForces(const Vector &DoFs, bool frozen) const {
    if ( frozen ) { //frozen internal variables
        return giveConductivityMatrix("secant") * DoFs;
    } else  {
        Vector pressureGrad(2 + 2 * friends.size() );
        pressureGrad [ 0 ] = ( DoFs [ 1 ] - DoFs [ 0 ] ) / length; //the real pressure gradient
        pressureGrad [ 1 ] = area; //area of the element face
        double elem_crack_opening;
        size_t m = 0;
        for ( auto &f: friends ) {
            elem_crack_opening = 0;
            for ( unsigned k = 0; k < f->giveIPNum(); k++ ) {
                elem_crack_opening += abs(f->giveIPValue("tempCrackOpening", k) );
            }
            pressureGrad [ 2 * m + 2 ] += elem_crack_opening / f->giveIPNum(); //average crack opening in friend mechanical element
            pressureGrad [ 2 * m + 3 ] = friendsweight [ m ]; //crack length in friend mechanical element
            m++;
        }
        Vector flux = stats [ 0 ]->giveStress(pressureGrad);
        Vector intf(2);
        intf [ 0 ] = flux [ 0 ] * area;
        intf [ 1 ] = -intf [ 0 ];
        return intf;
    }
};

//////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////
// 2D QUADRILATERAL TRANSPORT ELEMENT
TranspQuad :: TranspQuad(){
    ndim = 2;
    name = "Transport Quadrilateral";
}

//////////////////////////////////////////////////////////
void TranspQuad :: init(){    
    TransportElement :: init();
    ip_locs.resize(4);
    ip_weights.resize(4);
    stats.resize(4);
    double q = 1./pow(3.,0.5);
    ip_locs[0] = Point(-q,-q);
    ip_locs[1] = Point( q,-q);
    ip_locs[2] = Point( q, q);
    ip_locs[3] = Point(-q, q);
    Matrix phiGrad(2,4);
    for(unsigned k=0; k<4; k++){
        stats[k] = mat->giveNewMaterialStatus(this);
        ip_weights[k] = shapeFGrad(&(ip_locs[k]), phiGrad);
    }    
}

//////////////////////////////////////////////////////////
void TranspQuad :: readFromLine(istringstream &iss, NodeContainer *fullnodes, MaterialContainer *fullmatrs){
    unsigned num;
    nodes.resize(4);
    for(unsigned k=0; k<4; k++){
        iss >> num;
        nodes [ k ] = fullnodes->giveNode(num);
    }
    iss >> num;
    mat = fullmatrs->giveMaterial(num);    
}

//////////////////////////////////////////////////////////
void TranspQuad :: shapeF(const Point *x, Vector &phi) const{ 
    //x in natural coordinates
    phi[0] = 0.25*(1.+x->getX())*(1.+x->getY());
    phi[1] = 0.25*(1.-x->getX())*(1.+x->getY());
    phi[2] = 0.25*(1.-x->getX())*(1.-x->getY());
    phi[3] = 0.25*(1.+x->getX())*(1.-x->getY());
}

//////////////////////////////////////////////////////////
double TranspQuad :: shapeFGrad(const Point *x, Matrix &phiGrad) const{
    //x in natural coordinates
    phiGrad[0][0] = 0.25*(1.+x->getY());
    phiGrad[0][1] = -phiGrad[0][0];
    phiGrad[0][2] = -0.25*(1.-x->getY());
    phiGrad[0][3] = -phiGrad[0][2];
    phiGrad[1][0] = 0.25*(1.+x->getX());
    phiGrad[1][1] = 0.25*(1.-x->getX());
    phiGrad[1][2] = -phiGrad[1][1];
    phiGrad[1][3] = -phiGrad[1][0];

    //Jacobi Matrix
    Matrix Jac(2,2);
    Point n;
    for(unsigned i=0; i<4; i++){
        n = nodes[i]->givePoint();
        Jac[0][0] += phiGrad[0][i]*n.getX();
        Jac[0][1] += phiGrad[0][i]*n.getY();
        Jac[1][0] += phiGrad[1][i]*n.getX();
        Jac[1][1] += phiGrad[1][i]*n.getY();
    }
    double JacDet = Jac[0][0]*Jac[1][1] - Jac[0][1]*Jac[1][0];
    Matrix JacInv(2,2);
    JacInv[0][0] = Jac[1][1]/JacDet;
    JacInv[1][0] = -Jac[0][1]/JacDet;
    JacInv[0][1] = -Jac[1][0]/JacDet;
    JacInv[1][1] = Jac[0][0]/JacDet;
    
    //transorm to physical space
    phiGrad = JacInv*phiGrad;
    return JacDet;
}

//////////////////////////////////////////////////////////
Matrix TranspQuad :: giveConductivityMatrix(string matrixType) const {
    unsigned nDoFs = DoFids.size();
    Matrix K(nDoFs,nDoFs);
    Matrix phiGrad(ndim,nDoFs);
    Matrix D;
    for(unsigned i=0; i<stats.size(); i++){
        shapeFGrad(&ip_locs[i], phiGrad);
        D = stats [ i ] -> giveStiffnessTensor(matrixType, ndim);
        K += phiGrad.transpose() * D * (phiGrad * ip_weights[i]);
    }
    return K;
}

//////////////////////////////////////////////////////////
Matrix TranspQuad :: giveCapacityMatrix() const {
    unsigned nDoFs = DoFids.size();
    Matrix M(nDoFs,nDoFs);
    Vector phi(nDoFs);
    double c;
    for(unsigned i=0; i<stats.size(); i++){
        shapeF(&ip_locs[i], phi);
        c = stats [ i ] -> giveMassConstant();
        M += dyadicProduct(phi,phi)*(ip_weights[i]*c);
    }
    return M;
}

//////////////////////////////////////////////////////////
Vector TranspQuad :: giveInternalForces(const Vector &DoFs, bool frozen) const {
    if ( frozen ) { //frozen internal variables
        return giveConductivityMatrix("secant") * DoFs;
    } else  {    //evalution of internal variables
        Vector intF;
        Matrix B(ndim,DoFids.size());
        intF.resize(DoFids.size());
        for(unsigned i=0; i<stats.size(); i++){
            shapeFGrad(&ip_locs[i],B);         
            intF  -= B.transpose() * (stats[i]->giveStress(B*DoFs) * ip_weights[i]);
        }
        return intF;
    }
};
