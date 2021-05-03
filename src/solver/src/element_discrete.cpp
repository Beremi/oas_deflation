#include "element_discrete.h"
#include "element_container.h"
#include "boundary_condition.h"

//////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////
// RBSN ELEMENT
RigidBodyContact :: RigidBodyContact(const unsigned dim) {
    ndim = dim;
    numOfNodes = 2;
    nodes.resize(2);
    name = "LTCBEAM";
    vtk_cell_type = 3;
    shafunc = new Linear1DLineShapeF();
    inttype = new IntegrDiscrete1();
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

Matrix RigidBodyContact :: giveBMatrix(const Point *x) const {
    //Matrix B
    Matrix B = Matrix(ndim, 6 * ( ndim - 1 ) );
    Matrix Aa = giveAMatrix(nodes [ 0 ]->givePoint(), inttype->giveIPLocation(0) ) * ( -1. );
    Matrix Ab = giveAMatrix(nodes [ 1 ]->givePoint(), inttype->giveIPLocation(0) );
    for ( unsigned i = 0; i < ndim; i++ ) {
        for ( unsigned j = 0; j < 3 * ( ndim - 1 ); j++ ) {
            B [ i ] [ j ] = Aa [ i ] [ j ];
            B [ i ] [ j + 3 * ( ndim - 1 ) ] = Ab [ i ] [ j ];
        }
    }
    return matrix_multiply(R, B) / length;
}


//////////////////////////////////////////////////////////
void RigidBodyContact :: setIntegrationPointsAndWeights() {
    stats.resize(1);

    Point t;
    if ( ndim == 2 ) {
        // NOTE the following taken from few lines upper
        if ( !( vert.size() == 2 ) ) {
            cerr << "Error: exactly 2 vertices must be involved, " << vert.size() << " provided" << endl;
            exit(1);
        }
        inttype->setIPLocation(0, ( vert [ 0 ]->givePoint() + vert [ 1 ]->givePoint() ) / 2.);
        /////////////////////////////////////////////////////////
        t = vert [ 1 ]->givePoint() - vert [ 0 ]->givePoint();
        area = t.norm();
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
        if ( maxErr > 1e-6 ) {
            cerr << "Vertices are not coplanar!!! Coplanarity error: " << maxErr << endl;
            exit(1);
        }

        //JM: face normal vector made from first 3 vertices
        Point n = cross(vert [ 1 ]->givePoint() - vert [ 0 ]->givePoint(), vert [ 2 ]->givePoint() - vert [ 0 ]->givePoint() );
        n /= n.norm();

        //JM: Perpendicularity check of the beam and face directions
        //JM: normal of the face surface taken from first 3 vertices is (B - A) x (C - A)
        //JM: perpendicularity check: cross (beam, face)=>0
        Point prp = ( nodes [ 1 ]->givePoint() - nodes [ 0 ]->givePoint() ) * t;
        if ( prp.norm() > 1e-8 ) {
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
        Point centroid = Point(0.0, 0.0, 0.0);
        area = 0.0;
        double ai = 0.0;
        unsigned int j = 0;
        for ( unsigned int i = 0; i < vert.size(); i++ ) {
            j = i + 1;
            if ( i == vert.size() - 1 ) {
                j = 0;
            }
            //triangle area computed as a_i = norm(cross(AB, AC)) / 2
            ai = ( cross(vert [ i ]->givePoint() - avgPoint,   vert [ j ]->givePoint() - avgPoint) ).norm() / 2.;
            area += ai;
            //triangle cg_i is an average of simplex vertices, adding to CG coordinates multiplied by a_i weight
            centroid += ( avgPoint + vert [ i ]->givePoint() + vert [ j ]->givePoint() ) / 3.0 * ai;
        }
        centroid /= area;
        inttype->setIPLocation(0, centroid);

        //JM: Check if integration point is coplanar with face
        currErr = checkCoplanarity(vert [ 0 ]->givePoint(), vert [ 1 ]->givePoint(), vert [ 2 ]->givePoint(), inttype->giveIPLocation(0) );
        if ( abs(currErr) > 1e-6 ) {
            cerr << "Integration point is not coplanar with the face!!! Coplanarity error: " << currErr << endl;
            exit(1);
        }
    }

    normal = nodes [ 1 ]->givePoint() - nodes [ 0 ]->givePoint();
    length = normal.norm();
    normal = normal / length;
    if ( abs(normal * t) > 1e-8 ) {
        cout << vert [ 0 ]->givePoint().x << " " <<  vert [ 0 ]->givePoint().y <<  " X " << vert [ 1 ]->givePoint().x << " " <<  vert [ 1 ]->givePoint().y << endl;
        cout << nodes [ 0 ]->givePoint().x << " " <<  nodes [ 0 ]->givePoint().y <<  " X " << nodes [ 1 ]->givePoint().x << " " <<  nodes [ 1 ]->givePoint().y << endl;
        cerr << "Error: normal and contact vector are not parallel, error " << normal * t << " normal v." << normal.x << " " << normal.y << " contact v. " << t.x << " " << t.y << endl;
        // exit(1);
    }

    // Matrices according to habilitation of Jan Elias (2017, page 42): https://www.vutbr.cz/www_base/vutdisk.php?i=103116a130
    // Matrix R;
    if ( ndim == 2 ) {
        t1 = Point(-normal.y, normal.x);
        R = Matrix(2, 2);
        R [ 0 ] [ 0 ] = normal.x;
        R [ 0 ] [ 1 ] = normal.y;
        R [ 1 ] [ 0 ] = t1.x;
        R [ 1 ] [ 1 ] = t1.y;
    } else if ( ndim == 3 ) {
        // coordinate swap for tangential vector according to https://orbit.dtu.dk/files/126824972/onb_frisvad_jgt2012_v2.pdf
        Point arbit(sqrt(2.), -sqrt(3.), M_PI);
        if ( ( normal - arbit ).norm() < 1e-3 ) {
            t1 = cross(arbit, normal);
            t1.normalize();
            t2 = cross(normal, t1);
            t2.normalize();
        } else {
            // the following results in zeros in stiffness matrix in case of normal in direction of any of global base axes
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

    inttype->setIPWeight(0, length * area);  // NOTE JK: this works for stiffness matrix, since there is actually A * L, but not for the mass matrix, where volume should be taken (divide by dim) // JE: Agreed, independent integration of mass has been included for this element

    stats [ 0 ] = mat->giveNewMaterialStatus(this, 0);
}

//////////////////////////////////////////////////////////
void RigidBodyContact :: init() {
    Element :: init(); //calling base class method;

    checkNodeType();

    //create simplices
    for ( auto &v: vert ) {
        simplices.push_back(v->addElementToSimplex(this) );
    }

    //check that material is DisMechMat
    DisMechMaterial *p = dynamic_cast< DisMechMaterial * >( mat );
    if ( !p ) {
        cerr << "Error in " << name << ": material must be inherited from DisMechMaterial, " << mat->giveName() << " provided" << endl;
        exit(1);
    }
}

//////////////////////////////////////////////////////////
Matrix RigidBodyContact :: giveHMatrix(const Point *x) const {
    return Matrix(12, 12);
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
Vector RigidBodyContact :: giveContactStrainNT() const {
    return stats [ 0 ]->giveTempStrain();
};

//////////////////////////////////////////////////////////
Vector RigidBodyContact :: giveContactStrainXYZ() const {
    return this->R.transpose() * this->giveContactStrainNT();
};

//////////////////////////////////////////////////////////
Vector RigidBodyContact :: giveContactStressXYZ() {
    return this->R.transpose() * stats [ 0 ]->giveTempStress();
};

//////////////////////////////////////////////////////////
Vector RigidBodyContact :: transformToLocal(const Vector &DoFs) const {
    return this->R.transpose() * DoFs;
}

//////////////////////////////////////////////////////////
Vector RigidBodyContact :: transformToGlobal(const Vector &DoFs) const {
    return this->R * DoFs;
}

//////////////////////////////////////////////////////////
Vector RigidBodyContact :: giveVectorToNode(const unsigned &node_i, const unsigned &ip_id) const {
    Point distance = inttype->giveIPLocation(0) - nodes [ node_i ]->givePoint();
    Vector dst( ( double ) 0, ndim );
    for ( unsigned i = 0; i < ndim; i++ ) {
        if ( i == 0 ) {
            dst [ i ] = distance.getX();
        } else if ( i == 1 ) {
            dst [ i ] = distance.getY();
        } else if ( i == 2 ) {
            dst [ i ] = distance.getZ();
        }
    }
    return dst;
}

//////////////////////////////////////////////////////////
double RigidBodyContact :: giveVolume() const {
    return area * length / ndim;
};

//////////////////////////////////////////////////////////
double RigidBodyContact :: giveVolume(unsigned nodenum) const {
    if ( nodenum == 0 ) {
        return dot(vert [ 0 ]->givePoint() - nodes [ 0 ]->givePoint(), normal) * area / ndim;
    } else if ( nodenum == 1 ) {
        return -dot(vert [ 0 ]->givePoint() - nodes [ 1 ]->givePoint(), normal) * area / ndim;
    } else {
        cerr << "Error in " << name << ": attempting to reach node number different form 0 or 1." << endl;
        exit(1);
    }
};

//////////////////////////////////////////////////////////
Vector RigidBodyContact :: giveStrain(unsigned i, const Vector &DoFs) {
    //extract volumetric strains from simplices
    volumetricStrain = 0;
    unsigned validSnum = 0;
    for ( auto &s: simplices ) {
        if ( s->isValid() ) {
            validSnum++;
            volumetricStrain += s->giveVolumetricStrain();
        }
    }
    if ( validSnum > 0 ) {
        volumetricStrain /= validSnum;
    }

    return Element :: giveStrain(i, DoFs);
};

//////////////////////////////////////////////////////////
Matrix RigidBodyContact :: giveDampingMatrix() const {
    return giveStiffnessMatrix("elastic") * 1e-15;           //rough fix of zeros, here can be anything
}

//////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////
// COUPLED RBSN ELEMENT
RigidBodyContactCoupled :: RigidBodyContactCoupled(const unsigned dim) : RigidBodyContact(dim) {
    name = "LTCBEAMCoupled";
}

//////////////////////////////////////////////////////////
Vector RigidBodyContactCoupled :: giveStrain(unsigned i, const Vector &DoFs) {
    //extract pressure from simplices
    averagePressure = 0;
    unsigned validSnum = 0;
    for ( auto &s: simplices ) {
        if ( s->hasPressure() ) {
            validSnum++;
            averagePressure += s->givePressure();
        }
    }
    if ( validSnum > 0 ) {
        averagePressure /= validSnum;
    }

    return RigidBodyContact :: giveStrain(i, DoFs);
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
Matrix Truss :: giveBMatrix(const Point *x) const {
    //Matrix B
    Matrix B = Matrix(ndim, 2 * ndim);
    Matrix Aa = giveAMatrix(nodes [ 0 ]->givePoint(), inttype->giveIPLocation(0) ) * ( -1. );
    Matrix Ab = giveAMatrix(nodes [ 1 ]->givePoint(), inttype->giveIPLocation(0) );
    for ( unsigned i = 0; i < ndim; i++ ) {
        for ( unsigned j = 0; j < ndim; j++ ) {
            B [ i ] [ j ] = Aa [ i ] [ j ];
            B [ i ] [ j + ndim ] = Ab [ i ] [ j ];
        }
    }
    return matrix_multiply(R, B) / length;
}

//////////////////////////////////////////////////////////
Matrix Truss :: giveHMatrix(const Point *x) const {
    return Matrix(0, 0);
}
//////////////////////////////////////////////////////////
Vector Truss :: giveContactStrainNT(const Vector &DoFs) const {
    Vector strain = Bs [ 0 ] * DoFs;
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
    bound = false;
    name = "LTCTRSP";
    BolanderCapacityMatrix = false;
    shafunc = new Linear1DLineShapeF();
    inttype = new IntegrDiscrete1();
    vtk_cell_type = 3;
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
        if ( code.compare("BolanderCapacityMatrix") == 0 ) {
            BolanderCapacityMatrix = true;
        }
    }

    //  cout<< "Loaded 1D trsprt: "<<nodes.size()<<" nodes, "<<vert.size()<<" vertices"<<endl;
}



//////////////////////////////////////////////////////////
void Transp1D :: setIntegrationPointsAndWeights() {
    stats.resize(1);

    normal = nodes [ 1 ]->givePoint() - nodes [ 0 ]->givePoint();
    length = normal.norm();
    normal = normal / length;
    if ( length < 1e-8 ) {
        length = 1e-8;           //artificial increase of length in case of extremely short voronoi edge
    }
    Point t;
    if ( ndim == 2 ) {
        if ( !( vert.size() == 2 ) ) {
            cerr << "Error: exactly 2 vertices must be involved, " << vert.size() << " provided" << endl;
            exit(1);
        }

        inttype->setIPLocation(0, ( vert [ 0 ]->givePoint() + vert [ 1 ]->givePoint() ) / 2.);
        t = vert [ 1 ]->givePoint() - vert [ 0 ]->givePoint();
        area = t.norm();
        t = t / area;
    } else {
        //JM: Coplanarity check for vertices on the face
        //JM: checking coplanarity of every consecutive 4 nodes
        double maxErr = 0.0;
        double currErr = 0.0;
        //

        if (  vert.size() < 3  ) {
            cerr << name << " Error: three or more vertices are required, only " << vert.size() << " provided" << endl;
            exit(1);
        }

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
        if ( prp.norm() > 1e-8 ) {
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
        Point centroid = Point(0.0, 0.0, 0.0);
        area = 0.0;
        double ai = 0.0;
        unsigned int j = 0;
        for ( unsigned int i = 0; i < vert.size(); i++ ) {
            j = i + 1;
            if ( i == vert.size() - 1 ) {
                j = 0;
            }
            //triangle area computed as a_i = norm(cross(AB, AC)) / 2
            ai = ( cross(vert [ i ]->givePoint() - avgPoint,   vert [ j ]->givePoint() - avgPoint) ).norm() / 2.;
            area += ai;
            //triangle cg_i is an average of simplex vertices, adding to CG coordinates multiplied by a_i weight
            centroid += ( avgPoint + vert [ i ]->givePoint() + vert [ j ]->givePoint() ) / 3.0 * ai;
        }
        centroid /= area;
        inttype->setIPLocation(0, centroid);

        //JM: Check if integration point is coplanar with face
        currErr = checkCoplanarity(vert [ 0 ]->givePoint(), vert [ 1 ]->givePoint(), vert [ 2 ]->givePoint(), inttype->giveIPLocation(0) );
        if ( abs(currErr) > 1e-10 ) {
            cerr << "TRSPRT: Integration point is not coplanar with the face!!! Coplanarity error: " << currErr << endl;
            exit(1);
        }
    }

    if ( abs(normal * t) > 1e-5 ) {
        cout << vert [ 0 ]->givePoint().x << " " <<  vert [ 0 ]->givePoint().y <<  " X " << vert [ 1 ]->givePoint().x << " " <<  vert [ 1 ]->givePoint().y << endl;
        cout << nodes [ 0 ]->givePoint().x << " " <<  nodes [ 0 ]->givePoint().y <<  " X " << nodes [ 1 ]->givePoint().x << " " <<  nodes [ 1 ]->givePoint().y << endl;
        cerr << "TRSPRT: normal and contact vector are not parallel, error " << normal * t << endl;
        cout << " normal v.:";
        for ( unsigned p = 0; p < ndim; p++ ) {
            cout << "\t" << normal.giveCoord(p);
        }
        cout << endl;
        cout << " contact v.:";
        for ( unsigned p = 0; p < ndim; p++ ) {
            cout << "\t" << t.giveCoord(p);
        }
        cout << endl;
        //exit(1);
    }

    inttype->setIPWeight(0, length * area);   // NOTE JK: should not this be divided by dimension? Otherwise you integrate dim-times volume that is actually there (this is what caused problems in study on unstructured grid for FDM contribution)
    stats [ 0 ] = mat->giveNewMaterialStatus(this, 0);
}

//////////////////////////////////////////////////////////
void Transp1D :: init() {
    Element :: init(); //calling base class method;

    checkNodeType();

    //check that material is TrsprtMaterial
    TrsprtMaterial *p = dynamic_cast< TrsprtMaterial * >( mat );
    if ( !p ) {
        cerr << "Error in " << name << ": material must be inherited from TrsprtMaterial, " << mat->giveName() << " provided" << endl;
        exit(1);
    }
}

//////////////////////////////////////////////////////////
void Transp1D :: checkNodeType() const {
    for ( unsigned i = 0; i < 2; i++ ) {
        TrsNode *p = dynamic_cast< TrsNode * >( nodes [ i ] );
        if ( !p ) {
            cerr << "Error in " << name << ": nodes must be inherited from TrsNode, " << nodes [ i ]->giveName() << " provided" << endl;
            exit(1);
        }
    }
}

//////////////////////////////////////////////////////////
Matrix Transp1D :: giveBMatrix(const Point *x) const {
    ( void ) x;
    Matrix B(1, 2);
    B [ 0 ] [ 0 ] = -1. / length;
    B [ 0 ] [ 1 ] = 1. / length;
    return B;
}

//////////////////////////////////////////////////////////
Matrix Transp1D :: giveHMatrix(const Point *x) const {
    Matrix H(1, 2);
    double l1 = dot(* x - nodes [ 0 ]->givePoint(), normal);
    double l2 = dot(nodes [ 1 ]->givePoint() - * x, normal);
    H [ 0 ] [ 0 ] = l1 / length;
    H [ 0 ] [ 1 ] = l2 / length;
    return H;
}

//////////////////////////////////////////////////////////
Matrix Transp1D :: giveDampingMatrix() const {
    Matrix S(2, 2);
    double s = area * stats [ 0 ]->giveDampingConstant() * length /  ( 2. * ndim );

    S [ 0 ] [ 0 ] = S [ 1 ] [ 1 ] = s; //finite volume
    if ( BolanderCapacityMatrix ) { //from Bolander's papers
        S [ 0 ] [ 0 ] = S [ 1 ] [ 1 ] = 2. / 3. * s;
        S [ 1 ] [ 0 ] = S [ 0 ] [ 1 ] = s / 3.;
    }
    return S;
}


//////////////////////////////////////////////////////////
vector< double >Transp1D :: integrateLoad(BodyLoad *vl, double time) const {
    vector< double >load = Element :: integrateLoad(vl, time);
    for ( auto &p:load ) {
        p /= ndim;
    }
    return load;
}

//////////////////////////////////////////////////////////
double Transp1D :: giveVolume() const {
    return area * length / ndim;
};

//////////////////////////////////////////////////////////
double Transp1D :: giveVolume(unsigned nodenum) const {
    if ( nodenum == 0 ) {
        return dot(vert [ 0 ]->givePoint() - nodes [ 0 ]->givePoint(), normal) * area / ndim;
    } else if ( nodenum == 1 ) {
        return -dot(vert [ 0 ]->givePoint() - nodes [ 1 ]->givePoint(), normal) * area / ndim;
    } else {
        cerr << "Error in " << name << ": attempting to reach node number different form 0 or 1." << endl;
        exit(1);
    }
};

//////////////////////////////////////////////////////////
Vector Transp1D :: giveStrain(unsigned i, const Vector &DoFs) {
    averagePressure = ( DoFs [ 0 ] + DoFs [ 1 ] ) / 2.;
    return Element :: giveStrain(i, DoFs);
};


//////////////////////////////////////////////////////////
double Transp1D :: giveAveragePressure() {
    return averagePressure;
}


//////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////
// 1D TRANSPORT ELEMENT COUPLED WITH MECHANICS
//////////////////////////////////////////////////////////
void Transp1DCoupled :: init() {
    Transp1D :: init(); //calling base class method;
}

//////////////////////////////////////////////////////////
double Transp1DCoupled :: giveValue(string code) const {
    if ( code.compare("numOfFriends") == 0 ) {
        return friends.size();
    } else {
        return Transp1D :: giveValue(code);
    }
};

//////////////////////////////////////////////////////////
double Transp1DCoupled :: giveIPValue(string code, unsigned ipnum) const {
    if ( code.compare("numOfFriends") == 0 ) {
        return friends.size();
    } else {
        return Transp1D :: giveIPValue(code, ipnum);
    }
}

//////////////////////////////////////////////////////////
void Transp1DCoupled :: addNewFriend(RigidBodyContact *f, double weight) {
    friends.push_back(f);
    friendsweight.push_back(weight);
}

//////////////////////////////////////////////////////////
double Transp1DCoupled :: giveCrackOpeningInNeigborhood() {
    return crackInNeighborhood;
}

//////////////////////////////////////////////////////////
Vector Transp1DCoupled :: giveStrain(unsigned i, const Vector &DoFs) {
    crackInNeighborhood = 0;
    double elem_crack_opening;
    size_t m = 0;
    for ( auto &f: friends ) {
        elem_crack_opening = 0;
        for ( unsigned k = 0; k < f->giveNumIP(); k++ ) {
            elem_crack_opening += abs(f->giveIPValue("tempCrackOpening", k) );
        }
        crackInNeighborhood += pow(elem_crack_opening / f->giveNumIP(), 3) * friendsweight [ m ];
        m++;
    }
    return Transp1D :: giveStrain(i, DoFs);
};

//////////////////////////////////////////////////////////
void Transp1DCoupled :: collectInformationsFromNeigborhood() {
    findFriendsFromSimplices();
}

//////////////////////////////////////////////////////////
void Transp1DCoupled :: findFriendsFromSimplices() {
    unordered_set< RigidBodyContact * >allNeighbors;
    vector< RigidBodyContact * >simplexElems;
    Simplex *s;

    //collect all possible candidates
    for ( auto &n: nodes ) {
        s = n->giveSimplex();
        if ( s ) {
            simplexElems = s->giveElements();
            for ( auto &k: simplexElems ) {
                allNeighbors.insert(k);
            }
        }
    }
    //search
    double weight = 0;
    vector< Node * >verts;
    for ( auto &rbc:allNeighbors ) {
        verts = rbc->giveVertices();
        if ( find(verts.begin(), verts.end(), nodes [ 0 ]) != verts.end() && find(verts.begin(), verts.end(), nodes [ 1 ]) != verts.end() ) {
            if ( ndim == 2 ) {
                weight = rbc->giveArea();
            } else if ( ndim == 3 ) {
                weight = ( inttype->giveIPLocation(0) - ( rbc->giveNode(0)->givePoint() + rbc->giveNode(1)->givePoint() ) / 2. ).norm();
            }
            addNewFriend(rbc, weight);
        }
    }
}

//////////////////////////////////////////////////////////
Vector Transp1DCoupled :: giveInternalForces(const Vector &DoFs, bool frozen, double timeStep) {
    Vector intF = Element :: giveInternalForces(DoFs, frozen, timeStep);

    //Biot effect
    double volStrain = 0;
    Simplex *s0 = nodes [ 0 ]->giveSimplex();
    Simplex *s1 = nodes [ 1 ]->giveSimplex();
    if ( !s0 ) {
        if ( s1 ) {
            volStrain = s1->giveVolumetricStrain();
        }
    } else if ( !s1 ) {
        if ( s0 ) {
            volStrain = s0->giveVolumetricStrain();
        }
    } else {
        if ( ( s0->isValid() && s1->isValid() ) || ( !s0->isValid() && !s1->isValid() ) ) {
            volStrain = ( s0->giveVolumetricStrain() + s1->giveVolumetricStrain() ) / 2.;
        } else if ( s0->isValid() ) {
            volStrain = s0->giveVolumetricStrain();
        } else                                                                {
            volStrain = s1->giveVolumetricStrain();
        }
    }

    TrsprtCoupledMaterialStatus *cstat = static_cast< TrsprtCoupledMaterialStatus * >( stats [ 0 ] );
    double volumetricBiotPart = cstat->computeBiotEffect(volStrain, timeStep);
    for ( unsigned i = 0; i < 2; i++ ) {
        intF [ i ] -= volumetricBiotPart * giveVolume(i);
    }

    return intF;
}
