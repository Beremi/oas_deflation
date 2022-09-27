#include "element_discrete.h"
#include "element_container.h"
#include "boundary_condition.h"
#include "material_coulomb_friction.h"
#include "model.h"

using namespace std;

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
void RigidBodyContact :: giveValues(string code, Vector &result) const {
    if ( code.compare("damage") == 0 ) {
        stats [ 0 ]->giveValues(code, result);
    } else if ( code.compare("normal") == 0 ) {
        result.resize(ndim);
        for ( unsigned i = 0; i < ndim; i++ ) {
            result [ i ] = normal(i);
        }
    } else if ( code.compare("t1") == 0 ) {
        result.resize(ndim);
        for ( unsigned i = 0; i < ndim; i++ ) {
            result [ i ] = R(1, i);
        }
    } else if ( code.compare("t2") == 0 ) {
        result.resize(ndim);
        for ( unsigned i = 0; i < ndim; i++ ) {
            result [ i ] = R(1, i);
        }
    } else if ( code.compare("volume") == 0 ) {
        result.resize(1);
        result [ 0 ] = area * length / ndim;
    } else if ( code.compare("area") == 0 ) {
        result.resize(1);
        result [ 0 ] = area;
    } else if ( code.compare("length") == 0 ) {
        result.resize(1);
        result [ 0 ] = length;
    } else {
        MechanicalElement :: giveValues(code, result);
    }
}


/////////////////////////////////////////////////////////
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
    ( void ) x;
    //MyMatrix B
    Matrix B = Matrix :: Zero(ndim, 6 * ( ndim - 1 ) );
    Particle *a = static_cast< Particle * >( nodes [ 0 ] );
    Matrix Aa = a->giveRigidBodyMotionMatrix( inttype->giveIPLocationPointer(0) ) * ( -1. );
    a = static_cast< Particle * >( nodes [ 1 ] );
    Matrix Ab = a->giveRigidBodyMotionMatrix( inttype->giveIPLocationPointer(0) );
    for ( unsigned i = 0; i < ndim; i++ ) {
        for ( unsigned j = 0; j < 3 * ( ndim - 1 ); j++ ) {
            B(i, j) = Aa(i, j);
            B(i, j + 3 * ( ndim - 1 ) ) = Ab(i, j);
        }
    }
    return ( R * B ) / length;
}


//////////////////////////////////////////////////////////
void RigidBodyContact :: setIntegrationPointsAndWeights() {
    stats.resize(1);

    Point t(0, 0, 0);
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
        if ( maxErr > 1e-4 ) {
            cerr << "Vertices are not coplanar!!! Coplanarity error: " << maxErr << endl;
            exit(1);
        }

        //JM: face normal vector made from first 3 vertices
        Point n = ( vert [ 1 ]->givePoint() - vert [ 0 ]->givePoint() ).cross(vert [ 2 ]->givePoint() - vert [ 0 ]->givePoint() );
        n /= n.norm();

        //JM: Perpendicularity check of the beam and face directions
        //JM: normal of the face surface taken from first 3 vertices is (B - A) x (C - A)
        //JM: perpendicularity check: cross (beam, face)=>0
        double prp = ( nodes [ 1 ]->givePoint() - nodes [ 0 ]->givePoint() ).dot(t);
        if ( prp > 1e-8 ) {
            cerr << "Face surface is not perpendicular to beam direction!!! Error: " << prp << endl;
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
        perimeter = 0;
        double ai = 0.0;
        unsigned int j = 0;
        for ( unsigned int i = 0; i < vert.size(); i++ ) {
            j = i + 1;
            if ( i == vert.size() - 1 ) {
                j = 0;
            }
            //triangle area computed as a_i = norm(cross(AB, AC)) / 2
            ai = ( ( vert [ i ]->givePoint() - avgPoint ).cross(vert [ j ]->givePoint() - avgPoint) ).norm() / 2.;
            area += ai;
            perimeter += ( vert [ i ]->givePoint() - vert [ j ]->givePoint() ).norm();
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
    if ( abs(normal.dot(t) ) > 1e-8 ) {
        cout << vert [ 0 ]->givePoint().x() << " " <<  vert [ 0 ]->givePoint().y() <<  " X " << vert [ 1 ]->givePoint().x() << " " <<  vert [ 1 ]->givePoint().y() << endl;
        cout << nodes [ 0 ]->givePoint().x() << " " <<  nodes [ 0 ]->givePoint().y() <<  " X " << nodes [ 1 ]->givePoint().x() << " " <<  nodes [ 1 ]->givePoint().y() << endl;
        cerr << "Error: normal and contact vector are not parallel, error " << normal.dot(t) << " normal v." << normal.x() << " " << normal.y() << " contact v. " << t.x() << " " << t.y() << endl;
        // exit(1);
    }

    // Matrices according to habilitation of Jan Elias (2017, page 42): https://www.vutbr.cz/www_base/vutdisk.php?i=103116a130
    // MyMatrix R;
    if ( ndim == 2 ) {
        t1 = Point(-normal.y(), normal.x(), 0.0);
        R = Matrix :: Zero(2, 2);
        R(0, 0) = normal.x();
        R(0, 1) = normal.y();
        R(1, 0) = t1.x();
        R(1, 1) = t1.y();
    } else if ( ndim == 3 ) {
        // coordinate swap for tangential vector according to https://orbit.dtu.dk/files/126824972/onb_frisvad_jgt2012_v2.pdf
        Point arbit(sqrt(2.), -sqrt(3.), M_PI);
        if ( ( normal - arbit ).norm() < 1e-3 ) {
            t1 = arbit.cross(normal);
        } else {
            // the following results in zeros in stiffness matrix in case of normal in direction of any of global base axes
            if ( abs(normal.x() ) > 1e-3 ) {
                t1 = Point(-normal.y() / normal.x(), 1, 0);
            } else if ( abs(normal.y() ) > 1e-3 ) {
                t1 = Point(0, -normal.z() / normal.y(), 1);
            } else {
                t1 = Point(1, 0, -normal.x() / normal.z() );
            }
        }
        t1.normalize();
        t2 = normal.cross(t1);

        R = Matrix :: Zero(3, 3);
        R(0, 0) = normal.x();
        R(0, 1) = normal.y();
        R(0, 2) = normal.z();
        R(1, 0) = t1.x();
        R(1, 1) = t1.y();
        R(1, 2) = t1.z();
        R(2, 0) = t2.x();
        R(2, 1) = t2.y();
        R(2, 2) = t2.z();
    } else {
        cerr << "Error - RigidBodyContact: dimension " << ndim << "not implemented" << endl;
        exit(EXIT_FAILURE);
    }

    inttype->setIPWeight(0, length * area / ndim);

    stats [ 0 ] = mat->giveNewMaterialStatus(this, 0);
    volume = area * length / ndim;
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
    ( void ) x;
    return Matrix :: Zero(12, 12);  // NOTE JK: this should be based on ndim
}

//////////////////////////////////////////////////////////
Matrix RigidBodyContact :: giveMassMatrix() const {
    Matrix M = Matrix :: Zero( 6 * ( ndim - 1 ), 6 * ( ndim - 1 ) );
    DisMechMaterialStatus *mechstat = static_cast< DisMechMaterialStatus * >( stats [ 0 ] );
    double density = mechstat->giveDensity();
    double m0 = giveVolumeAssociatedWithNode(0) * density; ///mass
    double m1 = giveVolumeAssociatedWithNode(1) * density; ///mass
    if ( ndim == 2 ) {
        // Define points
        Point *A = nodes [ 0 ]->givePointPointer();
        Point *B = nodes [ 1 ]->givePointPointer();
        Point *C = vert [ 0 ]->givePointPointer();
        Point *D = vert [ 1 ]->givePointPointer();
        Point cg0 = ( ( * A ) + ( * C ) + ( * D ) ) / 3.;
        Point cg1 = ( ( * B ) + ( * C ) + ( * D ) ) / 3.;
        Point null(0, 0, 0);
        Point C_ = ( * C ) - ( * A );
        Point D_ = ( * D ) - ( * A );
        Point C__ = ( * C ) - ( * B );
        Point D__ = ( * D ) - ( * B );
        // MassMatrix
        M(0, 0) = M(1, 1) = m0;
        M(3, 3) = M(4, 4) = m1;
        M(2, 2) = density * ( triInertia2D(& null, & C_, & D_) );  // Inertia relative to the node A[0,0]
        M(5, 5) = density * ( triInertia2D(& null, & C__, & D__) );  // Inertia relative to the node B[0,0]
        M(0, 2) = M(2, 0) = -m0 * ( cg0.y() - A->y() );
        M(1, 2) = M(2, 1) = +m0 * ( cg0.x() - A->x() );
        M(3, 5) = M(5, 3) = -m1 * ( cg1.y() - B->y() );
        M(4, 5) = M(5, 4) = +m1 * ( cg1.x() - B->x() );
    } else if ( ndim == 3 )     {
        // Define points
        Point *A = nodes [ 0 ]->givePointPointer();
        Point *B = nodes [ 1 ]->givePointPointer();
        Point *FaceCentroid = inttype->giveIPLocationPointer(0);
        // Mass
        M(0, 0) = M(1, 1) = M(2, 2) = m0;
        M(6, 6) = M(7, 7) = M(8, 8) = m1;
        size_t n = vert.size();
        for ( unsigned i = 0; i < n; i++ ) {
            Point *C = vert [ i ]->givePointPointer();
            Point *D;
            if ( i == 0 ) {
                D = vert [ n - 1 ]->givePointPointer();
            } else {
                D = vert [ i - 1 ]->givePointPointer();
            }
            Point cg0 = ( ( * A ) + ( * C ) + ( * D ) + ( * FaceCentroid ) ) / 4.;
            Point cg1 = ( ( * B ) + ( * C ) + ( * D ) + ( * FaceCentroid ) ) / 4.;
            // Volume of tetrahedron
            double tetraVolume0 = tetraVolumeSigned(A, C, D, FaceCentroid);
            double tetraVolume1 = tetraVolumeSigned(B, C, D, FaceCentroid);
            //cout << tetraVolume0 << endl;
            // Inertia matrix relative to the centroid [0,0,0]
            Point A__ = ( * A ) - cg0;
            Point C__ = ( * C ) - cg0;
            Point D__ = ( * D ) - cg0;
            Point FaceCentroid__ = ( * FaceCentroid ) - cg0;
            Point B___ = ( * B ) - cg1;
            Point C___ = ( * C ) - cg1;
            Point D___ = ( * D ) - cg1;
            Point FaceCentroid___ = ( * FaceCentroid ) - cg1;
            Matrix Icg0 = tetraInertia3D(& A__, & C__, & D__, & FaceCentroid__);
            Matrix Icg1 = tetraInertia3D(& B___, & C___, & D___, & FaceCentroid___);

            // MassMatrix
            M(3, 3) += density * ( Icg0(0, 0) + tetraVolume0 * ( pow( ( cg0.y() - A->y() ), 2 ) + pow( ( cg0.z() - A->z() ), 2 ) ) );
            M(4, 4) += density * ( Icg0(1, 1) + tetraVolume0 * ( pow( ( cg0.x() - A->x() ), 2 ) + pow( ( cg0.z() - A->z() ), 2 ) ) );
            M(5, 5) += density * ( Icg0(2, 2) + tetraVolume0 * ( pow( ( cg0.x() - A->x() ), 2 ) + pow( ( cg0.y() - A->y() ), 2 ) ) );
            M(3, 4) = M(4, 3) += density * ( Icg0(0, 1) - tetraVolume0 * ( ( cg0.x() - A->x() ) * ( cg0.y() - A->y() ) ) );
            M(3, 5) = M(5, 3) += density * ( Icg0(0, 2) - tetraVolume0 * ( ( cg0.x() - A->x() ) * ( cg0.z() - A->z() ) ) );
            M(4, 5) = M(5, 4) += density * ( Icg0(1, 2) - tetraVolume0 * ( ( cg0.y() - A->y() ) * ( cg0.z() - A->z() ) ) );

            M(9, 9) += density * ( Icg1(0, 0) + tetraVolume1 * ( pow( ( cg1.y() - B->y() ), 2 ) + pow( ( cg1.z() - B->z() ), 2 ) ) );
            M(10, 10) += density * ( Icg1(1, 1) + tetraVolume1 * ( pow( ( cg1.x() - B->x() ), 2 ) + pow( ( cg1.z() - B->z() ), 2 ) ) );
            M(11, 11) += density * ( Icg1(2, 2) + tetraVolume1 * ( pow( ( cg1.x() - B->x() ), 2 ) + pow( ( cg1.y() - B->y() ), 2 ) ) );
            M(9, 10) = M(10, 9) += density * ( Icg1(0, 1) - tetraVolume1 * ( ( cg1.x() - B->x() ) * ( cg1.y() - B->y() ) ) );
            M(9, 11) = M(11, 9) += density * ( Icg1(0, 2) - tetraVolume1 * ( ( cg1.x() - B->x() ) * ( cg1.z() - B->z() ) ) );
            M(10, 11) = M(11, 10) += density * ( Icg1(1, 2) - tetraVolume1 * ( ( cg1.y() - B->y() ) * ( cg1.z() - B->z() ) ) );

            M(0, 4) = M(4, 0) += +( tetraVolume0 * density ) * ( cg0.z() - A->z() );
            M(0, 5) = M(5, 0) += -( tetraVolume0 * density ) * ( cg0.y() - A->y() );
            M(1, 3) = M(3, 1) += -( tetraVolume0 * density ) * ( cg0.z() - A->z() );
            M(1, 5) = M(5, 1) += +( tetraVolume0 * density ) * ( cg0.x() - A->x() );
            M(2, 3) = M(3, 2) += +( tetraVolume0 * density ) * ( cg0.y() - A->y() );
            M(2, 4) = M(4, 2) += -( tetraVolume0 * density ) * ( cg0.x() - A->x() );

            M(6, 10) = M(10, 6) += +( tetraVolume1 * density ) * ( cg0.z() - B->z() );
            M(6, 11) = M(11, 6) += -( tetraVolume1 * density ) * ( cg0.y() - B->y() );
            M(7, 9) = M(9, 7) += -( tetraVolume1 * density ) * ( cg0.z() - B->z() );
            M(7, 11) = M(11, 7) += +( tetraVolume1 * density ) * ( cg0.x() - B->x() );
            M(8, 9) = M(9, 8) += +( tetraVolume1 * density ) * ( cg0.y() - B->y() );
            M(8, 10) = M(10, 8) += -( tetraVolume1 * density ) * ( cg0.x() - B->x() );
        }
    }
    //    cout << M << endl;
    //    exit(1);
    return M;
}

//////////////////////////////////////////////////////////
Vector RigidBodyContact :: giveBoundingBox() {
    Vector bbox = giveFacetBoundingBox();

    Point *p;
    for ( auto &v:nodes ) {
        p = v->givePointPointer();
        for ( unsigned i = 0; i < ndim; i++ ) {
            bbox [ 2 * i ]   =  min(bbox [ 2 * i  ], ( * p ) [ i ]);
            bbox [ 2 * i + 1 ] =  max(bbox [ 2 * i + 1 ], ( * p ) [ i ]);
        }
    }
    return bbox;
}

//////////////////////////////////////////////////////////
Vector RigidBodyContact :: giveFacetBoundingBox() {
    Vector bbox = Vector :: Zero(ndim * 2);

    //init
    for ( unsigned i = 0; i < ndim; i++ ) {
        bbox [ 2 * i ]   =  std :: numeric_limits< double > :: infinity();
        bbox [ 2 * i + 1 ] = -std :: numeric_limits< double > :: infinity();
    }

    Point *p;
    for ( auto &v:vert ) {
        p = v->givePointPointer();
        for ( unsigned i = 0; i < ndim; i++ ) {
            bbox [ 2 * i ]   =  min(bbox [ 2 * i  ], ( * p ) [ i ]);
            bbox [ 2 * i + 1 ] =  max(bbox [ 2 * i + 1 ], ( * p ) [ i ]);
        }
    }
    return bbox;
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
Vector RigidBodyContact :: transformVectorToXYZ(Vector &result) const {
    return this->R.transpose() * result;
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
    ( void ) ip_id;
    Point distance = inttype->giveIPLocation(0) - nodes [ node_i ]->givePoint();
    Vector dst = Vector :: Zero(ndim);
    for ( unsigned i = 0; i < ndim; i++ ) {
        if ( i == 0 ) {
            dst [ i ] = distance.x();
        } else if ( i == 1 ) {
            dst [ i ] = distance.y();
        } else if ( i == 2 ) {
            dst [ i ] = distance.z();
        }
    }
    return dst;
}

//////////////////////////////////////////////////////////
double RigidBodyContact :: giveVolumeAssociatedWithNode(unsigned nodenum) const {
    if ( nodenum == 0 ) {
        return ( vert [ 0 ]->givePoint() - nodes [ 0 ]->givePoint() ).dot(normal) * area / ndim;
    } else if ( nodenum == 1 ) {
        return -( vert [ 0 ]->givePoint() - nodes [ 1 ]->givePoint() ).dot(normal) * area / ndim;
    } else {
        cerr << "Error in " << name << ": attempting to reach node number different form 0 or 1." << endl;
        exit(1);
    }
};

//////////////////////////////////////////////////////////
Vector RigidBodyContact :: giveStrain(unsigned i, const Vector &DoFs) {
    //extract volumetric strains from simplices
    double volumetricStrain = 0;
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
    stats [ 0 ]->setParameterValue("volumetric_strain", volumetricStrain);

    return Element :: giveStrain(i, DoFs);
};

//////////////////////////////////////////////////////////
Matrix RigidBodyContact :: giveStiffnessMatrix(string matrixType) const {
    return Element :: giveStiffnessMatrix(matrixType) * ndim; //ndim needs to be included here for discrete elements
}

//////////////////////////////////////////////////////////
Matrix RigidBodyContact :: giveDampingMatrix() const {
    return giveStiffnessMatrix("elastic") * 1e-15;           //rough fix of zeros, here can be anything
}

//////////////////////////////////////////////////////////
Vector RigidBodyContact :: giveInternalForces(const Vector &DoFs, bool frozen, double timeStep) {
    return Element :: giveInternalForces(DoFs, frozen, timeStep) * ndim; //ndim needs to be included here for discrete elements
}

//////////////////////////////////////////////////////////
Vector RigidBodyContact :: integrateLoad(BodyLoad *vl, double time) const {
    return Element :: integrateLoad(vl, time) / ndim;
}

//////////////////////////////////////////////////////////
Vector RigidBodyContact :: integrateInternalSources() {
    return Element :: integrateInternalSources() / ndim;
}






//////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////
void RigidBodyContact :: extrapolateIPValuesToNodes(string code, vector< Vector > &result, Vector &weights) const {
    Vector ipres;
    giveIPValues(code, 0, ipres);
    Vector A = giveVectorToNode(0, 0);
    Vector B = giveVectorToNode(1, 0);
    size_t d;

    weights.resize(2);
    weights [ 0 ] = giveVolumeAssociatedWithNode(0);
    weights [ 1 ] = giveVolumeAssociatedWithNode(1);

    if ( ipres.size() == 0 ) {   //empty answer
        result.resize(0);
    } else if ( ipres.size() == 1 ) {   //scalar times vector //needs to be checked, probably not theoretically correct
        result.resize(ndim);
        for ( d = 0; d < ndim; d++ ) {
            result [ d ].resize(2);
            result [ d ] [ 0 ] =  area * ipres [ 0 ] * abs(A [ d ]);
            result [ d ] [ 1 ] =  area * ipres [ 0 ] * abs(B [ d ]);
        }
    } else if ( ipres.size() == A.size() ) { //vector times vector of same length, symmetrization
        //transform result to xyz
        Vector ipresglobal = transformVectorToXYZ(ipres);

        //dyadic product
        unsigned k = A.size();
        result.resize( ( k * ( k - 1 ) ) / 2 + k);
        for ( d = 0; d < ( k * ( k - 1 ) ) / 2 + k; d++ ) {
            result [ d ].resize(2);
        }
        //diagonal
        for ( d = 0; d < k; d++ ) {
            result [ d ] [ 0 ] =  area * ipresglobal [ d ] * A [ d ];
            result [ d ] [ 1 ] =  -area * ipresglobal [ d ] * B [ d ];
        }
        //off diagonal
        if ( k == 2 ) {
            result [ 2 ] [ 0 ] =  area * ( ipresglobal [ 1 ] * A [ 0 ] + ipresglobal [ 0 ] * A [ 1 ] ) / 2.;
            result [ 2 ] [ 1 ] = -area * ( ipresglobal [ 1 ] * B [ 0 ] + ipresglobal [ 0 ] * B [ 1 ] ) / 2.;
        } else if ( k == 3 ) {
            result [ 3 ] [ 0 ] =  area * ( ipresglobal [ 1 ] * A [ 2 ] + ipresglobal [ 2 ] * A [ 1 ] ) / 2.;
            result [ 3 ] [ 1 ] = -area * ( ipresglobal [ 1 ] * B [ 2 ] + ipresglobal [ 2 ] * B [ 1 ] ) / 2.;
            result [ 4 ] [ 0 ] =  area * ( ipresglobal [ 2 ] * A [ 0 ] + ipresglobal [ 0 ] * A [ 2 ] ) / 2.;
            result [ 4 ] [ 1 ] = -area * ( ipresglobal [ 2 ] * B [ 0 ] + ipresglobal [ 0 ] * B [ 2 ] ) / 2.;
            result [ 5 ] [ 0 ] =  area * ( ipresglobal [ 1 ] * A [ 0 ] + ipresglobal [ 0 ] * A [ 1 ] ) / 2.;
            result [ 5 ] [ 1 ] = -area * ( ipresglobal [ 1 ] * B [ 0 ] + ipresglobal [ 0 ] * B [ 1 ] ) / 2.;
        } else {
            cerr << "Error in " << name << ": transformation of matrix of size " << k << " to vector not implemented" << endl;
            exit(1);
        }
    } else {
        cerr << "Error in " << name << ": dyadic product of vectors of different length in function extrapolateIPValuesToNodes" << endl;
        exit(1);
    }
}

//////////////////////////////////////////////////////////
// COUPLED RBSN ELEMENT
RigidBodyContactCoupled :: RigidBodyContactCoupled(const unsigned dim) : RigidBodyContact(dim) {
    name = "LTCBEAMCoupled";
}

void RigidBodyContactCoupled :: extractPressureFromSimplices() {
    double averagePressure = 0;
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

    stats [ 0 ]->setParameterValue("pressure", averagePressure);
}

//////////////////////////////////////////////////////////
Vector RigidBodyContactCoupled :: giveStrain(unsigned i, const Vector &DoFs) {
    // TODO put this into a separate fn to make it avaliable separately for derived methods
    //extract pressure from simplices
    this->extractPressureFromSimplices();

    return RigidBodyContact :: giveStrain(i, DoFs);
};

//////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////
// RBSN BOUNDARY ELEMENT
RigidBodyBoundary :: RigidBodyBoundary(const unsigned dim) : RigidBodyContact(dim) {
    name = "LTCBoundary";
    active = false;
    numOfNodes = 2; //will be in init to 1
    nodes.resize(2); //will be in init to 1
    name = "LTCBEAM";
    vtk_cell_type = 1; //point
}

//////////////////////////////////////////////////////////
void RigidBodyBoundary :: checkNodeType() const {
    // do nothing here, since the check is already performed during init() together setting correct order of nodes
};

//////////////////////////////////////////////////////////
void RigidBodyBoundary :: init() {
    //check that nodes are one particle and one auxNode
    if ( dynamic_cast< Particle * >( nodes [ 0 ] ) && dynamic_cast< AuxNode * >( nodes [ 1 ] ) ) {
        // this is fine, do nothing, just use it to check if particle and auxnode is there
    } else if ( dynamic_cast< Particle * >( nodes [ 1 ] ) && dynamic_cast< AuxNode * >( nodes [ 0 ] ) ) {
        std :: reverse( this->nodes.begin(), this->nodes.end() );
        std :: reverse( this->vert.begin(), this->vert.end() );
    } else {
        cerr << "Error in " << name << ": nodes must be inherited from Particle and AuxNode, " << nodes [ 0 ]->giveName() << "and " << nodes [ 1 ]->giveName() << " provided" << endl;
    }
    // init of parent class must be done after reverse of node, because the geometrical matrices are calculated in init()
    RigidBodyContact :: init();

    CoulombFrictionMaterial *cfm = dynamic_cast< CoulombFrictionMaterial * >( mat );
    if ( cfm ) {
        active = true;
    }

    numOfNodes = 1;
    nodes.resize(1);
}

//////////////////////////////////////////////////////////
void RigidBodyBoundary :: extrapolateIPValuesToNodes(string code, vector< Vector > &result, Vector &weights) const {
    Vector ipres;
    giveIPValues(code, 0, ipres);
    Vector A = giveVectorToNode(0, 0);
    size_t d;

    weights.resize(2);
    weights [ 0 ] = giveVolumeAssociatedWithNode(0);

    if ( ipres.size() == 0 ) {   //empty answer
        result.resize(0);
    } else if ( ipres.size() == 1 ) {   //scalar times vector //needs to be checked, probably not theoretically correct
        result.resize(ndim);
        for ( d = 0; d < ndim; d++ ) {
            result [ d ].resize(1);
            result [ d ] [ 0 ] =  area * ipres [ 0 ] * abs(A [ d ]);
        }
    } else if ( ipres.size() == A.size() ) { //vector times vector of same length, symmetrization
        //transform result to xyz
        Vector ipresglobal = transformVectorToXYZ(ipres);

        //dyadic product
        unsigned k = A.size();
        result.resize( ( k * ( k - 1 ) ) / 2 + k);
        for ( d = 0; d < ( k * ( k - 1 ) ) / 2 + k; d++ ) {
            result [ d ].resize(1);
        }
        //diagonal
        for ( d = 0; d < k; d++ ) {
            result [ d ] [ 0 ] =  area * ipresglobal [ d ] * A [ d ];
        }
        //off diagonal
        if ( k == 2 ) {
            result [ 2 ] [ 0 ] =  area * ( ipresglobal [ 1 ] * A [ 0 ] + ipresglobal [ 0 ] * A [ 1 ] ) / 2.;
        } else if ( k == 3 ) {
            result [ 3 ] [ 0 ] =  area * ( ipresglobal [ 1 ] * A [ 2 ] + ipresglobal [ 2 ] * A [ 1 ] ) / 2.;
            result [ 4 ] [ 0 ] =  area * ( ipresglobal [ 2 ] * A [ 0 ] + ipresglobal [ 0 ] * A [ 2 ] ) / 2.;
            result [ 5 ] [ 0 ] =  area * ( ipresglobal [ 1 ] * A [ 0 ] + ipresglobal [ 0 ] * A [ 1 ] ) / 2.;
        } else {
            cerr << "Error in " << name << ": transformation of matrix of size " << k << " to vector not implemented" << endl;
            exit(1);
        }
    } else {
        cerr << "Error in " << name << ": dyadic product of vectors of different length in function extrapolateIPValuesToNodes" << endl;
        exit(1);
    }
}

//////////////////////////////////////////////////////////
Vector RigidBodyBoundary :: giveStrain(unsigned i, const Vector &DoFs) {
    ( void ) i;
    ( void ) DoFs;
    Vector f_ext = masterModel->giveSolver()->giveNodalForces();
    double pressure = 0;
    for ( unsigned k = 0; k < ndim; k++ ) {
        pressure += f_ext [ DoFids [ k ] ] * normal [ k ];
    }
    pressure /= area;
    if ( active ) {
        stats [ 0 ]->setParameterValue("normal_stress", pressure);
        return RigidBodyContact :: giveStrain(i, DoFs);
    } else   {
        return Vector :: Zero( ( this->ndim - 1 ) * 3 );
    }
};

//////////////////////////////////////////////////////////
Matrix RigidBodyBoundary :: giveHMatrix(const Point *x) const {
    ( void ) x;
    return Matrix :: Zero( ( this->ndim - 1 ) * 3, ( this->ndim - 1 ) * 3 );
}

//////////////////////////////////////////////////////////
Matrix RigidBodyBoundary :: giveBMatrix(const Point *x) const {
    ( void ) x;
    // MyMatrix B = MyMatrix( ndim, 6 * ( ndim - 1 ) );
    Matrix B = Matrix :: Zero(ndim, 3 * ( ndim - 1 ) );
    Particle *a = static_cast< Particle * >( nodes [ 0 ] );
    Matrix Aa = a->giveRigidBodyMotionMatrix( inttype->giveIPLocationPointer(0) ) * ( -1. );
    for ( unsigned i = 0; i < ndim; i++ ) {
        for ( unsigned j = 0; j < 3 * ( ndim - 1 ); j++ ) {
            B(i, j) = Aa(i, j);
        }
    }
    return ( R * B ) / length;
}

//////////////////////////////////////////////////////////
Matrix RigidBodyBoundary :: giveStiffnessMatrix(std :: string matrixType) const {
    if ( active ) {
        return RigidBodyContact :: giveStiffnessMatrix(matrixType);
    } else {
        return Matrix :: Zero(6, 6);
    }
}

//////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////
// RBSN BOUNDARY ELEMENT COUPLED
RigidBodyBoundaryCoupled :: RigidBodyBoundaryCoupled(const unsigned dim) : RigidBodyContactCoupled(dim) {
    name = "LTCBoundaryCoupled";
}

//////////////////////////////////////////////////////////
void RigidBodyBoundaryCoupled :: checkNodeType() const {
    // do nothing here, since the check is already performed during init() together setting correct order of nodes
};

//////////////////////////////////////////////////////////
void RigidBodyBoundaryCoupled :: init() {
    //check that nodes are one particle and one auxNode
    if ( dynamic_cast< Particle * >( nodes [ 0 ] ) && dynamic_cast< AuxNode * >( nodes [ 1 ] ) ) {
        // this is fine, do nothing, just use it to check if particle and auxnode is there
    } else if ( dynamic_cast< Particle * >( nodes [ 1 ] ) && dynamic_cast< AuxNode * >( nodes [ 0 ] ) ) {
        std :: reverse( this->nodes.begin(), this->nodes.end() );
        std :: reverse( this->vert.begin(), this->vert.end() );
    } else {
        cerr << "Error in " << name << ": nodes must be inherited from Particle and AuxNode, " << nodes [ 0 ]->giveName() << "and " << nodes [ 1 ]->giveName() << " provided" << endl;
    }
    // init of parent class must be done after reverse of node, because the geometrical matrices are calculated in init()
    RigidBodyContact :: init();
}

//////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////
void RigidBodyBoundaryCoupled :: extrapolateIPValuesToNodes(string code, vector< Vector > &result, Vector &weights) const {
    Vector ipres;
    giveIPValues(code, 0, ipres);
    Vector A = giveVectorToNode(0, 0);
    size_t d;

    weights.resize(2);
    weights [ 0 ] = giveVolumeAssociatedWithNode(0);

    if ( ipres.size() == 0 ) {   //empty answer
        result.resize(0);
    } else if ( ipres.size() == 1 ) {   //scalar times vector //needs to be checked, probably not theoretically correct
        result.resize(ndim);
        for ( d = 0; d < ndim; d++ ) {
            result [ d ].resize(1);
            result [ d ] [ 0 ] =  area * ipres [ 0 ] * abs(A [ d ]);
        }
    } else if ( ipres.size() == ndim ) { //vector times vector of same length, symmetrization
        //transform result to xyz
        Vector ipresglobal = transformVectorToXYZ(ipres);

        //dyadic product
        unsigned k = A.size();
        result.resize( ( k * ( k - 1 ) ) / 2 + k);
        for ( d = 0; d < ( k * ( k - 1 ) ) / 2 + k; d++ ) {
            result [ d ].resize(1);
        }
        //diagonal
        for ( d = 0; d < k; d++ ) {
            result [ d ] [ 0 ] =  area * ipresglobal [ d ] * A [ d ];
        }
        //off diagonal
        if ( k == 2 ) {
            result [ 2 ] [ 0 ] =  area * ( ipresglobal [ 1 ] * A [ 0 ] + ipresglobal [ 0 ] * A [ 1 ] ) / 2.;
        } else if ( k == 3 ) {
            result [ 3 ] [ 0 ] =  area * ( ipresglobal [ 1 ] * A [ 2 ] + ipresglobal [ 2 ] * A [ 1 ] ) / 2.;
            result [ 4 ] [ 0 ] =  area * ( ipresglobal [ 2 ] * A [ 0 ] + ipresglobal [ 0 ] * A [ 2 ] ) / 2.;
            result [ 5 ] [ 0 ] =  area * ( ipresglobal [ 1 ] * A [ 0 ] + ipresglobal [ 0 ] * A [ 1 ] ) / 2.;
        } else {
            cerr << "Error in " << name << ": transformation of matrix of size " << k << " to vector not implemented" << endl;
            exit(1);
        }
    } else {
        cerr << "Error in " << name << ": dyadic product of vectors of different length in function extrapolateIPValuesToNodes" << endl;
        exit(1);
    }
}

// //////////////////////////////////////////////////////////
// MyVector RigidBodyBoundaryCoupled :: giveInternalForces(const MyVector &DoFs, bool frozen, double timeStep){
//     MyVector innttff = RigidBodyContactCoupled :: giveInternalForces(DoFs, frozen, timeStep);
//
//     // std::cout << "Boundary int F:" << '\t';
//     // for ( auto const &cc : innttff ) {
//     //     std::cout << cc << '\t';
//     // }
//     // std::cout << '\n';
//
//     // return MyVector((double)0, (this->ndim - 1) * 3);
//     return innttff;
// };

//////////////////////////////////////////////////////////
Vector RigidBodyBoundaryCoupled :: giveStrain(unsigned i, const Vector &DoFs) {
    // DONE  call only separate fn for update of pressure due to transport
    // MyVector dummy = RigidBodyContactCoupled :: giveStrain(i, DoFs);
    ( void ) i;
    ( void ) DoFs;
    this->extractPressureFromSimplices();
    // std::cout << "gstr DoFs size = " << DoFs.size() << '\n';
    return Vector :: Zero( ( this->ndim - 1 ) * 3 );
};

//////////////////////////////////////////////////////////
Matrix RigidBodyBoundaryCoupled :: giveHMatrix(const Point *x) const {
    ( void ) x;
    return Matrix :: Zero( ( this->ndim - 1 ) * 3, ( this->ndim - 1 ) * 3 );
}

//////////////////////////////////////////////////////////
Matrix RigidBodyBoundaryCoupled :: giveBMatrix(const Point *x) const {
    ( void ) x;
    // MyMatrix B = MyMatrix( ndim, 6 * ( ndim - 1 ) );
    Matrix B = Matrix :: Zero(ndim, 3 * ( ndim - 1 ) );
    Particle *a = static_cast< Particle * >( nodes [ 0 ] );
    Matrix Aa = a->giveRigidBodyMotionMatrix( inttype->giveIPLocationPointer(0) ) * ( -1. );
    for ( unsigned i = 0; i < ndim; i++ ) {
        for ( unsigned j = 0; j < 3 * ( ndim - 1 ); j++ ) {
            B(i, j) = Aa(i, j);
        }
    }
    return ( R * B ) / length;
}


//////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////
// TRUSS ELEMENT
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
    ( void ) x;
    //MyMatrix B
    Matrix B = Matrix :: Zero(ndim, 2 * ndim);
    Matrix Aa = Matrix :: Identity(ndim, ndim);
    for ( unsigned i = 0; i < ndim; i++ ) {
        for ( unsigned j = 0; j < ndim; j++ ) {
            B(i, j) = -Aa(i, j);
            B(i, j + ndim) = Aa(i, j);
        }
    }
    return ( R * B ) / length;
}

//////////////////////////////////////////////////////////
Matrix Truss :: giveHMatrix(const Point *x) const {
    ( void ) x;
    return Matrix(0, 0);
}
//////////////////////////////////////////////////////////
Vector Truss :: giveContactStrainNT(const Vector &DoFs) const {
    Vector strain = Bs [ 0 ] * DoFs;
    for ( size_t k = 1; k < ( size_t ) strain.size(); k++ ) {
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
    Point t = Point :: Zero();
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
        Point n = ( vert [ 1 ]->givePoint() - vert [ 0 ]->givePoint() ).cross(vert [ 2 ]->givePoint() - vert [ 0 ]->givePoint() );
        n /= n.norm();
        Point t2;
        if ( fabs(n.x() ) > fabs(n.z() ) ) {
            t2 = Point(-n.y(), n.x(), 0.0f);
        } else {
            t2 = Point(0.0f, -n.z(), n.y() );
        }
        t = t2.cross(n);
        t /= t.norm();

        //JM: Perpendicularity check of the beam and face directions
        //JM: normal of the face surface taken from first 3 vertices is (B - A) x (C - A)
        //JM: perpendicularity check: cross (beam, face)=>0
        double prp = ( nodes [ 1 ]->givePoint() - nodes [ 0 ]->givePoint() ).dot(t);
        if ( prp > 1e-8 ) {
            cerr << "TRSPRT: Face surface is not perpendicular to beam direction!!! Error: " << prp << endl;
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
            ai = ( ( vert [ i ]->givePoint() - avgPoint ).cross(vert [ j ]->givePoint() - avgPoint) ).norm() / 2.;
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

    if ( abs(normal.dot(t) ) > 1e-5 ) {
        cout << vert [ 0 ]->givePoint().x() << " " <<  vert [ 0 ]->givePoint().y() <<  " X " << vert [ 1 ]->givePoint().x() << " " <<  vert [ 1 ]->givePoint().y() << endl;
        cout << nodes [ 0 ]->givePoint().x() << " " <<  nodes [ 0 ]->givePoint().y() <<  " X " << nodes [ 1 ]->givePoint().x() << " " <<  nodes [ 1 ]->givePoint().y() << endl;
        cerr << "TRSPRT: normal and contact vector are not parallel, error " << normal.dot(t) << endl;
        cout << " normal v.:";
        for ( unsigned p = 0; p < ndim; p++ ) {
            cout << "\t" << normal(p);
        }
        cout << endl;
        cout << " contact v.:";
        for ( unsigned p = 0; p < ndim; p++ ) {
            cout << "\t" << t(p);
        }
        cout << endl;
        //exit(1);
    }

    inttype->setIPWeight(0, length * area / ndim);
    stats [ 0 ] = mat->giveNewMaterialStatus(this, 0);
    volume = area * length / ndim;
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
    Matrix B = Matrix :: Zero(1, 2);
    B(0, 0) = -1. / length;
    B(0, 1) = 1. / length;
    return B;
}

//////////////////////////////////////////////////////////
Matrix Transp1D :: giveHMatrix(const Point *x) const {
    ( void ) x;
    Matrix H = Matrix :: Zero(1, 2);
    //double l1 = dot(* x - nodes [ 0 ]->givePoint(), normal);
    //double l2 = dot(nodes [ 1 ]->givePoint() - * x, normal);
    //H(0, 0) = l1 / length;
    //H(0, 1) = l2 / length;
    H(0, 0) = giveVolumeAssociatedWithNode(0) / volume;
    H(0, 1) = giveVolumeAssociatedWithNode(1) / volume;


    return H;
}

//////////////////////////////////////////////////////////
Matrix Transp1D :: giveDampingMatrix() const {
    Matrix S = Matrix :: Zero(2, 2);
    double s = area * stats [ 0 ]->giveDampingTensor()(0, 0) * length /  ( 2. * ndim );

    S(0, 0) = S(1, 1) = s; //finite volume
    if ( BolanderCapacityMatrix ) { //from Bolander's papers
        S(0, 0) = S(1, 1) = 2. / 3. * s;
        S(1, 0) = S(0, 1) = s / 3.;
    }
    return S;
}

//////////////////////////////////////////////////////////
double Transp1D :: giveVolumeAssociatedWithNode(unsigned nodenum) const {
    if ( nodenum == 0 ) {
        return ( vert [ 0 ]->givePoint() - nodes [ 0 ]->givePoint() ).dot(normal) * area / ndim;
    } else if ( nodenum == 1 ) {
        return -( vert [ 0 ]->givePoint() - nodes [ 1 ]->givePoint() ).dot(normal) * area / ndim;
    } else {
        cerr << "Error in " << name << ": attempting to reach node number different form 0 or 1." << endl;
        exit(1);
    }
};

//////////////////////////////////////////////////////////
Vector Transp1D :: giveStrain(unsigned i, const Vector &DoFs) {
    double averagePressure = ( DoFs [ 0 ] * giveVolumeAssociatedWithNode(0) + DoFs [ 1 ] * giveVolumeAssociatedWithNode(1) ) / volume;
    stats [ 0 ]->setParameterValue("pressure", averagePressure);
    return Element :: giveStrain(i, DoFs);
};


//////////////////////////////////////////////////////////
Matrix Transp1D :: giveStiffnessMatrix(string matrixType) const {
    return Element :: giveStiffnessMatrix(matrixType) * ndim; //ndim needs to be included here for discrete elements
}


//////////////////////////////////////////////////////////
Vector Transp1D :: giveInternalForces(const Vector &DoFs, bool frozen, double timeStep) {
    //MyVector Q = Element :: giveInternalForces(DoFs, frozen, timeStep) * ndim; //ndim needs to be included here for discrete elements
    //for (auto p:Q) cout << " " << p;
    return Element :: giveInternalForces(DoFs, frozen, timeStep) * ndim; //ndim needs to be included here for discrete elements
}

//////////////////////////////////////////////////////////
Vector Transp1D :: integrateLoad(BodyLoad *vl, double time) const {
    return Element :: integrateLoad(vl, time) / ndim;
}

//////////////////////////////////////////////////////////
Vector Transp1D :: integrateInternalSources() {
    return Element :: integrateInternalSources() / ndim;
}

//////////////////////////////////////////////////////////
Vector Transp1D :: giveVectorToNode(const unsigned &node_i, const unsigned &ip_id) const {
    ( void ) ip_id;
    Point distance = inttype->giveIPLocation(0) - nodes [ node_i ]->givePoint();
    Vector dst = Vector :: Zero(ndim);
    for ( unsigned i = 0; i < ndim; i++ ) {
        if ( i == 0 ) {
            dst [ i ] = distance.x();
        } else if ( i == 1 ) {
            dst [ i ] = distance.y();
        } else if ( i == 2 ) {
            dst [ i ] = distance.z();
        }
    }
    return dst;
}

//////////////////////////////////////////////////////////
void Transp1D :: extrapolateIPValuesToNodes(string code, vector< Vector > &result, Vector &weights) const {
    Vector ipres;
    giveIPValues(code, 0, ipres);
    Vector A = giveVectorToNode(0, 0);
    Vector B = giveVectorToNode(1, 0);
    size_t d;

    weights.resize(2);
    weights [ 0 ] = giveVolumeAssociatedWithNode(0);
    weights [ 1 ] = giveVolumeAssociatedWithNode(1);

    if ( ipres.size() == 0 ) {   //empty answer
        result.resize(0);
    } else if ( ipres.size() == 1 ) {   //scalar times vector //needs to be checked, probably not theoretically correct
        result.resize(ndim);
        for ( d = 0; d < ndim; d++ ) {
            result [ d ].resize(2);
            result [ d ] [ 0 ] =  area * ipres [ 0 ] * abs(A [ d ]);
            result [ d ] [ 1 ] =  area * ipres [ 0 ] * abs(B [ d ]);
        }
    } else {
        cerr << "Error in " << name << ": dyadic product of vectors in function extrapolateIPValuesToNodes from transport 1D element" << endl;
        exit(1);
    }
}

//////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////
// 1D TRANSPORT ELEMENT COUPLED WITH MECHANICS
//////////////////////////////////////////////////////////
void Transp1DCoupled :: init() {
    Transp1D :: init(); //calling base class method;
}

//////////////////////////////////////////////////////////
void Transp1DCoupled :: giveValues(string code, Vector &result) const {
    if ( code.compare("numOfFriends") == 0 ) {
        result.resize(0);
        result [ 0 ] = friends.size();
    } else {
        Transp1D :: giveValues(code, result);
    }
};

//////////////////////////////////////////////////////////
void Transp1DCoupled :: addNewFriend(RigidBodyContact *f, double weight) {
    friends.push_back(f);
    friendsweight.push_back(weight);
}

//////////////////////////////////////////////////////////
Vector Transp1DCoupled :: giveStrain(unsigned i, const Vector &DoFs) {
    //crack opening
    double crackInNeighborhood = 0;
    double crackVolume = 0.;
    double elem_crack_opening;
    size_t m = 0;
    Vector res;
    for ( auto &f: friends ) {
        elem_crack_opening = 0.;
        for ( unsigned k = 0; k < f->giveNumIP(); k++ ) {
            f->giveIPValues("tempCrackOpening", k, res);
            elem_crack_opening += abs(res [ 0 ]);
        }
        crackInNeighborhood += pow(elem_crack_opening / f->giveNumIP(), 3) * friendsweight [ m ];
        if ( ndim == 3 ) {
            crackVolume += ( elem_crack_opening / f->giveNumIP() ) * f->giveArea() * length / f->givePerimeter();
        } else if ( ndim == 2 ) {
            crackVolume += ( elem_crack_opening / f->giveNumIP() ) * f->giveArea();
        }
        m++;
    }
    stats [ 0 ]->setParameterValue("crack_opening", crackInNeighborhood);
    stats [ 0 ]->setParameterValue("crack_volume", crackVolume);

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
    } else if ( s0 && s1 ) {
        if ( ( s0->isValid() && s1->isValid() ) || ( !s0->isValid() && !s1->isValid() ) ) {
            volStrain = ( s0->giveVolumetricStrain() + s1->giveVolumetricStrain() ) / 2.;
        } else if ( s0->isValid() ) {
            volStrain = s0->giveVolumetricStrain();
        } else {
            volStrain = s1->giveVolumetricStrain();
        }
    } else {
        cerr << "no simplex found" << endl;
        exit(1);
    }


    stats [ 0 ]->setParameterValue("volumetric_strain", volStrain); //3 is there to obtain mechanical volumetric strain

    return Transp1D :: giveStrain(i, DoFs);
};

//////////////////////////////////////////////////////////
void Transp1DCoupled :: collectInformationsFromNeigborhood() {
    findFriendsFromSimplices();
}

//////////////////////////////////////////////////////////
void Transp1DCoupled :: findFriendsFromSimplices() {
    friends.resize(0);
    friendsweight.resize(0);

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
                weight = 1.;//rbc->giveArea();
            } else if ( ndim == 3 ) {
                weight = ( inttype->giveIPLocation(0) - ( rbc->giveNode(0)->givePoint() + rbc->giveNode(1)->givePoint() ) / 2. ).norm();
            }
            addNewFriend(rbc, weight);
        }
    }
}
