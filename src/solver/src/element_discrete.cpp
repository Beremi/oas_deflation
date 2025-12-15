#include "element_discrete.h"
#include "element_container.h"
#include "boundary_condition.h"
#include "material_coulomb_friction.h"
#include "model.h"
#include "material_vectorial.h"

using namespace std;

//////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////
// RBSN ELEMENT
RigidBodyContact :: RigidBodyContact(const unsigned dim) : Element(dim) {
    physicalFields [ 0 ] = true; //does mechanics
    numOfNodes = 2;
    nodes.resize(2);
    name = "LTCBEAM";
    vtk_cell_type = 3;
    shafunc = new Linear1DLineShapeF();
    inttype = new IntegrDiscrete1();
    areIPLocsInNaturalCoords = false;    
    projectArea = true;
    intPoints = "centroid";
    userDefinedCentroid = false;
    ignoreNegativeAreas = false;
}


//////////////////////////////////////////////////////////
void RigidBodyContact :: giveIPValues(std :: string code, unsigned ipnum, Vector &result) const {
    if ( ipnum >= inttype->giveNumIP() ) {
        std :: cerr << name <<  " Error: intergration point number " << ipnum << " exceeds number of integration points" << std :: endl;
        exit(1);
    }
    if ( code.compare("location") == 0 ) {
        result.resize(ndim);
        Point *h = inttype->giveIPLocationPointer(ipnum);
        for ( unsigned k = 0; k < ndim; k++ ) {
            result [ k ] = ( * h ) [ k ];
        }
    } else {
        Element :: giveIPValues(code, ipnum, result);
    }
};


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
        if ( ndim < 3 ) {
            result.resize(0);
        } else {
            result.resize(ndim);
            for ( unsigned i = 0; i < ndim; i++ ) {
                result [ i ] = R(2, i);
            }
        }
    } else if ( code.compare("volume") == 0 ) {
        result.resize(1);
        result [ 0 ] = area * length / ndim;
    } else if ( code.compare("volumeA") == 0 ) {
        result.resize(1);
        result [ 0 ] = giveVolumeAssociatedWithNode(0);
    } else if ( code.compare("volumeB") == 0 ) {
        result.resize(1);
        result [ 0 ] = giveVolumeAssociatedWithNode(1);
    } else if ( code.compare("area") == 0 ) {
        result.resize(1);
        result [ 0 ] = area;
    } else if ( code.compare("length") == 0 ) {
        result.resize(1);
        result [ 0 ] = length;
    } else {
        Element :: giveValues(code, result);
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
    string param;
    while ( iss >> param ) {
        if ( param.compare("integration") == 0 ) {
            iss >> intPoints;
        } else if ( param.compare("centroid") == 0 && ndim == 3 ) {
            userDefinedCentroid = true;
            double x, y, z;
            iss >> x >> y >> z;
            userCentroid = Point(x, y, z);
        } else if ( param.compare("area_projection") == 0 ) {
            iss >> projectArea;
        } else if ( param.compare("ignore_negative_area") == 0 ) {
            ignoreNegativeAreas = true;
        }
    }
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

    //check that material is VectMechMat
    VectMechMaterial *p = dynamic_cast< VectMechMaterial * >( mat );
    if ( !p ) {
        cerr << "Error in " << name << ": material must be inherited from VectMechMaterial, " << mat->giveName() << " provided" << endl;
        exit(1);
    }
}

//////////////////////////////////////////////////////////

Matrix RigidBodyContact :: giveBMatrix(const Point *x) const {
    ( void ) x;
    //MyMatrix B
    Matrix B = Matrix :: Zero(ndim, 6 * ( ndim - 1 ) );
    Particle *a = static_cast< Particle * >( nodes [ 0 ] );
    Matrix Aa = a->giveRigidBodyMotionMatrix(x) * ( -1. );
    a = static_cast< Particle * >( nodes [ 1 ] );
    Matrix Ab = a->giveRigidBodyMotionMatrix(x);
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
    normal = nodes [ 1 ]->givePoint() - nodes [ 0 ]->givePoint();
    length = normal.norm();
    normal = normal / length;

    if ( ndim == 2 ) {
        if ( !( vert.size() == 2 ) ) {
            cerr << "Error: exactly 2 vertices must be involved, " << vert.size() << " provided" << endl;
            exit(1);
        }

        t1 = vert [ 1 ]->givePoint() - vert [ 0 ]->givePoint();
        area = t1.norm();
        Point faceNormal;
        if ( area < 1e-14 ) {
            faceNormal = Point(-t1 [ 1 ] / 1e-14, t1 [ 0 ] / 1e-14, 0);
        } else {
            faceNormal = Point(-t1 [ 1 ] / area, t1 [ 0 ] / area, 0);
        }
        if ( projectArea ) {
            area = abs( faceNormal.dot(normal) ) * area;
            t1 += normal * t1.dot(normal);
            t1 /= max(t1.norm(), 1e-14);
        }

        unsigned n;
        if ( intPoints.compare("centroid") == 0 ) {
            n = 1;
        } else if ( is_positive_integer(intPoints) ) {
            n = stoi(intPoints);
        } else {
            cerr << "RigidBodyContact Error: unknown value of intPoints parameter: " << intPoints << endl;
            exit(1);
        }

        centroid = ( vert [ 0 ]->givePoint() + vert [ 1 ]->givePoint() ) / 2.;
        IntegrDiscrete1 *it = static_cast< IntegrDiscrete1 * >( inttype );
        it->setNumIP(n);
        //Gauss integration, not used
        if ( false && n < 5 ) {
            Vector locs, weis;
            giveGaussIntegrationPointAndWeights(n, locs, weis);
            for ( unsigned i = 0; i < inttype->giveNumIP(); i++ ) {
                inttype->setIPLocation(i, centroid + t1 * area * locs [ i ] / 2.);
                inttype->setIPWeight(i, length * area * weis [ i ] / 2. / ndim);
            }
        } else {  //equidistant
            for ( unsigned i = 0; i < inttype->giveNumIP(); i++ ) {
                inttype->setIPLocation( i, centroid + t1 * area * ( ( i + 0.5 ) / n - 0.5 ) );
                inttype->setIPWeight( i, length * area / ( ndim * n ) );
            }
        }
    } else {
        Point avgPoint = Point(0.0, 0.0, 0.0);
        for ( unsigned int i = 0; i < vert.size(); i++ ) {
            avgPoint += vert [ i ]->givePoint();
        }
        avgPoint /= vert.size();

        if ( userDefinedCentroid ) {
            avgPoint = userCentroid;
        }

        double diff = 1e6;

        unsigned iter = 0;
        unsigned int j = 0;
        double ai = 0.0;
        while ( diff > 1e-10 ) {
            centroid = Point(0.0, 0.0, 0.0);
            area = 0.0;
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
            if ( area > 1e-20 ) {
                centroid /= area;
            } else {
                centroid = avgPoint;
            }
            iter++;
            diff = ( centroid - avgPoint ).norm();
            avgPoint = centroid;
        }

        if ( userDefinedCentroid ) {
            avgPoint = userCentroid;
            centroid = userCentroid;
        }

        Vector ais( vert.size() );     //projected areas of individual triangles
        Point ni;
        area = 0.0;
        perimeter = 0;
        j = 0;
        for ( unsigned int i = 0; i < vert.size(); i++ ) {
            j = i + 1;
            if ( i == vert.size() - 1 ) {
                j = 0;
            }
            ni = ( vert [ i ]->givePoint() - centroid ).cross(vert [ j ]->givePoint() - centroid);
            ais [ i ] = ni.norm() / 2;
            ni.normalize();
            ais [ i ] *= ni.dot(normal);
            if ( ais [ i ] < 1e-15 ) {
                if ( !ignoreNegativeAreas ) {
                    //cout << "RigidBodyContact Warning: negative area " << ais [ i ] << ", incorrect orientation of verices, corrected automatically" << endl;
                    ais [ i ] = std :: abs(ais [ i ]);
                } else {
                    ais [ i ] = 0.;
                }
            }
            area += ais [ i ];
            ni = ( vert [ i ]->givePoint() - vert [ j ]->givePoint() );   //just for periemeter
            ni = ni - normal * ( normal.dot(ni) );
            perimeter += ni.norm();
        }

        IntegrDiscrete1 *it = static_cast< IntegrDiscrete1 * >( inttype );
        if ( intPoints.compare("centroid") == 0 ) {
            it->setNumIP(1);
            inttype->setIPLocation(0, centroid);
            inttype->setIPWeight(0, length * area / ndim);
        } else if ( intPoints.compare("triangles") == 0 ) {
            it->setNumIP( vert.size() );
            j = 0;
            for ( unsigned int i = 0; i < vert.size(); i++ ) {
                j = i + 1;
                if ( i == vert.size() - 1 ) {
                    j = 0;
                }
                inttype->setIPWeight(i, length * ais [ i ] / ndim);
                inttype->setIPLocation(i, ( centroid + vert [ i ]->givePoint() + vert [ j ]->givePoint() ) / 3.0);
            }
        } else if ( intPoints.compare("triangles3") == 0 ) { //at medians of triangle
            it->setNumIP(vert.size() * 3);
            vector< Point >natcoords(3);
            double s = 1. / 6.;
            //medians of triangle in natural coordinate system
            natcoords [ 0 ] = Point(1. - 2 * s, s, s);
            natcoords [ 1 ] = Point(1. - 5 * s, s, 4 * s);
            natcoords [ 2 ] = Point(1. - 5 * s, 4 * s, s);

            j = 0;
            for ( unsigned int i = 0; i < vert.size(); i++ ) {
                j = i + 1;
                if ( i == vert.size() - 1 ) {
                    j = 0;
                }
                inttype->setIPWeight( 3 * i, length * ais [ i ] / ( 3 * ndim ) );
                inttype->setIPWeight( 3 * i + 1, length * ais [ i ] / ( 3 * ndim ) );
                inttype->setIPWeight( 3 * i + 2, length * ais [ i ] / ( 3 * ndim ) );
                inttype->setIPLocation( 3 * i, centroid * natcoords [ 0 ].x() + vert [ i ]->givePoint() * natcoords [ 0 ].y() + vert [ j ]->givePoint() * natcoords [ 0 ].z() );
                inttype->setIPLocation( 3 * i + 1, centroid * natcoords [ 1 ].x() + vert [ i ]->givePoint() * natcoords [ 1 ].y() + vert [ j ]->givePoint() * natcoords [ 1 ].z() );
                inttype->setIPLocation( 3 * i + 2, centroid * natcoords [ 2 ].x() + vert [ i ]->givePoint() * natcoords [ 2 ].y() + vert [ j ]->givePoint() * natcoords [ 2 ].z() );
            }
        } else {
            cerr << "RigidBodyContact Error: unknown value of intPoints parameter: " << intPoints << endl;
            exit(1);
        }
    }

    // Matrices according to habilitation of Jan Elias (2017, page 42): https://www.vutbr.cz/www_base/vutdisk.php?i=103116a130 integration triangles
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

    stats.resize( inttype->giveNumIP() );

    for ( unsigned i = 0; i < inttype->giveNumIP(); i++ ) {
        stats [ i ] = mat->giveNewMaterialStatus(this, i);
    }
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
}

//////////////////////////////////////////////////////////
Matrix RigidBodyContact :: giveHMatrix(const Point *x) const {
    ( void ) x;
    return Matrix :: Zero(12, 12);  // NOTE JK: this should be based on ndim
}

//////////////////////////////////////////////////////////
void RigidBodyContact :: computeMassMatrix() {
    massM = Matrix :: Zero( 6 * ( ndim - 1 ), 6 * ( ndim - 1 ) );
    VectMechMaterialStatus *mechstat = static_cast< VectMechMaterialStatus * >( stats [ 0 ] );
    double density = mechstat->giveDensity();
    double m0 = giveVolumeAssociatedWithNode(0) * density; ///mass
    double m1 = giveVolumeAssociatedWithNode(1) * density; ///mass
    if ( ndim == 2 ) {
        // Define points
        Point *A = nodes [ 0 ]->givePointPointer();
        Point *B = nodes [ 1 ]->givePointPointer();
        Point *C = vert [ 0 ]->givePointPointer();
        Point *D = vert [ 1 ]->givePointPointer();
        Point cg0 = ( nodes [ 0 ]->givePoint() + vert [ 0 ]->givePoint() + vert [ 1 ]->givePoint() ) / 3.;
        Point cg1 = ( nodes [ 1 ]->givePoint() + vert [ 0 ]->givePoint() + vert [ 1 ]->givePoint() ) / 3.;
        Point null(0, 0, 0);
        // MassMatrix
        massM(0, 0) = massM(1, 1) = m0;
        massM(3, 3) = massM(4, 4) = m1;
        Point C_ = ( * C ) - ( * A );
        Point D_ = ( * D ) - ( * A );
        massM(2, 2) = density * ( triInertia2D(& null, & C_, & D_) );  // Inertia relative to the node A[0,0]
        C_ = ( * C ) - ( * B );
        D_ = ( * D ) - ( * B );
        massM(5, 5) = density * ( triInertia2D(& null, & C_, & D_) );  // Inertia relative to the node B[0,0]
        massM(0, 2) = massM(2, 0) = -m0 * ( cg0.y() - A->y() );
        massM(1, 2) = massM(2, 1) = +m0 * ( cg0.x() - A->x() );
        massM(3, 5) = massM(5, 3) = -m1 * ( cg1.y() - B->y() );
        massM(4, 5) = massM(5, 4) = +m1 * ( cg1.x() - B->x() );
    } else if ( ndim == 3 ) {
        // Define points
        Point *A = nodes [ 0 ]->givePointPointer();
        // Mass
        massM(0, 0) = massM(1, 1) = massM(2, 2) = m0;
        massM(6, 6) = massM(7, 7) = massM(8, 8) = m1;
        Point *D, cg, A_, C_, D_, centroid_;
        Point null(0, 0, 0);
        Point *C = vert [ vert.size() - 1 ]->givePointPointer();
        double tetraVolume;
        for ( unsigned i = 0; i < vert.size(); i++ ) {
            D = C;
            C = vert [ i ]->givePointPointer();
            for ( unsigned k = 0; k < 2; k++ ) {
                A = nodes [ k ]->givePointPointer();
                cg = ( ( * A ) + ( * C ) + ( * D ) + centroid ) / 4.;
                tetraVolume = abs( tetraVolumeSigned(A, C, D, & centroid) );

                // Inertia matrix relative to the centroid [0,0,0]
                A_ = ( * A ) - cg;
                C_ = ( * C ) - cg;
                D_ = ( * D ) - cg;
                centroid_ = centroid - cg;
                Matrix I = tetraInertia3D(& A_, & C_, & D_, & centroid_);

                // MassMatrix
                massM(6 * k + 3, 6 * k + 3) += density * ( I(0, 0) + tetraVolume * ( pow( ( cg.y() - A->y() ), 2 ) + pow( ( cg.z() - A->z() ), 2 ) ) );
                massM(6 * k + 4, 6 * k + 4) += density * ( I(1, 1) + tetraVolume * ( pow( ( cg.x() - A->x() ), 2 ) + pow( ( cg.z() - A->z() ), 2 ) ) );
                massM(6 * k + 5, 6 * k + 5) += density * ( I(2, 2) + tetraVolume * ( pow( ( cg.x() - A->x() ), 2 ) + pow( ( cg.y() - A->y() ), 2 ) ) );
                massM(6 * k + 3, 6 * k + 4) += density * ( I(0, 1) - tetraVolume * ( ( cg.x() - A->x() ) * ( cg.y() - A->y() ) ) );
                massM(6 * k + 3, 6 * k + 5) += density * ( I(0, 2) - tetraVolume * ( ( cg.x() - A->x() ) * ( cg.z() - A->z() ) ) );
                massM(6 * k + 4, 6 * k + 5) += density * ( I(1, 2) - tetraVolume * ( ( cg.y() - A->y() ) * ( cg.z() - A->z() ) ) );

                massM(6 * k, 6 * k + 4)   += tetraVolume * density * ( cg.z() - A->z() );
                massM(6 * k, 6 * k + 5)   -= tetraVolume * density * ( cg.y() - A->y() );
                massM(6 * k + 1, 6 * k + 3) -= tetraVolume * density * ( cg.z() - A->z() );
                massM(6 * k + 1, 6 * k + 5) += tetraVolume * density * ( cg.x() - A->x() );
                massM(6 * k + 2, 6 * k + 3) += tetraVolume * density * ( cg.y() - A->y() );
                massM(6 * k + 2, 6 * k + 4) -= tetraVolume * density * ( cg.x() - A->x() );
            }
        }
        //symmetric
        for ( unsigned k = 0; k < 6; k++ ) {
            for ( unsigned l = max( k + 1, unsigned( 3 ) ); l < 6; l++ ) {
                massM(l, k) = massM(k, l);
                massM(l + 6, k + 6) = massM(k + 6, l + 6);
            }
        }
    }
}

//////////////////////////////////////////////////////////
Vector RigidBodyContact :: giveBoundingBox() const {
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
Vector RigidBodyContact :: giveFacetBoundingBox() const {
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
Vector RigidBodyContact :: transformToGlobal(const Vector &DoFs) const {
    return this->R.transpose() * DoFs;
}

//////////////////////////////////////////////////////////
Vector RigidBodyContact :: transformToLocal(const Vector &DoFs) const {
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
        if ( s->isUpdated() ) {
            validSnum++;
            volumetricStrain += s->giveVolumetricStrain();
        }
    }

    if ( validSnum > 0 ) {
        volumetricStrain /= validSnum;
    }
    if ( volumetricStrain == 0 ) {
        for ( auto &s: simplices ) {
            if ( s->isUpdated() ) {
                validSnum++;
            }
        }
    }
    for ( auto &s:stats ) {
        s->setParameterValue("volumetric_strain", volumetricStrain);
    }

    return Element :: giveStrain(i, DoFs);
};

//////////////////////////////////////////////////////////
Matrix RigidBodyContact :: giveStiffnessMatrix(string matrixType) const {
    return Element :: giveStiffnessMatrix(matrixType) * ndim; //ndim needs to be included here for discrete elements
}

//////////////////////////////////////////////////////////
void RigidBodyContact :: computeDampingMatrix() {
    dampC = giveStiffnessMatrix("elastic") * 1e-15;           //rough fix of zeros, here can be anything
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
bool RigidBodyContact :: isPointInside(Point *xn, const Point *x) const {
    if ( ndim == 3 ) {
        //https://stackoverflow.com/questions/25179693/how-to-check-whether-the-point-is-in-the-tetrahedron-or-not
        Matrix M = Matrix :: Zero(3, 3);
        Point X, Y, A;
        for (unsigned j = 0; j < 2; j++) {
            Point C = vert [ vert.size() - 1 ]->givePoint();
            A = nodes [ j ]->givePoint();
            X = centroid - A;
            M(0, 0) = X [ 0 ];
            M(1, 0) = X [ 1 ];
            M(2, 0) = X [ 2 ];
            for (unsigned i = 0; i < vert.size(); i++) {
                X = C - A;
                M(0, 1) = X [ 0 ];
                M(1, 1) = X [ 1 ];
                M(2, 1) = X [ 2 ];
                X = vert [ i ]->givePoint() - A;
                M(0, 2) = X [ 0 ];
                M(1, 2) = X [ 1 ];
                M(2, 2) = X [ 2 ];
                Y = M.colPivHouseholderQr().solve( ( * x ) - A);
                if ( Y [ 0 ] >= 0. && Y [ 1 ] >= 0. && Y [ 2 ] >= 0. && Y [ 0 ] + Y [ 1 ] + Y [ 2 ] <= 1. ) {
                    * xn = Point(j, 0, 0);
                    return true;
                }
                C = vert [ i ]->givePoint();
            }
        }
    }
    return false;
}

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
    } else if ( ipres.size() >= A.size() ) { //vector times vector of same length, symmetrization
        ipres.resize( A.size() );
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
double RigidBodyContact :: giveDissipatedEnergy() const {
    return Element::giveDissipatedEnergy()*ndim;
}

//////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////
// RBSN ELEMENT with rotational stiffness
RigidBodyContactWithRotationalStiffness :: RigidBodyContactWithRotationalStiffness(const unsigned dim) : RigidBodyContact(dim) {
    name = "RigidBodyContactWithRotationalStiffness";
}

//////////////////////////////////////////////////////////

Matrix RigidBodyContactWithRotationalStiffness :: giveBMatrix(const Point *x) const {
    Matrix B0 = RigidBodyContact :: giveBMatrix(x);
    Matrix C0;
    if ( ndim == 2 ) {
        C0 = Matrix :: Zero(1, 6 * ( ndim - 1 ) );
        C0(0, 2) = -1 / length;
        C0(0, 5) =  1. / length;
    } else if ( ndim == 3 ) {
        C0 = Matrix :: Zero(ndim, 6 * ( ndim - 1 ) );
        for ( unsigned i = 0; i < 3; i++ ) {
            C0(i, 3 + i) = -1 / length;
            C0(i, 9 + i) =  1 / length;
        }
        C0 = R * C0;
    }
    Matrix B = Matrix :: Zero( B0.rows() + C0.rows(), B0.cols() );
    B << B0, C0;
    return B;
}

//////////////////////////////////////////////////////////
void RigidBodyContactWithRotationalStiffness :: setIntegrationPointsAndWeights() {
    RigidBodyContact :: setIntegrationPointsAndWeights();
    if ( ndim == 3 ) {
        Inertia = 0;
    } else {
        Inertia = pow(area / inttype->giveNumIP(), 3) / 12.;
    }
}


//////////////////////////////////////////////////////////
Matrix RigidBodyContactWithRotationalStiffness :: giveStiffnessMatrix(string matrixType) const {
    return RigidBodyContact :: giveStiffnessMatrix(matrixType);
}

//////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////
// RBSN ELEMENT with rotational stiffness
RigidBodyContactWithDecoupledRotationsAndTranslations :: RigidBodyContactWithDecoupledRotationsAndTranslations(const unsigned dim) : RigidBodyContactWithRotationalStiffness(dim) {
    name = "RigidBodyContactWithDecoupledRotationsAndTranslations";
}

//////////////////////////////////////////////////////////
Matrix RigidBodyContactWithDecoupledRotationsAndTranslations :: giveBMatrix(const Point *x) const {
    Matrix B = RigidBodyContactWithRotationalStiffness :: giveBMatrix(x);
    //delete coupling terms
    if ( ndim == 2 ) {
        B(0, 2) = 0;
        B(1, 2) = 0;
        B(0, 5) = 0;
        B(1, 5) = 0;
    } else if ( ndim == 3 ) {
        for ( unsigned p = 0; p < 3; p++ ) {
            for ( unsigned q = 0; q < 3; q++ ) {
                B(p, 3 + q) = 0;
                B(p, 9 + q) = 0;
            }
        }
    }
    return B;
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
        for ( auto &s:stats ) {
            s->setParameterValue("normal_stress", pressure);
        }
        return RigidBodyContact :: giveStrain(i, DoFs);
    } else {
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
Vector Truss :: giveContactStrainNT() const {
    Vector strain = RigidBodyContact :: giveContactStrainNT();
    for ( size_t k = 1; k < ( size_t ) strain.size(); k++ ) {
        strain [ k ] = 0;                                  //only normal strain active in truss
    }
    return strain;
};

//////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////
// 1D TRANSPORT ELEMENT
DiscreteTrsprtElem :: DiscreteTrsprtElem(const unsigned dim) : Element(dim) {
    nodes.resize(2);
    bound = false;
    name = "LTCTRSP";
    BolanderCapacityMatrix = false;
    shafunc = new Linear1DLineShapeF();
    inttype = new IntegrDiscrete1();
    vtk_cell_type = 3;
    physicalFields [ 1 ] = true; //transport
    ignoreNegativeAreas = false;
    projectArea = true;
}

//////////////////////////////////////////////////////////
void DiscreteTrsprtElem :: readFromLine(istringstream &iss, NodeContainer *fullnodes, MaterialContainer *fullmatrs) {
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
void DiscreteTrsprtElem :: setIntegrationPointsAndWeights() {
    IntegrDiscrete1 *it = static_cast< IntegrDiscrete1 * >( inttype );
    it->setNumIP(1);
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

        //project area
        Point faceNormal = Point(-t.y(), t.x(), 0);
        area = abs( faceNormal.dot(normal) ) * area;
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
        //JM: perpendicularity check: cross (beam, face)=>0 integration triangles
        //double prp = ( nodes [ 1 ]->givePoint() - nodes [ 0 ]->givePoint() ).dot(t);
        //if ( prp > 1e-8 ) {
        //    cerr << "TRSPRT: Face surface is not perpendicular to beam direction!!! Error: " << prp << endl;
        //    //  exit(1);
        //}



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
        Point ni;
        for ( unsigned int i = 0; i < vert.size(); i++ ) {
            j = i + 1;
            if ( i == vert.size() - 1 ) {
                j = 0;
            }
            //triangle area computed as a_i = norm(cross(AB, AC)) / 2
            ni = ( vert [ i ]->givePoint() - avgPoint ).cross(vert [ j ]->givePoint() - avgPoint);
            ai = ni.norm() / 2.;
            ni.normalize();
            ai *= ni.dot(normal);
            if ( ai < 1e-15 ) {
                if ( !ignoreNegativeAreas ) {
                    //cout << "DiscreteTrsprtElem Warning: negative area " << ai << ", incorrect orientation of verices, corrected automatically" << endl;
                    ai = std :: abs(ai);
                } else {
                    ai = 0.;
                }
            }
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

    /*
     * if ( abs( normal.dot(t) ) > 1e-5 ) {
     *  cout << vert [ 0 ]->givePoint().x() << " " <<  vert [ 0 ]->givePoint().y() <<  " X " << vert [ 1 ]->givePoint().x() << " " <<  vert [ 1 ]->givePoint().y() << endl;
     *  cout << nodes [ 0 ]->givePoint().x() << " " <<  nodes [ 0 ]->givePoint().y() <<  " X " << nodes [ 1 ]->givePoint().x() << " " <<  nodes [ 1 ]->givePoint().y() << endl;
     *  cerr << "TRSPRT: normal and contact vector are not parallel, error " << normal.dot(t) << endl;
     *  cout << " normal v.:";
     *  for ( unsigned p = 0; p < ndim; p++ ) {
     *      cout << "\t" << normal(p);
     *  }
     *  cout << endl;
     *  cout << " contact v.:";
     *  for ( unsigned p = 0; p < ndim; p++ ) {
     *      cout << "\t" << t(p);
     *  }
     *  cout << endl;
     *  //exit(1);
     * }
     */

    inttype->setIPWeight(0, length * area / ndim);
    stats [ 0 ] = mat->giveNewMaterialStatus(this, 0);
    volume = area * length / ndim;
}

//////////////////////////////////////////////////////////
void DiscreteTrsprtElem :: init() {
    Element :: init(); //calling base class method;

    checkNodeAndMaterialType();
}

//////////////////////////////////////////////////////////
void DiscreteTrsprtElem :: checkNodeAndMaterialType() const {
    for ( unsigned i = 0; i < 2; i++ ) {
        TrsNode *p = dynamic_cast< TrsNode * >( nodes [ i ] );
        if ( !p ) {
            cerr << "Error in " << name << ": nodes must be inherited from TrsNode, " << nodes [ i ]->giveName() << " provided" << endl;
            exit(1);
        }
    }
    //check that material is TrsprtMaterial
    VectTrsprtMaterial *p = dynamic_cast< VectTrsprtMaterial * >( mat );
    if ( !p ) {
        cerr << "Error in " << name << ": material must be inherited from VectTrsprtMaterial, " << mat->giveName() << " provided" << endl;
        exit(1);
    }
}

//////////////////////////////////////////////////////////
Matrix DiscreteTrsprtElem :: giveBMatrix(const Point *x) const {
    ( void ) x;
    Matrix B = Matrix :: Zero(1, 2);
    B(0, 0) = -1. / length;
    B(0, 1) = 1. / length;
    return B;
}

//////////////////////////////////////////////////////////
Matrix DiscreteTrsprtElem :: giveHMatrix(const Point *x) const {
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
void DiscreteTrsprtElem :: giveValues(string code, Vector &result) const {
    if ( code.compare("normal") == 0 ) {
        result.resize(ndim);
        for ( unsigned i = 0; i < ndim; i++ ) {
            result [ i ] = normal(i);
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
        Element :: giveValues(code, result);
    }
}

//////////////////////////////////////////////////////////
void DiscreteTrsprtElem :: giveIPValues(std :: string code, unsigned ipnum, Vector &result) const {
    if ( ipnum >= inttype->giveNumIP() ) {
        std :: cerr << name <<  " Error: intergration point number " << ipnum << " exceeds number of integration points" << std :: endl;
        exit(1);
    }
    if ( code.compare("location") == 0 ) {
        result.resize(ndim);
        Point *h = inttype->giveIPLocationPointer(ipnum);
        for ( unsigned k = 0; k < ndim; k++ ) {
            result [ k ] = ( * h ) [ k ];
        }
    } else {
        Element :: giveIPValues(code, ipnum, result);
    }
};

//////////////////////////////////////////////////////////
void DiscreteTrsprtElem :: computeDampingMatrix() {
    dampC = Matrix :: Zero(2, 2);
    double s = area * stats [ 0 ]->giveDampingTensor()(0, 0) * length /  ( 2. * ndim );

    dampC(0, 0) = dampC(1, 1) = s; //finite volume
    if ( BolanderCapacityMatrix ) { //from Bolander's papers
        dampC(0, 0) = dampC(1, 1) = 2. / 3. * s;
        dampC(1, 0) = dampC(0, 1) = s / 3.;
    }
}

//////////////////////////////////////////////////////////
void DiscreteTrsprtElem :: computeMassMatrix() {
    massM = Matrix :: Zero(2, 2);
}

//////////////////////////////////////////////////////////
double DiscreteTrsprtElem :: giveVolumeAssociatedWithNode(unsigned nodenum) const {
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
Vector DiscreteTrsprtElem :: giveStrain(unsigned i, const Vector &DoFs) {
    double averagePressure = ( DoFs [ 0 ] * giveVolumeAssociatedWithNode(0) + DoFs [ 1 ] * giveVolumeAssociatedWithNode(1) ) / volume;
    stats [ 0 ]->setParameterValue("pressure", averagePressure);
    return Element :: giveStrain(i, DoFs);
};

//////////////////////////////////////////////////////////
Matrix DiscreteTrsprtElem :: giveStiffnessMatrix(string matrixType) const {
    return Element :: giveStiffnessMatrix(matrixType) * ndim; //ndim needs to be included here for discrete elements
}


//////////////////////////////////////////////////////////
Vector DiscreteTrsprtElem :: giveInternalForces(const Vector &DoFs, bool frozen, double timeStep) {
    //MyVector Q = Element :: giveInternalForces(DoFs, frozen, timeStep) * ndim; //ndim needs to be included here for discrete elements
    //for (auto p:Q) cout << " " << p;
    return Element :: giveInternalForces(DoFs, frozen, timeStep) * ndim; //ndim needs to be included here for discrete elements
}

//////////////////////////////////////////////////////////
Vector DiscreteTrsprtElem :: integrateLoad(BodyLoad *vl, double time) const {
    return Element :: integrateLoad(vl, time) / ndim;
}

//////////////////////////////////////////////////////////
Vector DiscreteTrsprtElem :: integrateInternalSources() {
    return Element :: integrateInternalSources() / ndim;
}

//////////////////////////////////////////////////////////
Vector DiscreteTrsprtElem :: giveVectorToNode(const unsigned &node_i, const unsigned &ip_id) const {
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
void DiscreteTrsprtElem :: extrapolateIPValuesToNodes(string code, vector< Vector > &result, Vector &weights) const {
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
// 1D HEAT CONDUCTION ELEMENT
DiscreteHeatConductionElem :: DiscreteHeatConductionElem(const unsigned dim) : DiscreteTrsprtElem(dim) {
    name = "Heat Conduction Element 1D";
    physicalFields [ 1 ] = false; //transport
    physicalFields [ 2 ] = true; //heat conduction
}

//////////////////////////////////////////////////////////
void DiscreteHeatConductionElem :: checkNodeAndMaterialType() const {
    for ( unsigned i = 0; i < 2; i++ ) {
        TempNode *p = dynamic_cast< TempNode * >( nodes [ i ] );
        if ( !p ) {
            cerr << "Error in " << name << ": nodes must be inherited from TempNode, " << nodes [ i ]->giveName() << " provided" << endl;
            exit(1);
        }
    }
    //check that material is TrsprtMaterial
    VectHeatConductionMaterial *p = dynamic_cast< VectHeatConductionMaterial * >( mat );
    if ( !p ) {
        cerr << "Error in " << name << ": material must be inherited from HeatConductionMaterial, " << mat->giveName() << " provided" << endl;
        exit(1);
    }
}

//////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////
// 1D TRANSPORT ELEMENT COUPLED WITH MECHANICS
//////////////////////////////////////////////////////////
void DiscreteTrsprtCoupledElem :: init() {
    DiscreteTrsprtElem :: init(); //calling base class method;
}

//////////////////////////////////////////////////////////
void DiscreteTrsprtCoupledElem :: giveValues(string code, Vector &result) const {
    if ( code.compare("numOfFriends") == 0 ) {
        result.resize(0);
        result [ 0 ] = friends.size();
    } else {
        DiscreteTrsprtElem :: giveValues(code, result);
    }
};

//////////////////////////////////////////////////////////
void DiscreteTrsprtCoupledElem :: addNewFriend(RigidBodyContact *f, double weight) {
    friends.push_back(f);
    friendsweight.push_back(weight);
}

//////////////////////////////////////////////////////////
Vector DiscreteTrsprtCoupledElem :: giveStrain(unsigned i, const Vector &DoFs) {
    //crack opening
    double crackInNeighborhood = 0;
    double crackVolume = 0.;
    double elem_crack_opening;
    size_t m = 0;
    Vector res;
    for ( auto &f: friends ) {
        elem_crack_opening = 0.;
        for ( unsigned k = 0; k < f->giveNumIP(); k++ ) {
            f->giveIPValues("crack_opening", k, res);
            if ( res.size() == 1 ) {
                elem_crack_opening += abs(res [ 0 ]);
            }
        }
        crackInNeighborhood += pow(elem_crack_opening / f->giveNumIP(), 3) * friendsweight [ m ];
        if ( ndim == 3 ) {
            crackVolume += ( elem_crack_opening / f->giveNumIP() ) * f->giveArea() * length / f->givePerimeter();
        } else if ( ndim == 2 ) {
            crackVolume += ( elem_crack_opening / f->giveNumIP() ) * f->giveArea();
        }
        m++;
    }
    stats [ 0 ]->setParameterValue("crack_param", crackInNeighborhood);
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


    stats [ 0 ]->setParameterValue("volumetric_strain", volStrain); //trace divided by dimension to obtain mechanical volumetric strain

    return DiscreteTrsprtElem :: giveStrain(i, DoFs);
};

//////////////////////////////////////////////////////////
void DiscreteTrsprtCoupledElem :: collectInformationsFromNeigborhood() {
    findFriendsFromSimplices();
}

//////////////////////////////////////////////////////////
void DiscreteTrsprtCoupledElem :: findFriendsFromSimplices() {
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


//////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////
// RigidBodyContactWithHeatConduction ELEMENT

RigidBodyContactWithHeatConduction :: RigidBodyContactWithHeatConduction(const unsigned dim) : RigidBodyContact(dim) {
    name = "RigidBodyContactWithHeatConduction";
    physicalFields [ 2 ] = true;
}

//////////////////////////////////////////////////////////
void RigidBodyContactWithHeatConduction :: giveValues(string code, Vector &result) const {
    RigidBodyContact :: giveValues(code, result);
}

//////////////////////////////////////////////////////////
void RigidBodyContactWithHeatConduction :: checkNodeType() const {
    //check that nodes are particles
    for ( unsigned i = 0; i < 2; i++ ) {
        ParticleWithTemperature *p = dynamic_cast< ParticleWithTemperature * >( nodes [ i ] );
        if ( !p ) {
            cerr << "Error in " << name << ": nodes must be inherited from ParticleWithTemperature, " << nodes [ i ]->giveName() << " provided" << endl;
            exit(1);
        }
    }
}

//////////////////////////////////////////////////////////

Matrix RigidBodyContactWithHeatConduction :: giveBMatrix(const Point *x) const {
    unsigned p = 3 * ( ndim - 1 );
    Matrix B = Matrix :: Zero(ndim + 1, 2 * p + 2);
    Matrix RB = RigidBodyContact :: giveBMatrix(x);
    //mechanics
    for ( unsigned i = 0; i < p; i++ ) {
        for ( unsigned j = 0; j < ndim; j++ ) {
            B(j, i) = RB(j, i);
            B(j, i + p + 1) = RB(j, i + p);
        }
    }
    //transport
    B(ndim, p) = -1. / length;
    ;
    B(ndim, 2 * p + 1) = 1. / length;
    ;
    return B;
}


//////////////////////////////////////////////////////////
Vector RigidBodyContactWithHeatConduction :: giveStrain(unsigned i, const Vector &DoFs) {
    unsigned p = 3 * ( ndim - 1 );
    stats [ 0 ]->setParameterValue("temperature", ( DoFs [ p ] + DoFs [ 2 * p + 1 ] ) / 2.);
    return RigidBodyContact :: giveStrain(i, DoFs);
}

//////////////////////////////////////////////////////////
void RigidBodyContactWithHeatConduction :: init() {
    RigidBodyContact :: init(); //calling base class method;
}

//////////////////////////////////////////////////////////
Matrix RigidBodyContactWithHeatConduction :: giveHMatrix(const Point *x) const {
    ( void ) x;
    return Matrix :: Zero(12, 12);  // NOTE JK: this should be based on ndim
}

//////////////////////////////////////////////////////////
void RigidBodyContactWithHeatConduction :: computeMassMatrix() {
    massM = Matrix :: Zero( 6 * ( ndim - 1 ), 6 * ( ndim - 1 ) );
}

//////////////////////////////////////////////////////////
Matrix RigidBodyContactWithHeatConduction :: giveStiffnessMatrix(string matrixType) const {
    return Element :: giveStiffnessMatrix(matrixType) * ndim; //ndim needs to be included here for discrete elements
}

//////////////////////////////////////////////////////////
void RigidBodyContactWithHeatConduction :: computeDampingMatrix() {
    dampC = giveStiffnessMatrix("elastic") * 1e-15;           //rough fix of zeros, here can be anything
}
