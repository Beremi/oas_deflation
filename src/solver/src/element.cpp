#include "element.h"
#include "element_container.h"
#include "boundary_condition.h"

//////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////
// BASIC ELEMENT - MASTER CLASS
Element :: ~Element() {
    for ( vector< MaterialStatus * > :: iterator e = stats.begin(); e != stats.end(); ++e ) {
        delete * e;
    }
}

// std :: string Element :: giveLineToSave(NodeContainer * nodes_all) const {
//   std :: string str = this->giveName() + "\t";
//   if ( this->nodes.size() > 2 ){
//     // for polygonal elems etc
//     str += to_string(this->nodes.size()) + "\t";
//   }
//   for ( auto const &node : this->nodes ){
//     str += to_string(nodes_all->giveNodeId(node)) + "\t";
//   }
//   return str;
// }

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

    setIntegrationPointsAndWeights();

    Bs.resize(ip_locs.size() );
    for ( k = 0; k < ip_locs.size(); k++ ) {
        Bs [ k ] = giveBMatrix(& ip_locs [ k ]);
    }
}

//////////////////////////////////////////////////////////
vector< unsigned >Element :: giveDoFsInDirection(unsigned dir) const {
    vector< unsigned >DoFinDir(nodes.size() );
    for ( unsigned i = 0; i < nodes.size(); i++ ) {
        DoFinDir [ i ] = nodes [ i ]->giveStartingDoF() + dir;
    }
    return DoFinDir;
}

//////////////////////////////////////////////////////////
void Element :: initMaterialStatuses() {
    unsigned num = 0;
    for ( vector< MaterialStatus * > :: iterator m = stats.begin(); m != stats.end(); ++m, num++ ) {
        ( * m )->setID(num);
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
Matrix Element :: giveStiffnessMatrix(string matrixType) const {
    unsigned nDoFs = DoFids.size();
    Matrix K(nDoFs, nDoFs);
    Matrix D(0, 0);
    for ( unsigned i = 0; i < stats.size(); i++ ) {
        D = stats [ i ]->giveStiffnessTensor(matrixType, ndim);
        K += Bs [ i ].transpose() * D * ( Bs [ i ] * ip_weights [ i ] );
    }
    return K;
}

//////////////////////////////////////////////////////////
Vector Element :: giveInternalForces(const Vector &DoFs, bool frozen) {
    Vector intF;
    intF.resize(DoFids.size() );
    Vector stress;
    for ( unsigned i = 0; i < stats.size(); i++ ) {
        if ( frozen ) {
            stress = stats [ i ]->giveStressWithFrozenIntVars(giveStrain(i, DoFs));  //frozen internal variables
        } else  {
            stress = stats [ i ]->giveStress(giveStrain(i, DoFs) ); //full evaluation of stress including change of state variables
        }
        intF  += Bs [ i ].transpose() * (  stress * ip_weights [ i ] );
    }
    return intF;
}

//////////////////////////////////////////////////////////
Matrix Element :: giveDampingMatrix() const {
    unsigned nDoFs = DoFids.size();
    Matrix M(nDoFs, nDoFs);
    double c;
    Matrix H;
    for ( unsigned i = 0; i < stats.size(); i++ ) {
        H = giveHMatrix(& ip_locs [ i ]);
        c = stats [ i ]->giveDampingConstant();
        M += matrix_multiply(H.transpose(), H) * ( ip_weights [ i ] * c );
    }
    return M;
}

//////////////////////////////////////////////////////////
Matrix Element :: giveMassMatrix() const {
    unsigned nDoFs = DoFids.size();
    Matrix M(nDoFs, nDoFs);
    double c;
    Matrix H;
    for ( unsigned i = 0; i < stats.size(); i++ ) {
        H = giveHMatrix(& ip_locs [ i ]);
        c = stats [ i ]->giveMassConstant();
        M += matrix_multiply(H.transpose(), H) * ( ip_weights [ i ] * c );
    }
    return M;
}

//////////////////////////////////////////////////////////
vector< double >Element :: integrateLoad(BodyLoad *vl, double time) const {
    unsigned nDoFs = DoFids.size();
    vector< double >load(nDoFs);
    double fvalue;
    unsigned dir = vl->giveDirection();
    Matrix H;
    for ( unsigned i = 0; i < stats.size(); i++ ) {
        H = giveHMatrix(& ip_locs [ i ]);
        fvalue = vl->giveValue(& ip_locs [ i ], time);
        for ( unsigned j = 0; j < nDoFs; j++ ) {
            load [ j ] += H [ dir ] [ j ] * fvalue * ip_weights [ i ];
        }
    }
    return load;
}

//////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////
// RBSN ELEMENT
RigidBodyContact :: RigidBodyContact(const unsigned dim) {
    ndim = dim;
    nodes.resize(2);
    name = "LTCBEAM";
    vtk_cell_type = 3;
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
    Matrix Aa = giveAMatrix(nodes [ 0 ]->givePoint(), ip_locs [ 0 ]) * ( -1. );
    Matrix Ab = giveAMatrix(nodes [ 1 ]->givePoint(), ip_locs [ 0 ]);
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
    ip_locs.resize(1);
    ip_weights.resize(1);
    stats.resize(1);


    Point t;
    if ( ndim == 2 ) {
        // NOTE the following taken from few lines upper
        if ( !( vert.size() == 2 ) ) {
            cerr << "Error: exactly 2 vertices must be involved, " << vert.size() << " provided" << endl;
            exit(1);
        }
        ip_locs [ 0 ] = ( vert [ 0 ]->givePoint() + vert [ 1 ]->givePoint() ) / 2.;
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
        if ( maxErr > 1e-10 ) {
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
            ai = ( cross(vert [ i ]->givePoint() - avgPoint,   vert [ j ]->givePoint() - avgPoint) ).norm()/2.;
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

    ip_weights [ 0 ] = length * area;  // NOTE JK: this works for stiffness matrix, since there is actually A * L, but not for the mass matrix, where volume should be taken (divide by dim) // JE: Agreed, independent integration of mass has been included for this element

    stats [ 0 ] = mat->giveNewMaterialStatus(this);
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
Vector RigidBodyContact :: giveContactStrainNT(const Vector &DoFs) const {
    return Bs [ 0 ] * DoFs;
};

//////////////////////////////////////////////////////////
Vector RigidBodyContact :: giveContactStrainXYZ(const Vector &DoFs) const {
    return this->R.transpose() * this->giveContactStrainNT(DoFs);
};

Vector RigidBodyContact :: giveContactStressXYZ(const Vector &DoFs) {
    return this->R.transpose() * this->giveMatStatus(0)->giveStressWithFrozenIntVars(this->giveContactStrainNT(DoFs));
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
    Point distance = this->ip_locs [ ip_id ] - this->nodes [ node_i ]->givePoint();
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
    Matrix Aa = giveAMatrix(nodes [ 0 ]->givePoint(), ip_locs [ 0 ]) * ( -1. );
    Matrix Ab = giveAMatrix(nodes [ 1 ]->givePoint(), ip_locs [ 0 ]);
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
    ip_locs.resize(1);
    ip_weights.resize(1);
    stats.resize(1);

    ip_locs [ 0 ] = ( nodes [ 0 ]->givePoint() + nodes [ 1 ]->givePoint() ) / 2.;

    normal = nodes [ 1 ]->givePoint() - nodes [ 0 ]->givePoint();
    length = normal.norm();

    Point t;
    if ( ndim == 2 ) {
        if ( !( vert.size() == 2 ) ) {
            cerr << "Error: exactly 2 vertices must be involved, " << vert.size() << " provided" << endl;
            exit(1);
        }

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
            ai = ( cross(vert [ i ]->givePoint() - avgPoint,   vert [ j ]->givePoint() - avgPoint) ).norm() / 2.;
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



    normal = normal / length;
    if ( abs(normal * t) > 1e-5 ) {
        cout << vert [ 0 ]->givePoint().x << " " <<  vert [ 0 ]->givePoint().y <<  " X " << vert [ 1 ]->givePoint().x << " " <<  vert [ 1 ]->givePoint().y << endl;
        cout << nodes [ 0 ]->givePoint().x << " " <<  nodes [ 0 ]->givePoint().y <<  " X " << nodes [ 1 ]->givePoint().x << " " <<  nodes [ 1 ]->givePoint().y << endl;
        cerr << "TRSPRT: normal and contact vector are not parallel, error " << normal * t << endl;
        cout << " normal v.:"; 
        for(unsigned p=0; p<ndim; p++ ) cout << "\t" << normal.giveCoord(p);
        cout << endl;
        cout << " contact v.:"; 
        for(unsigned p=0; p<ndim; p++ ) cout << "\t" << t.giveCoord(p);
        cout << endl;
        exit(1);
    }

    ip_weights [ 0 ] = length * area;  // NOTE JK: should not this be divided by dimension? Otherwise you integrate dim-times volume that is actually there (this is what caused problems in study on unstructured grid for FDM contribution)
    stats [ 0 ] = mat->giveNewMaterialStatus(this);
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
    Vector pressureGradPlain = Element :: giveStrain(i, DoFs);

    Vector strain(2);
    strain [ 0 ] = pressureGradPlain [ 0 ];
    strain [ 1 ] = (DoFs [ 0 ] + DoFs [ 1 ] ) / 2.;

    return strain;
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
    //ip point is centroid
    Point centroid = ip_locs [ 0 ];
    unsigned i = vert.size() - 1;
    for ( unsigned j = 0; j < vert.size(); i = j, j++ ) {   //circle around the face
        for ( unsigned k = 0; k < elemcont->giveSize(); k++ ) {
            RigidBodyContact *rbc = dynamic_cast< RigidBodyContact * >( elemcont->giveElement(k) );
            if ( rbc ) {
                if ( ( rbc->giveNode(0) == vert [ i ] && rbc->giveNode(1) == vert [ j ] ) || ( rbc->giveNode(0) == vert [ j ] && rbc->giveNode(1) == vert [ i ] ) ) {
                    friends.push_back(rbc);
                    friendsweight.push_back( ( centroid - ( vert [ i ]->givePoint() + vert [ j ]->givePoint() ) / 2. ).norm() );
                    break; //only one friend avalilable in 3D for given pair
                }
            }
        }
    }
}

//////////////////////////////////////////////////////////
double Transp1DCoupled :: giveValue(string code) const {
    if ( code.compare("numOfFriends") == 0 ) {
        return friends.size();
    } else  {
        return Transp1D :: giveValue(code);
    }
};

//////////////////////////////////////////////////////////
double Transp1DCoupled :: giveIPValue(string code, unsigned ipnum) const {
    if ( code.compare("numOfFriends") == 0 ) {
        return friends.size();
    } else  {
        return Transp1D :: giveIPValue(code, ipnum);
    }
}

//////////////////////////////////////////////////////////
void Transp1DCoupled :: findElementFriends(ElementContainer *elemcont) {
    if ( ndim == 2 ) {
        findFriends2D(elemcont);
    } else if ( ndim == 3 ) {
        findFriends3D(elemcont);
    }

    //collect nodes from friend elements
}

//////////////////////////////////////////////////////////
Vector Transp1DCoupled :: giveStrain(unsigned i, const Vector &DoFs) {
    Vector pressureGradPlain = Transp1D :: giveStrain(i, DoFs);

    Vector pressureGrad(2 + 2 * friends.size() );
    pressureGrad [ 0 ] = pressureGradPlain [ 0 ];
    pressureGrad [ 1 ] = pressureGradPlain [ 1 ];
    pressureGrad [ 2 ] = area;

    double elem_crack_opening;
    size_t m = 0;
    for ( auto &f: friends ) {
        elem_crack_opening = 0;
        for ( unsigned k = 0; k < f->giveIPNum(); k++ ) {
            elem_crack_opening += abs(f->giveIPValue("tempCrackOpening", k) );
        }
        pressureGrad [ 2 * m + 3 ] += elem_crack_opening / f->giveIPNum(); //average crack opening in friend mechanical element
        pressureGrad [ 2 * m + 4 ] = friendsweight [ m ]; //crack length in friend mechanical element
        m++;
    }

    return pressureGrad;
};

//////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////
// 2D QUADRILATERAL TRANSPORT ELEMENT
TrsprtQuad :: TrsprtQuad() {
    ndim = 2;
    name = "TrsprtQuad";
    vtk_cell_type = 9;
}

//////////////////////////////////////////////////////////
void TrsprtQuad :: setIntegrationPointsAndWeights() {
    unsigned nnodes = nodes.size();
    ip_locs.resize(nnodes);
    ip_weights.resize(nnodes);
    stats.resize(nnodes);
    double q = 1. / pow(3., 0.5);
    ip_locs [ 0 ] = Point( -q,  -q);
    ip_locs [ 1 ] = Point(  q,  -q);
    ip_locs [ 2 ] = Point(  q,   q);
    ip_locs [ 3 ] = Point( -q,   q);
    Matrix phiGrad(ndim, 4);
    for ( unsigned k = 0; k < nnodes; k++ ) {
        stats [ k ] = mat->giveNewMaterialStatus(this);
        ip_weights [ k ] = shapeFGrad(& ( ip_locs [ k ] ), phiGrad);
    }
};

//////////////////////////////////////////////////////////
void TrsprtQuad :: readFromLine(istringstream &iss, NodeContainer *fullnodes, MaterialContainer *fullmatrs) {
    unsigned num;
    nodes.resize(4);
    for ( unsigned k = 0; k < 4; k++ ) {
        iss >> num;
        nodes [ k ] = fullnodes->giveNode(num);
    }
    iss >> num;
    mat = fullmatrs->giveMaterial(num);
}

//////////////////////////////////////////////////////////
void TrsprtQuad :: shapeF(const Point *x, Vector &phi) const {
    //x in natural coordinates
    phi [ 0 ] = 0.25 * ( 1. - x->getX() ) * ( 1. - x->getY() );
    phi [ 1 ] = 0.25 * ( 1. + x->getX() ) * ( 1. - x->getY() );
    phi [ 2 ] = 0.25 * ( 1. + x->getX() ) * ( 1. + x->getY() );
    phi [ 3 ] = 0.25 * ( 1. - x->getX() ) * ( 1. + x->getY() );
}

//////////////////////////////////////////////////////////
double TrsprtQuad :: shapeFGrad(const Point *x, Matrix &phiGrad) const {
    //x in natural coordinates
    phiGrad [ 0 ] [ 0 ] = -0.25 * ( 1. - x->getY() );
    phiGrad [ 0 ] [ 1 ] = -phiGrad [ 0 ] [ 0 ];
    phiGrad [ 0 ] [ 2 ] = 0.25 * ( 1. + x->getY() );
    phiGrad [ 0 ] [ 3 ] = -phiGrad [ 0 ] [ 2 ];
    phiGrad [ 1 ] [ 0 ] = -0.25 * ( 1. - x->getX() );
    phiGrad [ 1 ] [ 1 ] = -0.25 * ( 1. + x->getX() );
    phiGrad [ 1 ] [ 2 ] = -phiGrad [ 1 ] [ 1 ];
    phiGrad [ 1 ] [ 3 ] = -phiGrad [ 1 ] [ 0 ];

    //Jacobi Matrix
    Matrix Jac(2, 2);
    Point n;
    for ( unsigned i = 0; i < 4; i++ ) {
        n = nodes [ i ]->givePoint();
        Jac [ 0 ] [ 0 ] += phiGrad [ 0 ] [ i ] * n.getX();
        Jac [ 0 ] [ 1 ] += phiGrad [ 0 ] [ i ] * n.getY();
        Jac [ 1 ] [ 0 ] += phiGrad [ 1 ] [ i ] * n.getX();
        Jac [ 1 ] [ 1 ] += phiGrad [ 1 ] [ i ] * n.getY();
    }
    double JacDet = Jac [ 0 ] [ 0 ] * Jac [ 1 ] [ 1 ] - Jac [ 0 ] [ 1 ] * Jac [ 1 ] [ 0 ];
    Matrix JacInv(2, 2);
    JacInv [ 0 ] [ 0 ] = Jac [ 1 ] [ 1 ] / JacDet;
    JacInv [ 0 ] [ 1 ] = -Jac [ 0 ] [ 1 ] / JacDet;
    JacInv [ 1 ] [ 0 ] = -Jac [ 1 ] [ 0 ] / JacDet;
    JacInv [ 1 ] [ 1 ] = Jac [ 0 ] [ 0 ] / JacDet;

    //transorm to physical space
    phiGrad = JacInv * phiGrad;
    return JacDet;
}

//////////////////////////////////////////////////////////
Matrix TrsprtQuad :: giveBMatrix(const Point *x) const {
    Matrix phiG(ndim, nodes.size() );
    shapeFGrad(x, phiG);
    return phiG;
}

//////////////////////////////////////////////////////////
Matrix TrsprtQuad :: giveHMatrix(const Point *x) const {
    Vector phi(DoFids.size() );
    shapeF(x, phi);
    Matrix H(1,DoFids.size());
    for(unsigned k=0; k<DoFids.size(); k++) H[0][k] = phi[k];
    return H;
}


//////////////////////////////////////////////////////////
Vector TrsprtQuad :: giveStrain(unsigned i, const Vector &DoFs){
    Vector pressureGradPlain = Element :: giveStrain(i, DoFs);

    Vector strain(pressureGradPlain.size()+1);
    for(unsigned k=0; k<pressureGradPlain.size(); k++) {
        strain [ k ] = pressureGradPlain [ k ];
    }

    //evaluate pressure at gauss point to account for nonlinearity
    strain[pressureGradPlain.size()] = matrix_vector_multiply(giveHMatrix(&ip_locs [ i ]), DoFs)[0];

    return strain;
}

//////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////
// 3D BRICK TRANSPORT ELEMENT
TrsprtBrick :: TrsprtBrick() {
    ndim = 3;
    name = "TrsprtBrick";
    vtk_cell_type = 12;
}

//////////////////////////////////////////////////////////
void TrsprtBrick :: setIntegrationPointsAndWeights() {
    unsigned nnodes = nodes.size();
    ip_locs.resize(nnodes);
    ip_weights.resize(nnodes);
    stats.resize(nnodes);
    double q = 1. / pow(3., 0.5);
    ip_locs [ 0 ] = Point( -q,  -q, -q);
    ip_locs [ 1 ] = Point(  q,  -q, -q);
    ip_locs [ 2 ] = Point(  q,   q, -q);
    ip_locs [ 3 ] = Point( -q,   q, -q);
    ip_locs [ 4 ] = Point( -q,  -q,  q);
    ip_locs [ 5 ] = Point(  q,  -q,  q);
    ip_locs [ 6 ] = Point(  q,   q,  q);
    ip_locs [ 7 ] = Point( -q,   q,  q);
    Matrix phiGrad(ndim, 8);
    for ( unsigned k = 0; k < nnodes; k++ ) {
        stats [ k ] = mat->giveNewMaterialStatus(this);
        ip_weights [ k ] = shapeFGrad(& ( ip_locs [ k ] ), phiGrad);
    }
};

//////////////////////////////////////////////////////////
void TrsprtBrick :: readFromLine(istringstream &iss, NodeContainer *fullnodes, MaterialContainer *fullmatrs) {
    unsigned num;
    nodes.resize(8);
    for ( unsigned k = 0; k < 8; k++ ) {
        iss >> num;
        nodes [ k ] = fullnodes->giveNode(num);
    }
    iss >> num;
    mat = fullmatrs->giveMaterial(num);
}

//////////////////////////////////////////////////////////
void TrsprtBrick :: shapeF(const Point *x, Vector &phi) const {
    //x in natural coordinates
    phi [ 0 ] = 0.125 * ( 1. - x->getX() ) * ( 1. - x->getY() )  * ( 1. - x->getZ() );
    phi [ 1 ] = 0.125 * ( 1. + x->getX() ) * ( 1. - x->getY() )  * ( 1. - x->getZ() );
    phi [ 2 ] = 0.125 * ( 1. + x->getX() ) * ( 1. + x->getY() )  * ( 1. - x->getZ() );
    phi [ 3 ] = 0.125 * ( 1. - x->getX() ) * ( 1. + x->getY() )  * ( 1. - x->getZ() );
    phi [ 4 ] = 0.125 * ( 1. - x->getX() ) * ( 1. - x->getY() )  * ( 1. + x->getZ() );
    phi [ 5 ] = 0.125 * ( 1. + x->getX() ) * ( 1. - x->getY() )  * ( 1. + x->getZ() );
    phi [ 6 ] = 0.125 * ( 1. + x->getX() ) * ( 1. + x->getY() )  * ( 1. + x->getZ() );
    phi [ 7 ] = 0.125 * ( 1. - x->getX() ) * ( 1. + x->getY() )  * ( 1. + x->getZ() );
}

//////////////////////////////////////////////////////////
double TrsprtBrick :: shapeFGrad(const Point *x, Matrix &phiGrad) const {
    //x in natural coordinates
    phiGrad [ 0 ] [ 0 ] = -0.125 * ( 1. - x->getY() )  * ( 1. - x->getZ() );
    phiGrad [ 0 ] [ 1 ] = -phiGrad [ 0 ] [ 0 ];
    phiGrad [ 0 ] [ 2 ] =  0.125 * ( 1. + x->getY() )  * ( 1. - x->getZ() );
    phiGrad [ 0 ] [ 3 ] = -phiGrad [ 0 ] [ 2 ];
    phiGrad [ 0 ] [ 4 ] = -0.125 * ( 1. - x->getY() )  * ( 1. + x->getZ() );
    phiGrad [ 0 ] [ 5 ] = -phiGrad [ 0 ] [ 4 ];
    phiGrad [ 0 ] [ 6 ] =  0.125 * ( 1. + x->getY() )  * ( 1. + x->getZ() );
    phiGrad [ 0 ] [ 7 ] = -phiGrad [ 0 ] [ 6 ];
    phiGrad [ 1 ] [ 0 ] = -0.125 * ( 1. - x->getX() )  * ( 1. - x->getZ() );
    phiGrad [ 1 ] [ 1 ] = -0.125 * ( 1. + x->getX() )  * ( 1. - x->getZ() );
    phiGrad [ 1 ] [ 2 ] = -phiGrad [ 1 ] [ 1 ];
    phiGrad [ 1 ] [ 3 ] = -phiGrad [ 1 ] [ 0 ];
    phiGrad [ 1 ] [ 4 ] = -0.125 * ( 1. - x->getX() )  * ( 1. + x->getZ() );
    phiGrad [ 1 ] [ 5 ] = -0.125 * ( 1. + x->getX() )  * ( 1. + x->getZ() );
    phiGrad [ 1 ] [ 6 ] = -phiGrad [ 1 ] [ 5 ];
    phiGrad [ 1 ] [ 7 ] = -phiGrad [ 1 ] [ 4 ];
    phiGrad [ 2 ] [ 0 ] = -0.125 * ( 1. - x->getX() )  * ( 1. - x->getY() );
    phiGrad [ 2 ] [ 1 ] = -0.125 * ( 1. + x->getX() )  * ( 1. - x->getY() );
    phiGrad [ 2 ] [ 2 ] = -0.125 * ( 1. + x->getX() )  * ( 1. + x->getY() );
    phiGrad [ 2 ] [ 3 ] = -0.125 * ( 1. - x->getX() )  * ( 1. + x->getY() );
    phiGrad [ 2 ] [ 4 ] = -phiGrad [ 2 ] [ 0 ];
    phiGrad [ 2 ] [ 5 ] = -phiGrad [ 2 ] [ 1 ];
    phiGrad [ 2 ] [ 6 ] = -phiGrad [ 2 ] [ 2 ];
    phiGrad [ 2 ] [ 7 ] = -phiGrad [ 2 ] [ 3 ];


    //Jacobi Matrix
    Matrix Jac(3, 3);
    Point n;
    for ( unsigned i = 0; i < 8; i++ ) {
        n = nodes [ i ]->givePoint();
        Jac [ 0 ] [ 0 ] += phiGrad [ 0 ] [ i ] * n.getX();
        Jac [ 0 ] [ 1 ] += phiGrad [ 0 ] [ i ] * n.getY();
        Jac [ 0 ] [ 2 ] += phiGrad [ 0 ] [ i ] * n.getZ();
        Jac [ 1 ] [ 0 ] += phiGrad [ 1 ] [ i ] * n.getX();
        Jac [ 1 ] [ 1 ] += phiGrad [ 1 ] [ i ] * n.getY();
        Jac [ 1 ] [ 2 ] += phiGrad [ 1 ] [ i ] * n.getZ();
        Jac [ 2 ] [ 0 ] += phiGrad [ 2 ] [ i ] * n.getX();
        Jac [ 2 ] [ 1 ] += phiGrad [ 2 ] [ i ] * n.getY();
        Jac [ 2 ] [ 2 ] += phiGrad [ 2 ] [ i ] * n.getZ();
    }
    //detrminant and inverse of 3x3 matrix
    double JacDet = Jac [ 0 ] [ 0 ] * ( Jac [ 1 ] [ 1 ] * Jac [ 2 ] [ 2 ] - Jac [ 2 ] [ 1 ] * Jac [ 1 ] [ 2 ]) 
                 -  Jac [ 0 ] [ 1 ] * ( Jac [ 1 ] [ 0 ] * Jac [ 2 ] [ 2 ] - Jac [ 2 ] [ 0 ] * Jac [ 1 ] [ 2 ]) 
                 +  Jac [ 0 ] [ 2 ] * ( Jac [ 1 ] [ 0 ] * Jac [ 2 ] [ 1 ] - Jac [ 2 ] [ 0 ] * Jac [ 1 ] [ 1 ]);

    Matrix JacInv(3, 3);
    JacInv [ 0 ] [ 0 ] =  (Jac [ 1 ] [ 1 ] * Jac [ 2 ] [ 2 ] - Jac [ 2 ] [ 1 ] * Jac [ 1 ] [ 2 ]) / JacDet;
    JacInv [ 0 ] [ 1 ] = -(Jac [ 0 ] [ 1 ] * Jac [ 2 ] [ 2 ] - Jac [ 2 ] [ 1 ] * Jac [ 0 ] [ 2 ]) / JacDet;
    JacInv [ 0 ] [ 2 ] =  (Jac [ 0 ] [ 1 ] * Jac [ 1 ] [ 2 ] - Jac [ 1 ] [ 1 ] * Jac [ 0 ] [ 2 ]) / JacDet;

    JacInv [ 1 ] [ 0 ] = -(Jac [ 1 ] [ 0 ] * Jac [ 2 ] [ 2 ] - Jac [ 2 ] [ 0 ] * Jac [ 1 ] [ 2 ]) / JacDet;
    JacInv [ 1 ] [ 1 ] =  (Jac [ 0 ] [ 0 ] * Jac [ 2 ] [ 2 ] - Jac [ 2 ] [ 0 ] * Jac [ 0 ] [ 2 ]) / JacDet;
    JacInv [ 1 ] [ 2 ] = -(Jac [ 0 ] [ 0 ] * Jac [ 1 ] [ 2 ] - Jac [ 1 ] [ 0 ] * Jac [ 0 ] [ 2 ]) / JacDet;

    JacInv [ 2 ] [ 0 ] =  (Jac [ 1 ] [ 0 ] * Jac [ 2 ] [ 1 ] - Jac [ 2 ] [ 0 ] * Jac [ 1 ] [ 1 ]) / JacDet;
    JacInv [ 2 ] [ 1 ] = -(Jac [ 0 ] [ 0 ] * Jac [ 2 ] [ 1 ] - Jac [ 2 ] [ 0 ] * Jac [ 0 ] [ 1 ]) / JacDet;
    JacInv [ 2 ] [ 2 ] =  (Jac [ 0 ] [ 0 ] * Jac [ 1 ] [ 1 ] - Jac [ 1 ] [ 0 ] * Jac [ 0 ] [ 1 ]) / JacDet;

    //transorm to physical space
    phiGrad = JacInv * phiGrad;
    return JacDet;
}

//////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////
// 2D QUADRILATERAL MECHANICAL ELEMENT
MechanicalQuad :: MechanicalQuad() {
    ndim = 2;
    name = "MechanicalQuad";
    vtk_cell_type = 9;
}

//////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////
void MechanicalQuad :: setIntegrationPointsAndWeights() {
    unsigned nnodes = nodes.size();
    ip_locs.resize(nnodes);
    ip_weights.resize(nnodes);
    stats.resize(nnodes);
    double q = 1. / pow(3., 0.5);
    ip_locs [ 0 ] = Point(-q, -q);
    ip_locs [ 1 ] = Point(q, -q);
    ip_locs [ 2 ] = Point(q, q);
    ip_locs [ 3 ] = Point(-q, q);
    Matrix phiGrad(ndim, 4);
    for ( unsigned k = 0; k < nnodes; k++ ) {
        stats [ k ] = mat->giveNewMaterialStatus(this);
        ip_weights [ k ] = shapeFGrad(& ( ip_locs [ k ] ), phiGrad);
    }
};

//////////////////////////////////////////////////////////
void MechanicalQuad :: readFromLine(istringstream &iss, NodeContainer *fullnodes, MaterialContainer *fullmatrs) {
    unsigned num;
    nodes.resize(4);
    for ( unsigned k = 0; k < 4; k++ ) {
        iss >> num;
        nodes [ k ] = fullnodes->giveNode(num);
    }
    iss >> num;
    mat = fullmatrs->giveMaterial(num);
}

//////////////////////////////////////////////////////////
void MechanicalQuad :: shapeF(const Point *x, Vector &phi) const {
    //x in natural coordinates
    phi [ 0 ] = 0.25 * ( 1. + x->getX() ) * ( 1. + x->getY() );
    phi [ 1 ] = 0.25 * ( 1. - x->getX() ) * ( 1. + x->getY() );
    phi [ 2 ] = 0.25 * ( 1. - x->getX() ) * ( 1. - x->getY() );
    phi [ 3 ] = 0.25 * ( 1. + x->getX() ) * ( 1. - x->getY() );
}

//////////////////////////////////////////////////////////
double MechanicalQuad :: shapeFGrad(const Point *x, Matrix &phiGrad) const {
    //x in natural coordinates
    phiGrad [ 0 ] [ 0 ] = 0.25 * ( 1. + x->getY() );
    phiGrad [ 0 ] [ 1 ] = -phiGrad [ 0 ] [ 0 ];
    phiGrad [ 0 ] [ 2 ] = -0.25 * ( 1. - x->getY() );
    phiGrad [ 0 ] [ 3 ] = -phiGrad [ 0 ] [ 2 ];
    phiGrad [ 1 ] [ 0 ] = 0.25 * ( 1. + x->getX() );
    phiGrad [ 1 ] [ 1 ] = 0.25 * ( 1. - x->getX() );
    phiGrad [ 1 ] [ 2 ] = -phiGrad [ 1 ] [ 1 ];
    phiGrad [ 1 ] [ 3 ] = -phiGrad [ 1 ] [ 0 ];

    //Jacobi Matrix
    Matrix Jac(2, 2);
    Point n;
    for ( unsigned i = 0; i < 4; i++ ) {
        n = nodes [ i ]->givePoint();
        Jac [ 0 ] [ 0 ] += phiGrad [ 0 ] [ i ] * n.getX();
        Jac [ 0 ] [ 1 ] += phiGrad [ 0 ] [ i ] * n.getY();
        Jac [ 1 ] [ 0 ] += phiGrad [ 1 ] [ i ] * n.getX();
        Jac [ 1 ] [ 1 ] += phiGrad [ 1 ] [ i ] * n.getY();
    }
    double JacDet = Jac [ 0 ] [ 0 ] * Jac [ 1 ] [ 1 ] - Jac [ 0 ] [ 1 ] * Jac [ 1 ] [ 0 ];
    Matrix JacInv(2, 2);
    JacInv [ 0 ] [ 0 ] = Jac [ 1 ] [ 1 ] / JacDet;
    JacInv [ 0 ] [ 1 ] = -Jac [ 0 ] [ 1 ] / JacDet;
    JacInv [ 1 ] [ 0 ] = -Jac [ 1 ] [ 0 ] / JacDet;
    JacInv [ 1 ] [ 1 ] = Jac [ 0 ] [ 0 ] / JacDet;

    //transorm to physical space
    phiGrad = JacInv * phiGrad;

    return JacDet;
}

//////////////////////////////////////////////////////////
Matrix MechanicalQuad :: giveBMatrix(const Point *x) const {
    Matrix phiG(ndim, nodes.size() );
    shapeFGrad(x, phiG);
    Matrix B(3, DoFids.size() );
    for ( unsigned i = 0; i < nodes.size(); i++ ) {
        B [ 0 ] [ 2 * i ]     =   B [ 2 ] [ 2 * i + 1 ] =   phiG [ 0 ] [ i ];
        B [ 1 ] [ 2 * i ]     =   B [ 0 ] [ 2 * i + 1 ] =   0.;
        B [ 1 ] [ 2 * i + 1 ]   =   B [ 2 ] [ 2 * i ]   =   phiG [ 1 ] [ i ];
    }
    return B;
}

//////////////////////////////////////////////////////////
Matrix MechanicalQuad :: giveHMatrix(const Point *x) const {
    Vector phi(nodes.size() );
    shapeF(x, phi);
    Matrix H(ndim,DoFids.size());
    for(unsigned i=0; i<ndim; i++){
        for(unsigned j=0; j<4; j++){ 
            H[i][ndim*j+i] = phi[j];
        }       
    }
    return H;
}

//////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////
// 2D QUADRILATERAL COSSERAT MECHANICAL ELEMENT
CosseratQuad :: CosseratQuad() {
    ndim = 2;
    name = "CosseratQuad";
}


//////////////////////////////////////////////////////////
Matrix CosseratQuad :: giveBMatrix(const Point *x) const {
    Matrix phiG(ndim, nodes.size() );
    shapeFGrad(x, phiG);
    Vector phi(nodes.size() );
    shapeF(x, phi);
    Matrix B(6, DoFids.size() );
    for ( unsigned i = 0; i < nodes.size(); i++ ) {
        B [ 0 ] [ 3 * i ]       =   B [ 2 ] [ 3 * i + 1 ] =   B [ 4 ] [ 3 * i + 2 ]   =   phiG [ 0 ] [ i ];
        B [ 1 ] [ 3 * i + 1 ]   =   B [ 3 ] [ 3 * i ]   =   B [ 5 ] [ 3 * i + 2 ]   =   phiG [ 1 ] [ i ];
        B [ 2 ] [ 3 * i + 2 ]   =   -phi [ i ];
        B [ 3 ] [ 3 * i + 2 ]   =   phi [ i ];
    }
    return B;
}

//////////////////////////////////////////////////////////
Matrix CosseratQuad :: giveHMatrix(const Point *x) const {
    Vector phi(nodes.size() );
    shapeF(x, phi);
    return Matrix(0,0);
}
