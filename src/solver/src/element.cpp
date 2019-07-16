#include "element.h"
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
    return 0;
};

//////////////////////////////////////////////////////////
double Element :: giveIPValue(string code, unsigned ipnum) const {
    if ( code.compare("x") == 0 ) {
        return ip_locs [ ipnum ].x;
    } else if ( code.compare("y") == 0 )       {
        return ip_locs [ ipnum ].y;
    } else  {
        return stats [ ipnum ]->giveValue(code);
    }
};

//////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////
// RBSN ELEMENT
RigidBodyContact :: RigidBodyContact(const unsigned dim) {
    ndim = dim;
    nodes.resize(2);
    ip_locs.resize(1);
    stats.resize(1);
    name = "RigidBodyContact";
}

//////////////////////////////////////////////////////////
double RigidBodyContact :: giveValue(string code) const {
    if ( code.compare("damage") == 0 ) {
        DisMechMaterialStatus *dmstats = static_cast< DisMechMaterialStatus * >( stats [ 0 ] );
        return dmstats->giveValue(code);
    } else   {
        return mechanicalElement :: giveValue(code);
    }
}

//////////////////////////////////////////////////////////
double RigidBodyContact :: giveIPValue(string code, unsigned ipnum) const {
  if ( code.compare("normal_x") == 0 )       {
      return normal.getX();
  } else if ( code.compare("normal_y") == 0 )       {
      return normal.getY();
  } else {
    return mechanicalElement :: giveIPValue(code, ipnum);
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
void RigidBodyContact :: init() {
    Element :: init(); //calling base class method;

    //check that nodes are particles
    for ( unsigned i = 0; i < 2; i++ ) {
        Particle *p = dynamic_cast< Particle * >( nodes [ i ] );
        if ( !p ) {
            cerr << "Error in " << name << ": nodes must be inherited from Particle, " << nodes [ i ]->giveName() << " provided" << endl;
            exit(1);
        }
    }
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
    } else   {
        //JM: Coplanarity check for vertices on the face
        //JM: checking coplanarity of every consecutive 4 nodes
        double maxErr = 0.0;
        double currErr = 0.0;
        //
        for (unsigned int i=0; i<vert.size()-3; i++){
          currErr = checkCoplanarity(  vert[i]->givePoint(), vert[i+1]->givePoint(), vert[i+2]->givePoint(), vert[i+3]->givePoint()  );
          if (abs(currErr) > maxErr){ maxErr = abs(currErr); }
        }
        //JM: also checking if the beam midpoint is coplanar with the face
        Point midPoint = (nodes [1]->givePoint() + nodes [0]->givePoint())/2.;
        currErr = checkCoplanarity(  vert[0]->givePoint(), vert[1]->givePoint(), vert[2]->givePoint(), midPoint );
        if (abs(currErr) > maxErr){ maxErr = abs(currErr); }
        //
        if (maxErr > 1e-12){
          cerr << "Vertices are not coplanar!!! Coplanarity error: " << maxErr << endl;
          exit(1);
        }

        //JM: face normal vector made from first 3 vertices
        //JM: coordinate swap for tangential vector according to https://orbit.dtu.dk/files/126824972/onb_frisvad_jgt2012_v2.pdf
        Point n = cross(vert[1]->givePoint()-vert[0]->givePoint() , vert[2]->givePoint()-vert[0]->givePoint());
        n /= n.norm();
        Point t2 ;
        if( fabs (n.x ) > fabs (n.z )) t2 = Point (-n.y , n.x , 0.0f );
        else t2 = Point (0.0f , -n.z , n.y );
        t = cross (t2 , n);
        t /= t.norm();


        //JM: Perpendicularity check of the beam and face directions
        //JM: normal of the face surface taken from first 3 vertices is (B - A) x (C - A)
        //JM: perpendicularity check: cross (beam, face)=>0
        Point prp = (nodes [1]->givePoint() - nodes [0]->givePoint())* t ;
        if (prp.norm() > 1e-8){
          cerr << "Face surface is not perpendicular to beam direction!!! Error: " << prp.norm() << endl;
        //  exit(1);
        }

        //JM: finding position of the SINGLE integration point -> center of gravity of the face polygon
        //JM: average point of the polygon for triangulation
        Point avgPoint = Point(0.0, 0.0, 0.0);
        for (unsigned int i=0; i<vert.size(); i++){ avgPoint += vert[i]->givePoint(); }
        avgPoint /= vert.size();

        //JM: integration point coordinates as an average of CGs of face triangles weighted by areas
        ip_locs[0] = Point(0.0, 0.0, 0.0);
        area = 0.0;
        double ai = 0.0;
        unsigned int j = 0;
        for (unsigned int i=0; i<vert.size(); i++){
          j=i+1;
          if (i==vert.size()-1){ j=0; }
          //triangle area computed as a_i = norm(cross(AB, AC)) / 2
          ai = (cross(vert[i]->givePoint() - avgPoint,   vert[j]->givePoint() - avgPoint) ).norm();
          area += ai;
          //triangle cg_i is an average of simplex vertices, adding to CG coordinates multiplied by a_i weight
          ip_locs[0] += ( avgPoint +vert[i]->givePoint() +vert[j]->givePoint() )/3.0 * ai;
        }
        ip_locs[0] /= area;

        //JM: Check if integration point is coplanar with face
        currErr = checkCoplanarity( vert[0]->givePoint(), vert[1]->givePoint(), vert[2]->givePoint(), ip_locs[0] );
        if (abs(currErr) > 1e-12){
          cerr << "Integration point is not coplanar with the face!!! Coplanarity error: " << currErr << endl;
          exit(1);
        }
    }

    stats [ 0 ] = mat->giveNewMaterialStatus(this);
    normal = nodes [ 1 ]->givePoint() - nodes [ 0 ]->givePoint();
    length = normal.norm();
    normal = normal / length;
    if ( abs(normal * t) > 1e-6 ) {
        cout << vert [ 0 ]->givePoint().x << " " <<  vert [ 0 ]->givePoint().y <<  " X " << vert [ 1 ]->givePoint().x << " " <<  vert [ 1 ]->givePoint().y << endl;
        cout << nodes [ 0 ]->givePoint().x << " " <<  nodes [ 0 ]->givePoint().y <<  " X " << nodes [ 1 ]->givePoint().x << " " <<  nodes [ 1 ]->givePoint().y << endl;
        cerr << "Error: normal and contact vector are not parallel, error " << normal * t << " normal v." << normal.x << " " << normal.y << " contact v. " << t.x << " " << t.y << endl;
        exit(1);
    }

    //Matrices according to habilitation of Jan Elias (2017, page 42): https://www.vutbr.cz/www_base/vutdisk.php?i=103116a130

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
    B = B / length;

    //Matrix R
    Matrix R;
    if ( ndim == 2 ) {
        Point t1 = Point(-normal.y, normal.x);
        R = Matrix(2, 2);
        R [ 0 ] [ 0 ] = normal.x;
        R [ 0 ] [ 1 ] = normal.y;
        R [ 1 ] [ 0 ] = t1.x;
        R [ 1 ] [ 1 ] = t1.y;
    } else if ( ndim == 3 )       {
        Point t1, t2;
        if ( abs(normal.x) > 1e-3 ) {
            t1 = Point(-normal.y / normal.x, 1, 0);
        } else if ( abs(normal.y) > 1e-3 ) {
            t1 = Point(0, -normal.z / normal.y, 1);
        } else                                                                       {
            t1 = Point(1, 0, -normal.x / normal.z);
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
    } else  {
        cerr << "Error - RigidBodyContact: dimension " << ndim << "not implemented" << endl;
        exit(0);
    }
    RB = R * B;
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
    } else if ( ndim == 2 )      {
        A [ 0 ] [ 0 ] = A [ 1 ] [ 1 ] = 1;
        A [ 0 ] [ 2 ] = a.y - x.y;
        A [ 1 ] [ 2 ] = x.x - a.x;
    } else  {
        cerr << "Error - RigidBodyContact: dimension " << ndim << "not implemented" << endl;
        exit(0);
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
    Matrix K = RB.transpose() * Alpha * RB;
    return K * ( length * area );
}

//////////////////////////////////////////////////////////
Matrix RigidBodyContact :: giveInertiaMatrix() const {
    Matrix S(12, 12);
    //TrsprtMaterialStatus *tstats = static_cast<TrsprtMaterialStatus *>(stats[0]);
    //TO BE DONE
    return S;
}

//////////////////////////////////////////////////////////
Vector RigidBodyContact :: giveInternalForces(const Vector &DoFs) const {
    Vector strainNT = RB * DoFs;
    DisMechMaterialStatus *dmstats = static_cast< DisMechMaterialStatus * >( stats [ 0 ] );
    Vector stressNT = dmstats->giveStress(strainNT);
    return RB.transpose() * ( stressNT * ( length * area ) );
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
        if ( num2 < 0 ) {
            bound = true;
            vert [ i ] = NULL;
        } else  {
            vert [ i ] = fullnodes->giveNode(num2);
        }
    }
    iss >> num;
    mat = fullmatrs->giveMaterial(num);
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

    Point v0, v1, t, normal;
    if ( ndim == 2 ) {
        if ( !( vert.size() == 2 ) ) {
            cerr << "Error: exactly 2 vertices must be involved, " << vert.size() << " provided" << endl;
            exit(1);
        }

        if ( bound ) {
            if ( vert [ 0 ] ) {
                v0 = vert [ 0 ]->givePoint();
            } else {
                v0 = ( nodes [ 1 ]->givePoint() + nodes [ 0 ]->givePoint() ) / 2.;
            }
            if ( vert [ 1 ] ) {
                v1 = vert [ 1 ]->givePoint();
            } else {
                v1 = ( nodes [ 1 ]->givePoint() + nodes [ 0 ]->givePoint() ) / 2.;
            }
        } else  {
            v0 = vert [ 0 ]->givePoint();
            v1 = vert [ 1 ]->givePoint();
        }
        ip_locs [ 0 ] = ( v0 + v1 ) / 2.;
        t = v1 - v0;
        area = t.norm();
        t = t / area;
    } else   {
      if ( !( vert.size() == 2 ) ) {
          cerr << "Error: exactly 2 vertices must be involved, " << vert.size() << " provided" << endl;
          exit(1);
      }

      //JM: coplanarity controll of nodes divided by the 1D trsprt link (v0,v1)
      for (unsigned int i = 0; i<vert.size(); i++){
        cout << "Vertex nr." << i << ": " << vert[i]->getPoint() << endl 
      }




      //Work in progress
      cerr << "Dimension " << ndim << " transport implementation is in progress. JM" << endl;

      exit(0);
    }

    stats [ 0 ] = mat->giveNewMaterialStatus(this);
    normal = nodes [ 1 ]->givePoint() - nodes [ 0 ]->givePoint();
    length = normal.norm();
    normal = normal / length;

    if ( abs(normal * t) > 1e-6 ) {
        cout << v0.x << " " <<  v0.y <<  " X " << v1.x << " " <<  v1.y << endl;
        cout << nodes [ 0 ]->givePoint().x << " " <<  nodes [ 0 ]->givePoint().y <<  " X " << nodes [ 1 ]->givePoint().x << " " <<  nodes [ 1 ]->givePoint().y << endl;
        cerr << "Error: normal and contact vector are not parallel, error " << normal * t << " normal v." << normal.x << " " << normal.y << " contact v. " << t.x << " " << t.y << endl;
        exit(1);
    }
}

//////////////////////////////////////////////////////////
Matrix Transp1D :: giveConductivityMatrix(string matrixType) const {
    Matrix C(2, 2);
    TrsprtMaterialStatus *tstats = static_cast< TrsprtMaterialStatus * >( stats [ 0 ] );
    double c = area * tstats->giveConductivity() / length;
    C [ 0 ] [ 0 ] = C [ 1 ] [ 1 ] = c;
    C [ 1 ] [ 0 ] = C [ 0 ] [ 1 ] = -c;
    return C;
}

//////////////////////////////////////////////////////////
Matrix Transp1D :: giveCapacityMatrix() const {
    Matrix S(2, 2);
    TrsprtMaterialStatus *tstats = static_cast< TrsprtMaterialStatus * >( stats [ 0 ] );
    double s = area * tstats->giveCapacity() * length / 6.;
    S [ 0 ] [ 0 ] = S [ 1 ] [ 1 ] = 2 * s;
    S [ 1 ] [ 0 ] = S [ 0 ] [ 1 ] = s;
    return S;
}

//////////////////////////////////////////////////////////
Vector Transp1D :: giveInternalForces(const Vector &DoFs) const {
    return giveConductivityMatrix("elastic") * DoFs;
};
