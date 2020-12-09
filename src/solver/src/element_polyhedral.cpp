#include "element_polyhedral.h"

//////////////////////////////////////////////////////////
double triArea2D(const Point a, const Point b, const Point c) { //points in counter clockwise direction
    return 0.5 * ( a.getX() * ( b.getY() - c.getY() ) + b.getX() * ( c.getY() - a.getY() ) + c.getX() * ( a.getY() - b.getY() ) );
}

//////////////////////////////////////////////////////////
Point triNormal3D(const Point *a, const Point *b, const Point *c) { //points
    Point normal = Point(0, 0, 0);
    normal.setX( ( b->getY() - a->getY() ) * ( c->getZ() - a->getZ() ) - ( b->getZ() - a->getZ() ) * ( c->getY() - a->getY() ) );
    normal.setY( ( b->getZ() - a->getZ() ) * ( c->getX() - a->getX() ) - ( b->getX() - a->getX() ) * ( c->getZ() - a->getZ() ) );
    normal.setZ( ( b->getX() - a->getX() ) * ( c->getY() - a->getY() ) - ( b->getY() - a->getY() ) * ( c->getX() - a->getX() ) );
    return normal / normal.norm();
}

//////////////////////////////////////////////////////////
double triArea3D(const Point *a, const Point *b, const Point *c) { //points
    return abs(0.5 * pow(pow( ( b->getY() - a->getY() ) * ( c->getZ() - a->getZ() ) - ( b->getZ() - a->getZ() ) * ( c->getY() - a->getY() ), 2 ) + pow( ( b->getZ() - a->getZ() ) * ( c->getX() - a->getX() ) - ( b->getX() - a->getX() ) * ( c->getZ() - a->getZ() ), 2 ) + pow( ( b->getX() - a->getX() ) * ( c->getY() - a->getY() ) - ( b->getY() - a->getY() ) * ( c->getX() - a->getX() ), 2 ), 0.5) );
}

//////////////////////////////////////////////////////////
double tetVolume3D(const Point *a, const Point *b, const Point *c, const Point *d) {
    return ( ( ( b->getY() - a->getY() ) * ( c->getZ() - a->getZ() ) - ( b->getZ() - a->getZ() ) * ( c->getY() - a->getY() ) ) * ( d->getX() - a->getX() ) + ( ( b->getZ() - a->getZ() ) * ( c->getX() - a->getX() ) - ( b->getX() - a->getX() ) * ( c->getZ() - a->getZ() ) ) * ( d->getY() - a->getY() ) + ( ( b->getX() - a->getX() ) * ( c->getY() - a->getY() ) - ( b->getY() - a->getY() ) * ( c->getX() - a->getX() ) ) * ( d->getZ() - a->getZ() ) ) / 6.;
}

//////////////////////////////////////////////////////////
// 2 D
//////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////
// TRANSPORT POLYHEDRAL ELEMENT
TranspPolygonal :: TranspPolygonal(const unsigned dim) {
    ndim = dim;
    name = "TranspPolygonal";
    ip_type = "quad";
}

//////////////////////////////////////////////////////////
void TranspPolygonal :: readFromLine(istringstream &iss, NodeContainer *fullnodes, MaterialContainer *fullmatrs) {
    unsigned num, num2;

    iss >> nnodes;
    nodes.resize(nnodes);
    for ( unsigned i = 0; i < nnodes; i++ ) {
        iss >> num2;
        nodes [ i ] = fullnodes->giveNode(num2);
    }
    iss >> num;
    mat = fullmatrs->giveMaterial(num);
}

//////////////////////////////////////////////////////////
void TranspPolygonal :: sort2D() {
    //estimate centroid
    Point cpoint;
    for ( const auto n: nodes ) {
        cpoint += n->givePoint();
    }
    cpoint /= nnodes;

    //compute angles
    vector< pair< double, unsigned > >angles;
    angles.resize(nnodes);
    for ( unsigned i = 0; i < nnodes; i++ ) {
        angles [ i ].second = i;
        angles [ i ].first = atan2(nodes [ i ]->givePoint().getY() - cpoint.getY(), nodes [ i ]->givePoint().getX() - cpoint.getX() );
    }

    //sort to have counterclockwise direction
    sort(angles.begin(), angles.end() );
    vector< Node * >newnodes;
    newnodes.resize(nnodes);
    for ( unsigned i = 0; i < nnodes; i++ ) {
        newnodes [ i ] = nodes [ angles [ i ].second ];
    }
    nodes = newnodes;


    volume = 0;
    double triVolume;
    faces.resize(nnodes);
    normals.resize(nnodes);
    surfaces.resize(nnodes);
    centroid = Point(0., 0., 0.);
    Point diff;
    //face X must start at node X and end at node X+1
    for ( unsigned i = 0; i < nnodes; i++ ) {
        faces [ i ].resize(2);
        faces [ i ] [ 0 ] = i;
        if(i==nnodes-1) faces [ i ] [ 1 ] = 0;
        else faces [ i ] [ 1 ] = i+1;
        diff = nodes [ faces [ i ] [ 1 ] ]->givePoint() - nodes [ faces [ i ] [ 0 ] ]->givePoint();
        surfaces [ i ] = diff.norm();
        normals [ i ] = Point(diff.y / surfaces [ i ], -diff.x / surfaces [ i ], 0);
        triVolume = triArea2D(nodes [ faces [ i ] [ 0 ] ]->givePoint(), nodes [ faces [ i ] [ 1 ] ]->givePoint(), cpoint);
        centroid += ( nodes [ faces [ i ] [ 0 ] ]->givePoint() + nodes [ faces [ i ] [ 1 ] ]->givePoint() + cpoint ) * triVolume;
        volume += triVolume;
    }
    centroid /= volume * 3.;
}

//////////////////////////////////////////////////////////
void TranspPolygonal :: WachspressShapeF(const Point *x, Vector &phi) const {
    Vector h(nnodes);
    Point oldNormal = normals [ nnodes - 1 ];
    for ( unsigned i = 0; i < nnodes; i++ ) {
        h [ i ] = abs(dot(nodes [ faces [ i ] [ 0 ] ]->givePoint() - *x, normals [ i ]) );
        phi [ i ] = abs(oldNormal.x * normals [ i ].y - oldNormal.y * normals [ i ].x);
        oldNormal = normals [ i ];
    }
    double oldh = h [ nnodes - 1 ];
    double sumW = 0;
    for ( unsigned i = 0; i < nnodes; i++ ) {
        phi [ i ] /= oldh * h [ i ];
        sumW += phi [ i ];
        oldh = h [ i ];
    }
    phi  /= sumW;
}

//////////////////////////////////////////////////////////
double TranspPolygonal :: WachspressShapeFGrad(const Point *x, Matrix &phiGrad) const {
    Vector phi(nnodes);
    WachspressShapeF(x, phi);
    Matrix R(nnodes, 2);

    Vector h(nnodes);
    for ( unsigned i = 0; i < nnodes; i++ ) {
        h [ i ] = abs(dot(nodes [ faces [ i ] [ 0 ] ]->givePoint() - *x, normals [ i ]) );
    }
    unsigned oldi = nnodes - 1;
    Vector phiR(2);
    phiR [ 0 ] = 0;
    phiR [ 1 ] = 0;
    for ( unsigned i = 0; i < nnodes; i++ ) {
        R [ i ] [ 0 ] = normals [ oldi ].getX() / h [ oldi ] + normals [ i ].getX() / h [ i ];
        R [ i ] [ 1 ] = normals [ oldi ].getY() / h [ oldi ] + normals [ i ].getY() / h [ i ];
        phiR [ 0 ] += R [ i ] [ 0 ] * phi [ i ];
        phiR [ 1 ] += R [ i ] [ 1 ] * phi [ i ];
        oldi = i;
    }
    for ( unsigned i = 0; i < nnodes; i++ ) {
        for ( unsigned j = 0; j < 2; j++ ) {
            phiGrad [ j ] [ i ] = phi [ i ] * ( R [ i ] [ j ] - phiR [ j ] );
        }
    }
    return 1;
}

//////////////////////////////////////////////////////////
void TranspPolygonal :: setIntegrationPointsAndWeights() {
    if ( ip_type.compare("quad") == 0 ) {
        //based on quadrilateral isoparametric elements
        ip_locs.resize(4 * nnodes);
        ip_weights.resize(4 * nnodes);
        stats.resize(4 * nnodes);
        Point a = ( nodes [ faces [ nnodes - 1 ] [ 0 ] ]->givePoint() + nodes [ faces [ nnodes - 1 ] [ 1 ] ]->givePoint() ) / 2;
        Point b, c, d, derxyr, derxys;
        double detJ;
        d = centroid;
        double r = 1. / sqrt(3.);

        for ( unsigned i = 0; i < nnodes; i++ ) {
            c = ( nodes [ faces [ i ] [ 0 ] ]->givePoint() + nodes [ faces [ i ] [ 1 ] ]->givePoint() ) / 2;
            b = nodes [ i ]->givePoint();
            for ( int k = -1; k < 2; k = k + 2 ) {
                for ( int l = -1; l < 2; l = l + 2 ) {
                    ip_locs [ 4 * i + ( k + 1 ) + ( l + 1 ) / 2 ] = ( a * ( ( 1 + k * r ) * ( 1 + l * r ) ) + b * ( ( 1 - k * r ) * ( 1 + l * r ) ) + c * ( ( 1 - k * r ) * ( 1 - l * r ) ) + d * ( ( 1 + k * r ) * ( 1 - l * r ) ) ) / 4.;
                    derxyr = ( a * ( 1 + l * r ) - b * ( 1 + l * r ) - c * ( 1 - l * r ) + d * ( 1 - l * r ) ) / 4.;
                    derxys = ( a * ( 1 + k * r ) + b * ( 1 - k * r ) - c * ( 1 - k * r ) - d * ( 1 + k * r ) ) / 4.;
                    detJ = derxyr.x * derxys.y - derxyr.y * derxys.x;
                    ip_weights [ 4 * i + ( k + 1 ) + ( l + 1 ) / 2 ] = detJ;
                    stats [ 4 * i + ( k + 1 ) + ( l + 1 ) / 2 ] = mat->giveNewMaterialStatus(this);
                }
            }
            a = c;
        }
    } else if ( ip_type.compare("tri") == 0 ) {
        //based on triangular isoparametric elements
        ip_locs.resize(3 * nnodes);
        ip_weights.resize(3 * nnodes);
        stats.resize(3 * nnodes);
        Point a, b;

        for ( unsigned i = 0; i < nnodes; i++ ) {
            a = nodes [ faces [ i ] [ 0 ] ]->givePoint();
            b = nodes [ faces [ i ] [ 1 ] ]->givePoint();
            double tarea = triArea2D(a, b, centroid);
            ip_locs [ 3 * i  ] = b * ( 1. - 1. / 6. - 1. / 6. ) + a / 6. + centroid / 6.;
            ip_locs [ 3 * i + 1 ] = b * ( 1. - 2. / 3. - 1. / 6. ) + a * 2. / 3. + centroid / 6.;
            ip_locs [ 3 * i + 2 ] = b * ( 1. - 1. / 6. - 2. / 3. ) + a / 6. + centroid * 2. / 3.;
            for ( unsigned t = 0; t < 3; t++ ) {
                ip_weights [ 3 * i + t ] = tarea / 3.;
                stats [ 3 * i + t ] = mat->giveNewMaterialStatus(this);
            }
        }
    } else {
        cerr << "Error in " << name << ": ip_type '" << ip_type << "' not implemented" << endl;
    }
}

//////////////////////////////////////////////////////////
void TranspPolygonal :: init() {

    //check that nodes are TrsNodes
    for ( const auto n: nodes ) {
        TrsNode *p = dynamic_cast< TrsNode * >( n );
        if ( !p ) {
            cerr << "Error in " << name << ": nodes must be inherited from TrsNode, " << n->giveName() << " provided" << endl;
            exit(1);
        }
    }
    //check that material is DisMechMat
    TrsprtMaterial *p = dynamic_cast< TrsprtMaterial * >( mat );
    if ( !p ) {
        cerr << "Error in " << name << ": material must be inherited from TrsprtMaterial, " << mat->giveName() << " provided" << endl;
        exit(1);
    }

    //reorder nodes before calling base calls initialization
    sort2D();

    Element :: init(); //calling base class method;
}

//////////////////////////////////////////////////////////
Matrix TranspPolygonal :: giveBMatrix(const Point *x) const  {
    Matrix B(ndim, nodes.size() );
    shapeFGrad(x, B);
    return B;
}

//////////////////////////////////////////////////////////
Matrix TranspPolygonal :: giveHMatrix(const Point *x) const {
    Vector X(nnodes);
    shapeF(x, X);
    Matrix H(1,nnodes);
    for(unsigned k=0; k<nnodes; k++) H[0][k] = X[k];
    return H;
}


//////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////
// TRANSPORT VIRTUAL POLYGONAL ELEMENT

TranspVirtPolygonal :: TranspVirtPolygonal(const unsigned dim) : TranspPolygonal(dim) {
    name = "TranspVirtPolygonal";
}

//////////////////////////////////////////////////////////
void TranspVirtPolygonal :: init() {
    TranspPolygonal :: init(); //calling base class method;

    double radius = pow(volume/M_PI, 0.5);
    Matrix D(nnodes,ndim+1);
    Point x;
    unsigned i,j,v;
    for( i=0; i<nnodes; i++){
        D[i][0] = 1.;
        x = nodes[i]->givePoint();
        for( v=0; v<ndim; v++){
            D[i][v+1] = (x.giveCoord(v) - centroid.giveCoord(v))/radius;
        }
    }

    Matrix B(ndim+1, nnodes);
     j = nnodes - 1;
    for (  i = 0; i < nnodes; i++ ) {
        B[0][i] = 1./nnodes;
        for( v=0; v<ndim; v++){
            B [ v+1 ] [ i ] = ( normals [ i ].giveCoord(v) * surfaces [ i ] + normals [ j ].giveCoord(v) * surfaces [ j ] );
        }
        j = i;
    }
    B /= 2.*radius;

    Matrix G = matrix_multiply(B,D);
    Matrix Gtilde = G;
    for( i=0; i<ndim+1; i++) Gtilde[0][i] = 0.;

    double Gdet = G[0][0]*G[1][1]*G[2][2] + G[0][2]*G[1][0]*G[2][1] + G[2][0]*G[0][1]*G[1][2] - G[0][0]*G[2][1]*G[1][2] - G[0][1]*G[1][0]*G[2][2] - G[0][2]*G[1][1]*G[2][0];
    Matrix Ginv(ndim+1, ndim+1);
    Ginv[0][0] = G[1][1]*G[2][2]-G[1][2]*G[2][1];
    Ginv[1][0] = -G[1][0]*G[2][2]-G[1][2]*G[2][0];
    Ginv[2][0] = G[1][0]*G[2][1]-G[1][1]*G[2][0];
    Ginv[0][1] = -G[0][1]*G[2][2]-G[0][2]*G[2][1];
    Ginv[1][1] = G[0][0]*G[2][2]-G[0][2]*G[2][0];
    Ginv[2][1] = -G[0][0]*G[2][1]-G[0][1]*G[2][0];
    Ginv[0][2] = G[0][1]*G[1][2]-G[0][2]*G[1][1];
    Ginv[1][2] = -G[0][0]*G[1][2]-G[0][2]*G[1][0];
    Ginv[2][2] = G[0][0]*G[1][1]-G[0][1]*G[1][0];
    Ginv /= Gdet;

    V1 = matrix_multiply(Ginv,B);
    V2 =matrix_multiply(D,V1) * (-1.);
    for(i=0; i<nnodes; i++) V2[i][i] += 1.;
    V1 = matrix_multiply(matrix_multiply(V1.transpose(),Gtilde),V1);

    Matrix H(ndim+1,ndim+1);
    Vector m(ndim+1);
    m[0] = 1;
    for ( size_t i = 0; i < ip_weights.size(); i++ ) {
        for( v=0; v<ndim; v++)  m[v+1] = (ip_locs[i].giveCoord(v)-centroid.giveCoord(v))/radius;
        H += dyadicProduct(m,m*ip_weights [ i ]);
    }

    double Hdet = H[0][0]*H[1][1]*H[2][2] + H[0][2]*H[1][0]*H[2][1] + H[2][0]*H[0][1]*H[1][2] - H[0][0]*H[2][1]*H[1][2] - H[0][1]*H[1][0]*H[2][2] - H[0][2]*H[1][1]*H[2][0];
    Matrix Hinv(ndim+1, ndim+1);
    Hinv[0][0] = H[1][1]*H[2][2]-H[1][2]*H[2][1];
    Hinv[1][0] = -H[1][0]*H[2][2]-H[1][2]*H[2][0];
    Hinv[2][0] = H[1][0]*H[2][1]-H[1][1]*H[2][0];
    Hinv[0][1] = -H[0][1]*H[2][2]-H[0][2]*H[2][1];
    Hinv[1][1] = H[0][0]*H[2][2]-H[0][2]*H[2][0];
    Hinv[2][1] = -H[0][0]*H[2][1]-H[0][1]*H[2][0];
    Hinv[0][2] = H[0][1]*H[1][2]-H[0][2]*H[1][1];
    Hinv[1][2] = -H[0][0]*H[1][2]-H[0][2]*H[1][0];
    Hinv[2][2] = H[0][0]*H[1][1]-H[0][1]*H[1][0];
    Hinv /= Hdet;

    W2 = matrix_multiply(Ginv,B);
    Matrix C = matrix_multiply(H,W2);
    W2 = matrix_multiply(D,W2)*(-1.);
    for(i=0; i<nnodes; i++) W2[i][i] += 1.;
    W1 = matrix_multiply(matrix_multiply(C.transpose(), Hinv), C);
}

//////////////////////////////////////////////////////////
Matrix TranspVirtPolygonal :: giveSteadyStateMatrix(string matrixType) const {
    Matrix C = TranspPolygonal :: giveSteadyStateMatrix(matrixType);
    double cond = 0;
    for ( size_t i = 0; i < ip_weights.size(); i++ ) {
        cond += ip_weights [ i ] * stats [ i ]->giveStiffnessTensor("elastic", ndim)[0][0];
    }
    cond /= volume;

    return V1 * cond + matrix_multiply(matrix_multiply(V2.transpose(), C), V2);
}

//////////////////////////////////////////////////////////
Matrix TranspVirtPolygonal :: giveMassMatrix() const {
    Matrix M = TranspPolygonal :: giveMassMatrix();
    double cap = 0;
    for ( size_t i = 0; i < ip_weights.size(); i++ ) {
        cap += ip_weights [ i ] * stats [ i ]->giveMassConstant();
    }
    cap /= volume;

    /*
    cout << "------------" << endl;
    M.print();
    cout << "------------" << endl;
    ( W1 * cap + matrix_multiply(matrix_multiply(W2.transpose(), M), W2) ).print();
    cout << "------------" << endl;
    ( (W1 + matrix_multiply(W2.transpose(), W2)*volume*volume )* cap ).print();
    exit(1);
    return M;
    */
    return W1 * cap + matrix_multiply(matrix_multiply(W2.transpose(), M), W2);
}

//////////////////////////////////////////////////////////
Vector TranspVirtPolygonal :: giveInternalForces(const Vector &DoFs, bool frozen){
    (void) frozen;
    //return Element::giveInternalForces(DoFs, frozen); //incorrect integration
    return  giveSteadyStateMatrix("elastic")*DoFs;  //using VEM integration, only elastic material!
}

//////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////
// TRANSPORT POLYGONAL ELEMENT CONSTRUCTED BY STATIC CONDENSATION OF ISOPARAMETRIC TRIANGLES

TranspCondensedPolygonal :: TranspCondensedPolygonal(const unsigned dim) : TranspPolygonal(dim) {
    name = "TranspCondensedPolygonal";
    ip_type = "tri";
}

//////////////////////////////////////////////////////////
void TranspCondensedPolygonal :: fullShapeF(const Point *x, Vector &phi) const {
    unsigned face = findFaceNumber(*x);
    double tarea = triArea2D(centroid, nodes [ faces [ face ] [ 0 ] ]->givePoint(), nodes [ faces [ face ] [ 1 ] ]->givePoint() );
    phi [ faces [ face ] [ 1 ] ] = triArea2D(*x, centroid, nodes [ faces [ face ] [ 0 ] ]->givePoint() ) / tarea;
    phi [ faces [ face ] [ 0 ] ] = triArea2D(*x, nodes [ faces [ face ] [ 1 ] ]->givePoint(), centroid) / tarea;
    phi [ nnodes ] = triArea2D(*x, nodes [ faces [ face ] [ 0 ] ]->givePoint(), nodes [ faces [ face ] [ 1 ] ]->givePoint() ) / tarea;
}

//////////////////////////////////////////////////////////
double TranspCondensedPolygonal :: fullShapeFGrad(const Point *x, Matrix &phiGrad) const {
    unsigned face = findFaceNumber(*x);
    phiGrad *= 0.;
    double tarea = triArea2D(centroid, nodes [ faces [ face ] [ 0 ] ]->givePoint(), nodes [ faces [ face ] [ 1 ] ]->givePoint() );
    phiGrad [ 0 ] [ faces [ face ] [ 1 ] ] = 0.5 * ( centroid.getY() - nodes [ faces [ face ] [ 0 ] ]->givePoint().getY() ) / tarea;
    phiGrad [ 1 ] [ faces [ face ] [ 1 ] ] = 0.5 * ( nodes [ faces [ face ] [ 0 ] ]->givePoint().getX() - centroid.getX() ) / tarea;
    phiGrad [ 0 ] [ faces [ face ] [ 0 ] ] = 0.5 * ( nodes [ faces [ face ] [ 1 ] ]->givePoint().getY() - centroid.getY() ) / tarea;
    phiGrad [ 1 ] [ faces [ face ] [ 0 ] ] = 0.5 * ( centroid.getX() - nodes [ faces [ face ] [ 1 ] ]->givePoint().getX() ) / tarea;
    phiGrad [ 0 ] [ nnodes ]               = 0.5 * ( nodes [ faces [ face ] [ 0 ] ]->givePoint().getY() - nodes [ faces [ face ] [ 1 ] ]->givePoint().getY() ) / tarea;
    phiGrad [ 1 ] [ nnodes ]               = 0.5 * ( nodes [ faces [ face ] [ 1 ] ]->givePoint().getX() - nodes [ faces [ face ] [ 0 ] ]->givePoint().getX() ) / tarea;
    return 1.;
}

//////////////////////////////////////////////////////////
unsigned TranspCondensedPolygonal :: findFaceNumber(Point x) const {
    double alpha = atan2(x.getY() - centroid.getY(), x.getX() - centroid.getX() );
    double a, b;
    unsigned face;
    for ( face = 0; face < nnodes; face++ ) {
        a = angles[faces [ face ] [ 0 ]];
        b = angles[faces [ face ] [ 1 ]];
        if (a<b && a<alpha && b>alpha) {
            return face;
        }
        if (a>b && (a<alpha || b>alpha) ){
           return face;
        }
    }
    cerr << "Error in TranspCondensedPolygonal: findFaceNumber should never go here" << endl;
    exit(1);
}

//////////////////////////////////////////////////////////
void TranspCondensedPolygonal :: shapeF(const Point *x, Vector &phi) const {
    Vector full(nnodes+1);
    fullShapeF(x, full);
    for ( unsigned i = 0; i < nnodes; i++ ) {
        phi [ i ] = full [ i ] + full [ nnodes ] * red2full [ i ];
    }
}

//////////////////////////////////////////////////////////
double TranspCondensedPolygonal :: shapeFGrad(const Point *x, Matrix &phiGrad) const {
    Matrix full(ndim, nnodes+1);
    fullShapeFGrad(x, full);
    for ( unsigned d = 0; d < ndim; d++ ) {
        for ( unsigned i = 0; i < nnodes; i++ ) {
            phiGrad [ d ] [ i ] = full [ d ] [ i ] + full [ d ] [ nnodes ] * red2full [ i ];
        }
    }
    return 1.;
}

//////////////////////////////////////////////////////////
void TranspCondensedPolygonal :: init() {
    red2full.resize(nnodes);

    sort2D();
    angles.resize(nnodes);
    for ( unsigned i = 0; i < nnodes; i++ ) {
        angles [ i ] = atan2(nodes [ i ]->givePoint().getY() - centroid.getY(), nodes [ i ]->givePoint().getX() - centroid.getX() );
    }

    TranspPolygonal :: init(); //calling base class method;


    //build transformation matrix allowing to calculate inner degree of freedom
    Matrix FullK(nnodes + 1, nnodes + 1);
    Matrix phiGrad(ndim, nnodes+1);
    for ( size_t i = 0; i < ip_weights.size(); i++ ) {
        fullShapeFGrad(& ip_locs [ i ], phiGrad);
        FullK += matrix_multiply(phiGrad.transpose(), phiGrad) * ip_weights [ i ];
    }

    for ( unsigned i = 0; i < nnodes; i++ ) {
        red2full [ i ] = -FullK [ i ] [ nnodes ] / FullK [ nnodes ] [ nnodes ];
    }

    //update of B matrices requested, because the previous calculation didn't have correct shape functions
    for ( unsigned k = 0; k < ip_locs.size(); k++ ) {
        Bs [ k ] = giveBMatrix(& ip_locs [ k ]);
    }
}

/*
 * //////////////////////////////////////////////////////////
 * // 3 D
 * //////////////////////////////////////////////////////////
 * //////////////////////////////////////////////////////////
 * // POLYHEDRAL FACE
 * PolyhedralFace :: PolyhedralFace(const unsigned dim) {
 *  ndim = 3;
 *  if (dim!=3){
 *      cerr << "Polyhedral Face error: this element can be used only in 3D models" << endl;
 *      exit(1);
 *  }
 *  name = "PolyhedralFace";
 * }
 *
 * //////////////////////////////////////////////////////////
 * void PolyhedralFace :: readFromLine(istringstream &iss, NodeContainer *fullnodes, MaterialContainer *fullmatrs){
 *  (void) fullmatrs;
 *
 *  unsigned num, num2;
 *
 *  iss >> num;
 *  nodes.resize(num);
 *  for ( unsigned i = 0; i < num; i++ ) {
 *      iss >> num2;
 *      nodes [ i ] = fullnodes->giveNode(num2);
 *  }
 * }
 *
 * //////////////////////////////////////////////////////////
 * void PolyhedralFace :: init(){
 * }
 */

/*
//////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////
// TRANSPORT POLYHEDRAL ELEMENT
TranspPolyhedral :: TranspPolyhedral(const unsigned dim)  : TranspPolygonal(2) {
    ndim = dim;
    name = "TranspPolyhedral";
    ip_type = "quad";
}

//////////////////////////////////////////////////////////
void TranspPolyhedral :: readFromLine(istringstream &iss, NodeContainer *fullnodes, MaterialContainer *fullmatrs) {
    unsigned num, num2;
    unordered_set< unsigned >nn;

    iss >> nfaces;
    faces.resize(nfaces);
    for ( unsigned i = 0; i < nfaces; i++ ) {
        iss >> num;
        faces [ i ].resize(num);
        for ( unsigned j = 0; j < num; j++ ) {
            iss >> num2;
            nn.insert(num2);
            faces [ i ] [ j ] = num2;
        }
    }
    iss >> num;
    mat = fullmatrs->giveMaterial(num);

    nnodes = nn.size();
    vector< unsigned >nv;
    nv.resize(nnodes);
    nodes.resize(nnodes);
    faceConnectivity.resize(nnodes);

    unsigned i = 0;
    for ( auto &n: nn ) {
        nv [ i ] = n;
        nodes [ i ] = fullnodes->giveNode(n);
        i++;
    }
    for ( unsigned f = 0; f < nfaces; f++ ) {
        for ( unsigned n = 0; n < faces [ f ].size(); n++ ) {
            i = 0;
            while ( nv [ i ] != faces [ f ] [ n ] && i < nv.size() ) {
                i++;
            }
            if ( i == nv.size() ) {
                cerr << "Inconzistence in Polyhedral input file" << endl;
                exit(1);
            }
            faces [ f ] [ n ] = i;
            faceConnectivity [ i ].push_back(f);
        }
    }
}


//////////////////////////////////////////////////////////
void TranspPolyhedral :: WachspressShapeF(const Point x, Matrix phi) const {
    //compute hf
    Vector h(nfaces);
    for ( unsigned i = 0; i < nfaces; i++ ) {
        h [ i ] = dot(nodes [ faces [ i ] [ 0 ] ]->givePoint() - x, normals [ i ]);
    }

    //compute wf
    double sumW = 0;
    phi *= 0;
    for ( unsigned i = 0; i < nnodes; i++ ) {
        phi[ i ] = 1;
        for ( unsigned j = 0; j < 3; j++ ) {
            w [0] [ i ] *= h [ faceConnectivity [ i ] [ j ] ];
        }
        phi[ i ] = determinants [ i ] / w [ i ];
        sumW += phi[ i ];
    }
    return phi / sumW;
}

//////////////////////////////////////////////////////////
double TranspPolyhedral :: WachspressShapeFGrad(const Point x, Matrix phiGrad) const {
    Vector phi(nnodes);
    WachspressShapeF(x, phi);
    Matrix R(nnodes, 3);

    Vector h(nfaces);
    for ( unsigned i = 0; i < nfaces; i++ ) {
        h [ i ] = dot(nodes [ faces [ i ] [ 0 ] ]->givePoint() - x, normals [ i ]);
    }

    unsigned n;
    Vector phiR(3);
    phiR [ 0 ] = 0;
    phiR [ 1 ] = 0;
    phiR [ 2 ] = 0;
    for ( unsigned i = 0; i < nnodes; i++ ) {
        R [ i ] [ 0 ] = 0;
        R [ i ] [ 1 ] = 0;
        R [ i ] [ 2 ] = 0;
        for ( unsigned j = 0; j < 3; j++ ) {
            n = faceConnectivity [ i ] [ 0 ];
            R [ i ] [ 0 ] += normals [ n ].getX() / h [ n ];
            R [ i ] [ 1 ] += normals [ n ].getY() / h [ n ];
            R [ i ] [ 2 ] += normals [ n ].getZ() / h [ n ];
        }
        phiR [ 0 ] += R [ i ] [ 0 ] * phi[ i ];
        phiR [ 1 ] += R [ i ] [ 1 ] * phi[ i ];
        phiR [ 2 ] += R [ i ] [ 2 ] * phi[ i ];
    }

    for ( unsigned i = 0; i < nnodes; i++ ) {
        for ( unsigned j = 0; j < 3; j++ ) {
            phiGrad [ j ] [ i ] = phi[ i ] * ( R [ i ] [ j ] - phiR [ j ] );
        }
    }

    return 1;
}

//////////////////////////////////////////////////////////
void TranspPolyhedral :: findIntegrationPoints() {
    if ( ip_type.compare("quad") == 0 ) {
        //based on quadrilateral isoparametric elements
        ip_locs.resize(4 * nnodes);
        ip_weights.resize(4 * nnodes);
        stats.resize(4 * nnodes);
        Point a = ( nodes [ faces [ nnodes - 1 ] [ 0 ] ]->givePoint() + nodes [ faces [ nnodes - 1 ] [ 1 ] ]->givePoint() ) / 2;
        Point b, c, d, derxyr, derxys;
        double detJ;
        d = centroid;
        double r = 1. / sqrt(3.);
        *
        for ( unsigned i = 0; i < nnodes; i++ ) {
            c = ( nodes [ faces [ i ] [ 0 ] ]->givePoint() + nodes [ faces [ i ] [ 1 ] ]->givePoint() ) / 2;
            b = nodes [ i ]->givePoint();
            for ( int k = -1; k < 2; k = k + 2 ) {
                for ( int l = -1; l < 2; l = l + 2 ) {
                    ip_locs [ 4 * i + ( k + 1 ) + ( l + 1 ) / 2 ] = ( a * ( ( 1 + k * r ) * ( 1 + l * r ) ) + b * ( ( 1 - k * r ) * ( 1 + l * r ) ) + c * ( ( 1 - k * r ) * ( 1 - l * r ) ) + d * ( ( 1 + k * r ) * ( 1 - l * r ) ) ) / 4.;
                    derxyr = ( a * ( 1 + l * r ) - b * ( 1 + l * r ) - c * ( 1 - l * r ) + d * ( 1 - l * r ) ) / 4.;
                    derxys = ( a * ( 1 + k * r ) + b * ( 1 - k * r ) - c * ( 1 - k * r ) - d * ( 1 + k * r ) ) / 4.;
                    detJ = derxyr.x * derxys.y - derxyr.y * derxys.x;
                    ip_weights [ 4 * i + ( k + 1 ) + ( l + 1 ) / 2 ] = detJ;
                    stats [ 4 * i + ( k + 1 ) + ( l + 1 ) / 2 ] = mat->giveNewMaterialStatus(this);
                }
            }
            a = c;
        }
    } else if ( ip_type.compare("tri") == 0 ) {
        //based on triangular isoparametric elements
        ip_locs.resize(3 * nnodes);
        ip_weights.resize(3 * nnodes);
        stats.resize(3 * nnodes);
        Point a, b;

        for ( unsigned i = 0; i < nnodes; i++ ) {
            a = nodes [ faces [ i ] [ 0 ] ]->givePoint();
            b = nodes [ faces [ i ] [ 1 ] ]->givePoint();
            double tarea = triArea2D(a, b, centroid);
            ip_locs [ 3 * i  ] = b * ( 1. - 1. / 6. - 1. / 6. ) + a / 6. + centroid / 6.;
            ip_locs [ 3 * i + 1 ] = b * ( 1. - 2. / 3. - 1. / 6. ) + a * 2. / 3. + centroid / 6.;
            ip_locs [ 3 * i + 2 ] = b * ( 1. - 1. / 6. - 2. / 3. ) + a / 6. + centroid * 2. / 3.;
            for ( unsigned t = 0; t < 3; t++ ) {
                ip_weights [ 3 * i + t ] = tarea / 3.;
                stats [ 3 * i + t ] = mat->giveNewMaterialStatus(this);
            }
        }
       } else {
        cerr << "Error in " << name << ": ip_type '" << ip_type << "' not implemented" << endl;
       }
}

//////////////////////////////////////////////////////////
void TranspPolyhedral :: init() {
    Element :: init(); //calling base class method;

    //check that nodes are TrsNodes
    for ( const auto n: nodes ) {
        TrsNode *p = dynamic_cast< TrsNode * >( n );
        if ( !p ) {
            cerr << "Error in " << name << ": nodes must be inherited from TrsNode, " << n->giveName() << " provided" << endl;
            exit(1);
        }
    }
    //check that material is TrsprtMaterial
    TrsprtMaterial *p = dynamic_cast< TrsprtMaterial * >( mat );
    if ( !p ) {
        cerr << "Error in " << name << ": material must be inherited from TrsprtMaterial, " << mat->giveName() << " provided" << endl;
        exit(1);
    }

    //centroid estimation
    Point estcentroid = Point(0, 0, 0);
    for ( auto &n: nodes ) {
        estcentroid += n->givePoint();
    }
    estcentroid /= nnodes;

    //centers, areas and volumes of faces
    surfaces.resize(nfaces);
    faceCenters.resize(nfaces);
    volumes.resize(nfaces);

    Point *last;
    Point *current;
    Point estcenter;
    unsigned i = 0;
    double fa, fv;
    volume = 0.;
    centroid = Point(0, 0, 0);
    for ( auto &f: faces ) {
        //estimate face center
        estcenter = Point(0, 0, 0);
        for ( auto &n: f ) {
            estcenter += nodes [ n ]->givePoint();
        }
        estcenter /= f.size();


        last = nodes [ f [ f.size() - 1 ] ]->givePointPointer();
        surfaces [ i ] = 0;
        volumes [ i ] = 0;
        faceCenters [ i ] = Point(0, 0, 0);
        for ( auto &n: f ) {
            current = nodes [ n ]->givePointPointer(); \
            fa = triArea3D(last, current, & estcenter);
            fv = -tetVolume3D(last, current, & estcenter, & estcentroid); //normal oriented outside
            surfaces [ i ] += fa;
            volumes [ i ] += fv;
            faceCenters [ i ] = ( * last + * current ) * fa;
            centroid = ( * last + * current + estcenter ) * fv;
            last = current;
        }
        faceCenters [ i ] /= 3. * surfaces [ i ];
        faceCenters [ i ] += estcenter / 3.;
        volume += volumes [ i ];
        i++;
    }
    centroid /= 4. * volume;
    centroid += estcentroid / 4;

    i = 0;
    //correct volumes of each face
    for ( auto &f: faces ) {
        last = nodes [ f [ f.size() - 1 ] ]->givePointPointer();
        volumes [ i ] = 0;
        for ( auto &n: f ) {
            current = nodes [ n ]->givePointPointer(); \
            volumes [ i ] = tetVolume3D(last, current, & faceCenters [ i ], & centroid);
            last = current;
        }
        i++;
    }

    //normals
    i = 0;
    for ( auto &f: faces ) {
        Point a = nodes [ f [ 1 ] ]->givePoint() - faceCenters [ i ];
        Point b = nodes [ f [ 0 ] ]->givePoint() - faceCenters [ i ];
        normals [ i ].setX(a.getY() * b.getZ() - a.getZ() * b.getY() );
        normals [ i ].setX(a.getZ() * b.getX() - a.getX() * b.getZ() );
        normals [ i ].setX(a.getX() * b.getY() - a.getY() * b.getX() );
        normals [ i ] /= normals [ i ].norm();
        i++;
    }

    //check face connectivity and compute determinants from normals
    i = 0;
    determinants.resize(nnodes);
    Point a, b, c;
    for ( auto &fc: faceConnectivity ) {
        if ( fc.size() != 3 ) {
            cerr << "Error in TranspPolyhedral: each node must have 3 connected faces, " << fc.size() << " found" << endl;
            exit(EXIT_FAILURE);
        }
        a = nodes [ fc [ 0 ] ]->givePoint();
        b = nodes [ fc [ 1 ] ]->givePoint();
        c = nodes [ fc [ 2 ] ]->givePoint();
        determinants [ i ] = a.getX() * ( b.getY() * c.getZ() - b.getZ() * c.getY() ) - b.getX() * ( a.getY() * c.getZ() - a.getZ() * c.getY() ) + c.getX() * ( a.getY() * b.getZ() - a.getZ() * b.getY() );
        if ( determinants [ i ] < 0 ) {
            determinants [ i ] *= -1;
            unsigned p = fc [ 1 ];
            fc [ 1 ] = fc [ 2 ];
            fc [ 2 ] = p;
        }
        i++;
    }

    findIntegrationPoints();
}
*/

/*
 * //////////////////////////////////////////////////////////
 * //////////////////////////////////////////////////////////
 * // TRANSPORT VIRTUAL POLYHEDRAL ELEMENT
 *
 * TranspVirtPolyhedral :: TranspVirtPolyhedral(const unsigned dim) : TranspPolyhedral(dim) {
 *  name = "TranspVirtPolyhedral";
 * }
 *
 * //////////////////////////////////////////////////////////
 * void TranspVirtPolyhedral :: init() {
 *  TranspPolyhedral :: init(); //calling base class method;
 *
 *  Matrix R(nnodes, ndim);
 *  unsigned j = nnodes - 1;
 *  for ( unsigned i = 0; i < nnodes; i++ ) {
 *      R [ i ] [ 0 ] = ( normals [ i ].x * surfaces [ i ] + normals [ j ].x * surfaces [ j ] ) / 2.;
 *      R [ i ] [ 1 ] = ( normals [ i ].y * surfaces [ i ] + normals [ j ].y * surfaces [ j ] ) / 2.;
 *      j = i;
 *  }
 *
 *  Matrix N(nnodes, ndim);
 *  Point x;
 *  for ( unsigned i = 0; i < nnodes; i++ ) {
 *      x = nodes [ i ]->givePoint();
 *      N [ i ] [ 0 ] = x.x;
 *      N [ i ] [ 1 ] = x.y;
 *  }
 *
 *  Matrix I(nnodes, nnodes);
 *  for ( unsigned i = 0; i < nnodes; i++ ) {
 *      I [ i ] [ i ] = 1.;
 *  }
 *
 *  Matrix P0(nnodes, nnodes);
 *  for ( unsigned i = 0; i < nnodes; i++ ) {
 *      for ( unsigned j = 0; j < nnodes; j++ ) {
 *          P0 [ i ] [ j ] = 1. / nnodes;
 *      }
 *  }
 *
 *  Matrix H = matrix_multiply(N, R.transpose() ) / volume;
 *  Matrix P =  H +  matrix_multiply(P0, I - H);
 *
 *  V1 = matrix_multiply(R, R.transpose() ) / volume;
 *  V2 = I - P;
 * }
 *
 * //////////////////////////////////////////////////////////
 * Matrix TranspVirtPolyhedral :: giveConductivityMatrix(string matrixType) const {
 *  Matrix C = TranspPolyhedral :: giveConductivityMatrix(matrixType);
 *  TrsprtMaterial *tmat = static_cast< TrsprtMaterial * >( mat );
 *  return V1 * tmat->giveConductivity() + matrix_multiply(matrix_multiply(V2.transpose(), C), V2);
 * }
 *
 *
 * //////////////////////////////////////////////////////////
 * //////////////////////////////////////////////////////////
 * // TRANSPORT POLYHEDRAL ELEMENT CONSTRUCTED BY STATIC CONDENSATION OF ISOPARAMETRIC TRIANGLES
 *
 * TranspCondensedPolyhedral :: TranspCondensedPolyhedral(const unsigned dim) : TranspPolyhedral(dim) {
 *  name = "TranspCondensedPolyhedral";
 *  ip_type = "tri";
 * }
 *
 * //////////////////////////////////////////////////////////
 * Vector TranspCondensedPolyhedral :: fullTriShapeF(Point x) const {
 *  unsigned face = findFaceNumber(x);
 *  Vector phi(0., nnodes + 1); //include the centroid
 *  double tarea = triArea2D(centroid, nodes [ faces [ face ] [ 0 ] ]->givePoint(), nodes [ faces [ face ] [ 1 ] ]->givePoint() );
 *  phi [ faces [ face ] [ 1 ] ] = triArea2D(x, centroid, nodes [ faces [ face ] [ 0 ] ]->givePoint() ) / tarea;
 *  phi [ faces [ face ] [ 0 ] ] = triArea2D(x, nodes [ faces [ face ] [ 1 ] ]->givePoint(), centroid) / tarea;
 *  phi [ nnodes ] = triArea2D(x, nodes [ faces [ face ] [ 0 ] ]->givePoint(), nodes [ faces [ face ] [ 1 ] ]->givePoint() ) / tarea;
 *  return phi;
 * }
 *
 * //////////////////////////////////////////////////////////
 * Matrix TranspCondensedPolyhedral :: fullTriShapeFGrad(Point x) const {
 *  unsigned face = findFaceNumber(x);
 *  Matrix phiGrad(ndim, nnodes + 1); //include the centroid
 *  double tarea = triArea2D(centroid, nodes [ faces [ face ] [ 0 ] ]->givePoint(), nodes [ faces [ face ] [ 1 ] ]->givePoint() );
 *  phiGrad [ 0 ] [ faces [ face ] [ 1 ] ] = 0.5 * ( centroid.getY() - nodes [ faces [ face ] [ 0 ] ]->givePoint().getY() ) / tarea;
 *  phiGrad [ 1 ] [ faces [ face ] [ 1 ] ] = 0.5 * ( nodes [ faces [ face ] [ 0 ] ]->givePoint().getX() - centroid.getX() ) / tarea;
 *  phiGrad [ 0 ] [ faces [ face ] [ 0 ] ] = 0.5 * ( nodes [ faces [ face ] [ 1 ] ]->givePoint().getY() - centroid.getY() ) / tarea;
 *  phiGrad [ 1 ] [ faces [ face ] [ 0 ] ] = 0.5 * ( centroid.getX() - nodes [ faces [ face ] [ 1 ] ]->givePoint().getX() ) / tarea;
 *  phiGrad [ 0 ] [ nnodes ]         = 0.5 * ( nodes [ faces [ face ] [ 0 ] ]->givePoint().getY() - nodes [ faces [ face ] [ 1 ] ]->givePoint().getY() ) / tarea;
 *  phiGrad [ 1 ] [ nnodes ]         = 0.5 * ( nodes [ faces [ face ] [ 1 ] ]->givePoint().getX() - nodes [ faces [ face ] [ 0 ] ]->givePoint().getX() ) / tarea;
 *  return phiGrad;
 * }
 *
 * //////////////////////////////////////////////////////////
 * unsigned TranspCondensedPolyhedral :: findFaceNumber(Point x) const {
 *  double alpha = atan2(x.getY() - centroid.getY(), x.getX() - centroid.getX() );
 *  unsigned face;
 *  for ( face = 0; face < nnodes; face++ ) {
 *      if ( alpha < angles [ faces [ face ] [ 1 ] ] &&  alpha > angles [ faces [ face ] [ 0 ] ] ) {
 *          return face;
 *      }
 *  }
 *  return nodeMaxAngle;
 * }
 *
 * //////////////////////////////////////////////////////////
 * Vector TranspCondensedPolyhedral :: condTriShapeF(Point x) const {
 *  Vector full = fullTriShapeF(x);
 *  Vector reduced(nnodes);
 *  for ( unsigned i = 0; i < nnodes; i++ ) {
 *      reduced [ i ] = full [ i ] + full [ nnodes ] * red2full [ i ];
 *  }
 *  return reduced;
 * }
 *
 * //////////////////////////////////////////////////////////
 * Matrix TranspCondensedPolyhedral :: condTriShapeFGrad(Point x) const {
 *  Matrix full = fullTriShapeFGrad(x);
 *  Matrix reduced(ndim, nnodes);
 *  for ( unsigned d = 0; d < ndim; d++ ) {
 *      for ( unsigned i = 0; i < nnodes; i++ ) {
 *          reduced [ d ] [ i ] = full [ d ] [ i ] + full [ d ] [ nnodes ] * red2full [ i ];
 *      }
 *  }
 *  return reduced;
 * }
 *
 * //////////////////////////////////////////////////////////
 * void TranspCondensedPolyhedral :: init() {
 *  TranspPolyhedral :: init(); //calling base class method;
 *  angles.resize(nnodes);
 *
 *  nodeMaxAngle = 0;
 *  double maxAngle = -2;
 *  for ( unsigned i = 0; i < nnodes; i++ ) {
 *      angles [ i ] = atan2(nodes [ i ]->givePoint().getY() - centroid.getY(), nodes [ i ]->givePoint().getX() - centroid.getX() );
 *      if (maxAngle<angles [ i ]){
 *          maxAngle = angles [ i ];
 *          nodeMaxAngle = i;
 *      }
 *  }
 *
 *  //build transformation matrix allowing to calculate inner degree of freedom
 *  Matrix FullK(nnodes + 1, nnodes + 1);
 *  Matrix phiGrad;
 *  for ( size_t i = 0; i < ip_weights.size(); i++ ) {
 *      phiGrad = fullTriShapeFGrad(ip_locs [ i ]);
 *      FullK += matrix_multiply(phiGrad.transpose(), phiGrad) * ip_weights [ i ];
 *  }
 *
 *  red2full.resize(nnodes);
 *  for ( unsigned i = 0; i < nnodes; i++ ) {
 *      red2full [ i ] = -FullK [ i ] [ nnodes ] / FullK [ nnodes ] [ nnodes ];
 *  }
 * }
 */
