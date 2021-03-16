#include "element_continuous.h"
#include "element_container.h"
#include "boundary_condition.h"

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
    name = "MechanicalQuad";
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
/*
//////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////
// 3D BRICK COSSERAT MECHANICAL ELEMENT
CosseratBrick :: CosseratBrick() {
    ndim = 3;
    name = "CosseratBrick";
}


//////////////////////////////////////////////////////////
Matrix CosseratBrick :: giveBMatrix(const Point *x) const {
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
Matrix CosseratBrick :: giveHMatrix(const Point *x) const {
    Vector phi(nodes.size() );
    shapeF(x, phi);
    return Matrix(0,0);
}
*/
