#include "shape_functions.h"

//////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////
// SHAPE FUNCTIONS - MASTER CLASS
//////////////////////////////////////////////////////////
void ShapeFunc :: init( vector< Node * > &nodes ){
    points.resize(nodes.size());
    unsigned k=0;
    for (auto &n: nodes){
        points[k] = n->givePointPointer();
        k++;
    }
}

//////////////////////////////////////////////////////////
void ShapeFunc :: init( vector< Point * > pp ){
    points = pp;
}

//////////////////////////////////////////////////////////
void ShapeFunc :: giveShapeFGrad(const Point *x, Matrix &phiGrad) const {
    giveShapeFGradNatural(x, phiGrad);
    Matrix JacobiMInverse(ndim, ndim);
    giveJacobiMInverse(x, JacobiMInverse);
    phiGrad = JacobiMInverse * phiGrad;
}

//////////////////////////////////////////////////////////
void ShapeFunc :: giveShapeFGrad(const Point *x, const Matrix &JacobiMInverse, Matrix &phiGrad) const {
    giveShapeFGradNatural(x, phiGrad);
    phiGrad = JacobiMInverse * phiGrad;
}

//////////////////////////////////////////////////////////
void ShapeFunc :: giveJacobiM(const Matrix &phiGradNat, Matrix &JacobiM) const {
    JacobiM *= 0;
    Point *n;
    for ( unsigned i = 0; i < points.size(); i++ ) {
        n = points[i];
        for ( unsigned v1 = 0; v1 < ndim; v1++ ) {
            for ( unsigned v2 = 0; v2 < ndim; v2++ ) {
                JacobiM [ v1 ] [ v2 ] += phiGradNat [ v1 ] [ i ] * n->giveCoord(v2);
            }
        }
    }
}

//////////////////////////////////////////////////////////
void ShapeFunc :: giveJacobiM(const Point *x, Matrix &JacobiM) const {
    Matrix phiGradNat( ndim, points.size() );
    giveShapeFGradNatural(x, phiGradNat);
    JacobiM *= 0;
    Point *n;
    for ( unsigned i = 0; i < points.size(); i++ ) {
        n = points[i];
        for ( unsigned v1 = 0; v1 < ndim; v1++ ) {
            for ( unsigned v2 = 0; v2 < ndim; v2++ ) {
                JacobiM [ v1 ] [ v2 ] += phiGradNat [ v1 ] [ i ] * n->giveCoord(v2);
            }
        }
    }
}

//////////////////////////////////////////////////////////
double ShapeFunc :: giveJacobian(const Point *x) const {
    Matrix JacobiM(ndim, ndim);
    giveJacobiM(x, JacobiM);
    return giveJacobian(JacobiM);
}

//////////////////////////////////////////////////////////
double ShapeFunc :: giveJacobian(const Matrix &JacobiM) const {
    if ( ndim == 1 ) {
        return JacobiM [ 0 ] [ 0 ];
    } else if ( ndim == 2 ) {
        return JacobiM [ 0 ] [ 0 ] * JacobiM [ 1 ] [ 1 ] - JacobiM [ 0 ] [ 1 ] * JacobiM [ 1 ] [ 0 ];
    } else if ( ndim == 3 ) {
        return JacobiM [ 0 ] [ 0 ] * JacobiM [ 1 ] [ 1 ] * JacobiM [ 2 ] [ 2 ] - JacobiM [ 0 ] [ 0 ] * JacobiM [ 1 ] [ 2 ] * JacobiM [ 2 ] [ 1 ] + JacobiM [ 0 ] [ 1 ] * JacobiM [ 1 ] [ 2 ] * JacobiM [ 2 ] [ 0 ] - JacobiM [ 0 ] [ 1 ] * JacobiM [ 1 ] [ 0 ] * JacobiM [ 2 ] [ 2 ] + JacobiM [ 0 ] [ 2 ] * JacobiM [ 1 ] [ 0 ] * JacobiM [ 2 ] [ 1 ] - JacobiM [ 0 ] [ 2 ] * JacobiM [ 1 ] [ 1 ] * JacobiM [ 2 ] [ 0 ];
    } else {
        return 0;
    }
}

//////////////////////////////////////////////////////////
void ShapeFunc :: giveJacobiMInverse(const Matrix &JacobiM, Matrix &JacobiMInverse) const {
    double JacDet = giveJacobian(JacobiM);
    if ( ndim == 1 ) {
        JacobiMInverse [ 0 ] [ 0 ] = 1. / JacobiM [ 0 ] [ 0 ];
    } else if ( ndim == 2 ) {
        JacobiMInverse [ 0 ] [ 0 ] = JacobiM [ 1 ] [ 1 ] / JacDet;
        JacobiMInverse [ 0 ] [ 1 ] = -JacobiM [ 0 ] [ 1 ] / JacDet;
        JacobiMInverse [ 1 ] [ 0 ] = -JacobiM [ 1 ] [ 0 ] / JacDet;
        JacobiMInverse [ 1 ] [ 1 ] = JacobiM [ 0 ] [ 0 ] / JacDet;
    } else if ( ndim == 3 ) {
        JacobiMInverse [ 0 ] [ 0 ] =  ( JacobiM [ 1 ] [ 1 ] * JacobiM [ 2 ] [ 2 ] - JacobiM [ 2 ] [ 1 ] * JacobiM [ 1 ] [ 2 ] ) / JacDet;
        JacobiMInverse [ 0 ] [ 1 ] = -( JacobiM [ 0 ] [ 1 ] * JacobiM [ 2 ] [ 2 ] - JacobiM [ 2 ] [ 1 ] * JacobiM [ 0 ] [ 2 ] ) / JacDet;
        JacobiMInverse [ 0 ] [ 2 ] =  ( JacobiM [ 0 ] [ 1 ] * JacobiM [ 1 ] [ 2 ] - JacobiM [ 1 ] [ 1 ] * JacobiM [ 0 ] [ 2 ] ) / JacDet;
        JacobiMInverse [ 1 ] [ 0 ] = -( JacobiM [ 1 ] [ 0 ] * JacobiM [ 2 ] [ 2 ] - JacobiM [ 2 ] [ 0 ] * JacobiM [ 1 ] [ 2 ] ) / JacDet;
        JacobiMInverse [ 1 ] [ 1 ] =  ( JacobiM [ 0 ] [ 0 ] * JacobiM [ 2 ] [ 2 ] - JacobiM [ 2 ] [ 0 ] * JacobiM [ 0 ] [ 2 ] ) / JacDet;
        JacobiMInverse [ 1 ] [ 2 ] = -( JacobiM [ 0 ] [ 0 ] * JacobiM [ 1 ] [ 2 ] - JacobiM [ 1 ] [ 0 ] * JacobiM [ 0 ] [ 2 ] ) / JacDet;
        JacobiMInverse [ 2 ] [ 0 ] =  ( JacobiM [ 1 ] [ 0 ] * JacobiM [ 2 ] [ 1 ] - JacobiM [ 2 ] [ 0 ] * JacobiM [ 1 ] [ 1 ] ) / JacDet;
        JacobiMInverse [ 2 ] [ 1 ] = -( JacobiM [ 0 ] [ 0 ] * JacobiM [ 2 ] [ 1 ] - JacobiM [ 2 ] [ 0 ] * JacobiM [ 0 ] [ 1 ] ) / JacDet;
        JacobiMInverse [ 2 ] [ 2 ] =  ( JacobiM [ 0 ] [ 0 ] * JacobiM [ 1 ] [ 1 ] - JacobiM [ 1 ] [ 0 ] * JacobiM [ 0 ] [ 1 ] ) / JacDet;
    }
}

//////////////////////////////////////////////////////////
void ShapeFunc :: giveJacobiMInverse(const Point *x, Matrix &JacobiMInverse) const {
    Matrix JacobiM(ndim, ndim);
    giveJacobiM(x, JacobiM);
    giveJacobiMInverse(JacobiM, JacobiMInverse);
}

//////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////
//  1D LINEAR SHAPE FUNCTIONS
//////////////////////////////////////////////////////////
void Linear1DLineShapeF :: giveShapeF(const Point *x, Vector &phi) const {
    //x in natural coordinates
    phi [ 0 ] = 0.5 * ( 1. - x->getX() );
    phi [ 1 ] = 0.5 * ( 1. + x->getX() );
}

//////////////////////////////////////////////////////////
void Linear1DLineShapeF :: giveShapeFGradNatural(const Point *x, Matrix &phiGrad) const {
    ( void ) x;
    //x in natural coordinates
    phiGrad [ 0 ] [ 0 ] = -0.5;
    phiGrad [ 0 ] [ 1 ] = 0.5;
}


//////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////
//  2D LINEAR SHAPE FUNCTIONS FOR QUAD
//////////////////////////////////////////////////////////
void Linear2DQuadShapeF :: giveShapeF(const Point *x, Vector &phi) const {
    //x in natural coordinates
    phi [ 0 ] = 0.25 * ( 1. - x->getX() ) * ( 1. - x->getY() );
    phi [ 1 ] = 0.25 * ( 1. + x->getX() ) * ( 1. - x->getY() );
    phi [ 2 ] = 0.25 * ( 1. + x->getX() ) * ( 1. + x->getY() );
    phi [ 3 ] = 0.25 * ( 1. - x->getX() ) * ( 1. + x->getY() );
}

//////////////////////////////////////////////////////////
void Linear2DQuadShapeF :: giveShapeFGradNatural(const Point *x, Matrix &phiGrad) const {
    //x in natural coordinates
    phiGrad [ 0 ] [ 0 ] = -0.25 * ( 1. - x->getY() );
    phiGrad [ 0 ] [ 1 ] = -phiGrad [ 0 ] [ 0 ];
    phiGrad [ 0 ] [ 2 ] = 0.25 * ( 1. + x->getY() );
    phiGrad [ 0 ] [ 3 ] = -phiGrad [ 0 ] [ 2 ];
    phiGrad [ 1 ] [ 0 ] = -0.25 * ( 1. - x->getX() );
    phiGrad [ 1 ] [ 1 ] = -0.25 * ( 1. + x->getX() );
    phiGrad [ 1 ] [ 2 ] = -phiGrad [ 1 ] [ 1 ];
    phiGrad [ 1 ] [ 3 ] = -phiGrad [ 1 ] [ 0 ];
}

//////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////
//  3D LINEAR SHAPE FUNCTIONS FOR BRICK
//////////////////////////////////////////////////////////
void Linear3DBrickShapeF :: giveShapeF(const Point *x, Vector &phi) const {
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
void Linear3DBrickShapeF :: giveShapeFGradNatural(const Point *x, Matrix &phiGrad) const {
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
}

//////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////
// 2D WACHSPRESS 
//////////////////////////////////////////////////////////
void Wachspress2DShapeF :: giveShapeF(const Point *x, Vector &phi) const {
    unsigned numOfNodes = points.size();
    Vector h(numOfNodes);
    Point oldNormal = normals [ numOfNodes - 1 ];
    for ( unsigned i = 0; i < numOfNodes; i++ ) {
        h [ i ] = abs( dot(* points [ faces [ i ] [ 0 ] ] - * x, normals [ i ]) );
        phi [ i ] = abs(oldNormal.x * normals [ i ].y - oldNormal.y * normals [ i ].x);
        oldNormal = normals [ i ];
    }
    double oldh = h [ numOfNodes - 1 ];
    double sumW = 0;
    for ( unsigned i = 0; i < numOfNodes; i++ ) {
        phi [ i ] /= oldh * h [ i ];
        sumW += phi [ i ];
        oldh = h [ i ];
    }
    phi  /= sumW;
};

//////////////////////////////////////////////////////////
void Wachspress2DShapeF :: giveShapeFGrad(const Point *x, Matrix &phiGrad) const {
    unsigned numOfNodes = points.size();
    Vector phi(numOfNodes);
    giveShapeF(x, phi);
    Matrix R(numOfNodes, 2);

    Vector h(numOfNodes);
    for ( unsigned i = 0; i < numOfNodes; i++ ) {
        h [ i ] = abs( dot(* points [ faces [ i ] [ 0 ] ] - * x, normals [ i ]) );
    }
    unsigned oldi = numOfNodes - 1;
    Vector phiR(2);
    phiR [ 0 ] = 0;
    phiR [ 1 ] = 0;
    for ( unsigned i = 0; i < numOfNodes; i++ ) {
        R [ i ] [ 0 ] = normals [ oldi ].getX() / h [ oldi ] + normals [ i ].getX() / h [ i ];
        R [ i ] [ 1 ] = normals [ oldi ].getY() / h [ oldi ] + normals [ i ].getY() / h [ i ];
        phiR [ 0 ] += R [ i ] [ 0 ] * phi [ i ];
        phiR [ 1 ] += R [ i ] [ 1 ] * phi [ i ];
        oldi = i;
    }
    for ( unsigned i = 0; i < numOfNodes; i++ ) {
        for ( unsigned j = 0; j < 2; j++ ) {
            phiGrad [ j ] [ i ] = phi [ i ] * ( R [ i ] [ j ] - phiR [ j ] );
        }
    }
};
