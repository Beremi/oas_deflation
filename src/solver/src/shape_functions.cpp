#include "shape_functions.h"
#include "integration.h"

using namespace std;

//////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////
// SHAPE FUNCTIONS - MASTER CLASS
//////////////////////////////////////////////////////////
void ShapeFunc :: init(vector< Node * > &nodes) {
    points.resize( nodes.size() );
    unsigned k = 0;
    for ( auto &n: nodes ) {
        points [ k ] = n->givePointPointer();
        k++;
    }
}

//////////////////////////////////////////////////////////
void ShapeFunc :: init(vector< Point * > &pp) {
    points = pp;
}

//////////////////////////////////////////////////////////
void ShapeFunc :: giveShapeFGrad(const Point *x, Matrix &phiGrad) const {
    giveShapeFGradNatural(x, phiGrad);
    Matrix JacobiMInverse = Matrix :: Zero(ndim, ndim);
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
    JacobiM.setZero();
    Point *n;
    for ( unsigned i = 0; i < points.size(); i++ ) {
        n = points [ i ];
        for ( unsigned v1 = 0; v1 < ndim; v1++ ) {
            for ( unsigned v2 = 0; v2 < ndim; v2++ ) {
                JacobiM(v1, v2) += phiGradNat(v1, i) * ( * n )(v2);
            }
        }
    }
}

//////////////////////////////////////////////////////////
void ShapeFunc :: giveJacobiM(const Point *x, Matrix &JacobiM) const {
    Matrix phiGradNat = Matrix :: Zero(ndim, points.size() );
    giveShapeFGradNatural(x, phiGradNat);
    JacobiM.setZero();
    Point *n;
    for ( unsigned i = 0; i < points.size(); i++ ) {
        n = points [ i ];
        for ( unsigned v1 = 0; v1 < ndim; v1++ ) {
            for ( unsigned v2 = 0; v2 < ndim; v2++ ) {
                JacobiM(v1, v2) += phiGradNat(v1, i) * ( * n )(v2);
            }
        }
    }
}

//////////////////////////////////////////////////////////
double ShapeFunc :: giveJacobian(const Point *x) const {
    Matrix JacobiM = Matrix :: Zero(ndim, ndim);
    giveJacobiM(x, JacobiM);
    return giveJacobian(JacobiM);
}

//////////////////////////////////////////////////////////
double ShapeFunc :: giveJacobian(const Matrix &JacobiM) const {
    if ( ndim == 1 ) {
        return JacobiM(0, 0);
    } else if ( ndim == 2 ) {
        return JacobiM(0, 0) * JacobiM(1, 1) - JacobiM(0, 1) * JacobiM(1, 0);
    } else if ( ndim == 3 ) {
        return JacobiM(0, 0) * JacobiM(1, 1) * JacobiM(2, 2) - JacobiM(0, 0) * JacobiM(1, 2) *
               JacobiM(2, 1) + JacobiM(0, 1) * JacobiM(1, 2) * JacobiM(2, 0) - JacobiM(0, 1) *
               JacobiM(1, 0) * JacobiM(2, 2) + JacobiM(0, 2) * JacobiM(1, 0) * JacobiM(2, 1) -
               JacobiM(0, 2) * JacobiM(1, 1) * JacobiM(2, 0);
    } else {
        return 0;
    }
}

//////////////////////////////////////////////////////////
void ShapeFunc :: giveJacobiMInverse(const Matrix &JacobiM, Matrix &JacobiMInverse) const {
    double JacDet = giveJacobian(JacobiM);
    if ( ndim == 1 ) {
        JacobiMInverse(0, 0) = 1. / JacobiM(0, 0);
    } else if ( ndim == 2 ) {
        JacobiMInverse(0, 0) = JacobiM(1, 1) / JacDet;
        JacobiMInverse(0, 1) = -JacobiM(0, 1) / JacDet;
        JacobiMInverse(1, 0) = -JacobiM(1, 0) / JacDet;
        JacobiMInverse(1, 1) = JacobiM(0, 0) / JacDet;
    } else if ( ndim == 3 ) {
        JacobiMInverse(0, 0) =  ( JacobiM(1, 1) * JacobiM(2, 2) - JacobiM(2, 1) * JacobiM(1, 2) ) / JacDet;
        JacobiMInverse(0, 1) = -( JacobiM(0, 1) * JacobiM(2, 2) - JacobiM(2, 1) * JacobiM(0, 2) ) / JacDet;
        JacobiMInverse(0, 2) =  ( JacobiM(0, 1) * JacobiM(1, 2) - JacobiM(1, 1) * JacobiM(0, 2) ) / JacDet;
        JacobiMInverse(1, 0) = -( JacobiM(1, 0) * JacobiM(2, 2) - JacobiM(2, 0) * JacobiM(1, 2) ) / JacDet;
        JacobiMInverse(1, 1) =  ( JacobiM(0, 0) * JacobiM(2, 2) - JacobiM(2, 0) * JacobiM(0, 2) ) / JacDet;
        JacobiMInverse(1, 2) = -( JacobiM(0, 0) * JacobiM(1, 2) - JacobiM(1, 0) * JacobiM(0, 2) ) / JacDet;
        JacobiMInverse(2, 0) =  ( JacobiM(1, 0) * JacobiM(2, 1) - JacobiM(2, 0) * JacobiM(1, 1) ) / JacDet;
        JacobiMInverse(2, 1) = -( JacobiM(0, 0) * JacobiM(2, 1) - JacobiM(2, 0) * JacobiM(0, 1) ) / JacDet;
        JacobiMInverse(2, 2) =  ( JacobiM(0, 0) * JacobiM(1, 1) - JacobiM(1, 0) * JacobiM(0, 1) ) / JacDet;
    }
}

//////////////////////////////////////////////////////////
void ShapeFunc :: giveJacobiMInverse(const Point *x, Matrix &JacobiMInverse) const {
    Matrix JacobiM = Matrix :: Zero(ndim, ndim);
    giveJacobiM(x, JacobiM);
    giveJacobiMInverse(JacobiM, JacobiMInverse);
}

//////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////
// 2D WACHSPRESS
//////////////////////////////////////////////////////////
void NullShapeF :: giveShapeF(const Point *x, Vector &phi) const {
    ( void ) x;
    ( void ) phi;
};

//////////////////////////////////////////////////////////
void NullShapeF :: giveShapeFGrad(const Point *x, Matrix &phiGrad) const {
    ( void ) x;
    ( void ) phiGrad;
};

//////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////
//  1D LINEAR SHAPE FUNCTIONS
//////////////////////////////////////////////////////////
void Linear1DLineShapeF :: giveShapeF(const Point *x, Vector &phi) const {
    //x in natural coordinates
    phi [ 0 ] = 0.5 * ( 1. - x->x() );
    phi [ 1 ] = 0.5 * ( 1. + x->x() );
}

//////////////////////////////////////////////////////////
void Linear1DLineShapeF :: giveShapeFGradNatural(const Point *x, Matrix &phiGrad) const {
    ( void ) x;
    //x in natural coordinates
    phiGrad(0, 0) = -0.5;
    phiGrad(0, 1) = 0.5;
}


//////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////
//  2D LINEAR SHAPE FUNCTIONS FOR QUAD
//////////////////////////////////////////////////////////
void Linear2DQuadShapeF :: giveShapeF(const Point *x, Vector &phi) const {
    //x in natural coordinates
    phi [ 0 ] = 0.25 * ( 1. - x->x() ) * ( 1. - x->y() );
    phi [ 1 ] = 0.25 * ( 1. + x->x() ) * ( 1. - x->y() );
    phi [ 2 ] = 0.25 * ( 1. + x->x() ) * ( 1. + x->y() );
    phi [ 3 ] = 0.25 * ( 1. - x->x() ) * ( 1. + x->y() );
}

//////////////////////////////////////////////////////////
void Linear2DQuadShapeF :: giveShapeFGradNatural(const Point *x, Matrix &phiGrad) const {
    //x in natural coordinates
    phiGrad(0, 0) = -0.25 * ( 1. - x->y() );
    phiGrad(0, 1) = -phiGrad(0, 0);
    phiGrad(0, 2) = 0.25 * ( 1. + x->y() );
    phiGrad(0, 3) = -phiGrad(0, 2);
    phiGrad(1, 0) = -0.25 * ( 1. - x->x() );
    phiGrad(1, 1) = -0.25 * ( 1. + x->x() );
    phiGrad(1, 2) = -phiGrad(1, 1);
    phiGrad(1, 3) = -phiGrad(1, 0);
}

//////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////
//  3D LINEAR SHAPE FUNCTIONS FOR BRICK
//////////////////////////////////////////////////////////
void Linear3DBrickShapeF :: giveShapeF(const Point *x, Vector &phi) const {
    //x in natural coordinates
    phi [ 0 ] = 0.125 * ( 1. - x->x() ) * ( 1. - x->y() )  * ( 1. - x->z() );
    phi [ 1 ] = 0.125 * ( 1. + x->x() ) * ( 1. - x->y() )  * ( 1. - x->z() );
    phi [ 2 ] = 0.125 * ( 1. + x->x() ) * ( 1. + x->y() )  * ( 1. - x->z() );
    phi [ 3 ] = 0.125 * ( 1. - x->x() ) * ( 1. + x->y() )  * ( 1. - x->z() );
    phi [ 4 ] = 0.125 * ( 1. - x->x() ) * ( 1. - x->y() )  * ( 1. + x->z() );
    phi [ 5 ] = 0.125 * ( 1. + x->x() ) * ( 1. - x->y() )  * ( 1. + x->z() );
    phi [ 6 ] = 0.125 * ( 1. + x->x() ) * ( 1. + x->y() )  * ( 1. + x->z() );
    phi [ 7 ] = 0.125 * ( 1. - x->x() ) * ( 1. + x->y() )  * ( 1. + x->z() );
}

//////////////////////////////////////////////////////////
void Linear3DBrickShapeF :: giveShapeFGradNatural(const Point *x, Matrix &phiGrad) const {
    //x in natural coordinates
    phiGrad(0, 0) = -0.125 * ( 1. - x->y() )  * ( 1. - x->z() );
    phiGrad(0, 1) = -phiGrad(0, 0);
    phiGrad(0, 2) =  0.125 * ( 1. + x->y() )  * ( 1. - x->z() );
    phiGrad(0, 3) = -phiGrad(0, 2);
    phiGrad(0, 4) = -0.125 * ( 1. - x->y() )  * ( 1. + x->z() );
    phiGrad(0, 5) = -phiGrad(0, 4);
    phiGrad(0, 6) =  0.125 * ( 1. + x->y() )  * ( 1. + x->z() );
    phiGrad(0, 7) = -phiGrad(0, 6);
    phiGrad(1, 0) = -0.125 * ( 1. - x->x() )  * ( 1. - x->z() );
    phiGrad(1, 1) = -0.125 * ( 1. + x->x() )  * ( 1. - x->z() );
    phiGrad(1, 2) = -phiGrad(1, 1);
    phiGrad(1, 3) = -phiGrad(1, 0);
    phiGrad(1, 4) = -0.125 * ( 1. - x->x() )  * ( 1. + x->z() );
    phiGrad(1, 5) = -0.125 * ( 1. + x->x() )  * ( 1. + x->z() );
    phiGrad(1, 6) = -phiGrad(1, 5);
    phiGrad(1, 7) = -phiGrad(1, 4);
    phiGrad(2, 0) = -0.125 * ( 1. - x->x() )  * ( 1. - x->y() );
    phiGrad(2, 1) = -0.125 * ( 1. + x->x() )  * ( 1. - x->y() );
    phiGrad(2, 2) = -0.125 * ( 1. + x->x() )  * ( 1. + x->y() );
    phiGrad(2, 3) = -0.125 * ( 1. - x->x() )  * ( 1. + x->y() );
    phiGrad(2, 4) = -phiGrad(2, 0);
    phiGrad(2, 5) = -phiGrad(2, 1);
    phiGrad(2, 6) = -phiGrad(2, 2);
    phiGrad(2, 7) = -phiGrad(2, 3);
}


//////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////
// 2D TRIANGLE LINEAR
//////////////////////////////////////////////////////////
void Linear2DTriShapeF :: giveShapeF(const Point *x, Vector &phi) const {
    phi [ 0 ] = 1. - x->x() - x->y();
    phi [ 1 ] = x->x();
    phi [ 2 ] = x->y();
};

//////////////////////////////////////////////////////////
void Linear2DTriShapeF :: giveShapeFGradNatural(const Point *x, Matrix &phiGrad) const {
    ( void ) x;
    phiGrad(0, 0) = -1.;
    phiGrad(1, 0) = -1.;
    phiGrad(0, 1) = 1.;
    phiGrad(1, 1) = 0.;
    phiGrad(0, 2) = 0.;
    phiGrad(1, 2) = 1.;
};

//////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////
//  3D LINEAR SHAPE FUNCTIONS FOR TETRAHEDRON
//////////////////////////////////////////////////////////
void Linear3DTetraShapeF :: giveShapeF(const Point *x, Vector &phi) const {
    phi [ 0 ] = 1. - x->x() - x->y() - x->z();
    phi [ 1 ] = x->x();
    phi [ 2 ] = x->y();
    phi [ 3 ] = x->z();
};

//////////////////////////////////////////////////////////
void Linear3DTetraShapeF :: giveShapeFGradNatural(const Point *x, Matrix &phiGrad) const {
    ( void ) x;
    phiGrad(0, 0) = -1.;
    phiGrad(1, 0) = -1.;
    phiGrad(2, 0) = -1.;
    phiGrad(0, 1) = 1.;
    phiGrad(1, 1) = 0.;
    phiGrad(2, 1) = 0.;
    phiGrad(0, 2) = 0.;
    phiGrad(1, 2) = 1.;
    phiGrad(2, 2) = 0.;
    phiGrad(0, 3) = 0.;
    phiGrad(1, 3) = 0.;
    phiGrad(2, 3) = 1.;
};

//////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////
// 2D LINEAR TRIANGULAR BASED SHAPE FUNCTIONS IN POLYGON
//////////////////////////////////////////////////////////
void Linear2DPolygonShapeF :: init(vector< Node * > &nodes) {
    ShapeFunc :: init(nodes);

    unsigned n = points.size();
    angles.resize( points.size() );
    for ( unsigned i = 0; i < points.size(); i++ ) {
        angles [ i ] = atan2(points [ i ]->y() - centroid.y(), points [ i ]->x() - centroid.x() );
    }
    triangles.resize( faces.size() );
    vector< Point * >trinodes;
    trinodes.resize(3);
    trinodes [ 2 ] = & centroid;
    for ( unsigned i = 0; i < faces.size(); i++ ) {
        trinodes [ 0 ] = points [ faces [ i ] [ 0 ] ];
        trinodes [ 1 ] = points [ faces [ i ] [ 1 ] ];
        triangles [ i ].init(trinodes);
    }


    Matrix FullK = Matrix :: Zero(n + 1, n + 1);
    Matrix phiFullGrad = Matrix :: Zero(2, n + 1);
    for ( unsigned i = 0; i < inttype->giveNumIP(); i++ ) {
        giveFullShapeFGrad(inttype->giveIPLocationPointer(i), phiFullGrad);
        FullK += ( phiFullGrad.transpose() * phiFullGrad ) * inttype->giveIPWeight(i);
    }

    red2full.resize(n);
    for ( unsigned i = 0; i < n; i++ ) {
        red2full [ i ] = -FullK(i, n) / FullK(n, n);
    }
}


//////////////////////////////////////////////////////////
unsigned Linear2DPolygonShapeF :: findFaceNumber(const Point *x) const {
    double alpha = atan2(x->y() - centroid.y(), x->x() - centroid.x() );
    double a, b;
    unsigned face;
    for ( face = 0; face < faces.size(); face++ ) {
        a = angles [ faces [ face ] [ 0 ] ];
        b = angles [ faces [ face ] [ 1 ] ];
        if ( a < b && a < alpha && b > alpha ) {
            return face;
        }
        if ( a > b && ( a < alpha || b > alpha ) ) {
            return face;
        }
    }
    cerr << "Error in Linear2DPolygonShapeF: findFaceNumber should never go here" << endl;
    exit(1);
}

//////////////////////////////////////////////////////////
void Linear2DPolygonShapeF :: setFacesCentroidAndIntegration(vector< vector< unsigned > > &f, Point c, IntegrationType *it) {
    faces = f;
    centroid = c;
    inttype = it;
}

//////////////////////////////////////////////////////////
void Linear2DPolygonShapeF :: giveFullShapeF(const Point *x, Vector &phi) const {
    unsigned t = findFaceNumber(x);
    Vector triphi = Vector :: Zero(3);
    triangles [ t ].giveShapeF(x, triphi);
    phi.setZero();
    phi [ faces [ t ] [ 0 ] ]  = triphi [ 0 ];
    phi [ faces [ t ] [ 1 ] ]  = triphi [ 1 ];
    phi [ faces.size() ] = triphi [ 2 ]; //centroid
}

//////////////////////////////////////////////////////////
void Linear2DPolygonShapeF :: giveFullShapeFGrad(const Point *x, Matrix &phiGrad) const {
    unsigned t = findFaceNumber(x);
    Matrix triphiGrad = Matrix :: Zero(2, 3);
    triangles [ t ].giveShapeFGrad(x, triphiGrad);
    phiGrad.setZero();
    for ( unsigned i = 0; i < 2; i++ ) {
        phiGrad(i, faces [ t ] [ 0 ])  = triphiGrad(i, 0);
        phiGrad(i, faces [ t ] [ 1 ])  = triphiGrad(i, 1);
        phiGrad(i, faces.size() ) = triphiGrad(i, 2);       //centroid
    }
};


//////////////////////////////////////////////////////////
void Linear2DPolygonShapeF :: giveShapeF(const Point *x, Vector &phi) const {
    unsigned t = findFaceNumber(x);
    Vector triphi = Vector :: Zero(3);
    triangles [ t ].giveShapeF(x, triphi);
    phi.setZero();
    phi [ faces [ t ] [ 0 ] ]  = triphi [ 0 ];
    phi [ faces [ t ] [ 1 ] ]  = triphi [ 1 ];
    for ( unsigned r = 0; r < red2full.size(); r++ ) {
        phi [ r ] += triphi [ 2 ] * red2full [ r ];
    }
}

//////////////////////////////////////////////////////////
void Linear2DPolygonShapeF :: giveShapeFGrad(const Point *x, Matrix &phiGrad) const {
    unsigned t = findFaceNumber(x);
    Matrix triphiGrad = Matrix :: Zero(2, 3);
    triangles [ t ].giveShapeFGrad(x, triphiGrad);
    phiGrad.setZero();
    for ( unsigned i = 0; i < 2; i++ ) {
        phiGrad(i, faces [ t ] [ 0 ])  = triphiGrad(i, 0);
        phiGrad(i, faces [ t ] [ 1 ])  = triphiGrad(i, 1);

        for ( unsigned r = 0; r < red2full.size(); r++ ) {
            phiGrad(i, r) += triphiGrad(i, 2) * red2full [ r ];
        }
    }
};



//////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////
// 2D WACHSPRESS
//////////////////////////////////////////////////////////
void Wachspress2DShapeF :: giveShapeF(const Point *x, Vector &phi) const {
    unsigned numOfNodes = points.size();
    Vector h = Vector :: Zero(numOfNodes);
    Point oldNormal = normals [ numOfNodes - 1 ];
    for ( unsigned i = 0; i < numOfNodes; i++ ) {
        h [ i ] = abs( ( * points [ faces [ i ] [ 0 ] ] - * x ).dot(normals [ i ]) );
        phi [ i ] = abs(oldNormal.x() * normals [ i ].y() - oldNormal.y() * normals [ i ].x() );
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
    Vector phi = Vector :: Zero(numOfNodes);
    giveShapeF(x, phi);
    Matrix R = Matrix :: Zero(numOfNodes, 2);

    Vector h = Vector :: Zero(numOfNodes);
    for ( unsigned i = 0; i < numOfNodes; i++ ) {
        h [ i ] = abs( ( * points [ faces [ i ] [ 0 ] ] - * x ).dot(normals [ i ]) );
    }
    unsigned oldi = numOfNodes - 1;
    Vector phiR = Vector :: Zero(2);
    phiR [ 0 ] = 0;
    phiR [ 1 ] = 0;
    for ( unsigned i = 0; i < numOfNodes; i++ ) {
        R(i, 0) = normals [ oldi ].x() / h [ oldi ] + normals [ i ].x() / h [ i ];
        R(i, 1) = normals [ oldi ].y() / h [ oldi ] + normals [ i ].y() / h [ i ];
        phiR [ 0 ] += R(i, 0) * phi [ i ];
        phiR [ 1 ] += R(i, 1) * phi [ i ];
        oldi = i;
    }
    for ( unsigned i = 0; i < numOfNodes; i++ ) {
        for ( unsigned j = 0; j < 2; j++ ) {
            phiGrad(j, i) = phi [ i ] * ( R(i, j) - phiR [ j ] );
        }
    }
};
