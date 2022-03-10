#ifndef _SHAPE_F_H
#define _SHAPE_F_H

#include "node_container.h"
class IntegrationType; //forward declaration

//////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////
// SHAPE FUNCTIONS - MASTER CLASS
class ShapeFunc
{
private:

protected:
    unsigned ndim;
    string name;
    bool is_natural; //boolean determining if shape functions are computed in natural or real coordinates
    vector< Point * >points;

    virtual void giveShapeFGradNatural(const Point *x, MyMatrix &phiGradNat) const { ( void ) x; ( void ) phiGradNat; }; //without Jacobi transformation, natural coordinates
    void giveShapeFGrad(const Point *x, const MyMatrix &JacobiMInverse, MyMatrix &phiGrad) const;
    void giveJacobiM(const MyMatrix &phiGradNat, MyMatrix &JacobiM) const;
    void giveJacobiMInverse(const MyMatrix &JacobiM, MyMatrix &JacobiMInverse) const;
    double giveJacobian(const MyMatrix &JacobiM) const;

public:
    ShapeFunc() { name = "basic shape functions"; is_natural = true; }
    virtual ~ShapeFunc() {};
    virtual void init(vector< Node * > &nodes);
    virtual void init(vector< Point * >points);
    unsigned giveDimension() const { return ndim; }
    virtual void giveShapeF(const Point *x, MyVector &phi) const { ( void ) x; ( void ) phi; };
    virtual void giveShapeFGrad(const Point *x, MyMatrix &phiGrad) const;
    void giveJacobiM(const Point *x, MyMatrix &JacobiM) const;
    void giveJacobiMInverse(const Point *x, MyMatrix &JacobiMInverse) const;
    double giveJacobian(const Point *x) const;
};

//////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////
// 1D LINEAR SHAPE FUNCTIONS
class Linear1DLineShapeF : public ShapeFunc
{
protected:
    virtual void giveShapeFGradNatural(const Point *x, MyMatrix &phiGradNat) const;

public:
    Linear1DLineShapeF() { name = "1D linear shape functions for line element"; ndim = 1;  is_natural = true; };
    virtual ~Linear1DLineShapeF() {};
    virtual void giveShapeF(const Point *x, MyVector &phi) const;
};

//////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////
// 2D LINEAR SHAPE FUNCTIONS IN QUAD
class Linear2DQuadShapeF : public ShapeFunc
{
protected:
    virtual void giveShapeFGradNatural(const Point *x, MyMatrix &phiGradNat) const;

public:
    Linear2DQuadShapeF() { name = "2D linear shape functions for Quad"; ndim = 2;  is_natural = true; };
    virtual ~Linear2DQuadShapeF() {};
    virtual void giveShapeF(const Point *x, MyVector &phi) const;
};

//////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////
// 2D LINEAR SHAPE FUNCTIONS IN TRIANGLE
class Linear2DTriShapeF : public ShapeFunc
{
    double area;
public:
    Linear2DTriShapeF() { name = "2D linear shape functions for Triangle"; ndim = 2;  is_natural = false; };
    virtual ~Linear2DTriShapeF() {};
    virtual void init(vector< Node * > &nodes);
    virtual void init(vector< Point * > &points);
    virtual void giveShapeF(const Point *x, MyVector &phi) const;
    virtual void giveShapeFGrad(const Point *x, MyMatrix &phiGrad) const;
    double giveArea()const { return area; };
};

//////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////
// 2D LINEAR TRIANGULAR BASED SHAPE FUNCTIONS IN POLYGON
class Linear2DPolygonShapeF : public ShapeFunc
{
protected:
    Point centroid;
    vector< vector< unsigned > >faces;
    vector< double >angles;
    vector< Linear2DTriShapeF >triangles;
    IntegrationType *inttype;
    MyVector red2full;

    virtual void giveFullShapeF(const Point *x, MyVector &phi) const;
    virtual void giveFullShapeFGrad(const Point *x, MyMatrix &phiGrad) const;
public:
    Linear2DPolygonShapeF() { name = "2D linear shape functions for Polygon"; ndim = 2;  is_natural = false; };
    virtual ~Linear2DPolygonShapeF() {};
    virtual void init(vector< Node * > &nodes);
    virtual void giveShapeF(const Point *x, MyVector &phi) const;
    virtual void giveShapeFGrad(const Point *x, MyMatrix &phiGrad) const;
    void setFacesCentroidAndIntegration(vector< vector< unsigned > > &f, Point c, IntegrationType *it);
    unsigned findFaceNumber(const Point *x) const;
};


//////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////
// 2D WACHSPRESS
class Wachspress2DShapeF : public ShapeFunc
{
protected:
    unsigned nfaces;
    vector< vector< unsigned > >faces;
    vector< Point >normals;

public:
    Wachspress2DShapeF() { name = "2D Wachspress shape functions"; ndim = 2;  is_natural = false; };
    virtual ~Wachspress2DShapeF() {};
    virtual void giveShapeF(const Point *x, MyVector &phi) const;
    virtual void giveShapeFGrad(const Point *x, MyMatrix &phiGrad) const;
    void setFacesAndNormals(vector< vector< unsigned > > &f, vector< Point >n) { normals = n; faces = f; nfaces = faces.size(); }
};


//////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////
// 3D LINEAR SHAPE FUNCTIONS IN BRICK
class Linear3DBrickShapeF : public ShapeFunc
{
protected:
    virtual void giveShapeFGradNatural(const Point *x, MyMatrix &phiGradNat) const;

public:
    Linear3DBrickShapeF() { name = "3D linear shape functions for Brick"; ndim = 3;  is_natural = true; };
    virtual ~Linear3DBrickShapeF() {};
    virtual void giveShapeF(const Point *x, MyVector &phi) const;
};

#endif  /* _SHAPE_F_H */
