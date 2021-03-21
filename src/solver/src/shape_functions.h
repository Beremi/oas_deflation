#ifndef _SHAPE_F_H
#define _SHAPE_F_H

#include "node_container.h"

//////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////
// SHAPE FUNCTIONS - MASTER CLASS
class ShapeFunc
{
private:

protected:
    unsigned ndim;
    string name;   
    
public:
    ShapeFunc() { name = "basic shape functions"; }
    virtual ~ShapeFunc(){};
    virtual void init(){};
    unsigned giveDimension() const { return ndim; }
    virtual void giveShapeF(const Point *x, Vector &phi ) const { (void) x; (void) phi; };
    virtual void giveShapeFGradNatural(const Point *x, Matrix &phiGradNat) const { ( void ) x; ( void ) phiGradNat;}; //without Jacobi transformation, natural coordinates
    void giveShapeFGrad(const Point *x, const vector<Node*> & nodes, Matrix &phiGrad) const;
    void giveShapeFGrad(const Point *x, const Matrix &JacobiMInverse, Matrix &phiGrad) const;
    void giveJacobiM(const Point *x, const vector<Node*> & nodes, const Matrix &phiGradNat, Matrix &JacobiM) const;
    void giveJacobiM(const Point *x, const vector<Node*> & nodes, Matrix &JacobiM) const;
    void giveJacobiMInverse(const Matrix &JacobiM, Matrix &JacobiMInverse) const;
    void giveJacobiMInverse(const Point *x, const vector<Node*> & nodes, Matrix &JacobiMInverse) const;
    double giveJacobian(const Point *x, const vector<Node*> & nodes) const;
    double giveJacobian(const Matrix &JacobiM) const;
};

//////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////
// 1D LINEAR SHAPE FUNCTIONS
class Linear1DLineShapeF : public ShapeFunc
{
public:
    Linear1DLineShapeF() { name = "1D linear shape functions for line element"; ndim = 1;};
    virtual ~Linear1DLineShapeF(){};
    virtual void giveShapeF(const Point *x, Vector &phi) const;
    virtual void giveShapeFGradNatural(const Point *x, Matrix &phiGradNat) const;
};

//////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////
// 2D LINEAR SHAPE FUNCTIONS
class Linear2DQuadShapeF : public ShapeFunc
{
public:
    Linear2DQuadShapeF() { name = "2D linear shape functions for Quad"; ndim = 2;};
    virtual ~Linear2DQuadShapeF(){};
    virtual void giveShapeF(const Point *x, Vector &phi) const;
    virtual void giveShapeFGradNatural(const Point *x, Matrix &phiGradNat) const;
};


//////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////
// 3D LINEAR SHAPE FUNCTIONS
class Linear3DBrickShapeF : public ShapeFunc
{
public:
    Linear3DBrickShapeF() { name = "3D linear shape functions for Brick"; ndim = 3;};
    virtual ~Linear3DBrickShapeF(){};
    virtual void giveShapeF(const Point *x, Vector &phi) const;
    virtual void giveShapeFGradNatural(const Point *x, Matrix &phiGradNat) const;
};

#endif  /* _SHAPE_F_H */
