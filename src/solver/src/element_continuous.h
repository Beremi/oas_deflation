#ifndef _ELEMENT_CONTINUOUS_H
#define _ELEMENT_CONTINUOUS_H

#include "element.h"

class ElementContainer; //forward declaration;
class BodyLoad; //forward declaration

//////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////
// 2D QUADRILATERAL TRANSPORT ELEMENT
class TrsprtQuad : public TransportElement
{
protected:
public:
    TrsprtQuad();
    virtual ~TrsprtQuad() {};
    virtual MyMatrix giveBMatrix(const Point *x) const;
    virtual MyMatrix giveHMatrix(const Point *x) const;
    virtual MyVector giveStrain(unsigned i, const MyVector &DoFs);
};

//////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////
// 3D BRICK TRANSPORT ELEMENT
class TrsprtBrick : public TrsprtQuad
{
protected:

public:
    TrsprtBrick();
    virtual ~TrsprtBrick() {};
};


//////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////
// 3D BRICK TRANSPORT + TEMPERATURE  ELEMENT
class TrsprtTemprtrCoupledBrick : public TrsprtBrick
{
protected:

public:
    TrsprtTemprtrCoupledBrick();
    virtual ~TrsprtTemprtrCoupledBrick() {};
    virtual MyMatrix giveBMatrix(const Point *x) const;
    virtual MyMatrix giveHMatrix(const Point *x) const;
    virtual MyVector giveStrain(unsigned i, const MyVector &DoFs);
};
//////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////
// 2D QUADRILATERAL MECHANICAL ELEMENT
class MechanicalQuad : public MechanicalElement
{
protected:
public:
    MechanicalQuad();
    virtual ~MechanicalQuad() {};
    virtual MyMatrix giveBMatrix(const Point *x) const;
    virtual MyMatrix giveHMatrix(const Point *x) const;
};


//////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////
// 3D BRICK MECHANICAL ELEMENT
class MechanicalBrick : public MechanicalQuad
{
protected:
public:
    MechanicalBrick();
    virtual ~MechanicalBrick() {};
    virtual MyMatrix giveBMatrix(const Point *x) const;
};


//////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////
// 2D QUADRILATERAL COSSERAT MECHANICAL ELEMENT
class CosseratQuad : public MechanicalQuad
{
protected:

public:
    CosseratQuad();
    virtual ~CosseratQuad() {};
    virtual MyMatrix giveBMatrix(const Point *x) const;
    virtual MyMatrix giveHMatrix(const Point *x) const;
};


//////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////
// 3D BRICK COSSERAT MECHANICAL ELEMENT
class CosseratBrick : public MechanicalBrick
{
protected:

public:
    CosseratBrick();
    virtual ~CosseratBrick() {};
    virtual MyMatrix giveBMatrix(const Point *x) const;
    virtual MyMatrix giveHMatrix(const Point *x) const;
};

//////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////
// 2D QUADRILATERAL COSSERAT COUPLED MECHANICAL-TRANSPORT ELEMENT
class CoupledCosseratQuad : public CosseratQuad
{
protected:

public:
    CoupledCosseratQuad();
    virtual ~CoupledCosseratQuad() {};
    virtual void init();
    virtual MyMatrix giveBMatrix(const Point *x) const;
    virtual MyMatrix giveHMatrix(const Point *x) const;
};

//////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////
// 3D BRICK COSSERAT COUPLED MECHANICAL-TRANSPORT ELEMENT
class CoupledCosseratBrick : public CosseratBrick
{
protected:

public:
    CoupledCosseratBrick();
    virtual ~CoupledCosseratBrick() {};
    virtual MyMatrix giveBMatrix(const Point *x) const;
    virtual MyMatrix giveHMatrix(const Point *x) const;
    virtual MyMatrix giveDampingMatrix() const;
    virtual MyVector giveStrain(unsigned i, const MyVector &DoFs);
};

//////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////
// 3D BRICK COSSERAT COUPLED MECHANICAL-TRANSPORT ELEMENT
class CoupledCosseratBrickWithDependentUpperZLayer : public CoupledCosseratBrick
{
protected:
    bool bindlayers;

public:
    CoupledCosseratBrickWithDependentUpperZLayer();
    virtual ~CoupledCosseratBrickWithDependentUpperZLayer() {};
    virtual void setIntegrationPointsAndWeights();
    virtual MyMatrix giveStiffnessMatrix(std :: string matrixType) const;
    virtual MyVector giveInternalForces(const MyVector &DoFs, bool frozen, double timeStep);
    virtual MyVector integrateInternalSources();
    virtual MyMatrix giveDampingMatrix() const;
    virtual MyMatrix giveMassMatrix() const;
    virtual MyVector giveStrain(unsigned i, const MyVector &DoFs);
};
#endif  /* _ELEMENT_STRUCT_H */
