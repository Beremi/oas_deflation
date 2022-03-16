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
    virtual Matrix giveBMatrix(const Point *x) const;
    virtual Matrix giveHMatrix(const Point *x) const;
    virtual Vector giveStrain(unsigned i, const Vector &DoFs);
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
    virtual Matrix giveBMatrix(const Point *x) const;
    virtual Matrix giveHMatrix(const Point *x) const;
    virtual Vector giveStrain(unsigned i, const Vector &DoFs);
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
    virtual Matrix giveBMatrix(const Point *x) const;
    virtual Matrix giveHMatrix(const Point *x) const;
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
    virtual Matrix giveBMatrix(const Point *x) const;
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
    virtual Matrix giveBMatrix(const Point *x) const;
    virtual Matrix giveHMatrix(const Point *x) const;
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
    virtual Matrix giveBMatrix(const Point *x) const;
    virtual Matrix giveHMatrix(const Point *x) const;
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
    virtual Matrix giveBMatrix(const Point *x) const;
    virtual Matrix giveHMatrix(const Point *x) const;
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
    virtual Matrix giveBMatrix(const Point *x) const;
    virtual Matrix giveHMatrix(const Point *x) const;
    virtual Matrix giveDampingMatrix() const;
    virtual Vector giveStrain(unsigned i, const Vector &DoFs);
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
    virtual Matrix giveStiffnessMatrix(std :: string matrixType) const;
    virtual Vector giveInternalForces(const Vector &DoFs, bool frozen, double timeStep);
    virtual Vector integrateInternalSources();
    virtual Matrix giveDampingMatrix() const;
    virtual Matrix giveMassMatrix() const;
    virtual Vector giveStrain(unsigned i, const Vector &DoFs);
};
#endif  /* _ELEMENT_CONTINUOUS_H */
