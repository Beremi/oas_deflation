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
    ~TrsprtQuad() {};
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
    ~TrsprtBrick() {};
};


//////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////
// 3D BRICK TRANSPORT + TEMPERATURE  ELEMENT
class TrsprtTemprtrCoupledBrick : public TrsprtBrick
{
protected:

public:
    TrsprtTemprtrCoupledBrick();
    ~TrsprtTemprtrCoupledBrick() {};
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
    ~MechanicalQuad() {};
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
    ~MechanicalBrick() {};
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
    ~CosseratQuad() {};
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
    ~CosseratBrick() {};
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
    ~CoupledCosseratQuad() {};
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
    vector< double >IPpressures;
    vector< double >IPvolstrains;

public:
    CoupledCosseratBrick();
    ~CoupledCosseratBrick() {};
    virtual Matrix giveBMatrix(const Point *x) const;
    virtual Matrix giveHMatrix(const Point *x) const;
    virtual Matrix giveDampingMatrix() const;
    virtual Vector giveStrain(unsigned i, const Vector &DoFs);
    double givePressureAtINtegrationPoint(unsigned i) const { return IPpressures [ i ]; };
    virtual Vector giveInternalForces(const Vector &DoFs, bool frozen, double timeStep);
};
#endif  /* _ELEMENT_STRUCT_H */
