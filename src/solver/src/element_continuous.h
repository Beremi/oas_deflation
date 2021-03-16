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
    virtual void setIntegrationPointsAndWeights();
public:
    TrsprtQuad();
    ~TrsprtQuad() {};
    void readFromLine(istringstream &iss, NodeContainer *fullnodes, MaterialContainer *fullmatrs);
    virtual void shapeF(const Point *x, Vector &phi) const;
    virtual double shapeFGrad(const Point *x, Matrix &phiGrad) const;
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
    virtual void setIntegrationPointsAndWeights();

public:
    TrsprtBrick();
    ~TrsprtBrick() {};
    void readFromLine(istringstream &iss, NodeContainer *fullnodes, MaterialContainer *fullmatrs);
    virtual void shapeF(const Point *x, Vector &phi) const;
    virtual double shapeFGrad(const Point *x, Matrix &phiGrad) const;
};

//////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////
// 2D QUADRILATERAL MECHANICAL ELEMENT
class MechanicalQuad : public TrsprtQuad
{
protected:
public:
    MechanicalQuad();
    ~MechanicalQuad() {};
    virtual Matrix giveBMatrix(const Point *x) const;
    virtual Matrix giveHMatrix(const Point *x) const;
};

/*
//////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////
// 3D BRICK MECHANICAL ELEMENT
class MechanicalBrick : public MechanicalQuad
{
protected:
    virtual void setIntegrationPointsAndWeights();

public:
    MechanicalBrick();
    ~MechanicalBrick() {};
    void readFromLine(istringstream &iss, NodeContainer *fullnodes, MaterialContainer *fullmatrs);
    virtual void shapeF(const Point *x, Vector &phi) const;
    virtual double shapeFGrad(const Point *x, Matrix &phiGrad) const;
    virtual Matrix giveBMatrix(const Point *x) const;
    virtual Matrix giveHMatrix(const Point *x) const;
};
*/

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

/*
//////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////
// 3D BRICK COSSERAT MECHANICAL ELEMENT
class CosseratBrick : public MechanicalQuad
{
protected:

public:
    CosseratBrick();
    ~CosseratBrick() {};
    virtual Matrix giveBMatrix(const Point *x) const;
    virtual Matrix giveHMatrix(const Point *x) const;
};
*/
#endif  /* _ELEMENT_STRUCT_H */
