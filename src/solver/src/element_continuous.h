#ifndef _ELEMENT_CONTINUOUS_H
#define _ELEMENT_CONTINUOUS_H

#include "element.h"

class ElementContainer; //forward declaration;
class BodyLoad; //forward declaration

//////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////
// 2D TRIANGULAR TRANSPORT ELEMENT
class TrsprtTriangle : public Element
{
protected:
public:
    TrsprtTriangle();
    virtual ~TrsprtTriangle() {};
    virtual Matrix giveBMatrix(const Point *x) const;
    virtual Matrix giveHMatrix(const Point *x) const;
    virtual Vector giveStrain(unsigned i, const Vector &DoFs);
};

//////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////
// 2D QUADRILATERAL TRANSPORT ELEMENT
class TrsprtQuad : public TrsprtTriangle
{
protected:
public:
    TrsprtQuad();
    virtual ~TrsprtQuad() {};
};

//////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////
// 3D TETRA TRANSPORT ELEMENT
class TrsprtTetra : public TrsprtTriangle
{
protected:

public:
    TrsprtTetra();
    virtual ~TrsprtTetra() {};
    virtual void init();
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
// 2D TRIANGULAR MECHANICAL ELEMENT
class MechanicalTriangle : public Element
{
protected:
    Matrix averageVolumeB;
    bool b_bar_integration_split;
    void applyAverageVolumeB(Matrix &B, const Matrix phiG) const;
    virtual void computeAverageBVolumeMatrix();
    virtual void setIntegrationPointsAndWeights();

public:
    MechanicalTriangle();
    virtual ~MechanicalTriangle() {};
    virtual Matrix giveBMatrix(const Point *x) const;
    virtual Matrix giveHMatrix(const Point *x) const;
    virtual void readFromLine(std :: istringstream &iss, NodeContainer *fullnodes, MaterialContainer *fullmatrs);
};

//////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////
// 2D QUADRILATERAL MECHANICAL ELEMENT
class MechanicalQuad : public MechanicalTriangle
{
protected:
public:
    MechanicalQuad();
    virtual ~MechanicalQuad() {};
};

//////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////
// 3D TETRAHEDRAL MECHANICAL ELEMENT
class MechanicalTetra : public MechanicalTriangle
{
protected:
public:
    MechanicalTetra();
    virtual ~MechanicalTetra() {};
    virtual Matrix giveBMatrix(const Point *x) const;
    virtual void init();
};

//////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////
// 3D BRICK MECHANICAL ELEMENT
class MechanicalBrick : public MechanicalTetra
{
protected:
public:
    MechanicalBrick();
    virtual ~MechanicalBrick() {};
    virtual void init();
};


//////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////
// 2D QUADRILATERAL COSSERAT MECHANICAL ELEMENT
class CosseratQuad : public MechanicalQuad
{
protected:
    virtual void computeMassMatrix();
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
class CoupledCosseratTransportQuad : public CosseratQuad
{
protected:

public:
    CoupledCosseratTransportQuad();
    virtual ~CoupledCosseratTransportQuad() {};
    virtual void init();
    virtual Matrix giveBMatrix(const Point *x) const;
    virtual Matrix giveHMatrix(const Point *x) const;
};

//////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////
// 3D BRICK COSSERAT COUPLED MECHANICAL-TRANSPORT ELEMENT
class CoupledCosseratTransportBrick : public CosseratBrick
{
protected:
    virtual void computeDampingMatrix();
public:
    CoupledCosseratTransportBrick();
    virtual ~CoupledCosseratTransportBrick() {};
    virtual Matrix giveBMatrix(const Point *x) const;
    virtual Matrix giveHMatrix(const Point *x) const;
    virtual Vector giveStrain(unsigned i, const Vector &DoFs);
};

//////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////
// 3D BRICK COSSERAT COUPLED MECHANICAL-TRANSPORT ELEMENT
class CoupledCosseratBrickWithDependentUpperZLayer : public CoupledCosseratTransportBrick
{
protected:
    bool bindlayers;
    virtual void computeMassMatrix();
    virtual void computeDampingMatrix();
public:
    CoupledCosseratBrickWithDependentUpperZLayer();
    virtual ~CoupledCosseratBrickWithDependentUpperZLayer() {};
    virtual void setIntegrationPointsAndWeights();
    virtual Matrix giveStiffnessMatrix(std :: string matrixType) const;
    virtual Vector giveInternalForces(const Vector &DoFs, bool frozen, double timeStep);
    virtual Vector integrateInternalSources();
    virtual Vector giveStrain(unsigned i, const Vector &DoFs);
};
#endif  /* _ELEMENT_CONTINUOUS_H */
