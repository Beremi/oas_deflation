#ifndef _MAT_RVE_H
#define _MAT_RVE_H

#include "material.h"
#include "globals.h"


class Model; //forward declaraion
class Node; //forward declaraion
class Transp1D; //forward declaraion
class BoundaryCondition; //forward declaraion
class PieceWiseLinearFunction; //forward declaraion



//////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////
// GENERAL RVE MATERIAL

class RVEMaterial;
class RVEMaterialStatus : public MaterialStatus
{
protected:
    Model *RVE;

    fs :: path inputfile;

    //setup for volumetric average
    PieceWiseLinearFunction *volumAverFunc;

    virtual void generateVolumetricAverageBC() {};
public:
    RVEMaterialStatus(RVEMaterial *m, Element *e, fs :: path masterfile);
    virtual ~RVEMaterialStatus();
    virtual void init();
    virtual void update();
};

//////////////////////////////////////////////////////////
class RVEMaterial : public Material
{
protected:
    fs :: path inputfile;


public:
    RVEMaterial() { name = "generic RVE material"; };
    virtual ~RVEMaterial() {};
    virtual void readFromLine(istringstream &iss);
    virtual MaterialStatus *giveNewMaterialStatus(Element *e);
};


//////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////
// DISCRETE RVE MATERIAL FOR BOTH TRANSPORT AND MECHANICAL HOMOGENIZATION

class DiscreteRVEMaterial;
class DiscreteRVEMaterialStatus : public RVEMaterialStatus
{
protected:
    virtual void generateVolumetricAverageBC();
    void calculateCentroid();

    Point centroid;
    bool active_mechanics, active_transport;
    vector< vector< Vector > >mechProjectors;

public:
    DiscreteRVEMaterialStatus(DiscreteRVEMaterial *m, Element *e, fs :: path masterfile);
    virtual ~DiscreteRVEMaterialStatus() {};
    virtual void init();
    virtual Vector giveStress(const Vector &strain);//terminology from mechanics, it returns flux
    virtual Matrix giveStiffnessTensor(string type, unsigned dimension) const;
    virtual double giveMassConstant() const;
};

//////////////////////////////////////////////////////////
class DiscreteRVEMaterial : public RVEMaterial
{
protected:

public:
    DiscreteRVEMaterial() { name = "transport RVE material"; };
    virtual ~DiscreteRVEMaterial() {};
    virtual MaterialStatus *giveNewMaterialStatus(Element *e);
};

#endif /* _MAT_RVE_H */
