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
    virtual Vector giveStressWithFrozenIntVars(const Vector &strain);
    virtual Matrix giveStiffnessTensor(string type, unsigned dimension) const;
    virtual double giveDampingConstant() const;
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

//////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////
// DISCRETE RVE MATERIAL FOR TRANSPORT PRECOMPUTED

class DiscreteRVEMaterialPrecomputed;
class DiscreteRVEMaterialPrecomputedStatus : public DiscreteRVEMaterialStatus
{
protected:
    double temp_nonlin;
    bool is_master_status;

public:
    DiscreteRVEMaterialPrecomputedStatus(DiscreteRVEMaterialPrecomputed *m, Element *e, fs :: path masterfile);
    virtual ~DiscreteRVEMaterialPrecomputedStatus() {};
    virtual void init();
    virtual Vector giveStress(const Vector &strain);//terminology from mechanics, it returns flux
    virtual Vector giveStressWithFrozenIntVars(const Vector &strain);
    virtual Matrix giveStiffnessTensor(string type, unsigned dimension) const;
    virtual double giveDampingConstant() const;
    virtual void update(); 
};

//////////////////////////////////////////////////////////
class DiscreteRVEMaterialPrecomputed : public DiscreteRVEMaterial
{
protected:
    TrsprtMaterialStatus* masterStatus;
    TrsprtMaterial* masterMaterial;
    Matrix conductivity; // precomputed form one integration point, then used everywhere
    double capacity;
public:
    DiscreteRVEMaterialPrecomputed() { name = "transport RVE precomputed material"; conductivity = Matrix(0,0); };
    virtual ~DiscreteRVEMaterialPrecomputed() {};
    virtual MaterialStatus *giveNewMaterialStatus(Element *e);
    Matrix givePrecomputedConductivity() const { return conductivity;};
    double givePrecomputedCapacity() const { return capacity;};
    TrsprtMaterialStatus* giveMasterStatus() { return masterStatus;};
    TrsprtMaterial* giveMasterMaterial() { return masterMaterial;};
    void setPrecomputedConductivityAndCapacityAndMasterMaterial(Matrix lam, double c, TrsprtMaterialStatus* masterS,  TrsprtMaterial* masterM) {conductivity=lam; capacity = c; masterStatus = masterS; masterMaterial = masterM;} ;    
};

#endif /* _MAT_RVE_H */
