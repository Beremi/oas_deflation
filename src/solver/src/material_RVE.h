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
    unsigned ndim;
    fs :: path inputfile;

    //setup for volumetric average
    PieceWiseLinearFunction *volumAverFunc;

    virtual void generateRandomFixedBC() {};
    virtual void generateVolumetricAverageBC() {};
public:
    RVEMaterialStatus(RVEMaterial *m, Element *e, unsigned ipnum, fs :: path masterfile);
    virtual ~RVEMaterialStatus();
    virtual void init();
    virtual void update();
    Model *giveWholeRVE() { return RVE; };
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
    virtual MaterialStatus *giveNewMaterialStatus(Element *e, unsigned ipnum);
    fs :: path givePathToInputFile() const { return inputfile; };
};


//////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////
// DISCRETE TRANSPORT RVE

class DiscreteTransportRVEMaterial;
class DiscreteTransportRVEMaterialStatus : public RVEMaterialStatus
{
protected:
    vector< double >aux_params; //macroscopic pressure coeff
    vector< double >orig_mater_params; //material parameters of the original model

    virtual void generateRandomFixedBC();
    virtual void generateVolumetricAverageBC();
    virtual void applyEigenStrains();
    virtual void collectStresses();
    virtual unsigned giveStrainSize(unsigned ndim) const;

public:
    DiscreteTransportRVEMaterialStatus(RVEMaterial *m, Element *e, unsigned ipnum, fs :: path masterfile);
    virtual ~DiscreteTransportRVEMaterialStatus() {};
    virtual void init();
    virtual Vector giveStress(const Vector &strain, double timeStep);//terminology from mechanics, it returns flux
    virtual Vector giveStressWithFrozenIntVars(const Vector &strain, double timeStep);
    virtual Matrix giveStiffnessTensor(string type, unsigned dimension) const;
    virtual double giveDampingConstant() const;
    virtual unsigned giveStrainSize() const { return giveStrainSize(ndim); };
};

//////////////////////////////////////////////////////////
class DiscreteTransportRVEMaterial : public RVEMaterial
{
protected:

public:
    DiscreteTransportRVEMaterial() { name = "transport RVE material"; };
    virtual ~DiscreteTransportRVEMaterial() {};
    virtual MaterialStatus *giveNewMaterialStatus(Element *e, unsigned ipnum);
};

//////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////
// DISCRETE MECHANICAL RVE

class DiscreteMechanicalRVEMaterial;
class DiscreteMechanicalRVEMaterialStatus : public DiscreteTransportRVEMaterialStatus
{
protected:
    virtual void applyEigenStrains();
    virtual void collectStresses();
    void calculateCentroid();
    virtual unsigned giveStrainSize(unsigned ndim) const;

    Point centroid;
    vector< vector< Vector > >projectors;

public:
    DiscreteMechanicalRVEMaterialStatus(RVEMaterial *m, Element *e, unsigned ipnum, fs :: path masterfile);
    virtual ~DiscreteMechanicalRVEMaterialStatus() {};
    virtual void init();
    virtual Matrix giveStiffnessTensor(string type, unsigned dimension) const;
    virtual unsigned giveStrainSize() const { return giveStrainSize(ndim); };
};

//////////////////////////////////////////////////////////
class DiscreteMechanicalRVEMaterial : public RVEMaterial
{
protected:

public:
    DiscreteMechanicalRVEMaterial() { name = "mechanical RVE material"; };
    virtual ~DiscreteMechanicalRVEMaterial() {};
    virtual MaterialStatus *giveNewMaterialStatus(Element *e, unsigned ipnum);
};

//////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////
// DISCRETE COUPLED RVE MATERIAL

class DiscreteCoupledRVEMaterial;
class DiscreteCoupledRVEMaterialStatus : public MaterialStatus
{
protected:
    fs :: path inputfileM;
    fs :: path inputfileT;
    DiscreteMechanicalRVEMaterialStatus *mechRVEstat;
    DiscreteTransportRVEMaterialStatus *trspRVEstat;

    void findFriends();
public:
    DiscreteCoupledRVEMaterialStatus(Material *m, Element *e, unsigned ipnum, fs :: path masterfileM, fs :: path masterfileT);
    virtual ~DiscreteCoupledRVEMaterialStatus();
    virtual void init();
    virtual void update();

    virtual Vector giveStress(const Vector &strain, double timeStep);
    virtual Vector giveStressWithFrozenIntVars(const Vector &strain, double timeStep);
    virtual double giveValue(string code);
    virtual Matrix giveStiffnessTensor(string type, unsigned dimension) const;
    virtual void setEigenStrain(Vector &x);
    virtual std :: string giveLineToSave() const;
};

//////////////////////////////////////////////////////////
class DiscreteCoupledRVEMaterial : public Material
{
protected:
    DiscreteMechanicalRVEMaterial *mechRVEmat;
    DiscreteTransportRVEMaterial *trspRVEmat;
public:
    DiscreteCoupledRVEMaterial() { name = "coupled RVE material"; };
    virtual ~DiscreteCoupledRVEMaterial();
    virtual MaterialStatus *giveNewMaterialStatus(Element *e, unsigned ipnum);
    virtual void readFromLine(istringstream &iss);
    DiscreteMechanicalRVEMaterial *giveMechanicalRVEmat() { return mechRVEmat; }
    DiscreteTransportRVEMaterial *giveTransportRVEmat() { return trspRVEmat; }
};

//////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////
// DISCRETE RVE MATERIAL FOR TRANSPORT PRECOMPUTED

class DiscreteTransportRVEMaterialPrecomputed;
class DiscreteTransportRVEMaterialPrecomputedStatus : public DiscreteTransportRVEMaterialStatus
{
protected:
    double temp_nonlin;
    bool is_master_status;

public:
    DiscreteTransportRVEMaterialPrecomputedStatus(RVEMaterial *m, Element *e, unsigned ipnum, fs :: path masterfile);
    virtual ~DiscreteTransportRVEMaterialPrecomputedStatus() {};
    virtual void init();
    virtual Vector giveStress(const Vector &strain, double timeStep);//terminology from mechanics, it returns flux
    virtual Vector giveStressWithFrozenIntVars(const Vector &strain, double timeStep);
    virtual Matrix giveStiffnessTensor(string type, unsigned dimension) const;
    virtual double giveDampingConstant() const;
    virtual void update();
};

//////////////////////////////////////////////////////////
class DiscreteTransportRVEMaterialPrecomputed : public DiscreteTransportRVEMaterial
{
protected:
    TrsprtMaterialStatus *masterStatus;
    TrsprtMaterial *masterMaterial;
    Matrix conductivity; // precomputed form one integration point, then used everywhere
    double capacity;
public:
    DiscreteTransportRVEMaterialPrecomputed() { name = "transport RVE precomputed material"; conductivity = Matrix(0, 0); };
    virtual ~DiscreteTransportRVEMaterialPrecomputed() {};
    virtual MaterialStatus *giveNewMaterialStatus(Element *e, unsigned ipnum);
    Matrix givePrecomputedConductivity() const { return conductivity; };
    double givePrecomputedCapacity() const { return capacity; };
    TrsprtMaterialStatus *giveMasterStatus() { return masterStatus; };
    TrsprtMaterial *giveMasterMaterial() { return masterMaterial; };
    void setPrecomputedConductivityAndCapacityAndMasterMaterial(Matrix lam, double c, TrsprtMaterialStatus *masterS,  TrsprtMaterial *masterM);
};




/*
 * //////////////////////////////////////////////////////////
 * //////////////////////////////////////////////////////////
 * // DISCRETE RVE MATERIAL FOR BOTH TRANSPORT AND MECHANICAL HOMOGENIZATION
 *
 * class DiscreteRVEMaterial;
 * class DiscreteRVEMaterialStatus : public RVEMaterialStatus
 * {
 * protected:
 *  virtual void generateRandomFixedBC();
 *  virtual void generateVolumetricAverageBC();
 *  void calculateCentroid();
 *
 *  Point centroid;
 *  bool active_mechanics, active_transport;
 *  vector< vector< Vector > >projectors;
 *
 * public:
 *  DiscreteRVEMaterialStatus(DiscreteRVEMaterial *m, Element *e, unsigned ipnum, fs :: path masterfile);
 *  virtual ~DiscreteRVEMaterialStatus() {};
 *  virtual void init();
 *  virtual Vector giveStress(const Vector &strain, double timeStep);//terminology from mechanics, it returns flux
 *  virtual Vector giveStressWithFrozenIntVars(const Vector &strain, double timeStep);
 *  virtual Matrix giveStiffnessTensor(string type, unsigned dimension) const;
 *  virtual double giveDampingConstant() const;
 * };
 *
 * //////////////////////////////////////////////////////////
 * class DiscreteRVEMaterial : public RVEMaterial
 * {
 * protected:
 *
 * public:
 *  DiscreteRVEMaterial() { name = "transport RVE material"; };
 *  virtual ~DiscreteRVEMaterial() {};
 *  virtual MaterialStatus *giveNewMaterialStatus(Element *e, unsigned ipnum);
 * };
 */
#endif /* _MAT_RVE_H */
