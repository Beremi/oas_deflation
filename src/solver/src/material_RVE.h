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
    unsigned giveNumOfDimensions() const { return ndim; };
};

//////////////////////////////////////////////////////////
class RVEMaterial : public Material
{
protected:
    fs :: path inputfile;


public:
    RVEMaterial() { name = "generic RVE material";};
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
    vector< double >orig_mater_params; //material parameters of the original model
    double macro_pressure;
    double temp_nonlin;
    bool is_master_status;
    bool is_precomputed;

    virtual void generateRandomFixedBC();
    virtual void generateVolumetricAverageBC();
    virtual void applyEigenStrains();
    virtual void collectStresses();
    virtual unsigned giveStrainSize(unsigned rdim) const;
    virtual Vector giveStressPrecomputed(const Vector &strain, double timeStep);
    virtual Matrix giveStiffnessTensorPrecomputed(string type, unsigned dimension) const;
    virtual Matrix giveDampingTensorPrecomputed() const;

public:
    DiscreteTransportRVEMaterialStatus(RVEMaterial *m, Element *e, unsigned ipnum, fs :: path masterfile);
    virtual ~DiscreteTransportRVEMaterialStatus() {};
    virtual void init();
    virtual Vector giveStress(const Vector &strain, double timeStep);//terminology from mechanics, it returns flux
    virtual Vector giveStressWithFrozenIntVars(const Vector &strain, double timeStep);
    virtual Matrix giveStiffnessTensor(string type, unsigned dimension) const;
    virtual Matrix giveDampingTensor() const;
    virtual unsigned giveStrainSize() const { return giveStrainSize(ndim); }; 
    virtual void update();   
    void setFromPrecomputedToFullModel();
    virtual void setParameterValue(string code, double value);
    void setToPrecomputed(){ is_precomputed = true; };
    bool isPrecomputed() const {return is_precomputed;};
    void setToMasterStatus(){ is_master_status = true; };
};

//////////////////////////////////////////////////////////
class DiscreteTransportRVEMaterial : public RVEMaterial
{
protected:
    TrsprtMaterialStatus *masterStatus;
    TrsprtMaterial *masterMaterial;
    Matrix conductivity; // precomputed form one integration point, then used everywhere
    double capacity;
    unsigned ndim;

public:
    DiscreteTransportRVEMaterial() { name = "transport RVE material";  conductivity = Matrix(0, 0); };
    virtual ~DiscreteTransportRVEMaterial() {};
    virtual MaterialStatus *giveNewMaterialStatus(Element *e, unsigned ipnum);
    Matrix givePrecomputedConductivity() const { return conductivity; };
    double givePrecomputedCapacity() const { return capacity; };
    TrsprtMaterialStatus *giveMasterStatus() { return masterStatus; };
    TrsprtMaterial *giveMasterMaterial() { return masterMaterial; };
    void setPrecomputedConductivityAndCapacityAndMasterMaterial(Matrix lam, double c, TrsprtMaterialStatus *masterS,  TrsprtMaterial *masterM, unsigned dimension); 
    unsigned giveDimension(){return ndim;}; 
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
    virtual unsigned giveStrainSize(unsigned rdim) const;    
    bool checkOttosenCriterion();

    Point calculateCentroid();
    vector< vector< Vector > > calculateProjectors(const Point centroid);
    virtual Vector giveStressPrecomputed(const Vector &strain, double timeStep);

public:
    DiscreteMechanicalRVEMaterialStatus(RVEMaterial *m, Element *e, unsigned ipnum, fs :: path masterfile);
    virtual ~DiscreteMechanicalRVEMaterialStatus() {};
    virtual void init();
    virtual Vector giveStress(const Vector &strain, double timeStep);//terminology from mechanics, it returns flux
    virtual Vector giveStressWithFrozenIntVars(const Vector &strain, double timeStep);
    virtual Matrix giveStiffnessTensor(string type, unsigned dimension) const;
    virtual Matrix giveDampingTensor() const;
    virtual Matrix giveInertiaTensor() const;
    virtual unsigned giveStrainSize() const { return giveStrainSize(ndim); };
    void setFromPrecomputedToFullModel();
    virtual void setToPrecomputed(){ is_precomputed = true; };
};

//////////////////////////////////////////////////////////
class DiscreteMechanicalRVEMaterial : public RVEMaterial
{
protected:
    Matrix precompElastic, precompDamping, precompInertia; 
    Point centroid;
    vector< vector< Vector > >projectors;   
    unsigned ndim;

public:
    DiscreteMechanicalRVEMaterial() { name = "mechanical RVE material"; precompElastic=Matrix(0,0);};
    virtual ~DiscreteMechanicalRVEMaterial() {};
    virtual MaterialStatus *giveNewMaterialStatus(Element *e, unsigned ipnum);
    void setPrecomputedElasticTensor(Matrix ela){precompElastic = ela;};
    void setPrecomputedDampingTensor(Matrix dam){precompDamping = dam;};
    void setPrecomputedInertiaTensor(Matrix ine){precompInertia = ine;};
    void setCentroidProjectorsAndDimensions(Point c, vector< vector< Vector > >p, unsigned d);
    Matrix givePrecomputedElasticTensor() const { return precompElastic; };
    Matrix givePrecomputedDampingTensor() const { return precompDamping; };
    Matrix givePrecomputedInertiaTensor() const { return precompInertia; };
    Point giveCentroid(){return centroid;};
    vector< vector< Vector > > * giveProjectors(){return &projectors;};
    unsigned giveDimension(){return ndim;}; 
};
 	
//////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////
// DISCRETE COUPLED RVE MATERIAL

class DiscreteCoupledRVEMaterial;
class DiscreteCoupledRVEMaterialStatus : public RVEMaterialStatus
{
protected:
    fs :: path inputfileM;
    fs :: path inputfileT;
    DiscreteMechanicalRVEMaterialStatus *mechRVEstat;
    DiscreteTransportRVEMaterialStatus *trspRVEstat;

    bool is_master_status;
    bool is_precomputed;

    double temp_volumetricStrain, volumetricStrain, volStrainRate;
    double temp_pressure, pressure, pressureRate;

    void findFriends();
public:
    DiscreteCoupledRVEMaterialStatus(RVEMaterial *m, Element *e, unsigned ipnum, fs :: path masterfileM, fs :: path masterfileT);
    virtual ~DiscreteCoupledRVEMaterialStatus();
    virtual void init();
    virtual void update();

    virtual Vector giveStress(const Vector &strain, double timeStep);
    virtual Vector giveStressWithFrozenIntVars(const Vector &strain, double timeStep);
    virtual double giveValue(string code);
    virtual Matrix giveStiffnessTensor(string type, unsigned dimension) const;
    virtual Matrix giveDampingTensor() const;
    virtual Matrix giveInertiaTensor() const;
    virtual void setEigenStrain(Vector &x);
    virtual std :: string giveLineToSave() const;
    virtual double computeBiotEffect() const;
    virtual void setParameterValue(string code, double value);
    virtual void setFromPrecomputedToFullModel();
    virtual void setToPrecomputed(){ is_precomputed = true; mechRVEstat->setToPrecomputed(); trspRVEstat->setToPrecomputed();};
};

//////////////////////////////////////////////////////////
class DiscreteCoupledRVEMaterial : public RVEMaterial
{
protected:
    DiscreteMechanicalRVEMaterial *mechRVEmat;
    DiscreteTransportRVEMaterial *trspRVEmat;
    Matrix precompElastic, precompDamping, precompInertia; 
    double biotCoeff;   

public:
    DiscreteCoupledRVEMaterial() { name = "coupled RVE material";  precompElastic=Matrix(0,0);};
    virtual ~DiscreteCoupledRVEMaterial();
    virtual MaterialStatus *giveNewMaterialStatus(Element *e, unsigned ipnum);
    virtual void readFromLine(istringstream &iss);
    DiscreteMechanicalRVEMaterial *giveMechanicalRVEmat() { return mechRVEmat; }
    DiscreteTransportRVEMaterial *giveTransportRVEmat() { return trspRVEmat; }
    double giveBiotCoefficient() const { return biotCoeff; };
    void setPrecomputedElasticDampingAndInertiaTensors(Matrix ela, Matrix dam, Matrix ine);
    Matrix givePrecomputedElasticTensor() const { return precompElastic; };
    Matrix givePrecomputedDampingTensor() const { return precompDamping; };
    Matrix givePrecomputedInertiaTensor() const { return precompInertia; };
};

#endif /* _MAT_RVE_H */
