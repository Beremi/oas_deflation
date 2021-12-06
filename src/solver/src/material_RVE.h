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
    
    Matrix transf;   //transformation matrix
    Matrix axDirs;  //directional matrix

    //setup for volumetric average
    PieceWiseLinearFunction *volumAverFunc;

    virtual void generateRandomFixedBC() {};
    virtual void generateVolumetricAverageBC() {};
public:
    RVEMaterialStatus(RVEMaterial *m, Element *e, unsigned ipnum, fs :: path masterfile, unsigned ndim);
    virtual ~RVEMaterialStatus();
    virtual void init();
    virtual void update();
    Model *giveWholeRVE() { return RVE; };
    virtual void setReferenceSystemDirections(Matrix r) {axDirs = r;};
    Matrix giveTransformationMatrix() {return transf;};
    Matrix giveReferenceSystemDirections() {return axDirs;};
};

//////////////////////////////////////////////////////////
class RVEMaterial : public Material
{
protected:
    fs :: path inputfile;
    unsigned ndim;
    bool nonlinear;  

public:
    RVEMaterial() { name = "generic RVE material"; nonlinear = true; };
    virtual ~RVEMaterial() {};
    virtual void readFromLine(istringstream &iss);
    virtual MaterialStatus *giveNewMaterialStatus(Element *e, unsigned ipnum);
    fs :: path givePathToInputFile() const { return inputfile; };
    unsigned giveNumOfDimensions() const { return ndim; };
    void setNumOfDimensions(unsigned dim) { ndim = dim; };
    void setPathToInputFolder(string f){ inputfile = GlobPaths :: BASEDIR  / f; };
    void enforceLinearity() {nonlinear = false;};
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

    Vector local_strain, local_stress;

    virtual void generateRandomFixedBC();
    virtual void generateVolumetricAverageBC();
    virtual void applyEigenStrains();
    virtual void collectStresses();
    virtual unsigned giveStrainSize(unsigned rdim) const;
    virtual Vector giveStressPrecomputed(const Vector &strain, double timeStep);
    virtual Matrix giveStiffnessTensorLocal(string type, unsigned dimension) const;
    virtual Matrix giveStiffnessTensorPrecomputedLocal(string type, unsigned dimension) const;
    virtual Matrix giveDampingTensorPrecomputed() const;

    virtual void transformStrain();
    virtual void transformStress();
    virtual void calculateTransformationMatrix();

public:
    DiscreteTransportRVEMaterialStatus(RVEMaterial *m, Element *e, unsigned ipnum, fs :: path masterfile, unsigned ndim);
    virtual ~DiscreteTransportRVEMaterialStatus() {};
    virtual void init();
    virtual Vector giveStress(const Vector &strain, double timeStep);//terminology from mechanics, it returns flux
    virtual Vector giveStressWithFrozenIntVars(const Vector &strain, double timeStep);
    virtual Matrix giveStiffnessTensor(string type, unsigned dimension) const;
    virtual Matrix giveDampingTensor() const;
    virtual unsigned giveStrainSize() const;
    virtual void update();
    void setFromPrecomputedToFullModel();
    virtual void setParameterValue(string code, double value);
    void setToPrecomputed() { is_precomputed = true; };
    bool isPrecomputed() const { return is_precomputed; };
    void setToMasterStatus() { is_master_status = true; };
};

//////////////////////////////////////////////////////////
class DiscreteTransportRVEMaterial : public RVEMaterial
{
protected:
    TrsprtMaterialStatus *masterStatus;
    TrsprtMaterial *masterMaterial;
    Matrix conductivity; // precomputed form one integration point, then used everywhere
    double capacity;

public:
    DiscreteTransportRVEMaterial() { name = "transport RVE material";  conductivity = Matrix(0, 0); };
    virtual ~DiscreteTransportRVEMaterial() {};
    virtual MaterialStatus *giveNewMaterialStatus(Element *e, unsigned ipnum);
    Matrix givePrecomputedConductivity() const { return conductivity; };
    double givePrecomputedCapacity() const { return capacity; };
    TrsprtMaterialStatus *giveMasterStatus() { return masterStatus; };
    TrsprtMaterial *giveMasterMaterial() { return masterMaterial; };
    void setPrecomputedConductivityAndCapacityAndMasterMaterial(Matrix lam, double c, TrsprtMaterialStatus *masterS,  TrsprtMaterial *masterM);
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
    vector< vector< Vector > >calculateProjectors(const Point centroid);
    virtual Vector giveStressPrecomputed(const Vector &strain, double timeStep);
    virtual Matrix giveStiffnessTensorLocal(string type, unsigned dimension) const;

    virtual void transformStrain();
    virtual void transformStress();
    virtual void calculateTransformationMatrix();

public:
    DiscreteMechanicalRVEMaterialStatus(RVEMaterial *m, Element *e, unsigned ipnum, fs :: path masterfile, unsigned ndim);
    virtual ~DiscreteMechanicalRVEMaterialStatus() {};
    virtual void init();
    virtual Vector giveStress(const Vector &strain, double timeStep);//terminology from mechanics, it returns flux
    virtual Vector giveStressWithFrozenIntVars(const Vector &strain, double timeStep);
    virtual Matrix giveStiffnessTensor(string type, unsigned dimension) const;
    virtual Matrix giveDampingTensor() const;
    virtual Matrix giveInertiaTensor() const;
    virtual unsigned giveStrainSize() const;
    virtual double giveCrackVolume() const;
    void setFromPrecomputedToFullModel();
    virtual void setToPrecomputed() { is_precomputed = true; };
};

//////////////////////////////////////////////////////////
class DiscreteMechanicalRVEMaterial : public RVEMaterial
{
protected: 
    Matrix precompElastic, precompDamping, precompInertia;
    Point centroid;
    vector< vector< Vector > >projectors;

public:
    DiscreteMechanicalRVEMaterial() { name = "mechanical RVE material"; precompElastic = Matrix(0, 0);};
    virtual ~DiscreteMechanicalRVEMaterial() {};
    virtual MaterialStatus *giveNewMaterialStatus(Element *e, unsigned ipnum);
    void setPrecomputedElasticTensor(Matrix ela) { precompElastic = ela; };
    void setPrecomputedDampingTensor(Matrix dam) { precompDamping = dam; };
    void setPrecomputedInertiaTensor(Matrix ine) { precompInertia = ine; };
    void setCentroidAndProjectors(Point c, vector< vector< Vector > >p);
    Matrix givePrecomputedElasticTensor() const { return precompElastic; };
    Matrix givePrecomputedDampingTensor() const { return precompDamping; };
    Matrix givePrecomputedInertiaTensor() const { return precompInertia; };
    Point giveCentroid() { return centroid; };
    vector< vector< Vector > > *giveProjectors() { return & projectors; };
    bool isNonlinear() const { return nonlinear; }; 
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
    double temp_crackVolume, crackVolume, crackVolumeRate;

    void findFriends();
    void updateRateVariables(double timeStep);
public:
    DiscreteCoupledRVEMaterialStatus(RVEMaterial *m, Element *e, unsigned ipnum, fs :: path masterfileM, fs :: path masterfileT, unsigned ndim);
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
    virtual Vector giveInternalSource() const;
    virtual void setParameterValue(string code, double value);
    virtual void setFromPrecomputedToFullModel();
    virtual void setToPrecomputed() { is_precomputed = true; mechRVEstat->setToPrecomputed(); trspRVEstat->setToPrecomputed(); };
    virtual void setReferenceSystemDirections(Matrix r);
};

//////////////////////////////////////////////////////////
class DiscreteCoupledRVEMaterial : public RVEMaterial
{
protected:
    DiscreteMechanicalRVEMaterial *mechRVEmat;
    DiscreteTransportRVEMaterial *trspRVEmat;
    Matrix precompElastic, precompDamping, precompInertia;
    double biotCoeff, PUCVolume;
    DiscreteTrsprtCoupledMaterial *masterMaterial;

public:
    DiscreteCoupledRVEMaterial() { name = "coupled RVE material";  precompElastic = Matrix(0, 0); produceInternalSources = true; PUCVolume = 0; };
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
    DiscreteTrsprtCoupledMaterial *giveMasterMaterial() { return masterMaterial; };
    void setMasterMaterial(DiscreteTrsprtCoupledMaterial *masterM) { masterMaterial = masterM; };
    void setPUCVolume(double vol) { PUCVolume = vol; };
    double givePUCVolume() const { return PUCVolume; };
};

#endif /* _MAT_RVE_H */
