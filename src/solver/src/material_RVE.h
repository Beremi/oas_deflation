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

    MyMatrix transf;   //transformation matrix
    MyMatrix axDirs;  //directional matrix

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
    virtual void setReferenceSystemDirections(MyMatrix r) { axDirs = r; };
    MyMatrix giveTransformationMatrix() { return transf; };
    MyMatrix giveReferenceSystemDirections() { return axDirs; };
};

//////////////////////////////////////////////////////////
class RVEMaterial : public Material
{
protected:
    fs :: path inputfile;
    unsigned ndim;
    bool nonlinear;

    bool elastic_sol_is_Voigt;  //distinguish whether the solution in initial precomputed state is really solved elastically or using Voigt constraint

public:
    RVEMaterial() { name = "generic RVE material"; nonlinear = true; elastic_sol_is_Voigt = false; };
    virtual ~RVEMaterial() {};
    virtual void readFromLine(std :: istringstream &iss);
    virtual MaterialStatus *giveNewMaterialStatus(Element *e, unsigned ipnum);
    fs :: path givePathToInputFile() const { return inputfile; };
    unsigned giveNumOfDimensions() const { return ndim; };
    void setNumOfDimensions(unsigned dim) { ndim = dim; };
    void setPathToInputFolder(std :: string f) { inputfile = GlobPaths :: BASEDIR  / f; };
    void enforceLinearity() { nonlinear = false; };
    bool isElasticSolutionVoigt() const { return elastic_sol_is_Voigt; };
};


//////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////
// DISCRETE TRANSPORT RVE

class DiscreteTransportRVEMaterial;
class DiscreteTransportRVEMaterialStatus : public RVEMaterialStatus
{
protected:
    std :: vector< double >orig_mater_params; //material parameters of the original model
    double macro_pressure;
    double temp_nonlin;
    bool is_master_status;
    bool is_precomputed;

    MyVector local_strain, local_stress;

    virtual void generateRandomFixedBC();
    virtual void generateVolumetricAverageBC();
    virtual void applyEigenStrains();
    virtual void collectStresses();
    virtual unsigned giveStrainSize(unsigned rdim) const;
    virtual MyVector giveStressPrecomputed(const MyVector &strain, double timeStep);
    virtual MyMatrix giveStiffnessTensorLocal(std :: string type, unsigned dimension) const;
    virtual MyMatrix giveStiffnessTensorPrecomputedLocal(std :: string type, unsigned dimension) const;
    virtual MyMatrix giveDampingTensorPrecomputed() const;

    virtual void transformStrain();
    virtual void transformStress();
    virtual void calculateTransformationMatrix();

public:
    DiscreteTransportRVEMaterialStatus(RVEMaterial *m, Element *e, unsigned ipnum, fs :: path masterfile, unsigned ndim);
    virtual ~DiscreteTransportRVEMaterialStatus() {};
    virtual void init();
    virtual MyVector giveStress(const MyVector &strain, double timeStep);//terminology from mechanics, it returns flux
    virtual MyVector giveStressWithFrozenIntVars(const MyVector &strain, double timeStep);
    virtual MyMatrix giveStiffnessTensor(std :: string type, unsigned dimension) const;
    virtual MyMatrix giveDampingTensor() const;
    virtual unsigned giveStrainSize() const;
    virtual void update();
    void setFromPrecomputedToFullModel();
    virtual void setParameterValue(std :: string code, double value);
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
    MyMatrix conductivity; // precomputed form one integration point, then used everywhere
    double capacity;

public:
    DiscreteTransportRVEMaterial() { name = "transport RVE material";  conductivity = MyMatrix(0, 0); };
    virtual ~DiscreteTransportRVEMaterial() {};
    virtual MaterialStatus *giveNewMaterialStatus(Element *e, unsigned ipnum);
    MyMatrix givePrecomputedConductivity() const { return conductivity; };
    double givePrecomputedCapacity() const { return capacity; };
    TrsprtMaterialStatus *giveMasterStatus() { return masterStatus; };
    TrsprtMaterial *giveMasterMaterial() { return masterMaterial; };
    void setPrecomputedConductivityAndCapacityAndMasterMaterial(MyMatrix lam, double c, TrsprtMaterialStatus *masterS,  TrsprtMaterial *masterM);
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
    std :: vector< std :: vector< MyVector > >calculateProjectors(const Point centroid);
    virtual MyVector giveStressPrecomputed(const MyVector &strain, double timeStep);
    virtual MyMatrix giveStiffnessTensorLocal(std :: string type, unsigned dimension) const;

    virtual void transformStrain();
    virtual void transformStress();
    virtual void calculateTransformationMatrix();

public:
    DiscreteMechanicalRVEMaterialStatus(RVEMaterial *m, Element *e, unsigned ipnum, fs :: path masterfile, unsigned ndim);
    virtual ~DiscreteMechanicalRVEMaterialStatus() {};
    virtual void init();
    virtual MyVector giveStress(const MyVector &strain, double timeStep);//terminology from mechanics, it returns flux
    virtual MyVector giveStressWithFrozenIntVars(const MyVector &strain, double timeStep);
    virtual MyMatrix giveStiffnessTensor(std :: string type, unsigned dimension) const;
    virtual MyMatrix giveDampingTensor() const;
    virtual MyMatrix giveInertiaTensor() const;
    virtual unsigned giveStrainSize() const;
    virtual double giveCrackVolume() const;
    void setFromPrecomputedToFullModel();
    virtual void setToPrecomputed() { is_precomputed = true; };
};

//////////////////////////////////////////////////////////
class DiscreteMechanicalRVEMaterial : public RVEMaterial
{
protected:
    MyMatrix precompElastic, precompDamping, precompInertia;
    Point centroid;
    std :: vector< std :: vector< MyVector > >projectors;

public:
    DiscreteMechanicalRVEMaterial() { name = "mechanical RVE material"; precompElastic = MyMatrix(0, 0); };
    virtual ~DiscreteMechanicalRVEMaterial() {};
    virtual MaterialStatus *giveNewMaterialStatus(Element *e, unsigned ipnum);
    void setPrecomputedElasticTensor(MyMatrix ela) { precompElastic = ela; };
    void setPrecomputedDampingTensor(MyMatrix dam) { precompDamping = dam; };
    void setPrecomputedInertiaTensor(MyMatrix ine) { precompInertia = ine; };
    void setCentroidAndProjectors(Point c, std :: vector< std :: vector< MyVector > >p);
    MyMatrix givePrecomputedElasticTensor() const { return precompElastic; };
    MyMatrix givePrecomputedDampingTensor() const { return precompDamping; };
    MyMatrix givePrecomputedInertiaTensor() const { return precompInertia; };
    Point giveCentroid() { return centroid; };
    std :: vector< std :: vector< MyVector > > *giveProjectors() { return & projectors; };
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

    virtual MyVector giveStress(const MyVector &strain, double timeStep);
    virtual MyVector giveStressWithFrozenIntVars(const MyVector &strain, double timeStep);
    virtual void giveValues(std :: string code, MyVector & result) const;
    virtual MyMatrix giveStiffnessTensor(std :: string type, unsigned dimension) const;
    virtual MyMatrix giveDampingTensor() const;
    virtual MyMatrix giveInertiaTensor() const;
    virtual void setEigenStrain(MyVector &x);
    virtual std :: string giveLineToSave() const;
    virtual MyVector giveInternalSource() const;
    virtual void setParameterValue(std :: string code, double value);
    virtual void setFromPrecomputedToFullModel();
    virtual void setToPrecomputed() { is_precomputed = true; mechRVEstat->setToPrecomputed(); trspRVEstat->setToPrecomputed(); };
    virtual void setReferenceSystemDirections(MyMatrix r);
};

//////////////////////////////////////////////////////////
class DiscreteCoupledRVEMaterial : public RVEMaterial
{
protected:
    DiscreteMechanicalRVEMaterial *mechRVEmat;
    DiscreteTransportRVEMaterial *trspRVEmat;
    MyMatrix precompElastic, precompDamping, precompInertia;
    double biotCoeff, PUCVolume;
    DiscreteTrsprtCoupledMaterial *masterMaterial;

public:
    DiscreteCoupledRVEMaterial() { name = "coupled RVE material";  precompElastic = MyMatrix(0, 0); produceInternalSources = true; PUCVolume = 0; };
    virtual ~DiscreteCoupledRVEMaterial();
    virtual MaterialStatus *giveNewMaterialStatus(Element *e, unsigned ipnum);
    virtual void readFromLine(std :: istringstream &iss);
    DiscreteMechanicalRVEMaterial *giveMechanicalRVEmat() { return mechRVEmat; }
    DiscreteTransportRVEMaterial *giveTransportRVEmat() { return trspRVEmat; }
    double giveBiotCoefficient() const { return biotCoeff; };
    void setPrecomputedElasticDampingAndInertiaTensors(MyMatrix ela, MyMatrix dam, MyMatrix ine);
    MyMatrix givePrecomputedElasticTensor() const { return precompElastic; };
    MyMatrix givePrecomputedDampingTensor() const { return precompDamping; };
    MyMatrix givePrecomputedInertiaTensor() const { return precompInertia; };
    DiscreteTrsprtCoupledMaterial *giveMasterMaterial() { return masterMaterial; };
    void setMasterMaterial(DiscreteTrsprtCoupledMaterial *masterM) { masterMaterial = masterM; };
    void setPUCVolume(double vol) { PUCVolume = vol; };
    double givePUCVolume() const { return PUCVolume; };
};

#endif /* _MAT_RVE_H */
