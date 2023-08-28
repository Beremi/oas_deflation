#ifndef _MATERIAL_RVE_H
#define _MATERIAL_RVE_H

#include "material_vectorial.h"
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


    bool is_precomputed;

    virtual void generateRandomFixedBC() {};
    virtual void generateVolumetricAverageBC() {};
public:
    RVEMaterialStatus(RVEMaterial *m, Element *e, unsigned ipnum, fs :: path masterfile, unsigned ndim);
    virtual ~RVEMaterialStatus();
    virtual void init();
    virtual void update();
    Model *giveWholeRVE() { return RVE; };
    virtual void setReferenceSystemDirections(Matrix r) { axDirs = r; };
    Matrix giveTransformationMatrix() { return transf; };
    Matrix giveReferenceSystemDirections() { return axDirs; };
};

//////////////////////////////////////////////////////////
class RVEMaterial : public Material
{
protected:
    fs :: path inputfile;
    bool nonlinear;

    bool elastic_sol_is_Voigt;  //distinguish whether the solution in initial precomputed state is really solved elastically or using Voigt constraint
    bool start_from_precomputed;

public:
    RVEMaterial(unsigned dimension) : Material(dimension)  { name = "generic RVE material"; nonlinear = true; elastic_sol_is_Voigt = false; start_from_precomputed = true; };
    virtual ~RVEMaterial() {};
    virtual void readFromLine(std :: istringstream &iss);
    virtual MaterialStatus *giveNewMaterialStatus(Element *e, unsigned ipnum);
    fs :: path givePathToInputFile() const { return inputfile; };
    void setPathToInputFolder(std :: string f) { inputfile = GlobPaths :: BASEDIR  / f; };
    void enforceLinearity() { nonlinear = false; };
    bool isElasticSolutionVoigt() const { return elastic_sol_is_Voigt; };
    bool shouldStartFromPrecomputed() const { return start_from_precomputed; };
    void setStartFromPrecomputed(bool s) { start_from_precomputed = s; };
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

    Vector local_strain, local_stress;

    virtual void generateRandomFixedBC();
    virtual void generateVolumetricAverageBC();
    virtual void applyEigenStrains();
    virtual void collectStresses();
    virtual Vector giveStressPrecomputed(const Vector &strain, double timeStep);
    virtual Matrix giveStiffnessTensorLocal(std :: string type) const;
    Matrix giveStiffnessTensorLocalExact(std :: string type);
    virtual Matrix giveStiffnessTensorPrecomputedLocal(std :: string type) const;
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
    virtual Matrix giveStiffnessTensor(std :: string type) const;
    virtual Matrix giveDampingTensor() const;
    virtual void update();
    void setFromPrecomputedToFullModel();
    virtual void setParameterValue(std :: string code, double value);
    void setToPrecomputed() { is_precomputed = true; };
    bool isPrecomputed() const { return is_precomputed; };
    void setToMasterStatus() { is_master_status = true; };
    virtual bool giveValues(std :: string code, Vector &result) const;
};

//////////////////////////////////////////////////////////
class DiscreteTransportRVEMaterial : public RVEMaterial
{
protected:
    VectTrsprtMaterialStatus *masterStatus;
    VectTrsprtMaterial *masterMaterial;
    Matrix conductivity; // precomputed form one integration point, then used everywhere
    double capacity;

public:
    DiscreteTransportRVEMaterial(unsigned dimension);
    virtual ~DiscreteTransportRVEMaterial() {};
    virtual MaterialStatus *giveNewMaterialStatus(Element *e, unsigned ipnum);
    Matrix givePrecomputedConductivity() const { return conductivity; };
    double givePrecomputedCapacity() const { return capacity; };
    VectTrsprtMaterialStatus *giveMasterStatus() { return masterStatus; };
    VectTrsprtMaterial *giveMasterMaterial() { return masterMaterial; };
    void setPrecomputedConductivity(Matrix lam);
    void setMasterMaterial(VectTrsprtMaterialStatus *masterS,  VectTrsprtMaterial *masterM);
    void setPrecomputedCapacity(double c);
};

//////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////
// DISCRETE MECHANICAL COSSERAT RVE

class DiscreteMechanicalRVEMaterial;
class DiscreteMechanicalRVEMaterialStatus : public DiscreteTransportRVEMaterialStatus
{
protected:
    virtual void applyEigenStrains();
    virtual void collectStresses();
    //virtual unsigned giveStrainSize(unsigned rdim) const;
    bool checkOttosenCriterion();

    Point calculateCentroid();
    std :: vector< std :: vector< Vector > >calculateProjectors(const Point centroid);
    virtual Vector giveStressPrecomputed(const Vector &strain, double timeStep);
    virtual Matrix giveStiffnessTensorLocal(std :: string type) const;

    virtual void transformStrain();
    virtual void transformStress();
    virtual void calculateTransformationMatrix();

public:
    DiscreteMechanicalRVEMaterialStatus(RVEMaterial *m, Element *e, unsigned ipnum, fs :: path masterfile, unsigned ndim);
    virtual ~DiscreteMechanicalRVEMaterialStatus() {};
    virtual void init();
    virtual Vector giveStress(const Vector &strain, double timeStep);//terminology from mechanics, it returns flux
    virtual Vector giveStressWithFrozenIntVars(const Vector &strain, double timeStep);
    virtual Matrix giveStiffnessTensor(std :: string type) const;
    virtual Matrix giveDampingTensor() const;
    virtual Matrix giveInertiaTensor() const;
    //virtual unsigned giveStrainSize() const;
    virtual double giveCrackVolume() const;
    void setFromPrecomputedToFullModel();
    virtual void setToPrecomputed() { is_precomputed = true; };
    virtual bool giveValues(std :: string code, Vector &result) const;
};

//////////////////////////////////////////////////////////
class DiscreteMechanicalRVEMaterial : public DiscreteTransportRVEMaterial
{
protected:
    Matrix precompElastic, precompDamping, precompInertia;
    Point centroid;
    std :: vector< std :: vector< Vector > >projectors;
    bool project_curvature;

public:
    DiscreteMechanicalRVEMaterial(unsigned dimension);
    virtual ~DiscreteMechanicalRVEMaterial() {};
    virtual MaterialStatus *giveNewMaterialStatus(Element *e, unsigned ipnum);
    void setPrecomputedElasticTensor(Matrix ela) { precompElastic = ela; };
    void setPrecomputedDampingTensor(Matrix dam) { precompDamping = dam; };
    void setPrecomputedInertiaTensor(Matrix ine) { precompInertia = ine; };
    void setCentroidAndProjectors(Point c, std :: vector< std :: vector< Vector > >p);
    Matrix givePrecomputedElasticTensor() const { return precompElastic; };
    Matrix givePrecomputedDampingTensor() const { return precompDamping; };
    Matrix givePrecomputedInertiaTensor() const { return precompInertia; };
    Point giveCentroid() { return centroid; };
    std :: vector< std :: vector< Vector > > *giveProjectors() { return & projectors; };
    bool isNonlinear() const { return nonlinear; };
    bool projectCurvature() const {return project_curvature;};
    virtual void readFromLine(std :: istringstream &iss);
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
    virtual bool giveValues(std :: string code, Vector &result) const;
    virtual Matrix giveStiffnessTensor(std :: string type) const;
    virtual Matrix giveDampingTensor() const;
    virtual Matrix giveInertiaTensor() const;
    virtual void setEigenStrain(Vector &x);
    virtual std :: string giveLineToSave() const;
    virtual Vector giveInternalSource() const;
    virtual void setParameterValue(std :: string code, double value);
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
    VectTrsprtCoupledMaterial *masterMaterial;

public:
    DiscreteCoupledRVEMaterial(unsigned dimension);
    virtual ~DiscreteCoupledRVEMaterial();
    virtual MaterialStatus *giveNewMaterialStatus(Element *e, unsigned ipnum);
    virtual void readFromLine(std :: istringstream &iss);
    DiscreteMechanicalRVEMaterial *giveMechanicalRVEmat() { return mechRVEmat; }
    DiscreteTransportRVEMaterial *giveTransportRVEmat() { return trspRVEmat; }
    double giveBiotCoefficient() const { return biotCoeff; };
    void setPrecomputedDampingAndInertiaTensors(Matrix dam, Matrix ine);
    void setPrecomputedElasticTensor(Matrix ela);
    Matrix givePrecomputedElasticTensor() const { return precompElastic; };
    Matrix givePrecomputedDampingTensor() const { return precompDamping; };
    Matrix givePrecomputedInertiaTensor() const { return precompInertia; };
    VectTrsprtCoupledMaterial *giveMasterMaterial() { return masterMaterial; };
    void setMasterMaterial(VectTrsprtCoupledMaterial *masterM) { masterMaterial = masterM; };
    void setPUCVolume(double vol) { PUCVolume = vol; };
    double givePUCVolume() const { return PUCVolume; };
};

#endif /* _MATERIAL_RVE_H */
