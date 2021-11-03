#ifndef _MATERIAL_H
#define _MATERIAL_H

#include "linear_algebra.h"
#include <vector>
#include <iostream>
#include <fstream>

//////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////
//BASIC MATERIAL

class Element; //forward declaration

class Material;
class MaterialStatus
{
private:

public:
    MaterialStatus(Material *m, Element *e, unsigned ipnum) { name = "basic mat. status"; mat = m; element = e; idx = ipnum; };
    MaterialStatus(Material *m) { name = "basic mat. status"; mat = m; };
    virtual ~MaterialStatus() {};
    string whoAmI() { return name; }
    string giveName() { return name; }
    Material *giveMaterial() { return mat; };
    virtual void init() {};
    virtual void update();
    virtual Vector giveStress(const Vector &strain, double timeStep) { return giveStressWithFrozenIntVars(strain, timeStep); };
    virtual Vector giveStressWithFrozenIntVars(const Vector &strain, double timeStep) { ( void ) strain; ( void ) timeStep; return Vector(0); };
    virtual double giveValue(string code) const { ( void ) code; return 0; };
    virtual Vector giveTempStress() const { return temp_stress; };
    virtual Vector giveUpdatedStress() const { return updt_stress; };
    virtual Vector giveTempStrain() const { return temp_strain; };
    virtual Vector giveUpdatedStrain() const { return updt_strain; };
    virtual Matrix giveStiffnessTensor(string type, unsigned dimension) const { ( void ) dimension; ( void ) type; return Matrix(0, 0); };
    virtual Matrix giveMassTensor() const { return Matrix(0, 0); };
    virtual Matrix giveDampingTensor() const { return Matrix(0, 0); };
    virtual void setEigenStrain(Vector &x);
    //virtual void setID(unsigned i) { idx = i; };
    virtual std :: string giveLineToSave() const { return "no internal variables to export, you need to implement this possibility for " + this->name; }
    virtual void readFromLine(istringstream &iss) {
        ( void ) iss;
        std :: cout << "no internal variables to read, you need to implement this possibility for " << this->name << '\n';
    };
    virtual Vector giveInternalSource()const { return Vector(0); };
    virtual bool isElastic(const bool &now = false) const;
    virtual void setParameterValue(string code, double value) { ( void ) code; ( void ) value; };

protected:
    Vector addEigenStrain(const Vector &totalStrain) const;
    Element *element;
    string name;
    Material *mat;
    Vector eigenstrain;
    Vector updt_strain, temp_strain, updt_stress, temp_stress;
    unsigned idx;
};


//////////////////////////////////////////////////////////
class Material
{
private:
    unsigned id;  // to be able to save element
protected:
    bool produceInternalSources;
public:
    Material() { name = "basic material"; produceInternalSources = false; };
    virtual ~Material() {};
    virtual void readFromLine(istringstream &iss) { ( void ) iss; };
    virtual MaterialStatus *giveNewMaterialStatus(Element *e, unsigned ipnum) { MaterialStatus *newStatus = new MaterialStatus(this, e, ipnum); return newStatus; };
    string whoAmI() { return name; }
    string giveName() { return name; }
    unsigned giveId() { return id; }
    void setId(const unsigned &i) { this->id = i; }
    virtual void init() { };
    bool isProducingInternalSources()const { return produceInternalSources; }
protected:
    string name;
};

//////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////
// TRANSPORT MATERIAL

class TrsprtMaterial;
class TrsprtMaterialStatus : public MaterialStatus
{
protected:
    double effConductivity, temp_effConductivity;
    double avgPressure;
public:
    TrsprtMaterialStatus(TrsprtMaterial *m, Element *e, unsigned ipnum);
    virtual ~TrsprtMaterialStatus() {};
    virtual Vector giveStress(const Vector &strain, double timeStep); //terminology from mechanics, it returns flux
    virtual Vector giveStressWithFrozenIntVars(const Vector &strain, double timeStep);
    virtual Matrix giveStiffnessTensor(string type, unsigned dimension) const;
    virtual Matrix giveDampingTensor() const;
    virtual double giveValue(string code) const;
    virtual double giveEffectiveConductivity(string type) const;
    virtual double updateEffectiveConductivity() const;
    virtual double calculatePressureDependentPermeability(double pressure) const;
    virtual bool isElastic(const bool &now = false) const;
    virtual void setParameterValue(string code, double value);
};

//////////////////////////////////////////////////////////
class TrsprtMaterial : public Material
{
protected:
    double permeability, viscosity, capacity, density, a, m;
public:
    TrsprtMaterial() { name = "transport material"; a = -1.; m = 0; };
    ~TrsprtMaterial() {};
    double giveCapacity() const { return capacity; };
    double giveDensity() const { return density; };
    double givePermeability() const { return permeability; };
    double giveViscosity() const { return viscosity; };
    double giveParamA() const { return a; };
    double giveParamM() const { return m; };
    void setPermeability(double new_p) { permeability = new_p; };
    void setParamA(double new_a) { a = new_a; };
    void readFromLine(istringstream &iss);
    MaterialStatus *giveNewMaterialStatus(Element *e, unsigned ipnum);
};

//////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////
// DISCRETE TRANSPORT MATERIAL

class DiscreteTrsprtMaterial;
class DiscreteTrsprtMaterialStatus : public TrsprtMaterialStatus
{
protected:
public:
    DiscreteTrsprtMaterialStatus(TrsprtMaterial *m, Element *e, unsigned ipnum);
    virtual ~DiscreteTrsprtMaterialStatus() {};
    virtual Matrix giveStiffnessTensor(string type, unsigned dimension) const;
};

//////////////////////////////////////////////////////////
class DiscreteTrsprtMaterial : public TrsprtMaterial
{
protected:

public:
    DiscreteTrsprtMaterial() { name = "discrete transport material"; };
    ~DiscreteTrsprtMaterial() {};
    MaterialStatus *giveNewMaterialStatus(Element *e, unsigned ipnum);
};

//////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////
// DISCRETE TRANSPORT MATERIAL WITH COUPLING TERMS

class DiscreteTrsprtCoupledMaterial;
class DiscreteTrsprtCoupledMaterialStatus : public DiscreteTrsprtMaterialStatus
{
protected:
    double tempVolumetricStrain, volumetricStrain, volStrainRate, crackParam, tempCrackVolume, crackVolume, crackVolumeRate, pressure, pressureRate;
public:
    DiscreteTrsprtCoupledMaterialStatus(TrsprtMaterial *m, Element *e, unsigned ipnum);
    virtual ~DiscreteTrsprtCoupledMaterialStatus() {};
    virtual Vector giveStress(const Vector &strain, double timeStep);
    virtual Vector giveStressWithFrozenIntVars(const Vector &strain, double timeStep);
    virtual double giveEffectiveConductivity(string type) const;
    virtual void update();
    virtual void init() { volumetricStrain = 0; };
    virtual Vector giveInternalSource() const;
    virtual double giveValue(string code) const;
    virtual void setParameterValue(string code, double value);
    virtual double updateEffectiveConductivity() const;
};

//////////////////////////////////////////////////////////
class DiscreteTrsprtCoupledMaterial : public DiscreteTrsprtMaterial
{
private:
    double crack_turtuosity, biotCoeff, refP, Kw; 
public:
    DiscreteTrsprtCoupledMaterial() { name = "coupled transport material";  produceInternalSources = true; refP = 0.; Kw = 2.15e9;};
    ~DiscreteTrsprtCoupledMaterial() {};
    void readFromLine(istringstream &iss);
    MaterialStatus *giveNewMaterialStatus(Element *e, unsigned ipnum);
    double giveTurtuosity() { return crack_turtuosity; };
    double giveBiotCoeff() const { return biotCoeff; };
    double giveKw() const {return Kw;};
    double giveReferencePressure() const {return refP; };
};

//////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////
// TENSORIAL MECHANICAL MATERIAL

class ElasticMechMaterial;
class ElasticMechMaterialStatus : public MaterialStatus
{
protected:

public:
    ElasticMechMaterialStatus(ElasticMechMaterial *m, Element *e, unsigned ipnum);
    virtual ~ElasticMechMaterialStatus() {};
    virtual Matrix giveStiffnessTensor(string type, unsigned dim) const;
    virtual Vector giveStress(const Vector &strain, double timeStep);
    virtual Vector giveStressWithFrozenIntVars(const Vector &strain, double timeStep);
    virtual double giveMassConstant() const;
};

//////////////////////////////////////////////////////////
class ElasticMechMaterial : public Material
{
protected:
    double E, nu, density;
    bool planeStress;

public:
    ElasticMechMaterial() { name = "elastic tensorial mechanical material"; planeStress = true; };
    ~ElasticMechMaterial() {};
    virtual void readFromLine(istringstream &iss);
    virtual MaterialStatus *giveNewMaterialStatus(Element *e, unsigned ipnum);
    double giveElasticModulus() const { return E; }
    double givePoissonsRatio() const { return nu; }
    double giveDensity() const { return density; };
    bool isPlaneStress() { return planeStress; };
};

//////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////
// COSSERAT MECHANICAL MATERIAL

class CosseratMechMaterial;
class CosseratMechMaterialStatus : public ElasticMechMaterialStatus
{
protected:

public:
    CosseratMechMaterialStatus(CosseratMechMaterial *m, Element *e, unsigned ipnum);
    virtual ~CosseratMechMaterialStatus() {};
    virtual Matrix giveStiffnessTensor(string type, unsigned dim) const;
    virtual Vector giveStress(const Vector &strain, double timeStep);
    virtual Vector giveStressWithFrozenIntVars(const Vector &strain, double timeStep);
};

//////////////////////////////////////////////////////////
class CosseratMechMaterial : public ElasticMechMaterial
{
protected:
    double lc, muc;

public:
    CosseratMechMaterial() { name = "elastic Cosserat mechanical material"; planeStress = true; };
    ~CosseratMechMaterial() {};
    virtual void readFromLine(istringstream &iss);
    virtual MaterialStatus *giveNewMaterialStatus(Element *e, unsigned ipnum);
    double giveCharacteristicLength() const { return lc; }
    double giveCosseratShearParam() const { return muc; }
};

//////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////
// DISCRETE MECHANICAL LINEAR MATERIAL

class DisMechMaterial;
class DisMechMaterialStatus : public MaterialStatus
{
protected:

public:
    DisMechMaterialStatus(DisMechMaterial *m, Element *e, unsigned ipnum);
    virtual ~DisMechMaterialStatus() {};
    virtual Matrix giveStiffnessTensor(string type, unsigned dim) const;
    virtual Vector giveStress(const Vector &strain, double timeStep);
    virtual Vector giveStressWithFrozenIntVars(const Vector &strain, double timeStep);
    double giveDensity() const;
    virtual bool isElastic(const bool &now = false) const { ( void ) now; return true; };
    virtual double giveValue(string code) const;
};

//////////////////////////////////////////////////////////
class DisMechMaterial : public Material
{
protected:
    double E0, alpha, density;

public:
    DisMechMaterial() { name = "discrete mechanical material"; };
    ~DisMechMaterial() {};
    double giveDensity() { return density; };
    virtual void readFromLine(istringstream &iss);
    virtual MaterialStatus *giveNewMaterialStatus(Element *e, unsigned ipnum);
    double giveAlpha() const { return alpha; }
    double giveE0() const { return E0; }
};

#endif /* _MATERIAL_H */
