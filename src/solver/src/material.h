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
    MaterialStatus(Material *m, Element *e) { name = "basic mat. status"; mat = m; element = e; };
    MaterialStatus(Material *m) { name = "basic mat. status"; mat = m; };
    virtual ~MaterialStatus() {};
    string whoAmI() { return name; }
    Material *giveMaterial() { return mat; };
    virtual void init() {};
    virtual void update();
    virtual Vector giveStress(const Vector &strain) { return giveStressWithFrozenIntVars(strain); };
    virtual Vector giveStressWithFrozenIntVars(const Vector &strain) const { ( void ) strain; return Vector(0); };
    virtual double giveValue(string code) const { ( void ) code; return 0; };
    virtual Vector giveTempStress() const { return temp_stress; };
    virtual Vector giveUpdatedStress() const { return updt_stress; };
    virtual Vector giveTempStrain() const { return temp_strain; };
    virtual Vector giveUpdatedStrain() const { return updt_strain; };
    virtual Matrix giveStiffnessTensor(string type, unsigned dimension) const { ( void ) dimension; return Matrix(0, 0); };
    virtual double giveMassConstant() const { return 0; };
    virtual void setEigenStrain(Vector &x) { eigenstrain = x; };
    void setID(unsigned i) { idx = i; };
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

public:
    Material() { name = "basic material"; };
    virtual ~Material() {};
    virtual void readFromLine(istringstream &iss) { ( void ) iss; };
    virtual MaterialStatus *giveNewMaterialStatus(Element *e) { MaterialStatus *newStatus = new MaterialStatus(this, e); return newStatus; };
    string giveName() { return name; }
    virtual void init() {};
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
    double effConductivity;
public:
    TrsprtMaterialStatus(TrsprtMaterial *m, Element *e);
    virtual ~TrsprtMaterialStatus() {};
    virtual Vector giveStress(const Vector &strain); //terminology from mechanics, it returns flux
    virtual Vector giveStressWithFrozenIntVars(const Vector &strain) const;
    virtual double giveEffectiveConductivity(string type) const;
    virtual Matrix giveStiffnessTensor(string type, unsigned dimension) const;
    virtual double giveMassConstant() const;
};

//////////////////////////////////////////////////////////
class TrsprtMaterial : public Material
{
protected:
    double permeability, viscosity, capacity, density;
public:
    TrsprtMaterial() { name = "transport material"; };
    ~TrsprtMaterial() {};
    double giveCapacity() { return capacity; };
    double giveDensity() { return density; };
    double givePermeability() { return permeability; };
    double giveViscosity() { return viscosity; };
    void readFromLine(istringstream &iss);
    MaterialStatus *giveNewMaterialStatus(Element *e);
};

//////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////
// TRANSPORT MATERIAL WITH COUPLING TERMS

class TrsprtCoupledMaterial;
class TrsprtCoupledMaterialStatus : public TrsprtMaterialStatus
{
protected:
    double temp_effConductivity;
public:
    TrsprtCoupledMaterialStatus(TrsprtMaterial *m, Element *e);
    virtual ~TrsprtCoupledMaterialStatus() {};
    virtual Vector giveStress(const Vector &strain);
    virtual Vector giveStressWithFrozenIntVars(const Vector &strain) const;
    virtual double giveEffectiveConductivity(string type) const;
    virtual void update();
};

//////////////////////////////////////////////////////////
class TrsprtCoupledMaterial : public TrsprtMaterial
{
private:
    double crack_turtuosity;
public:
    TrsprtCoupledMaterial() { name = "coupled transport material"; };
    ~TrsprtCoupledMaterial() {};
    void readFromLine(istringstream &iss);
    MaterialStatus *giveNewMaterialStatus(Element *e);
    double giveTurtuosity() { return crack_turtuosity; };
};

//////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////
// TENSORIAL MECHANICAL MATERIAL

class ElasticMechMaterial;
class ElasticMechMaterialStatus : public MaterialStatus
{
protected:

public:
    ElasticMechMaterialStatus(ElasticMechMaterial *m, Element *e);
    virtual ~ElasticMechMaterialStatus() {};
    virtual Matrix giveStiffnessTensor(string type, unsigned dim) const;
    virtual Vector giveStress(const Vector &strain);
    virtual Vector giveStressWithFrozenIntVars(const Vector &strain) const;
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
    virtual MaterialStatus *giveNewMaterialStatus(Element *e);
    double giveElasticModulus() const { return E; }
    double givePoissonsRatio() const { return nu; }
    double giveDensity() { return density; };
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
    CosseratMechMaterialStatus(CosseratMechMaterial *m, Element *e);
    virtual ~CosseratMechMaterialStatus() {};
    virtual Matrix giveStiffnessTensor(string type, unsigned dim) const;
    virtual Vector giveStress(const Vector &strain);
    virtual Vector giveStressWithFrozenIntVars(const Vector &strain) const;
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
    virtual MaterialStatus *giveNewMaterialStatus(Element *e);
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
    DisMechMaterialStatus(DisMechMaterial *m, Element *e);
    virtual ~DisMechMaterialStatus() {};
    virtual Matrix giveStiffnessTensor(string type, unsigned dim) const;
    virtual Vector giveStress(const Vector &strain);
    virtual Vector giveStressWithFrozenIntVars(const Vector &strain) const;
    double giveDensity() const;
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
    virtual MaterialStatus *giveNewMaterialStatus(Element *e);
    double giveAlpha() const { return alpha; }
    double giveE0() const { return E0; }
};

#endif /* _MATERIAL_H */
