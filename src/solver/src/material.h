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
    virtual void update() {};
    virtual Vector giveStress(const Vector &strain) { ( void ) strain; return Vector(0); };
    virtual double giveValue(string code) const { ( void ) code; return 0; };
protected:
    Element *element;
    string name;
    Material *mat;
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
    virtual double giveEffectiveConductivity(string type) const { return effConductivity; };
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
// DISCRETE MECHANICAL LINEAR MATERIAL

class DisMechMaterial;
class DisMechMaterialStatus : public MaterialStatus
{
protected:

public:
    DisMechMaterialStatus(DisMechMaterial *m, Element *e);
    virtual ~DisMechMaterialStatus() {};
    Vector giveElasticNormalShearStiffness() const;
    virtual Vector giveNormalShearStiffness(string type) const { ( void ) type; return giveElasticNormalShearStiffness(); }   //only elastic
    virtual Vector giveNormalShearStiffness(string type, const Vector &strain) { giveStress(strain); return giveNormalShearStiffness(type); }
    virtual Vector giveStress(const Vector &strain);
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
