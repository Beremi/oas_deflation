#ifndef _MATERIAL_CONTINUOUS_H
#define _MATERIAL_CONTINUOUS_H

#include "linalg.h"
#include "material.h"
#include <vector>
#include <iostream>
#include <fstream>


//////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////
// TRANSPORT MATERIAL

class TensTrsprtMaterial;
class TensTrsprtMaterialStatus : public MaterialStatus
{
protected:
    double effConductivity, temp_effConductivity;
    double temp_pressure;
public:
    TensTrsprtMaterialStatus(TensTrsprtMaterial *m, Element *e, unsigned ipnum);
    virtual ~TensTrsprtMaterialStatus() {};
    virtual Vector giveStress(const Vector &strain, double timeStep); //terminology from mechanics, it returns flux
    virtual Vector giveStressWithFrozenIntVars(const Vector &strain, double timeStep);
    virtual Matrix giveStiffnessTensor(std :: string type) const;
    virtual Matrix giveDampingTensor() const;
    virtual bool giveValues(std :: string code, Vector &result) const;
    virtual double giveEffectiveConductivity(std :: string type) const;
    virtual double updateEffectiveConductivity() const;
    virtual double calculatePressureDependentPermeability(double pressure) const;
    virtual bool isElastic(const bool &now = false) const;
    virtual void setParameterValue(std :: string code, double value);
};

//////////////////////////////////////////////////////////
class TensTrsprtMaterial : public Material
{
protected:
    double permeability, viscosity, capacity, density, a, m;
public:
    TensTrsprtMaterial(unsigned dimension) : Material(dimension) { name = "transport material"; a = -1.; m = 0; };
    ~TensTrsprtMaterial() {};
    double giveCapacity() const { return capacity; };
    double giveDensity() const { return density; };
    double givePermeability() const { return permeability; };
    double giveViscosity() const { return viscosity; };
    double giveParamA() const { return a; };
    double giveParamM() const { return m; };
    void setPermeability(double new_p) { permeability = new_p; };
    void setParamA(double new_a) { a = new_a; };
    void readFromLine(std :: istringstream &iss);
    MaterialStatus *giveNewMaterialStatus(Element *e, unsigned ipnum);
    virtual void init(MaterialContainer *matcont) { Material :: init(matcont); strainsize = dim; }
};

//////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////
// HEAT CONDUCTION MATERIAL

class TensHeatConductionMaterial;
class TensHeatConductionMaterialStatus : public MaterialStatus
{
protected:
    double effConductivity, temp_effConductivity;
    double temp_temperature;
public:
    TensHeatConductionMaterialStatus(TensHeatConductionMaterial *m, Element *e, unsigned ipnum);
    virtual ~TensHeatConductionMaterialStatus() {};
    virtual Vector giveStress(const Vector &strain, double timeStep); //terminology from mechanics, it returns flux
    virtual Vector giveStressWithFrozenIntVars(const Vector &strain, double timeStep);
    virtual Matrix giveStiffnessTensor(std :: string type) const;
    virtual Matrix giveDampingTensor() const;
    virtual bool giveValues(std :: string code, Vector &result) const;
    virtual void setParameterValue(std :: string code, double value);
};

//////////////////////////////////////////////////////////
class TensHeatConductionMaterial : public Material
{
protected:
    double conductivity, capacity, density;
public:
    TensHeatConductionMaterial(unsigned dimension) : Material(dimension) { name = "heat conduction material"; };
    ~TensHeatConductionMaterial() {};
    double giveCapacity() const { return capacity; };
    double giveDensity() const { return density; };
    double giveConductivity() const { return conductivity; };
    void readFromLine(std :: istringstream &iss);
    MaterialStatus *giveNewMaterialStatus(Element *e, unsigned ipnum);
    virtual void init(MaterialContainer *matcont) { Material :: init(matcont); strainsize = dim; }
};

//////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////
// TENSORIAL MECHANICAL MATERIAL

class TensMechMaterial;
class TensMechMaterialStatus : public MaterialStatus
{
protected:

public:
    TensMechMaterialStatus(TensMechMaterial *m, Element *e, unsigned ipnum);
    virtual ~TensMechMaterialStatus() {};
    virtual Matrix giveStiffnessTensor(std :: string type) const;
    virtual Vector giveStress(const Vector &strain, double timeStep);
    virtual Vector giveStressWithFrozenIntVars(const Vector &strain, double timeStep);
    virtual double giveMassConstant() const;
    virtual bool giveValues(std :: string code, Vector &result) const;
    virtual void update();
};

//////////////////////////////////////////////////////////
class TensMechMaterial : public Material
{
protected:
    double E, nu, density;
    bool planeStress;

public:
    TensMechMaterial(unsigned dimension) : Material(dimension) { name = "elastic tensorial mechanical material"; planeStress = true; };
    ~TensMechMaterial() {};
    virtual void readFromLine(std :: istringstream &iss);
    virtual MaterialStatus *giveNewMaterialStatus(Element *e, unsigned ipnum);
    double giveElasticModulus() const { return E; }
    double givePoissonsRatio() const { return nu; }
    double giveDensity() const { return density; };
    bool isPlaneStress() { return planeStress; };
    virtual void init(MaterialContainer *matcont) { Material :: init(matcont); strainsize = ( dim - 1 ) * 3; }
};

//////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////
// COSSERAT MECHANICAL MATERIAL

class TensCosseratMechMaterial;
class TensCosseratMechMaterialStatus : public TensMechMaterialStatus
{
protected:

public:
    TensCosseratMechMaterialStatus(TensCosseratMechMaterial *m, Element *e, unsigned ipnum);
    virtual ~TensCosseratMechMaterialStatus() {};
    virtual Matrix giveStiffnessTensor(std :: string type) const;
    virtual Vector giveStress(const Vector &strain, double timeStep);
    virtual Vector giveStressWithFrozenIntVars(const Vector &strain, double timeStep);
};

//////////////////////////////////////////////////////////
class TensCosseratMechMaterial : public TensMechMaterial
{
protected:
    double lc, muc;

public:
    TensCosseratMechMaterial(unsigned dimension) : TensMechMaterial(dimension) { name = "elastic Cosserat mechanical material"; planeStress = true; };
    ~TensCosseratMechMaterial() {};
    virtual void readFromLine(std :: istringstream &iss);
    virtual MaterialStatus *giveNewMaterialStatus(Element *e, unsigned ipnum);
    double giveCharacteristicLength() const { return lc; }
    double giveCosseratShearParam() const { return muc; }
    virtual void init(MaterialContainer *matcont) { TensMechMaterial :: init(matcont); strainsize = dim * dim; }
};

#endif /* _MATERIAL_CONTINUOUS_H */
