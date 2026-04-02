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
    double effConductivity;
    double temp_pressure;
public:
    TensTrsprtMaterialStatus(TensTrsprtMaterial *m, Element *e, unsigned ipnum);
    virtual ~TensTrsprtMaterialStatus() {};
    virtual void computeStress(double timeStep);  //terminology from mechanics, it returns flux
    virtual void computeStressWithFrozenIntVars(double timeStep);
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
    double permeability; //[m2]
    double viscosity; //[Pa s]
    double capacity; //[1/Pa = m s2/kg]
    double density; //[kg/m3]
    double a, m;
public:
    TensTrsprtMaterial(unsigned dimension) : Material(dimension) { name = "transport material"; a = -1.; m = 0; capacity=-1;};
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
    double effConductivity;
    double temp_temperature;
public:
    TensHeatConductionMaterialStatus(TensHeatConductionMaterial *m, Element *e, unsigned ipnum);
    virtual ~TensHeatConductionMaterialStatus() {};
    virtual void computeStress(double timeStep);  //terminology from mechanics, it returns flux
    virtual void computeStressWithFrozenIntVars(double timeStep);
    virtual Matrix giveStiffnessTensor(std :: string type) const;
    virtual Matrix giveDampingTensor() const;
    virtual bool giveValues(std :: string code, Vector &result) const;
    virtual void setParameterValue(std :: string code, double value);
    virtual double updateEffectiveConductivity() const;    
};

//////////////////////////////////////////////////////////
class TensHeatConductionMaterial : public Material
{
protected:
    double conductivity; //[J/s/m/K]
    double capacity; //[J/K/Kg]
    double density; //[kg/m3]
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
    Matrix giveElasticStiffnessTensor3D() const;
public:
    TensMechMaterialStatus(TensMechMaterial *m, Element *e, unsigned ipnum);
    virtual ~TensMechMaterialStatus() {};
    virtual Matrix giveStiffnessTensor(std :: string type) const;
    virtual void computeStress(double timeStep);
    virtual void computeStressWithFrozenIntVars(double timeStep);
    virtual double giveMassConstant() const;
    virtual bool giveValues(std :: string code, Vector &result) const;
    virtual void update();
    virtual Matrix giveMassTensor() const;
    virtual Matrix giveDampingTensor() const;
};

//////////////////////////////////////////////////////////
class TensMechMaterial : public Material
{
protected:
    double E, nu, density;
    bool planeStress;

public:
    TensMechMaterial(unsigned dimension) : Material(dimension) { name = "elastic tensorial mechanical material"; planeStress = true; strainsize = ( dim - 1 ) * 3; };
    ~TensMechMaterial() {};
    virtual void readFromLine(std :: istringstream &iss);
    virtual MaterialStatus *giveNewMaterialStatus(Element *e, unsigned ipnum);
    double giveElasticModulus() const { return E; }
    double giveShearModulus() const { return E / ( 2. + 2. * nu ); }
    double givePoissonsRatio() const { return nu; }
    double giveDensity() const { return density; };
    bool isPlaneStress() { return planeStress; };
    virtual void init(MaterialContainer *matcont) { Material :: init(matcont);  }
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
    virtual void computeStress(double timeStep);
    virtual void computeStressWithFrozenIntVars(double timeStep);
};

//////////////////////////////////////////////////////////
class TensCosseratMechMaterial : public TensMechMaterial
{
protected:
    double lc, muc, psize;

public:
    TensCosseratMechMaterial(unsigned dimension);
    ~TensCosseratMechMaterial() {};
    virtual void readFromLine(std :: istringstream &iss);
    virtual MaterialStatus *giveNewMaterialStatus(Element *e, unsigned ipnum);
    double giveCharacteristicLength() const { return lc; }
    double giveCosseratShearParam() const { return muc; }
    double giveParticleSize() const { return psize; }
    virtual void init(MaterialContainer *matcont) { TensMechMaterial :: init(matcont); }
};

#endif /* _MATERIAL_CONTINUOUS_H */
