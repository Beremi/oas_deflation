#ifndef _MATERIAL_THERMOMECHANICAL_H
#define _MATERIAL_THERMOMECHANICAL_H

#include "material.h"

class MaterialContainer; //forward declaration

//////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////
// DISCRETE COUPLED MECHANICS AND HEAT CONDUCTION MATERIAL

class ThermoMechanicalMaterial;
class ThermoMechanicalMaterialStatus : public CoupledMaterialStatus
{
protected:
    double temperature;
    void addTemperatureEffectToMechanics();
public:
    ThermoMechanicalMaterialStatus(ThermoMechanicalMaterial *m, Element *e, unsigned ipnum);
    virtual ~ThermoMechanicalMaterialStatus() {};
    bool giveValues(std :: string code, Vector &result) const;
    virtual void setParameterValue(std :: string code, double value);
    virtual Vector giveStress(const Vector &strain, double timeStep);
    virtual Vector giveStressWithFrozenIntVars(const Vector &strain, double timeStep);
};

//////////////////////////////////////////////////////////
class ThermoMechanicalMaterial : public CoupledMaterial
{
protected:
    double tec;  //thermal expansion ceofficient
    double initialTemp;
public:
    ThermoMechanicalMaterial(unsigned dimension) : CoupledMaterial(dimension) { name = "generic discrete coupled  material"; initialTemp = 0; };
    ~ThermoMechanicalMaterial() {};
    void init(MaterialContainer *matcont);
    double giveThermalExpansionCoeff() const { return tec; };
    void setThermalExpansionCoeff(double value) { tec = value; };
    double giveInitialTemperature() const { return initialTemp; };
    virtual void readFromLine(std :: istringstream &iss);
    virtual MaterialStatus *giveNewMaterialStatus(Element *e, unsigned ipnum);
};

#endif /* _MATERIAL_THERMOMECHANICAL_H */
