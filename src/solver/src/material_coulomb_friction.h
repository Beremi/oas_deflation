#ifndef _COULOMB_FRICTION_MATERIAL_H
#define _COULOMB_FRICTION_MATERIAL_H

#include "material.h"

//////////////////////////////////////////////////////////
// COULOMB FRICTION MATERIAL

class CoulombFrictionMaterial;
class CoulombFrictionMaterialStatus : public DisMechMaterialStatus
{
private:
    double normalStress;
public:
    CoulombFrictionMaterialStatus(CoulombFrictionMaterial *m, Element *e, unsigned ipnum);
    virtual ~CoulombFrictionMaterialStatus() {};
    virtual void init();
    virtual void update();
    virtual void resetTemporaryVariables();
    virtual Matrix giveStiffnessTensor(std :: string type, unsigned dim) const;
    virtual Vector giveStress(const Vector &strain, double timeStep);
    virtual Vector giveStressWithFrozenIntVars(const Vector &strain, double timeStep);
    virtual void setParameterValue(std :: string code, double value);
    virtual void giveValues(std :: string code, Vector &result) const;
};


class CoulombFrictionMaterial : public DisMechMaterial
{
private:
    double friction_angle, init_stiffness;
public:
    CoulombFrictionMaterial() { name = "Coulomb friction material"; };
    virtual ~CoulombFrictionMaterial() {};
    virtual void readFromLine(std :: istringstream &iss);
    virtual CoulombFrictionMaterialStatus *giveNewMaterialStatus(Element *e, unsigned ipnum);
    double giveFrictionAngle() const { return friction_angle; };
    double giveInitialStiffness() const { return init_stiffness; };
    virtual void init();
};

#endif /* _COULOMB_FRICTION_MATERIAL_H */
