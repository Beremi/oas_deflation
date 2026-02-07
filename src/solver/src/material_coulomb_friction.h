#ifndef _COULOMB_FRICTION_MATERIAL_H
#define _COULOMB_FRICTION_MATERIAL_H

#include "material.h"
#include "material_vectorial.h"

//////////////////////////////////////////////////////////
// COULOMB FRICTION MATERIAL

class CoulombFrictionMaterial;
class CoulombFrictionMaterialStatus : public VectMechMaterialStatus
{
private:
    double normalStress;
public:
    CoulombFrictionMaterialStatus(CoulombFrictionMaterial *m, Element *e, unsigned ipnum);
    virtual ~CoulombFrictionMaterialStatus() {};
    virtual void init();
    virtual void update();
    virtual void resetTemporaryVariables();
    virtual Matrix giveStiffnessTensor(std :: string type) const;
    virtual void computeStress( double timeStep);
    virtual void computeStressWithFrozenIntVars( double timeStep);
    virtual void setParameterValue(std :: string code, double value);
    virtual bool giveValues(std :: string code, Vector &result) const;
};


class CoulombFrictionMaterial : public VectMechMaterial
{
private:
    double friction_angle, init_stiffness;
public:
    CoulombFrictionMaterial(unsigned dimension) : VectMechMaterial(dimension) { name = "Coulomb friction material"; };
    virtual ~CoulombFrictionMaterial() {};
    virtual void readFromLine(std :: istringstream &iss);
    virtual CoulombFrictionMaterialStatus *giveNewMaterialStatus(Element *e, unsigned ipnum);
    double giveFrictionAngle() const { return friction_angle; };
    double giveInitialStiffness() const { return init_stiffness; };
    virtual void init(MaterialContainer *matcont);
};

#endif /* _COULOMB_FRICTION_MATERIAL_H */
