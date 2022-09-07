#ifndef _MATERIAL_FIBER_H
#define _MATERIAL_FIBER_H

#include "material.h"

//////////////////////////////////////////////////////////
// MATERIAL FIBER

class FiberMaterial;
class FiberMaterialStatus : public MaterialStatus
{
private:
    double crack_opening;

public:
    FiberMaterialStatus(FiberMaterial *m, Element *e, unsigned ipnum);
    virtual ~FiberMaterialStatus() {};
    virtual void init();
    virtual void update();
    virtual void resetTemporaryVariables();
    virtual Matrix giveStiffnessTensor(std :: string type, unsigned dim) const;
    virtual Vector giveStress(const Vector &strain, double timeStep);
    virtual Vector giveStressWithFrozenIntVars(const Vector &strain, double timeStep);
    virtual void giveValues(std :: string code, Vector &result) const;
    virtual void setParameterValue(std :: string code, double value);
};


class FiberMaterial : public Material
{
private:
    double Ef, Gd, tau0, betaf, ft;
public:
    FiberMaterial() { name = "fiber material"; };
    virtual ~FiberMaterial() {};
    virtual void readFromLine(std :: istringstream &iss);
    virtual MaterialStatus *giveNewMaterialStatus(Element *e, unsigned ipnum);
    double giveEf() { return Ef; }
    double giveGd() { return Gd; }
    double giveTau0() { return tau0; }
    double giveBetaf() { return betaf; }
    double giveFt() { return ft; }
    virtual void init();
};
#endif /* _LDPM_MATERIAL_H */
