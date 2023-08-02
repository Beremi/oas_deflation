#ifndef _MATERIAL_FIBER_H
#define _MATERIAL_FIBER_H

#include "material_tensorial.h"

//////////////////////////////////////////////////////////
// MATERIAL FIBER

class FiberMaterial;
class FiberMaterialStatus : public TensMechMaterialStatus
{
private:
    double crack_opening, temp_crack_opening, incrementOfCrack, maxCrackOpening;
    Vector crackOpeningVector;

    double right_pullout, temp_rightPullout, left_pullout, temp_leftPullout;

    double bridgingForce, temp_bridgingForce, rightForce, leftForce;
    Vector stress;

    double spallingLength, temp_spallingLength, deflectionAngle, fiberForce;
    Vector w, nf;

    int debondedFiber, temp_debondedFiber, pulloutOfFiber, temp_pulloutOfFiber, ruptureOfFiber, temp_ruptureOfFiber, closingCrack, temp_closingCrack;

    double df, rightLe, leftLe;
    Point fiberNormal;
    Vector fiberNormalLocal;

    double contactLength;
    Point contactNormal;

    double inclineAngle;

    double right_F0, left_F0, limit_rightPullout, limit_leftPullout;

    double A, B, C, v0, v_Fmax, Fmax, AA, BB, CC, v21, v22, v2Bonded_Fmax, v2Debonded_Fmax;
    double Le1, F01, v1, vd1, Le2, F02, v2, vd2;
    double bridgingForce1, bridgingForce2;

public:
    FiberMaterialStatus(FiberMaterial *m, Element *e, unsigned ipnum);
    virtual ~FiberMaterialStatus() {};
    virtual void init();
    virtual void update();
    virtual void resetTemporaryVariables();
    virtual Matrix giveStiffnessTensor(std :: string type) const;
    virtual Vector giveStress(const Vector &strain, double timeStep);
    virtual Vector giveStressWithFrozenIntVars(const Vector &strain, double timeStep);
    virtual bool giveValues(std :: string code, Vector &result) const;
    virtual void setParameterValue(std :: string code, double value);
};


class FiberMaterial : public TensMechMaterial
{
private:
    double Ef, Gd, tau0, betaf, ft, Ksn, Ksp, Krup;
public:
    FiberMaterial(unsigned dimension);
    virtual ~FiberMaterial() {};
    virtual void readFromLine(std :: istringstream &iss);
    virtual MaterialStatus *giveNewMaterialStatus(Element *e, unsigned ipnum);
    double giveEf() { return Ef; }
    double giveGd() { return Gd; }
    double giveTau0() { return tau0; }
    double giveBetaf() { return betaf; }
    double giveFt() { return ft; }
    double giveKsn() { return Ksn; }
    double giveKsp() { return Ksp; }
    double giveKrup() { return Krup; }
    virtual void init(MaterialContainer *matcont);
};

double bridgingForce_bonded(double v, double vd, double df, double Ef, double tau0, double Gd);
double derivF_bonded(double v, double vd, double df, double Ef, double tau0, double Gd);
double bridgingForce_debonded(double v, double vd, double Le, double F0, double df, double betaf);
double derivF_debonded(double v, double vd, double Le, double F0, double df, double betaf);
double bridgingForce_unloading(double deltaX, double deltaY, double x);
double derivF_unloading(double deltaX, double deltaY);

#endif /* _LDPM_MATERIAL_H */
