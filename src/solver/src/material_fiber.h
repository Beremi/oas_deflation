#ifndef _MATERIAL_FIBER_H
#define _MATERIAL_FIBER_H

#include "material_tensorial.h"

//////////////////////////////////////////////////////////
// MATERIAL FIBER

class FiberMaterial;
class FiberMaterialStatus : public TensMechMaterialStatus
{
private:
    double crackOpening, temp_crackOpening, maxCrackOpening, incrementOfCrack;
    Vector crackOpeningVector;
    
    double rightPullout, temp_rightPullout, leftPullout, temp_leftPullout;
    
    double bridgingForce, temp_bridgingForce, fiberForce, temp_fiberForce; 
    double rightForce, leftForce;
    
    double inclineAngle, deflectionAngle, spallingLength, temp_spallingLength; 
    Vector w, mf;
    
    int debondedFiber, temp_debondedFiber, pulledOutFiber, temp_pulledOutFiber, rupturedFiber, temp_rupturedFiber, closingCrack, temp_closingCrack;
    
    double df, rightLe, leftLe, fiberLength;
    
    Point fiberNormal, contactNormal;
    Vector fiberNormalLocal;
    double contactLength;

    double right_F0, left_F0;
    double criticalRightPullout, criticalLeftPullout;
    
    double Le1, F01, v1, vd1, Le2, F02, v2, vd2, limitPullout;
    double bridgingForce1, bridgingForce2;
    
    double A1, B1, C1, v1_Fmax, Fmax; 
    double A2, B2, C2, v2_Fmax1, v2_Fmax2, v2_pullingOutFmax, v2_debondingFmax;
    
    
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

//////////////////////////////////////////////////////////
// ------- EQUATIONS FOR FIBER CONSTITUTIVE LAW ------- //
//////////////////////////////////////////////////////////
//------------------------------------------------------//
//                   DEBONDING STAGE                    //
//------------------------------------------------------//
// -------- ORIGINAL EQUATIONS FROM LITERATURE -------- // with zero crack, there is non-zero bridging force in fiber due to the fracture energy Gd
//double bridgingForceDebonding( double v, double vd, double df, double Ef, double tau0, double Gd );
//double derivBridgingForceDebonding( double v, double vd, double df, double Ef, double tau0, double Gd );
// ---------------- MODIFIED EQUATIONS ---------------- // linear distribution of Gd
double bridgingForceDebonding( double v, double vd, double df, double Ef, double tau0, double Gd );
double derivBridgingForceDebonding( double v, double vd, double df, double Ef, double tau0, double Gd );
//------------------------------------------------------//
//                  PULLING OUT STAGE                   //
//------------------------------------------------------//
double bridgingForcePullingOut( double v, double vd, double Le, double F0, double df, double betaf );
double derivBridgingForcePullingOut( double v, double vd, double Le, double F0, double df, double betaf );
//------------------------------------------------------//
//                   UNLOADING STAGE                    //
//------------------------------------------------------//
double bridgingForceUnloading( double deltaX, double deltaY, double x );
double derivBridgingForceUnloading( double deltaX, double deltaY );

#endif /* _LDPM_MATERIAL_H */
