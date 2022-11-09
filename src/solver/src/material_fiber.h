#ifndef _MATERIAL_FIBER_H
#define _MATERIAL_FIBER_H

#include "material.h"

//////////////////////////////////////////////////////////
// MATERIAL FIBER

class FiberMaterial;
class FiberMaterialStatus : public ElasticMechMaterialStatus
{
private:
    double crack_opening;
    double temp_crack_opening;
    double incrementOfCrack;
    
    double contactLength; 
    Point contactNormal;
    double df, rightLe, leftLe; 
    Point fiberNormal;
    double inclineAngle; 
    
    double right_pullout;
    double temp_rightPullout;
    double left_pullout; 
    double temp_leftPullout; 
    
    double bridgingForce;
    double temp_bridgingForce;
    double rightForce, leftForce;
    
    double limit_rightPullout, limit_leftPullout, right_F0, left_F0;
    
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


class FiberMaterial : public ElasticMechMaterial
{
private:
    double Ef, Gd, tau0, betaf, ft, f;
public:
    FiberMaterial(); //{ name = "fiber material"; };
    virtual ~FiberMaterial() {};
    virtual void readFromLine(std :: istringstream &iss);
    virtual MaterialStatus *giveNewMaterialStatus(Element *e, unsigned ipnum);
    double giveEf() { return Ef; }
    double giveGd() { return Gd; }
    double giveTau0() { return tau0; }
    double giveBetaf() { return betaf; }
    double giveFt() { return ft; }
    double giveF() { return f; }
    virtual void init();
};

double bridgingForce_bonded( double v, double df, double Ef, double tau0, double Gd ); 
double derivF_bonded( double v, double df, double Ef, double tau0, double Gd ); 
double bridgingForce_debonded( double v, double vd, double Le, double F0, double df, double betaf );
double derivF_debonded( double v, double vd, double Le, double F0, double df, double betaf ); 
double bridgingForce_unloading( double deltaX, double deltaY, double x ); 
double derivF_unloading( double deltaX, double deltaY ); 

#endif /* _LDPM_MATERIAL_H */
