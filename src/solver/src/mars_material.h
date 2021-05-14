#ifndef _MARS_MATERIAL_H
#define _MARS_MATERIAL_H

#include "material.h"

//////////////////////////////////////////////////////////
// CUSATIS/MARS MATERIAL

class MarsMaterial;
class MarsMaterialStatus : public DisMechMaterialStatus
{
private:
    double omega0, maxEpsT, maxEpsN, temp_maxEpsT, temp_maxEpsN;
    // Vector temp_strain;
    double damage, temp_damage;
    double Kt, Ks, L, nt;
    double RAND_H;
    double temp_crackOpening;
    double volumetricStrain;

    double giveS0tension(double omega) const;
    double giveS0compression(double omega) const;
    void computeOmega0();
    void computeDamage(Vector strain);
    void computeKsAnsKt();
public:
    MarsMaterialStatus(MarsMaterial *m, Element *e, unsigned ipnum);
    virtual ~MarsMaterialStatus() {};
    virtual void init();
    virtual void update();
    virtual Matrix giveStiffnessTensor(string type, unsigned dim) const;
    virtual Vector giveStress(const Vector &strain, double timeStep);
    virtual Vector giveStressWithFrozenIntVars(const Vector &strain, double timeStep);
    virtual double giveValue(string code) const;
    virtual std :: string giveLineToSave() const;
    virtual void readFromLine(istringstream &iss);
    virtual bool isElastic(const bool &now = false) const;
    virtual void setParameterValue(string code, double value);
};


class MarsMaterial : public DisMechMaterial
{
private:
    double ft, Gt;
    double fs, Gs, fc, Kc, beta, mu, nc;
    double Lcrs, Lcrt;
public:
    MarsMaterial() { name = "Mars material"; };
    virtual ~MarsMaterial() {};
    virtual void readFromLine(istringstream &iss);
    virtual MaterialStatus *giveNewMaterialStatus(Element *e, unsigned ipnum);
    double giveFt() { return ft; }
    double giveGt() { return Gt; }
    double giveFs() { return fs; }
    double giveGs() { return Gs; }
    double giveFc() { return fc; }
    double giveBeta() { return beta; }
    double giveMu() { return mu; }
    double giveNc() { return nc; }
    double giveLcrs() { return Lcrs; }
    double giveLcrt() { return Lcrt; }
    double giveKc() { return Kc; }

    virtual void init();
};



//////////////////////////////////////////////////////////
//COUPLED MARS MATERIAL
class CoupledMarsMaterial;
class CoupledMarsMaterialStatus : public MarsMaterialStatus
{
private:
    void updateStressByBiotEffect(double timeStep);
    double avgPressure;

public:
    CoupledMarsMaterialStatus(MarsMaterial *m, Element *e, unsigned ipnum);
    ~CoupledMarsMaterialStatus() {};
    virtual Vector giveStress(const Vector &strain, double timeStep);
    virtual Vector giveStressWithFrozenIntVars(const Vector &strain, double timeStep);
    virtual double giveValue(string code) const;
    virtual void init();
    virtual void update();
    virtual void setParameterValue(string code, double value);
};

class CoupledMarsMaterial : public MarsMaterial
{
private:
    double biotCoeff;
public:
    CoupledMarsMaterial() { name = "Coupled Mars material"; };
    virtual ~CoupledMarsMaterial() {};
    virtual void readFromLine(istringstream &iss);
    virtual MaterialStatus *giveNewMaterialStatus(Element *e, unsigned ipnum);
    virtual void init();
    double giveBiotCoeff() const { return biotCoeff; };
};

#endif /* _MARS_MATERIAL_H */
