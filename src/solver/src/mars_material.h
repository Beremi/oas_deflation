#ifndef _MARS_MATERIAL_H
#define _MARS_MATERIAL_H

#include "material.h"

//////////////////////////////////////////////////////////
// CUSATIS/MARS MATERIAL 2007

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

    double crackOpening;
public:
    MarsMaterialStatus(MarsMaterial *m, Element *e, unsigned ipnum);
    virtual ~MarsMaterialStatus() {};
    virtual void init();
    virtual void update();
    virtual void resetTemporaryVariables();
    virtual Matrix giveStiffnessTensor(string type, unsigned dim) const;
    virtual Vector giveStress(const Vector &strain, double timeStep);
    virtual Vector giveStressWithFrozenIntVars(const Vector &strain, double timeStep);
    virtual void giveValues(string code, Vector &result) const;
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
    double damage_residuum = 0.0;
    double stress_residuum_fraction = 0.0;
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
    double giveDamageResiduum() { return damage_residuum; }
    double giveStressResiduum() { return ft * stress_residuum_fraction; }

    virtual void init();
};

//////////////////////////////////////////////////////////
// LDPM MATERIAL (2011)

class LDPMMaterial;
class LDPMMaterialStatus : public MarsMaterialStatus
{
private:

public:
    LDPMMaterialStatus(LDPMMaterial *m, Element *e, unsigned ipnum);
    virtual ~LDPMMaterialStatus() {};
    void computeDamage(Vector strain);
    virtual Vector giveStress(const Vector &strain, double timeStep);
};


class LDPMMaterial : public MarsMaterial
{
private:

public:
    LDPMMaterial() { name = "LDPM material"; };
    virtual ~LDPMMaterial() {};
    virtual void readFromLine(istringstream &iss);
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
    virtual void giveValues(string code, Vector &result) const;
    virtual void init();
    virtual void update();
    virtual void resetTemporaryVariables();
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
