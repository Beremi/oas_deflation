#ifndef _CSL_MATERIAL_H
#define _CSL_MATERIAL_H

#include "material_vectorial.h"

//////////////////////////////////////////////////////////
// CSL MATERIAL 2007

class CSLMaterial;
class CSLMaterialStatus : public VectMechMaterialStatus
{
protected:
    double omega0, maxEpsT, maxEpsN, temp_maxEpsT, temp_maxEpsN;
    // MyVector temp_strain;
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
    CSLMaterialStatus(CSLMaterial *m, Element *e, unsigned ipnum);
    virtual ~CSLMaterialStatus() {};
    virtual void init();
    virtual void update();
    virtual void resetTemporaryVariables();
    virtual Matrix giveStiffnessTensor(std :: string type) const;
    virtual Vector giveStress(const Vector &strain, double timeStep);
    virtual Vector giveStressWithFrozenIntVars(const Vector &strain, double timeStep);
    virtual bool giveValues(std :: string code, Vector &result) const;
    virtual std :: string giveLineToSave() const;
    virtual void readFromLine(std :: istringstream &iss);
    virtual bool isElastic(const bool &now = false) const;
    virtual void setParameterValue(std :: string code, double value);
    Vector giveCrackOpeningVector() const;
};

//////////////////////////////////////////////////////////
class CSLMaterial : public VectMechMaterial
{
protected:
    double ft, Gt;
    double fs, Gs, fc, Kc, beta, mu, nc, lam0;
    double Lcrs, Lcrt;
    double damage_residuum = 0.0;
    double stress_residuum_fraction = 0.0;
public:
    CSLMaterial(unsigned dimension) : VectMechMaterial(dimension) { name = "CSL material"; lam0 = 1e10; }; //confinement removed
    virtual ~CSLMaterial() {};
    virtual void readFromLine(std :: istringstream &iss);
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
    double giveLam0() { return lam0; }
    double giveDamageResiduum() { return damage_residuum; }
    double giveStressResiduum() { return ft * stress_residuum_fraction; }

    virtual void init(MaterialContainer *matcont);
};

//////////////////////////////////////////////////////////
// CSL MATERIAL with update by tensorial stress

class CSLMaterialWithTensorialStressUpdate;
class CSLMaterialWithTensorialStressUpdateStatus : public CSLMaterialStatus
{
public:
    CSLMaterialWithTensorialStressUpdateStatus(CSLMaterialWithTensorialStressUpdate *m, Element *e, unsigned ipnum);
    virtual ~CSLMaterialWithTensorialStressUpdateStatus() {};
    virtual bool giveValues(std :: string code, Vector &result) const;
    virtual Vector giveStress(const Vector &strain, double timeStep);
    virtual Vector giveStressWithFrozenIntVars(const Vector &strain, double timeStep);
private:
    Vector giveEigenStrainFromTensorialStress();
};

//////////////////////////////////////////////////////////
class CSLMaterialWithTensorialStressUpdate : public CSLMaterial
{
private:
    double poisson;
    std :: vector< Vector >tensstress;
public:
    CSLMaterialWithTensorialStressUpdate(unsigned dimension);
    virtual ~CSLMaterialWithTensorialStressUpdate() {};
    virtual void readFromLine(std :: istringstream &iss);
    virtual MaterialStatus *giveNewMaterialStatus(Element *e, unsigned ipnum);
    double givePoissonNumber() { return poisson; }
    virtual void init(MaterialContainer *matcont);
    virtual void prepareForStressEvaluation(ElementContainer *elems);
    Vector giveAveragePrincipalStress(unsigned Anode, unsigned Bnode);
};


//////////////////////////////////////////////////////////
//COUPLED CSL MATERIAL
class CoupledCSLMaterial;
class CoupledCSLMaterialStatus : public CSLMaterialStatus
{
private:
    void updateStressByBiotEffect(double timeStep);
    double avgPressure;

public:
    CoupledCSLMaterialStatus(CSLMaterial *m, Element *e, unsigned ipnum);
    ~CoupledCSLMaterialStatus() {};
    virtual Vector giveStress(const Vector &strain, double timeStep);
    virtual Vector giveStressWithFrozenIntVars(const Vector &strain, double timeStep);
    virtual bool giveValues(std :: string code, Vector &result) const;
    virtual void init();
    virtual void update();
    virtual void resetTemporaryVariables();
    virtual void setParameterValue(std :: string code, double value);
};

//////////////////////////////////////////////////////////
class CoupledCSLMaterial : public CSLMaterial
{
private:
    double biotCoeff;
public:
    CoupledCSLMaterial(unsigned dimension) : CSLMaterial(dimension) { name = "Coupled CSL material"; };
    virtual ~CoupledCSLMaterial() {};
    virtual void readFromLine(std :: istringstream &iss);
    virtual MaterialStatus *giveNewMaterialStatus(Element *e, unsigned ipnum);
    virtual void init(MaterialContainer *matcont);
    double giveBiotCoeff() const { return biotCoeff; };
};


#endif /* _CSL_MATERIAL_H */
