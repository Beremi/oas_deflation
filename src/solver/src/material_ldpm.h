#ifndef _LDPM_MATERIAL_H
#define _LDPM_MATERIAL_H

#include "material_vectorial.h"

//////////////////////////////////////////////////////////
// LDPM MATERIAL 2011

class LDPMMaterial;
class LDPMMaterialStatus : public VectMechMaterialStatus
{
private:
    double maxEpsT, maxEpsN, temp_maxEpsT, temp_maxEpsN;
    double Kt, Ks, L, nt;
    double RAND_H;
    double crackOpening, temp_crackOpening;
    double volumetricStrain;
    double virtual_damage;
    Vector updt_mech_strain; //last strain without eigenstrain
    Vector temp_mech_strain; //current strain without eigenstrain

    double giveStrengthLimit(double omega);
    Vector giveTension(const Vector &strain, const Vector strain_prev, const Vector stress_prev);
    Vector giveCompression(const Vector &strain, const Vector strain_prev, const Vector stress_prev);
    Vector passZero(const Vector &strain);
    Vector passThroughZero(const Vector &strain);
    void giveVirtualDamage();

public:
    LDPMMaterialStatus(LDPMMaterial *m, Element *e, unsigned ipnum);
    virtual ~LDPMMaterialStatus() {};
    virtual void init();
    virtual bool giveValues(std :: string code, Vector &result) const;
    virtual Vector giveStress(const Vector &strain, double timeStep);
    virtual Vector giveStressWithFrozenIntVars(const Vector &strain, double timeStep);
    virtual void update();
    virtual void resetTemporaryVariables();
    virtual Matrix giveStiffnessTensor(std :: string type) const;
    virtual std :: string giveLineToSave() const;
    virtual void setParameterValue(std :: string code, double value);
    virtual void readFromLine(std :: istringstream &iss);
    virtual bool isElastic(const bool &now = false) const;
    virtual void initializeStressAndStrainVector(unsigned num);
};


class LDPMMaterial : public VectMechMaterial
{
private:
    double ft, Gt;
    double nt, kt, beta;
    double fc, fc0, Ed, Hc0, Hc1, Kc0, Kc1, Kc2, Kc3;
    double fs, fs0, Et, mu0, muinf;
    double damage_residuum;
    double stress_residuum_fraction = 0.0;
public:
    LDPMMaterial(unsigned dimension) : VectMechMaterial(dimension) { name = "LDPM material"; };
    virtual ~LDPMMaterial() {};
    virtual void init(MaterialContainer *matcont);
    virtual void readFromLine(std :: istringstream &iss);
    virtual MaterialStatus *giveNewMaterialStatus(Element *e, unsigned ipnum);
    double giveFt() { return ft; }
    double giveGt() { return Gt; }
    double givent() { return nt; }
    double givekt() { return kt; }
    double giveBeta() { return beta; }
    double giveFc() { return fc; }
    double giveFc0() { return fc0; }
    double giveEd() { return Ed; }
    double giveHc0() { return Hc0; }
    double giveHc1() { return Hc1; }
    double giveKc0() { return Kc0; }
    double giveKc1() { return Kc1; }
    double giveKc2() { return Kc2; }
    double giveKc3() { return Kc3; }
    double giveFs() { return fs; }
    double giveFs0() { return fs0; }
    double giveEt() { return Et; }
    double giveMu0() { return mu0; }
    double giveMuinf() { return muinf; }
    double giveDamageResiduum() { return damage_residuum; }
    double giveStressResiduum() { return ft * stress_residuum_fraction; }
};
#endif /* _LDPM_MATERIAL_H */
