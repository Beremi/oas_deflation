#ifndef _FATIGUE_MATERIAL_H
#define _FATIGUE_MATERIAL_H

#include "material.h"

//////////////////////////////////////////////////////////
// DAMAGE DRIVEN BY CUMMULATIVE SHEAR SLIP
//////////////////////////////////////////////////////////
// MATERIAL ACCORDING TO ALGORITHM FROM https://doi.org/10.1016/j.ijfatigue.2018.04.020

class FatigueShearMaterial;
class FatigueShearMaterialStatus : public DisMechMaterialStatus
{
private:
    Point slip; ///< slip
    Point sPi; ///< irreversible slip
    Point alphaKin;  ///< kinematic hardening variable
    double damageShear; ///< damage in tangential direction
    double zIso;  ///< isotropic hardening variable
    double tang_stiff;  ///< consistent algorithmic (= tangent) shear stiffness

    Point temp_sPi, temp_slip, temp_alphaKin;
    Point stressT, temp_stressT;
    double temp_damageShear, temp_zIso; ///<temporary variables

    double damageShear_set_from_the_outside = 0.0;

    double prev_damageShear, prev_zIso;
    Point prev_sPi, prev_alphaKin, prev_stressT, prev_slip;

    double Ynext;
    double lambda;
    double temp_lambda;

    double strain_slip_multiplier;
    double regularization_multiplier_area;

    bool checkReturnMap, useAnaliticalLambda, newIter, bisectionMeth;

    // the following are all densities ()
    double energy_PL, energy_D, energy_Kin, energy_Iso;
    double work_tot;

    void print() const;

    double coup_dam;
    bool comp_dam;
public:
    FatigueShearMaterialStatus(FatigueShearMaterial *m, Element *e, unsigned ipnum);
    virtual ~FatigueShearMaterialStatus() {};
    void init();
    virtual void update();
    virtual Matrix giveStiffnessTensor(string type, unsigned dim) const;
    virtual Vector giveStress(const Vector &strain, double timeStep);
    virtual Vector giveStressWithFrozenIntVars(const Vector &strain, double timeStep);
    virtual double giveValue(string code) const;
    double isDamageCoupled() const { return coup_dam; }
protected:
    void setDamage(const double &new_damage) {
        if ( new_damage > this->damageShear_set_from_the_outside ) {
            this->damageShear_set_from_the_outside = new_damage;
        }
    }
};


class FatigueShearMaterial : public DisMechMaterial
{
private:
    double tauBar; ///< reversibility limit
    double Kin;  ///< isotropic hardening modulus
    double gamma;  ///< kinematic hardening modulus
    double S;  ///< damage strength
    double c, r;  // parameters controling the damage acumullation, c >= 1.0
    double mC, mT;  ///< parameters controling the pressure sensitivity (under Compression or Tension)
    bool use_slip, check_retturn_mapping, analytical_lambda, newIterOn, bisecOn;
    double coup_dam;
    bool comp_dam;
    double comp_thresh = 0.0;
public:
    FatigueShearMaterial() { name = "Fatigue Shear material"; };
    ~FatigueShearMaterial() {};
    void readFromLine(istringstream &iss);
    MaterialStatus *giveNewMaterialStatus(Element *e, unsigned ipnum);
    double giveTauBar() const { return tauBar; }
    double giveKin() const { return Kin; }
    double giveGamma() const { return gamma; }
    double giveS() const { return S; }
    double giveC() const { return c; }
    double giveR() const { return r; }
    double giveMC() const { return mC; }
    double giveMT() const { return mT; }
    virtual void init();

    bool useSlip() const { return use_slip; }
    bool checkReturnMap() const { return check_retturn_mapping; }
    bool analyticalLambda() const { return analytical_lambda; }
    bool newIterativeApproachOn() const { return newIterOn; }
    bool bisectionMethOn() const { return bisecOn; }
    double isDamageCoupled() const { return coup_dam; }
    bool isCompressiveDamageOff() const { return comp_dam; }
    double giveCompressiveThreshold() const { return comp_thresh; }
};

//////////////////////////////////////////////////////////
// NORMAL DIRECTION: TENISON - DAMAGE, COMPRESSION - PLASTICITY
// according to paper http://www.eccm-ecfd2018.org/admin/files/filePaper/p924.pdf

class DamagePlasticMaterial;
class DamagePlasticMaterialStatus : public DisMechMaterialStatus
{
private:
    double epsN; ///< slip
    double damage; ///< damage in tangential direction
    double epsNP; ///< irreversible slip
    double alphaN;  ///< kinematic hardening variable
    double zN;  ///< isotropic hardening variable (for normal direction can be different than the tangential)
    double rN;

    double temp_epsN, temp_damage, temp_epsNP, temp_alphaN, temp_zN, temp_rN; ///<temporary variables
    double temp_stressN, stressN;

    double damageNormal_set_from_the_outside = 0.0;

    double strain_displ_multiplier;

    double prev_damage, prev_zN, prev_epsNP, prev_alphaN, prev_stressN, prev_epsN;

    double energy_PL, energy_D, energy_Kin, energy_Iso, work_tot;
    double temp_Y, prev_Y, Y_next;

    bool symmetric;  ///> if true, symmetric behavior tension/compression is applied

    double Kt;    ///>  initial slope of the softening curve (when Gt - fracture energy - used)

    void print() const;
public:
    DamagePlasticMaterialStatus(DamagePlasticMaterial *m, Element *e, unsigned ipnum);
    virtual ~DamagePlasticMaterialStatus() {};
    void init();
    virtual void update();
    virtual Matrix giveStiffnessTensor(string type, unsigned dim) const;
    virtual Vector giveStress(const Vector &strain, double timeStep);
    virtual Vector giveStressWithFrozenIntVars(const Vector &strain, double timeStep);
    virtual double giveValue(string code) const;
protected:
    void setDamage(const double &new_damage) {
        if ( new_damage > this->damageNormal_set_from_the_outside ) {
            this->damageNormal_set_from_the_outside = new_damage;
        }
    }
};


class DamagePlasticMaterial : public DisMechMaterial
{
private:
    double fc; ///< compressive plastic yielding stress
    double ft;  ///< tensile elastic limit
    double KinN;  ///< isotropic hardening modulus
    double gammaN;  ///< kinematic hardening modulus
    ///< tensile part can be prescribed by Ad or Gt:
    double Ad;  ///< brittlenes of damage evolution
    double Gt;  ///< fracture energy
    double Kt;  ///< initial slope of the softening curve
    double m;  ///< hardening parameter
    bool use_displ;  ///< whether to use absolute values of displacements instead of strains

    bool sym;

public:
    DamagePlasticMaterial() { name = "Damage Plastic material"; };
    ~DamagePlasticMaterial() {};
    void readFromLine(istringstream &iss);
    MaterialStatus *giveNewMaterialStatus(Element *e, unsigned ipnum);
    double giveYieldStress() const { return fc; }
    double giveTensileStrength() const { return ft; }
    double giveElasticLimit() const { return ft / giveE0(); }
    double giveGammaN() const { return gammaN; }
    double giveKinN() const { return KinN; }
    double giveAd() const { return Ad; }
    double giveGt() const { return Gt; }
    double giveKt() const { return Kt; }
    double giveM() const { return m; }
    virtual void init();

    bool useDispl() const { return use_displ; }
    bool isSym() const { return sym; }
};

///////////////////////////////////////////////////////////
// constitutive law for normal and tangential direction together
class FatigueMaterial;
class FatigueMaterialStatus : public FatigueShearMaterialStatus, public DamagePlasticMaterialStatus
{
private:
    double coupled_damage;
public:
    FatigueMaterialStatus(FatigueMaterial *m, Element *e, unsigned ipnum);
    virtual ~FatigueMaterialStatus() {};
    void init();
    virtual void update();
    virtual Matrix giveStiffnessTensor(string type, unsigned dim) const;
    virtual Vector giveStress(const Vector &strain, double timeStep);
    virtual Vector giveStressWithFrozenIntVars(const Vector &strain, double timeStep);
    virtual double giveValue(string code) const;
};


class FatigueMaterial : public FatigueShearMaterial, public DamagePlasticMaterial
{
private:

public:
    FatigueMaterial() { FatigueShearMaterial :: name = "Fatigue material"; DamagePlasticMaterial :: name = "Fatigue material"; };
    ~FatigueMaterial() {};
    void readFromLine(istringstream &iss);
    MaterialStatus *giveNewMaterialStatus(Element *e, unsigned ipnum);
    virtual void init();
};


//////////////////////////////////////////////////////////
// MATERIAL FOR FATIGUE DAMAGE OF CONCRETE
//////////////////////////////////////////////////////////
// ACCORDING TO Damage model for fatigue loading of concrete, ALLICHE, 2004, doi.org/10.1016/j.ijfatigue.2004.02.006
// see also Classification and evaluation of phenomenological numerical models for concrete fatigue behavior under compression, A. Baktheer, R. Chudoba 2019 doi.org/10.1016/j.conbuildmat.2019.06.022

class AllicheMaterial;
class AllicheMaterialStatus : public DisMechMaterialStatus
{
private:
    Point sigma; ///< stress
    Point temp_Y, Y, Y_plus; ///< energy release rate and its positive part
    Point Y_prev;
    Point temp_eps, eps, eps_plus, temp_eps_plus, eps_plus_prev; ///< strain and its positive part
    Point damage, temp_damage; ///< damage

    // Point shear_eps_cur, shear_eps_prev;
    // double temp_damage, damage;

public:
    AllicheMaterialStatus(AllicheMaterial *m, Element *e, unsigned ipnum);
    virtual ~AllicheMaterialStatus() {};
    void init();
    virtual void update();
    virtual Matrix giveStiffnessTensor(string type, unsigned dim) const;
    virtual Vector giveStress(const Vector &strain, double timeStep);
    virtual Vector giveStressWithFrozenIntVars(const Vector &strain, double timeStep);
    // virtual Vector giveStressWrong(const Vector &strain);
    void calculateDamage(const Vector &strain);
    virtual double giveValue(string code) const;
};


class AllicheMaterial : public DisMechMaterial
{
private:
    double lambda, mu;   // Lame parameters
    double C0, C1;
    double K;
    double n;
    double g;  // constant relevant to damage-induced residual stress
    double alphaDam, betaDam;  // parameters of damage degradation moduli
public:
    AllicheMaterial() { name = "Alliche material"; };
    ~AllicheMaterial() {};
    void readFromLine(istringstream &iss);
    MaterialStatus *giveNewMaterialStatus(Element *e, unsigned ipnum);
    double giveLambda() const { return lambda; }
    double giveMu() const { return mu; }
    double giveC0() const { return C0; }
    double giveC1() const { return C1; }
    double giveK() const { return K; }
    double giveN() const { return n; }
    double giveG() const { return g; }
    double giveAlphaDam() const { return alphaDam; }
    double giveBetaDam() const { return betaDam; }
    virtual void init();
};


//////////////////////////////////////////////////////////
// MATERIAL FOR FATIGUE DAMAGE OF CONCRETE
//////////////////////////////////////////////////////////
// ACCORDING TO Continuum damage mechanics for hysteresis and fatigue of quasi-brittle materials and structures
// R. Desmorat et al. 2007 10.1002/nag.532
// see also Classification and evaluation of phenomenological numerical models for concrete fatigue behavior under compression, A. Baktheer, R. Chudoba 2019 doi.org/10.1016/j.conbuildmat.2019.06.022

class DesmoratMaterial;
class DesmoratMaterialStatus : public DisMechMaterialStatus
{
private:
    Point temp_sigma, sigma; ///< stress
    double temp_Y, Y; ///< energy release rate
    Point epsT; ///< strain
    double epsN;
    Point temp_epsPi, epsPi; ///< irreversible strain
    double temp_damage, damage; ///< damage
    double temp_zIso, zIso;
    Point temp_alphaKin, alphaKin;

public:
    DesmoratMaterialStatus(DesmoratMaterial *m, Element *e, unsigned ipnum);
    virtual ~DesmoratMaterialStatus() {};
    void init();
    virtual void update();
    virtual Matrix giveStiffnessTensor(string type, unsigned dim) const;
    virtual Vector giveStress(const Vector &strain, double timeStep);
    virtual Vector giveStressWithFrozenIntVars(const Vector &strain, double timeStep);
    virtual double giveValue(string code) const;
};


class DesmoratMaterial : public DisMechMaterial
{
private:
    double E2;
    double Sigma0;
    double K;
    double gamma;
    double S;
public:
    DesmoratMaterial() { name = "Desmorat material"; };
    ~DesmoratMaterial() {};
    void readFromLine(istringstream &iss);
    MaterialStatus *giveNewMaterialStatus(Element *e, unsigned ipnum);
    double giveE2() const { return E0 * alpha; }
    double giveSigma0() const { return Sigma0; }
    double giveK() const { return K; }
    double giveGamma() const { return gamma; }
    double giveS() const { return S; }
    virtual void init();
};

#endif /* _FATIGUE_MATERIAL_H */
