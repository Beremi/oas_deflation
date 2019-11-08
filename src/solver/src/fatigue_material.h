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
    Point stressT;
    double temp_damageShear, temp_zIso; ///<temporary variables

    double prev_damageShear, prev_zIso;
    Point prev_sPi, prev_alphaKin;

    double Ynext;
    double lambda;
    double temp_lambda;

    double strain_slip_multiplier;
    double regularization_multiplier_area;

    bool checkReturnMap;

    void print() const;
public:
    FatigueShearMaterialStatus(FatigueShearMaterial *m, Element *e);
    virtual ~FatigueShearMaterialStatus() {};
    void init();
    virtual void update();
    virtual Vector giveNormalShearStiffness(string type) const;
    virtual Vector giveStress(const Vector &strain);
    virtual double giveValue(string code) const;
};


class FatigueShearMaterial : public DisMechMaterial
{
private:
    double tauBar; ///< reversibility limit
    double Kin;  ///< isotropic hardening modulus
    double gamma;  ///< kinematic hardening modulus
    double S;  ///< damage strength
    double c, r;  // parameters controling the damage acumullation, c >= 1.0
    double m;  ///< parameter controling the pressure sensitivity, TODO rename to "a" due to coupling with normal direction
    bool use_slip, check_retturn_mapping;
public:
    FatigueShearMaterial() { name = "Fatigue Shear material"; };
    ~FatigueShearMaterial() {};
    void readFromLine(istringstream &iss);
    MaterialStatus *giveNewMaterialStatus(Element *e);
    double giveTauBar() const { return tauBar; }
    double giveKin() const { return Kin; }
    double giveGamma() const { return gamma; }
    double giveS() const { return S; }
    double giveC() const { return c; }
    double giveR() const { return r; }
    double giveM() const { return m; }
    virtual void init();

    bool useSlip() const { return use_slip; }
    bool checkReturnMap() const { return check_retturn_mapping; }
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
    double stressN;

    double strain_displ_multiplier;

    void print() const;
public:
    DamagePlasticMaterialStatus(DamagePlasticMaterial *m, Element *e);
    virtual ~DamagePlasticMaterialStatus() {};
    void init();
    virtual void update();
    virtual Vector giveNormalShearStiffness(string type) const;
    virtual Vector giveStress(const Vector &strain);
    virtual double giveValue(string code) const;
};


class DamagePlasticMaterial : public DisMechMaterial
{
private:
    double fc; ///< compressive plastic yielding stress
    double ft;  ///< tensile elastic limit
    double KinN;  ///< isotropic hardening modulus
    double gammaN;  ///< kinematic hardening modulus
    double Ad;  ///< brittlenes of damage evolution
    double m;  ///< hardening parameter
    bool use_displ;  ///< whether to use absolute values of displacements instead of strains

public:
    DamagePlasticMaterial() { name = "Damage Plastic material"; };
    ~DamagePlasticMaterial() {};
    void readFromLine(istringstream &iss);
    MaterialStatus *giveNewMaterialStatus(Element *e);
    double giveYieldStress() const { return fc; }
    double giveElasticLimit() const { return ft/giveE0(); }
    double giveGammaN() const { return gammaN; }
    double giveKinN() const { return KinN; }
    double giveAd() const { return Ad; }
    double giveM() const { return m; }
    virtual void init();

    bool useDispl() const { return use_displ; }
};

///////////////////////////////////////////////////////////
// constitutive law for normal and tangential direction together
class FatigueMaterial;
class FatigueMaterialStatus : public FatigueShearMaterialStatus, public DamagePlasticMaterialStatus
{
private:

public:
    FatigueMaterialStatus(FatigueMaterial *m, Element *e);
    virtual ~FatigueMaterialStatus() {};
    void init();
    virtual void update();
    virtual Vector giveNormalShearStiffness(string type) const;
    virtual Vector giveStress(const Vector &strain);
    virtual double giveValue(string code) const;
};


class FatigueMaterial : public FatigueShearMaterial, public DamagePlasticMaterial
{
private:

public:
    FatigueMaterial() { FatigueShearMaterial :: name = "Fatigue material"; DamagePlasticMaterial :: name = "Fatigue material"; };
    ~FatigueMaterial() {};
    void readFromLine(istringstream &iss);
    MaterialStatus *giveNewMaterialStatus(Element *e);
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
    AllicheMaterialStatus(AllicheMaterial *m, Element *e);
    virtual ~AllicheMaterialStatus() {};
    void init();
    virtual void update();
    virtual Vector giveNormalShearStiffness(string type) const;
    virtual Vector giveStress(const Vector &strain);
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
    MaterialStatus *giveNewMaterialStatus(Element *e);
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
    DesmoratMaterialStatus(DesmoratMaterial *m, Element *e);
    virtual ~DesmoratMaterialStatus() {};
    void init();
    virtual void update();
    virtual Vector giveNormalShearStiffness(string type) const;
    virtual Vector giveStress(const Vector &strain);
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
    MaterialStatus *giveNewMaterialStatus(Element *e);
    double giveE2() const { return E0 * alpha; }
    double giveSigma0() const { return Sigma0; }
    double giveK() const { return K; }
    double giveGamma() const { return gamma; }
    double giveS() const { return S; }
    virtual void init();
};

#endif /* _FATIGUE_MATERIAL_H */
