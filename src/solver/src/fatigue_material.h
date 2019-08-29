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


#endif /* _FATIGUE_MATERIAL_H */
