#ifndef _FATIGUE_MATERIAL_H
#define _FATIGUE_MATERIAL_H

#include "material.h"

//////////////////////////////////////////////////////////
// DAMAGE DRIVEN BY CUMMULATIVE SHEAR SLIP

class FatigueShearMaterial;
class FatigueShearMaterialStatus : public DisMechMaterialStatus
{
private:
    double slip_cur, slip_prev;
    double maxEpsT, temp_maxEpsT;
    double damageShear, temp_damageShear;
    double temp_stiffMultip, stiffMultip;
    double sPi, temp_sPi;  // cumulative slip
    // double tauPi, temp_tauPi ;
    double tauPiTrial, tauTildaPiTrial;  //  sliding bond stress
    double alphaKin, temp_alphaKin;  // kinematic hardening variable
    double zIso, temp_zIso;  // isotropic hardening variable
    double temp_EAlg, EAlg;  // tangential shear stiffness
    double RAND_H;

    double computeTrials(const double &stressN, const FatigueShearMaterial *m);  // returns f_trial (negative f_trial -> elastic step)
    void computeShearStifness(const Vector &strain);
    void print() const;
public:
    FatigueShearMaterialStatus(FatigueShearMaterial *m, Element *e);
    virtual ~FatigueShearMaterialStatus() {};
    void init();
    virtual void update();
    virtual Vector giveSecantNormalShearStiffness() const;
    virtual Vector giveSecantNormalShearStiffness(const Vector &strain);
    virtual Vector giveTangentNormalShearStiffness() const;
    virtual Vector giveTangentNormalShearStiffness(const Vector &strain);
    virtual Vector giveStress(const Vector &strain);
    virtual double giveValue(string code) const;
};


class FatigueShearMaterial : public DisMechMaterial
{
private:
    // double Eb;  // bond stiffness (~shear stifness? Eb = E_0*alpha...?)
    double tauBar; // reversibility limit
    double Kin;  // isotropic hardening modulus
    double gamma;  // kinematic hardening modulus
    double S;  // damage strength
    double c, r;  // parameters controling the damage acumullation, c >= 1.0
    double m;  // parameter controling the pressure sensitivity
public:
    FatigueShearMaterial() { name = "Fatigue Shear material"; };
    ~FatigueShearMaterial() {};
    void readFromLine(istringstream &iss);
    MaterialStatus *giveNewMaterialStatus(Element *e);
    double giveEb() const { return alpha * E0; }
    double giveTauBar() const { return tauBar; }
    double giveKin() const { return Kin; }
    double giveGamma() const { return gamma; }
    double giveS() const { return S; }
    double giveC() const { return c; }
    double giveR() const { return r; }
    double giveM() const { return m; }
    virtual void init();
};


#endif /* _FATIGUE_MATERIAL_H */
