#ifndef _FATIGUE_MATERIAL_H
#define _FATIGUE_MATERIAL_H

#include "material.h"

//////////////////////////////////////////////////////////
// DAMAGE DRIVEN BY CUMMULATIVE SHEAR SLIP

class FatigueShearMaterial;
class FatigueShearMaterialStatus : public DisMechMaterialStatus
{
private:
    double slip; ///< slip
    double damageShear; ///< damage in tangential direction
    double sPi; ///< irreversible slip
    double alphaKin;  ///< kinematic hardening variable
    double zIso;  ///< isotropic hardening variable
    double tang_stiff;  ///< consistent algorithmic (= tangent) shear stiffness 

    double temp_sPi, temp_damageShear, temp_slip, temp_alphaKin, temp_zIso; ///<temporary variables
    
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
    double m;  ///< parameter controling the pressure sensitivity
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


#endif /* _FATIGUE_MATERIAL_H */
