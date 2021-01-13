#ifndef _MATERIAL_MISC_H
#define _MATERIAL_MISC_H
// file containing micelanous materials
#include "material.h"


//////////////////////////////////////////////////////////
// BRITTLE MATERIAL

class BrittleMaterial;
class BrittleMaterialStatus : public DisMechMaterialStatus
{
protected:
    /*
     * comp_recovery: full stiffness if loading goes back to negative values after previously reaching tensile strength
     * eqiv: calculate in equivalent space (normal + shear transformed into equivalent single - normal - axis)
     */
    bool damage, temp_damage, comp_recovery;
    double RAND_H;
    double L;
    double temp_crackOpening, temp_normal_strain;
private:
    void computeDamage(const Vector &strain);
public:
    BrittleMaterialStatus(BrittleMaterial *m, Element *e);
    virtual ~BrittleMaterialStatus() {};
    void init();
    virtual void update();
    virtual Matrix giveStiffnessTensor(string type, unsigned dim) const;
    virtual Vector giveStress(const Vector &strain);
    virtual Vector giveStressWithFrozenIntVars(const Vector &strain) const;
    virtual double giveValue(string code) const;
};


class BrittleMaterial : public DisMechMaterial
    // TODO do only elasto brittle - fully elastic in compression and brittle in tension/shear with possíbility to calc in equiv. space
{
protected:
    double ft, fs;
    // bool compression_recovery;  // NOTE compression recovery true
public:
    BrittleMaterial() { name = "Brittle material"; };
    ~BrittleMaterial() {};
    void readFromLine(istringstream &iss);
    MaterialStatus *giveNewMaterialStatus(Element *e);
    double giveFt() { return ft; }
    double giveFs() { return fs; }

    virtual void init();
};

// //////////////////////////////////////////////////////////
// // PLASTO - BRITTLE MATERIAL
//
// class PlastoBrittleMaterial;
// class PlastoBrittleMaterialStatus : public BrittleMaterialStatus
// {
// private:
//     double plast_strain;
//     double RAND_H;
//     double temp_crackOpening;
//
//     void calcPlastStrain();
// public:
//     PlastoBrittleMaterialStatus(PlastoBrittleMaterial *m, Element *e);
//     virtual ~PlastoBrittleMaterialStatus() {};
//     void init();
//     virtual void update();
//     virtual Vector giveNormalShearStiffness(string type) const;
//     virtual Vector giveStress(const Vector &strain);
//     virtual Vector giveStressWithFrozenIntVars(const Vector &strain) const;
//     virtual double giveValue(string code) const;
// };
//
//
// class PlastoBrittleMaterial : public BrittleMaterial
//   // TODO do plasto brittle - can be useful in any way? (well, can be implemented in future)
//   // plastic strain in compression and brittle in tension/shear with possíbility to calc in equiv. space
// {
// private:
//     double fc;
// public:
//     PlastoBrittleMaterial() { name = "PlastoBrittle material"; };
//     ~PlastoBrittleMaterial() {};
//     void readFromLine(istringstream &iss);
//     MaterialStatus *giveNewMaterialStatus(Element *e);
//     double giveFc() { return fc; }
//
//     virtual void init();
// };

//////////////////////////////////////////////////////////
// CONTACT SHEAR MATERIAL
/*
 * Material (and material status) to mimic contact behavior - no transfer in tension, elastic in compression and if friction prescribed, then resisting in friction according to pressure
 */

class ContactMaterial;
class ContactMaterialStatus : public DisMechMaterialStatus
{
private:
    double temp_normal_strain;
public:
    ContactMaterialStatus(ContactMaterial *m, Element *e);
    virtual ~ContactMaterialStatus() {};
    void init();
    virtual void update();
    virtual Matrix giveStiffnessTensor(string type, unsigned dim) const;
    virtual Vector giveStress(const Vector &strain);
    virtual Vector giveStressWithFrozenIntVars(const Vector &strain) const;
    // virtual double giveValue(string code) const;
};


class ContactMaterial : public DisMechMaterial
    // TODO do only elasto brittle - fully elastic in compression and brittle in tension/shear with possíbility to calc in equiv. space
{
private:
    double friction_coef;  // friction coefficient
public:
    ContactMaterial() { name = "Contact material"; };
    ~ContactMaterial() {};
    void readFromLine(istringstream &iss);
    MaterialStatus *giveNewMaterialStatus(Element *e);
    double giveFrictionCoef() const { return friction_coef; };

    virtual void init();
};

#endif
