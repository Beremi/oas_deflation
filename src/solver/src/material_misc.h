#ifndef _MATERIAL_MISC_H
#define _MATERIAL_MISC_H
// file containing micelanous materials
#include "material_vectorial.h"


//////////////////////////////////////////////////////////
// BRITTLE MATERIAL

class BrittleMaterial;
class BrittleMaterialStatus : public VectMechMaterialStatus
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
    BrittleMaterialStatus(BrittleMaterial *m, Element *e, unsigned ipnum);
    virtual ~BrittleMaterialStatus() {};
    void init();
    virtual void update();
    virtual Matrix giveStiffnessTensor(std :: string type) const;
    virtual Vector giveStress(const Vector &strain, double timeStep);
    virtual Vector giveStressWithFrozenIntVars(const Vector &strain, double timeStep);
    virtual bool giveValues(std :: string code, Vector &result) const;
};


class BrittleMaterial : public VectMechMaterial
    // TODO do only elasto brittle - fully elastic in compression and brittle in tension/shear with possíbility to calc in equiv. space
{
protected:
    double ft, fs;
    // bool compression_recovery;  // NOTE compression recovery true
public:
    BrittleMaterial(unsigned dimension) : VectMechMaterial(dimension)  { name = "Brittle material"; };
    ~BrittleMaterial() {};
    void readFromLine(std :: istringstream &iss);
    MaterialStatus *giveNewMaterialStatus(Element *e, unsigned ipnum);
    double giveFt() { return ft; }
    double giveFs() { return fs; }

    virtual void init(MaterialContainer *matcont);
};

//////////////////////////////////////////////////////////
// CONTACT SHEAR MATERIAL
/*
 * Material (and material status) to mimic contact behavior - no transfer in tension, elastic in compression and if friction prescribed, then resisting in friction according to pressure
 */

class ContactMaterial;
class ContactMaterialStatus : public VectMechMaterialStatus
{
private:
    double temp_normal_strain;
public:
    ContactMaterialStatus(ContactMaterial *m, Element *e, unsigned ipnum);
    virtual ~ContactMaterialStatus() {};
    void init();
    virtual void update();
    virtual Matrix giveStiffnessTensor(std :: string type) const;
    virtual Vector giveStress(const Vector &strain, double timeStep);
    virtual Vector giveStressWithFrozenIntVars(const Vector &strain, double timeStep);
    //virtual bool giveValues(string code, MyVector &result) const;
};


class ContactMaterial : public VectMechMaterial
    // TODO do only elasto brittle - fully elastic in compression and brittle in tension/shear with possíbility to calc in equiv. space
{
private:
    double friction_coef;  // friction coefficient
public:
    ContactMaterial(unsigned dimension) : VectMechMaterial(dimension)  { name = "Contact material"; };
    ~ContactMaterial() {};
    void readFromLine(std :: istringstream &iss);
    MaterialStatus *giveNewMaterialStatus(Element *e, unsigned ipnum);
    double giveFrictionCoef() const { return friction_coef; };

    virtual void init(MaterialContainer *matcont);
};

#endif /* _MATERIAL_MISC_H */
