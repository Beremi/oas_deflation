#ifndef _MATERIAL_PLASTICITY_H
#define _MATERIAL_PLASTICITY_H

#include "material_tensorial.h"

//////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////
// VON MISES PLASTICITY

class VonMisesPlasticMaterial;
class VonMisesPlasticMaterialStatus : public TensMechMaterialStatus
{
protected:
    Vector backstress, plasticstrain;
    Vector temp_backstress, temp_plasticstrain;
    double sigmay, temp_sigmay;

public:
    VonMisesPlasticMaterialStatus(VonMisesPlasticMaterial *m, Element *e, unsigned ipnum);
    virtual ~VonMisesPlasticMaterialStatus() {};
    virtual Matrix giveStiffnessTensor(std :: string type) const;
    virtual Vector giveStress(const Vector &strain, double timeStep);
    virtual Vector giveStressWithFrozenIntVars(const Vector &strain, double timeStep);
    virtual bool giveValues(std :: string code, Vector &result) const;
    virtual void update();
};

//////////////////////////////////////////////////////////
class VonMisesPlasticMaterial : public TensMechMaterial
{
protected:
    double H, beta, sigma0;    

public:
    VonMisesPlasticMaterial(unsigned dimension);
    ~VonMisesPlasticMaterial() {};
    virtual void readFromLine(std :: istringstream &iss);
    virtual MaterialStatus *giveNewMaterialStatus(Element *e, unsigned ipnum);
    double giveHardeningModulus() const { return H; }
    double giveBetaRatio() const { return beta; }
    double giveSigma0() const { return sigma0; };
    virtual void init(MaterialContainer *matcont) { Material :: init(matcont);  }
};

#endif /* _MATERIAL_PLASTICITY_H */
