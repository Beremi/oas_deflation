#ifndef _MAT_RVE_H
#define _MAT_RVE_H

#include "material.h"
#include "globals.h"

class Model; //forward declaraion

//////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////
// TRANSPORT RVE MATERIAL

class TrsprtRVEMaterial;
class TrsprtRVEMaterialStatus : public TrsprtMaterialStatus
{
protected:
    Model *model;

    fs :: path inputfile;
public:
    TrsprtRVEMaterialStatus(TrsprtRVEMaterial *m, Element *e, fs :: path masterfile);
    virtual ~TrsprtRVEMaterialStatus();
    virtual void init();
    virtual Vector giveStress(const Vector &strain);//terminology from mechanics, it returns flux
    virtual double giveEffectiveConductivity(string type) const;
    virtual Matrix giveStiffnessTensor(string type, unsigned dimension) const;
    virtual double giveMassConstant() const;
};

//////////////////////////////////////////////////////////
class TrsprtRVEMaterial : public TrsprtMaterial
{
protected:
    fs :: path inputfile;
public:
    TrsprtRVEMaterial() { name = "transport RVE material";};
    ~TrsprtRVEMaterial() {};
    void readFromLine(istringstream &iss);
    MaterialStatus *giveNewMaterialStatus(Element *e);
};

#endif /* _NODE_C_H */
