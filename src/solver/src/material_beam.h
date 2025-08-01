#ifndef _BEAM_MATERIAL_H
#define _BEAM_MATERIAL_H

#include "material_tensorial.h"

class CrossSection; //forward declaration

//////////////////////////////////////////////////////////
// BEAM MATERIAL

class BeamMaterial;
class BeamMaterialStatus : public TensMechMaterialStatus
{
protected:
    CrossSection* CS;    
public:
    BeamMaterialStatus(BeamMaterial *m, Element *e, unsigned ipnum);
    virtual ~BeamMaterialStatus() {};
    virtual void init();
    virtual void update();
    void setCrossSection(CrossSection* X);
    CrossSection* giveCrossSection() { return CS; }
    virtual void resetTemporaryVariables();
    virtual Matrix giveStiffnessTensor(std :: string type) const;
    virtual Vector giveStress(const Vector &strain, double timeStep);
    virtual Vector giveStressWithFrozenIntVars(const Vector &strain, double timeStep);
    virtual bool giveValues(std :: string code, Vector &result) const;
    virtual void readFromLine(std :: istringstream &iss);
    virtual Matrix giveMassTensor() const ;
};

//////////////////////////////////////////////////////////
class BeamMaterial : public TensMechMaterial
{
protected:

public:
    BeamMaterial() : TensMechMaterial(3) { name = "Beam material";};
    virtual ~BeamMaterial() {};
    virtual void readFromLine(std :: istringstream &iss);
    virtual MaterialStatus *giveNewMaterialStatus(Element *e, unsigned ipnum);
    virtual void init(MaterialContainer *matcont);
};


#endif /* _BEAM_MATERIAL_H */
