#ifndef _MATERIAL_H
#define _MATERIAL_H

#include "linalg.h"
#include <vector>
#include <iostream>
#include <fstream>


class Element; //forward declaration
class MaterialContainer; //forward declaration

//////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////
//BASIC MATERIAL

class Material;
class MaterialStatus
{
private:
    std :: vector< MaterialStatus * >matStatComponents;
public:
    MaterialStatus(Material *m, Element *e, unsigned ipnum);
    MaterialStatus(Material *m) { name = "generic mat. status"; mat = m; };
    virtual ~MaterialStatus();
    std :: string giveName() { return name; }
    Material *giveMaterial() { return mat; };
    virtual void init() {};
    virtual void update();
    virtual void resetTemporaryVariables();  ///> if step reset applied, reset temprary variables to last converged state
    virtual bool giveValues(std :: string code, Vector &result) const;
    virtual Vector giveStress(const Vector &strain, double timeStep) { return giveStressWithFrozenIntVars(strain, timeStep); };
    virtual Vector giveStressWithFrozenIntVars(const Vector &strain, double timeStep) { ( void ) strain; ( void ) timeStep; return Vector(0); };
    virtual Vector giveTempStress() const { return temp_stress; };
    virtual Vector giveUpdatedStress() const { return updt_stress; };
    virtual Vector giveTempStrain() const { return temp_strain; };
    virtual Vector giveUpdatedStrain() const { return updt_strain; };
    virtual Matrix giveStiffnessTensor(std :: string type) const { ( void ) type; return Matrix(0, 0); };
    virtual Matrix giveMassTensor() const { return Matrix(0, 0); };
    virtual Matrix giveDampingTensor() const { return Matrix(0, 0); };
    virtual void setEigenStrain(Vector &x);
    //virtual void setID(unsigned i) { idx = i; };
    virtual std :: string giveLineToSave() const { return "no internal variables to export, you need to implement this possibility for " + this->name; }
    virtual void readFromLine(std :: istringstream &iss);
    virtual Vector giveInternalSource()const { return Vector(0); };
    virtual bool isElastic(const bool &now = false) const;
    virtual void setParameterValue(std :: string code, double value) { ( void ) code; ( void ) value; };
    virtual void initializeStressAndStrainVector(unsigned num);
    double giveDissipatedEnergyDensity() const {return dissipEnergyDensity;};  
    virtual double giveEnergyDissipationIncrement() const;
    virtual void computeEnergyDensities();
protected:
    Vector addEigenStrain(const Vector &totalStrain) const;
    Element *element;
    std :: string name;
    Material *mat;
    Vector eigenstrain;
    Vector updt_strain, temp_strain, updt_stress, temp_stress;
    double totalEnergyDensity, strainEnergyDensity, dissipEnergyDensity, updt_dissip_energy;
    unsigned idx;
};

class ElementContainer; //forward declaration

//////////////////////////////////////////////////////////
class Material
{
private:
protected:
    bool produceInternalSources;
    unsigned idx;  // to be able to save element
    unsigned dim;
    unsigned strainsize;
    std :: vector< Material * >matComponents;
    bool dampMatUpdate, massMatUpdate;
    std :: string name;

public:
    Material(unsigned dimension) { name = "basic material"; produceInternalSources = false; dim = dimension; dampMatUpdate = false; massMatUpdate = false; };
    virtual ~Material();
    virtual void readFromLine(std :: istringstream &iss) { ( void ) iss; };
    virtual MaterialStatus *giveNewMaterialStatus(Element *e, unsigned ipnum) { MaterialStatus *newStatus = new MaterialStatus(this, e, ipnum); return newStatus; };
    std :: string whoAmI() { return name; }
    std :: string giveName() { return name; }
    unsigned giveId() { return idx; }
    unsigned giveDimension() const { return dim; }
    unsigned giveStrainSize() const { return strainsize; }
    void setId(const unsigned &i) { this->idx = i; }
    virtual void init(MaterialContainer *matcont) { ( void ) matcont; };
    bool isProducingInternalSources()const { return produceInternalSources; }
    virtual void prepareForStressEvaluation(ElementContainer *elems) { ( void ) elems; };
    bool requiresMassMatrixUpdate() const { return massMatUpdate; };
    bool requiresDampingsMatrixUpdate() const { return dampMatUpdate; };
    virtual bool requestTetrahedralBackgroundMesh()const{return false;} //for volumetric strain
};



//////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////
// COUPLED MATERIAL

class CoupledMaterial;
class CoupledMaterialStatus : public MaterialStatus
{
protected:
    std :: vector< MaterialStatus * >stats;
public:
    CoupledMaterialStatus(Material *m, Element *e, unsigned ipnum);
    virtual ~CoupledMaterialStatus();
    virtual Matrix giveStiffnessTensor(std :: string type) const;
    virtual Vector giveStress(const Vector &strain, double timeStep);
    virtual Vector giveStressWithFrozenIntVars(const Vector &strain, double timeStep);
    virtual bool giveValues(std :: string code, Vector &result) const;
    virtual void update();
    virtual void setParameterValue(std :: string code, double value);
};

//////////////////////////////////////////////////////////
class CoupledMaterial : public Material
{
protected:
    unsigned nmats;
    std :: vector< Material * >mats;
    std :: vector< unsigned >matnums;
public:
    CoupledMaterial(unsigned dimension) : Material(dimension) { name = "generic discrete coupled  material"; };
    ~CoupledMaterial();
    virtual void readFromLine(std :: istringstream &iss);
    virtual void init(MaterialContainer *matcont);
    virtual MaterialStatus *giveNewMaterialStatus(Element *e, unsigned ipnum);
    std :: vector< Material * >giveMaterials() const { return mats; };
    Material *giveMaterial(unsigned i) const;
};



template< typename T >
std :: string to_string_sci(const T a_value)
{
    std :: ostringstream out;
    out << std :: scientific;
    out << a_value;
    return out.str();
}


template< typename T >
std :: string to_string_with_precision(const T a_value, const int n = 6)
{
    std :: ostringstream out;
    out.precision(n);
    out << std :: fixed << a_value;
    return out.str();
}
template< typename T >int sgn(T &val) {
    // NOTE this returns 1 for val = 0 (this is an intention, do not repair it!!)
    return ( T(0) <= val ) - ( val < T(0) );
}


#endif /* _MATERIAL_H */
