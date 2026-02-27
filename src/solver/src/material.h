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
    virtual void computeStress(double timeStep) { computeStressWithFrozenIntVars(timeStep); };
    virtual void computeStressWithFrozenIntVars(double timeStep) { ( void ) timeStep; };
    virtual Vector giveTempStress() const { return temp_stress; };
    virtual Vector giveUpdatedStress() const { return updt_stress; };
    virtual void setTotalTempStrain(Vector str) {  temp_strain_total = str; temp_strain = str; };
    virtual Vector giveTempStrain() const { return temp_strain; };
    virtual Vector giveUpdatedStrain() const { return updt_strain; };
    virtual Vector giveTotalTempStrain() const { return temp_strain_total; };
    virtual Vector giveTotalUpdatedStrain() const { return updt_strain_total; };
    virtual Matrix giveStiffnessTensor(std :: string type) const { ( void ) type; return Matrix(0, 0); };
    virtual Matrix giveMassTensor() const { return Matrix(0, 0); };
    virtual Matrix giveDampingTensor() const { return Matrix(0, 0); };
    virtual void addToEigenStrain(const Vector &x);
    virtual void addToEigenVolumetricStrain(double x) { ( void ) x; };
    virtual void removeEigenStrain() { eigenstrain.setZero(); };
    Vector giveEigenStrain() const {return eigenstrain;};    
    //virtual void setID(unsigned i) { idx = i; };
    virtual std :: string giveLineToSave() const { return "no internal variables to export, you need to implement this possibility for " + this->name; }
    virtual void readFromLine(std :: istringstream &iss);
    virtual Vector giveInternalSource()const { return Vector(0); };
    virtual bool isElastic(const bool &now = false) const;
    virtual void setParameterValue(std :: string code, double value) { ( void ) code; ( void ) value; };
    virtual void initializeStressAndStrainVector();
    double giveDissipatedEnergyDensity() const { return dissipEnergyDensity; };
    virtual double giveEnergyDissipationIncrement() const;
    virtual void computeEnergyDensities();
    virtual MaterialStatus * giveMechanicalMaterialStatus() { return nullptr; };
    virtual MaterialStatus * giveTransportMaterialStatus() { return nullptr; };
    virtual MaterialStatus * giveHeatConductionMaterialStatus() { return nullptr; };
protected:
    virtual Vector addEigenStrain(const Vector &totalStrain) const;
    virtual double addEigenVolumetricStrain(double x) const { return x; };
    Element *element;
    std :: string name;
    Material *mat;
    Vector eigenstrain;
    Vector updt_strain, temp_strain, updt_stress, temp_stress;
    Vector updt_strain_total, temp_strain_total;
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
    virtual bool requestTetrahedralBackgroundMesh()const { return false; } //for volumetric strain
    virtual Material * giveMechanicalMaterial() { return nullptr; };
    virtual Material * giveTransportMaterial() { return nullptr; };
    virtual Material * giveHeatConductionMaterial() { return nullptr; };
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
    virtual void init();
    virtual Matrix giveStiffnessTensor(std :: string type) const;
    virtual void computeStress(double timeStep);
    virtual void computeStressWithFrozenIntVars(double timeStep);
    virtual void setTotalTempStrain(Vector str);
    virtual bool giveValues(std :: string code, Vector &result) const;
    virtual void update();
    virtual void setParameterValue(std :: string code, double value);
    virtual MaterialStatus * giveMechanicalMaterialStatus();
    virtual MaterialStatus * giveTransportMaterialStatus();
    virtual MaterialStatus * giveHeatConductionMaterialStatus();
    virtual void initializeStressAndStrainVector();
    virtual void addToEigenVolumetricStrain(double x);
    virtual void addToEigenStrain(const Vector &x);
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
    virtual Material * giveMechanicalMaterial();
    virtual Material * giveTransportMaterial();
    virtual Material * giveHeatConductionMaterial();
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
