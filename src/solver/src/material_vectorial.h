#ifndef _MATERIAL_Vect_H
#define _MATERIAL_Vect_H

#include "linalg.h"
#include "material_tensorial.h"
#include <vector>
#include <iostream>
#include <fstream>
#include "material.h"


class Element; //forward declaration

//////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////
// Vect TRANSPORT MATERIAL

class VectTrsprtMaterial;
class VectTrsprtMaterialStatus : public TensTrsprtMaterialStatus
{
protected:
public:
    VectTrsprtMaterialStatus(VectTrsprtMaterial *m, Element *e, unsigned ipnum);
    virtual ~VectTrsprtMaterialStatus() {};
    virtual Matrix giveStiffnessTensor(std :: string type) const;
    virtual MaterialStatus * giveTransportMaterialStatus() { return this; };
};

//////////////////////////////////////////////////////////
class VectTrsprtMaterial : public TensTrsprtMaterial
{
protected:

public:
    VectTrsprtMaterial(unsigned dimension) : TensTrsprtMaterial(dimension) { name = "Vect transport material"; };
    ~VectTrsprtMaterial() {};
    virtual void init(MaterialContainer *matcont) { TensTrsprtMaterial :: init(matcont); strainsize = 1; }
    MaterialStatus *giveNewMaterialStatus(Element *e, unsigned ipnum);
    virtual Material * giveTransportMaterial() { return this; };
};

//////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////
// Vect HEAT CONDUCTION MATERIAL

class VectHeatConductionMaterial;
class VectHeatConductionMaterialStatus : public TensHeatConductionMaterialStatus
{
protected:
    double relativeNormalStrain = 0.;
public:
    VectHeatConductionMaterialStatus(VectHeatConductionMaterial *m, Element *e, unsigned ipnum);
    virtual ~VectHeatConductionMaterialStatus() {};
    virtual Matrix giveStiffnessTensor(std :: string type) const;
    virtual MaterialStatus * giveHeatConductionMaterialStatus() { return this; };
    virtual double updateEffectiveConductivity() const;
    virtual void setParameterValue(std :: string code, double value);  
    virtual Vector giveInternalSource() const;
};

//////////////////////////////////////////////////////////
class VectHeatConductionMaterial : public TensHeatConductionMaterial
{
protected:
    double eta_c = 0.;
    double eta_s = 1.;    
public:
    VectHeatConductionMaterial(unsigned dimension) : TensHeatConductionMaterial(dimension) { name = "Vect heat conduction material"; };
    ~VectHeatConductionMaterial() {};
    void readFromLine(std :: istringstream &iss);    
    virtual void init(MaterialContainer *matcont) { TensHeatConductionMaterial :: init(matcont); strainsize = 1; }
    MaterialStatus *giveNewMaterialStatus(Element *e, unsigned ipnum);
    virtual Material * giveHeatConductionMaterial() { return this; };
    double giveEtaS() const {return eta_s;};
    double giveEtaC() const {return eta_c;};    
};


//////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////
// Vect TRANSPORT MATERIAL WITH COUPLING TERMS

class VectTrsprtCoupledMaterial;
class VectTrsprtCoupledMaterialStatus : public VectTrsprtMaterialStatus
{
protected:
    double temp_volumetricStrain, volumetricStrain, volStrainRate, crackParam, temp_crackVolume, crackVolume, crackVolumeRate, pressure, pressureRate;
    double temperature;
public:
    VectTrsprtCoupledMaterialStatus(VectTrsprtCoupledMaterial *m, Element *e, unsigned ipnum);
    virtual ~VectTrsprtCoupledMaterialStatus() {};
    virtual void computeStress(double timeStep);
    virtual void computeStressWithFrozenIntVars(double timeStep);
    virtual double giveEffectiveConductivity(std :: string type) const;
    virtual void update();
    virtual void resetTemporaryVariables();
    virtual void init() { volumetricStrain = 0; };
    virtual Vector giveInternalSource() const;
    virtual bool giveValues(std :: string code, Vector &result) const;
    virtual void setParameterValue(std :: string code, double value);
    virtual double updateEffectiveConductivity() const;
    void updateRateVariables(double timeStep);
};

//////////////////////////////////////////////////////////
class VectTrsprtCoupledMaterial : public VectTrsprtMaterial
{
private:
    double crack_turtuosity, biotCoeff, refP, Kw;
public:
    VectTrsprtCoupledMaterial(unsigned dimension) : VectTrsprtMaterial(dimension) { name = "coupled transport material";  produceInternalSources = true; refP = 0.; Kw = 2.15e9; };
    ~VectTrsprtCoupledMaterial() {};
    void readFromLine(std :: istringstream &iss);
    MaterialStatus *giveNewMaterialStatus(Element *e, unsigned ipnum);
    double giveTurtuosity() { return crack_turtuosity; };
    double giveBiotCoeff() const { return biotCoeff; };
    double giveKw() const { return Kw; };
    double giveReferencePressure() const { return refP; };
};

//////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////
// Vect MECHANICAL LINEAR MATERIAL

class VectMechMaterial;
class VectMechMaterialStatus : public MaterialStatus
{
protected:
    double normalEnergyDensity, shearEnergyDensity;
    double eigenVolumetricStrain;
    double temp_volumetricStrain, volumetricStrain;
    double temp_volumetricStrain_total, volumetricStrain_total;

    virtual double addEigenVolumetricStrain(double x) const { return x + eigenVolumetricStrain; };
public:
    VectMechMaterialStatus(VectMechMaterial *m, Element *e, unsigned ipnum);
    virtual ~VectMechMaterialStatus() {};
    virtual Matrix giveStiffnessTensor(std :: string type) const;
    virtual void computeStress(double timeStep);
    virtual void computeStressWithFrozenIntVars(double timeStep);
    double giveDensity() const;
    virtual bool isElastic(const bool &now = false) const { ( void ) now; return true; };
    virtual bool giveValues(std :: string code, Vector &result) const;
    virtual double giveNormalCrackOpening()const { return 0; };
    virtual void update();
    virtual MaterialStatus * giveMechanicalMaterialStatus() { return this; };
    virtual void removeEigenStrain() { MaterialStatus :: removeEigenStrain(); eigenVolumetricStrain = 0; };
    virtual void addToEigenVolumetricStrain(double x);
    virtual void setParameterValue(std :: string code, double value);
};

//////////////////////////////////////////////////////////
class VectMechMaterial : public Material
{
protected:
    double E0; //[Pa]
    double alpha; //[-]
    double density; //[kg/m3] 
public:
    VectMechMaterial(unsigned dimension);
    ~VectMechMaterial() {};
    double giveDensity() { return density; };
    virtual void readFromLine(std :: istringstream &iss);
    virtual MaterialStatus *giveNewMaterialStatus(Element *e, unsigned ipnum);
    double giveAlpha() const { return alpha; }
    double giveE0() const { return E0; }
    virtual void init(MaterialContainer *matcont) { Material :: init(matcont); strainsize = dim; }
    virtual Material * giveMechanicalMaterial() { return this; };
};

//////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////
// Vect MECHANICAL LINEAR MATERIAL WITH VOLUMETRIC-DEVIATORIC SPLIT

class VectMechVolDevSplitMaterial;
class VectMechVolDevSplitMaterialStatus : public VectMechMaterialStatus
{
protected:
    double temp_volumetricStrain;
public:
    VectMechVolDevSplitMaterialStatus(VectMechVolDevSplitMaterial *m, Element *e, unsigned ipnum);
    virtual ~VectMechVolDevSplitMaterialStatus() {};
    virtual Matrix giveStiffnessTensor(std :: string type) const;
    virtual void computeStressWithFrozenIntVars(double timeStep);
    double giveDensity() const;
    virtual bool isElastic(const bool &now = false) const { ( void ) now; return true; };
    virtual bool giveValues(std :: string code, Vector &result) const;
    virtual void setParameterValue(std :: string code, double value);
};

//////////////////////////////////////////////////////////
class VectMechVolDevSplitMaterial : public VectMechMaterial
{
protected:
public:
    VectMechVolDevSplitMaterial(unsigned dimension);
    ~VectMechVolDevSplitMaterial() {};
    virtual void readFromLine(std :: istringstream &iss);
    virtual MaterialStatus *giveNewMaterialStatus(Element *e, unsigned ipnum);
    virtual void init(MaterialContainer *matcont) { Material :: init(matcont); strainsize = dim; }
};

//////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////
// Vect MECHANICAL LINEAR MATERIAL with ROTATIONAL STIFFNESS

class VectMechMaterialWithRotationalStiffness;
class VectMechMaterialWithRotationalStiffnessStatus : public VectMechMaterialStatus
{
protected:
public:
    VectMechMaterialWithRotationalStiffnessStatus(VectMechMaterialWithRotationalStiffness *m, Element *e, unsigned ipnum);
    virtual ~VectMechMaterialWithRotationalStiffnessStatus() {};
    virtual Matrix giveStiffnessTensor(std :: string type) const;
    virtual void computeStressWithFrozenIntVars(double timeStep);
    double giveDensity() const;
    virtual bool isElastic(const bool &now = false) const { ( void ) now; return true; };
    virtual bool giveValues(std :: string code, Vector &result) const;
};

//////////////////////////////////////////////////////////
class VectMechMaterialWithRotationalStiffness : public VectMechMaterial
{
protected:
    double beta;
public:
    VectMechMaterialWithRotationalStiffness(unsigned dimension) : VectMechMaterial(dimension) { name = "Vect mechanical material with rotational stiffness"; };
    ~VectMechMaterialWithRotationalStiffness() {};
    virtual void readFromLine(std :: istringstream &iss);
    virtual MaterialStatus *giveNewMaterialStatus(Element *e, unsigned ipnum);
    double giveBeta() const { return beta; }
    virtual void init(MaterialContainer *matcont) { Material :: init(matcont); strainsize = dim + ( ( dim == 3 ) ? 3 : 1 ); }
};

#endif /* _MATERIAL_Vect_H */
