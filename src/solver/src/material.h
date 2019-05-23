#ifndef _MATERIAL_H
#define _MATERIAL_H

#include "linear_algebra.h"
#include <vector>
#include <iostream>
#include <fstream>

//////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////
//BASIC MATERIAL


class Material;
class MaterialStatus {
private:
    
public:
    MaterialStatus(){name="basic mat. status";};
    MaterialStatus(Material *m){name="basic mat. status"; mat=m;};
    virtual ~MaterialStatus(){};    
    string whoAmI(){return name;}
    Material* giveMaterial(){return mat;};
protected:
    string name;
    Material *mat;
};


//////////////////////////////////////////////////////////
class Material {
private:
    
public:
    Material(){name="basic material";};
    virtual ~Material(){};
    virtual void readFromLine(istringstream &iss){}; 
    virtual MaterialStatus* giveNewMaterialStatus(){MaterialStatus* newStatus = new MaterialStatus(this); return newStatus;}; 
    string giveName(){return name;}
protected:
    string name;
};

//////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////
// TRANSPORT MATERIAL

class TrsprtMaterial;
class TrsprtMaterialStatus: public MaterialStatus {
public:
    TrsprtMaterialStatus(TrsprtMaterial *m);
    virtual ~TrsprtMaterialStatus(){};  
    double giveConductivity() const; 
    double giveCapacity() const;
};

//////////////////////////////////////////////////////////
class TrsprtMaterial: public Material {
private:
    double conductivity, capacity;
public:
    TrsprtMaterial(){name="transport material";};
    ~TrsprtMaterial(){};
    double giveConductivity(){return conductivity;}; 
    double giveCapacity(){return capacity;};
    void readFromLine(istringstream &iss);
    MaterialStatus* giveNewMaterialStatus();
};

//////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////
// DISCRETE MECHANICAL LINEAR MATERIAL

class DisMechMaterial;
class DisMechMaterialStatus: public MaterialStatus {
private:
    DisMechMaterialStatus(){};
public:
    DisMechMaterialStatus(DisMechMaterial *m);
    virtual ~DisMechMaterialStatus(){};  
    vector<double> giveNormalShearStiffness() const; 
    double giveDensity() const;
};

//////////////////////////////////////////////////////////
class DisMechMaterial: public Material {
private:
    double E0, alpha, density;
    
public:
    DisMechMaterial(){name="discrete mechanical material";};
    ~DisMechMaterial(){};
    virtual vector<double> giveNormalShearStiffness(){vector<double> c; c.resize(2); c[0] = E0; c[1] = alpha*E0; return c;} 
    double giveDensity(){return density;};
    virtual void readFromLine(istringstream &iss);
    virtual MaterialStatus* giveNewMaterialStatus();
};

//////////////////////////////////////////////////////////
// CUSATIS/MARS MATERIAL
/*
class MarsMaterial;
class MarsMaterialStatus: public DisMechMaterial {
private:
    double epsEQ, epsT, epsN, epsmaxEQ, epsmaxT, epsmaxN;
    double omega, S0, chi, K0, K1;
    double damage, temp_damage;
    double RAND_H;
public:
    MarsMaterialStatus(MarsMaterialMaterial *m);
    virtual ~MarsMaterialStatus(){};  
};


class MarsMaterial: public Material {
private:

public:
    MarsMaterial(){name="discrete mechanical material";};
    ~MarsMaterial(){};
    vector<double> giveNormalShearStiffness(){vector<double> c; c.resize(2); c[0] = E0*(1-temp_damage); c[1] = alpha*E0*(1-temp_damage); return c;} 
    void readFromLine(istringstream &iss);
    MaterialStatus* giveNewMaterialStatus();
};
*/

/*
//////////////////////////////////////////////////////////
//COUPLED DISCRETE MECHANICAL - TRANSPORT MATERIAL

class DisMechLinMat;
class DisMechLinMatStatus: public MaterialStatus {
private:

public:    
    DisMechLinMatStatus();
    DisMechLinMatStatus(DisMechLinMat *m);
    ~DisMechLinMatStatus(){};
};

class DisMechLinMat: public Material {
private:
    double E0, alpha, density;
public:    
    DisMechLinMat();
    ~DisMechLinMat(){};
    void readFromLine(istringstream &iss);
    vector<double> giveStress(vector<double> strain) const;
    MaterialStatus* giveNewMaterialStatus(); 
};
*/

#endif /* _MATERIAL_H */
