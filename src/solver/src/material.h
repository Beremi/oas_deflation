#ifndef _MATERIAL_H
#define _MATERIAL_H

#include "linear_algebra.h"
#include <vector>
#include <iostream>
#include <fstream>

//////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////
//BASIC MATERIAL

class Element; //forward declaration

class Material;
class MaterialStatus {
private:
    
public:
    MaterialStatus(Material *m, Element *e){name="basic mat. status"; mat = m; element = e;};
    MaterialStatus(Material *m){name="basic mat. status"; mat=m;};
    virtual ~MaterialStatus(){};    
    string whoAmI(){return name;}
    Material* giveMaterial(){return mat;};
    virtual void init(){};
    virtual void update(){};
    virtual double giveValue(string code)const{return 0;};
protected:
    Element *element;
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
    virtual MaterialStatus* giveNewMaterialStatus(Element *e){MaterialStatus* newStatus = new MaterialStatus(this, e); return newStatus;}; 
    string giveName(){return name;}
    virtual void init(){};
protected:
    string name;
};

//////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////
// TRANSPORT MATERIAL

class TrsprtMaterial;
class TrsprtMaterialStatus: public MaterialStatus {
public:
    TrsprtMaterialStatus(TrsprtMaterial *m, Element *e);
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
    MaterialStatus* giveNewMaterialStatus(Element *e);
};

//////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////
// DISCRETE MECHANICAL LINEAR MATERIAL

class DisMechMaterial;
class DisMechMaterialStatus: public MaterialStatus {
protected:

public:
    DisMechMaterialStatus(DisMechMaterial *m, Element *e);
    virtual ~DisMechMaterialStatus(){};  
    Vector giveElasticNormalShearStiffness() const; 
    virtual Vector giveSecantNormalShearStiffness()const{return giveElasticNormalShearStiffness();}   //only elastic
    virtual Vector giveSecantNormalShearStiffness(const Vector &strain)const{return giveElasticNormalShearStiffness();} //only elastic 
    virtual Vector giveStress(const Vector &strain);
    double giveDensity() const;
};

//////////////////////////////////////////////////////////
class DisMechMaterial: public Material {
protected:
    double E0, alpha, density;
    
public:
    DisMechMaterial(){name="discrete mechanical material";};
    ~DisMechMaterial(){};
    double giveDensity(){return density;};
    virtual void readFromLine(istringstream &iss);
    virtual MaterialStatus* giveNewMaterialStatus(Element *e);
    double giveAlpha()const{return alpha;}
    double giveE0()const{return E0;}
};

#endif /* _MATERIAL_H */
