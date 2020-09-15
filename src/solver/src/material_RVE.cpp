#include "material_RVE.h"
#include "model.h"
//////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////
// transport RVE status
TrsprtRVEMaterialStatus :: TrsprtRVEMaterialStatus(TrsprtRVEMaterial *m, Element *e, fs :: path masterfile):TrsprtMaterialStatus(m, e){
    name = "transport RVE mat. status";
    inputfile = masterfile;
    model = new Model(false);
}

//////////////////////////////////////////////////////////
TrsprtRVEMaterialStatus :: ~TrsprtRVEMaterialStatus() {
    delete model;
}

//////////////////////////////////////////////////////////
void TrsprtRVEMaterialStatus :: init(){
    model -> readFromFile(inputfile.string());
}

////////////////////////////////////////////////////////// 
Vector TrsprtRVEMaterialStatus :: giveStress(const Vector &strain){
return Vector(0);
}

//////////////////////////////////////////////////////////
double TrsprtRVEMaterialStatus :: giveEffectiveConductivity(string type) const {
return 0;
}

//////////////////////////////////////////////////////////
Matrix TrsprtRVEMaterialStatus :: giveStiffnessTensor(string type, unsigned dimension) const{
return Matrix(0,0);
}

//////////////////////////////////////////////////////////
double TrsprtRVEMaterialStatus :: giveMassConstant() const{
return 0;
};

//////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////
// transport RVE material
void TrsprtRVEMaterial :: readFromLine(istringstream &iss){
    string filename;
    iss >> filename;
    inputfile = GlobPaths :: BASEDIR  / filename;
}

//////////////////////////////////////////////////////////
MaterialStatus* TrsprtRVEMaterial :: giveNewMaterialStatus(Element *e){
    TrsprtRVEMaterialStatus *newstat = new TrsprtRVEMaterialStatus(this, e, inputfile);
    return newstat;
}
