#include "material_RVE.h"
#include "model.h"
//////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////
// transport RVE status
TrsprtRVEMaterialStatus :: TrsprtRVEMaterialStatus(TrsprtRVEMaterial *m, Element *e, fs :: path masterfile):TrsprtMaterialStatus(m, e){
    name = "transport RVE mat. status";
    inputfile = masterfile;
    RVE = new Model(false);
}

//////////////////////////////////////////////////////////
TrsprtRVEMaterialStatus :: ~TrsprtRVEMaterialStatus() {
    delete RVE;
}

//////////////////////////////////////////////////////////
void TrsprtRVEMaterialStatus :: init(){
    RVE -> readFromFile(inputfile.string());
    RVE-> init();
}

////////////////////////////////////////////////////////// 
Vector TrsprtRVEMaterialStatus :: giveStress(const Vector &strain){
    RVE->solve();

    ElementContainer *elems = RVE->giveElements();
    unsigned ndim = elems->giveElement(0)->giveDimension();
    Vector stress(ndim);
    Vector DoFs = RVE -> giveSolver() -> giveDoFValues();
    Matrix lambdaEff(ndim,ndim);    
    Transp1D *e;
    double volume = 0;
    Point normal;
    Vector n(ndim);
    for(unsigned i=0; i<elems->giveSize(); i++){
        e = dynamic_cast<Transp1D*> (elems->giveElement(i));
        if(e){
            normal = e->giveNormal();
            n[0] = normal.getX();
            n[1] = normal.getY();
            if(ndim==3) n[2] = normal.getZ();
            lambdaEff += dyadicProduct(n,n)*(e->giveLength()*e->giveArea());
            volume += e->giveVolume();       
            cout << DoFs[e->giveNode(1)->giveStartingDoF()] << " " << DoFs[e->giveNode(0)->giveStartingDoF()] << endl;
            //stress -= n*(e->giveArea()*(DoFs[e->giveNode(1)->giveStartingDoF()]-DoFs[e->giveNode(0)->giveStartingDoF()]));
        }
    }
  
    for(unsigned i=0; i<ndim; i++){
        for(unsigned j=0; j<ndim; j++) stress[i] -= strain[j]*lambdaEff[j][i];        
    }
    
    return stress/volume;
}

//////////////////////////////////////////////////////////
double TrsprtRVEMaterialStatus :: giveEffectiveConductivity(string type) const {
    return 0;
}

//////////////////////////////////////////////////////////
Matrix TrsprtRVEMaterialStatus :: giveStiffnessTensor(string type, unsigned dimension) const{
    ElementContainer *elems = RVE->giveElements();
    unsigned ndim = elems->giveElement(0)->giveDimension();
    Matrix lambdaEff(ndim,ndim);
    Transp1D *e;
    double volume = 0;
    Point normal;
    Vector n(ndim);
    for(unsigned i=0; i<elems->giveSize(); i++){
        e = dynamic_cast<Transp1D*> (elems->giveElement(i));
        if(e){
            normal = e->giveNormal();
            n[0] = normal.getX();
            n[1] = normal.getY();
            if(ndim==3) n[2] = normal.getZ();
            lambdaEff += dyadicProduct(n,n)*(e->giveLength()*e->giveArea());
            volume += e->giveVolume();            
        }
    }
    return lambdaEff/volume;
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
