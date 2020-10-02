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
    genereteMicroSources();
    genereteVolumetricAverageBC();
    RVE-> init();
    updateMicroSources();
}

/////////////////////////////////.///////////////////////// 
Vector TrsprtRVEMaterialStatus :: giveStress(const Vector &strain){

    unsigned ndim = RVE->giveDimension();

    //set sources and BC
    for(unsigned i=0; i<ndim; i++){ 
        MSfunctions[i]->setYValue(-strain[i],0);
    }
    volumAverFunc->setYValue(0.,0);

    cout << "strain " << strain[0] << " " << strain[1] << endl;

    //solve
    RVE->resetTime();
    RVE->solve();

    //collect results
    ElementContainer *elems = RVE->giveElements();
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
            lambdaEff += dyadicProduct(n,n)*(e->giveLength()*e->giveArea()*e->giveMatStatus(0)->giveStiffnessTensor("secant",ndim)[0][0]);
            volume += e->giveVolume();  
            stress -= (e->giveArea()*e->giveMatStatus(0)->giveStiffnessTensor("secant",ndim)[0][0])*(n*(DoFs[e->giveNode(1)->giveStartingDoF()]-DoFs[e->giveNode(0)->giveStartingDoF()]));            
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
    unsigned ndim = RVE->giveDimension();
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
            lambdaEff += dyadicProduct(n,n)*(e->giveLength()*e->giveArea()*e->giveMatStatus(0)->giveStiffnessTensor(type,ndim)[0][0]);
            volume += e->giveVolume();            
        }
    }
    return lambdaEff/volume;
}

//////////////////////////////////////////////////////////
double TrsprtRVEMaterialStatus :: giveMassConstant() const{
    ElementContainer *elems = RVE->giveElements();
    Transp1D *e;
    double volume = 0;
    double mass = 0;
    for(unsigned i=0; i<elems->giveSize(); i++){
        e = dynamic_cast<Transp1D*> (elems->giveElement(i));
        if(e){
            mass += e->giveVolume()*e->giveMatStatus(0)->giveMassConstant();
            volume += e->giveVolume();            
        }
    }
    return mass/volume;
};

//////////////////////////////////////////////////////////
void TrsprtRVEMaterialStatus :: genereteVolumetricAverageBC() {

        NodeContainer *nodes = RVE->giveNodes();
        BCContainer *bconds = RVE->giveBC();
        FunctionContainer* funcs = RVE->giveFunctions();
        ElementContainer *elems = RVE->giveElements(); 
        ConstraintContainer *constrs = RVE->giveConstraints(); 

        TrsDoF *tn = new TrsDoF(1);
        nodes->addNode(tn);

        VolumetricAverage *va;
        vector< Node * >vm;
        for ( unsigned n = 0; n < nodes->giveSize(); n++ ) {
            if ( nodes->giveNode(n)->doesTransport() && ( dynamic_cast< TrsDoF * >( nodes->giveNode(n) ) == nullptr ) ) {
                vm.push_back(nodes->giveNode(n) );
            }
        }
        vector< unsigned >dirs;
        dirs.resize(vm.size() );
        va = new VolumetricAverage(vm, dirs, tn, 0, elems, constrs);
        constrs->addConstraint(va);

        unsigned funsize = RVE->giveFunctions()->giveSize();
        vector< double >x, y;
        x.resize(1, 0);
        y.resize(1, 0);
        volumAverFunc = new PieceWiseLinearFunction(x, y);
        RVE->giveFunctions()->addFunction(volumAverFunc);    

        BoundaryCondition *bc;
        vector< int >dBC, nBC;
        dBC.resize(1, -1);
        nBC.resize(1, -1);
        dBC [ 0 ] = funsize;
        bc = new BoundaryCondition(tn, dBC, nBC);
        bconds->addBoundaryCondition(bc);
}

//////////////////////////////////////////////////////////
void TrsprtRVEMaterialStatus :: genereteMicroSources() {

    unsigned funsize = RVE->giveFunctions()->giveSize();
    unsigned ndim = RVE->giveDimension();
    vector< double >x, y;
    x.resize(1, 0);
    y.resize(1, 0);
    MSfunctions.resize(ndim);
    for(unsigned i=0; i<ndim; i++){
        PieceWiseLinearFunction *newf = new PieceWiseLinearFunction(x, y);
        RVE->giveFunctions()->addFunction(newf);    
        MSfunctions[i] = newf;
    }

    BoundaryCondition *bc;
    vector< int >dBC, nBC;
    dBC.resize(ndim, -1);
    nBC.resize(ndim, -1);
    for(unsigned i=0; i<ndim; i++){    
        nBC [ i ] = funsize+i;
    }

    NodeContainer *nodes = RVE->giveNodes();
    for ( unsigned n = 0; n < nodes->giveSize(); n++ ) {
        if ( nodes->giveNode(n)->doesTransport() && ( dynamic_cast< TrsDoF * >( nodes->giveNode(n) ) == nullptr ) ) {
            MSnodes.push_back(nodes->giveNode(n) );
        }
    }
    MSbc.resize(MSnodes.size());
    BCContainer *bconds = RVE->giveBC();
    for ( unsigned n = 0; n < MSnodes.size(); n++ ) {
        bc = new BoundaryCondition(MSnodes[n], dBC, nBC);
        bconds->addBoundaryCondition(bc);
        MSbc[n] = bc;
    }

    MSelems.resize(MSnodes.size());
    MSorder.resize(MSnodes.size());
    unsigned r;
    Transp1D *et;
    ElementContainer *elems = RVE->giveElements();    
    for ( unsigned e = 0; e < elems->giveSize(); e++ ) {
        et = dynamic_cast< Transp1D * >( elems->giveElement(e) );
        if ( et ) {
            for ( unsigned p = 0; p < 2; p++ ) {
                r = ( std :: find(MSnodes.begin(), MSnodes.end(), et->giveNode(p) ) - MSnodes.begin() );
                MSelems[r].push_back(et);
                MSorder[r].push_back(p);
            }
        }
    }
}

//////////////////////////////////////////////////////////
void TrsprtRVEMaterialStatus :: updateMicroSources() {
    unsigned ndim = RVE->giveDimension();
    vector <double > m(ndim);
    Point normal;
    double factor;
    Transp1D * elem;
    for(unsigned i=0; i<MSnodes.size(); i++){
        fill(m.begin(), m.end(), 0);
        for(unsigned e=0; e<MSelems[i].size(); e++){    
            elem = MSelems[i][e];
            normal = elem->giveNormal();
            factor = elem->giveArea()*elem->giveMatStatus(0)->giveStiffnessTensor("secant",ndim)[0][0];
            if(MSorder[i][e]==1) factor *= -1;
            m[0] += normal.getX()*factor;
            m[1] += normal.getY()*factor;
            if(ndim==3) m[2] += normal.getZ()*factor;
        }   
        MSbc[i]->setMultipliers(m);     
    }
}

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
