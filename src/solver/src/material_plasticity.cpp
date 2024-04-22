#include "material_plasticity.h"


using namespace std;


//////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////
// VON MISES PLASTICITY
//////////////////////////////////////////////////////////
VonMisesPlasticMaterialStatus :: VonMisesPlasticMaterialStatus(VonMisesPlasticMaterial *m, Element *e, unsigned ipnum) : TensMechMaterialStatus(m, e, ipnum) {
    name = "tensorial mechanical mat. status";
    unsigned strainsize = m->giveStrainSize();
    temp_plasticstrain = Vector::Zero(strainsize);
    plasticstrain = Vector::Zero(strainsize);
    temp_backstress = Vector::Zero(strainsize);
    backstress = Vector::Zero(strainsize);
    sigmay = 0;
    temp_sigmay = 0;
}

//////////////////////////////////////////////////////////
Vector VonMisesPlasticMaterialStatus :: giveStress(const Vector &strain, double timeStep) {
    temp_strain = addEigenStrain(strain);
    //HERE COMES THE CONSTITUTIVE ROUTINES

    VonMisesPlasticMaterial *vmpm = static_cast<VonMisesPlasticMaterial*>(mat);
    vmpm->giveHardeningModulus();



    temp_stress = VonMisesPlasticMaterialStatus :: giveStressWithFrozenIntVars(strain, timeStep);
    return temp_stress;
};

//////////////////////////////////////////////////////////
Vector VonMisesPlasticMaterialStatus :: giveStressWithFrozenIntVars(const Vector &strain, double timeStep) {
    ( void ) timeStep;
    temp_strain = addEigenStrain(strain);
    temp_stress = giveStiffnessTensor("elastic") * (temp_strain - temp_plasticstrain);
    return temp_stress;
};

//////////////////////////////////////////////////////////
Matrix VonMisesPlasticMaterialStatus :: giveStiffnessTensor(string type) const {
    ( void ) type;
    return TensMechMaterialStatus::giveStiffnessTensor(type);
};

//////////////////////////////////////////////////////////
void VonMisesPlasticMaterialStatus :: update() {
    TensMechMaterialStatus :: update();
    backstress = temp_backstress;
    plasticstrain = temp_plasticstrain;
    sigmay = temp_sigmay;
}

//////////////////////////////////////////////////////////
bool VonMisesPlasticMaterialStatus :: giveValues(string code, Vector &result) const {
    if ( code.compare("plastic_strain") == 0 || code.compare("plasticstrain") == 0 ) {
        unsigned size = temp_plasticstrain.size();
        result.resize(size);
        for ( unsigned p = 0; p < size; p++ ) {
            result [ p ] = temp_plasticstrain [ p ];
        }
        return true;
    } else if ( code.compare("backstress") == 0 || code.compare("back_stress") == 0 ) {
        unsigned size = temp_backstress.size();
        result.resize(size);
        for ( unsigned p = 0; p < size; p++ ) {
            result [ p ] = temp_backstress [ p ];
        }
        return true;
    } else if ( code.compare("sigmay") == 0 ) {
        result.resize(1);
        result [ 0 ] = temp_sigmay;
        return true;
    } else {
        return TensMechMaterialStatus :: giveValues(code, result);
    }
}

//////////////////////////////////////////////////////////
VonMisesPlasticMaterial::VonMisesPlasticMaterial(unsigned dimension) : TensMechMaterial(dimension) { 
    name = "von Mises plastic material"; 
    strainsize = ( dim - 1 ) * 3; 
    sigma0 = 0;
    H = 0;
    beta = 0;
};

//////////////////////////////////////////////////////////
void VonMisesPlasticMaterial :: readFromLine(istringstream &iss) {

    TensMechMaterial :: readFromLine(iss); //read elastic parameters

    iss.clear(); // clear string stream
    iss.seekg(0, iss.beg); //reset position in string stream

    string param;
    bool bs0 = false;

    while (  iss >> param ) {
        if ( param.compare("sigma0") == 0 ) {
            bs0 = true;
            iss >> sigma0;
        } else if ( param.compare("H") == 0 ) {
            iss >> H;
        } else if ( param.compare("beta") == 0 ) {
            iss >> beta;
        }
    }
    if ( !bs0 ) {
        cerr << name << ": material parameter 'sigma0' was not specified" << endl;
        exit(EXIT_FAILURE);
    }
};

//////////////////////////////////////////////////////////
MaterialStatus *VonMisesPlasticMaterial :: giveNewMaterialStatus(Element *e, unsigned ipnum) {
    VonMisesPlasticMaterialStatus *newStatus = new VonMisesPlasticMaterialStatus(this, e, ipnum);
    return newStatus;
};

