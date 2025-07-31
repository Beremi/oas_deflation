#include "material_thermomechanical.h"
#include "material_container.h"
#include "material_vectorial.h"

using namespace std;

//////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////
// DISCRETE COUPLED MECHANICS AND HEAT CONDUCTION MATERIAL
ThermoMechanicalMaterialStatus :: ThermoMechanicalMaterialStatus(ThermoMechanicalMaterial *m, Element *e, unsigned ipnum) : CoupledMaterialStatus(m, e, ipnum) {};

//////////////////////////////////////////////////////////
bool ThermoMechanicalMaterialStatus :: giveValues(std :: string code, Vector &result) const {
    return CoupledMaterialStatus :: giveValues(code, result);
}

//////////////////////////////////////////////////////////
Vector ThermoMechanicalMaterialStatus :: giveStress(const Vector &strain, double timeStep) {
    temp_strain = strain;
    Vector mstrain, mstress;

    addTemperatureEffectToMechanics();

    unsigned k = 0;
    unsigned h, i;
    for ( auto &s:stats ) {
        h = s->giveMaterial()->giveStrainSize();
        mstrain.resize(h);
        for ( i = 0; i < h; i++ ) {
            mstrain [ i ] = temp_strain [ k + i ];
        }
        mstress = s->giveStress(mstrain, timeStep);
        for ( i = 0; i < h; i++ ) {
            temp_stress [ k + i ] = mstress [ i ];
        }
        k += h;
    }
    return temp_stress;
}


//////////////////////////////////////////////////////////
void ThermoMechanicalMaterialStatus :: setParameterValue(std :: string code, double value) {
    CoupledMaterialStatus :: setParameterValue(code, value);
    if ( code.compare("temperature") == 0 ) {
        temperature = value;
    }
}

//////////////////////////////////////////////////////////
Vector ThermoMechanicalMaterialStatus :: giveStressWithFrozenIntVars(const Vector &strain, double timeStep) {
    temp_strain = strain;
    Vector mstrain, mstress;

    addTemperatureEffectToMechanics();

    unsigned k = 0;
    unsigned h, i;
    for ( auto &s:stats ) {
        h = s->giveMaterial()->giveStrainSize();
        mstrain.resize(h);
        for ( i = 0; i < h; i++ ) {
            mstrain [ i ] = temp_strain [ k + i ];
        }
        mstress = s->giveStressWithFrozenIntVars(mstrain, timeStep);
        for ( i = 0; i < h; i++ ) {
            temp_stress [ k + i ] = mstress [ i ];
        }
        k += h;
    }
    return temp_stress;
}

//////////////////////////////////////////////////////////
void ThermoMechanicalMaterial :: init(MaterialContainer *matcont) {
    CoupledMaterial :: init(matcont);
    if ( mats.size() != 2 ) {
        cerr << name << "Error: two materials are reqested, " << mats.size() << " provided." << endl;
        exit(1);
    }
    VectMechMaterial *vmm = dynamic_cast< VectMechMaterial * >( mats [ 0 ] );
    VectTrsprtMaterial *vtm = dynamic_cast< VectTrsprtMaterial * >( mats [ 1 ] );
    if ( !vmm || !vtm ) {
        cerr << name << "Error: first material must be inherited from VectMechMaterial, the second from VectTrsprtMaterial." << endl;
        exit(1);
    }
};

//////////////////////////////////////////////////////////
void ThermoMechanicalMaterial :: readFromLine(std :: istringstream &iss) {
    CoupledMaterial :: readFromLine(iss);

    //iss.clear(); // clear string stream
    //iss.seekg(0, iss.beg); //reset position in string stream

    //VectMechMaterial *vmm = static_cast< VectMechMaterial * >( mats [ 0 ] );

    string param;
    bool btec = false;
    bool binitemp = false;
    while (  iss >> param ) {
        if ( param.compare("therm_exp_coeff") == 0 ) {
            btec = true;
            iss >> tec;
        } else if ( param.compare("initial_temperature") == 0 ) {
            binitemp = true;
            iss >> initialTemp;
        }
    }
    if ( !btec ) {
        cerr << name << ": material parameter 'therm_exp_coeff' was not specified" << endl;
        exit(EXIT_FAILURE);
    }
    if ( !binitemp ) {
        cerr << name << ": material parameter 'initial_temperature' was not specified, considered ZERO" << endl;
    }
};

//////////////////////////////////////////////////////////
MaterialStatus *ThermoMechanicalMaterial :: giveNewMaterialStatus(Element *e, unsigned ipnum) {
    ThermoMechanicalMaterialStatus *newStatus = new ThermoMechanicalMaterialStatus(this, e, ipnum);
    return newStatus;
};

//////////////////////////////////////////////////////////
void ThermoMechanicalMaterialStatus ::  addTemperatureEffectToMechanics() {
    ThermoMechanicalMaterial *tmm = static_cast< ThermoMechanicalMaterial * >( mat );
    Vector eig = Vector :: Zero( mat->giveDimension() );
    eig [ 0 ] = ( temperature - tmm->giveInitialTemperature() ) * tmm->giveThermalExpansionCoeff();
    stats [ 0 ]->setEigenStrain(eig);
};
