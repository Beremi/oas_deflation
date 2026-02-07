#include "material_beam.h"
#include "element_container.h"
#include "cross_section.h"
#include "function.h"
#include "model.h"

using namespace std;

//////////////////////////////////////////////////////////
// BEAM MATERIAL STATUS
//////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////

BeamMaterialStatus :: BeamMaterialStatus(BeamMaterial *m, Element *e, unsigned ipnum) : TensMechMaterialStatus(m, e, ipnum) {
    name = "Beam mat. status";
}

//////////////////////////////////////////////////////////
void BeamMaterialStatus :: setCrossSection(CrossSection *X) {
    CS = X;
}

//////////////////////////////////////////////////////////
bool BeamMaterialStatus :: giveValues(string code, Vector &result) const {
    if ( code.compare("normalforce") == 0 ) {
        result.resize(1);
        result [ 0 ] = temp_stress [ 0 ];
        return true;        
    } else if ( code.compare("shearforceY") == 0 ) {
        result.resize(1);
        result [ 0 ] = temp_stress [ 1 ];
        return true;
    } else if ( code.compare("shearforceZ") == 0 ) {
        result.resize(1);
        result [ 0 ] = temp_stress [ 2 ];
        return true;
    } else if ( code.compare("torque") == 0 ) {
        result.resize(1);
        result [ 0 ] = temp_stress [ 3 ];
        return true;
    } else if ( code.compare("bendingmomentY") == 0 ) {
        result.resize(1);
        result [ 0 ] = temp_stress [ 4 ];
        return true;
    } else if ( code.compare("bendingmomentZ") == 0 ) {
        result.resize(1);
        result [ 0 ] = temp_stress [ 5 ];
        return true;
    } else {
        return TensMechMaterialStatus :: giveValues(code, result);
    }
}

//////////////////////////////////////////////////////////
void BeamMaterialStatus :: init() {
}

//////////////////////////////////////////////////////////
void BeamMaterialStatus :: update() {
    TensMechMaterialStatus :: update();
}

//////////////////////////////////////////////////////////
void BeamMaterialStatus :: resetTemporaryVariables() {
    TensMechMaterialStatus :: resetTemporaryVariables();
}


//////////////////////////////////////////////////////////
Matrix BeamMaterialStatus :: giveStiffnessTensor(string type) const {
    ( void ) type;
    Matrix D = Matrix :: Zero(6, 6);
    BeamMaterial *tm = static_cast< BeamMaterial * >( mat );
    D(0, 0) = tm->giveElasticModulus() * CS->giveArea();
    D(1, 1) = CS->giveKappaY() * CS->giveArea() * tm->giveShearModulus();
    D(2, 2) = CS->giveKappaZ() * CS->giveArea() * tm->giveShearModulus();
    D(3, 3) = tm->giveShearModulus() * CS->giveJ();
    D(4, 4) = tm->giveElasticModulus() * CS->giveIy();
    D(5, 5) = tm->giveElasticModulus() * CS->giveIz();
    return D;
}

//////////////////////////////////////////////////////////
void BeamMaterialStatus :: computeStress( double timeStep) {
    //computes internal forces in local reference system (N, Vy, Vz, T, My, Mz)
    BeamMaterialStatus :: computeStressWithFrozenIntVars( timeStep);
}

//////////////////////////////////////////////////////////
void BeamMaterialStatus :: computeStressWithFrozenIntVars( double timeStep) {
    //computes internal forces in local reference system (N, Vy, Vz, T, My, Mz)
    ( void ) timeStep;
    temp_stress = giveStiffnessTensor("elastic") * temp_strain;
}

//////////////////////////////////////////////////////////
void BeamMaterialStatus :: readFromLine(istringstream &iss) {
    TensMechMaterialStatus :: readFromLine(iss);
}

/////////////////////////////////////////////////////////
Matrix BeamMaterialStatus :: giveMassTensor() const {
    //possible issue with rotational inertia, should it be disregarded?
    TensMechMaterial *tm = static_cast< TensMechMaterial * >( mat );
    Matrix c = Matrix :: Zero(6, 6);
    c(0, 0) = c(1, 1) = c(2, 2) = CS->giveArea();
    c(3, 3) = CS->giveJ();
    c(4, 4) = CS->giveIy();
    c(5, 5) = CS->giveIz();
    c *= tm->giveDensity();
    return c;
}

//////////////////////////////////////////////////////////
// BEAM MATERIAL
//////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////

//////////////////////////////////////////////////////////
void BeamMaterial :: readFromLine(istringstream &iss) {
    TensMechMaterial :: readFromLine(iss); //read elastic parameters
};

//////////////////////////////////////////////////////////
MaterialStatus *BeamMaterial :: giveNewMaterialStatus(Element *e, unsigned ipnum) {
    BeamMaterialStatus *newStatus = new BeamMaterialStatus(this, e, ipnum); //needs to be deleted manually
    return newStatus;
};


//////////////////////////////////////////////////////////
void BeamMaterial :: init(MaterialContainer *matcont) {
    TensMechMaterial :: init(matcont);
};

//////////////////////////////////////////////////////////
// BEAM MATERIAL WITH PLASTICITY IN NORMAL DIRECTION STATUS
//////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////

NormalPlasticBeamMaterialStatus :: NormalPlasticBeamMaterialStatus(NormalPlasticBeamMaterial *m, Element *e, unsigned ipnum) : BeamMaterialStatus(m, e, ipnum) {
    name = "Beam mat. w. normal plasticity status";
}

//////////////////////////////////////////////////////////
bool NormalPlasticBeamMaterialStatus :: giveValues(string code, Vector &result) const {
    if ( code.compare("normal_plastic_strain") == 0 ) {
        result.resize(1);
        result [ 0 ] = normalPlasticStrain;
        return true;
    } else {
        return BeamMaterialStatus :: giveValues(code, result);
    }
}

//////////////////////////////////////////////////////////
void NormalPlasticBeamMaterialStatus :: update() {
     normalPlasticStrain = temp_normalPlasticStrain;
     cumulPlasticStrain = temp_cumulPlasticStrain;
     BeamMaterialStatus :: update();
}

//////////////////////////////////////////////////////////
void NormalPlasticBeamMaterialStatus :: resetTemporaryVariables() {
    temp_normalPlasticStrain = normalPlasticStrain;
    temp_cumulPlasticStrain = cumulPlasticStrain;     
    TensMechMaterialStatus :: resetTemporaryVariables();
}

//////////////////////////////////////////////////////////
Matrix NormalPlasticBeamMaterialStatus :: giveStiffnessTensor(string type) const {
    ( void ) type;
    return BeamMaterialStatus :: giveStiffnessTensor(type);
}

//////////////////////////////////////////////////////////
void NormalPlasticBeamMaterialStatus :: computeStress( double timeStep) {
    //computes internal forces in local reference system (N, Vy, Vz, T, My, Mz)    
    computeStressWithFrozenIntVars(timeStep);
    NormalPlasticBeamMaterial *tm = static_cast< NormalPlasticBeamMaterial * >( mat );         
    double E = tm->giveElasticModulus();
    double normstress = temp_stress[0]/CS->giveArea();
    temp_cumulPlasticStrain = cumulPlasticStrain;
    temp_normalPlasticStrain = normalPlasticStrain; 
    double limitValue = tm->givePlasticLimit(temp_cumulPlasticStrain);
    double plstrdiff = (abs(normstress)-limitValue)/E;
    while(plstrdiff>1e-6){
        temp_cumulPlasticStrain += plstrdiff;         
        temp_normalPlasticStrain += copysign(plstrdiff,normstress);        
        normstress = copysign(limitValue,normstress);
        limitValue = tm->givePlasticLimit(temp_cumulPlasticStrain);
        plstrdiff = (abs(normstress)-limitValue)/E;
        
    }
    temp_stress[0] = normstress*CS->giveArea();
}

//////////////////////////////////////////////////////////
void NormalPlasticBeamMaterialStatus :: computeStressWithFrozenIntVars( double timeStep) {
    //computes internal forces in local reference system (N, Vy, Vz, T, My, Mz)
    ( void ) timeStep;    
    Vector shiftedstrain = temp_strain;
    shiftedstrain[0] -= normalPlasticStrain;
    temp_stress = BeamMaterialStatus::giveStiffnessTensor("elastic") * shiftedstrain;
}

//////////////////////////////////////////////////////////
void NormalPlasticBeamMaterialStatus :: readFromLine(istringstream &iss) {
     BeamMaterialStatus :: readFromLine(iss);
}

//////////////////////////////////////////////////////////
// BEAM MATERIAL WITH PLASTICITY IN NORMAL DIRECTION
//////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////

//////////////////////////////////////////////////////////
void NormalPlasticBeamMaterial :: readFromLine(istringstream &iss) {
    BeamMaterial :: readFromLine(iss); //read elastic parameters
    
    iss.clear(); // clear string stream
    iss.seekg(0, iss.beg); //reset position in string stream

    string param;
    bool bplastic_envelope = false;
    while (  iss >> param ) {
        if ( param.compare("plastic_envelope") == 0 ) {
            bplastic_envelope = true;
            iss >> plasticEnvelopeFuncNum;
        }
    }
    if ( !bplastic_envelope ) {
        cerr << name << ": material parameter 'plastic_envelope' was not specified" << endl;
        exit(EXIT_FAILURE);
    }
};

//////////////////////////////////////////////////////////
MaterialStatus *NormalPlasticBeamMaterial :: giveNewMaterialStatus(Element *e, unsigned ipnum) {
    NormalPlasticBeamMaterialStatus *newStatus = new NormalPlasticBeamMaterialStatus(this, e, ipnum); //needs to be deleted manually
    return newStatus;
};


//////////////////////////////////////////////////////////
void NormalPlasticBeamMaterial :: init(MaterialContainer *matcont) {
    BeamMaterial :: init(matcont);
    plasticEnvelope = matcont->giveModel()->giveFunctions()->giveFunction(plasticEnvelopeFuncNum);    
};

//////////////////////////////////////////////////////////
double NormalPlasticBeamMaterial :: givePlasticLimit(double cumstrain) const{
  return plasticEnvelope->giveY(cumstrain);
}
