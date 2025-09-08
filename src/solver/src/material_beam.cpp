#include "material_beam.h"
#include "element_container.h"
#include "cross_section.h"

using namespace std;

//////////////////////////////////////////////////////////
// BEAM MATERIAL STATUS

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
Vector BeamMaterialStatus :: giveStress(const Vector &strain, double timeStep) {
    //computes internal forces in local reference system (N, Vy, Vz, T, My, Mz)
    return BeamMaterialStatus :: giveStressWithFrozenIntVars(strain, timeStep);
}

//////////////////////////////////////////////////////////
Vector BeamMaterialStatus :: giveStressWithFrozenIntVars(const Vector &strain, double timeStep) {
    //computes internal forces in local reference system (N, Vy, Vz, T, My, Mz)
    ( void ) timeStep;
    temp_stress = giveStiffnessTensor("elastic") * strain;
    return temp_stress;
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
