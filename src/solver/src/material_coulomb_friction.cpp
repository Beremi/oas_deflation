#include "material_coulomb_friction.h"
#include "element_discrete.h"
#include "model.h"

using namespace std;

//////////////////////////////////////////////////////////
// COULOMB FRICTION MATERIAL STATUS

CoulombFrictionMaterialStatus :: CoulombFrictionMaterialStatus(CoulombFrictionMaterial *m, Element *e, unsigned ipnum) : DisMechMaterialStatus(m, e, ipnum) {
    name = "Coulomb friction mat. status";
    normalStress = 0.;
}

//////////////////////////////////////////////////////////
void CoulombFrictionMaterialStatus :: giveValues(string code, Vector &result) const {
    DisMechMaterialStatus :: giveValues(code, result);
}

//////////////////////////////////////////////////////////
void CoulombFrictionMaterialStatus :: init() {
}

//////////////////////////////////////////////////////////
void CoulombFrictionMaterialStatus :: update() {
    DisMechMaterialStatus :: update();
}

//////////////////////////////////////////////////////////
void CoulombFrictionMaterialStatus :: resetTemporaryVariables() {
    DisMechMaterialStatus :: resetTemporaryVariables();
}

//////////////////////////////////////////////////////////
Matrix CoulombFrictionMaterialStatus :: giveStiffnessTensor(string type, unsigned dim) const {
    ( void ) type;
    RigidBodyBoundary *rb = static_cast< RigidBodyBoundary * >( element );
    CoulombFrictionMaterial *m = static_cast< CoulombFrictionMaterial * >( mat );
    Matrix D = Matrix :: Zero(dim, dim);
    //double strain_norm = temp_strain.norm();
    //if(strain_norm<1e-10)  strain_norm = 1e-10;
    //double stress_norm = temp_stress.norm();    
    double stiff = 0.;
    if (normalStress<=-1) stiff = m->giveInitialStiffness();
    //if (normalStress<=-1) stiff = max(stiffX,0.);
    for ( unsigned i = 1; i < dim; i++ ) {        
        D(i, i) =  stiff * rb->giveLength();
    }
    return D;
}

//////////////////////////////////////////////////////////
Vector CoulombFrictionMaterialStatus :: giveStress(const Vector &strain, double timeStep) {
    return CoulombFrictionMaterialStatus :: giveStressWithFrozenIntVars(strain, timeStep);
}

//////////////////////////////////////////////////////////
Vector CoulombFrictionMaterialStatus :: giveStressWithFrozenIntVars(const Vector &strain, double timeStep) {
    (void) timeStep;
    CoulombFrictionMaterial *m = static_cast< CoulombFrictionMaterial * >( mat );
    RigidBodyBoundary *rb = static_cast< RigidBodyBoundary * >( element );
    temp_strain = strain * rb->giveLength();
    temp_strain[0] = 0.;
    Vector strain_inc = temp_strain - updt_strain;
    temp_stress = updt_stress + strain_inc * m->giveInitialStiffness();
    double eff_stress = temp_stress.norm();
    if (normalStress>-1) temp_stress *= 0.;
    else {
        double max_stress = -normalStress*m->giveFrictionAngle();
        if (eff_stress>max_stress) temp_stress *= max_stress/eff_stress;
    }
    return temp_stress;
}


//////////////////////////////////////////////////////////
void CoulombFrictionMaterialStatus :: setParameterValue(string code, double value) {
    if ( code.compare("normal_stress") == 0 ) {
        normalStress = value;
    } else {
        DisMechMaterialStatus :: setParameterValue(code, value);
    }
}

//////////////////////////////////////////////////////////
// COULOMB FRICTION MATERIAL

//////////////////////////////////////////////////////////
void CoulombFrictionMaterial :: readFromLine(istringstream &iss) {
    iss.clear(); // clear string stream
    iss.seekg(0, iss.beg); //reset position in string stream

    // initialize all values to zero (NOTE probably no ned in linux, but in windows necessary)
    init_stiffness = density = friction_angle = 0;
    bool bis, bfa;
    bis = bfa = false;
    string param;

    while ( !iss.eof() ) {
        iss >> param;
        if ( param.compare("initStiffness") == 0 ) {
            bis = true;
            iss >> init_stiffness;
        } else if ( param.compare("frictionAngle") == 0 ) {
            bfa = true;
            iss >> friction_angle;
        } else if ( param.compare("density") == 0 ) {
            iss >> density;
        }
    }
    if ( !bis ) {
        cerr << name << ": material parameter 'initStiffness' was not specified" << endl;
        exit(EXIT_FAILURE);
    }
    ;
    if ( !bfa ) {
        cerr << name << ": material parameter 'frictionAngle' was not specified" << endl;
        exit(EXIT_FAILURE);
    }
    ;
};

//////////////////////////////////////////////////////////
CoulombFrictionMaterialStatus *CoulombFrictionMaterial :: giveNewMaterialStatus(Element *e, unsigned ipnum) {
    CoulombFrictionMaterialStatus *newStatus = new CoulombFrictionMaterialStatus(this, e, ipnum); //needs to be deleted manually
    return newStatus;
};


//////////////////////////////////////////////////////////
void CoulombFrictionMaterial :: init() {    
};

