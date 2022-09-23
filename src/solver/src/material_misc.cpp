#include "material_misc.h"
#include "element_discrete.h"

using namespace std;

BrittleMaterialStatus :: BrittleMaterialStatus(BrittleMaterial *m, Element *e, unsigned ipnum) : DisMechMaterialStatus(m, e, ipnum) {
    name = "BRITTLE mat. status";
    RAND_H = 1.0;
}

//////////////////////////////////////////////////////////
void BrittleMaterialStatus :: giveValues(string code, Vector &result) const {
    if ( code.rfind("damage", 0) == 0 || code.rfind("damageN", 0) == 0 || code.rfind("damageT", 0) == 0 ) {
        result.resize(1);
        result [ 0 ] = temp_damage;
    } else {
        DisMechMaterialStatus :: giveValues(code, result);
    }
}



//////////////////////////////////////////////////////////
void BrittleMaterialStatus :: init() {
    damage = temp_damage = false;
    temp_normal_strain = 0;

    RigidBodyContact *rbc = dynamic_cast< RigidBodyContact * >( element );
    if ( !rbc ) {
        cerr << "Material " << name << " can be used only for RigidBodyContact elements" << endl;
        exit(EXIT_FAILURE);
    }
    L = rbc->giveLength();
}


//////////////////////////////////////////////////////////
void BrittleMaterialStatus :: computeDamage(const Vector &strain) {
    BrittleMaterial *m = static_cast< BrittleMaterial * >( mat );
    temp_normal_strain = strain [ 0 ];
    double epsT = 0;
    for ( unsigned ei = 1; ei < strain.size(); ei++ ) {
        epsT += pow(strain [ ei ], 2);
    }
    epsT = sqrt(epsT);

    if ( !damage ) {
        if ( temp_normal_strain > m->giveFt() / m->giveE0() ) {
            temp_damage = true;
        }
        if ( epsT > m->giveFs() / ( m->giveE0() * m->giveAlpha() ) ) {
            temp_damage = true;
        }
    }
}


void BrittleMaterialStatus :: update() {
    damage = temp_damage;
}

//////////////////////////////////////////////////////////
Matrix BrittleMaterialStatus :: giveStiffnessTensor(string type, unsigned dim) const {
    Matrix stiff = DisMechMaterialStatus :: giveStiffnessTensor(type, dim);
    if ( type.compare("elastic") == 0 ) {
        return stiff;
    } else if ( type.compare("secant") == 0 || type.compare("unloading") == 0 || type.compare("tangent") == 0 ) {
        if ( temp_normal_strain > 0 ) {
            stiff(0, 0) *= ( temp_damage ? 1e-10 : 1 );
        }
        for ( size_t i = 1; i < dim; i++ ) {
            stiff(i, i) *= ( temp_damage ? 1e-10 : 1 );
        }
        return stiff;
    } else {
        cerr << "Error: BrittleMaterialStatus does not provide '" << type << "' stiffness";
        exit(1);
    };
}

//////////////////////////////////////////////////////////
Vector BrittleMaterialStatus :: giveStress(const Vector &strain, double timeStep) {
    computeDamage(strain);
    if ( damage ) {
        Matrix stiff = giveStiffnessTensor("secant", strain.size() );
        return stiff * strain;
    } else {
        return DisMechMaterialStatus :: giveStress(strain, timeStep);
    }
}


//////////////////////////////////////////////////////////
Vector BrittleMaterialStatus :: giveStressWithFrozenIntVars(const Vector &strain, double timeStep) {
    ( void ) timeStep;
    ( void ) strain;
    return Vector(0);  //TOTO: FIX
}

//////////////////////////////////////////////////////////
// BRITTLE MATERIAL

//////////////////////////////////////////////////////////
void BrittleMaterial :: readFromLine(istringstream &iss) {
    DisMechMaterial :: readFromLine(iss); //read elastic parameters

    iss.clear(); // clear string stream
    iss.seekg(0, iss.beg); //reset position in string stream

    ft = fs = 0;

    string param;
    bool bft;
    bft = false;

    while ( !iss.eof() ) {
        iss >> param;
        if ( param.compare("ft") == 0 ) {
            bft = true;
            iss >> ft;
        } else if ( param.compare("fs") == 0 ) {
            iss >> fs;
        }
    }
    if ( !bft ) {
        cerr << name << ": material parameter 'ft' was not specified" << endl;
        exit(EXIT_FAILURE);
    }
    ;
};

//////////////////////////////////////////////////////////
MaterialStatus *BrittleMaterial :: giveNewMaterialStatus(Element *e, unsigned ipnum) {
    BrittleMaterialStatus *newStatus = new BrittleMaterialStatus(this, e, ipnum); //needs to be deleted manually
    return newStatus;
};


//////////////////////////////////////////////////////////
void BrittleMaterial :: init() {
    // if variables not specified on the input, use default multipliers
    fs = ( fs == 0 ) ? 3 * ft : fs;
};

///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////


ContactMaterialStatus :: ContactMaterialStatus(ContactMaterial *m, Element *e, unsigned ipnum) : DisMechMaterialStatus(m, e, ipnum) {
    name = "CONTACT mat. status";
}


//////////////////////////////////////////////////////////
void ContactMaterialStatus :: init() {
    if ( !dynamic_cast< RigidBodyContact * >( element ) ) {
        cerr << "Material " << name << " can be used only for RigidBodyContact elements" << endl;
        exit(EXIT_FAILURE);
    }
}


//////////////////////////////////////////////////////////
void ContactMaterialStatus :: update() {
    temp_normal_strain = 0;
}

//////////////////////////////////////////////////////////
Matrix ContactMaterialStatus :: giveStiffnessTensor(string type, unsigned dim) const {
    Matrix stiff = DisMechMaterialStatus :: giveStiffnessTensor(type, dim);
    if ( type.compare("elastic") == 0 ) {
        return stiff;
        // this is unfortunately wrong, cannot be done this way, shear stress will change in abruptly
    } else if ( type.compare("secant") == 0 || type.compare("unloading") == 0 || type.compare("tangent") == 0 ) {
        if ( temp_normal_strain > 0 ) {
            stiff = stiff * 1e-10;
        }
        return stiff;
    } else {
        cerr << "Error: ContactMaterialStatus does not provide '" << type << "' stiffness";
        exit(1);
    };
}

//////////////////////////////////////////////////////////
Vector ContactMaterialStatus :: giveStress(const Vector &strain, double timeStep) {
    ( void ) timeStep;
    temp_normal_strain = strain [ 0 ];
    Vector stress = Vector :: Zero(strain.size() );
    if ( temp_normal_strain < 0 ) {
        ContactMaterial *m = static_cast< ContactMaterial * >( mat );
        stress [ 0 ] = strain [ 0 ] * m->giveE0();
        stress [ 1 ] = stress [ 0 ] * m->giveFrictionCoef();
        if ( strain.size() > 2 ) {
            stress [ 2 ] = stress [ 0 ] * m->giveFrictionCoef();
        }
    }
    return stress;
}


//////////////////////////////////////////////////////////
Vector ContactMaterialStatus :: giveStressWithFrozenIntVars(const Vector &strain, double timeStep) {
    ( void ) timeStep;
    ( void ) strain;
    return Vector(0);  //TOTO: FIX
}

//////////////////////////////////////////////////////////
// BRITTLE MATERIAL

//////////////////////////////////////////////////////////
void ContactMaterial :: readFromLine(istringstream &iss) {
    DisMechMaterial :: readFromLine(iss); //read elastic parameters

    iss.clear(); // clear string stream
    iss.seekg(0, iss.beg); //reset position in string stream

    friction_coef = 0;

    string param;

    while ( !iss.eof() ) {
        iss >> param;
        if ( param.compare("friction_coef") == 0 ) {
            iss >> friction_coef;
        }
    }
}

//////////////////////////////////////////////////////////
MaterialStatus *ContactMaterial :: giveNewMaterialStatus(Element *e, unsigned ipnum) {
    ContactMaterialStatus *newStatus = new ContactMaterialStatus(this, e, ipnum); //needs to be deleted manually
    return newStatus;
}


//////////////////////////////////////////////////////////
void ContactMaterial :: init() {
    // if variables not specified on the input, use default multipliers
}
