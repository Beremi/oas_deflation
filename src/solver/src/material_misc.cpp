#include "material_misc.h"
#include "element.h"


BrittleMaterialStatus :: BrittleMaterialStatus(BrittleMaterial *m, Element *e) : DisMechMaterialStatus(m, e) {
    name = "BRITTLE mat. status";
    RAND_H = 1.0;
}

//////////////////////////////////////////////////////////
double BrittleMaterialStatus :: giveValue(string code) const {
    if ( code.rfind( "damage", 0 ) == 0 || code.rfind( "damageN", 0 ) == 0 || code.rfind( "damageT", 0 ) == 0) {
        return temp_damage;

    } else {
        return DisMechMaterialStatus :: giveValue(code);
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
    temp_normal_strain = strain[0];
    double epsT = 0;
    for ( unsigned ei = 1; ei < strain.size(); ei++){
      epsT += pow(strain[ ei ], 2);
    }
    epsT = sqrt(epsT);

    if ( !damage ){
      if ( temp_normal_strain > m->giveFt() / this->giveElasticNormalShearStiffness()[ 0 ] ) {
        temp_damage = true;
      }
      if ( epsT > m->giveFs() / this->giveElasticNormalShearStiffness()[ 1 ] ) {
        temp_damage = true;
      }
    }
}


void BrittleMaterialStatus :: update() {
    damage = temp_damage;
}

//////////////////////////////////////////////////////////
Vector BrittleMaterialStatus :: giveNormalShearStiffness(string type) const {
    Vector stiff = giveElasticNormalShearStiffness();
    if ( type.compare("elastic") == 0 ) {
        return stiff;
    } else if ( type.compare("secant") == 0 || type.compare("unloading") == 0 || type.compare("tangent") == 0 ) {
        stiff[ 1 ] = stiff[ 1 ] * ( temp_damage ? 1e-10 : 1 );
        if ( temp_normal_strain > 0 ) stiff[ 0 ] = stiff[ 0 ] * ( temp_damage ? 1e-10 : 1 );
        return stiff;
    } else {
        cerr << "Error: MarsMaterialStatus does not provide '" << type << "' stiffness";
        exit(1);
    };
}

//////////////////////////////////////////////////////////
Vector BrittleMaterialStatus :: giveStress(const Vector &strain) {
    computeDamage(strain);
    if ( damage ) {
      Vector stress, stiff;
      stress.resize(strain.size());
      stiff = this->giveNormalShearStiffness("secant");
      for ( unsigned i = 0; i < strain.size(); i++ ){
        if ( i < 2 ){
          stress[ i ] = strain[ i ] * stiff[ i ];
        } else {
          stress[ i ] = strain[ i ] * stiff[ 1 ];
        }
      }
      return stress;
    } else {
      return DisMechMaterialStatus :: giveStress(strain);
    }
}



//////////////////////////////////////////////////////////
// CUSATIS MATERIAL

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
MaterialStatus *BrittleMaterial :: giveNewMaterialStatus(Element *e) {
    BrittleMaterialStatus *newStatus = new BrittleMaterialStatus(this, e); //needs to be deleted manually
    return newStatus;
};


//////////////////////////////////////////////////////////
void BrittleMaterial :: init() {
    // if variables not specified on the input, use default multipliers
    fs = ( fs == 0 ) ? 3 * ft : fs;
};
