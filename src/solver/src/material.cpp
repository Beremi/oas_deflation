#include "material.h"
#include "element.h"

//////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////
// TRANSPORT MATERIAL

TrsprtMaterialStatus :: TrsprtMaterialStatus(TrsprtMaterial *m, Element *e) : MaterialStatus(m, e) {
    name = "transport mat. status";
    effConductivity = m->giveDensity() * m->givePermeability() / m->giveViscosity();
}


//////////////////////////////////////////////////////////
Vector TrsprtMaterialStatus :: giveStress(const Vector &strain) {
    return -effConductivity * strain;
};

//////////////////////////////////////////////////////////
void TrsprtMaterial :: readFromLine(istringstream &iss) {
    string param;
    bool bcapacity, bpermeability, bviscosity, bdensity;
    bcapacity = bpermeability = bviscosity = bdensity = false;

    while ( !iss.eof() ) {
        iss >> param;
        if ( param.compare("capacity") == 0 ) {
            bcapacity = true;
            iss >> capacity;
        } else if ( param.compare("viscosity") == 0 ) {
            bviscosity = true;
            iss >> viscosity;
        } else if ( param.compare("permeability") == 0 ) {
            bpermeability = true;
            iss >> permeability;
        } else if ( param.compare("density") == 0 ) {
            bdensity = true;
            iss >> density;
        }
    }
    if ( !bcapacity ) {
        cerr << name << ": material parameter 'capacity' was not specified" << endl;
        exit(EXIT_FAILURE);
    }
    if ( !bviscosity ) {
        cerr << name << ": material parameter 'viscosity' was not specified" << endl;
        exit(EXIT_FAILURE);
    }
    if ( !bpermeability ) {
        cerr << name << ": material parameter 'permeability' was not specified" << endl;
        exit(EXIT_FAILURE);
    }
    if ( !bdensity ) {
        cerr << name << ": material parameter 'density' was not specified" << endl;
        exit(EXIT_FAILURE);
    }
};

//////////////////////////////////////////////////////////
MaterialStatus *TrsprtMaterial :: giveNewMaterialStatus(Element *e) {
    TrsprtMaterialStatus *newStatus = new TrsprtMaterialStatus(this, e); //needs to be deleted manually
    return newStatus;
};

//////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////
// COUPLED TRANSPORT MATERIAL
//////////////////////////////////////////////////////////

TrsprtCoupledMaterialStatus :: TrsprtCoupledMaterialStatus(TrsprtMaterial *m, Element *e) : TrsprtMaterialStatus(m, e) {
    name = "coupled transport mat. status";
    temp_effConductivity = effConductivity;
}

//////////////////////////////////////////////////////////
double TrsprtCoupledMaterialStatus :: giveEffectiveConductivity(string type) const {
    if ( type.compare("elastic") == 0  ) {
        TrsprtMaterial *m = static_cast< TrsprtMaterial * >( mat );
        return m->giveDensity() * m->givePermeability() / m->giveViscosity();
    } else if ( type.compare("secant") == 0 || type.compare("unloading") == 0 || type.compare("tangent") == 0 ) {
        return temp_effConductivity;
    } else {
        cerr << "Error: TrsprtCoupledMaterialStatus does not provide '" << type << "' stiffness";
        exit(1);
    };
}

//////////////////////////////////////////////////////////
void TrsprtCoupledMaterialStatus :: update() {
    effConductivity = temp_effConductivity;
};

//////////////////////////////////////////////////////////
Vector TrsprtCoupledMaterialStatus :: giveStress(const Vector &strain) {
    //first strain term is pressure gradient, second strain face are, then it is alway crack opening and crack length coulples
    TrsprtCoupledMaterial *m = static_cast< TrsprtCoupledMaterial * >( mat );
    Vector flux(1);
    double crackParam = 0;
    for ( size_t i = 2; i < strain.size(); i += 2 ) {
        crackParam += pow(strain [ i ], 3) * strain [ i + 1 ];
    }
    temp_effConductivity = m->giveDensity() * m->givePermeability() / m->giveViscosity() + m->giveTurtuosity() / ( 12. * m->giveViscosity() * strain [ 1 ] ) * crackParam;
    flux [ 0 ] = -temp_effConductivity * strain [ 0 ];
    return flux;
};

//////////////////////////////////////////////////////////
void TrsprtCoupledMaterial :: readFromLine(istringstream &iss) {
    TrsprtMaterial :: readFromLine(iss);

    iss.clear(); // clear string stream
    iss.seekg(0, iss.beg); //reset position in string stream

    string param;
    bool bturtuosity;
    bturtuosity = false;

    while ( !iss.eof() ) {
        iss >> param;
        if ( param.compare("crack_turtuosity") == 0 ) {
            bturtuosity = true;
            iss >> crack_turtuosity;
        }
    }
    if ( !bturtuosity ) {
        cerr << name << ": material parameter 'crack_turtuosity' was not specified" << endl;
        exit(EXIT_FAILURE);
    }
};

//////////////////////////////////////////////////////////
MaterialStatus *TrsprtCoupledMaterial :: giveNewMaterialStatus(Element *e) {
    TrsprtCoupledMaterialStatus *newStatus = new TrsprtCoupledMaterialStatus(this, e); //needs to be deleted manually
    return newStatus;
};

//////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////
// DISCRETE MECHANICAL ELASTIC MATERIAL

DisMechMaterialStatus :: DisMechMaterialStatus(DisMechMaterial *m, Element *e) : MaterialStatus(m, e) {
    name = "discrete mechanical mat. status";
    mat = m;
}

//////////////////////////////////////////////////////////
Vector DisMechMaterialStatus :: giveElasticNormalShearStiffness() const {
    DisMechMaterial *m = static_cast< DisMechMaterial * >( mat );
    Vector c(2);
    c [ 0 ] = m->giveE0();
    c [ 1 ] = m->giveAlpha() * m->giveE0();
    return c;
}

//////////////////////////////////////////////////////////
double DisMechMaterialStatus :: giveDensity() const {
    DisMechMaterial *tmat = static_cast< DisMechMaterial * >( mat );
    return tmat->giveDensity();
}

//////////////////////////////////////////////////////////
Vector DisMechMaterialStatus :: giveStress(const Vector &strain) {
    Vector stiff = giveElasticNormalShearStiffness();
    Vector stress(strain.size() );
    stress [ 0 ] = stiff [ 0 ] * strain [ 0 ];
    for ( unsigned i = 1; i < strain.size(); i++ ) {
        stress [ i ] = stiff [ 1 ] * strain [ i ];
    }
    return stress;
};

//////////////////////////////////////////////////////////
void DisMechMaterial :: readFromLine(istringstream &iss) {
    string param;

    bool bE0, balpha, bdensity;
    bE0 = balpha = bdensity = false;

    while ( !iss.eof() ) {
        iss >> param;
        if ( param.compare("E0") == 0 ) {
            bE0 = true;
            iss >> E0;
        } else if ( param.compare("alpha") == 0 ) {
            balpha = true;
            iss >> alpha;
        } else if ( param.compare("density") == 0 ) {
            bdensity = true;
            iss >> density;
        }
    }
    if ( !bE0 ) {
        cerr << name << ": material parameter 'E0' was not specified" << endl;
        exit(EXIT_FAILURE);
    }
    ;
    if ( !balpha ) {
        cerr << name << ": material parameter 'alpha' was not specified" << endl;
        exit(EXIT_FAILURE);
    }
    ;
    if ( !bdensity ) {
        cerr << name << ": material parameter 'density' was not specified" << endl;
        exit(EXIT_FAILURE);
    }
    ;
};

//////////////////////////////////////////////////////////
MaterialStatus *DisMechMaterial :: giveNewMaterialStatus(Element *e) {
    DisMechMaterialStatus *newStatus = new DisMechMaterialStatus(this, e); //needs to be deleted manually
    return newStatus;
};
