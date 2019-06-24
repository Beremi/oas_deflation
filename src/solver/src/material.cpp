#include "material.h"
#include "element.h"

//////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////
// TRANSPORT MATERIAL

TrsprtMaterialStatus :: TrsprtMaterialStatus(TrsprtMaterial *m, Element *e) : MaterialStatus(m, e) {
    name = "transport mat. status";
}

//////////////////////////////////////////////////////////
double TrsprtMaterialStatus :: giveConductivity() const {
    TrsprtMaterial *tmat = static_cast< TrsprtMaterial * >( mat );
    return tmat->giveConductivity();
}

//////////////////////////////////////////////////////////
double TrsprtMaterialStatus :: giveCapacity() const {
    TrsprtMaterial *tmat = static_cast< TrsprtMaterial * >( mat );
    return tmat->giveCapacity();
}

//////////////////////////////////////////////////////////
void TrsprtMaterial :: readFromLine(istringstream &iss) {
    string param;
    bool bcapacity, bconductivity;
    bcapacity = bconductivity = false;

    while ( !iss.eof() ) {
        iss >> param;
        if ( param.compare("capacity") == 0 ) {
            bcapacity = true;
            iss >> capacity;
        } else if ( param.compare("conductivity") == 0 )    {
            bconductivity = true;
            iss >> conductivity;
        }
    }
    if ( !bcapacity ) {
        cerr << name << ": material parameter 'capacity' was not specified" << endl;
        exit(0);
    }
    ;
    if ( !bconductivity ) {
        cerr << name << ": material parameter 'conductivity' was not specified" << endl;
        exit(0);
    }
    ;
};

//////////////////////////////////////////////////////////
MaterialStatus *TrsprtMaterial :: giveNewMaterialStatus(Element *e) {
    TrsprtMaterialStatus *newStatus = new TrsprtMaterialStatus(this, e); //needs to be deleted manually
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
        } else if ( param.compare("alpha") == 0 )    {
            balpha = true;
            iss >> alpha;
        } else if ( param.compare("density") == 0 )    {
            bdensity = true;
            iss >> density;
        }
    }
    if ( !bE0 ) {
        cerr << name << ": material parameter 'E0' was not specified" << endl;
        exit(0);
    }
    ;
    if ( !balpha ) {
        cerr << name << ": material parameter 'alpha' was not specified" << endl;
        exit(0);
    }
    ;
    if ( !bdensity ) {
        cerr << name << ": material parameter 'density' was not specified" << endl;
        exit(0);
    }
    ;
};

//////////////////////////////////////////////////////////
MaterialStatus *DisMechMaterial :: giveNewMaterialStatus(Element *e) {
    DisMechMaterialStatus *newStatus = new DisMechMaterialStatus(this, e); //needs to be deleted manually
    return newStatus;
};
