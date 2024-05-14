#include "material.h"
#include "element.h"
#include "material_tensorial.h"


using namespace std;


//////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////
// TRANSPORT MATERIAL

TensTrsprtMaterialStatus :: TensTrsprtMaterialStatus(TensTrsprtMaterial *m, Element *e, unsigned ipnum) : MaterialStatus(m, e, ipnum) {
    name = "transport mat. status";
    temp_pressure = 0;
    effConductivity = updateEffectiveConductivity();
}


//////////////////////////////////////////////////////////
Vector TensTrsprtMaterialStatus :: giveStress(const Vector &strain, double timeStep) {
    effConductivity = updateEffectiveConductivity(); //nonlinear effecto of pressure
    return giveStressWithFrozenIntVars(strain, timeStep);
};

//////////////////////////////////////////////////////////
Vector TensTrsprtMaterialStatus :: giveStressWithFrozenIntVars(const Vector &strain, double timeStep) {
    ( void ) timeStep;
    temp_strain = strain;
    temp_stress = -effConductivity *addEigenStrain(temp_strain); //DO NOT update here effConductivity, it is used for RVE material
    return temp_stress;
};

//////////////////////////////////////////////////////////
Matrix TensTrsprtMaterialStatus :: giveStiffnessTensor(string type) const {
    ( void ) type;
    unsigned ss = mat->giveStrainSize();
    Matrix T = Matrix :: Zero(ss, ss);
    for ( unsigned i = 0; i < ss; i++ ) {
        T(i, i) = -giveEffectiveConductivity(type);
    }
    return T;
};

//////////////////////////////////////////////////////////
Matrix TensTrsprtMaterialStatus :: giveDampingTensor() const {
    TensTrsprtMaterial *m = static_cast< TensTrsprtMaterial * >( mat );
    Matrix M = Matrix :: Zero(1, 1);
    M(0, 0) = -m->giveCapacity() * m->giveDensity();
    return M;
}

//////////////////////////////////////////////////////////
double TensTrsprtMaterialStatus :: giveEffectiveConductivity(string type) const {
    if ( type.compare("elastic") == 0 ) {
        TensTrsprtMaterial *tmat = static_cast< TensTrsprtMaterial * >( mat );
        return calculatePressureDependentPermeability(0.) * tmat->giveDensity() / tmat->giveViscosity();
    } else {
        return effConductivity;
    }
}

//////////////////////////////////////////////////////////
double TensTrsprtMaterialStatus :: calculatePressureDependentPermeability(double pressure) const {
    TensTrsprtMaterial *tmat = static_cast< TensTrsprtMaterial * >( mat );
    if ( tmat->giveParamA() < 0 ) {
        return tmat->givePermeability();
    } else {
        double m = tmat->giveParamM();
        double saturation = pow(1. + pow(pressure / tmat->giveParamA(), 1. / ( 1. - m ) ), -m);
        return tmat->givePermeability() * pow(saturation, 0.5) * pow(1. - pow(1. - pow(saturation, 1. / m), m), 2.);
    }
}

//////////////////////////////////////////////////////////
bool TensTrsprtMaterialStatus :: giveValues(string code, Vector &result) const {
    if ( code.compare("permeability") == 0 ) {
        TensTrsprtMaterial *m = dynamic_cast< TensTrsprtMaterial * >( mat );
        result.resize(1);
        result [ 0 ] =  m->givePermeability();
        return true;
    } else if ( code.compare("flux") == 0 ) {
        result.resize(1);
        result [ 0 ] =  abs(temp_stress [ 0 ]);
        return true;
    } else if ( code.compare("pressure_gradient") == 0 ) {
        result.resize(1);
        result [ 0 ] = abs(temp_strain [ 0 ]);
        return true;
    } else {
        return MaterialStatus :: giveValues(code, result);
    }
}

//////////////////////////////////////////////////////////
double TensTrsprtMaterialStatus :: updateEffectiveConductivity() const {
    TensTrsprtMaterial *tmat = static_cast< TensTrsprtMaterial * >( mat );
    return calculatePressureDependentPermeability(temp_pressure) * tmat->giveDensity() / tmat->giveViscosity();
}

//////////////////////////////////////////////////////////
bool TensTrsprtMaterialStatus :: isElastic(const bool &now) const {
    ( void ) now;
    return true; //this is not true, discuss with Pepa
}

//////////////////////////////////////////////////////////
void TensTrsprtMaterialStatus :: setParameterValue(string code, double value) {
    if ( code.compare("pressure") == 0 ) {
        temp_pressure = value;
    } else {
        MaterialStatus :: setParameterValue(code, value);
    }
}

//////////////////////////////////////////////////////////
void TensTrsprtMaterial :: readFromLine(istringstream &iss) {
    string param;
    bool bcapacity, bpermeability, bviscosity, bdensity;
    bcapacity = bpermeability = bviscosity = bdensity = false;

    while (  iss >> param ) {
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
        } else if ( param.compare("a") == 0 ) {
            iss >> a; //negative "a" means linear material behavior
        } else if ( param.compare("m") == 0 ) {
            iss >> m;
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
MaterialStatus *TensTrsprtMaterial :: giveNewMaterialStatus(Element *e, unsigned ipnum) {
    TensTrsprtMaterialStatus *newStatus = new TensTrsprtMaterialStatus(this, e, ipnum); //needs to be deleted manually
    return newStatus;
};

//////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////
// HEAT CONDUCTION MATERIAL

TensHeatConductionMaterialStatus :: TensHeatConductionMaterialStatus(TensHeatConductionMaterial *m, Element *e, unsigned ipnum) : MaterialStatus(m, e, ipnum) {
    name = "heat conduction mat. status";
    temp_temperature = 0;
}


//////////////////////////////////////////////////////////
Vector TensHeatConductionMaterialStatus :: giveStress(const Vector &strain, double timeStep) {
    return giveStressWithFrozenIntVars(strain, timeStep);
};

//////////////////////////////////////////////////////////
Vector TensHeatConductionMaterialStatus :: giveStressWithFrozenIntVars(const Vector &strain, double timeStep) {
    ( void ) timeStep;
    TensHeatConductionMaterial *m = static_cast< TensHeatConductionMaterial * >( mat );
    temp_strain = strain;
    temp_stress = m->giveConductivity() * addEigenStrain(temp_strain); //DO NOT update here effConductivity, it is used for RVE material
    return temp_stress;
};

//////////////////////////////////////////////////////////
Matrix TensHeatConductionMaterialStatus :: giveStiffnessTensor(string type) const {
    ( void ) type;
    TensHeatConductionMaterial *m = static_cast< TensHeatConductionMaterial * >( mat );
    unsigned ss = mat->giveStrainSize();
    Matrix T = Matrix :: Zero(ss, ss);
    for ( unsigned i = 0; i < ss; i++ ) {
        T(i, i) = m->giveConductivity();
    }
    return T;
};

//////////////////////////////////////////////////////////
Matrix TensHeatConductionMaterialStatus :: giveDampingTensor() const {
    TensHeatConductionMaterial *m = static_cast< TensHeatConductionMaterial * >( mat );
    Matrix M = Matrix :: Zero(1, 1);
    M(0, 0) = m->giveCapacity() * m->giveDensity();
    return M;
}

//////////////////////////////////////////////////////////
bool TensHeatConductionMaterialStatus :: giveValues(string code, Vector &result) const {
    if ( code.compare("heat_flux") == 0 ) {
        result.resize(1);
        result [ 0 ] =  abs(temp_stress [ 0 ]);
        return true;
    } else if ( code.compare("temperature_gradient") == 0 ) {
        result.resize(1);
        result [ 0 ] = abs(temp_strain [ 0 ]);
        return true;
    } else {
        return MaterialStatus :: giveValues(code, result);
    }
}

//////////////////////////////////////////////////////////
void TensHeatConductionMaterialStatus :: setParameterValue(string code, double value) {
    if ( code.compare("temperature") == 0 ) {
        temp_temperature = value;
    } else {
        MaterialStatus :: setParameterValue(code, value);
    }
}

//////////////////////////////////////////////////////////
void TensHeatConductionMaterial :: readFromLine(istringstream &iss) {
    string param;
    bool bcapacity, bconductivity, bdensity;
    bcapacity = bconductivity = bdensity = false;

    while (  iss >> param ) {
        if ( param.compare("capacity") == 0 ) {
            bcapacity = true;
            iss >> capacity;
        } else if ( param.compare("conductivity") == 0 ) {
            bconductivity = true;
            iss >> conductivity;
        } else if ( param.compare("density") == 0 ) {
            bdensity = true;
            iss >> density;
        }
    }
    if ( !bcapacity ) {
        cerr << name << ": material parameter 'capacity' was not specified" << endl;
        exit(EXIT_FAILURE);
    }
    if ( !bconductivity ) {
        cerr << name << ": material parameter 'conductivity' was not specified" << endl;
        exit(EXIT_FAILURE);
    }
    if ( !bdensity ) {
        cerr << name << ": material parameter 'density' was not specified" << endl;
        exit(EXIT_FAILURE);
    }
};

//////////////////////////////////////////////////////////
MaterialStatus *TensHeatConductionMaterial :: giveNewMaterialStatus(Element *e, unsigned ipnum) {
    TensHeatConductionMaterialStatus *newStatus = new TensHeatConductionMaterialStatus(this, e, ipnum); //needs to be deleted manually
    return newStatus;
};

//////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////
// ELASTIC TENSORIAL MECHANICAL MATERIAL

Vector TensMechMaterialStatus :: giveStress(const Vector &strain, double timeStep) {
    return TensMechMaterialStatus :: giveStressWithFrozenIntVars(strain, timeStep);
};

//////////////////////////////////////////////////////////
Matrix TensMechMaterialStatus :: giveDampingTensor() const {
    //TensMechMaterial *m = static_cast< TensMechMaterial * >( mat );
    Matrix M = Matrix :: Zero(1, 1);
    //M(0, 0) = m->giveDampingConstant();
    M(0, 0) = 0;
    return M;
}

//////////////////////////////////////////////////////////
double TensMechMaterialStatus :: giveMassConstant() const {
    TensMechMaterial *material = static_cast< TensMechMaterial * >( mat );
    return material->giveDensity();
}

//////////////////////////////////////////////////////////
Vector TensMechMaterialStatus :: giveStressWithFrozenIntVars(const Vector &strain, double timeStep) {
    ( void ) timeStep;
    temp_strain = addEigenStrain(strain);
    temp_stress = giveStiffnessTensor("elastic") * temp_strain;
    return temp_stress;
};

//////////////////////////////////////////////////////////
Matrix TensMechMaterialStatus :: giveStiffnessTensor(string type) const {
    ( void ) type;
    unsigned size = 1;
    unsigned dimension = mat->giveDimension();
    if ( dimension == 1 ) {
        size = 1;
    } else if ( dimension == 2 ) {
        size = 3;
    } else if ( dimension == 3 ) {
        return giveElasticStiffnessTensor3D();
    } else {
        cerr << name << ": unsupported dimension " << dimension << endl;
        exit(1);
    }
    Matrix D = Matrix :: Zero(size, size);
    TensMechMaterial *m = static_cast< TensMechMaterial * >( mat );
    if ( dimension == 1 ) {
        D(0, 0) = m->giveElasticModulus();
    } else if ( dimension == 2 ) {
        if ( m->isPlaneStress() ) {
            double factor = m->giveElasticModulus() / ( 1. - pow(m->givePoissonsRatio(), 2) );
            D(0, 0) = D(1, 1) = factor;
            D(0, 1) = D(1, 0) = m->givePoissonsRatio() * factor;
            D(2, 2) = ( 1. - m->givePoissonsRatio() ) / 2. * factor;
        } else {  //plane strain
            double factor = m->giveElasticModulus() / (( 1. - 2. * m->givePoissonsRatio() ) * ( 1. + m->givePoissonsRatio() ));
            D(0, 0) = D(1, 1) = factor * ( 1. - m->givePoissonsRatio() );
            D(0, 1) = D(1, 0) = m->givePoissonsRatio() * factor;
            D(2, 2) = ( 1. - 2. * m->givePoissonsRatio() ) / 2. * factor;
        }
    } else {
        cerr << name << " error: dimension " << dimension << " not implemented" << endl;
        exit(1);
    }
    return D;
};

//////////////////////////////////////////////////////////

Matrix  TensMechMaterialStatus :: giveElasticStiffnessTensor3D() const{
    Matrix D = Matrix :: Zero(6, 6);
    TensMechMaterial *m = static_cast< TensMechMaterial * >( mat );
    double factor = m->giveElasticModulus() / ( 1. - 2. * m->givePoissonsRatio() ) / ( 1. + m->givePoissonsRatio() );
    D(0, 0) = D(1, 1) = D(2, 2) = factor * ( 1. - m->givePoissonsRatio() );
    D(0, 1) = D(1, 0) = D(0, 2) = D(2, 0) = D(2, 1) = D(1, 2) = m->givePoissonsRatio() * factor;
    D(3, 3) = D(4, 4) = D(5, 5) = ( 1. - 2. * m->givePoissonsRatio() ) / 2. * factor;
    return D;
}

//////////////////////////////////////////////////////////
void TensMechMaterialStatus :: update() {
    MaterialStatus :: update();
}

//////////////////////////////////////////////////////////
Matrix TensMechMaterialStatus :: giveMassTensor() const {
    TensMechMaterial *m = static_cast< TensMechMaterial * >( mat );
    unsigned dimension = mat->giveDimension();
    return Eigen :: MatrixXd :: Identity(dimension, dimension) * m->giveDensity();
}

//////////////////////////////////////////////////////////
TensMechMaterialStatus :: TensMechMaterialStatus(TensMechMaterial *m, Element *e, unsigned ipnum) : MaterialStatus(m, e, ipnum) {
    name = "tensorial mechanical mat. status";
}

//////////////////////////////////////////////////////////
bool TensMechMaterialStatus :: giveValues(string code, Vector &result) const {
    if ( code.compare("stress") == 0 || code.compare("stresses") == 0 ) {
        unsigned size = temp_stress.size();
        result.resize(size);
        for ( unsigned p = 0; p < size; p++ ) {
            result [ p ] = temp_stress [ p ];
        }
        return true;
    } else if ( code.compare("strain") == 0 || code.compare("strains") == 0 ) {
        unsigned size = temp_strain.size();
        result.resize(size);
        for ( unsigned p = 0; p < size; p++ ) {
            result [ p ] = temp_strain [ p ];
        }
        return true;
    } else {
        return MaterialStatus :: giveValues(code, result);
    }
}

//////////////////////////////////////////////////////////
void TensMechMaterial :: readFromLine(istringstream &iss) {
    string param;
    bool bE, bnu, bdensity;
    bE = bnu = bdensity = false;

    while (  iss >> param ) {
        if ( param.compare("E") == 0 ) {
            bE = true;
            iss >> E;
        } else if ( param.compare("nu") == 0 ) {
            bnu = true;
            iss >> nu;
        } else if ( param.compare("density") == 0 ) {
            bdensity = true;
            iss >> density;
        } else if ( param.compare("planeStrain") == 0 ) {
            planeStress = false;
        }
    }
    if ( !bE ) {
        cerr << name << ": material parameter 'E' was not specified" << endl;
        exit(EXIT_FAILURE);
    }
    if ( !bnu ) {
        cerr << name << ": material parameter 'nu' was not specified" << endl;
        exit(EXIT_FAILURE);
    }
    if ( !bdensity ) {
        cerr << name << ": material parameter 'density' was not specified" << endl;
        exit(EXIT_FAILURE);
    }
};

//////////////////////////////////////////////////////////
MaterialStatus *TensMechMaterial :: giveNewMaterialStatus(Element *e, unsigned ipnum) {
    TensMechMaterialStatus *newStatus = new TensMechMaterialStatus(this, e, ipnum); //needs to be deleted manually
    return newStatus;
};

//////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////
// ELASTIC COSSERAT MECHANICAL MATERIAL

TensCosseratMechMaterial::TensCosseratMechMaterial(unsigned dimension) : TensMechMaterial(dimension) {
    name = "elastic Cosserat mechanical material";
    if (dim==2){
        planeStress = true;
        strainsize = 6;
    }else if (dim==3){
        strainsize = 18;        
    }
}

//////////////////////////////////////////////////////////
Vector TensCosseratMechMaterialStatus ::  giveStress(const Vector &strain, double timeStep) {
    return TensCosseratMechMaterialStatus ::  giveStressWithFrozenIntVars(strain, timeStep);
};

//////////////////////////////////////////////////////////
Vector TensCosseratMechMaterialStatus ::  giveStressWithFrozenIntVars(const Vector &strain, double timeStep) {
    ( void ) timeStep;
    temp_strain = addEigenStrain(strain);
    temp_stress = giveStiffnessTensor("elastic") *  temp_strain;
    return temp_stress;
};

//////////////////////////////////////////////////////////
Matrix TensCosseratMechMaterialStatus :: giveStiffnessTensor(string type) const {
    ( void ) type;
    unsigned size;
    unsigned dimension = mat->giveDimension();
    if ( dimension == 1 ) {
        size = 1;
    } else if ( dimension == 2 ) {
        size = 6;
    } else if ( dimension == 3 ) {
        size = 18;
    } else {
        cerr << name << ": unsupported dimension " << dimension << endl;
        exit(1);
    }
    Matrix D = Matrix :: Zero(size, size);
    TensCosseratMechMaterial *m = static_cast< TensCosseratMechMaterial * >( mat );
    double lammeL = m->giveElasticModulus() * m->givePoissonsRatio() / ( ( 1. + m->givePoissonsRatio() ) * ( 1. - 2 * m->givePoissonsRatio() ) );
    double lammeM = m->giveElasticModulus() / ( 2. * ( 1. + m->givePoissonsRatio() ) );
    if ( dimension == 1 ) {
        D(0, 0) = m->giveElasticModulus();
        D(1, 1) = 4. * m->giveCosseratShearParam() * pow(m->giveCharacteristicLength(), 2);
    } else if ( dimension == 2 ) {
        if ( m->isPlaneStress() ) {
            cerr << name << " error: dimension " << dimension << ", Cosserat plane stress not implemented" << endl;
            exit(1);
        } else {  //plane strain
            D(0, 0) = D(1, 1) = lammeL + 2. * lammeM;
            D(0, 1) = D(1, 0) = lammeL;
            D(2, 2) = D(3, 3) = lammeM + m->giveCosseratShearParam();
            D(2, 3) = D(3, 2) = lammeM - m->giveCosseratShearParam();
            D(4, 4) = D(5, 5) = lammeM * 4. * pow(m->giveCharacteristicLength(), 2);
        }
    } else {
            D(0, 0) = D(1, 1)  = D(2, 2) = lammeL + 2. * lammeM;
            D(0, 1) = D(1, 0) = D(0, 2) = D(2, 0) = D(2, 1) = D(1, 2) = lammeL;
            D(3, 3) = D(4, 4) = D(5, 5) = D(6, 6) = D(7, 7) = D(8, 8) = lammeM + m->giveCosseratShearParam();
            D(4, 3) = D(3, 4) = D(6, 5) = D(5, 6) = D(8, 7) = D(7, 8) = lammeM - m->giveCosseratShearParam();
            D(9, 9) = D(10, 10) = D(11, 11) = D(12, 12) = D(13, 13) = D(14, 14) = D(15, 15) = D(16, 16) = D(17, 17) = lammeM * 4. * pow(m->giveCharacteristicLength(), 2);
    }
    return D;
};

//////////////////////////////////////////////////////////
TensCosseratMechMaterialStatus :: TensCosseratMechMaterialStatus(TensCosseratMechMaterial *m, Element *e, unsigned ipnum) : TensMechMaterialStatus(m, e, ipnum) {
    name = "tensorial cosserat mat. status";
}

//////////////////////////////////////////////////////////
void TensCosseratMechMaterial :: readFromLine(istringstream &iss) {
    TensMechMaterial :: readFromLine(iss);

    iss.clear(); // clear string stream
    iss.seekg(0, iss.beg); //reset position in string stream

    string param;
    bool blc, bmuc;
    blc = bmuc = false;

    while (  iss >> param ) {
        if ( param.compare("lc") == 0 ) {
            blc = true;
            iss >> lc;
        } else if ( param.compare("muc") == 0 ) {
            bmuc = true;
            iss >> muc;
        }
    }
    if ( !blc ) {
        cerr << name << ": material parameter 'lc' was not specified" << endl;
        exit(EXIT_FAILURE);
    }
    if ( !bmuc ) {
        cerr << name << ": material parameter 'muc' was not specified" << endl;
        exit(EXIT_FAILURE);
    }
};

//////////////////////////////////////////////////////////
MaterialStatus *TensCosseratMechMaterial :: giveNewMaterialStatus(Element *e, unsigned ipnum) {
    TensCosseratMechMaterialStatus *newStatus = new TensCosseratMechMaterialStatus(this, e, ipnum);
    return newStatus;
};
