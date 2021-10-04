#include "material.h"
#include "element.h"
#include "element_discrete.h"

//////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////
// BASIC MATERIAL
Vector MaterialStatus :: addEigenStrain(const Vector &totalStrain) const {
    if ( eigenstrain.size() > 0 ) {
        if ( eigenstrain.size() != totalStrain.size() ) {
            cerr << "Material status error: cannot apply eigenstrain of size " << eigenstrain.size() << " to  total strain of size " << totalStrain.size() << endl;
            exit(1);
        }
        Vector activeStrain = totalStrain - eigenstrain;
        return activeStrain;
    } else {
        return totalStrain;
    }
}

//////////////////////////////////////////////////////////
void MaterialStatus :: update() {
    updt_strain = temp_strain;
    updt_stress = temp_stress;
}


//////////////////////////////////////////////////////////
void MaterialStatus :: setEigenStrain(Vector &x) {
    eigenstrain = x;
}

//////////////////////////////////////////////////////////
bool MaterialStatus :: isElastic(const bool &now) const {
    ( void ) now;
    if ( this->name != "basic mat. status" ) {
        std :: cout << "using elastic check for base class MaterialStatus, if this is not a desire, you need to implement method \'isElastic\' for " << this->name << '\n';
    }
    return true;
}


//////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////
// TRANSPORT MATERIAL

TrsprtMaterialStatus :: TrsprtMaterialStatus(TrsprtMaterial *m, Element *e, unsigned ipnum) : MaterialStatus(m, e, ipnum) {
    name = "transport mat. status";
    avgPressure = 0;
    effConductivity = updateEffectiveConductivity();
}


//////////////////////////////////////////////////////////
Vector TrsprtMaterialStatus :: giveStress(const Vector &strain, double timeStep) {
    effConductivity = updateEffectiveConductivity(); //nonlinear effecto of pressure
    return giveStressWithFrozenIntVars(strain, timeStep);
};

//////////////////////////////////////////////////////////
Vector TrsprtMaterialStatus :: giveStressWithFrozenIntVars(const Vector &strain, double timeStep) {
    ( void ) timeStep;
    temp_strain = strain;
    temp_stress = -effConductivity *addEigenStrain(temp_strain); //DO NOT update here effConductivity, it is used for RVE material
    return temp_stress;
};

//////////////////////////////////////////////////////////
Matrix TrsprtMaterialStatus :: giveStiffnessTensor(string type, unsigned dimension) const {
    ( void ) type;
    Matrix T(dimension, dimension);
    for ( unsigned i = 0; i < dimension; i++ ) {
        T [ i ] [ i ] = - giveEffectiveConductivity(type);
    }
    return T;
};

//////////////////////////////////////////////////////////
Matrix TrsprtMaterialStatus :: giveDampingTensor() const {
    TrsprtMaterial *m = static_cast< TrsprtMaterial * >( mat );
    Matrix M(1, 1);
    M [ 0 ] [ 0 ] = -m->giveCapacity() * m->giveDensity();
    return M;
}

//////////////////////////////////////////////////////////
double TrsprtMaterialStatus :: giveEffectiveConductivity(string type) const {
    if ( type.compare("elastic") == 0 ) {
        TrsprtMaterial *tmat = static_cast< TrsprtMaterial * >( mat );
        return calculatePressureDependentPermeability(0.) * tmat->giveDensity() / tmat->giveViscosity();
    } else {
        return effConductivity;
    }
}

//////////////////////////////////////////////////////////
double TrsprtMaterialStatus :: calculatePressureDependentPermeability(double pressure) const {
    TrsprtMaterial *tmat = static_cast< TrsprtMaterial * >( mat );
    if ( tmat->giveParamA() < 0 ) {
        return tmat->givePermeability();
    } else {
        double m = tmat->giveParamM();
        double saturation = pow(1. + pow(pressure / tmat->giveParamA(), 1. / ( 1. - m ) ), -m);
        return tmat->givePermeability() * pow(saturation, 0.5) * pow(1. - pow(1. - pow(saturation, 1. / m), m), 2.);
    }
}

//////////////////////////////////////////////////////////
double TrsprtMaterialStatus :: giveValue(string code) const {
    if ( code.compare("permeability") == 0 ) {
        TrsprtMaterial *m = dynamic_cast< TrsprtMaterial * >( mat );
        return m->givePermeability();
    } else if ( code.compare("flux") == 0 ) {
        return abs(temp_stress [ 0 ]);
    } else if ( code.compare("pressure_gradient") == 0 ) {
        return abs(temp_strain [ 0 ]);
    } else {
        return MaterialStatus :: giveValue(code);
    }
}

//////////////////////////////////////////////////////////
double TrsprtMaterialStatus :: updateEffectiveConductivity() const {
    TrsprtMaterial *tmat = static_cast< TrsprtMaterial * >( mat );
    return calculatePressureDependentPermeability(avgPressure) * tmat->giveDensity() / tmat->giveViscosity();
}

//////////////////////////////////////////////////////////
bool TrsprtMaterialStatus :: isElastic(const bool &now) const {
    ( void ) now;
    return true; //this is not true, discuss with Pepa
}

//////////////////////////////////////////////////////////
void TrsprtMaterialStatus :: setParameterValue(string code, double value) {
    if ( code.compare("pressure") == 0 ) {
        avgPressure = value;
    } else {
        MaterialStatus :: setParameterValue(code, value);
    }
}



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
MaterialStatus *TrsprtMaterial :: giveNewMaterialStatus(Element *e, unsigned ipnum) {
    TrsprtMaterialStatus *newStatus = new TrsprtMaterialStatus(this, e, ipnum); //needs to be deleted manually
    return newStatus;
};

//////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////
// DISCRETE TRANSPORT MATERIAL

DiscreteTrsprtMaterialStatus :: DiscreteTrsprtMaterialStatus(TrsprtMaterial *m, Element *e, unsigned ipnum) : TrsprtMaterialStatus(m, e, ipnum) {
    name = "transport mat. status";
}

//////////////////////////////////////////////////////////
Matrix DiscreteTrsprtMaterialStatus :: giveStiffnessTensor(string type, unsigned dimension) const {
    ( void ) type;
    ( void ) dimension;
    Matrix T(1, 1); //discrete material, only one direction in any dimension
    T [ 0 ] [ 0 ] = - giveEffectiveConductivity(type);
    return T;
};

//////////////////////////////////////////////////////////
MaterialStatus *DiscreteTrsprtMaterial :: giveNewMaterialStatus(Element *e, unsigned ipnum) {
    DiscreteTrsprtMaterialStatus *newStatus = new DiscreteTrsprtMaterialStatus(this, e, ipnum); //needs to be deleted manually
    return newStatus;
};

//////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////
// ELASTIC TENSORIAL MECHANICAL MATERIAL

Vector ElasticMechMaterialStatus :: giveStress(const Vector &strain, double timeStep) {
    return ElasticMechMaterialStatus :: giveStressWithFrozenIntVars(strain, timeStep);
};


//////////////////////////////////////////////////////////
double ElasticMechMaterialStatus :: giveMassConstant() const {
    ElasticMechMaterial *material = static_cast< ElasticMechMaterial * >( mat );
    return material->giveDensity();
}

//////////////////////////////////////////////////////////
Vector ElasticMechMaterialStatus :: giveStressWithFrozenIntVars(const Vector &strain, double timeStep) {
    ( void ) timeStep;
    temp_strain = strain;
    unsigned dim = 0;
    if ( strain.size() == 1 ) {
        dim = 1;
    } else if ( strain.size() == 3 ) {
        dim = 2;
    } else if ( strain.size() == 6 ) {
        dim = 3;
    }
    temp_stress = giveStiffnessTensor("elastic", dim) * addEigenStrain(strain);
    return temp_stress;
};

//////////////////////////////////////////////////////////
Matrix ElasticMechMaterialStatus :: giveStiffnessTensor(string type, unsigned dimension) const {
    ( void ) type;
    unsigned size = 1;
    if ( dimension == 1 ) {
        size = 1;
    } else if ( dimension == 2 ) {
        size = 3;
    } else if ( dimension == 3 ) {
        size = 6;
    } else {
        cerr << name << ": unsupported dimension " << dimension << endl;
        exit(1);
    }
    Matrix D(size, size);
    ElasticMechMaterial *m = static_cast< ElasticMechMaterial * >( mat );
    if ( dimension == 1 ) {
        D [ 0 ] [ 0 ] = m->giveElasticModulus();
    } else if ( dimension == 2 ) {
        if ( m->isPlaneStress() ) {
            double factor = m->giveElasticModulus() / ( 1. - pow(m->givePoissonsRatio(), 2) );
            D [ 0 ] [ 0 ] = D [ 1 ] [ 1 ] = factor;
            D [ 0 ] [ 1 ] = D [ 1 ] [ 0 ] = m->givePoissonsRatio() * factor;
            D [ 2 ] [ 2 ] = ( 1. - m->givePoissonsRatio() ) / 2. * factor;
        } else {  //plane strain
            double factor = m->giveElasticModulus() / ( 1. - 2. * m->givePoissonsRatio() ) / ( 1. + m->givePoissonsRatio() );
            D [ 0 ] [ 0 ] = D [ 1 ] [ 1 ] = factor * ( 1. - m->givePoissonsRatio() );
            D [ 0 ] [ 1 ] = D [ 1 ] [ 0 ] = m->givePoissonsRatio() * factor;
            D [ 2 ] [ 2 ] = ( 1. - 2. * m->givePoissonsRatio() ) / 2. * factor;
        }
    } else if ( dimension == 3 ) {
        double factor = m->giveElasticModulus() / ( 1. - 2. * m->givePoissonsRatio() ) / ( 1. + m->givePoissonsRatio() );
        D [ 0 ] [ 0 ] = D [ 1 ] [ 1 ] = D [ 2 ] [ 2 ] = factor * ( 1. - m->givePoissonsRatio() );
        D [ 0 ] [ 1 ] = D [ 1 ] [ 0 ] = D [ 0 ] [ 2 ] = D [ 2 ] [ 0 ] = D [ 2 ] [ 1 ] = D [ 1 ] [ 2 ] = m->givePoissonsRatio() * factor;
        D [ 3 ] [ 3 ] = D [ 4 ] [ 4 ] = D [ 5 ] [ 5 ] = ( 1. - 2. * m->givePoissonsRatio() ) / 2. * factor;
    } else {
        cerr << name << " error: dimension " << dimension << " not implemented" << endl;
        exit(1);
    }
    return D;
};

//////////////////////////////////////////////////////////
ElasticMechMaterialStatus :: ElasticMechMaterialStatus(ElasticMechMaterial *m, Element *e, unsigned ipnum) : MaterialStatus(m, e, ipnum) {
    name = "tensorial mechanical mat. status";
}

//////////////////////////////////////////////////////////
void ElasticMechMaterial :: readFromLine(istringstream &iss) {
    string param;
    bool bE, bnu, bdensity;
    bE = bnu = bdensity = false;

    while ( !iss.eof() ) {
        iss >> param;
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
MaterialStatus *ElasticMechMaterial :: giveNewMaterialStatus(Element *e, unsigned ipnum) {
    ElasticMechMaterialStatus *newStatus = new ElasticMechMaterialStatus(this, e, ipnum); //needs to be deleted manually
    return newStatus;
};

//////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////
// ELASTIC COSSERAT MECHANICAL MATERIAL

Vector CosseratMechMaterialStatus ::  giveStress(const Vector &strain, double timeStep) {
    return CosseratMechMaterialStatus ::  giveStressWithFrozenIntVars(strain, timeStep);
};

//////////////////////////////////////////////////////////
Vector CosseratMechMaterialStatus ::  giveStressWithFrozenIntVars(const Vector &strain, double timeStep) {
    ( void ) timeStep;
    temp_strain = strain;
    unsigned dim;
    if ( strain.size() == 2 ) {
        dim = 1;
    } else if ( strain.size() == 6 ) {
        dim = 2;
    } else if ( strain.size() == 18 ) {
        dim = 3;
    } else {
        cerr << name << " error: unsupported dimension" << endl;
        exit(1);
    }
    temp_stress = giveStiffnessTensor("elastic", dim) * addEigenStrain(strain);
    return temp_stress;
};

//////////////////////////////////////////////////////////
Matrix CosseratMechMaterialStatus :: giveStiffnessTensor(string type, unsigned dimension) const {
    ( void ) type;
    unsigned size;
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
    Matrix D(size, size);
    CosseratMechMaterial *m = static_cast< CosseratMechMaterial * >( mat );
    double lammeL = m->giveElasticModulus() * m->givePoissonsRatio() / ( ( 1. + m->givePoissonsRatio() ) * ( 1. - 2 * m->givePoissonsRatio() ) );
    double lammeM = m->giveElasticModulus() / ( 2. * ( 1. + m->givePoissonsRatio() ) );
    if ( dimension == 1 ) {
        D [ 0 ] [ 0 ] = m->giveElasticModulus();
        D [ 1 ] [ 1 ] = 4. * m->giveCosseratShearParam() * pow(m->giveCharacteristicLength(), 2);
    } else if ( dimension == 2 ) {
        if ( m->isPlaneStress() ) {
            cerr << name << " error: dimension " << dimension << ", Cosserat plane stress not implemented" << endl;
            exit(1);
        } else {  //plane strain ACTUALLY NOT IMPLEMENTED YET
            D [ 0 ] [ 0 ] = D [ 1 ] [ 1 ] = lammeL + 2. * lammeM;
            D [ 0 ] [ 1 ] = D [ 1 ] [ 0 ] = lammeL;
            D [ 2 ] [ 2 ] = D [ 3 ] [ 3 ] = lammeM + m->giveCosseratShearParam();
            D [ 2 ] [ 3 ] = D [ 3 ] [ 2 ] = lammeM - m->giveCosseratShearParam();
            D [ 4 ] [ 4 ] = D [ 5 ] [ 5 ] = lammeM * 4. * m->giveCharacteristicLength();
        }
    } else {
        cerr << name << " error: dimension " << dimension << " not implemented" << endl;
        exit(1);
    }
    return D;
};

//////////////////////////////////////////////////////////
CosseratMechMaterialStatus :: CosseratMechMaterialStatus(CosseratMechMaterial *m, Element *e, unsigned ipnum) : ElasticMechMaterialStatus(m, e, ipnum) {
    name = "tensorial cosserat mat. status";
}

//////////////////////////////////////////////////////////
void CosseratMechMaterial :: readFromLine(istringstream &iss) {
    ElasticMechMaterial :: readFromLine(iss);

    iss.clear(); // clear string stream
    iss.seekg(0, iss.beg); //reset position in string stream

    string param;
    bool blc, bmuc;
    blc = bmuc = false;

    while ( !iss.eof() ) {
        iss >> param;
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
MaterialStatus *CosseratMechMaterial :: giveNewMaterialStatus(Element *e, unsigned ipnum) {
    CosseratMechMaterialStatus *newStatus = new CosseratMechMaterialStatus(this, e, ipnum);
    return newStatus;
};


//////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////
// COUPLED TRANSPORT MATERIAL
//////////////////////////////////////////////////////////

DiscreteTrsprtCoupledMaterialStatus :: DiscreteTrsprtCoupledMaterialStatus(TrsprtMaterial *m, Element *e, unsigned ipnum) : DiscreteTrsprtMaterialStatus(m, e, ipnum) {
    name = "discrete coupled transport mat. status";
}

//////////////////////////////////////////////////////////
double DiscreteTrsprtCoupledMaterialStatus :: giveEffectiveConductivity(string type) const {
    if ( type.compare("elastic") == 0  ) {
        TrsprtMaterial *tmat = static_cast< TrsprtMaterial * >( mat );
        return calculatePressureDependentPermeability(0.) * tmat->giveDensity() / tmat->giveViscosity();
    } else if ( type.compare("secant") == 0 || type.compare("unloading") == 0 || type.compare("tangent") == 0 ) {
        return updateEffectiveConductivity();
    } else {
        cerr << "Error: DiscreteTrsprtCoupledMaterialStatus does not provide '" << type << "' stiffness";
        exit(1);
    };
}


//////////////////////////////////////////////////////////
double DiscreteTrsprtCoupledMaterialStatus :: updateEffectiveConductivity() const {
    DiscreteTrsprtCoupledMaterial *tmat = static_cast< DiscreteTrsprtCoupledMaterial * >( mat );
    Transp1DCoupled *tc = static_cast< Transp1DCoupled * >( element );  
    return (TrsprtMaterialStatus::updateEffectiveConductivity()) + tmat->giveTurtuosity() * tmat->giveDensity() / ( 12. * tmat->giveViscosity() * tc->giveArea() ) * crackParam;
}

//////////////////////////////////////////////////////////
Vector DiscreteTrsprtCoupledMaterialStatus :: giveStress(const Vector &strain, double timeStep) {
    if ( timeStep > 0 ) {
        volStrainRate = ( tempVolumetricStrain - volumetricStrain ) / timeStep;
    } else {
        volStrainRate = 0;
    }

    effConductivity = updateEffectiveConductivity();
    temp_strain = strain;
    temp_stress = -effConductivity *addEigenStrain(temp_strain);
    return temp_stress;
};


//////////////////////////////////////////////////////////
Vector DiscreteTrsprtCoupledMaterialStatus :: giveInternalSource() const {
    Vector ints(1);
    DiscreteTrsprtCoupledMaterial *m = static_cast< DiscreteTrsprtCoupledMaterial * >( mat );
    ints [ 0 ] = -m->giveBiotCoeff() * m->giveDensity() * 3. * volStrainRate; //Biot coeff times volumetric strain rate
    return ints;
}


//////////////////////////////////////////////////////////
void DiscreteTrsprtCoupledMaterialStatus :: setParameterValue(string code, double value) {
    if ( code.compare("volumetric_strain") == 0 ) {
        tempVolumetricStrain = value;
    } else if ( code.compare("crack_opening") == 0 ) {
        crackParam = value;
    } else {
        TrsprtMaterialStatus :: setParameterValue(code, value);
    }
}


//////////////////////////////////////////////////////////
Vector DiscreteTrsprtCoupledMaterialStatus :: giveStressWithFrozenIntVars(const Vector &strain, double timeStep) {
    ( void ) timeStep;
    temp_strain.resize(1);
    temp_strain [ 0 ] = strain [ 0 ];
    temp_stress = -effConductivity *addEigenStrain(strain);

    return temp_stress;
};


//////////////////////////////////////////////////////////
void DiscreteTrsprtCoupledMaterialStatus ::  update() {
    TrsprtMaterialStatus :: update();
    volumetricStrain = tempVolumetricStrain;
}

//////////////////////////////////////////////////////////
double DiscreteTrsprtCoupledMaterialStatus ::  giveValue(string code) const {
    if ( code.compare("volumetric_strain") == 0 ) {
        return volumetricStrain;
    } else {
        return TrsprtMaterialStatus :: giveValue(code);
    }
}

//////////////////////////////////////////////////////////
void DiscreteTrsprtCoupledMaterial :: readFromLine(istringstream &iss) {
    TrsprtMaterial :: readFromLine(iss);

    iss.clear(); // clear string stream
    iss.seekg(0, iss.beg); //reset position in string stream

    string param;
    bool bturtuosity, bbiot;
    bturtuosity = bbiot = false;

    while ( !iss.eof() ) {
        iss >> param;
        if ( param.compare("crack_turtuosity") == 0 ) {
            bturtuosity = true;
            iss >> crack_turtuosity;
        } else if ( param.compare("biot_coeff") == 0 ) {
            bbiot = true;
            iss >> biotCoeff;
        }
    }
    if ( !bturtuosity ) {
        cerr << name << ": material parameter 'crack_turtuosity' was not specified" << endl;
        exit(EXIT_FAILURE);
    }
    if ( !bbiot ) {
        cerr << name << ": material parameter 'biot_coeff' was not specified" << endl;
        exit(EXIT_FAILURE);
    }
};

//////////////////////////////////////////////////////////
MaterialStatus *DiscreteTrsprtCoupledMaterial :: giveNewMaterialStatus(Element *e, unsigned ipnum) {
    DiscreteTrsprtCoupledMaterialStatus *newStatus = new DiscreteTrsprtCoupledMaterialStatus(this, e, ipnum); //needs to be deleted manually
    return newStatus;
};

//////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////
// DISCRETE MECHANICAL ELASTIC MATERIAL

DisMechMaterialStatus :: DisMechMaterialStatus(DisMechMaterial *m, Element *e, unsigned ipnum) : MaterialStatus(m, e, ipnum) {
    name = "discrete mechanical mat. status";
    mat = m;
}

//////////////////////////////////////////////////////////
Matrix DisMechMaterialStatus :: giveStiffnessTensor(string type, unsigned dimension) const {
    ( void ) type;
    DisMechMaterial *m = static_cast< DisMechMaterial * >( mat );
    Matrix D(dimension, dimension);
    D [ 0 ] [ 0 ] = m->giveE0();
    for ( size_t i = 1; i < dimension; i++ ) {
        D [ i ] [ i ] =  m->giveAlpha() * m->giveE0();
    }
    return D;
}

//////////////////////////////////////////////////////////
double DisMechMaterialStatus :: giveDensity() const {
    DisMechMaterial *tmat = static_cast< DisMechMaterial * >( mat );
    return tmat->giveDensity();
}

//////////////////////////////////////////////////////////
Vector DisMechMaterialStatus ::  giveStress(const Vector &strain, double timeStep) {
    return DisMechMaterialStatus :: giveStressWithFrozenIntVars(strain, timeStep);
};

//////////////////////////////////////////////////////////
Vector DisMechMaterialStatus ::  giveStressWithFrozenIntVars(const Vector &strain, double timeStep) {
    ( void ) timeStep;
    temp_strain = strain;
    DisMechMaterial *m = static_cast< DisMechMaterial * >( mat );
    temp_stress.resize(strain.size() );
    Vector activeStrain = addEigenStrain(strain);
    temp_stress [ 0 ] = m->giveE0() * activeStrain [ 0 ];
    for ( unsigned i = 1; i < strain.size(); i++ ) {
        temp_stress [ i ] = m->giveAlpha() * m->giveE0() * activeStrain [ i ];
    }
    return temp_stress;
};

//////////////////////////////////////////////////////////
double DisMechMaterialStatus ::  giveValue(string code) const {
    if ( code.compare("s_N") == 0 ) {
        return temp_stress [ 0 ];
    } else if ( code.compare("s_M") == 0 ) {
        return temp_stress [ 1 ];
    } else if ( code.compare("s_L") == 0 ) {
        return temp_stress [ 2 ];
    } else if ( code.compare("e_N") == 0 ) {
        return temp_strain [ 0 ];
    } else if ( code.compare("e_M") == 0 ) {
        return temp_strain [ 1 ];
    } else if ( code.compare("e_L") == 0 ) {
        return temp_strain [ 2 ];
    } else {
        return MaterialStatus :: giveValue(code);
    }
}

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
MaterialStatus *DisMechMaterial :: giveNewMaterialStatus(Element *e, unsigned ipnum) {
    DisMechMaterialStatus *newStatus = new DisMechMaterialStatus(this, e, ipnum); //needs to be deleted manually
    return newStatus;
};
