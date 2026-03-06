#include "material_vectorial.h"
#include "element_discrete.h"
#include "material_container.h"

using namespace std;

//////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////
// VECTORIAL TRANSPORT MATERIAL
//////////////////////////////////////////////////////////

VectTrsprtMaterialStatus :: VectTrsprtMaterialStatus(VectTrsprtMaterial *m, Element *e, unsigned ipnum) : TensTrsprtMaterialStatus(m, e, ipnum) {
    name = "transport mat. status";
}

//////////////////////////////////////////////////////////
Matrix VectTrsprtMaterialStatus :: giveStiffnessTensor(string type) const {
    ( void ) type;
    Matrix T = Matrix :: Zero(1, 1); //discrete material, only one direction in any dimension
    T(0, 0) = -giveEffectiveConductivity(type);
    return T;
};

//////////////////////////////////////////////////////////
MaterialStatus *VectTrsprtMaterial :: giveNewMaterialStatus(Element *e, unsigned ipnum) {
    VectTrsprtMaterialStatus *newStatus = new VectTrsprtMaterialStatus(this, e, ipnum); //needs to be deleted manually
    return newStatus;
};

//////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////
// VECTORIAL HEAT CONDUCTION MATERIAL MATERIAL
//////////////////////////////////////////////////////////

VectHeatConductionMaterialStatus :: VectHeatConductionMaterialStatus(VectHeatConductionMaterial *m, Element *e, unsigned ipnum) : TensHeatConductionMaterialStatus(m, e, ipnum) {
    name = "vectorial heat conduction mat. status";
}

//////////////////////////////////////////////////////////
Matrix VectHeatConductionMaterialStatus :: giveStiffnessTensor(string type) const {
    ( void ) type;
    VectHeatConductionMaterial *hcm = static_cast< VectHeatConductionMaterial * >( mat );
    Matrix T = Matrix :: Zero(1, 1); //discrete material, only one direction in any dimension
    T(0, 0) = hcm->giveConductivity();
    return T;
};

//////////////////////////////////////////////////////////
MaterialStatus *VectHeatConductionMaterial :: giveNewMaterialStatus(Element *e, unsigned ipnum) {
    VectHeatConductionMaterialStatus *newStatus = new VectHeatConductionMaterialStatus(this, e, ipnum); //needs to be deleted manually
    return newStatus;
};


//////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////
// COUPLED TRANSPORT MATERIAL
//////////////////////////////////////////////////////////

VectTrsprtCoupledMaterialStatus :: VectTrsprtCoupledMaterialStatus(VectTrsprtCoupledMaterial *m, Element *e, unsigned ipnum) : VectTrsprtMaterialStatus(m, e, ipnum) {
    name = "discrete coupled transport mat. status";
    crackParam = 0.;
    temperature = 0;
    temp_volumetricStrain = volumetricStrain = volStrainRate = 0.;
    temp_crackVolume = crackVolume = crackVolumeRate = 0.;
    pressure = pressureRate = 0.;
}

//////////////////////////////////////////////////////////
double VectTrsprtCoupledMaterialStatus :: giveEffectiveConductivity(string type) const {
    if ( type.compare("elastic") == 0  ) {
        TensTrsprtMaterial *tmat = static_cast< TensTrsprtMaterial * >( mat );
        return calculatePressureDependentPermeability(0.) * tmat->giveDensity() / tmat->giveViscosity();
    } else if ( type.compare("secant") == 0 || type.compare("unloading") == 0 || type.compare("tangent") == 0 ) {
        return updateEffectiveConductivity();
    } else {
        cerr << "Error: VectTrsprtCoupledMaterialStatus does not provide '" << type << "' stiffness";
        exit(1);
    };
}


//////////////////////////////////////////////////////////
double VectTrsprtCoupledMaterialStatus :: updateEffectiveConductivity() const {
    VectTrsprtCoupledMaterial *tmat = static_cast< VectTrsprtCoupledMaterial * >( mat );
    DiscreteTrsprtCoupledElem *tc = static_cast< DiscreteTrsprtCoupledElem * >( element );
    return ( TensTrsprtMaterialStatus :: updateEffectiveConductivity() ) + tmat->giveTurtuosity() * tmat->giveDensity() / ( 12. * tmat->giveViscosity() * tc->giveArea() ) * crackParam;
}

//////////////////////////////////////////////////////////
void VectTrsprtCoupledMaterialStatus ::  updateRateVariables(double timeStep) {
    if ( timeStep > 0 ) {
        volStrainRate = ( temp_volumetricStrain - volumetricStrain ) / timeStep;
        crackVolumeRate = ( temp_crackVolume - crackVolume ) / timeStep;
        pressureRate = ( temp_pressure - pressure ) / timeStep;
    } else {
        volStrainRate = 0.;
        crackVolumeRate = 0.;
        pressureRate = 0.;
    }
}

//////////////////////////////////////////////////////////
void VectTrsprtCoupledMaterialStatus :: computeStress(double timeStep) {
    updateRateVariables(timeStep);
    effConductivity = updateEffectiveConductivity();
    computeStressWithFrozenIntVars(timeStep);
};


//////////////////////////////////////////////////////////
Vector VectTrsprtCoupledMaterialStatus :: giveInternalSource() const {
    Vector ints = Vector :: Zero(1);
    VectTrsprtCoupledMaterial *m = static_cast< VectTrsprtCoupledMaterial * >( mat );


    ints [ 0 ]  = -m->giveBiotCoeff() *  3. * volStrainRate; //Biot coeff times volumetric strain rate
    if ( crackVolumeRate > 0 || pressureRate > 0 ) {
        DiscreteTrsprtElem *trs = static_cast< DiscreteTrsprtElem * >( element );
        double vol = trs->giveVolume();
        if ( temp_crackVolume > 0 ) {
            ints [ 0 ] -= temp_crackVolume * pressureRate / ( vol * m->giveKw() );
        }
        ints [ 0 ] -= crackVolumeRate / vol * ( 1. - m->giveBiotCoeff() +  ( temp_pressure - m->giveReferencePressure() ) / m->giveKw() );
    }
    return ints * m->giveDensity();
}


//////////////////////////////////////////////////////////
void VectTrsprtCoupledMaterialStatus :: setParameterValue(string code, double value) {
    if ( code.compare("volumetric_strain") == 0 ) {
        temp_volumetricStrain = value;
    } else if ( code.compare("crack_param") == 0 ) {
        crackParam = value;
    } else if ( code.compare("crack_volume") == 0 ) {
        temp_crackVolume = value;
    } else if ( code.compare("temperature") == 0 ) {
        temperature = value;
    } else {
        VectTrsprtMaterialStatus :: setParameterValue(code, value);
    }
}


//////////////////////////////////////////////////////////
void VectTrsprtCoupledMaterialStatus :: computeStressWithFrozenIntVars(double timeStep) {
    updateRateVariables(timeStep);
    computeConstitutiveStrain();
    temp_stress = -effConductivity * temp_strain;
};



//////////////////////////////////////////////////////////
void VectTrsprtCoupledMaterialStatus ::  update() {
    TensTrsprtMaterialStatus :: update();
    volumetricStrain = temp_volumetricStrain;
    crackVolume = temp_crackVolume;
    pressure = temp_pressure;
}

//////////////////////////////////////////////////////////
void VectTrsprtCoupledMaterialStatus ::  resetTemporaryVariables() {
    TensTrsprtMaterialStatus :: resetTemporaryVariables();
    temp_volumetricStrain = volumetricStrain;
    temp_crackVolume = crackVolume;
    temp_pressure = pressure;
}



//////////////////////////////////////////////////////////
bool VectTrsprtCoupledMaterialStatus ::  giveValues(string code, Vector &result) const {
    if ( code.compare("volumetric_strain") == 0 ) {
        result.resize(1);
        result [ 0 ] = volumetricStrain;
        return true;
    } else if ( code.compare("crack_volume") == 0 ) {
        result.resize(1);
        result [ 0 ] = temp_crackVolume;
        return true;
    } else if ( code.compare("crack_volume_rate") == 0 ) {
        result.resize(1);
        result [ 0 ] = crackVolumeRate;
        return true;
    } else if ( code.compare("temperature") == 0 ) {
        result.resize(1);
        result [ 0 ] = temperature;
        return true;
    } else if ( code.compare("rel_crack_volume_rate") == 0 ) {
        DiscreteTrsprtElem *trs = static_cast< DiscreteTrsprtElem * >( element );
        double vol = trs->giveVolume();
        result.resize(1);
        result [ 0 ] = crackVolumeRate / vol;
        return true;
    } else {
        return TensTrsprtMaterialStatus :: giveValues(code, result);
    }
}

//////////////////////////////////////////////////////////
void VectTrsprtCoupledMaterial :: readFromLine(istringstream &iss) {
    TensTrsprtMaterial :: readFromLine(iss);

    iss.clear(); // clear string stream
    iss.seekg(0, iss.beg); //reset position in string stream

    string param;
    bool bturtuosity, bbiot;
    bturtuosity = bbiot = false;

    while (  iss >> param ) {
        if ( param.compare("crack_turtuosity") == 0 ) {
            bturtuosity = true;
            iss >> crack_turtuosity;
        } else if ( param.compare("biot_coeff") == 0 ) {
            bbiot = true;
            iss >> biotCoeff;
        } else if ( param.compare("reference_pressure") == 0 ) {
            iss >> refP;
        } else if ( param.compare("Kw") == 0 ) {
            iss >> Kw;
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
MaterialStatus *VectTrsprtCoupledMaterial :: giveNewMaterialStatus(Element *e, unsigned ipnum) {
    VectTrsprtCoupledMaterialStatus *newStatus = new VectTrsprtCoupledMaterialStatus(this, e, ipnum); //needs to be deleted manually
    return newStatus;
};

//////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////
// VECTORIAL MECHANICAL ELASTIC MATERIAL
//////////////////////////////////////////////////////////

VectMechMaterialStatus :: VectMechMaterialStatus(VectMechMaterial *m, Element *e, unsigned ipnum) : MaterialStatus(m, e, ipnum) {
    name = "discrete mechanical mat. status";
    mat = m;
    normalEnergyDensity = 0;
    shearEnergyDensity = 0;
    temp_volumetricStrain = 0;
    volumetricStrain = 0;
    temp_volumetricStrain_total = 0;
    volumetricStrain_total = 0;
    eigenVolumetricStrain = 0;
}

//////////////////////////////////////////////////////////
Matrix VectMechMaterialStatus :: giveStiffnessTensor(string type) const {
    ( void ) type;
    unsigned ss = mat->giveStrainSize();
    VectMechMaterial *m = static_cast< VectMechMaterial * >( mat );
    Matrix D = Matrix :: Zero(ss, ss);
    D(0, 0) = m->giveE0();
    for ( size_t i = 1; i < ss; i++ ) {
        D(i, i) =  m->giveAlpha() * m->giveE0();
    }
    return D;
}

//////////////////////////////////////////////////////////
double VectMechMaterialStatus :: giveDensity() const {
    VectMechMaterial *tmat = static_cast< VectMechMaterial * >( mat );
    return tmat->giveDensity();
}

//////////////////////////////////////////////////////////
void VectMechMaterialStatus ::  computeStress(double timeStep) {
    computeStressWithFrozenIntVars(timeStep);
};

//////////////////////////////////////////////////////////
void VectMechMaterialStatus :: setParameterValue(string code, double value) {
    if ( code.compare("volumetric_strain") == 0 ) {
        temp_volumetricStrain_total = value;
    } else {
        MaterialStatus :: setParameterValue(code, value);
    }
}

//////////////////////////////////////////////////////////
void VectMechMaterialStatus ::  update() {
    normalEnergyDensity += ( temp_strain [ 0 ] - updt_strain [ 0 ] ) * ( updt_stress [ 0 ] + temp_stress [ 0 ] ) / 2.;
    shearEnergyDensity += ( temp_strain [ 1 ] - updt_strain [ 1 ] ) * ( updt_stress [ 1 ] + temp_stress [ 1 ] ) / 2.;
    if ( mat->giveDimension() == 3 ) {
        shearEnergyDensity += ( temp_strain [ 2 ] - updt_strain [ 2 ] ) * ( updt_stress [ 2 ] + temp_stress [ 2 ] ) / 2.;
    }

    volumetricStrain = temp_volumetricStrain;
    volumetricStrain_total = temp_volumetricStrain_total;

    MaterialStatus :: update();
}

//////////////////////////////////////////////////////////
void VectMechMaterialStatus ::  computeStressWithFrozenIntVars(double timeStep) {
    ( void ) timeStep;
    computeConstitutiveStrain();
    VectMechMaterial *m = static_cast< VectMechMaterial * >( mat );
    temp_stress.resize( temp_strain.size() );
    temp_stress [ 0 ] = m->giveE0() * temp_strain [ 0 ];
    for ( unsigned i = 1; i < temp_strain.size(); i++ ) {
        temp_stress [ i ] = m->giveAlpha() * m->giveE0() * temp_strain [ i ];
    }
};

//////////////////////////////////////////////////////////
bool VectMechMaterialStatus ::  giveValues(string code, Vector &result) const {
    if ( code.compare("stress") == 0 || code.compare("stresses") == 0 || code.compare("solid_stress") == 0 ) {
        unsigned size = element->giveDimension();
        result.resize(size);
        if ( size > temp_stress.size() ) {
            size = temp_stress.size();
        }
        for ( unsigned p = 0; p < size; p++ ) {
            result [ p ] = temp_stress [ p ];
        }
        return true;
    } else if ( code.compare("strain") == 0 || code.compare("strains") == 0 ) {
        unsigned size = element->giveDimension();
        result.resize(size);
        if ( size > temp_strain.size() ) {
            size = temp_strain.size();
        }
        for ( unsigned p = 0; p < size; p++ ) {
            result [ p ] = temp_strain [ p ];
        }
        return true;
    } else if ( code.compare("volumetric_strain") == 0 ) {
        result.resize(1);
        result [ 0 ] = temp_volumetricStrain;
        return true;
    } else if ( code.compare("total_volumetric_strain") == 0 ) {
        result.resize(1);
        result [ 0 ] = temp_volumetricStrain_total;
        return true;
    } else if ( code.compare("E0") == 0 ) {
        result.resize(1);
        VectMechMaterial *m = static_cast< VectMechMaterial * >( mat );
        result [ 0 ] = m->giveE0();
        return true;
    } else if ( code.compare("s_N") == 0 ) {
        result.resize(1);
        return true;

        result [ 0 ] = temp_stress [ 0 ];
    } else if ( code.compare("s_M") == 0 ) {
        result.resize(1);
        result [ 0 ] = temp_stress [ 1 ];
        return true;
    } else if ( code.compare("s_L") == 0 ) {
        result.resize(1);
        result [ 0 ] = temp_stress [ 2 ];
        return true;
    } else if ( code.compare("e_N") == 0 ) {
        result.resize(1);
        result [ 0 ] = temp_strain [ 0 ];
        return true;
    } else if ( code.compare("e_M") == 0 ) {
        result.resize(1);
        result [ 0 ] = temp_strain [ 1 ];
        return true;
    } else if ( code.compare("e_L") == 0 ) {
        result.resize(1);
        result [ 0 ] = temp_strain [ 2 ];
        return true;
    } else {
        return MaterialStatus :: giveValues(code, result);
    }
}


//////////////////////////////////////////////////////////
void VectMechMaterialStatus :: addToEigenVolumetricStrain(double x) {
    eigenVolumetricStrain += x;
}

//////////////////////////////////////////////////////////
VectMechMaterial :: VectMechMaterial(unsigned dimension) : Material(dimension) {
    name = "Vect mechanical material";
}

//////////////////////////////////////////////////////////
void VectMechMaterial :: readFromLine(istringstream &iss) {
    string param;

    bool bE0, balpha, bdensity;
    bE0 = balpha = bdensity = false;

    while (  iss >> param ) {
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
MaterialStatus *VectMechMaterial :: giveNewMaterialStatus(Element *e, unsigned ipnum) {
    VectMechMaterialStatus *newStatus = new VectMechMaterialStatus(this, e, ipnum); //needs to be deleted manually
    return newStatus;
};

//////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////
// VECTORIAL MECHANICAL ELASTIC MATERIAL WITH VOLUMETRIC-DEVIATORIC SPLIT
//////////////////////////////////////////////////////////

VectMechVolDevSplitMaterialStatus :: VectMechVolDevSplitMaterialStatus(VectMechVolDevSplitMaterial *m, Element *e, unsigned ipnum) : VectMechMaterialStatus(m, e, ipnum) {
    name = "discrete mechanical mat. status with volumetric-deviatoric split";
    temp_volumetricStrain = 0;
}

//////////////////////////////////////////////////////////
Matrix VectMechVolDevSplitMaterialStatus :: giveStiffnessTensor(string type) const {
    ( void ) type;
    unsigned ss = mat->giveStrainSize();
    VectMechVolDevSplitMaterial *m = static_cast< VectMechVolDevSplitMaterial * >( mat );
    Matrix D = Matrix :: Zero(ss, ss);
    D(0, 0) = m->giveE0();
    for ( size_t i = 1; i < ss; i++ ) {
        D(i, i) =  m->giveE0();
    }
    return D;
}

//////////////////////////////////////////////////////////
void VectMechVolDevSplitMaterialStatus ::  computeStressWithFrozenIntVars(double timeStep) {
    ( void ) timeStep;
    computeConstitutiveStrain();
    VectMechVolDevSplitMaterial *m = static_cast< VectMechVolDevSplitMaterial * >( mat );
    temp_stress.resize( temp_strain.size() );
    double ED = m->giveE0();
    double EV = m->giveAlpha();

    temp_stress [ 0 ] = ED * temp_strain [ 0 ] + ( EV - ED ) * temp_volumetricStrain;
    for ( unsigned i = 1; i < temp_strain.size(); i++ ) {
        temp_stress [ i ] = ED * temp_strain [ i ];
    }
};

//////////////////////////////////////////////////////////
void VectMechVolDevSplitMaterialStatus :: setParameterValue(string code, double value) {
    if ( code.compare("volumetric_strain") == 0 ) {
        temp_volumetricStrain = value;
    } else {
        MaterialStatus :: setParameterValue(code, value);
    }
}

//////////////////////////////////////////////////////////
bool VectMechVolDevSplitMaterialStatus ::  giveValues(string code, Vector &result) const {
    if ( code.compare("volumetric_strain") == 0 ) {
        result.resize(1);
        result [ 0 ] = temp_volumetricStrain;
        return true;
    } else {
        return VectMechMaterialStatus :: giveValues(code, result);
    }
}


//////////////////////////////////////////////////////////
VectMechVolDevSplitMaterial :: VectMechVolDevSplitMaterial(unsigned dimension) : VectMechMaterial(dimension) {
    name = "Vect mechanical material with volumetric-deviatoric split";
}

//////////////////////////////////////////////////////////
void VectMechVolDevSplitMaterial :: readFromLine(istringstream &iss) {
    string param;

    bool bED, bEV, bdensity;
    bED = bEV = bdensity = false;

    while (  iss >> param ) {
        if ( param.compare("ED") == 0 ) {
            bED = true;
            iss >> E0;
        } else if ( param.compare("EV") == 0 ) {
            bEV = true;
            iss >> alpha;
        } else if ( param.compare("density") == 0 ) {
            bdensity = true;
            iss >> density;
        }
    }
    if ( !bED ) {
        cerr << name << ": material parameter 'ED' was not specified" << endl;
        exit(EXIT_FAILURE);
    }
    ;
    if ( !bEV ) {
        cerr << name << ": material parameter 'EV' was not specified" << endl;
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
MaterialStatus *VectMechVolDevSplitMaterial :: giveNewMaterialStatus(Element *e, unsigned ipnum) {
    VectMechVolDevSplitMaterialStatus *newStatus = new VectMechVolDevSplitMaterialStatus(this, e, ipnum); //needs to be deleted manually
    return newStatus;
};

//////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////
// VECTORIAL MECHANICAL ELASTIC MATERIAL WITH ROTATIONAL STIFFNESS
//////////////////////////////////////////////////////////

VectMechMaterialWithRotationalStiffnessStatus :: VectMechMaterialWithRotationalStiffnessStatus(VectMechMaterialWithRotationalStiffness *m, Element *e, unsigned ipnum) : VectMechMaterialStatus(m, e, ipnum) {
    name = "discrete mechanical mat. status with rotational stiffness";
    mat = m;
}

//////////////////////////////////////////////////////////
Matrix VectMechMaterialWithRotationalStiffnessStatus :: giveStiffnessTensor(string type) const {
    ( void ) type;
    unsigned ss = mat->giveStrainSize();
    VectMechMaterialWithRotationalStiffness *m = static_cast< VectMechMaterialWithRotationalStiffness * >( mat );
    RigidBodyContactWithRotationalStiffness *rbcr = static_cast< RigidBodyContactWithRotationalStiffness * >( element );
    double A = rbcr->giveArea();
    double I = rbcr->giveMomentOfInertia();
    Matrix D = Matrix :: Zero(ss, ss);
    D(0, 0) = m->giveE0();
    size_t i = 1;
    for ( ; i < m->giveDimension(); i++ ) {
        D(i, i) =  m->giveAlpha() * m->giveE0();
    }
    for ( ; i < ss; i++ ) {
        D(i, i) =  m->giveBeta() * m->giveE0() * I  * rbcr->giveNumIP()  / A;
    }
    return D;
}

//////////////////////////////////////////////////////////
void VectMechMaterialWithRotationalStiffnessStatus ::  computeStressWithFrozenIntVars(double timeStep) {
    ( void ) timeStep;
    computeConstitutiveStrain();
    VectMechMaterialWithRotationalStiffness *m = static_cast< VectMechMaterialWithRotationalStiffness * >( mat );
    temp_stress.resize( temp_strain.size() );
    temp_stress [ 0 ] = m->giveE0() * temp_strain [ 0 ];
    unsigned dim = m->giveDimension();
    RigidBodyContactWithRotationalStiffness *rbcr = static_cast< RigidBodyContactWithRotationalStiffness * >( element );
    double A = rbcr->giveArea();
    double I = rbcr->giveMomentOfInertia();
    size_t i = 1;
    for ( ; i < dim; i++ ) {
        temp_stress [ i ] = m->giveAlpha() * m->giveE0() * temp_strain [ i ];
    }
    for ( ; i < ( size_t ) temp_strain.size(); i++ ) {
        temp_stress [ i ] =  m->giveBeta() * m->giveE0() * temp_strain [ i ] * I *rbcr->giveNumIP() / A;
    }
}

//////////////////////////////////////////////////////////
bool VectMechMaterialWithRotationalStiffnessStatus ::  giveValues(string code, Vector &result) const {
    return MaterialStatus :: giveValues(code, result);
}

//////////////////////////////////////////////////////////
void VectMechMaterialWithRotationalStiffness :: readFromLine(istringstream &iss) {
    VectMechMaterial :: readFromLine(iss);
    iss.clear(); // clear string stream
    iss.seekg(0, iss.beg); //reset position in string stream

    string param;
    bool bbeta;
    bbeta = false;

    while (  iss >> param ) {
        if ( param.compare("beta") == 0 ) {
            bbeta = true;
            iss >> beta;
        }
    }
    if ( !bbeta ) {
        cerr << name << ": material parameter 'beta' was not specified" << endl;
        exit(EXIT_FAILURE);
    }
};

//////////////////////////////////////////////////////////
MaterialStatus *VectMechMaterialWithRotationalStiffness :: giveNewMaterialStatus(Element *e, unsigned ipnum) {
    VectMechMaterialWithRotationalStiffnessStatus *newStatus = new VectMechMaterialWithRotationalStiffnessStatus(this, e, ipnum); //needs to be deleted manually
    return newStatus;
};
