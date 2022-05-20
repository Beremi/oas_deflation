#include "material_mars.h"
#include "element_discrete.h"

using namespace std;

//////////////////////////////////////////////////////////
// CUSATIS MATERIAL STATUS

MarsMaterialStatus :: MarsMaterialStatus(MarsMaterial *m, Element *e, unsigned ipnum) : DisMechMaterialStatus(m, e, ipnum) {
    name = "MARS mat. status";
    RAND_H = 1.0;
}

//////////////////////////////////////////////////////////
void MarsMaterialStatus :: giveValues(string code, Vector &result) const {
    if ( code.compare("tempCrackOpening") == 0 || code.compare("crack_opening") == 0 ) {
        result.resize(1);
        result [ 0 ] = temp_crackOpening;
    } else if ( code.compare("volumetric_strain") == 0 ) {
        result.resize(1);
        result [ 0 ] = volumetricStrain;
    } else if ( code.rfind("damage", 0) == 0 || code.rfind("damageN", 0) == 0 || code.rfind("damageT", 0) == 0 ) {
        result.resize(1);
        result [ 0 ] = temp_damage;
    } else if ( code.compare("ft") == 0 ) {
        MarsMaterial *m = static_cast< MarsMaterial * >( mat );
        result.resize(1);
        result [ 0 ] = m->giveFt();
    } else if ( code.compare("Gt") == 0 ) {
        MarsMaterial *m = static_cast< MarsMaterial * >( mat );
        result.resize(1);
        result [ 0 ] = m->giveGt();
    } else if ( code.compare("fs") == 0 ) {
        MarsMaterial *m = static_cast< MarsMaterial * >( mat );
        result.resize(1);
        result [ 0 ] = m->giveFs();
    } else if ( code.compare("Gs") == 0 ) {
        MarsMaterial *m = static_cast< MarsMaterial * >( mat );
        result.resize(1);
        result [ 0 ] = m->giveGs();
    } else if ( ( code.compare("strainN") == 0 ) ) {
        result.resize(1);
        result [ 0 ] = temp_strain [ 0 ];
    } else if ( ( code.compare("stressN") == 0 ) ) {
        result.resize(1);
        result [ 0 ] = temp_stress [ 0 ];
    } else {
        DisMechMaterialStatus :: giveValues(code, result);
    }
}

//////////////////////////////////////////////////////////
void MarsMaterialStatus :: init() {
    computeOmega0();
    maxEpsT = 0;
    maxEpsN = 0;
    damage = 0;
    temp_damage = 0;
    temp_maxEpsN = 0;
    temp_maxEpsT = 0;
    temp_crackOpening = 0;
    volumetricStrain = 0;


    crackOpening = 0;

    RigidBodyContact *rbc = dynamic_cast< RigidBodyContact * >( element );
    if ( rbc ) {
        L = rbc->giveLength();
    } else {
        cerr << "Material " << name << " can be used only for RigidBodyContact elements" << endl;
        exit(EXIT_FAILURE);
    }

    MarsMaterial *m = static_cast< MarsMaterial * >( mat );
    Ks = 2 * m->giveAlpha() * m->giveE0() / ( m->giveLcrs() / L - 1 );
    Kt = 2 * m->giveE0() / ( m->giveLcrt() / L - 1 );
    nt = log( Kt / ( Kt - Ks ) ) / log(1 - 2 * omega0 / M_PI);

    if ( Ks < 0 || Kt < 0 ) {
        cerr << "Error " << name << ": snap back occured" << endl;
        exit(1);
    }

}

//////////////////////////////////////////////////////////
double MarsMaterialStatus :: giveS0tension(double omega) const {
    MarsMaterial *m = static_cast< MarsMaterial * >( mat );

    double s = sin(omega);
    double s2 = s * s;
    double c2 = cos(omega) * cos(omega);
    double ft = m->giveFt() * RAND_H;
    double fs = m->giveFs() * RAND_H;

    double sa = .5 * ft * ( pow(fs / ( m->giveMu() * ft ), 2) - 1. );


    if ( omega == atan( sqrt( m->giveAlpha() ) / m->giveMu() ) ) {
        // for this anle, the later equation is undetermined, but hyperbola eq. gives this (see Two Scale Study - Cusatis 2007 doi.org/10.1016/j.engfracmech.2006.01.021)
        return .5 * ( ft + 2 * sa ) * ft / ( ( ft + sa ) * s );
    } else {
        return ( -( ft + sa ) * s + sqrt(pow( ( ft + sa ) * s, 2) + ( m->giveAlpha() * ( c2 / pow(m->giveMu(), 2) ) - s2 ) * ( ft + 2 * sa ) * ft) ) / ( m->giveAlpha() * ( c2 / pow(m->giveMu(), 2) ) - s2 );
    }
}

//////////////////////////////////////////////////////////
double MarsMaterialStatus :: giveS0compression(double omega) const {
    MarsMaterial *m = static_cast< MarsMaterial * >( mat );

    double fc = m->giveFc() * RAND_H;

    return fc / sqrt( pow(sin(omega), 2) + ( m->giveAlpha() * pow(cos(omega), 2) ) / m->giveBeta() );
}

//////////////////////////////////////////////////////////
void MarsMaterialStatus :: computeOmega0() {
    double o1 = 0.1;
    double o2 = -M_PI;
    omega0 = ( o1 + o2 ) / 2;
    double err = giveS0tension(omega0) - giveS0compression(omega0);
    unsigned itmax = 1000;
    unsigned it = 0;
    while ( fabs(err) > 1e-5 && it < itmax ) {
        if ( err * o1 < 0. ) {
            o1 = omega0;
        } else {
            o2 = omega0;
        }
        omega0 = ( o1 + o2 ) / 2.;
        err = giveS0tension(omega0) - giveS0compression(omega0);
        it++;
    }
    if ( it == itmax ) {
        cerr << "Mars Material Error: omega0 does not converging, error " << fabs(err) << endl;
        exit(1);
    }

    if ( omega0 > 0.0 ) {
        omega0 -= 2. * M_PI;
    }
    if ( omega0 < -0.5 * M_PI ) {
        omega0 = -M_PI - omega0;
    }
}

//////////////////////////////////////////////////////////
void MarsMaterialStatus :: computeDamage(Vector strain) {
    MarsMaterial *m = static_cast< MarsMaterial * >( mat );
    double epsN = strain [ 0 ];
    double epsT;
    if ( strain.size() == 2 ) {
        epsT = abs(strain [ 1 ]);                //2D
    } else {
        epsT = sqrt( pow(strain [ 1 ], 2) + pow(strain [ 2 ], 2) );   //3D
    }
    double epsEQ = sqrt( pow(epsN, 2) + m->giveAlpha() * pow(epsT, 2) );   //equivalent strain

    if ( epsEQ > 0 && damage < 1.0 ) {
        double omega, S0, chi, K0, strEQ;
        if ( epsT > 0 ) {
            omega = atan( epsN / ( sqrt( m->giveAlpha() ) * epsT ) );
        } else if ( epsN > 0 ) {
            omega = 0.5 * M_PI;
        } else {
            omega = -0.5 * M_PI;
        }

        if ( omega < omega0 ) { //compression
            S0 = giveS0compression(omega);
            chi = epsEQ;
            K0 = m->giveKc() * ( 1. - pow( ( omega + 0.5 * M_PI ) / ( omega0 + 0.5 * M_PI ), m->giveNc() ) );
        } else { //tension-shear
            double emax, f;
            temp_maxEpsN = max(maxEpsN, epsN);
            temp_maxEpsT = max(maxEpsT, epsT);
            emax = sqrt( pow(temp_maxEpsN, 2) + m->giveAlpha() * pow(temp_maxEpsT, 2) );
            S0 = giveS0tension(omega);

            //confinement not applied
            //if (fabs(confiningStrain) / cLaw.lambda < 1. || confiningStrain > 0.) f = 1.;
            //else f = 1. / (1. - confiningStrain / cLaw.lambda);
            f = 1;


            K0 = -f * Kt * ( 1. - pow( ( omega - 0.5 * M_PI ) / ( omega0 - 0.5 * M_PI ), nt) );
            if ( omega < 0.0 ) {
                chi = epsEQ * omega / omega0 + emax * ( 1. - omega / omega0 );
            } else {
                chi = emax;
            }
        }
        if ( chi - S0 / m->giveE0() > 0 ) {
            strEQ = S0 * exp( K0 / S0 * ( chi - S0 / m->giveE0() ) );
        } else {
            strEQ = S0;
        }

        temp_damage = 1. - strEQ / ( m->giveE0() * epsEQ );
    } else {
        temp_damage = 0.0;
    }

    // no recovery
    if ( temp_damage < damage ) {
        temp_damage = damage;
    }

    //temp_damage = min( temp_damage, m->giveMaxDamage() ); //dangerous, better switched off

    // compression recovery
    //if (temp_damage < damage && epsN>0)  temp_damage = damage;
    //if (temp_damage < damageC && epsN<0) temp_damage = damageC;
}


//////////////////////////////////////////////////////////
void MarsMaterialStatus :: update() {
    DisMechMaterialStatus :: update();
    damage = temp_damage;
    maxEpsN = temp_maxEpsN;
    maxEpsT = temp_maxEpsT;

    crackOpening = temp_crackOpening;
}

//////////////////////////////////////////////////////////
void MarsMaterialStatus :: resetTemporaryVariables() {
    DisMechMaterialStatus :: resetTemporaryVariables();
    temp_damage = damage;
    temp_maxEpsN = maxEpsN;
    temp_maxEpsT = maxEpsT;
    temp_crackOpening = crackOpening;
}



//////////////////////////////////////////////////////////
Matrix MarsMaterialStatus :: giveStiffnessTensor(string type, unsigned dim) const {
    Matrix stiff = DisMechMaterialStatus :: giveStiffnessTensor(type, dim);
    if ( type.compare("elastic") == 0 ) {
        return stiff;
    } else if ( type.compare("secant") == 0 ) {
        MarsMaterial *m = static_cast< MarsMaterial * >( mat );
        if ( m->giveDamageResiduum() > 0.0 ) {
            return stiff * fmax(1 - temp_damage, m->giveDamageResiduum() );
        } else if ( m->giveStressResiduum() > 0.0 ) {
            // TODO finish this JK
            // QUESTION is this performed before update?
            if ( temp_damage > damage ) { // if damage increases, apply residuum
                double sN, sT;
                sN = temp_stress [ 0 ];
                sT = 0.0;
                for ( unsigned i = 1; i < temp_stress.size(); i++ ) {
                    sT += pow(temp_stress [ i ], 2);
                }
                double strs = sqrt( pow(sN, 2) + ( sT / m->giveAlpha() ) );
                if ( strs  < m->giveStressResiduum() ) {
                    double epsN, epsT;
                    epsN = temp_strain [ 0 ];
                    epsT = 0.0;
                    for ( unsigned i = 1; i < temp_strain.size(); i++ ) {
                        epsT += pow(temp_strain [ i ], 2);
                    }
                    double epsEQ = sqrt( pow(epsN, 2) + epsT * m->giveAlpha() );
                    return stiff * ( 1 - m->giveStressResiduum() / ( m->giveE0() * epsEQ ) );
                } else {
                    return stiff * ( 1 - temp_damage );
                }
            } else {
                return stiff * ( 1 - temp_damage );
            }
        } else {
            return stiff * ( 1 - temp_damage );
        }
    } else if ( type.compare("unloading") == 0 ) {
        return stiff * ( 1 - temp_damage );
    } else if ( type.compare("tangent") == 0 ) {
        return stiff * ( 1 - temp_damage );                                   //not implemented, used unloading
    } else {
        cerr << "Error: MarsMaterialStatus does not provide '" << type << "' stiffness";
        exit(1);
    };
}

//////////////////////////////////////////////////////////
Vector MarsMaterialStatus :: giveStress(const Vector &strain, double timeStep) {
    computeDamage( addEigenStrain(strain) );
    return MarsMaterialStatus :: giveStressWithFrozenIntVars(strain, timeStep);
}

//////////////////////////////////////////////////////////
Vector MarsMaterialStatus :: giveStressWithFrozenIntVars(const Vector &strain, double timeStep) {
    temp_stress = DisMechMaterialStatus :: giveStressWithFrozenIntVars(strain, timeStep) * ( 1. - temp_damage );   //without eigen strain, it will be applied later
    if ( temp_strain [ 0 ] > 0 ) {
        temp_crackOpening = ( L * temp_damage ) * temp_strain [ 0 ]; //normal opening only
        //temp_crackOpening = l2_norm( ( L * temp_damage ) * temp_strain );  //total opening
    } else {
        temp_crackOpening = 0;
    }
    return temp_stress;
}

//////////////////////////////////////////////////////////
std :: string MarsMaterialStatus :: giveLineToSave() const {
    return "damage " + to_string_sci(this->damage) + " maxEpsN " + to_string_sci(this->maxEpsN) + " maxEpsT " + to_string_sci(this->maxEpsT);
}

//////////////////////////////////////////////////////////
void MarsMaterialStatus :: setParameterValue(string code, double value) {
    if ( code.compare("volumetric_strain") == 0 ) {
        volumetricStrain = value;
    } else {
        DisMechMaterialStatus :: setParameterValue(code, value);
    }
}

//////////////////////////////////////////////////////////
void MarsMaterialStatus :: readFromLine(istringstream &iss) {
    std :: string param;
    while ( !iss.eof() ) {
        iss >> param;
        if ( param.compare("damage") == 0 ) {
            iss >> this->damage;
            temp_damage = damage;
        } else if ( param.compare("maxEpsN") == 0 ) {
            iss >> this->maxEpsN;
            temp_maxEpsN = maxEpsN;
        } else if ( param.compare("maxEpsT") == 0 ) {
            iss >> this->maxEpsT;
            temp_maxEpsT = maxEpsT;
        }
    }
}

//////////////////////////////////////////////////////////
bool MarsMaterialStatus :: isElastic(const bool &now) const {
    if ( now && this->temp_damage != 0.0 ) {
        return false;
    } else if ( this->damage != 0.0 ) {
        return false;
    }
    return true;
}


//////////////////////////////////////////////////////////
// CUSATIS MATERIAL

//////////////////////////////////////////////////////////
void MarsMaterial :: readFromLine(istringstream &iss) {
    DisMechMaterial :: readFromLine(iss); //read elastic parameters

    iss.clear(); // clear string stream
    iss.seekg(0, iss.beg); //reset position in string stream

    // initialize all values to zero (NOTE probably no ned in linux, but in windows necessary)
    fs = Gs = fc = Kc = 0;

    string param;
    bool bft, bGt;
    bft = bGt = false;

    while ( !iss.eof() ) {
        iss >> param;
        if ( param.compare("Gt") == 0 ) {
            bGt = true;
            iss >> Gt;
        } else if ( param.compare("ft") == 0 ) {
            bft = true;
            iss >> ft;
        } else if ( param.compare("fs") == 0 ) {
            iss >> fs;
        } else if ( param.compare("Gs") == 0 ) {
            iss >> Gs;
        } else if ( param.compare("fc") == 0 ) {
            iss >> fc;
        } else if ( param.compare("Kc") == 0 ) {
            iss >> Kc;
        } else if ( param.compare("damage_residuum") == 0 ) {
            iss >> damage_residuum;
        } else if ( param.compare("stress_residuum_fraction") == 0 ) {
            iss >> stress_residuum_fraction;
        }
    }
    if ( !bft ) {
        cerr << name << ": material parameter 'ft' was not specified" << endl;
        exit(EXIT_FAILURE);
    }
    ;
    if ( !bGt ) {
        cerr << name << ": material parameter 'Gt' was not specified" << endl;
        exit(EXIT_FAILURE);
    }
    ;
};

//////////////////////////////////////////////////////////
MaterialStatus *MarsMaterial :: giveNewMaterialStatus(Element *e, unsigned ipnum) {
    MarsMaterialStatus *newStatus = new MarsMaterialStatus(this, e, ipnum); //needs to be deleted manually
    return newStatus;
};


//////////////////////////////////////////////////////////
void MarsMaterial :: init() {
    // if variables not specified on the input, use default multipliers
    fs = ( fs == 0 ) ? 3 * ft : fs;
    Gs = ( Gs == 0 ) ? 16 * Gt : Gs;
    fc = ( fc == 0 ) ? 16 * ft : fc;
    Kc = ( Kc == 0 ) ? 0.26 * E0 : Kc;
    beta = 1.;
    mu = 0.2;
    nc = 2;

    Lcrt = 2 * E0 * Gt / pow(ft, 2);
    Lcrs = 2 * alpha * E0 * Gs / pow(fs, 2);
};

//////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////
// COUPLED CUSATIS MATERIAL STATUS
//////////////////////////////////////////////////////////
CoupledMarsMaterialStatus :: CoupledMarsMaterialStatus(MarsMaterial *m, Element *e, unsigned ipnum) : MarsMaterialStatus(m, e, ipnum) {
    name = "Coupled MARS mat. status";
}

//////////////////////////////////////////////////////////
void CoupledMarsMaterialStatus :: init() {
    MarsMaterialStatus :: init();
    avgPressure = 0;
}


//////////////////////////////////////////////////////////
void CoupledMarsMaterialStatus :: setParameterValue(string code, double value) {
    if ( code.compare("pressure") == 0 ) {
        avgPressure = value;
    } else {
        MarsMaterialStatus :: setParameterValue(code, value);
    }
}

//////////////////////////////////////////////////////////
void CoupledMarsMaterialStatus :: giveValues(string code, Vector &result) const {
    if ( code.compare("pressure") == 0 || code.compare("avg_pressure") == 0 ) {
        result.resize(1);
        result [ 0 ] = avgPressure;
    } else if ( code.compare("solid_stress") == 0 ) {
        MarsMaterialStatus :: giveValues("stress", result); //standard stress including Biot's effect
        CoupledMarsMaterial *m = static_cast< CoupledMarsMaterial * >( mat );
        result [ 0 ] += m->giveBiotCoeff() * avgPressure; //stress without Biot's effect
    }
    else return MarsMaterialStatus :: giveValues(code, result);
}



//////////////////////////////////////////////////////////
void CoupledMarsMaterialStatus :: updateStressByBiotEffect(double timeStep) {
    ( void ) timeStep;
    CoupledMarsMaterial *m = static_cast< CoupledMarsMaterial * >( mat );
    temp_stress [ 0 ] -= m->giveBiotCoeff() * avgPressure;
}

//////////////////////////////////////////////////////////
Vector CoupledMarsMaterialStatus :: giveStress(const Vector &strain, double timeStep) {
    MarsMaterialStatus :: giveStress(strain, timeStep);
    updateStressByBiotEffect(timeStep);
    return temp_stress;
}

//////////////////////////////////////////////////////////
Vector CoupledMarsMaterialStatus :: giveStressWithFrozenIntVars(const Vector &strain, double timeStep) {
    MarsMaterialStatus :: giveStressWithFrozenIntVars(strain, timeStep);
    updateStressByBiotEffect(timeStep);
    return temp_stress;
}

//////////////////////////////////////////////////////////
void CoupledMarsMaterialStatus :: update() {
    MarsMaterialStatus :: update();
}

//////////////////////////////////////////////////////////
void CoupledMarsMaterialStatus :: resetTemporaryVariables() {
    MarsMaterialStatus :: resetTemporaryVariables();
}

//////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////
// COUPLED CUSATIS MATERIAL
//////////////////////////////////////////////////////////
void CoupledMarsMaterial :: init() {
    MarsMaterial :: init();
}

//////////////////////////////////////////////////////////
MaterialStatus *CoupledMarsMaterial :: giveNewMaterialStatus(Element *e, unsigned ipnum) {
    CoupledMarsMaterialStatus *newStatus = new CoupledMarsMaterialStatus(this, e, ipnum);
    return newStatus;
};

//////////////////////////////////////////////////////////
void CoupledMarsMaterial :: readFromLine(istringstream &iss) {
    MarsMaterial :: readFromLine(iss);

    iss.clear(); // clear string stream
    iss.seekg(0, iss.beg); //reset position in string stream

    string param;
    bool bbiot = false;

    while ( !iss.eof() ) {
        iss >> param;
        if ( param.compare("biot_coeff") == 0 ) {
            bbiot = true;
            iss >> biotCoeff;
        }
    }
    if ( !bbiot ) {
        cerr << name << ": material parameter 'biot_coeff' was not specified" << endl;
        exit(EXIT_FAILURE);
    }
}
