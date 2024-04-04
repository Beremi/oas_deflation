#include "material_ldpm.h"
#include "element_discrete.h"
#include "element_ldpm.h"

using namespace std;

//////////////////////////////////////////////////////////
// LDPM MATERIAL STATUS

LDPMMaterialStatus :: LDPMMaterialStatus(LDPMMaterial *m, Element *e, unsigned ipnum) : VectMechMaterialStatus(m, e, ipnum) {
    name = "LDPM mat. status";
    RAND_H = 1.0;
}

//////////////////////////////////////////////////////////
void LDPMMaterialStatus :: init() {
    maxEpsT = 0;
    maxEpsN = 0;
    temp_maxEpsN = 0;
    temp_maxEpsT = 0;
    temp_volumetricStrain = 0;
    volumetricStrain = 0;
    crackOpening = 0;
    temp_crackOpening = 0;
    virtual_damage = 0;

    RigidBodyContact *rbc = dynamic_cast< RigidBodyContact * >( element );
    LDPMTetra *tet = dynamic_cast< LDPMTetra * >( element );
    if ( rbc ) {
        L = rbc->giveLength();
    } else if ( tet ) {
        L = tet->giveLength(idx);
    } else {
        cerr << "Material " << name << " can be used only for RigidBodyContact or LDMPTetra elements" << endl;
        exit(EXIT_FAILURE);
    }
}

//////////////////////////////////////////////////////////
void LDPMMaterialStatus :: initializeStressAndStrainVector(unsigned num) {
    updt_mech_strain = Vector :: Zero(num);
    VectMechMaterialStatus :: initializeStressAndStrainVector(num);
}

//////////////////////////////////////////////////////////
bool LDPMMaterialStatus :: giveValues(string code, Vector &result) const {
    if ( code.compare("tempCrackOpening") == 0 || code.compare("crack_opening") == 0 ) {
        result.resize(1);
        result [ 0 ] = temp_crackOpening;
        return true;
    } else if ( code.compare("volumetric_strain") == 0 ) {
        result.resize(1);
        result [ 0 ] = temp_volumetricStrain;
        return true;
    } else if ( code.rfind("virtual_damage") == 0 ) {
        result.resize(1);
        result [ 0 ] = virtual_damage;
        return true;
    } else if ( code.compare("ft") == 0 ) {
        LDPMMaterial *m = static_cast< LDPMMaterial * >( mat );
        result.resize(1);
        result [ 0 ] = m->giveFt();
        return true;
    } else if ( code.compare("Gt") == 0 ) {
        LDPMMaterial *m = static_cast< LDPMMaterial * >( mat );
        result.resize(1);
        result [ 0 ] = m->giveGt();
        return true;
    } else if ( code.compare("fs") == 0 ) {
        LDPMMaterial *m = static_cast< LDPMMaterial * >( mat );
        result.resize(1);
        result [ 0 ] = m->giveFs();
        return true;
    } else if ( ( code.compare("strainN") == 0 ) ) {
        result.resize(1);
        result [ 0 ] = temp_strain [ 0 ];
        return true;
    } else if ( ( code.compare("stressN") == 0 ) ) {
        result.resize(1);
        result [ 0 ] = temp_stress [ 0 ];
        return true;
    } else {
        return VectMechMaterialStatus :: giveValues(code, result);
    }
}

//////////////////////////////////////////////////////////
double LDPMMaterialStatus :: giveStrengthLimit(double omega) {
    LDPMMaterial *m = static_cast< LDPMMaterial * >( mat );
    double sigma0;
    double rst = m->giveFs() / m->giveFt();
    double s = sin(omega);
    double s2 = s * s;
    double ac_rst = m->giveAlpha() * cos(omega) * cos(omega) / ( rst * rst );

    if ( omega > ( M_PI / 2 - 1e-6 ) ) {
        sigma0 = m->giveFt();
    } else {
        sigma0 = m->giveFt() * ( -s + sqrt(s2 + 4 * ac_rst) ) / ( 2 * ac_rst );
    }

    return sigma0;
}

//////////////////////////////////////////////////////////
Vector LDPMMaterialStatus :: giveTension(const Vector &strain, Vector strain_prev, Vector stress_prev) {
    LDPMMaterial *m = static_cast< LDPMMaterial * >( mat );

    // new strains & strains + stresses from previous step
    double epsN = strain [ 0 ], epsN_prev = strain_prev [ 0 ], strN_prev = stress_prev [ 0 ];
    double epsT, epsT_prev, strT_prev;

    if ( strain.size() == 2 ) {     //2D
        epsT = abs(strain [ 1 ]);
        epsT_prev = abs(strain_prev [ 1 ]);
        strT_prev = stress_prev [ 1 ];
    } else {    //3D
        epsT = sqrt(pow(strain [ 1 ], 2) + pow(strain [ 2 ], 2) );
        epsT_prev = sqrt(pow(strain_prev [ 1 ], 2) + pow(strain_prev [ 2 ], 2) );
        strT_prev = sqrt(pow(stress_prev [ 1 ], 2) + pow(stress_prev [ 2 ], 2) );
    }
    double epsEff = sqrt(pow(epsN, 2) + m->giveAlpha() * pow(epsT, 2) );        // effective strains
    double epsEff_prev;

    if ( epsN_prev < 1e-18 ) {
        double epsEff_prev_tmp = sqrt( pow(epsN_prev, 2) + m->giveAlpha() * pow(epsT_prev, 2) );
        epsEff_prev = sqrt( pow(epsEff_prev_tmp, 2) - pow(epsN_prev, 2) );
    } else {
        epsEff_prev = sqrt(pow(epsN_prev, 2) + m->giveAlpha() * pow(epsT_prev, 2) );
    }

    double strEff_prev = sqrt(pow(strN_prev, 2) + pow(strT_prev, 2) / m->giveAlpha() );

    // new max strains
    temp_maxEpsN = max(maxEpsN, epsN);
    temp_maxEpsT = max(maxEpsT, epsT);
    double temp_maxEpsEff = sqrt(pow(temp_maxEpsN, 2) + m->giveAlpha() * pow(temp_maxEpsT, 2) );        // max effective strains

    // elasticity
    double dEps = epsEff - epsEff_prev;
    double dStrElastic = m->giveE0() * dEps;

    // softening
    double omega, Lt, Ht, H0, str0, eps0, strBt;
    if ( epsT == 0 ) {
        omega = 0.5 * M_PI;
    } else {
        omega = atan(epsN / ( sqrt(m->giveAlpha() ) * epsT ) );
    }
    Lt = 2 * m->giveE0() * m->giveGt() / pow(m->giveFt(), 2);
    Ht = 2 * m->giveE0() / ( ( Lt / L ) - 1 );
    H0 = Ht * pow( 2 * omega / M_PI, m->givent() );
    str0 = giveStrengthLimit(omega);
    eps0 = str0 / m->giveE0();
    strBt = str0 * exp(-H0 * ( temp_maxEpsEff - eps0 ) / str0);

    // reloading
    double epsTr;
    epsTr = m->givekt() * ( temp_maxEpsEff - strBt / m->giveE0() );
    if ( epsEff_prev < epsTr ) {
        dStrElastic = m->giveE0() * max(epsEff - epsTr, 0.0);
    }

    // effective stress
    double strElastic = strEff_prev + dStrElastic;
    double strEff = max( 0.0, min(strElastic, strBt) );
    Vector intStress = Vector :: Zero(strain.size() );           // vector to collect stress

    if ( epsEff > 10e-20 ) {
        intStress [ 0 ] = strEff * strain [ 0 ] / epsEff;
        intStress [ 1 ] = m->giveAlpha() * strEff * strain [ 1 ] / epsEff;
        if ( strain.size() == 3 ) {
            intStress [ 2 ] = m->giveAlpha() * strEff * strain [ 2 ] / epsEff;
        }
    }
    return intStress;
}

//////////////////////////////////////////////////////////
double LDPMMaterialStatus :: giveSigmaBCDiff(double relt, double *sigmaBC) {
    double eN = updt_strain [ 0 ] + relt * ( temp_strain [ 0 ] - updt_strain [ 0 ] );
    double eV = volumetricStrain + relt * ( temp_volumetricStrain - volumetricStrain );

    LDPMMaterial *m = static_cast< LDPMMaterial * >( mat );
    double epsC0 = m->giveFc0() / m->giveE0();
    double epsC1 = m->giveKc0() * epsC0;
    double deviatoricStrain = eV - eN;
    double epsDV = eV + m->giveBeta() * deviatoricStrain;
    double epsV0 = m->giveKc3() * epsC0;
    double rDV = -abs(deviatoricStrain) / ( eV - epsV0 );
    double Hc = ( m->giveHc0() - m->giveHc1() ) / ( 1 + m->giveKc2() * max(0., rDV - m->giveKc1() ) ) + m->giveHc1();
    double sigmaC1 = m->giveFc0() + ( epsC1 - epsC0 ) * Hc;
    ( * sigmaBC ) = sigmaC1 * exp( ( -epsDV - epsC1 ) * Hc / sigmaC1);

    double deDdt = deVdt - deNdt;
    //absolut value for rDV
    double drDVdt = -( deDdt * ( eV - epsV0 ) - abs(deviatoricStrain) * deVdt ) / pow(eV - epsV0, 2);
    if ( deviatoricStrain < 0 ) {
        drDVdt = -( -deDdt * ( eV - epsV0 ) - abs(deviatoricStrain) * deVdt ) / pow(eV - epsV0, 2);
    }
    //max in Hc
    double dHcdt = 0;
    if ( rDV - m->giveKc1() > 0 ) {
        dHcdt = ( m->giveHc0() - m->giveHc1() ) / pow(1 + m->giveKc2() * max(0., rDV - m->giveKc1() ), 2) * m->giveKc2() * drDVdt;
    }
    double dsc1dt = ( epsC1 - epsC0 ) * dHcdt;
    double deDVdt  = deVdt + m->giveBeta() * deDdt;
    double dsBCdt = dsc1dt * exp( ( -epsDV - epsC1 ) * Hc / sigmaC1) + sigmaC1 * exp( ( -epsDV - epsC1 ) * Hc / sigmaC1) * ( -deDVdt * Hc * sigmaC1 + ( -epsDV - epsC1 ) * dHcdt * sigmaC1 - ( -epsDV - epsC1 ) * Hc * dsc1dt ) / pow(sigmaC1, 2);
    //cout << std::setprecision(15) << endl;

    //cout << temp_strain[0] << "\t" << rDV << "\t" << deviatoricStrain << "\t" << Hc << "\t" << sigmaBC << endl;
    //cout << 1 << "\t" << drDVdt/deNdt << "\t" << deDdt/deNdt << "\t" << dHcdt/deNdt <<  "\t" << dsBCdt/deNdt << endl;
    
    return dsBCdt / deNdt;
}

//////////////////////////////////////////////////////////
Vector LDPMMaterialStatus :: giveCompression(const Vector &strain, Vector strain_prev, Vector stress_prev) {
    LDPMMaterial *m = static_cast< LDPMMaterial * >( mat );

    temp_maxEpsN = max(maxEpsN, strain [ 0 ]);

    // NORMAL STRESS
    double deviatoricStrain, epsT, rDV, epsDV, epsV0, strNElastic, dEpsN;
    double sigmaBC, epsC0, Hc, epsC1;

    epsC0 = m->giveFc0() / m->giveE0();
    epsC1 = m->giveKc0() * epsC0;
    deviatoricStrain = temp_volumetricStrain - strain [ 0 ];
    epsDV = temp_volumetricStrain + m->giveBeta() * deviatoricStrain;


    dEpsN = strain [ 0 ] - strain_prev [ 0 ];
    // densified normal modulus and elastic step w/ or w/out change in normal modulus
    double ENc, dEps1, dEps2;
    if ( stress_prev [ 0 ] <=  -m->giveFc0() ) {
        ENc = m->giveEd();
    } else {
        ENc = m->giveE0();
    }

    strNElastic = stress_prev [ 0 ] + dEpsN * ENc;

    if ( ( m->giveFc0() + stress_prev [ 0 ] ) * ( m->giveFc0() + strNElastic ) < 0 ) {
        // Fc0 (+), str (-) -> if E0 changes to Ed, result negative
        dEps1 = ( -m->giveFc0() - stress_prev [ 0 ] ) / ENc;
        dEps2 = dEpsN - dEps1;
        if ( ENc == m->giveEd() ) {
            strNElastic = stress_prev [ 0 ] + dEps1 * ENc + dEps2 * m->giveE0();
        } else {
            strNElastic = stress_prev [ 0 ] + dEps1 * m->giveE0() + dEps2 * ENc;
        }
    } else {
        strNElastic = stress_prev [ 0 ] + ENc * dEpsN;
    }

    if ( strain.size() == 2 ) {     //2D
        epsT = abs(strain [ 1 ]);
    } else {    //3D
        epsT = sqrt(pow(strain [ 1 ], 2) + pow(strain [ 2 ], 2) );
    }

    /*
     * if ( epsT == 0 ) {
     *  omega = 0.5 * M_PI;
     * } else {
     *  omega = atan( strain [ 0 ] / ( sqrt( m->giveAlpha() ) * epsT ) );
     * }
     *
     * str0 = giveStrengthLimit(omega);
     * epsV0 = m->giveKc3() * str0 / m->giveE0();
     */
    epsV0 = m->giveKc3() * epsC0;

    if ( temp_volumetricStrain <= 0. ) {
        rDV = -abs(deviatoricStrain) / ( temp_volumetricStrain - epsV0 );
    } else {
        rDV = abs(deviatoricStrain) / epsV0;
    }

    Hc = ( m->giveHc0() - m->giveHc1() ) / ( 1 + m->giveKc2() * max(0., rDV - m->giveKc1() ) ) + m->giveHc1();

    // inelastic boundary
    if ( epsDV >= 0 ) {
        sigmaBC = m->giveFc0();
    } else if ( 0 <= -epsDV && -epsDV <= epsC1 ) {
        sigmaBC = m->giveFc0() + max(0., -epsDV - epsC0) * Hc;
    } else {
        double sigmaC1 = m->giveFc0() + ( epsC1 - epsC0 ) * Hc;
        sigmaBC = sigmaC1 * exp( ( -epsDV - epsC1 ) * Hc / sigmaC1);
        //check evolution of sigma bc derivative in time
        if ( deNdt > 0 ) {
            double transSigmaBC;
            double d0 = giveSigmaBCDiff(0, & transSigmaBC);
            double d1 = giveSigmaBCDiff(1, & transSigmaBC);
            if ( d1 < -m->giveEd() && d0 > -m->giveEd() ) {          
                //Newton method
                double transt = ( -m->giveEd() - d0 ) / ( d1 - d0 );
                double err = giveSigmaBCDiff(transt, & transSigmaBC) + m->giveEd();
                unsigned itmax = 20;
                unsigned it = 0;
                double diff;
                double Ed = m->giveEd();
                double shiftt = 1.;
                while ( abs(shiftt) > 1e-6 && it < itmax ) {
                    diff = ( ( giveSigmaBCDiff(transt + 1e-8, & transSigmaBC) + Ed ) - ( giveSigmaBCDiff(transt, & transSigmaBC) + Ed ) ) / 1e-8;
                    shiftt = err / diff;
                    if (shiftt>0) transt -= min(shiftt, transt);
                    else transt -= max(shiftt, transt-1);
                    err = giveSigmaBCDiff(transt, & transSigmaBC) + Ed;
                    it++;                
                }
                if ( it == itmax ) {
                    //bisection method
                    it = 0;
                    itmax = 60;
                    double t0 = 0.;
                    double t1 = 1.;
                    while ( t1-t0 > 1e-6 && it < itmax ) {
                        transt = (t0+t1)/2.;                    
                        err = giveSigmaBCDiff(transt, & transSigmaBC);
                        if (err > -m->giveEd()) t0=transt;
                        else t1=transt;                        
                        it++;                
                    }
                }
                if ( it == itmax ) {
                    cerr << "LDPM Material Error: transitional time not found, error " << fabs(err) << endl;
                    exit(1);
                }
                


                //found transitional time, from that time the evolution of bc should go with slope Ed
                sigmaBC = transSigmaBC + m->giveEd() * ( 1 - transt ) * deNdt;
            }
        }
    }

    Vector intStress = Vector :: Zero(strain.size() );           // vector to collect stress
    intStress [ 0 ] = min(0., max(-sigmaBC, strNElastic) );

    // SHEAR
    double dEpsM, dEpsL, strMElastic, strLElastic, strTElastic, strBs, strT;
    dEpsM = strain [ 1 ] - strain_prev [ 1 ];
    strMElastic = stress_prev [ 1 ] + m->giveEt() * dEpsM;
    if ( strain.size() == 2 ) {     //2D
        epsT = abs(strain [ 1 ]);
        strTElastic = abs(strMElastic);
    } else {    //3D
        epsT = sqrt(pow(strain [ 1 ], 2) + pow(strain [ 2 ], 2) );
        dEpsL = strain [ 2 ] - strain_prev [ 2 ];
        strLElastic = stress_prev [ 2 ] + m->giveEt() * dEpsL;
        strTElastic = sqrt(pow(strMElastic, 2) + pow(strLElastic, 2) );
    }

    temp_maxEpsT = max(maxEpsT, epsT);
    double dMu = m->giveMu0() - m->giveMuinf();
    strBs = m->giveFs() + dMu * m->giveFs0() - m->giveMuinf() * intStress [ 0 ] - dMu * m->giveFs0() * exp( intStress [ 0 ] / m->giveFs0() );
    strT = min( strBs, max(0.0, strTElastic) );

    if ( strT == 0 ) {
        intStress [ 1 ] = 0;
        if ( strain.size() == 3 ) {
            intStress [ 2 ] = 0;
        }
    } else {
        intStress [ 1 ] = strT * strMElastic / strTElastic;
        if ( strain.size() == 3 ) {
            intStress [ 2 ] = strT * strLElastic / strTElastic;
        }
    }

    return intStress;
}

//////////////////////////////////////////////////////////
Vector LDPMMaterialStatus :: passZero(const Vector &strain) {
    Vector intStrain = Vector :: Zero(strain.size() );     // itermediary strains when passing 0
    Vector intStress = Vector :: Zero(strain.size() );     // itermediary stresses when passing 0

    intStrain [ 0 ] = 0;
    intStrain [ 1 ] = updt_strain [ 1 ] + ( strain [ 1 ] - updt_strain [ 1 ] ) * ( -updt_strain [ 0 ] ) / ( strain [ 0 ] - updt_strain [ 0 ] );
    if ( strain.size() == 3 ) {
        intStrain [ 2 ] = updt_strain [ 2 ] + ( strain [ 2 ] - updt_strain [ 2 ] ) * ( -updt_strain [ 0 ] ) / ( strain [ 0 ] - updt_strain [ 0 ] );
    }

    if ( updt_strain [ 0 ] < 0 ) {  // passing from compression to tension
        intStress = giveCompression(intStrain, updt_strain, updt_stress);
        temp_stress = giveTension(strain, intStrain, intStress);
    } else if ( updt_strain [ 0 ] > 0 ) {  // passing from tension to compression
        intStress = giveTension(intStrain, updt_strain, updt_stress);
        temp_stress = giveCompression(strain, intStrain, intStress);
    }

    return temp_stress;
}

//////////////////////////////////////////////////////////
Vector LDPMMaterialStatus :: giveStress(const Vector &strain, double timeStep) {
    temp_strain = strain;
    temp_mech_strain = addEigenStrain(strain);

    //compute differentiations
    deVdt = ( temp_volumetricStrain - volumetricStrain ) / timeStep;
    deNdt = ( temp_strain [ 0 ] - updt_strain [ 0 ] ) / timeStep;

    LDPMMaterial *m = static_cast< LDPMMaterial * >( mat );
    double epsNState = temp_mech_strain [ 0 ] * updt_mech_strain [ 0 ];  // gives information about the evolution of normal strains
    if ( epsNState < 0 ) {  // change of sign of EpsN
        temp_stress = passZero(temp_mech_strain);
    } else {
        //if ( min(temp_mech_strain [ 0 ], updt_mech_strain[0])<= 0 ) {  // normal evolution in compression
        if ( temp_mech_strain [ 0 ] > 0 || updt_mech_strain [ 0 ] > 0 ) {     // normal evolution in tension
            temp_stress = giveTension(temp_mech_strain, updt_mech_strain, updt_stress);
        } else {
            temp_stress = giveCompression(temp_mech_strain, updt_mech_strain, updt_stress);
        }
    }

    if ( temp_stress [ 0 ] > 0 ) {
        temp_crackOpening = ( temp_mech_strain [ 0 ] - ( temp_stress [ 0 ] / m->giveE0() ) ) * L;
    } else {
        temp_crackOpening = 0;
    }

    giveVirtualDamage();

    return temp_stress;
}

//////////////////////////////////////////////////////////
void LDPMMaterialStatus :: giveVirtualDamage() {
    LDPMMaterial *m = static_cast< LDPMMaterial * >( mat );
    double temp_epsEff, temp_strEff;

    double epsN = temp_mech_strain [ 0 ], strN = temp_stress [ 0 ];
    double epsT, strT;

    if ( temp_mech_strain.size() == 2 ) {     //2D
        epsT = abs(temp_mech_strain [ 1 ]);
        strT = abs(temp_stress [ 1 ]);
    } else {    //3D
        epsT = sqrt(pow(temp_mech_strain [ 1 ], 2) + pow(temp_mech_strain [ 2 ], 2) );
        strT = sqrt(pow(temp_stress [ 1 ], 2) + pow(temp_stress [ 2 ], 2) );
    }

    temp_epsEff = sqrt(pow(epsN, 2) + m->giveAlpha() * pow(epsT, 2) );        // effective strains
    temp_strEff = sqrt(pow(strN, 2) + pow(strT, 2) / m->giveAlpha() );         // effective stress

    double temp_E;
    if ( epsN < -m->giveFc0() ) {
        temp_E = m->giveEd();
    } else {
        temp_E = m->giveE0();
    }

    virtual_damage = 1 - temp_strEff / ( temp_E * temp_epsEff );

    if ( virtual_damage < 1e-10 ) {
        virtual_damage = 0.;
    }
}

//////////////////////////////////////////////////////////
void LDPMMaterialStatus :: update() {
    VectMechMaterialStatus :: update();
    maxEpsN = temp_maxEpsN;
    maxEpsT = temp_maxEpsT;

    crackOpening = temp_crackOpening;
    updt_mech_strain = temp_mech_strain;
    volumetricStrain = temp_volumetricStrain;
}

//////////////////////////////////////////////////////////
void LDPMMaterialStatus :: resetTemporaryVariables() {
    VectMechMaterialStatus :: resetTemporaryVariables();
    temp_maxEpsN = maxEpsN;
    temp_maxEpsT = maxEpsT;
    temp_crackOpening = crackOpening;
}

//////////////////////////////////////////////////////////
Matrix LDPMMaterialStatus :: giveStiffnessTensor(string type) const {
    Matrix stiff = VectMechMaterialStatus :: giveStiffnessTensor(type);
    /*if ( type.compare("elastic") == 0 ) {
     *  return stiff;
     * } else if ( type.compare("secant") == 0 ) {
     *  return stiff;
     * } else if ( type.compare("unloading") == 0 ) {
     *  return stiff;
     * } else if ( type.compare("tangent") == 0 ) {
     *  return stiff;       //not implemented, used unloading
     * } else {
     *  cerr << "Error: LDPMMaterialStatus does not provide '" << type << "' stiffness";
     *  exit(1);
     * };
     */

    if ( type.compare("elastic") == 0 ) {
        return stiff;
    } else if ( type.compare("secant") == 0 ) {
        //LDPMMaterial *m = static_cast< LDPMMaterial * >( mat );
        //stiff(0,0) *= max( 1 - virtual_damage, m->giveDamageResiduum() );
        return stiff;
    } else {
        cerr << "Error: LDPMMaterialStatus does not provide '" << type << "' stiffness";
        exit(1);
    };
}

//////////////////////////////////////////////////////////
std :: string LDPMMaterialStatus :: giveLineToSave() const {
    return "maxEpsN " + to_string_sci(this->maxEpsN) + " maxEpsT " + to_string_sci(this->maxEpsT);
}

//////////////////////////////////////////////////////////
void LDPMMaterialStatus :: setParameterValue(string code, double value) {
    if ( code.compare("volumetric_strain") == 0 ) {
        temp_volumetricStrain = value;     // volumetric change, cela stopa matice, proto *3
    } else {
        VectMechMaterialStatus :: setParameterValue(code, value);
    }
}

//////////////////////////////////////////////////////////
void LDPMMaterialStatus :: readFromLine(istringstream &iss) {
    std :: string param;
    while (  iss >> param ) {
        if ( param.compare("maxEpsN") == 0 ) {
            iss >> this->maxEpsN;
            temp_maxEpsN = maxEpsN;
        } else if ( param.compare("maxEpsT") == 0 ) {
            iss >> this->maxEpsT;
            temp_maxEpsT = maxEpsT;
        }
    }
}

//////////////////////////////////////////////////////////
Vector LDPMMaterialStatus :: giveStressWithFrozenIntVars(const Vector &strain, double timeStep) {
    ( void ) timeStep;
    LDPMMaterial *m = static_cast< LDPMMaterial * >( mat );

    temp_strain = strain;
    temp_mech_strain = addEigenStrain(strain);

    // absolute prediction
    // double strN_tmp = m->giveE0() * temp_mech_strain [ 0 ];
    // if ( temp_stress [ 0 ] < m->giveFc0() ) {
    //     double eps1 = temp_mech_strain [ 0 ] * m->giveFc0() / strN_tmp;
    //     double eps2 = temp_mech_strain [ 0 ] * ( strN_tmp - m->giveFc0() ) / strN_tmp;
    //     temp_stress [ 0 ] = -eps1 *m->giveE0() - eps2 * m->giveEd();
    // } else {
    //     temp_stress [ 0 ] = strN_tmp;
    // }
    // temp_stress [ 1 ] = m->giveEt() * temp_mech_strain [ 1 ];
    // if ( strain.size() == 3 ) {
    //     temp_stress [ 2 ] = m->giveEt() * temp_mech_strain [ 2 ];
    // }

    // incremental prediction
    double strN_tmp = updt_stress [ 0 ] + m->giveE0() * ( temp_mech_strain [ 0 ] - updt_mech_strain [ 0 ] );

    // if ( strN_tmp < - m->giveFc0() ) { // Fc0 is positive
    //     double eps1 = abs( temp_mech_strain [ 0 ] * m->giveFc0() / strN_tmp );
    //     double eps2 = abs( temp_mech_strain [ 0 ] * ( strN_tmp + m->giveFc0() ) / strN_tmp );
    //     temp_stress [ 0 ] = - eps1 * m->giveE0() - eps2 * m->giveEd();
    // } else {
    //     temp_stress [ 0 ] = strN_tmp;
    // }

    temp_stress [ 0 ] = strN_tmp;

    temp_stress [ 1 ] = updt_stress [ 1 ] + m->giveEt() * ( temp_mech_strain [ 1 ] - updt_mech_strain [ 1 ] );
    if ( strain.size() == 3 ) {
        temp_stress [ 2 ] = updt_stress [ 2 ] + m->giveEt() * ( temp_mech_strain [ 2 ] - updt_mech_strain [ 2 ] );
    }

    return temp_stress;
}

//////////////////////////////////////////////////////////
bool LDPMMaterialStatus :: isElastic(const bool &now) const {
    ( void ) now;
    // if ( now && this->virtual_damage != 0.0 ) {
    //     return false;
    // } else if ( this->damage != 0.0 ) {
    //     return false;
    // }
    return true;
}

//////////////////////////////////////////////////////////
// LDPM MATERIAL

//////////////////////////////////////////////////////////
void LDPMMaterial :: readFromLine(istringstream &iss) {
    VectMechMaterial :: readFromLine(iss); //read elastic parameters

    iss.clear(); // clear string stream
    iss.seekg(0, iss.beg); //reset position in string stream

    // initialize all values to zero (NOTE probably no need in linux, but in windows necessary)
    nt = kt = beta = fc = Ed = Hc0 = Hc1 = Kc0 = Kc1 = Kc2 = Kc3 = fs = mu0 = muinf = Et = -1;

    string param;
    bool bft, bGt, bfc0, bfs0;
    bft = bGt = bfc0 = bfs0 = false;

    while (  iss >> param ) {
        if ( param.compare("Gt") == 0 ) {
            bGt = true;
            iss >> Gt;
        } else if ( param.compare("ft") == 0 ) {
            bft = true;
            iss >> ft;
        } else if ( param.compare("nt") == 0 ) {
            iss >> nt;
        } else if ( param.compare("kt") == 0 ) {
            iss >> kt;
        } else if ( param.compare("beta") == 0 ) {
            iss >> beta;
        } else if ( param.compare("fc") == 0 ) { //NO USE???
            iss >> fc;
        } else if ( param.compare("fc0") == 0 ) {
            bfc0 = true;
            iss >> fc0;
        } else if ( param.compare("Ed") == 0 ) {
            iss >> Ed;
        } else if ( param.compare("Hc0") == 0 ) {
            iss >> Hc0;
        } else if ( param.compare("Hc1") == 0 ) {
            iss >> Hc1;
        } else if ( param.compare("Kc0") == 0 ) {
            iss >> Kc0;
        } else if ( param.compare("Kc1") == 0 ) {
            iss >> Kc1;
        } else if ( param.compare("Kc2") == 0 ) {
            iss >> Kc2;
        } else if ( param.compare("Kc3") == 0 ) {
            iss >> Kc3;
        } else if ( param.compare("fs") == 0 ) {
            iss >> fs;
        } else if ( param.compare("fs0") == 0 ) { //sigma N0, transitional stress
            bfs0 = true;
            iss >> fs0;
        } else if ( param.compare("Et") == 0 ) {
            iss >> Et;
        } else if ( param.compare("mu0") == 0 ) {
            iss >> mu0;
        } else if ( param.compare("muinf") == 0 ) {
            iss >> muinf;
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
    if ( !bfc0 ) {
        cerr << name << ": material parameter 'fc0' was not specified" << endl;
        exit(EXIT_FAILURE);
    }
    ;
    if ( !bfs0 ) {
        cerr << name << ": material parameter 'fs0' was not specified" << endl;
        exit(EXIT_FAILURE);
    }
    ;
};

//////////////////////////////////////////////////////////
MaterialStatus *LDPMMaterial :: giveNewMaterialStatus(Element *e, unsigned ipnum) {
    LDPMMaterialStatus *newStatus = new LDPMMaterialStatus(this, e, ipnum); //needs to be deleted manually
    return newStatus;
};

//////////////////////////////////////////////////////////
void LDPMMaterial :: init(MaterialContainer *matcont) {
    VectMechMaterial :: init(matcont);
    ;

    // if variables not specified on the input, use default multipliers
    nt = ( nt == -1 ) ? 0.2 : nt;
    kt = ( kt == -1 ) ? 0.5 : kt;
    beta = ( beta == -1 ) ? 0. : beta;

    fc = ( fc == -1 ) ? 16 * ft : fc;
    Ed = ( Ed == -1 ) ? 2 * E0 : Ed;
    Hc0 = ( Hc0 == -1 ) ? 0.6 * E0 : Hc0;
    Hc1 = ( Hc1 == -1 ) ? 0.1 * E0 : Hc1;
    Kc0 = ( Kc0 == -1 ) ? 4 : Kc0;
    Kc1 = ( Kc1 == -1 ) ? 1 : Kc1;
    Kc2 = ( Kc2 == -1 ) ? 10 : Kc2;
    Kc3 = ( Kc3 == -1 ) ? 0.1 : Kc3;

    fs = ( fs == -1 ) ? 3 * ft : fs;
    Et = ( Et == -1 ) ? alpha * E0 : Et;
    mu0 = ( mu0 == -1 ) ? 0.1 : mu0;
    muinf = ( muinf == -1 ) ? 0.0125 : muinf;

    damage_residuum = ( damage_residuum == -1 ) ? 0. : damage_residuum;
};
