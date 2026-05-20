#include "material_csl.h"
#include "element_discrete.h"
#include "element_ldpm.h"
#include "element_container.h"
#include <algorithm>
#include <cmath>
#include <cstring>


using namespace std;

namespace {
constexpr double CSL_TINY = 1e-14;

void hashDoubleCSL(std :: uint64_t &hash, double value) {
    std :: uint64_t bits = 0;
    std :: memcpy(&bits, &value, sizeof(double) );
    hash ^= bits;
    hash *= 1099511628211ULL;
}

double clampOmegaForCSL(double omega) {
    return std :: max(-0.5 * M_PI + 1e-10, std :: min(0.5 * M_PI - 1e-10, omega) );
}
}

//////////////////////////////////////////////////////////
// CSL MATERIAL STATUS

CSLMaterialStatus :: CSLMaterialStatus(CSLMaterial *m, Element *e, unsigned ipnum) : VectMechMaterialStatus(m, e, ipnum) {
    name = "CSL mat. status";
    RAND_H = 1.0;
}

//////////////////////////////////////////////////////////
bool CSLMaterialStatus :: giveValues(string code, Vector &result) const {
    if ( code.compare("tempCrackOpening") == 0 || code.compare("crack_opening") == 0 || code.compare("normal_crack_opening") == 0 ) {
        result.resize(1);
        result [ 0 ] = temp_crackOpening;
        return true;
    } else if ( code.rfind("damage", 0) == 0 || code.rfind("damageN", 0) == 0 || code.rfind("damageT", 0) == 0 ) {
        result.resize(1);
        result [ 0 ] = temp_damage;
        return true;
    } else if ( code.compare("ft") == 0 ) {
        CSLMaterial *m = static_cast< CSLMaterial * >( mat );
        result.resize(1);
        result [ 0 ] = m->giveFt();
        return true;
    } else if ( code.compare("Gt") == 0 ) {
        CSLMaterial *m = static_cast< CSLMaterial * >( mat );
        result.resize(1);
        result [ 0 ] = m->giveGt();
        return true;
    } else if ( code.compare("fs") == 0 ) {
        CSLMaterial *m = static_cast< CSLMaterial * >( mat );
        result.resize(1);
        result [ 0 ] = m->giveFs();
        return true;
    } else if ( code.compare("Gs") == 0 ) {
        CSLMaterial *m = static_cast< CSLMaterial * >( mat );
        result.resize(1);
        result [ 0 ] = m->giveGs();
        return true;
    } else if ( ( code.compare("strainN") == 0 ) ) {
        result.resize(1);
        result [ 0 ] = temp_strain [ 0 ];
        return true;
    } else if ( ( code.compare("stressN") == 0 ) ) {
        result.resize(1);
        result [ 0 ] = temp_stress [ 0 ];
        return true;
    } else if ( code.compare("normal_dissipation_density") == 0 ) {
        result.resize(1);
        result [ 0 ] = normalEnergyDensity - 0.5 * temp_stress [ 0 ] * temp_strain [ 0 ];
        return true;
    } else if ( code.compare("shear_dissipation_density") == 0 ) {
        VectMechMaterial *m = static_cast< VectMechMaterial * >( mat );
        result.resize(1);
        if ( m->giveDimension() == 2 ) {
            result [ 0 ] = shearEnergyDensity  - 0.5 * temp_stress [ 1 ] * temp_strain [ 1 ];
        } else  {
            result [ 0 ] = shearEnergyDensity - 0.5 * ( temp_stress [ 1 ] * temp_strain [ 1 ] + temp_stress [ 2 ] * temp_strain [ 2 ] );
        }
        return true;
    } else {
        return VectMechMaterialStatus :: giveValues(code, result);
    }
}

//////////////////////////////////////////////////////////
void CSLMaterialStatus :: init() {
    computeOmega0();
    maxEpsT = 0;
    maxEpsN = 0;
    damage = 0;
    temp_damage = 0;
    temp_maxEpsN = 0;
    temp_maxEpsT = 0;
    temp_crackOpening = 0;


    crackOpening = 0;

    RigidBodyContact *rbc = dynamic_cast< RigidBodyContact * >( element );
    LDPMTetra *tet = dynamic_cast< LDPMTetra * >( element );
    MaterialTestElement *mte = dynamic_cast< MaterialTestElement * >( element );
    if ( rbc ) {
        L = rbc->giveLength();
    } else if ( tet ) {
        L = tet->giveLength(idx);
    } else if ( mte ) {
        L = 1.;
    } else {
        cerr << "Material " << name << " can be used only for RigidBodyContact, LDPMTetra, or MaterialTestElement elements" << endl;
        exit(EXIT_FAILURE);
    }

    CSLMaterial *m = static_cast< CSLMaterial * >( mat );
    Ks = 2 * m->giveAlphaForDamage() * m->giveE0() / ( m->giveLcrs() / L - 1 );
    Kt = 2 * m->giveE0() / ( m->giveLcrt() / L - 1 );
    nt = log( Kt / ( Kt - Ks ) ) / log(1 - 2 * omega0 / M_PI);

    if ( Ks < 0 || Kt < 0 ) {
        cerr << "Error in " << name << ": snap back occured" << endl;
        //exit(1);
    }

    if ( Kt <= Ks ) {
        cerr << "Error in " << name << ": Ks (" << Ks << ") is greater than Kt (" << Kt << ")" << L << endl;
        //exit(1);
    }
}

//////////////////////////////////////////////////////////
double CSLMaterialStatus :: giveS0tension(double omega) const {
    CSLMaterial *m = static_cast< CSLMaterial * >( mat );

    double s = sin(omega);
    double s2 = s * s;
    double c2 = cos(omega) * cos(omega);
    double ft = m->giveFt() * RAND_H;
    double fs = m->giveFs() * RAND_H;

    double sa = .5 * ft * ( pow(fs / ( m->giveMu() * ft ), 2) - 1. );


    if ( omega == atan( sqrt( m->giveAlphaForDamage() ) / m->giveMu() ) ) {
        // for this anle, the later equation is undetermined, but hyperbola eq. gives this (see Two Scale Study - Cusatis 2007 doi.org/10.1016/j.engfracmech.2006.01.021)
        return .5 * ( ft + 2 * sa ) * ft / ( ( ft + sa ) * s );
    } else {
        return ( -( ft + sa ) * s + sqrt(pow( ( ft + sa ) * s, 2) + ( m->giveAlphaForDamage() * ( c2 / pow(m->giveMu(), 2) ) - s2 ) * ( ft + 2 * sa ) * ft) ) / ( m->giveAlphaForDamage() * ( c2 / pow(m->giveMu(), 2) ) - s2 );
    }
}

//////////////////////////////////////////////////////////
double CSLMaterialStatus :: giveS0compression(double omega) const {
    CSLMaterial *m = static_cast< CSLMaterial * >( mat );

    double fc = m->giveFc() * RAND_H;

    return fc / sqrt( pow(sin(omega), 2) + ( m->giveAlphaForDamage() * pow(cos(omega), 2) ) / m->giveBeta() );
}

//////////////////////////////////////////////////////////
void CSLMaterialStatus :: computeOmega0() {
    double step = 0.01;
    double minomega = 0;
    double minerr = 1e100;
    double err;
    for ( omega0 = -M_PI / 2.; omega0 < 0; omega0 += step ) {
        err = giveS0tension(omega0) - giveS0compression(omega0);
        if ( minerr > fabs(err) ) {
            minomega = omega0;
            minerr = fabs(err);
        }
    }
    omega0 = minomega;

    unsigned itmax = 1000;
    unsigned it = 0;
    double diff;
    while ( fabs(err) > 1e-5 && it < itmax ) {
        diff = ( ( giveS0tension(omega0 + 1e-8) - giveS0compression(omega0 + 1e-8) ) - ( giveS0tension(omega0) - giveS0compression(omega0) ) ) / 1e-8;
        omega0 -= err / diff;
        err = giveS0tension(omega0) - giveS0compression(omega0);
        it++;
    }
    if ( it == itmax ) {
        cerr << "CSL Material Error: omega0 does not converging, error " << fabs(err) << endl;
        exit(1);
    }
}

//////////////////////////////////////////////////////////
void CSLMaterialStatus :: computeDamage(Vector strain) {
    CSLMaterial *m = static_cast< CSLMaterial * >( mat );
    double epsN = strain [ 0 ];
    double epsT;
    if ( strain.size() == 2 ) {
        epsT = abs(strain [ 1 ]);                //2D
    } else {
        epsT = sqrt( pow(strain [ 1 ], 2) + pow(strain [ 2 ], 2) );          //3D
    }
    double epsEQ = sqrt( pow(epsN, 2) + m->giveAlphaForDamage() * pow(epsT, 2) );          //equivalent strain

    if ( epsEQ > 0 && damage < 1.0 ) {
        double omega, S0, chi, K0, strEQ;
        if ( epsT > 0 ) {
            omega = atan( epsN / ( sqrt( m->giveAlphaForDamage() ) * epsT ) );
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
            double emax;
            temp_maxEpsN = max(maxEpsN, epsN);
            temp_maxEpsT = max(maxEpsT, epsT);
            emax = sqrt( pow(temp_maxEpsN, 2) + m->giveAlphaForDamage() * pow(temp_maxEpsT, 2) );
            S0 = giveS0tension(omega);

            unsigned dim = element->giveDimension();
            double flam = 1.;
            temp_volumetricStrain = addEigenVolumetricStrain(temp_volumetricStrain_total);
            if ( m->giveLam0() > 0 ) {
                flam = 1. / ( 1. + max(-( temp_volumetricStrain * dim - epsN ) / ( dim * m->giveLam0() ), 0.) );                  //projection of the trace perpendicularly to the connection
            }
            K0 = -flam * Kt * ( 1. - pow( ( omega - 0.5 * M_PI ) / ( omega0 - 0.5 * M_PI ), nt) );

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
Vector CSLMaterialStatus :: computeDamageGradient(Vector strain) const {
    Vector gradient = Vector :: Zero(strain.size() );
    CSLMaterial *m = static_cast< CSLMaterial * >( mat );

    const double epsN = strain [ 0 ];
    double epsT2 = 0.;
    for ( unsigned i = 1; i < strain.size(); i++ ) {
        epsT2 += strain [ i ] * strain [ i ];
    }
    const double epsT = sqrt(epsT2);
    const double alpha = m->giveAlphaForDamage();
    const double sqrtAlpha = sqrt(alpha);
    const double epsEQ = sqrt(epsN * epsN + alpha * epsT2);

    if ( epsEQ <= CSL_TINY || damage >= 1.0 ) {
        return gradient;
    }

    Vector depsEQ = Vector :: Zero(strain.size() );
    depsEQ [ 0 ] = epsN / epsEQ;
    for ( unsigned i = 1; i < strain.size(); i++ ) {
        depsEQ [ i ] = alpha * strain [ i ] / epsEQ;
    }

    double omega;
    Vector domega = Vector :: Zero(strain.size() );
    if ( epsT > CSL_TINY ) {
        omega = atan(epsN / ( sqrtAlpha * epsT ) );
        const double epsEQ2 = epsEQ * epsEQ;
        domega [ 0 ] = sqrtAlpha * epsT / epsEQ2;
        for ( unsigned i = 1; i < strain.size(); i++ ) {
            domega [ i ] = -sqrtAlpha * epsN * strain [ i ] / ( epsEQ2 * epsT );
        }
    } else if ( epsN > 0 ) {
        omega = 0.5 * M_PI;
    } else {
        omega = -0.5 * M_PI;
    }

    auto omegaDerivative = [] (auto fn, double x) {
        const double h = 1e-7;
        const double xMinus = clampOmegaForCSL(x - h);
        const double xPlus = clampOmegaForCSL(x + h);
        if ( fabs(xPlus - xMinus) <= CSL_TINY ) {
            return 0.;
        }
        return ( fn(xPlus) - fn(xMinus) ) / ( xPlus - xMinus );
    };

    double S0 = 0.;
    double K0 = 0.;
    double chi = 0.;
    Vector dS0 = Vector :: Zero(strain.size() );
    Vector dK0 = Vector :: Zero(strain.size() );
    Vector dchi = Vector :: Zero(strain.size() );

    if ( omega < omega0 ) { // compression
        S0 = giveS0compression(omega);
        const double dS0dOmega = omegaDerivative([this] (double x) { return giveS0compression(x); }, omega);
        dS0 = dS0dOmega * domega;

        chi = epsEQ;
        dchi = depsEQ;

        const double denom = omega0 + 0.5 * M_PI;
        const double ratio = ( omega + 0.5 * M_PI ) / denom;
        const double ratioPow = pow(std :: max(ratio, 0.), m->giveNc() );
        K0 = m->giveKc() * ( 1. - ratioPow );
        double dK0dOmega = 0.;
        if ( ratio > CSL_TINY ) {
            dK0dOmega = -m->giveKc() * m->giveNc() * pow(ratio, m->giveNc() - 1.) / denom;
        }
        dK0 = dK0dOmega * domega;
    } else { // tension-shear
        S0 = giveS0tension(omega);
        const double dS0dOmega = omegaDerivative([this] (double x) { return giveS0tension(x); }, omega);
        dS0 = dS0dOmega * domega;

        const double trialMaxEpsN = std :: max(maxEpsN, epsN);
        const double trialMaxEpsT = std :: max(maxEpsT, epsT);
        const double emax = sqrt(trialMaxEpsN * trialMaxEpsN + alpha * trialMaxEpsT * trialMaxEpsT);
        Vector demax = Vector :: Zero(strain.size() );
        if ( emax > CSL_TINY ) {
            if ( epsN > maxEpsN ) {
                demax [ 0 ] = trialMaxEpsN / emax;
            }
            if ( epsT > maxEpsT && epsT > CSL_TINY ) {
                for ( unsigned i = 1; i < strain.size(); i++ ) {
                    demax [ i ] = alpha * trialMaxEpsT * strain [ i ] / ( emax * epsT );
                }
            }
        }

        double flam = 1.;
        Vector dflam = Vector :: Zero(strain.size() );
        if ( m->giveLam0() > 0 ) {
            const unsigned dim = element->giveDimension();
            const double tempVolumetricStrain = addEigenVolumetricStrain(temp_volumetricStrain_total);
            const double confinement = ( epsN - tempVolumetricStrain * dim ) / ( dim * m->giveLam0() );
            if ( confinement > 0. ) {
                flam = 1. / ( 1. + confinement );
                dflam [ 0 ] = -1. / ( dim * m->giveLam0() * pow(1. + confinement, 2) );
            }
        }

        const double denom = omega0 - 0.5 * M_PI;
        const double ratio = ( omega - 0.5 * M_PI ) / denom;
        const double safeRatio = std :: max(ratio, 0.);
        const double ratioPow = pow(safeRatio, nt);
        K0 = -flam * Kt * ( 1. - ratioPow );
        double dK0dOmega = 0.;
        if ( safeRatio > CSL_TINY ) {
            dK0dOmega = flam * Kt * nt * pow(safeRatio, nt - 1.) / denom;
        }
        dK0 = dK0dOmega * domega - Kt * ( 1. - ratioPow ) * dflam;

        if ( omega < 0. ) {
            chi = epsEQ * omega / omega0 + emax * ( 1. - omega / omega0 );
            dchi = ( omega / omega0 ) * depsEQ
                    + ( ( epsEQ - emax ) / omega0 ) * domega
                    + ( 1. - omega / omega0 ) * demax;
        } else {
            chi = emax;
            dchi = demax;
        }
    }

    if ( S0 <= CSL_TINY || epsEQ <= CSL_TINY ) {
        return gradient;
    }

    double strEQ = S0;
    Vector dStrEQ = dS0;
    if ( chi - S0 / m->giveE0() > 0 ) {
        const double exponent = K0 / S0 * ( chi - S0 / m->giveE0() );
        strEQ = S0 * exp(exponent);
        dStrEQ = strEQ * (
                      dS0 / S0
                      + ( chi / S0 - 1. / m->giveE0() ) * dK0
                      + ( K0 / S0 ) * dchi
                      - ( K0 * chi / ( S0 * S0 ) ) * dS0
                  );
    }

    const double rawDamage = 1. - strEQ / ( m->giveE0() * epsEQ );
    if ( rawDamage <= damage ) {
        return gradient;
    }

    gradient = -dStrEQ / ( m->giveE0() * epsEQ )
               + ( strEQ / ( m->giveE0() * epsEQ * epsEQ ) ) * depsEQ;
    return gradient;
}

//////////////////////////////////////////////////////////
double CSLMaterialStatus :: giveEnergyDissipationIncrement() const {
    return ( updt_stress.dot(temp_strain) - temp_stress.dot(updt_strain) ) / 2;
}

//////////////////////////////////////////////////////////
void CSLMaterialStatus :: computeEnergyDensities() {
    VectMechMaterialStatus :: computeEnergyDensities();
}

//////////////////////////////////////////////////////////
void CSLMaterialStatus :: update() {
    VectMechMaterialStatus :: update();
    damage = temp_damage;
    maxEpsN = temp_maxEpsN;
    maxEpsT = temp_maxEpsT;

    crackOpening = temp_crackOpening;
}

//////////////////////////////////////////////////////////
void CSLMaterialStatus :: resetTemporaryVariables() {
    VectMechMaterialStatus :: resetTemporaryVariables();
    temp_damage = damage;
    temp_maxEpsN = maxEpsN;
    temp_maxEpsT = maxEpsT;
    temp_crackOpening = crackOpening;
}



//////////////////////////////////////////////////////////
std :: unique_ptr< MaterialStatus > CSLMaterialStatus :: cloneState() const {
    std :: unique_ptr< CSLMaterialStatus > copy = std :: make_unique< CSLMaterialStatus >( *this );
    copy->clearSnapshotComponentPointers();
    return copy;
}

//////////////////////////////////////////////////////////
void CSLMaterialStatus :: restoreStateFrom(const MaterialStatus &other) {
    const CSLMaterialStatus *otherCSL = dynamic_cast< const CSLMaterialStatus * >( &other );
    if ( !otherCSL ) {
        std :: cerr << "Material status snapshot restore type mismatch for " << name << std :: endl;
        exit(1);
    }
    *this = *otherCSL;
    clearSnapshotComponentPointers();
}

//////////////////////////////////////////////////////////
std :: uint64_t CSLMaterialStatus :: stateHash() const {
    std :: uint64_t hash = MaterialStatus :: stateHash();
    hashDoubleCSL(hash, normalEnergyDensity);
    hashDoubleCSL(hash, shearEnergyDensity);
    hashDoubleCSL(hash, eigenVolumetricStrain);
    hashDoubleCSL(hash, temp_volumetricStrain);
    hashDoubleCSL(hash, volumetricStrain);
    hashDoubleCSL(hash, temp_volumetricStrain_total);
    hashDoubleCSL(hash, volumetricStrain_total);
    hashDoubleCSL(hash, omega0);
    hashDoubleCSL(hash, maxEpsT);
    hashDoubleCSL(hash, maxEpsN);
    hashDoubleCSL(hash, temp_maxEpsT);
    hashDoubleCSL(hash, temp_maxEpsN);
    hashDoubleCSL(hash, damage);
    hashDoubleCSL(hash, temp_damage);
    hashDoubleCSL(hash, Kt);
    hashDoubleCSL(hash, Ks);
    hashDoubleCSL(hash, L);
    hashDoubleCSL(hash, nt);
    hashDoubleCSL(hash, RAND_H);
    hashDoubleCSL(hash, crackOpening);
    hashDoubleCSL(hash, temp_crackOpening);
    return hash;
}

//////////////////////////////////////////////////////////
Matrix CSLMaterialStatus :: giveStiffnessTensor(string type) const {
    Matrix stiff = VectMechMaterialStatus :: giveStiffnessTensor(type);
    if ( type.compare("elastic") == 0 ) {
        return stiff;
    } else if ( type.compare("secant") == 0 ) {
        CSLMaterial *m = static_cast< CSLMaterial * >( mat );
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
                double strs = sqrt( pow(sN, 2) + ( sT / m->giveAlphaForDamage() ) );
                if ( strs  < m->giveStressResiduum() ) {
                    double epsN, epsT;
                    epsN = temp_strain [ 0 ];
                    epsT = 0.0;
                    for ( unsigned i = 1; i < temp_strain.size(); i++ ) {
                        epsT += pow(temp_strain [ i ], 2);
                    }
                    double epsEQ = sqrt( pow(epsN, 2) + epsT * m->giveAlphaForDamage() );
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
    } else if ( type.compare("archived_csl_damage_tangent") == 0 ) {
        // Archived experiment: locally consistent CSL active-damage tangent.
        // It closes finite-difference diagnostics but failed TS-N65 strict phase closure.
        Matrix tangent = stiff * ( 1 - temp_damage );
        if ( temp_damage > damage ) {
            const Vector elasticStress = stiff * temp_strain;
            const Vector damageGradient = computeDamageGradient(temp_strain);
            tangent -= elasticStress * damageGradient.transpose();
        }
        return tangent;
    } else if ( type.compare("consistent") == 0 ) {
        // Reference actual-state tangent for damage-growth diagnostics.
        const double eps = 1e-8;
        Matrix consistent = Matrix :: Zero(stiff.rows(), stiff.cols() );
        const Vector baseStress = temp_stress;
        const Vector baseStrain = temp_strain_total;
        for ( unsigned i = 0; i < static_cast< unsigned >( baseStrain.size() ); i++ ) {
            std :: unique_ptr< MaterialStatus > perturbed = cloneState();
            Vector perturbedStrain = baseStrain;
            perturbedStrain [ i ] += eps;
            perturbed->setTotalTempStrain(perturbedStrain);
            perturbed->computeStress(1.);
            consistent.col(i) = ( perturbed->giveTempStress() - baseStress ) / eps;
        }
        return consistent;
    } else {
        cerr << "Error: CSLMaterialStatus does not provide '" << type << "' stiffness";
        exit(1);
    };
}

//////////////////////////////////////////////////////////
void CSLMaterialStatus :: computeStress(double timeStep) {
    computeConstitutiveStrain();
    computeDamage( temp_strain );
    CSLMaterialStatus :: computeStressWithFrozenIntVars(timeStep);
}

//////////////////////////////////////////////////////////
void CSLMaterialStatus :: computeStressWithFrozenIntVars(double timeStep) {
    computeConstitutiveStrain();
    VectMechMaterialStatus :: computeStressWithFrozenIntVars(timeStep);
    temp_stress *=   1. - temp_damage;
    if ( temp_strain [ 0 ] > 0 ) {
        temp_crackOpening = ( L * temp_damage ) * temp_strain [ 0 ]; //normal opening only
        //temp_crackOpening = l2_norm( ( L * temp_damage ) * temp_strain );  //total opening
    } else {
        temp_crackOpening = 0;
    }
}

//////////////////////////////////////////////////////////
Vector CSLMaterialStatus :: giveCrackOpeningVector() const {
    return ( L * temp_damage ) * temp_strain;
}
//////////////////////////////////////////////////////////
std :: string CSLMaterialStatus :: giveLineToSave() const {
    return "damage " + to_string_sci(this->damage) + " maxEpsN " + to_string_sci(this->maxEpsN) + " maxEpsT " + to_string_sci(this->maxEpsT);
}

//////////////////////////////////////////////////////////
void CSLMaterialStatus :: setParameterValue(string code, double value) {
    VectMechMaterialStatus :: setParameterValue(code, value);
}

//////////////////////////////////////////////////////////
void CSLMaterialStatus :: readFromLine(istringstream &iss) {
    VectMechMaterialStatus :: readFromLine(iss);

    iss.clear(); // clear string stream
    iss.seekg(0, iss.beg); //reset position in string stream

    std :: string param;
    while (  iss >> param ) {
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
bool CSLMaterialStatus :: isElastic(const bool &now) const {
    if ( now && this->temp_damage != 0.0 ) {
        return false;
    } else if ( this->damage != 0.0 ) {
        return false;
    }
    return true;
}


//////////////////////////////////////////////////////////
// CSL MATERIAL

//////////////////////////////////////////////////////////
void CSLMaterial :: readFromLine(istringstream &iss) {
    VectMechMaterial :: readFromLine(iss); //read elastic parameters

    iss.clear(); // clear string stream
    iss.seekg(0, iss.beg); //reset position in string stream

    // initialize all values to zero (NOTE probably no ned in linux, but in windows necessary)
    fs = Gs = fc = Kc = 0;

    string param;
    bool bft, bGt;
    bft = bGt = false;

    while (  iss >> param ) {
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
        } else if ( param.compare("lambda0") == 0 ) {
            iss >> lam0;
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
    peakTensileStrain = ft/E0;
};

//////////////////////////////////////////////////////////
MaterialStatus *CSLMaterial :: giveNewMaterialStatus(Element *e, unsigned ipnum) {
    CSLMaterialStatus *newStatus = new CSLMaterialStatus(this, e, ipnum); //needs to be deleted manually
    return newStatus;
};


//////////////////////////////////////////////////////////
void CSLMaterial :: init(MaterialContainer *matcont) {
    VectMechMaterial :: init(matcont);

    // if variables not specified on the input, use default multipliers
    fs = ( fs == 0 ) ? 3 * ft : fs;
    Gs = ( Gs == 0 ) ? 16 * Gt : Gs;
    fc = ( fc == 0 ) ? 16 * ft : fc;
    Kc = ( Kc == 0 ) ? 0.26 * E0 : Kc;
    beta = 1.;
    mu = 0.2;
    nc = 2;

    Lcrt = 2 * E0 * Gt / pow(ft, 2);
    Lcrs = 2 * giveAlphaForDamage() * E0 * Gs / pow(fs, 2);
};


//////////////////////////////////////////////////////////
// CSL MATERIAL STATUS WITH TANSORIAL STIFFNESS UPDATE

CSLMaterialWithTensorialStressUpdateStatus :: CSLMaterialWithTensorialStressUpdateStatus(CSLMaterialWithTensorialStressUpdate *m, Element *e, unsigned ipnum) : CSLMaterialStatus(m, e, ipnum) {
    name = "CSL mat. with tensorial stress update status";
}

//////////////////////////////////////////////////////////
std :: unique_ptr< MaterialStatus > CSLMaterialWithTensorialStressUpdateStatus :: cloneState() const {
    std :: unique_ptr< CSLMaterialWithTensorialStressUpdateStatus > copy = std :: make_unique< CSLMaterialWithTensorialStressUpdateStatus >( *this );
    copy->clearSnapshotComponentPointers();
    return copy;
}

//////////////////////////////////////////////////////////
void CSLMaterialWithTensorialStressUpdateStatus :: restoreStateFrom(const MaterialStatus &other) {
    const CSLMaterialWithTensorialStressUpdateStatus *otherCSL = dynamic_cast< const CSLMaterialWithTensorialStressUpdateStatus * >( &other );
    if ( !otherCSL ) {
        std :: cerr << "Material status snapshot restore type mismatch for " << name << std :: endl;
        exit(1);
    }
    *this = *otherCSL;
    clearSnapshotComponentPointers();
}

//////////////////////////////////////////////////////////
bool CSLMaterialWithTensorialStressUpdateStatus :: giveValues(string code, Vector &result) const {
    //if ( code.compare("tempCrackOpening") == 0 || code.compare("crack_opening") == 0 ) {
    return CSLMaterialStatus :: giveValues(code, result);
    //}
}

//////////////////////////////////////////////////////////
Vector CSLMaterialWithTensorialStressUpdateStatus :: giveEigenStrainFromTensorialStress() const {
    CSLMaterialWithTensorialStressUpdate *m = static_cast< CSLMaterialWithTensorialStressUpdate * >( mat );
    Vector ts = m->giveAveragePrincipalStress( element->giveNode(0)->giveID(), element->giveNode(1)->giveID() );

    Vector eigenvalues;
    vector< Vector >eigenvectors;
    LinalgEigenSolver(ts, eigenvalues, eigenvectors);
    unsigned dim = m->giveDimension();
    Vector augstress(dim);

    Matrix dirs(dim, dim);
    for ( unsigned i = 0; i < dim; i++ ) {
        dirs.row(i) = eigenvectors [ i ];
    }
    Vector res(dim);
    if ( dim == 2 ) {
        res [ 0 ] = eigenvalues [ 1 ];
        res [ 1 ] = eigenvalues [ 0 ];
    } else if ( dim == 3 ) {
        res [ 0 ] = eigenvalues [ 1 ] + eigenvalues [ 2 ];
        res [ 1 ] = eigenvalues [ 0 ] + eigenvalues [ 2 ];
        res [ 2 ] = eigenvalues [ 0 ] + eigenvalues [ 1 ];
    }
    RigidBodyContact *rbc = static_cast< RigidBodyContact * >( element );
    Vector p = ( ( dirs.transpose() * res.asDiagonal() * dirs ) * rbc->giveNormal() ) * ( m->givePoissonNumber() / ( max(1. - temp_damage, 0.00001) * m->giveE0() ) );
    Vector q(dim);
    q [ 0 ] = rbc->giveNormal().dot(p);
    q [ 1 ] = rbc->giveT1().dot(p);
    if ( dim == 3 ) {
        q [ 2 ] = rbc->giveT2().dot(p);
    }
    return q;
}


//////////////////////////////////////////////////////////
void CSLMaterialWithTensorialStressUpdateStatus :: computeStress(double timeStep) {
    temp_strain += giveEigenStrainFromTensorialStress();
    CSLMaterialStatus :: computeStress(timeStep);
}

//////////////////////////////////////////////////////////
void CSLMaterialWithTensorialStressUpdateStatus :: computeStressWithFrozenIntVars(double timeStep) {
    temp_strain += giveEigenStrainFromTensorialStress();
    CSLMaterialStatus :: computeStressWithFrozenIntVars(timeStep);
}

//////////////////////////////////////////////////////////
// CSL MATERIAL WITH TENSORIAL STRESS UPDATE
//////////////////////////////////////////////////////////

CSLMaterialWithTensorialStressUpdate :: CSLMaterialWithTensorialStressUpdate(unsigned dimension, CSLMaterialWithTensorialStressUpdate *mm) : CSLMaterial(dimension) {
    name = "CSL material with tensorial stress update";
    master_material = mm;
};

//////////////////////////////////////////////////////////
void CSLMaterialWithTensorialStressUpdate :: prepareForStressEvaluation(ElementContainer *elems) {
    if ( not master_material ) {
        tensstress = elems->computePrincipalStresses();
    }
    //else it is already computed by master
}


//////////////////////////////////////////////////////////
Vector CSLMaterialWithTensorialStressUpdate :: giveAveragePrincipalStress(unsigned Anode, unsigned Bnode) const {
    if ( not master_material ) {
        return ( tensstress [ Anode ] + tensstress [ Bnode ] ) / 2.;
    } else {
        return master_material->giveAveragePrincipalStress(Anode, Bnode);
    }
}

//////////////////////////////////////////////////////////
void CSLMaterialWithTensorialStressUpdate :: readFromLine(istringstream &iss) {
    CSLMaterial :: readFromLine(iss); //read elastic parameters

    iss.clear(); // clear string stream
    iss.seekg(0, iss.beg); //reset position in string stream



    string param;
    bool bnu;
    bnu = false;

    while (  iss >> param ) {
        if ( param.compare("nu") == 0 ) {
            bnu = true;
            iss >> poisson;
        }
    }
    if ( !bnu ) {
        cerr << name << ": material parameter 'nu' was not specified" << endl;
        exit(EXIT_FAILURE);
    }
};

//////////////////////////////////////////////////////////
MaterialStatus *CSLMaterialWithTensorialStressUpdate :: giveNewMaterialStatus(Element *e, unsigned ipnum) {
    CSLMaterialWithTensorialStressUpdateStatus *newStatus = new CSLMaterialWithTensorialStressUpdateStatus(this, e, ipnum); //needs to be deleted manually
    return newStatus;
};


//////////////////////////////////////////////////////////
void CSLMaterialWithTensorialStressUpdate :: init(MaterialContainer *matcont) {
    alphaForDamage = alpha;
    alpha = 1;
    CSLMaterial :: init(matcont);
};



//////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////
// COUPLED CSL MATERIAL STATUS
//////////////////////////////////////////////////////////
CoupledCSLMaterialStatus :: CoupledCSLMaterialStatus(CSLMaterial *m, Element *e, unsigned ipnum) : CSLMaterialStatus(m, e, ipnum) {
    name = "Coupled CSL mat. status";
}

//////////////////////////////////////////////////////////
std :: unique_ptr< MaterialStatus > CoupledCSLMaterialStatus :: cloneState() const {
    std :: unique_ptr< CoupledCSLMaterialStatus > copy = std :: make_unique< CoupledCSLMaterialStatus >( *this );
    copy->clearSnapshotComponentPointers();
    return copy;
}

//////////////////////////////////////////////////////////
void CoupledCSLMaterialStatus :: restoreStateFrom(const MaterialStatus &other) {
    const CoupledCSLMaterialStatus *otherCSL = dynamic_cast< const CoupledCSLMaterialStatus * >( &other );
    if ( !otherCSL ) {
        std :: cerr << "Material status snapshot restore type mismatch for " << name << std :: endl;
        exit(1);
    }
    *this = *otherCSL;
    clearSnapshotComponentPointers();
}

//////////////////////////////////////////////////////////
std :: uint64_t CoupledCSLMaterialStatus :: stateHash() const {
    std :: uint64_t hash = CSLMaterialStatus :: stateHash();
    hashDoubleCSL(hash, avgPressure);
    return hash;
}

//////////////////////////////////////////////////////////
void CoupledCSLMaterialStatus :: init() {
    CSLMaterialStatus :: init();
    avgPressure = 0;
}


//////////////////////////////////////////////////////////
void CoupledCSLMaterialStatus :: setParameterValue(string code, double value) {
    if ( code.compare("pressure") == 0 ) {
        avgPressure = value;
    } else {
        CSLMaterialStatus :: setParameterValue(code, value);
    }
}

//////////////////////////////////////////////////////////
bool CoupledCSLMaterialStatus :: giveValues(string code, Vector &result) const {
    if ( code.compare("pressure") == 0 || code.compare("avg_pressure") == 0 ) {
        result.resize(1);
        result [ 0 ] = avgPressure;
        return true;
    } else if ( code.compare("solid_stress") == 0 ) {
        CSLMaterialStatus :: giveValues("stress", result); //standard stress including Biot's effect
        CoupledCSLMaterial *m = static_cast< CoupledCSLMaterial * >( mat );
        result [ 0 ] += m->giveBiotCoeff() * avgPressure; //stress without Biot's effect
        return true;
    } else {
        return CSLMaterialStatus :: giveValues(code, result);
    }
}



//////////////////////////////////////////////////////////
void CoupledCSLMaterialStatus :: updateStressByBiotEffect(double timeStep) {
    ( void ) timeStep;
    CoupledCSLMaterial *m = static_cast< CoupledCSLMaterial * >( mat );
    temp_stress [ 0 ] -= m->giveBiotCoeff() * avgPressure;
}

//////////////////////////////////////////////////////////
void CoupledCSLMaterialStatus :: computeStress(double timeStep) {
    CSLMaterialStatus :: computeStress(timeStep);
    updateStressByBiotEffect(timeStep);
}

//////////////////////////////////////////////////////////
void CoupledCSLMaterialStatus :: computeStressWithFrozenIntVars(double timeStep) {
    CSLMaterialStatus :: computeStressWithFrozenIntVars(timeStep);
    updateStressByBiotEffect(timeStep);
}

//////////////////////////////////////////////////////////
void CoupledCSLMaterialStatus :: update() {
    CSLMaterialStatus :: update();
}

//////////////////////////////////////////////////////////
void CoupledCSLMaterialStatus :: resetTemporaryVariables() {
    CSLMaterialStatus :: resetTemporaryVariables();
}

//////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////
// COUPLED CSL MATERIAL
//////////////////////////////////////////////////////////
void CoupledCSLMaterial :: init(MaterialContainer *matcont) {
    CSLMaterial :: init(matcont);
}

//////////////////////////////////////////////////////////
MaterialStatus *CoupledCSLMaterial :: giveNewMaterialStatus(Element *e, unsigned ipnum) {
    CoupledCSLMaterialStatus *newStatus = new CoupledCSLMaterialStatus(this, e, ipnum);
    return newStatus;
};

//////////////////////////////////////////////////////////
void CoupledCSLMaterial :: readFromLine(istringstream &iss) {
    CSLMaterial :: readFromLine(iss);

    iss.clear(); // clear string stream
    iss.seekg(0, iss.beg); //reset position in string stream

    string param;
    bool bbiot = false;

    while (  iss >> param ) {
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
