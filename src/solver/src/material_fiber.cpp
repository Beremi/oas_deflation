#include "material_fiber.h"
#include "element_fiber.h"
#include "element_discrete.h"
#include "material_csl.h"

using namespace std;

//////////////////////////////////////////////////////////
// FIBER MATERIAL STATUS

FiberMaterialStatus :: FiberMaterialStatus(FiberMaterial *m, Element *e, unsigned ipnum) : TensMechMaterialStatus(m, e, ipnum) {
    name = "fiber mat. status";
    crack_opening = 0;
    temp_crack_opening = 0;
    incrementOfCrack = 0;
    maxCrackOpening = 0;
    crackOpeningVector = Vector :: Zero( contactNormal.size() );

    right_pullout = 0;
    temp_rightPullout = 0;
    left_pullout = 0;
    temp_leftPullout = 0;

    bridgingForce = 0;
    temp_bridgingForce = 0;
    rightForce = leftForce = 0;
    temp_stress = Vector :: Zero( contactNormal.size() );

    spallingLength = 0;
    w = nf = Vector :: Zero( contactNormal.size() );
    deflectionAngle = 0;
    fiberForce = 0;

    debondedFiber = 0;
    temp_debondedFiber = 0;
    pulloutOfFiber = 0;
    temp_pulloutOfFiber = 0;
    ruptureOfFiber = 0;
    temp_ruptureOfFiber = 0;
    closingCrack = 0;
    temp_closingCrack = 0;
}

//////////////////////////////////////////////////////////
bool FiberMaterialStatus :: giveValues(string code, Vector &result) const {
    if ( code.compare("tempLeftPullout") == 0 || code.compare("left_pullout") == 0 ) {
        result.resize(1);
        result [ 0 ] = temp_leftPullout;
        return true;
    } else if ( code.compare("tempRightPullout") == 0 || code.compare("right_pullout") == 0 ) {
        result.resize(1);
        result [ 0 ] = temp_rightPullout;
        return true;
    } else if ( code.compare("tempBridgingForce") == 0 || code.compare("bridging_force") == 0 ) {
        result.resize(1);
        result [ 0 ] = temp_bridgingForce;
        return true;
    } else {
        return TensMechMaterialStatus :: giveValues(code, result);
    }
}

//////////////////////////////////////////////////////////
void FiberMaterialStatus :: init() {
    //  temp_stress = Vector :: Zero(contactNormal.size() );
    //  crackOpeningVector = Vector :: Zero(contactNormal.size() );

    Fiber *fibElement = static_cast< Fiber * >( element );
    df = fibElement->giveDiameter();
    rightLe = fibElement->giveRightLength(idx);
    leftLe = fibElement->giveLeftLength(idx);

    fiberNormal = fibElement->giveDirVector();

    RigidBodyContact *mechElement = static_cast< RigidBodyContact * >( fibElement->giveRBContact(idx) );
    contactLength = mechElement->giveLength();
    contactNormal = mechElement->giveNormal();

    fiberNormalLocal = mechElement->transformToLocal(fiberNormal);

    FiberMaterial *fibMaterial = static_cast< FiberMaterial * >( mat );
    double Ef = fibMaterial->giveEf();
    double tau0 = fibMaterial->giveTau0();
    double Gd = fibMaterial->giveGd();
    double betaf = fibMaterial->giveBetaf();

    right_F0 = M_PI * df * tau0 * rightLe;
    left_F0 = M_PI * df * tau0 * leftLe;
    limit_rightPullout = 2. * tau0 * pow(rightLe, 2) / ( Ef * df ) + sqrt(8. * Gd * pow(rightLe, 2) / ( Ef * df ) );
    limit_leftPullout = 2. * tau0 * pow(leftLe, 2) / ( Ef * df ) + sqrt(8. * Gd * pow(leftLe, 2) / ( Ef * df ) );

    if ( rightLe > leftLe ) {
        Le1 = leftLe;
        F01 = left_F0;
        vd1 = limit_leftPullout;
        Le2 = rightLe;
        F02 = right_F0;
        vd2 = limit_rightPullout;
    } else if ( rightLe <= leftLe ) {
        Le1 = rightLe;
        F01 = right_F0;
        vd1 = limit_rightPullout;
        Le2 = leftLe;
        F02 = left_F0;
        vd2 = limit_leftPullout;
    }
    // pullout v0 -> debonded curve reaches F=0
    v0 = Le1 + vd1;
    if ( betaf > 0 ) {
        A = -betaf;
        B = Le1 * betaf - df + 2. * betaf * vd1;
        C = Le1 * df - Le1 * betaf * vd1 + vd1 * df - betaf * pow(vd1, 2);
        v_Fmax = -B / ( 2. * A ); // pullout where parabole reaches the max. bridging force
        Fmax = bridgingForce_debonded(v_Fmax, vd1, Le1, F01, df, betaf);

        AA = -betaf;
        BB = Le2 * betaf - df + 2. * betaf * vd2;
        CC = Le2 * df - Le2 * betaf * vd2 + vd2 * df - betaf * pow(vd2, 2) - Fmax * Le2 * df / F02;
        v21 = ( -BB + sqrt(pow(BB, 2) - 4. * AA * CC) ) / ( 2. * AA );
        v22 = ( -BB - sqrt(pow(BB, 2) - 4. * AA * CC) ) / ( 2. * AA );
        if ( v21 <= v22 ) {
            v2Debonded_Fmax = v21;
        } else {
            v2Debonded_Fmax = v22;
        }
    }
}

//////////////////////////////////////////////////////////
void FiberMaterialStatus :: update() {
    crack_opening = temp_crack_opening;

    if ( crack_opening > maxCrackOpening ) {
        maxCrackOpening = crack_opening;
    } else {
        maxCrackOpening = maxCrackOpening;
    }

    right_pullout = temp_rightPullout;
    left_pullout = temp_leftPullout;
    bridgingForce = temp_bridgingForce;
    debondedFiber = temp_debondedFiber;
    pulloutOfFiber = temp_pulloutOfFiber;
    ruptureOfFiber = temp_ruptureOfFiber;
    closingCrack = temp_closingCrack;

    MaterialStatus :: update();
}

//////////////////////////////////////////////////////////
void FiberMaterialStatus :: resetTemporaryVariables() {
    MaterialStatus :: resetTemporaryVariables();
}

//////////////////////////////////////////////////////////
double bridgingForce_bonded(double v, double vd, double df, double Ef, double tau0, double Gd) {
    return sqrt(0.5 * pow(M_PI, 2) * Ef * pow(df, 3) * ( tau0 * v + Gd * v / vd ) );
}
double derivF_bonded(double v, double vd, double df, double Ef, double tau0, double Gd) {
    return 0.5 / sqrt(0.5 * pow(M_PI, 2) * Ef * pow(df, 3) * ( tau0 * v + Gd * v / vd ) ) * ( 0.5 * pow(M_PI, 2) * Ef * pow(df, 3) * ( tau0 + Gd / vd ) );
}
double bridgingForce_debonded(double v, double vd, double Le, double F0, double df, double betaf) {
    return F0 * ( 1 - ( v - vd ) / Le ) * ( 1 +  betaf * ( v - vd ) / df );
}
double derivF_debonded(double v, double vd, double Le, double F0, double df, double betaf) {
    return F0 / ( Le * df ) * ( Le * betaf - df - 2. * v * betaf + 2. * vd * betaf );
}
double bridgingForce_unloading(double deltaX, double deltaY, double x) {   // linear curve to zero
    return deltaY / deltaX * x;
}
double derivF_unloading(double deltaX, double deltaY) {
    return deltaY / deltaX;
}

//////////////////////////////////////////////////////////
Matrix FiberMaterialStatus :: giveStiffnessTensor(string type) const {
    ( void ) type;
    unsigned ss = mat->giveStrainSize();

    Matrix stiffness = Matrix :: Zero(ss, ss);

    if ( temp_bridgingForce < 1e-12 ) {
        double alpha1, alpha2;
        alpha1 = temp_bridgingForce / ( temp_rightPullout / contactLength );
        alpha2 = temp_bridgingForce / ( temp_leftPullout / contactLength );
        stiffness(0, 0) = stiffness(1, 1) = stiffness(2, 2) = alpha1 * alpha2 / ( alpha1 + alpha2 );
    }

    //cout << "FiberMaterialStatus::giveStiffnessTensor" << endl;
    cout.flush();
    return stiffness;
}

//////////////////////////////////////////////////////////
Vector FiberMaterialStatus :: giveStress(const Vector &strain, double timeStep) {
    ( void ) strain;
    ( void ) timeStep;

    FiberMaterial *fibMaterial = static_cast< FiberMaterial * >( mat );
    double Ef = fibMaterial->giveEf();
    double Gd = fibMaterial->giveGd();
    double tau0 = fibMaterial->giveTau0();
    double betaf = fibMaterial->giveBetaf();
    double ft = fibMaterial->giveFt();
    double Ksn = fibMaterial->giveKsn();   // snubbing coefficient
    double Ksp = fibMaterial->giveKsp();   // spalling coefficient
    double Krup = fibMaterial->giveKrup();   // rupture coefficien

    Fiber *fibElement = static_cast< Fiber * >( element );
    CSLMaterialStatus *CSLMatStatus = static_cast< CSLMaterialStatus * >( fibElement->giveRBContact(idx)->giveMatStatus(0) );
    crackOpeningVector = CSLMatStatus->giveCrackOpeningVector();

    if ( crackOpeningVector [ 0 ] < 0 ) {
        crackOpeningVector [ 0 ] = 0;
    }

    temp_crack_opening = crackOpeningVector.norm();
    incrementOfCrack = temp_crack_opening - crack_opening;

    int iteration = 0;
    int limit_iteration = 150;
    double tolerance = 1e-8;
    // ---------------------- FIBER PULLOUTS --------------------------------------------------
    temp_rightPullout = right_pullout + incrementOfCrack / 2.;
    temp_leftPullout = left_pullout + incrementOfCrack / 2.;
    // ---------------------- BRIDGING FORCE --------------------------------------------------
    if ( temp_crack_opening < 1e-12 ) { // fiber does not bridge the crack
        temp_bridgingForce = 0;
        temp_rightPullout = temp_leftPullout = 0;
    } else if ( pulloutOfFiber == 1 or ruptureOfFiber == 1 ) { // fiber was pulled out or ruptured in the previous step
        temp_bridgingForce = 0;
        temp_rightPullout = right_pullout;
        temp_leftPullout = left_pullout;
        if ( pulloutOfFiber == 1 ) {
            temp_pulloutOfFiber = 1;
        } else if ( ruptureOfFiber == 1 ) {
            temp_ruptureOfFiber = 1;
        }
    } else if ( incrementOfCrack > 0 ) { // crack is opening
        if ( closingCrack == 1 and maxCrackOpening >= temp_crack_opening ) { // crack was closing in past and now is opening again
            while ( iteration <= limit_iteration ) {
                rightForce = bridgingForce_unloading(right_pullout, bridgingForce, temp_rightPullout);
                leftForce = bridgingForce_unloading(left_pullout, bridgingForce, temp_leftPullout);
                if ( rightForce <= 0 or leftForce <= 0 ) { // check if the crack is not fully closed
                    temp_bridgingForce = 0;
                    temp_rightPullout = temp_leftPullout = 0;
                    temp_crack_opening = 0;
                    cout << "TOTO BY NEMALO NASTAT UZ VOBEC !!!" << endl;
                    exit(1);
                    break;
                }
                if ( abs(rightForce - leftForce) <= tolerance ) {
                    temp_bridgingForce = rightForce;
                    break;
                } else if ( abs(rightForce - leftForce) > tolerance ) {
                    temp_leftPullout = temp_leftPullout - ( bridgingForce_unloading(left_pullout, bridgingForce, temp_leftPullout) - bridgingForce_unloading(right_pullout, bridgingForce, temp_crack_opening - temp_leftPullout) ) / ( derivF_unloading(left_pullout, bridgingForce) + derivF_unloading(right_pullout, bridgingForce) );
                    temp_rightPullout = temp_crack_opening - temp_leftPullout;
                }
                iteration += 1;
            }
            if ( iteration > limit_iteration ) {
                cout << " WARNING1: Fiber force equilibrium was not found ! Limit iteration was reached ! " << endl;
                exit(1);
            }
        } else if ( temp_rightPullout <= limit_rightPullout and temp_leftPullout <= limit_leftPullout ) { // both sides of the fiber are still bonded
            temp_bridgingForce = bridgingForce_bonded(temp_rightPullout, limit_rightPullout, df, Ef, tau0, Gd);
        } else if ( abs(rightLe - leftLe) < 1e-12 ) { // both sides of the fiber are equal and fully debonded
            temp_bridgingForce = bridgingForce_debonded(temp_rightPullout, limit_rightPullout, rightLe, right_F0, df, betaf);
            if ( temp_rightPullout >= v0 or temp_bridgingForce <= 0 ) { // check if the fiber is not fully pulled out
                temp_bridgingForce = 0;
                temp_pulloutOfFiber = 1;
                if ( temp_rightPullout < v0 ) { // for betaf << 0
                    temp_rightPullout = temp_crack_opening / 2.;
                    temp_leftPullout = temp_crack_opening / 2.;
                } else {
                    temp_rightPullout = temp_leftPullout = v0;
                }
            }
        } else if ( rightLe != leftLe ) { // shorter side of the fiber (1) is debonded, longer (2) is unloading
            if ( rightLe < leftLe ) {
                v1 = temp_rightPullout;
                v2 = temp_leftPullout;
            } else if ( rightLe > leftLe ) {
                v1 = temp_leftPullout;
                v2 = temp_rightPullout;
            }

            if ( debondedFiber == 0 ) {
                while ( iteration <= limit_iteration ) { // Newton iteration v_new = v_old - f(v_old)/f'(v_old)
                    bridgingForce1 = bridgingForce_debonded(v1, vd1, Le1, F01, df, betaf);
                    bridgingForce2 = bridgingForce_unloading(vd1, bridgingForce_bonded(vd1, vd1, df, Ef, tau0, Gd), v2);
                    if ( v1 >= v0 or bridgingForce1 <= 0 ) { // check if the shorter side of the fiber is not fully pulled out
                        temp_bridgingForce = 0;
                        temp_pulloutOfFiber = 1;
                        if ( v1 < v0 ) { // for betaf << 0
                            v1 = temp_crack_opening;
                            v2 = 0;
                        } else {
                            v1 = v0;
                            v2 = 0;
                        }
                        break;
                    } else if ( abs(bridgingForce1 - bridgingForce2) <= tolerance ) {
                        temp_bridgingForce = bridgingForce1;
                        if ( temp_bridgingForce > bridgingForce_bonded(vd1, vd1, df, Ef, tau0, Gd) ) {   // longer (2) is bonded
                            temp_debondedFiber = 1;
                        }
                        break;
                    } else if ( abs(bridgingForce1 - bridgingForce2) > tolerance ) {
                        v1 = v1 - ( bridgingForce_debonded(v1, vd1, Le1, F01, df, betaf) - bridgingForce_unloading(vd1, bridgingForce_bonded(vd1, vd1, df, Ef, tau0, Gd), temp_crack_opening - v1) ) / ( derivF_debonded(v1, vd1, Le1, F01, df, betaf) + derivF_unloading(vd1, bridgingForce_bonded(vd1, vd1, df, Ef, tau0, Gd) ) );
                        v2 = temp_crack_opening - v1;
                    }
                    iteration += 1;
                }
            }

            if ( debondedFiber == 1 or temp_debondedFiber == 1 ) {
                while ( iteration <= limit_iteration ) { // Newton iteration v_new = v_old - f(v_old)/f'(v_old)
                    bridgingForce1 = bridgingForce_debonded(v1, vd1, Le1, F01, df, betaf);
                    bridgingForce2 = bridgingForce_bonded(v2, vd2, df, Ef, tau0, Gd);
                    if ( abs(bridgingForce1 - bridgingForce2) <= tolerance ) {
                        temp_bridgingForce = bridgingForce1;
                        temp_debondedFiber = 1;
                        if ( temp_bridgingForce > bridgingForce_bonded(vd2, vd2, df, Ef, tau0, Gd) ) {
                            temp_debondedFiber = 2;
                        } else if ( v1 > v_Fmax ) {
                            temp_debondedFiber = 3;
                        }
                        break;
                    } else if ( abs(bridgingForce1 - bridgingForce2) > tolerance ) {
                        v1 = v1 - ( bridgingForce_debonded(v1, vd1, Le1, F01, df, betaf) - bridgingForce_bonded(temp_crack_opening - v1, vd2, df, Ef, tau0, Gd) ) / ( derivF_debonded(v1, vd1, Le1, F01, df, betaf) + derivF_bonded(temp_crack_opening - v1, vd2, df, Ef, tau0, Gd) );
                        v2 = temp_crack_opening - v1;
                    }
                    iteration += 1;
                }
            }

            if ( debondedFiber == 2 or temp_debondedFiber == 2 ) {
                while ( iteration <= limit_iteration ) { // Newton iteration v_new = v_old - f(v_old)/f'(v_old)
                    bridgingForce1 = bridgingForce_debonded(v1, vd1, Le1, F01, df, betaf);
                    bridgingForce2 = bridgingForce_debonded(v2, vd2, Le2, F02, df, betaf);
                    if ( abs(bridgingForce1 - bridgingForce2) <= tolerance ) {
                        temp_bridgingForce = bridgingForce1;
                        temp_debondedFiber = 2;
                        if ( v1 > v_Fmax ) {
                            temp_debondedFiber = 4;
                        }
                        break;
                    } else if ( abs(bridgingForce1 - bridgingForce2) > tolerance ) {
                        v1 = v1 - ( bridgingForce_debonded(v1, vd1, Le1, F01, df, betaf) - bridgingForce_debonded(temp_crack_opening - v1, vd2, Le2, F02, df, betaf) ) / ( derivF_debonded(v1, vd1, Le1, F01, df, betaf) + derivF_debonded(temp_crack_opening - v1, vd2, Le2, F02, df, betaf) );
                        v2 = temp_crack_opening - v1;
                    }
                    iteration += 1;
                }
            }

            if ( debondedFiber == 3 or temp_debondedFiber == 3 ) {
                while ( iteration <= limit_iteration ) { // Newton iteration v_new = v_old - f(v_old)/f'(v_old)
                    v2Bonded_Fmax = pow(Fmax, 2) * 2. / ( pow(M_PI, 2) * Ef * pow(df, 3) * ( tau0 + Gd / vd2 ) );
                    bridgingForce1 = bridgingForce_debonded(v1, vd1, Le1, F01, df, betaf);
                    bridgingForce2 = bridgingForce_unloading(v2Bonded_Fmax, Fmax, v2);
                    if ( bridgingForce1 <= 0 ) { // check if the shorter side of the fiber is not fully pulled out
                        temp_bridgingForce = 0;
                        v1 = v0;
                        v2 = 0;
                        temp_pulloutOfFiber = 1;
                        break;
                    } else if ( abs(bridgingForce1 - bridgingForce2) <= tolerance ) {
                        temp_bridgingForce = bridgingForce1;
                        temp_debondedFiber = 3;
                        break;
                    } else if ( abs(bridgingForce1 - bridgingForce2) > tolerance ) {
                        v1 = v1 - ( bridgingForce_debonded(v1, vd1, Le1, F01, df, betaf) - bridgingForce_unloading(v2Bonded_Fmax, Fmax, temp_crack_opening - v1) ) / ( derivF_debonded(v1, vd1, Le1, F01, df, betaf) + derivF_unloading(v2Bonded_Fmax, Fmax) );
                        v2 = temp_crack_opening - v1;
                    }
                    iteration += 1;
                }
            }

            if ( debondedFiber == 4 or temp_debondedFiber == 4 ) {
                while ( iteration <= limit_iteration ) { // Newton iteration v_new = v_old - f(v_old)/f'(v_old)
                    bridgingForce1 = bridgingForce_debonded(v1, vd1, Le1, F01, df, betaf);
                    bridgingForce2 = bridgingForce_unloading(v2Debonded_Fmax, Fmax, v2);
                    if ( v1 >= v0 or bridgingForce1 <= 0 ) { // check if the shorter side of the fiber is not fully pulled out
                        temp_bridgingForce = 0;
                        v1 = v0;
                        v2 = 0;
                        temp_pulloutOfFiber = 1;
                        break;
                    } else if ( abs(bridgingForce1 - bridgingForce2) <= tolerance ) {
                        temp_bridgingForce = bridgingForce1;
                        temp_debondedFiber = 4;
                        break;
                    } else if ( abs(bridgingForce1 - bridgingForce2) > tolerance ) {
                        v1 = v1 - ( bridgingForce_debonded(v1, vd1, Le1, F01, df, betaf) - bridgingForce_unloading(v2Debonded_Fmax, Fmax, temp_crack_opening - v1) ) / ( derivF_debonded(v1, vd1, Le1, F01, df, betaf) + derivF_unloading(v2Debonded_Fmax, Fmax) );
                        v2 = temp_crack_opening - v1;
                    }
                    iteration += 1;
                }
            }

            if ( iteration > limit_iteration ) {
                cout << " WARNING2: Fiber force equilibrium was not found ! Limit iteration was reached ! " << endl;
                exit(1);
            } else if ( rightLe < leftLe ) {
                temp_rightPullout = v1;
                temp_leftPullout = v2;
            } else if ( rightLe > leftLe ) {
                temp_leftPullout = v1;
                temp_rightPullout = v2;
            }
        }
    } else if ( incrementOfCrack <= 0 ) { // crack is closing or stayed as opened as in the previous step
        temp_closingCrack = 1;
        while ( iteration <= limit_iteration ) {
            rightForce = bridgingForce_unloading(right_pullout, bridgingForce, temp_rightPullout);
            leftForce = bridgingForce_unloading(left_pullout, bridgingForce, temp_leftPullout);
            if ( abs(rightForce - leftForce) <= tolerance ) {
                temp_bridgingForce = rightForce;
                if ( temp_bridgingForce < 0 ) { // check if the crack is not fully closed
                    temp_bridgingForce = 0;
                    temp_rightPullout = temp_leftPullout = 0;
                    temp_crack_opening = 0;
                    cout << "TOTO BY NEMALO NASTAT !!!" << endl;
                    exit(1);
                }
                break;
            } else if ( abs(rightForce - leftForce) > tolerance ) {
                temp_leftPullout = temp_leftPullout - ( bridgingForce_unloading(left_pullout, bridgingForce, temp_leftPullout) - bridgingForce_unloading(right_pullout, bridgingForce, temp_crack_opening - temp_leftPullout) ) / ( derivF_unloading(left_pullout, bridgingForce) + derivF_unloading(right_pullout, bridgingForce) );
                temp_rightPullout = temp_crack_opening - temp_leftPullout;
            }
            iteration += 1;
        }
        if ( iteration > limit_iteration ) {
            cout << " WARNING3: Fiber force equilibrium was not found ! Limit iteration was reached ! " << endl;
            exit(1);
        }
    }
    // -------------------------- STRESS VECTOR ---------------------------------------------------

    inclineAngle = acos( ( fiberNormal ).dot(contactNormal) );
    if ( inclineAngle > M_PI / 2. and temp_bridgingForce > 0 ) {
        contactNormal = contactNormal * ( -1. );
        inclineAngle = acos( ( fiberNormal ).dot(contactNormal) );
    }

    if ( temp_bridgingForce < 1e-12 ) {
        temp_stress = Vector :: Zero( contactNormal.size() );
    } else if ( inclineAngle < 1e-12 ) { // WITHOUT SNUBBING MODEL - NO MICROEFFECT AT EXIT POINTS
        temp_stress = temp_bridgingForce * crackOpeningVector / temp_crack_opening;
    } else if ( inclineAngle != 0 ) { // SNUBBING MODEL - CONSIDERING MICROEFFECT AT EXIT POINTS
        iteration = 0;
        while ( iteration <= limit_iteration ) {
            if ( fiberNormalLocal [ 0 ] < 0 ) {
                fiberNormalLocal = fiberNormalLocal * ( -1. );
            }

            w = crackOpeningVector + 2. * spallingLength * fiberNormalLocal;

            nf = w / w.norm();

            deflectionAngle = acos(fiberNormalLocal.dot(nf) );

            fiberForce = temp_bridgingForce * exp(Ksn * deflectionAngle);

            temp_stress = fiberForce * nf;

            temp_spallingLength = ( temp_stress(0) * sin(inclineAngle / 2.) ) / ( Ksp * ft * df * pow(cos(inclineAngle / 2.), 2) );

            if ( ( abs(temp_spallingLength - spallingLength) ) < 1e-12 ) { // ABSOLUTE ERROR OF SPALLING LENGTH
                temp_bridgingForce = fiberForce;
                break;
            } else {
                spallingLength = temp_spallingLength;
            }
            iteration += 1;
        }

        if ( iteration > limit_iteration ) {
            cout << " FIBER ERROR: LIMIT ITERATION FOR SNUBBING MODEL WAS REACHED ! " << endl;
            cout << " temp_bridgingForce " << temp_bridgingForce << endl;
            cout << " fiberNormal " << fiberNormal << endl;
            cout << " contactNormal " << contactNormal << endl;
            cout << " fiberNormalLocal " << fiberNormalLocal << endl;
            cout << " crackOpeningVector " << crackOpeningVector << endl;
            cout << " " << endl;
            cout << " inclineAngle " << inclineAngle << endl;
            cout << " spallingLength " << spallingLength << endl;
            cout << " w " << w << endl;
            cout << " nf " << nf << endl;
            cout << " deflectionAngle " << deflectionAngle << endl;
            cout << " fiberForce " << fiberForce << endl;
            cout << " temp_stress " << temp_stress << endl;

            exit(1);
        }
    }
    // -------------------------- RUPTURE CONDITION ---------------------------------------------------
    double fiberStress = 4. * temp_bridgingForce / M_PI / pow(df, 2);
    double limitStress = ft * exp(-Krup * deflectionAngle);
    if ( fiberStress > limitStress ) {
        temp_bridgingForce = 0;
        temp_stress = Vector :: Zero( contactNormal.size() );
        temp_ruptureOfFiber = 1;
    }

    //cout << "FiberMaterialStatus::giveStress" << endl;
    cout.flush();
    return temp_stress;
}

//////////////////////////////////////////////////////////
Vector FiberMaterialStatus :: giveStressWithFrozenIntVars(const Vector &strain, double timeStep) {
    ( void ) strain;
    ( void ) timeStep;

    temp_crack_opening = strain.norm() * contactLength;

    Vector stressWithFrozenIntVars = giveStiffnessTensor("elastic") * temp_crack_opening;

    //cout << "FiberMaterialStatus::giveStressWithFrozenIntVars" << endl;
    cout.flush();
    return stressWithFrozenIntVars;
}

//////////////////////////////////////////////////////////
void FiberMaterialStatus :: setParameterValue(string code, double value) {
    if ( code.compare("crack_opening") == 0 ) {
        crack_opening = value;
    } else {
        MaterialStatus :: setParameterValue(code, value);
    }
}

//////////////////////////////////////////////////////////
// FIBER MATERIAL

FiberMaterial :: FiberMaterial(unsigned dimension) : TensMechMaterial(dimension) {
    name = "fiber material";
    // if variables not specified on the input, use default multipliers
    Ef = Gd = tau0 = betaf = ft = Ksn = Ksp = Krup = 0;
}

//////////////////////////////////////////////////////////
void FiberMaterial :: readFromLine(istringstream &iss) {
    Material :: readFromLine(iss); //read elastic parameters

    iss.clear(); // clear string stream
    iss.seekg(0, iss.beg); //reset position in string stream

    // initialize all values to zero (NOTE probably no need in linux, but in windows necessary)
    Ef = Gd = tau0 = betaf = ft = Ksn = Ksp = Krup = 0;

    string param;
    bool bft, bGd, btau0, bbetaf, bEf, bKsn, bKsp, bKrup;
    bft = bGd = btau0 = bbetaf = bEf = bKsn = bKsp = bKrup = false;

    while (  iss >> param ) {
        if ( param.compare("Ef") == 0 ) {
            bEf = true;
            iss >> Ef;
        } else if ( param.compare("Gd") == 0 ) {
            bGd = true;
            iss >> Gd;
        } else if ( param.compare("tau0") == 0 ) {
            btau0 = true;
            iss >> tau0;
        } else if ( param.compare("ft") == 0 ) {
            bft = true;
            iss >> ft;
        } else if ( param.compare("betaf") == 0 ) {
            bbetaf = true;
            iss >> betaf;
        } else if ( param.compare("Ksn") == 0 ) {
            bKsn = true;
            iss >> Ksn;
        } else if ( param.compare("Ksp") == 0 ) {
            bKsp = true;
            iss >> Ksp;
        } else if ( param.compare("Krup") == 0 ) {
            bKrup = true;
            iss >> Krup;
        }
    }
    if ( !bEf ) {
        cerr << name << ": material parameter 'Ef' was not specified" << endl;
        exit(EXIT_FAILURE);
    }
    ;
    if ( !bGd ) {
        cerr << name << ": material parameter 'Gd' was not specified" << endl;
        exit(EXIT_FAILURE);
    }
    ;
    if ( !btau0 ) {
        cerr << name << ": material parameter 'tau0' was not specified" << endl;
        exit(EXIT_FAILURE);
    }
    ;
    if ( !bft ) {
        cerr << name << ": material parameter 'ft' was not specified" << endl;
        exit(EXIT_FAILURE);
    }
    ;
    if ( !bbetaf ) {
        cerr << name << ": material parameter 'betaf' was not specified" << endl;
        exit(EXIT_FAILURE);
    }
    ;
    if ( !bKsn ) {
        cerr << name << ": material parameter 'Ksn' was not specified" << endl;
        exit(EXIT_FAILURE);
    }
    ;
    if ( !bKsp ) {
        cerr << name << ": material parameter 'Ksp' was not specified" << endl;
        exit(EXIT_FAILURE);
    }
    ;
    if ( !bKrup ) {
        cerr << name << ": material parameter 'Krup' was not specified" << endl;
        exit(EXIT_FAILURE);
    }
    ;
};

//////////////////////////////////////////////////////////
MaterialStatus *FiberMaterial :: giveNewMaterialStatus(Element *e, unsigned ipnum) {
    FiberMaterialStatus *newStatus = new FiberMaterialStatus(this, e, ipnum); //needs to be deleted manually
    return newStatus;
};

//////////////////////////////////////////////////////////
void FiberMaterial :: init(MaterialContainer *matcont) {
    TensMechMaterial :: init(matcont);
}
