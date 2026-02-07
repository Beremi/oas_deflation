#include "material_fiber.h"
#include "element_fiber.h"
#include "element_discrete.h"
#include "material_csl.h"

using namespace std;

//////////////////////////////////////////////////////////
// FIBER MATERIAL STATUS

FiberMaterialStatus :: FiberMaterialStatus(FiberMaterial *m, Element *e, unsigned ipnum) : TensMechMaterialStatus(m, e, ipnum) {
    name = "fiber mat. status";
    crackOpening = temp_crackOpening = maxCrackOpening = 0;
    incrementOfCrack = 0;
    crackOpeningVector = Vector :: Zero(contactNormal.size() );

    rightPullout = temp_rightPullout = 0;
    leftPullout = temp_leftPullout = 0;

    bridgingForce = temp_bridgingForce = 0;
    fiberForce = temp_fiberForce = 0;
    rightForce = leftForce = 0;
    temp_stress = Vector :: Zero(contactNormal.size() );

    deflectionAngle = spallingLength = temp_spallingLength = 0;
    w = mf = Vector :: Zero(contactNormal.size() );

    debondedFiber = temp_debondedFiber = 0;
    pulledOutFiber = temp_pulledOutFiber = 0;
    rupturedFiber = temp_rupturedFiber = 0;
    closingCrack = temp_closingCrack = 0;
}

//////////////////////////////////////////////////////////
bool FiberMaterialStatus :: giveValues(string code, Vector &result) const {
    if ( code.compare("templeftpullout") == 0 || code.compare("leftpullout") == 0 ) {
        result.resize(1);
        result [ 0 ] = temp_leftPullout;
        return true;
    } else if ( code.compare("temprightpullout") == 0 || code.compare("rightpullout") == 0 ) {
        result.resize(1);
        result [ 0 ] = temp_rightPullout;
        return true;
    } else if ( code.compare("tempbridgingforce") == 0 || code.compare("bridgingforce") == 0 ) {
        result.resize(1);
        result [ 0 ] = temp_bridgingForce;
        return true;
    } else if ( code.compare("tempfiberforce") == 0 || code.compare("fiberforce") == 0 ) {
        result.resize(1);
        result [ 0 ] = temp_fiberForce;
        return true;
    } else if ( code.compare("tempCrackOpening") == 0 || code.compare("crack_opening") == 0 ) {
        result.resize(1);
        result [ 0 ] = temp_crackOpening;
        return true;
    } else {
        return TensMechMaterialStatus :: giveValues(code, result);
    }
}

//////////////////////////////////////////////////////////
void FiberMaterialStatus :: init() {
    FiberMaterial *fibMaterial = static_cast< FiberMaterial * >( mat );
    double Ef = fibMaterial->giveEf();
    double tau0 = fibMaterial->giveTau0();
    double Gd = fibMaterial->giveGd();
    double betaf = fibMaterial->giveBetaf();

    Fiber *fibElement = static_cast< Fiber * >( element );
    df = fibElement->giveDiameter();
    rightLe = fibElement->giveRightLength(idx);
    leftLe = fibElement->giveLeftLength(idx);
    fiberLength = fibElement->giveLength();
    fiberNormal = fibElement->giveDirVector();

    RigidBodyContact *mechElement = static_cast< RigidBodyContact * >( fibElement->giveRBContact(idx) );
    fiberNormalLocal = mechElement->transformToLocal(fiberNormal);
    contactNormal = mechElement->giveNormal();
    contactLength = mechElement->giveLength();

    right_F0 = M_PI * df * tau0 * rightLe;
    left_F0 = M_PI * df * tau0 * leftLe;
    criticalRightPullout = 2. * tau0 * pow(rightLe, 2) / ( Ef * df ) + sqrt( 8. * Gd * pow(rightLe, 2) / ( Ef * df ) );
    criticalLeftPullout = 2. * tau0 * pow(leftLe, 2) / ( Ef * df ) + sqrt( 8. * Gd * pow(leftLe, 2) / ( Ef * df ) );

    if ( abs(rightLe - leftLe) > 1e-12 and rightLe > leftLe ) {
        Le1 = leftLe;
        F01 = left_F0;
        vd1 = criticalLeftPullout;
        Le2 = rightLe;
        F02 = right_F0;
        vd2 = criticalRightPullout;
        limitPullout = Le1 + vd1; // at least one side of the fiber is fully pulled out of matrix
    } else {
        Le1 = rightLe;
        F01 = right_F0;
        vd1 = criticalRightPullout;
        Le2 = leftLe;
        F02 = left_F0;
        vd2 = criticalLeftPullout;
        limitPullout = Le1 + vd1;
    }

    if ( betaf > 0 ) {
        A1 = -betaf;
        B1 = Le1 * betaf - df + 2. * betaf * vd1;
        C1 = Le1 * df - Le1 * betaf * vd1 + vd1 * df - betaf * pow(vd1, 2);
        v1_Fmax = -B1 / ( 2. * A1 );
        Fmax = bridgingForcePullingOut(v1_Fmax, vd1, Le1, F01, df, betaf);

        A2 = -betaf;
        B2 = Le2 * betaf - df + 2. * betaf * vd2;
        C2 = Le2 * df - Le2 * betaf * vd2 + vd2 * df - betaf * pow(vd2, 2) - Fmax * Le2 * df / F02;
        v2_Fmax1 = ( -B2 + sqrt(pow(B2, 2) - 4. * A2 * C2) ) / ( 2. * A2 );
        v2_Fmax2 = ( -B2 - sqrt(pow(B2, 2) - 4. * A2 * C2) ) / ( 2. * A2 );
        if ( v2_Fmax1 < v2_Fmax2 ) {
            v2_pullingOutFmax = v2_Fmax1;
        } else {
            v2_pullingOutFmax = v2_Fmax2;
        }

        v2_debondingFmax = pow(Fmax, 2) * 2. / ( pow(M_PI, 2) * Ef * pow(df, 3) * ( tau0 + Gd / vd2 ) );
    }

    inclineAngle = acos( ( fiberNormal ).dot(contactNormal) );
    if ( inclineAngle > M_PI / 2. ) {
        inclineAngle = M_PI - inclineAngle;
    }
}

//////////////////////////////////////////////////////////
void FiberMaterialStatus :: update() {
    crackOpening = temp_crackOpening;
    if ( crackOpening > maxCrackOpening ) {
        maxCrackOpening = crackOpening;
    } else {
        maxCrackOpening = maxCrackOpening;
    }

    rightPullout = temp_rightPullout;
    leftPullout = temp_leftPullout;

    bridgingForce = temp_bridgingForce;
    fiberForce = temp_fiberForce;

    debondedFiber = temp_debondedFiber;
    pulledOutFiber = temp_pulledOutFiber;
    rupturedFiber = temp_rupturedFiber;
    closingCrack = temp_closingCrack;

    MaterialStatus :: update();
}

//////////////////////////////////////////////////////////
void FiberMaterialStatus :: resetTemporaryVariables() {
    MaterialStatus :: resetTemporaryVariables();
}

//////////////////////////////////////////////////////////
// ------- EQUATIONS FOR FIBER CONSTITUTIVE LAW ------- //
//////////////////////////////////////////////////////////
//------------------------------------------------------//
//                   DEBONDING STAGE                    //
//------------------------------------------------------//
// -------- ORIGINAL EQUATIONS FROM LITERATURE -------- // with zero crack, there is non-zero bridging force in fiber due to the fracture energy Gd
//double bridgingForceDebonding( double v, double vd, double df, double Ef, double tau0, double Gd ) {
//    return sqrt( 0.5 * pow(M_PI,2) * Ef * pow(df,3) * ( tau0 * v + Gd ) );
//}
//double derivBridgingForceDebonding( double v, double vd, double df, double Ef, double tau0, double Gd ) {
//    return 0.5 / sqrt( 0.5 * pow(M_PI,2) * Ef * pow(df,3) * ( tau0 * v + Gd ) ) * ( 0.5 * pow(M_PI,2) * Ef * pow(df,3) * tau0 );
//}
// ---------------- MODIFIED EQUATIONS ---------------- // linear distribution of Gd
double bridgingForceDebonding(double v, double vd, double df, double Ef, double tau0, double Gd) {
    return sqrt( 0.5 * pow(M_PI, 2) * Ef * pow(df, 3) * ( tau0 * v + Gd * v / vd ) );
}
double derivBridgingForceDebonding(double v, double vd, double df, double Ef, double tau0, double Gd) {
    return 0.5 / sqrt( 0.5 * pow(M_PI, 2) * Ef * pow(df, 3) * ( tau0 * v + Gd * v / vd ) ) * ( 0.5 * pow(M_PI, 2) * Ef * pow(df, 3) * ( tau0 + Gd / vd ) );
}
//------------------------------------------------------//
//                  PULLING-OUT STAGE                   //
//------------------------------------------------------//
double bridgingForcePullingOut(double v, double vd, double Le, double F0, double df, double betaf) {
    return F0 * ( 1 - ( v - vd ) / Le ) * ( 1 +  betaf * ( v - vd ) / df );
}
double derivBridgingForcePullingOut(double v, double vd, double Le, double F0, double df, double betaf) {
    return F0 / ( Le * df ) * ( Le * betaf - df - 2. * v * betaf + 2. * vd * betaf );
}
//------------------------------------------------------//
//                   UNLOADING STAGE                    //
//------------------------------------------------------//
double bridgingForceUnloading(double deltaX, double deltaY, double x) {   // linear curve to zero
    return deltaY / deltaX * x;
}
double derivBridgingForceUnloading(double deltaX, double deltaY) {
    return deltaY / deltaX;
}

//////////////////////////////////////////////////////////
Matrix FiberMaterialStatus :: giveStiffnessTensor(string type) const {
    ( void ) type;

    Matrix stiffness = Matrix :: Zero(3, 3);

    if ( temp_fiberForce > 1e-12 ) {
        Vector alpha1 = temp_stress / ( temp_leftPullout / leftLe );

        Vector alpha2 = temp_stress / ( temp_rightPullout / rightLe );

        stiffness(0, 0) = alpha1 [ 0 ] * alpha2 [ 0 ] / ( alpha1 [ 0 ] + alpha2 [ 0 ] );
        stiffness(1, 1) = alpha1 [ 1 ] * alpha2 [ 1 ] / ( alpha1 [ 1 ] + alpha2 [ 1 ] );
        stiffness(2, 2) = alpha1 [ 2 ] * alpha2 [ 2 ] / ( alpha1 [ 2 ] + alpha2 [ 2 ] );
    }

    //cout << "FiberMaterialStatus::giveStiffnessTensor" << endl;
    cout.flush();
    return stiffness;
}

//////////////////////////////////////////////////////////
void FiberMaterialStatus :: computeStress( double timeStep) {
    ( void ) timeStep;

    FiberMaterial *fibMaterial = static_cast< FiberMaterial * >( mat );
    double Ef = fibMaterial->giveEf();
    double Gd = fibMaterial->giveGd();
    double tau0 = fibMaterial->giveTau0();
    double betaf = fibMaterial->giveBetaf();
    double ft = fibMaterial->giveFt();
    double Ksn = fibMaterial->giveKsn();   // snubbing coefficient
    double Ksp = fibMaterial->giveKsp();   // spalling coefficient
    double Krup = fibMaterial->giveKrup();   // rupture coefficient

    Fiber *fibElement = static_cast< Fiber * >( element );
    CSLMaterialStatus *CSLMatStatus = static_cast< CSLMaterialStatus * >( fibElement->giveRBContact(idx)->giveMatStatus(0) );

    crackOpeningVector = CSLMatStatus->giveCrackOpeningVector();
    if ( crackOpeningVector [ 0 ] < 0 ) {
        crackOpeningVector [ 0 ] = 0;
    }
    temp_crackOpening = crackOpeningVector.norm();
    incrementOfCrack = temp_crackOpening - crackOpening;

    int iteration = 0;
    int limIteration = 150;

    // --------------- PULLOUTS OF THE FIBER -----------------
    temp_rightPullout = rightPullout + incrementOfCrack / 2.;
    temp_leftPullout = leftPullout + incrementOfCrack / 2.;

    // ----------------- BRIDGING FORCE ---------------------- NOTE: considering fiber perpendicular to the contact/crack face
    if ( pulledOutFiber == 1 or rupturedFiber == 1 ) { // fiber was pulled out of matrix or ruptured in the previous step
        temp_bridgingForce = 0;
        temp_rightPullout = rightPullout;
        temp_leftPullout = leftPullout;
        if ( pulledOutFiber == 1 ) {
            temp_pulledOutFiber = 1;
        } else if ( rupturedFiber == 1 ) {
            temp_rupturedFiber = 1;
        }
    } else if ( temp_crackOpening <= 1e-12 ) { // fiber does not bridge the crack
        temp_rightPullout = temp_leftPullout = 0;
        temp_bridgingForce = 0;
    } else if ( incrementOfCrack > 0 ) { // crack is opening
        if ( closingCrack == 1 and maxCrackOpening >= temp_crackOpening ) { // crack was closing in the past and now is opening again, but the crack does not reach the max. opening from the past
            while ( iteration <= limIteration ) { // Newton iteration v_new = v_old - f(v_old)/f'(v_old)
                rightForce = bridgingForceUnloading(rightPullout, bridgingForce, temp_rightPullout);
                leftForce = bridgingForceUnloading(leftPullout, bridgingForce, temp_leftPullout);
                if ( abs(rightForce - leftForce) < 1e-9 ) {
                    temp_bridgingForce = rightForce;
                    break;
                } else {
                    temp_leftPullout = temp_leftPullout - ( bridgingForceUnloading(leftPullout, bridgingForce, temp_leftPullout) - bridgingForceUnloading(rightPullout, bridgingForce, temp_crackOpening - temp_leftPullout) ) / ( derivBridgingForceUnloading(leftPullout, bridgingForce) + derivBridgingForceUnloading(rightPullout, bridgingForce) );
                    temp_rightPullout = temp_crackOpening - temp_leftPullout;
                }
                iteration += 1;
            }
            if ( iteration > limIteration ) {
                cout << " FIBER ERROR (1): Force equilibrium was not found ! Limit iteration was reached ! " << endl;
                exit(1);
            }
        } else if ( temp_rightPullout <= criticalRightPullout and temp_leftPullout <= criticalLeftPullout ) { // both sides of the fiber are in debonding stage
            temp_bridgingForce = bridgingForceDebonding(temp_rightPullout, criticalRightPullout, df, Ef, tau0, Gd);
        } else if ( abs(rightLe - leftLe) <= 1e-12 ) { // both sides of the fiber are equal and fully debonded
            temp_bridgingForce = bridgingForcePullingOut(temp_rightPullout, criticalRightPullout, rightLe, right_F0, df, betaf);
            if ( temp_rightPullout >= limitPullout or temp_bridgingForce <= 0 ) { // check if the fiber is not fully pulled out of the matrix
                temp_bridgingForce = 0;
                temp_pulledOutFiber = 1;
                if ( temp_rightPullout < limitPullout ) { // for case with betaf << 0
                    temp_rightPullout = temp_crackOpening / 2.;
                    temp_leftPullout = temp_crackOpening / 2.;
                } else {
                    temp_rightPullout = temp_leftPullout = limitPullout;
                }
            }
        } else if ( rightLe != leftLe ) { // shorter side of the fiber (1) is pulling-out of the matrix, longer (2) is unloading
            if ( rightLe < leftLe ) {
                v1 = temp_rightPullout;
                v2 = temp_leftPullout;
            } else if ( rightLe > leftLe ) {
                v1 = temp_leftPullout;
                v2 = temp_rightPullout;
            }
            if ( debondedFiber == 0 ) {
                while ( iteration <= limIteration ) { // Newton iteration v_new = v_old - f(v_old)/f'(v_old)
                    bridgingForce1 = bridgingForcePullingOut(v1, vd1, Le1, F01, df, betaf);
                    bridgingForce2 = bridgingForceUnloading(vd1, bridgingForceDebonding(vd1, vd1, df, Ef, tau0, Gd), v2);
                    if ( v1 >= limitPullout or bridgingForce1 <= 0 ) { // check if the shorter side of the fiber is not fully pulled out of the matrix
                        temp_bridgingForce = 0;
                        temp_pulledOutFiber = 1;
                        v2 = 0;
                        if ( v1 < limitPullout ) { // for case with betaf << 0
                            v1 = temp_crackOpening;
                        } else {
                            v1 = limitPullout;
                        }
                        break;
                    } else if ( abs(bridgingForce1 - bridgingForce2) < 1e-9 ) {
                        temp_bridgingForce = bridgingForce1;
                        if ( temp_bridgingForce > bridgingForceDebonding(vd1, vd1, df, Ef, tau0, Gd) ) {   // longer side of the fiber (2) is in debonding stage
                            temp_debondedFiber = 1;
                        }
                        break;
                    } else {
                        v1 = v1 - ( bridgingForcePullingOut(v1, vd1, Le1, F01, df, betaf) - bridgingForceUnloading(vd1, bridgingForceDebonding(vd1, vd1, df, Ef, tau0, Gd), temp_crackOpening - v1) ) / ( derivBridgingForcePullingOut(v1, vd1, Le1, F01, df, betaf) + derivBridgingForceUnloading( vd1, bridgingForceDebonding(vd1, vd1, df, Ef, tau0, Gd) ) );
                        v2 = temp_crackOpening - v1;
                    }
                    iteration += 1;
                }
            }
            if ( debondedFiber == 1 or temp_debondedFiber == 1 ) {
                while ( iteration <= limIteration ) { // Newton iteration v_new = v_old - f(v_old)/f'(v_old)
                    bridgingForce1 = bridgingForcePullingOut(v1, vd1, Le1, F01, df, betaf);
                    bridgingForce2 = bridgingForceDebonding(v2, vd2, df, Ef, tau0, Gd);
                    if ( abs(bridgingForce1 - bridgingForce2) < 1e-9 ) {
                        temp_bridgingForce = bridgingForce1;
                        temp_debondedFiber = 1;
                        if ( temp_bridgingForce > bridgingForceDebonding(vd2, vd2, df, Ef, tau0, Gd) ) {   // longer side of the fiber (2) is in pulling-out stage
                            temp_debondedFiber = 2;
                        } else if ( v1 > v1_Fmax ) {
                            temp_debondedFiber = 3;
                        }
                        break;
                    } else {
                        v1 = v1 - ( bridgingForcePullingOut(v1, vd1, Le1, F01, df, betaf) - bridgingForceDebonding(temp_crackOpening - v1, vd2, df, Ef, tau0, Gd) ) / ( derivBridgingForcePullingOut(v1, vd1, Le1, F01, df, betaf) + derivBridgingForceDebonding(temp_crackOpening - v1, vd2, df, Ef, tau0, Gd) );
                        v2 = temp_crackOpening - v1;
                    }
                    iteration += 1;
                }
            }
            if ( debondedFiber == 2 or temp_debondedFiber == 2 ) {
                while ( iteration <= limIteration ) { // Newton iteration v_new = v_old - f(v_old)/f'(v_old)
                    bridgingForce1 = bridgingForcePullingOut(v1, vd1, Le1, F01, df, betaf);
                    bridgingForce2 = bridgingForcePullingOut(v2, vd2, Le2, F02, df, betaf);
                    if ( abs(bridgingForce1 - bridgingForce2) < 1e-9 ) {
                        temp_bridgingForce = bridgingForce1;
                        temp_debondedFiber = 2;
                        if ( v1 > v1_Fmax ) {
                            temp_debondedFiber = 4;
                        }
                        break;
                    } else {
                        v1 = v1 - ( bridgingForcePullingOut(v1, vd1, Le1, F01, df, betaf) - bridgingForcePullingOut(temp_crackOpening - v1, vd2, Le2, F02, df, betaf) ) / ( derivBridgingForcePullingOut(v1, vd1, Le1, F01, df, betaf) + derivBridgingForcePullingOut(temp_crackOpening - v1, vd2, Le2, F02, df, betaf) );
                        v2 = temp_crackOpening - v1;
                    }
                    iteration += 1;
                }
            }
            if ( debondedFiber == 3 or temp_debondedFiber == 3 ) {
                while ( iteration <= limIteration ) { // Newton iteration v_new = v_old - f(v_old)/f'(v_old)
                    bridgingForce1 = bridgingForcePullingOut(v1, vd1, Le1, F01, df, betaf);
                    bridgingForce2 = bridgingForceUnloading(v2_debondingFmax, Fmax, v2);
                    if ( bridgingForce1 <= 0 ) { // check if the shorter side of the fiber is not fully pulled out of the matrix
                        temp_bridgingForce = 0;
                        v1 = limitPullout;
                        v2 = 0;
                        temp_pulledOutFiber = 1;
                        break;
                    } else if ( abs(bridgingForce1 - bridgingForce2) < 1e-9 ) {
                        temp_bridgingForce = bridgingForce1;
                        temp_debondedFiber = 3;
                        break;
                    } else {
                        v1 = v1 - ( bridgingForcePullingOut(v1, vd1, Le1, F01, df, betaf) - bridgingForceUnloading(v2_debondingFmax, Fmax, temp_crackOpening - v1) ) / ( derivBridgingForcePullingOut(v1, vd1, Le1, F01, df, betaf) + derivBridgingForceUnloading(v2_debondingFmax, Fmax) );
                        v2 = temp_crackOpening - v1;
                    }
                    iteration += 1;
                }
            }
            if ( debondedFiber == 4 or temp_debondedFiber == 4 ) {
                while ( iteration <= limIteration ) { // Newton iteration v_new = v_old - f(v_old)/f'(v_old)
                    bridgingForce1 = bridgingForcePullingOut(v1, vd1, Le1, F01, df, betaf);
                    bridgingForce2 = bridgingForceUnloading(v2_pullingOutFmax, Fmax, v2);
                    if ( v1 >= limitPullout or bridgingForce1 <= 0 ) { // check if the shorter side of the fiber is not fully pulled out of the matrix
                        temp_bridgingForce = 0;
                        v1 = limitPullout;
                        v2 = 0;
                        temp_pulledOutFiber = 1;
                        break;
                    } else if ( abs(bridgingForce1 - bridgingForce2) < 1e-9 ) {
                        temp_bridgingForce = bridgingForce1;
                        temp_debondedFiber = 4;
                        break;
                    } else {
                        v1 = v1 - ( bridgingForcePullingOut(v1, vd1, Le1, F01, df, betaf) - bridgingForceUnloading(v2_pullingOutFmax, Fmax, temp_crackOpening - v1) ) / ( derivBridgingForcePullingOut(v1, vd1, Le1, F01, df, betaf) + derivBridgingForceUnloading(v2_pullingOutFmax, Fmax) );
                        v2 = temp_crackOpening - v1;
                    }
                    iteration += 1;
                }
            }
            if ( iteration > limIteration ) {
                cout << " FIBER ERROR (2): Force equilibrium was not found ! Limit iteration was reached ! " << endl;
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
        while ( iteration <= limIteration ) { // Newton iteration v_new = v_old - f(v_old)/f'(v_old)
            rightForce = bridgingForceUnloading(rightPullout, bridgingForce, temp_rightPullout);
            leftForce = bridgingForceUnloading(leftPullout, bridgingForce, temp_leftPullout);
            if ( abs(rightForce - leftForce) < 1e-9 ) {
                temp_bridgingForce = rightForce;
                if ( temp_bridgingForce < 0 ) { // check if the crack is not fully closed
                    temp_bridgingForce = 0;
                    temp_rightPullout = temp_leftPullout = 0;
                    temp_crackOpening = 0;
                    cout << " FIBER ERROR: Crack should not be fully closed ! " << endl;
                    exit(1);
                }
                break;
            } else {
                temp_leftPullout = temp_leftPullout - ( bridgingForceUnloading(leftPullout, bridgingForce, temp_leftPullout) - bridgingForceUnloading(rightPullout, bridgingForce, temp_crackOpening - temp_leftPullout) ) / ( derivBridgingForceUnloading(leftPullout, bridgingForce) + derivBridgingForceUnloading(rightPullout, bridgingForce) );
                temp_rightPullout = temp_crackOpening - temp_leftPullout;
            }
            iteration += 1;
        }
        if ( iteration > limIteration ) {
            cout << " FIBER ERROR (3): Force equilibrium was not found ! Limit iteration was reached ! " << endl;
            exit(1);
        }
    }

    // ------------------- STRESS VECTOR --------------------- NOTE: fiberForce is modified bridgingForce due to incline of the fiber
    if ( temp_bridgingForce <= 1e-12 ) {
        temp_fiberForce = 0;
        temp_stress = Vector :: Zero(temp_strain.size() );
    } else if ( inclineAngle <= 1e-12 or crackOpeningVector [ 0 ] <= 1e-12 ) { // without SNUBBING MODEL - no microeffect at exit points
        temp_fiberForce = temp_bridgingForce;
        temp_stress = temp_fiberForce * crackOpeningVector / temp_crackOpening;
    } else if ( inclineAngle > 1e-12 ) { // with SNUBBING MODEL - considering microeffect at exit points
        iteration = 0;
        if ( fiberNormalLocal [ 0 ] < 0 ) {
            fiberNormalLocal = fiberNormalLocal * ( -1. );
        }
        while ( iteration <= limIteration ) {
            w = crackOpeningVector + 2. * spallingLength * fiberNormalLocal;
            mf = w / w.norm();
            deflectionAngle = acos( fiberNormalLocal.dot(mf) );
            temp_fiberForce = temp_bridgingForce * exp(Ksn * deflectionAngle);
            temp_stress = temp_fiberForce * mf;
            temp_spallingLength = ( temp_stress(0) * sin(inclineAngle / 2.) ) / ( Ksp * ft * df * pow(cos(inclineAngle / 2.), 2) );
            if ( ( abs(temp_spallingLength - spallingLength) ) / temp_spallingLength < 1e-6 ) { // relative error of the spalling length
                break;
            } else {
                spallingLength = temp_spallingLength;
            }
            iteration += 1;
        }
        if ( iteration > limIteration ) {
            cout << " FIBER ERROR (4): Limit iteration for snubbing model was reached ! " << endl;
            exit(1);
        }
    }

    // -------------------------- RUPTURE CONDITION ---------------------------------------------------
    double temp_fiberStress = 4. * temp_fiberForce / M_PI / pow(df, 2);
    double limitFiberStress = ft * exp(-Krup * deflectionAngle);
    if ( temp_fiberStress > limitFiberStress ) {
        temp_bridgingForce = temp_fiberForce = 0;
        temp_stress = Vector :: Zero(temp_strain.size() );
        temp_rupturedFiber = 1;
    }
}

//////////////////////////////////////////////////////////
void FiberMaterialStatus :: computeStressWithFrozenIntVars( double timeStep) {
    ( void ) timeStep;

    //temp_crackOpening = strain.norm() * contactLength;
    //Vector stressWithFrozenIntVars = giveStiffnessTensor("elastic") * temp_crackOpening;

    Fiber *fibElement = static_cast< Fiber * >( element );
    CSLMaterialStatus *CSLMatStatus = static_cast< CSLMaterialStatus * >( fibElement->giveRBContact(idx)->giveMatStatus(0) );

    crackOpeningVector = CSLMatStatus->giveCrackOpeningVector();
    if ( crackOpeningVector [ 0 ] < 0 ) {
        crackOpeningVector [ 0 ] = 0;
    }
    temp_stress = giveStiffnessTensor("elastic") * crackOpeningVector;
}

//////////////////////////////////////////////////////////
void FiberMaterialStatus :: setParameterValue(string code, double value) {
    if ( code.compare("crackOpening") == 0 ) {
        crackOpening = value;
    } else {
        MaterialStatus :: setParameterValue(code, value);
    }
}

//////////////////////////////////////////////////////////
// FIBER MATERIAL

FiberMaterial :: FiberMaterial(unsigned dimension) : TensMechMaterial(dimension) {
    name = "fiber material";
    strainsize = dimension;

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
