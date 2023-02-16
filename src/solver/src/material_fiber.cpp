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

    right_pullout = 0;
    temp_rightPullout = 0;
    left_pullout = 0;
    temp_leftPullout = 0;

    bridgingForce = 0;
    temp_bridgingForce = 0;
    rightForce = leftForce = 0;
}

//////////////////////////////////////////////////////////
bool FiberMaterialStatus :: giveValues(string code, Vector &result) const {
    return TensMechMaterialStatus :: giveValues(code, result);
}

//////////////////////////////////////////////////////////
void FiberMaterialStatus :: init() {
    Fiber *fib = static_cast< Fiber * >( element );
    RigidBodyContact *mechElement = static_cast< RigidBodyContact * >( fib->giveRBContact(idx) );
    contactLength = mechElement->giveLength();
    contactNormal = mechElement->giveNormal();

    temp_stress = Vector :: Zero(contactNormal.size() );
    crackOpeningVector = Vector :: Zero(contactNormal.size() );

    FiberMaterial *fibermat = static_cast< FiberMaterial * >( mat );
    double Ef = fibermat->giveEf();
    double tau0 = fibermat->giveTau0();
    double Gd = fibermat->giveGd();

    Fiber *fiberelem = static_cast< Fiber * >( element );
    df = fiberelem->giveDiameter();
    rightLe = fiberelem->giveRightLength(idx);
    leftLe = fiberelem->giveLeftLength(idx);
    fiberNormal = fiberelem->giveDirVector();

    inclineAngle = acos( ( fiberNormal ).dot(contactNormal) );

    right_F0 = M_PI * df * tau0 * rightLe;
    left_F0 = M_PI * df * tau0 * leftLe;
    limit_rightPullout = 2. * tau0 * pow(rightLe, 2) / ( Ef * df ) + sqrt( 8. * Gd * pow(rightLe, 2) / ( Ef * df ) );
    limit_leftPullout = 2. * tau0 * pow(leftLe, 2) / ( Ef * df ) + sqrt( 8. * Gd * pow(leftLe, 2) / ( Ef * df ) );
}

//////////////////////////////////////////////////////////
void FiberMaterialStatus :: update() {
    crack_opening = temp_crack_opening;
    right_pullout = temp_rightPullout;
    left_pullout = temp_leftPullout;
    bridgingForce = temp_bridgingForce;

    MaterialStatus :: update();
}

//////////////////////////////////////////////////////////
void FiberMaterialStatus :: resetTemporaryVariables() {
    MaterialStatus :: resetTemporaryVariables();
}

//////////////////////////////////////////////////////////
double bridgingForce_bonded(double v, double df, double Ef, double tau0, double Gd) {
    return sqrt( 0.5 * pow(M_PI, 2) * Ef * pow(df, 3) * ( tau0 * v + Gd ) );
}
double derivF_bonded(double v, double df, double Ef, double tau0, double Gd) {
    return 0.5 / sqrt( 0.5 * pow(M_PI, 2) * Ef * pow(df, 3) * ( tau0 * v + Gd ) ) * ( 0.5 * pow(M_PI, 2) * Ef * pow(df, 3) * tau0 );
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

    if ( temp_rightPullout > 0 and temp_leftPullout > 0 ) {
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

    FiberMaterial *fibermat = static_cast< FiberMaterial * >( mat );
    double Ef = fibermat->giveEf();
    double Gd = fibermat->giveGd();
    double tau0 = fibermat->giveTau0();
    double betaf = fibermat->giveBetaf();
    double ft = fibermat->giveFt();
    double Ksn = fibermat->giveKsn();   // snubbing coefficient
    double Ksp = fibermat->giveKsp();   // spalling coefficient
    double Krup = fibermat->giveKrup();   // rupture coefficient
    // ---------------------------------------------------------------------
    Fiber *fib = static_cast< Fiber * >( element );
    CSLMaterialStatus *status = static_cast< CSLMaterialStatus * >( fib->giveRBContact(idx)->giveMatStatus(0) );
    crackOpeningVector = status->giveCrackOpeningVector();
    temp_crack_opening = crackOpeningVector.norm();
    incrementOfCrack = temp_crack_opening - crack_opening;
    // ---------------------------------------------------------------------
    double v0; // pullout for the extreme of functions
    double Le1, F01, v1, vd1, bridgingForce1; // shorter fiber, will be debonded
    double Le2, F02, v2, vd2, bridgingForce2; // longer fiber, will be unloaded
    Le1 = F01 = v1 = vd1 = bridgingForce1 = Le2 = F02 = v2 = vd2 = bridgingForce2 = 0;
    // ---------------------------------------------------------------------
    int iteration = 0;
    double limit_iteration = 100;
    double tolerance = 1e-8;
    // ---------------------------------------------------------------------
    if ( incrementOfCrack > 0 ) {
        temp_rightPullout = right_pullout + incrementOfCrack / 2.;
        temp_leftPullout = left_pullout + incrementOfCrack / 2.;
        if ( temp_rightPullout <= limit_rightPullout && temp_leftPullout <= limit_leftPullout ) {
            temp_bridgingForce = bridgingForce_bonded(temp_rightPullout, df, Ef, tau0, Gd);
        } else if ( rightLe == leftLe ) {
            temp_bridgingForce = bridgingForce_debonded(temp_rightPullout, limit_rightPullout, rightLe, right_F0, df, betaf);
            if ( temp_bridgingForce <= 0 ) {
                temp_bridgingForce = 0;
            } else if ( betaf < 0 ) {
                v0 = ( rightLe * betaf - df + 2. * limit_rightPullout * betaf ) / ( 2. * betaf );
                if ( temp_rightPullout >= v0 ) {
                    temp_bridgingForce = 0;
                }
            }
        } else if ( rightLe != leftLe ) {
            if ( rightLe > leftLe ) {
                Le1 = leftLe;
                F01 = left_F0;
                v1 = temp_leftPullout;
                vd1 = limit_leftPullout;
                Le2 = rightLe;
                F02 = right_F0;
                v2 = temp_rightPullout;
                vd2 = limit_rightPullout;
            } else if ( rightLe < leftLe ) {
                Le1 = rightLe;
                F01 = right_F0;
                v1 = temp_rightPullout;
                vd1 = limit_rightPullout;
                Le2 = leftLe;
                F02 = left_F0;
                v2 = temp_leftPullout;
                vd2 = limit_leftPullout;
            }
            while ( iteration <= limit_iteration ) {
                // Newton iteration v1 = v0 - f(v0)/f'(v0)
                bridgingForce1 = bridgingForce_debonded(v1, vd1, Le1, F01, df, betaf);
                bridgingForce2 = bridgingForce_unloading(vd1, bridgingForce_bonded(vd1, df, Ef, tau0, Gd), v2);
                if ( abs(bridgingForce1 - bridgingForce2) > tolerance ) {
                    v1 = v1 - ( bridgingForce_debonded(v1, vd1, Le1, F01, df, betaf) - bridgingForce_unloading(vd1, bridgingForce_bonded(vd1, df, Ef, tau0, Gd), temp_crack_opening - v1) ) / ( derivF_debonded(v1, vd1, Le1, F01, df, betaf) + derivF_unloading( vd1, bridgingForce_bonded(vd1, df, Ef, tau0, Gd) ) );
                    v2 = temp_crack_opening - v1;
                } else if ( abs(bridgingForce1 - bridgingForce2) <= tolerance ) {
                    temp_bridgingForce = bridgingForce1;
                    if ( bridgingForce1 <= 0 or bridgingForce2 <= 0 ) {
                        temp_bridgingForce = 0;
                        v1 = temp_crack_opening;
                        v2 = 0;
                    } else if ( betaf < 0 ) {
                        v0 = ( rightLe * betaf - df + 2. * limit_rightPullout * betaf ) / ( 2. * betaf );
                        if ( v1 >= v0 ) {
                            temp_bridgingForce = 0;
                            v1 = temp_crack_opening;
                            v2 = 0;
                        }
                    } else if ( betaf > 0 ) {
                        if ( v2 > vd1 ) {
                            cout << " v2 > vd1 " << endl;
                            while ( iteration <= limit_iteration ) {
                                // Newton iteration v1 = v0 - f(v0)/f'(v0)
                                bridgingForce1 = bridgingForce_debonded(v1, vd1, Le1, F01, df, betaf);
                                bridgingForce2 = bridgingForce_bonded(v2, df, Ef, tau0, Gd);
                                if ( abs(bridgingForce1 - bridgingForce2) > tolerance ) {
                                    v1 = v1 - ( bridgingForce_debonded(v1, vd1, Le1, F01, df, betaf) - bridgingForce_bonded(temp_crack_opening - v1, df, Ef, tau0, Gd) ) / ( derivF_debonded(v1, vd1, Le1, F01, df, betaf) + derivF_bonded(temp_crack_opening - v1, df, Ef, tau0, Gd) );
                                    v2 = temp_crack_opening - v1;
                                } else if ( abs(bridgingForce1 - bridgingForce2) <= tolerance ) {
                                    temp_bridgingForce = bridgingForce1;
                                    if ( v2 < vd1 ) {
                                        cout << " ERROR pullout v2! We need v2 > vd1.  " << endl;
                                    } else if ( v2 > vd2 ) {
                                        cout << " v2 > vd2 " << endl;
                                        while ( iteration <= limit_iteration ) {
                                            // Newton iteration v1 = v0 - f(v0)/f'(v0)
                                            bridgingForce1 = bridgingForce_debonded(v1, vd1, Le1, F01, df, betaf);
                                            bridgingForce2 = bridgingForce_debonded(v2, vd2, Le2, F02, df, betaf);
                                            if ( abs(bridgingForce1 - bridgingForce2) > tolerance ) {
                                                v1 = v1 - ( bridgingForce_debonded(v1, vd1, Le1, F01, df, betaf) - bridgingForce_debonded(temp_crack_opening - v1, vd2, Le2, F02, df, betaf) ) / ( derivF_debonded(v1, vd1, Le1, F01, df, betaf) + derivF_debonded(temp_crack_opening - v1, vd1, Le1, F01, df, betaf) );
                                                v2 = temp_crack_opening - v1;
                                            } else if ( abs(bridgingForce1 - bridgingForce2) <= tolerance ) {
                                                temp_bridgingForce = bridgingForce1;
                                                cout << " Both sides of the fibre are debonding. " << endl;
                                                break;
                                            }
                                            iteration += 1;
                                        }
                                    }
                                    break;
                                }
                                iteration += 1;
                            }
                        } else {}
                    }
                    break;
                }
                iteration += 1;
            }
            if ( abs(bridgingForce1 - bridgingForce2) > tolerance ) {
                cout << " WARNING: Limit iteration was reached ! " << endl;
            }
            if ( rightLe > leftLe ) {
                temp_leftPullout = v1;
                temp_rightPullout = v2;
            } else if ( rightLe < leftLe ) {
                temp_rightPullout = v1;
                temp_leftPullout = v2;
            }
        }
    } else if ( incrementOfCrack <= 0 ) {
        temp_rightPullout = right_pullout + incrementOfCrack / 2.;
        temp_leftPullout = left_pullout + incrementOfCrack / 2.;
        while ( iteration <= limit_iteration ) {
            // Newton iteration v1 = v0 - f(v0)/f'(v0)
            rightForce = bridgingForce_unloading(right_pullout, bridgingForce, temp_rightPullout);
            leftForce = bridgingForce_unloading(left_pullout, bridgingForce, temp_leftPullout);
            if ( abs(rightForce - leftForce) > tolerance ) {
                temp_leftPullout = temp_leftPullout - ( bridgingForce_unloading(left_pullout, bridgingForce, temp_leftPullout) - bridgingForce_unloading(right_pullout, bridgingForce, temp_crack_opening - temp_leftPullout) ) / ( derivF_unloading(left_pullout, bridgingForce) + derivF_unloading(right_pullout, bridgingForce) );
                temp_rightPullout = temp_crack_opening - temp_leftPullout;
            } else if ( abs(rightForce - leftForce) <= tolerance ) {
                temp_bridgingForce = rightForce;
                break;
            }
            iteration += 1;
        }
        if ( rightForce <= 0 or leftForce <= 0 ) {
            temp_bridgingForce = 0;
            temp_rightPullout = 0;
            temp_leftPullout = 0;
        }
    }




    // microeffect of spalling - excel hotovy





    // RUPTURE CONDITION
    double deflectionAngle = 0; // budu znat z microeffectu
    double fiberStress = temp_bridgingForce;
    double limitStress = ft * exp(-Krup * deflectionAngle);
    if ( fiberStress > limitStress ) {
        temp_bridgingForce = 0;
    } else if ( fiberStress <= limitStress ) {
        temp_bridgingForce = temp_bridgingForce;
    }

    // FINAL STRESS VECTOR
    if ( temp_crack_opening == 0 ) {
        temp_stress = Vector :: Zero(strain.size() );
    } else {
        Vector w = crackOpeningVector / temp_crack_opening; // tohle vyjmout do microeffectu a pridat prispevok spalling length sf
        temp_stress = temp_bridgingForce * w;
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

    // initialize all values to zero (NOTE probably no ned in linux, but in windows necessary)
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
