#include "material_fiber.h"
#include "element_fiber.h"
#include "element_discrete.h"

using namespace std;

//////////////////////////////////////////////////////////
// FIBER MATERIAL STATUS

FiberMaterialStatus :: FiberMaterialStatus(FiberMaterial *m, Element *e, unsigned ipnum) : MaterialStatus(m, e, ipnum) {
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
void FiberMaterialStatus :: giveValues(string code, Vector &result) const {
    MaterialStatus :: giveValues(code, result);
}

//////////////////////////////////////////////////////////
void FiberMaterialStatus :: init() { 

    RigidBodyContact* mechElement = static_cast < RigidBodyContact* > ( element );
    contactLength = mechElement -> giveLength();

    FiberMaterial* fibermat = static_cast < FiberMaterial* > ( mat );    
    double Ef = fibermat -> giveEf();    
    double tau0 = fibermat -> giveTau0();
    double Gd = fibermat -> giveGd();
    
    Fiber* fiberelem = static_cast < Fiber* > ( element );
    df = fiberelem -> giveDiameter();
    rightLe = 0.006;
    //fiberelem -> giveRightLength( idx );
    leftLe = 0.004;
    //fiberelem -> giveLeftLength( idx ); 
    
    right_F0 = M_PI * df * tau0 * rightLe;
    left_F0 = M_PI * df * tau0 * leftLe;
    limit_rightPullout = 2. * tau0 * pow(rightLe,2) / ( Ef * df ) + sqrt( 8. * Gd * pow(rightLe,2) / ( Ef * df ) );
    limit_leftPullout = 2. * tau0 * pow(leftLe,2) / ( Ef * df ) + sqrt( 8. * Gd * pow(leftLe,2) / ( Ef * df ) );
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
double bridgingForce_bonded( double v, double df, double Ef, double tau0, double Gd ) {
    return sqrt( 0.5 * pow(M_PI,2) * Ef * pow(df,3) * ( tau0 * v + Gd ) );
}
double derivF_bonded( double v, double df, double Ef, double tau0, double Gd ) {
    return 0.5 / sqrt( 0.5 * pow(M_PI,2) * Ef * pow(df,3) * ( tau0 * v + Gd ) ) * ( 0.5 * pow(M_PI,2) * Ef * pow(df,3) * tau0 );
}
double bridgingForce_debonded( double v, double vd, double Le, double F0, double df, double betaf ) {
    return F0 * ( 1 - ( v - vd ) / Le ) * ( 1 +  betaf * ( v - vd ) / df );
}
double derivF_debonded( double v, double vd, double Le, double F0, double df, double betaf ) {
    return F0 / ( Le * df ) * ( Le * betaf - df - 2. * v * betaf + 2. * vd * betaf ); 
}
double bridgingForce_unloading( double deltaX, double deltaY, double x ) { // linear curve to zero 
    return deltaY / deltaX * x; 
}
double derivF_unloading( double deltaX, double deltaY ) {
    return deltaY / deltaX; 
}

//////////////////////////////////////////////////////////
Matrix FiberMaterialStatus :: giveStiffnessTensor(string type, unsigned dim) const {
    //TODO
    ( void ) type;
    ( void ) dim;
    cout << "FiberMaterialStatus::giveStiffnessTensor" << endl;
    cout.flush();

    return Matrix :: Zero(dim, dim);
}

//////////////////////////////////////////////////////////
Vector FiberMaterialStatus :: giveStress(const Vector &strain, double timeStep) {
    //TODO
    ( void ) strain;
    ( void ) timeStep;
        
    FiberMaterial* fibermat = static_cast < FiberMaterial* > ( mat );    
    double Ef = fibermat -> giveEf();    
    double Gd = fibermat -> giveGd();
    double tau0 = fibermat -> giveTau0();
    double betaf = fibermat -> giveBetaf();
// ---------------------------------------------------------------------     
    temp_crack_opening = 2. * 0.0010027;
    //strain[0] * contactLength;
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
            temp_bridgingForce = bridgingForce_bonded( temp_rightPullout, df, Ef, tau0, Gd );
        } else if ( rightLe == leftLe ) {
            temp_bridgingForce = bridgingForce_debonded( temp_rightPullout, limit_rightPullout, rightLe, right_F0, df, betaf );
            if ( temp_bridgingForce <=0 ) {
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
                bridgingForce1 = bridgingForce_debonded( v1, vd1, Le1, F01, df, betaf ); 
                bridgingForce2 = bridgingForce_unloading( vd1, bridgingForce_bonded( vd1, df, Ef, tau0, Gd ), v2 );
                if ( abs(bridgingForce1 - bridgingForce2) > tolerance ) {
                    v1 = v1 - ( bridgingForce_debonded( v1, vd1, Le1, F01, df, betaf ) - bridgingForce_unloading( vd1, bridgingForce_bonded( vd1, df, Ef, tau0, Gd ), temp_crack_opening - v1 ) ) / ( derivF_debonded( v1, vd1, Le1, F01, df, betaf ) + derivF_unloading( vd1, bridgingForce_bonded( vd1, df, Ef, tau0, Gd ) ) ); 
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
                                bridgingForce1 = bridgingForce_debonded( v1, vd1, Le1, F01, df, betaf ); 
                                bridgingForce2 = bridgingForce_bonded( v2, df, Ef, tau0, Gd ); 
                                if ( abs(bridgingForce1 - bridgingForce2) > tolerance ) {
                                    v1 = v1 - ( bridgingForce_debonded( v1, vd1, Le1, F01, df, betaf ) - bridgingForce_bonded( temp_crack_opening - v1, df, Ef, tau0, Gd ) ) / ( derivF_debonded( v1, vd1, Le1, F01, df, betaf ) + derivF_bonded( temp_crack_opening - v1, df, Ef, tau0, Gd ) );
                                    v2 = temp_crack_opening - v1;
                                } else if ( abs(bridgingForce1 - bridgingForce2) <= tolerance ) {
                                    temp_bridgingForce = bridgingForce1;
                                    if ( v2 < vd1 ) {
                                        cout << " ERROR pullout v2! We need v2 > vd1.  " << endl;
                                    } else if ( v2 > vd2 ) {
                                    cout << " v2 > vd2 " << endl; 
                                        while ( iteration <= limit_iteration ) {
                                        // Newton iteration v1 = v0 - f(v0)/f'(v0) 
                                            bridgingForce1 = bridgingForce_debonded( v1, vd1, Le1, F01, df, betaf ); 
                                            bridgingForce2 = bridgingForce_debonded( v2, vd2, Le2, F02, df, betaf ); 
                                            if ( abs(bridgingForce1 - bridgingForce2) > tolerance ) {
                                                v1 = v1 - ( bridgingForce_debonded( v1, vd1, Le1, F01, df, betaf ) - bridgingForce_debonded( temp_crack_opening - v1, vd2, Le2, F02, df, betaf ) ) / ( derivF_debonded( v1, vd1, Le1, F01, df, betaf ) + derivF_debonded( temp_crack_opening - v1, vd1, Le1, F01, df, betaf ) );
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
                        } else {
                        }                      
                    }
                    break; 
                }
            iteration += 1; 
            }
            if ( abs(bridgingForce1 - bridgingForce2) > tolerance ) {
                cout << " WARNING: Limit iteration was reached ! " << endl;  
            }
        }
        if ( rightLe > leftLe ) {
                temp_leftPullout = v1; 
                temp_rightPullout = v2; 
            } else if ( rightLe < leftLe ) {
                temp_rightPullout = v1; 
                temp_leftPullout = v2; 
            }
    } else if ( incrementOfCrack <= 0 ) {
        temp_rightPullout = right_pullout + incrementOfCrack / 2.; 
        temp_leftPullout = left_pullout + incrementOfCrack / 2.;
        while ( iteration <= limit_iteration ) {
        // Newton iteration v1 = v0 - f(v0)/f'(v0) 
            rightForce = bridgingForce_unloading( right_pullout, bridgingForce, temp_rightPullout ); 
            leftForce = bridgingForce_unloading( left_pullout, bridgingForce, temp_leftPullout ); 
            if ( abs(rightForce - leftForce) > tolerance ) {
                temp_leftPullout = temp_leftPullout - ( bridgingForce_unloading( left_pullout, bridgingForce, temp_leftPullout ) - bridgingForce_unloading( right_pullout, bridgingForce, temp_crack_opening - temp_leftPullout ) ) / ( derivF_unloading( left_pullout, bridgingForce ) + derivF_unloading( right_pullout, bridgingForce ) ); 
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

    cout << "Forces " << temp_bridgingForce << endl;
    cout << temp_rightPullout << endl;
    cout << temp_leftPullout << endl;
    cout << "Soucet pulloutu " << temp_rightPullout + temp_leftPullout << endl;
    cout << "Crack width " << temp_crack_opening << endl;
    cout << "Iterations " << iteration << endl;

    
    cout << "FiberMaterialStatus::giveStress" << endl;
    cout.flush();
    return Vector :: Zero(strain.size() );
}

//////////////////////////////////////////////////////////
Vector FiberMaterialStatus :: giveStressWithFrozenIntVars(const Vector &strain, double timeStep) {
    //TODO
    ( void ) strain;
    ( void ) timeStep;
    cout << "FiberMaterialStatus::giveStressWithFrozenIntVars" << endl;
    cout.flush();
    return Vector :: Zero(strain.size() );
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

FiberMaterial :: FiberMaterial(){
    name = "fiber material";
    // if variables not specified on the input, use default multipliers
    Ef = Gd = tau0 = betaf = ft =0;
}

//////////////////////////////////////////////////////////
void FiberMaterial :: readFromLine(istringstream &iss) {
    Material :: readFromLine(iss); //read elastic parameters

    iss.clear(); // clear string stream
    iss.seekg(0, iss.beg); //reset position in string stream

    // initialize all values to zero (NOTE probably no ned in linux, but in windows necessary)
    Ef = Gd = tau0 = betaf = ft = 0;

    string param;
    bool bft, bGd, btau0, bbetaf, bEf;
    bft = bGd = btau0 = bbetaf = bEf =  false;

    while ( !iss.eof() ) {
        iss >> param;
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
    };
};

//////////////////////////////////////////////////////////
MaterialStatus *FiberMaterial :: giveNewMaterialStatus(Element *e, unsigned ipnum) {
    FiberMaterialStatus *newStatus = new FiberMaterialStatus(this, e, ipnum); //needs to be deleted manually
    return newStatus;
};

//////////////////////////////////////////////////////////
void FiberMaterial :: init() {};
