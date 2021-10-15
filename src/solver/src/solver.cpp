#include "solver.h"
#include "adaptivity.h"
#define EPS2displ 1e-20
#define EPS2press 1e-12

//////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////
// BASIC SOLVER CLASS
Solver :: Solver() {
    name = "basic solver";
}


void Solver :: setContainers(ElementContainer *e, NodeContainer *n, FunctionContainer *functions) {
    elems = e;
    nodes = n;
    funcs = functions;
}

//////////////////////////////////////////////////////////
Solver *Solver :: readFromFile(const string filename) {
    string param, paramA, line;
    ifstream inputfile(filename.c_str() );
    if ( inputfile.is_open() ) {
        while ( getline(inputfile >> std :: ws, line) ) {
            if ( line.empty() ) {
                continue;
            }
            if ( line.at(0) == '#' ) {
                continue;
            }
            istringstream iss(line);
            iss >> param >> paramA;
            break;
        }
        inputfile.close();
    }
    if ( paramA.compare("Adaptive") == 0 || paramA.compare("adaptive") == 0 ) {
        if ( param.compare("SteadyStateNonLinearSolver") == 0 ) {
            Solver *newsolver = new AdaptiveSolver< SteadyStateNonLinearSolver >;
            newsolver->readFromFile(filename);
            cout << "Input file '" <<  filename << "' succesfully loaded; " << newsolver->name << " found" << endl;
            return newsolver;
        } else if ( param.compare("SteadyStateLinearSolver") == 0 ) {
            Solver *newsolver = new AdaptiveSolver< SteadyStateLinearSolver >;
            newsolver->readFromFile(filename);
            cout << "Input file '" <<  filename << "' succesfully loaded; " << newsolver->name << " found" << endl;
            return newsolver;
        } else {
            cerr << "Error: Solver " << param << " is not implemented for adaptivity" << endl;
            exit(EXIT_FAILURE);
        }
    } else {
        if ( param.compare("SteadyStateLinearSolver") == 0 ) {
            SteadyStateLinearSolver *newsolver = new SteadyStateLinearSolver();
            newsolver->readFromFile(filename);
            cout << "Input file '" <<  filename << "' succesfully loaded; " << newsolver->name << " found" << endl;
            return newsolver;
        } else if ( param.compare("SteadyStateNonLinearSolver") == 0 ) {
            SteadyStateNonLinearSolver *newsolver = new SteadyStateNonLinearSolver();
            newsolver->readFromFile(filename);
            cout << "Input file '" <<  filename << "' succesfully loaded; " << newsolver->name << " found" << endl;
            return newsolver;
        } else if ( param.compare("TransientLinearTransportSolver") == 0 ) {
            TransientLinearTransportSolver *newsolver = new TransientLinearTransportSolver();
            newsolver->readFromFile(filename);
            cout << "Input file '" <<  filename << "' succesfully loaded; " << newsolver->name << " found" << endl;
            return newsolver;
        } else if ( param.compare("TransientNonLinearTransportSolver") == 0 ) {
            TransientNonLinearTransportSolver *newsolver = new TransientNonLinearTransportSolver();
            newsolver->readFromFile(filename);
            cout << "Input file '" <<  filename << "' succesfully loaded; " << newsolver->name << " found" << endl;
            return newsolver;
        } else if ( param.compare("TransientLinearMechanicalSolver") == 0 ) {
            TransientLinearMechanicalSolver *newsolver = new TransientLinearMechanicalSolver();
            newsolver->readFromFile(filename);
            cout << "Input file '" <<  filename << "' succesfully loaded; " << newsolver->name << " found" << endl;
            return newsolver;
        } else if ( param.compare("TransientNonLinearMechanicalSolver") == 0 ) {
            TransientNonLinearMechanicalSolver *newsolver = new TransientNonLinearMechanicalSolver();
            newsolver->readFromFile(filename);
            cout << "Input file '" <<  filename << "' succesfully loaded; " << newsolver->name << " found" << endl;
            return newsolver;
        } else {
            cerr << "Error: Solver " << param << " is not implemented" << endl;
            exit(EXIT_FAILURE);
        }
    }
}

//////////////////////////////////////////////////////////
void Solver :: setNextStepTime() {
    double nextExtremeTime = funcs->giveTimeOfNextExtreme(time);
    // NOTE 1/4 of time step added to prevent next step extremely short
    if ( nextExtremeTime < time + 1.25 * dt ) {
        dt = nextExtremeTime - time;
        time = nextExtremeTime;
    } else {
        time += dt;
    }
}

//////////////////////////////////////////////////////////
void Solver :: runBeforeEachStep() {
    step += 1;
    setNextStepTime();

    load_old = load; //copy old load to be used in generalized alpha method
    load *= 0;  //clear nodal load
    ddr *= 0; //clear step contribution;
}

//////////////////////////////////////////////////////////
void Solver :: runAfterEachStep() {
    r = trial_r;
    f_int_old = f_int;
    f_ext_old = f_ext;

    if ( time + dt > termination_time ) {
        dt = termination_time - time;
    }
    if ( dt < 1e-15 ) {
        terminated = true;
    }
    elems->updateMaterialStatuses();
}

//////////////////////////////////////////////////////////
void Solver :: setTime(double t) {
    time = t;
    if ( time - 1e-15 > termination_time ) {
        time = termination_time;
        terminated = true;
    } else {
        terminated = false;
    }
    dt = initdt;
}

//////////////////////////////////////////////////////////
void Solver :: init(string init_r_file, string init_v_file, const bool initial) {
    ( void ) init_r_file;
    ( void ) init_v_file;

    elems->readMatStatsFromFile(this->init_time, this->init_step);
    if ( initial ) {
        step = init_step;
        time = init_time;
    }

    terminated = false;

    freeDoFnum = nodes->giveNumFreeDoFs();
    fixedDoFnum = ( nodes->giveTotalNumDoFs() - freeDoFnum );
    totalDoFnum = freeDoFnum + fixedDoFnum;

    load = Vector(totalDoFnum);
    r = Vector(totalDoFnum);
    f_ext = Vector(totalDoFnum);
    f_int = Vector(totalDoFnum);
    f_dam = Vector(totalDoFnum);
    f_acc = Vector(totalDoFnum);
    trial_r = Vector(totalDoFnum);
    pbc = Vector(fixedDoFnum);
    f = Vector(freeDoFnum - nodes->giveNumConstrDoFs() );
    residuals = Vector(totalDoFnum);
    ddr = Vector(freeDoFnum - nodes->giveNumConstrDoFs() );
    full_ddr = Vector(totalDoFnum);
    f_int_old = Vector(totalDoFnum);
    f_ext_old = Vector(totalDoFnum);
}

//////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////
// STEADY STATE LINEAR SOLVER
SteadyStateLinearSolver :: SteadyStateLinearSolver() {
    name = "SteadyStateLinearSolver";
    conj_grad_precision = 1e-13;
    conj_grad_relative_maxit = 0.85;
}

//////////////////////////////////////////////////////////
SteadyStateLinearSolver :: ~SteadyStateLinearSolver() {}

//////////////////////////////////////////////////////////
void SteadyStateLinearSolver :: prepareSystemMatricesAndInitialField(string init_r_file, string init_v_file, const bool initial) {
    Solver :: init(init_r_file, init_v_file, initial);

    //nodes->addRHS_nodalLoad(load, 0); //to correctly account for abrupt initial change of BC
    //nodes->updateDirrichletBC(r, 0); //to correctly account for abrupt initial change of BC

    //initial conditions
    if ( init_r_file.compare("") != 0 ) {
        r = nodes->readInitialConditions(init_r_file);
        computeInternalExternalForces(r, load, false, -1); //to activate initial conditions at elements
        elems->updateMaterialStatuses();
    }
    elems->prepareStiffnessMatrix(K);
    elems->updateStiffnessMatrix(K, "elastic");
}


//////////////////////////////////////////////////////////
void SteadyStateLinearSolver :: init(string init_r_file, string init_v_file, const bool initial) {
    prepareSystemMatricesAndInitialField(init_r_file, init_v_file, initial);
    computeKeff();
}

//////////////////////////////////////////////////////////
void SteadyStateLinearSolver :: computeKeff() {
    Keff = K;
    if ( nodes->giveConstraints()->isActive() ) {
        nodes->giveConstraints()->transformToConstraintSpace(Keff);
    }
}

//////////////////////////////////////////////////////////
Solver *SteadyStateLinearSolver :: readFromFile(const string filename) {
    string param, line;
    bool bdt, bttime;
    bdt = bttime = false;
    ifstream inputfile(filename.c_str() );
    if ( inputfile.is_open() ) {
        while ( getline(inputfile >> std :: ws, line) ) {
            if ( line.empty() ) {
                continue;
            }
            if ( line.at(0) == '#' ) {
                continue;
            }
            istringstream iss(line);
            iss >> param;
            if ( param.compare("time_step") == 0 ) {
                bdt = true;
                iss >> initdt;
                dt = initdt;
            } else if ( param.compare("total_time") == 0 ) {
                bttime = true;
                iss >> termination_time;
            } else if ( param.compare("conj_grad_precission") == 0 || param.compare("conj_grad_precision") == 0 ) {
                iss >> conj_grad_precision;
            } else if ( param.compare("conj_grad_relative_maxit") == 0 ) {
                iss >> conj_grad_relative_maxit;
            } else if ( param.compare("init_time") == 0 ) {
                iss >> this->init_time;
            } else if ( param.compare("init_step") == 0 ) {
                iss >> this->init_step;
            } else if ( param.compare("symsolver_type") == 0 ) {
                iss >> symsolver_type;
            }
        }
        inputfile.close();
    }
    if ( !bdt ) {
        cerr << name << ": solver parameter 'time_step' was not specified" << endl;
        exit(EXIT_FAILURE);
    }
    ;
    if ( !bttime ) {
        cerr << name << ": solver parameter 'total_time' was not specified" << endl;
        exit(EXIT_FAILURE);
    }

    return this;
};

//////////////////////////////////////////////////////////
void SteadyStateLinearSolver :: solve() {
    nodes->addRHS_nodalLoad(load, time);  //add nodal load
    nodes->updateDirrichletBC(trial_r, time); //give prescribed DoFs
    updateFieldVariables();
    computeForcesAtIntegrationTime(true);

    //solve linear system
    nodes->giveReducedForceArray(residuals, f);
    if ( LinalgSymmetricSolver(Keff, ddr, f, ddr, conj_grad_precision, conj_grad_relative_maxit, symsolver_type) == false ) {
        terminated = true;
        cerr << "Conjugate gradients did not converge during initialization of solver" << endl;
        return;
    }

    /*
     * cout << "----- K ----" << endl;
     * K.print();
     * cout << "----- d ----" << endl;
     * for(auto p:ddr ) cout << p<< endl;
     * cout << "----- f ----" << endl;
     * for(auto p:f ) cout << p<< endl;
     */
    // if ( ConjGrad(K, ddr, f, ddr, conj_grad_precision, conj_grad_relative_maxit) == false ) {
    // if (terminated == true);
    //     cerr << "Conjugate gradients did not converge" << endl;
    // }

    updateFieldVariables(); //calculate master fields at the step end
    computeForcesAtStepEnd(false); //to obtain the actual stress, fluxes, ...
}


//////////////////////////////////////////////////////////
void SteadyStateLinearSolver :: computeInternalExternalForces(const Vector &rr, const Vector &ll, const bool frozen, double timeStep) {
    nodes->updateSimplexVolumetricStrains(rr); //this line computes volumetric strain in simplices
    elems->integrateInternalForces(rr, f_int, frozen, timeStep);
    nodes->updateExternalForcesByReactions(f_int, ll, f_dam, f_acc, f_ext);     //give prescribed DoFs
    residuals = f_ext - f_int;
}

//////////////////////////////////////////////////////////
void SteadyStateLinearSolver :: runBeforeEachStep() {
    Solver :: runBeforeEachStep();
    cout << "######### Solving step " << step << " at time " << time << "; time step " << dt << " #########" << endl;
}

//////////////////////////////////////////////////////////
void SteadyStateLinearSolver :: updateFieldVariables() {
    nodes->giveFullDoFArray(ddr, full_ddr);
    trial_r += full_ddr;
}

//////////////////////////////////////////////////////////
void SteadyStateLinearSolver :: runAfterEachStep() {
    Solver :: runAfterEachStep();
}

//////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////
// STEADY STATE NONLINEAR SOLVER

SteadyStateNonLinearSolver :: SteadyStateNonLinearSolver() {
    name = "SteadyStateNonLinearSolver";
    idc = nullptr;

    W_ext_oldM = 0;
    W_int_oldM = 0;
    W_ext_oldT = 0;
    W_int_oldT = 0;
}

//////////////////////////////////////////////////////////
SteadyStateNonLinearSolver :: ~SteadyStateNonLinearSolver() {
    if ( idc ) {
        delete idc;
    }
}

//////////////////////////////////////////////////////////
void SteadyStateNonLinearSolver :: init(string init_r_file, string init_v_file, const bool initial) {
    SteadyStateLinearSolver :: init(init_r_file, init_v_file, initial);

    if ( idc ) {
        idc->init(nodes, funcs, initial);   //indirect displacement control
        ddf = Vector(freeDoFnum - nodes->giveNumConstrDoFs() );
        full_ddf = Vector(totalDoFnum);
        f_last_iter = Vector(freeDoFnum - nodes->giveNumConstrDoFs() );
        if ( initial ) {
            idc_time = 0;
            idc_dt = 1e-6;
            idc_time_converged = 0;
        } else {
            // JK: during solver initialization after geometry update, idc_time must be set to time of previously converged step
            idc_time = idc_time_converged;
        }
    }
}

//////////////////////////////////////////////////////////
Solver *SteadyStateNonLinearSolver :: readFromFile(const string filename) {
    SteadyStateLinearSolver :: readFromFile(filename);

    maxIt = 30;
    enlargeIt = shortenIt = 0;
    disErr = resErr = eneErr = 1e-5;
    limitEneErr = limitResErr = limitDisErr = 0;
    step_increase = 1.25;
    step_decrease = 0.8;
    critical_step_decrease = 0.5;
    stiffnessMatrixUpdate = 1e3;
    dampingMatrixUpdate = -1;
    massMatrixUpdate = -1;

    string param, line;
    dtmax = dtmin = dt;
    bool bdtmin = false;
    bool bdtmax = false;
    bool ben = false;
    bool bsh = false;
    unsigned helpuint;
    double valueIN;
    ifstream inputfile(filename.c_str() );
    if ( inputfile.is_open() ) {
        while ( getline(inputfile >> std :: ws, line) ) {
            if ( line.empty() ) {
                continue;
            }
            if ( line.at(0) == '#' ) {
                continue;
            }
            istringstream iss(line);
            iss >> param;
            if ( param.compare("max_time_step") == 0 ) {
                bdtmax = true;
                iss >> dtmax;
            } else if ( param.compare("min_time_step") == 0 ) {
                bdtmin = true;
                iss >> dtmin;
            } else if ( param.compare("tolerance") == 0 ) {
                iss >> disErr;
                resErr = eneErr = disErr;
            } else if ( param.compare("stiffness_matrix_update") == 0 ) {
                iss >> valueIN;
                stiffnessMatrixUpdate = int( valueIN + 0.5 );
            } else if ( param.compare("mass_matrix_update") == 0 ) {
                iss >> valueIN;
                massMatrixUpdate = int( valueIN + 0.5 );
            } else if ( param.compare("damping_matrix_update") == 0 ) {
                iss >> valueIN;
                dampingMatrixUpdate = int( valueIN + 0.5 );
            } else if ( param.compare("limit_tolerance") == 0 ) {
                iss >> valueIN;
                limitEneErr = limitResErr = limitDisErr = valueIN;
            } else if ( param.compare("maxIt") == 0 ) {
                iss >> maxIt;
                if ( maxIt < 1 ) {
                    std :: cerr << "number of itteration cannot be smaller than " << maxIt << "!!!,";
                    maxIt = 30;
                    std :: cout << " setting to default value: maxIt = " << maxIt << '\n';
                } else if ( maxIt < 3 ) {
                    std :: cout << "solver parameter maxIt set to " << maxIt << ", be carefull with such a small number" << '\n';
                }
            } else if ( param.compare("enlargeIt") == 0 ) {
                iss >> enlargeIt;
                ben = true;
            } else if ( param.compare("shortenIt") == 0 ) {
                iss >> shortenIt;
                bsh = true;
            } else if ( param.compare("step_increase") == 0 ) {
                iss >> valueIN;
                if ( valueIN < 1 ) {
                    std :: cerr << "step_increase cannot be smaller than 1! leaving default value " << step_increase << '\n';
                } else {
                    step_increase = valueIN;
                }
            } else if ( param.compare("step_decrease") == 0 ) {
                iss >> valueIN;
                if ( valueIN > 1 ) {
                    std :: cerr << "step_decrease cannot be greater than 1! leaving default value " << step_decrease << '\n';
                } else {
                    step_decrease = valueIN;
                }
            } else if ( param.compare("critical_step_decrease") == 0 ) {
                iss >> valueIN;
                if ( valueIN > 1 ) {
                    std :: cerr << "critical_step_decrease cannot be greater than 1! leaving default value " << critical_step_decrease << '\n';
                } else {
                    critical_step_decrease = valueIN;
                }
            } else if ( param.compare("indirect_displacement_control") == 0 ) {
                iss >> helpuint;
                if ( !idc ) {
                    idc = new IndirectDC();
                }
                idc->readFromStream(helpuint, inputfile);
            }
        }
        inputfile.close();
    }
    if ( !bdtmin && !bdtmax ) {
        cout << "fixed time step used" << endl;
    } else if ( bdtmin || bdtmax ) {
        cout << "adaptive time step used" << endl;
        if ( dtmin > dtmax ) {
            std :: cerr << "dtmin > dtmax, swapping values" << '\n';
            double dm = dtmax;
            dtmax = dtmin;
            dtmin = dm;
        }
        if ( dt < dtmin ) {
            std :: cerr << "dt < dtmin, setting dtmin = dt" << '\n';
            dtmin = dt;
        } else if ( dt > dtmax ) {
            std :: cerr << "dt > dtmax, setting dtmax = dt" << '\n';
            dtmax = dt;
        }
    } else {
        cout << endl;
    }
    if ( enlargeIt > shortenIt ) {
        std :: cerr << "cannot set number of iterations for step enlargement higher than number of iterations for step shortening, setting back to default values" << '\n';
        ben = false;
        bsh = false;
    }
    if ( !ben ) {
        enlargeIt = maxIt / 3;
    }
    if ( !bsh ) {
        shortenIt = maxIt / 2;
    }

    return this;
};



//////////////////////////////////////////////////////////
bool SteadyStateNonLinearSolver :: updateSystemMatrices(string matrixType, unsigned iteration) {
    if ( stiffnessMatrixUpdate == 0 || ( stiffnessMatrixUpdate > 0 && iteration % stiffnessMatrixUpdate == 0 ) ) {
        elems->updateStiffnessMatrix(K, matrixType);
        return true;
    } else {
        return false;
    }
}


//////////////////////////////////////////////////////////
void SteadyStateNonLinearSolver :: evaluateErrors(double *displa_error, double *energy_error, double *residu_error) {
    vector< bool >mechDoFs  = nodes->giveMechDoFsIndicator();   //these fields should be on input
    vector< bool >transpDoFs = nodes->giveTranspDoFsIndicator();
    double f_extM = 0;
    double f_intM = 0;
    double f_damM = 0;
    double f_accM = 0;
    double residualM = 0;
    double full_ddrM = 0;
    double trial_rM = 0;
    double energyM = 0;
    double f_extT = 0;
    double f_intT = 0;
    double f_damT = 0;
    double f_accT = 0;
    double residualT = 0;
    double full_ddrT = 0;
    double trial_rT = 0;
    double energyT = 0;
    W_intM = W_int_oldM;
    W_extM = W_ext_oldM;
    W_intT = W_int_oldT;
    W_extT = W_ext_oldT;
    W_kinM = 0;
    W_kinT = 0;

    for ( unsigned i = 0; i < totalDoFnum; i++ ) {
        if ( mechDoFs [ i ] ) {
            residualM += pow(residuals [ i ], 2);
            f_extM += pow(f_ext [ i ], 2);
            f_intM += pow(f_int [ i ], 2);
            f_damM += pow(f_dam [ i ], 2);
            f_accM += pow(f_acc [ i ], 2);
            full_ddrM += pow(full_ddr [ i ], 2);
            trial_rM += pow(trial_r [ i ], 2);
            energyM += abs(residuals [ i ] * full_ddr [ i ]);
            W_intM += abs(0.5 * ( f_int [ i ] + f_int_old [ i ] ) * ( r [ i ] - trial_r [ i ] ) );
            W_extM += abs(0.5 * ( f_ext [ i ] + f_ext_old [ i ] ) * ( r [ i ] - trial_r [ i ] ) );
            W_kinM += 0.; //TODO: correct kinetic energy
        }
        if ( transpDoFs [ i ] ) {
            residualT += pow(residuals [ i ], 2);
            f_extT += pow(f_ext [ i ], 2);
            f_intT += pow(f_int [ i ], 2);
            f_damT += pow(f_dam [ i ], 2);
            f_accT += pow(f_acc [ i ], 2);
            full_ddrT += pow(full_ddr [ i ], 2);
            trial_rT += pow(trial_r [ i ], 2);
            energyT += abs(residuals [ i ] * full_ddr [ i ]);
            W_intT += abs(0.5 * ( f_int [ i ] + f_int_old [ i ] ) * ( r [ i ] - trial_r [ i ] ) );
            W_extT += abs(0.5 * ( f_ext [ i ] + f_ext_old [ i ] ) * ( r [ i ] - trial_r [ i ] ) );
            W_kinT += 0.; //TODO: correct kinetic energy
        }
    }
    * residu_error = sqrt(residualM / max(max(max(f_extM, f_intM), max(f_damM, f_accM) ), EPS2displ) + residualT / max(max(max(f_extT, f_intT), max(f_damT, f_accT) ), EPS2press) );
    * displa_error = sqrt(full_ddrM / max(trial_rM, EPS2displ) + full_ddrT / max(trial_rT, EPS2press) );
    * energy_error = energyM / max(max(max(W_extM, W_intM), W_kinM), EPS2displ) + energyT / max(max(max(W_extT, W_intT), W_kinT), EPS2press);
}


//////////////////////////////////////////////////////////
void SteadyStateNonLinearSolver :: solve() {
    double load_mult;
    bool converged = false;
    bool restarted = false;
    double displa_error = 0;
    double energy_error = 0;
    double residu_error = 0;

    while ( !converged ) {
        //setup loading

        if ( !idc ) {
            nodes->addRHS_nodalLoad(load, time); //add nodal load
            nodes->updateDirrichletBC(trial_r, time); //give prescribed DoFs
            updateFieldVariables();
            computeForcesAtIntegrationTime(true);
        }

        unsigned it = 0;
        while ( !converged && it < maxIt ) {
            if ( updateSystemMatrices("secant", it) ) {
                computeKeff();                                    //only if required
            }
            nodes->giveReducedForceArray(residuals, f);

            if ( idc ) {      //indirect displacement control
                f_last_iter = f;
                load *= 0.;
                nodes->addRHS_nodalLoad(load, idc_time + idc_dt); //add nodal load
                nodes->updateDirrichletBC(trial_r, idc_time + idc_dt); //give prescribed DoFs
                computeForcesAtIntegrationTime(true);
                nodes->giveReducedForceArray(residuals, f);

                if ( LinalgSymmetricSolver(Keff, ddr, f_last_iter, ddr, conj_grad_precision, conj_grad_relative_maxit, symsolver_type) == false ) {
                    terminated = true;
                    cerr << "Conjugate gradients did not converge" << endl;
                    return;
                }
                if ( LinalgSymmetricSolver(Keff, ddf, f - f_last_iter, ddf, conj_grad_precision, conj_grad_relative_maxit, symsolver_type) == false ) {
                    terminated = true;
                    cerr << "Conjugate gradients did not converge" << endl;
                    return;
                }
                //nodes->giveFullDoFArray(ddr, full_ddr);
                nodes->giveFullDoFArray(ddf, full_ddf);
                //load_mult = idc->giveMultiplierCorrection(trial_r, full_ddr, full_ddf, time);
                load_mult = idc->giveMultiplierCorrection(trial_r, full_ddf, time);
                ddr += load_mult * ddf;
                idc_time += idc_dt * load_mult;

                load *= 0;
                nodes->addRHS_nodalLoad(load, idc_time); //add nodal load
                nodes->updateDirrichletBC(trial_r, idc_time); //give prescribed DoFs
            } else {         //direct controll
                if ( LinalgSymmetricSolver(Keff, ddr, f, ddr, conj_grad_precision, conj_grad_relative_maxit, symsolver_type) == false ) {
                    terminated = true;
                    cerr << "Conjugate gradients did not converge" << endl;
                    return;
                }
            }

            //update DoFs
            updateFieldVariables();
            //compute residuals
            computeForcesAtIntegrationTime(false); //to obtain the actual stress, fluxes, ...

            //compute and print errors
            evaluateErrors(& displa_error, & energy_error, & residu_error);
            if ( it == 0 ) {
                displa_error = 0;                        //error in displacement change, only from second iteration
            }
            cout << setw(6) << it << setw(15) << residu_error;
            if ( it == 0 ) {
                cout << setw(15) << "---";
            } else {
                cout << setw(15) << displa_error;
            }
            cout << setw(15) << energy_error << endl;

            if ( std :: isnan(residu_error) || std :: isnan(displa_error) || std :: isnan(energy_error) ) {
                cerr << "calculating with NaN in ";
                if ( std :: isnan(residu_error) ) {
                    std :: cerr << "\nresidua ";
                }
                if ( std :: isnan(displa_error) ) {
                    std :: cerr << "\ndisplacements ";
                }
                if ( std :: isnan(energy_error) ) {
                    std :: cerr << "\nenergies ";
                }
                std :: cerr << "- exit" << '\n';
                terminated = true;
                return;
            }

            if ( displa_error > disErr || residu_error > resErr || energy_error > eneErr ) {
                converged = false;
            } else {
                converged = true;
            }
            it++;
            //if (it>10) exit(1);
        }
        if ( !converged && dt > dtmin ) {
            time -= dt;
            dt = fmax(dt * critical_step_decrease, dtmin);
            time += dt;
            cerr << "Restarting step, timestep = " << dt << ", time = " << time << endl;
            restarted = true;
            trial_r = r;
            f_int = f_int_old;
            f_ext = f_ext_old;
            load *= 0;
            if ( idc ) {
                idc_time = idc_time_converged;
            }
        } else if ( !converged ) {
            if ( displa_error < limitDisErr && residu_error < limitResErr && energy_error < limitEneErr ) {
                std :: cerr << "tolerance increased in this step" << '\n';
                converged = true;
            } else {
                std :: cerr << "Error: Nonlinear static solver did not converge to the solution" << endl;
                terminated = true;
                return;
                // exit(1);
            }
        } else if ( ( !restarted ) && converged && it < enlargeIt && dt < dtmax ) {
            dt = fmin(dt * step_increase, dtmax);
            std :: cout << "enlarging step, timestep = " << dt << '\n';
        } else if ( converged && it > shortenIt && dt > dtmin ) {
            dt = fmax(dt * step_decrease, dtmin);
            std :: cout << "shortening step, timestep = " << dt << '\n';
        }
    }

    computeForcesAtStepEnd(false); //to obtain the actual stress, fluxes, ...
}

//////////////////////////////////////////////////////////
void SteadyStateNonLinearSolver :: runBeforeEachStep() {
    SteadyStateLinearSolver :: runBeforeEachStep();

    cout <<  scientific; //cout << setprecision(8);
    cout << "----------------------------------------------------" << endl;
    cout << setw(6) << "iter." << setw(15) << "residual" << setw(15) << "displacement" << setw(15) << "energy error" << endl;
    cout << setw(6) << " " << setw(15) << resErr << setw(15) << disErr << setw(15) << eneErr << endl;
    cout << "----------------------------------------------------" << endl;
}

//////////////////////////////////////////////////////////
void SteadyStateNonLinearSolver :: printAllVectors() {
    std :: cout << "trial_r\t r\t full_ddr\t ddr\t f_int\t f_ext\t load\t f_reduced\t" << '\n';
    for ( unsigned i = 0; i < trial_r.size(); i++ ) {
        std :: cout << trial_r [ i ] << '\t';
        std :: cout << r [ i ] << '\t';
        std :: cout << full_ddr [ i ] << '\t';
        if ( i < ddr.size() ) {
            std :: cout << ddr [ i ] << '\t';
        } else {
            std :: cout << 'H' << '\t';
        }
        std :: cout << f_int [ i ] << '\t';
        std :: cout << f_ext [ i ] << '\t';
        std :: cout << load [ i ] << '\t';
        if ( i < f.size() ) {
            std :: cout << f [ i ] << '\n';
        } else {
            std :: cout << 'H' << '\n';
        }
    }
}

//////////////////////////////////////////////////////////
void SteadyStateNonLinearSolver :: runAfterEachStep() {
    if ( !terminated ) {
        SteadyStateLinearSolver :: runAfterEachStep();

        W_int_oldM = W_intM;
        W_ext_oldM = W_extM;
        W_int_oldT = W_intT;
        W_ext_oldT = W_extT;
        cout << "----------------------------------------------------" << endl;

        if ( idc ) {
            idc_time_converged = idc_time;
        }
    }
}

//////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////
// TRANSIENT LINEAR TRANSPORT SOLVER
TransientLinearTransportSolver :: TransientLinearTransportSolver() {
    name = "TransientLinearTransportSolver";
    applySpectralRadius(0.8);
}

//////////////////////////////////////////////////////////
void TransientLinearTransportSolver :: applySpectralRadius(double rhoinfty) {
    //set up the generalized-alpha method according to Chung and Hulbert 1993
    if ( abs(rhoinfty - 0.5) > 0.5 ) {
        cerr << "Error in solver: spectral radius must be inside interval [0,1]" << endl;
        exit(1);
    }

    alpha_m = ( 2. * rhoinfty - 1. ) / ( 1. + rhoinfty );
    alpha_f = rhoinfty / ( 1. + rhoinfty );
    gamma = 1. / 2. - alpha_m + alpha_f;
    beta = 1.; //not used
}

//////////////////////////////////////////////////////////
TransientLinearTransportSolver :: ~TransientLinearTransportSolver() {}

//////////////////////////////////////////////////////////
void TransientLinearTransportSolver :: prepareSystemMatricesAndInitialField(string init_r_file, string init_v_file, const bool initial) {
    SteadyStateLinearSolver :: prepareSystemMatricesAndInitialField(init_r_file, init_v_file, initial);

    v_old = Vector(totalDoFnum);
    elems->prepareDampingMatrix(C);
    elems->updateDampingMatrix(C);
}

//////////////////////////////////////////////////////////
void TransientLinearTransportSolver :: init(string init_r_file, string init_v_file, const bool initial) {
    prepareSystemMatricesAndInitialField(init_r_file, init_v_file, initial);
    computeKeff();

    //compute initial presure derivative
    nodes->giveReducedForceArray(residuals, f);
    CoordinateIndexedSparseMatrix Cred(C);
    if ( nodes->giveConstraints()->isActive() ) {
        nodes->giveConstraints()->transformToConstraintSpace(Cred);
    }
    if ( LinalgSymmetricSolver(Cred, ddr, f,  ddr, conj_grad_precision, conj_grad_relative_maxit, symsolver_type) == false ) {
        terminated = true;
        cerr << "Conjugate gradients did not converge during initialization of solver" << endl;
        exit(1);
    }
    v = Vector(totalDoFnum);
    nodes->giveFullDoFArray(ddr, v);
}

//////////////////////////////////////////////////////////
void TransientLinearTransportSolver :: computeForcesAtIntegrationTime(const bool frozen) {
    elems->integrateDampingForces(v * ( 1. - alpha_m ) +  v_old * alpha_m, f_dam);
    computeInternalExternalForces(r * alpha_f + trial_r * ( 1. - alpha_f ), load_old * alpha_f + load * ( 1. - alpha_f ), frozen, dt * (1.-alpha_f) );
    residuals -= f_dam;
}

//////////////////////////////////////////////////////////
void TransientLinearTransportSolver :: computeForcesAtStepEnd(const bool frozen) {
    elems->integrateDampingForces(v, f_dam);
    computeInternalExternalForces(trial_r, load, frozen, dt);
    residuals -= f_dam;
}

//////////////////////////////////////////////////////////
void TransientLinearTransportSolver :: computeKeff() {
    Keff = C * ( ( 1. - alpha_m ) / ( dt * gamma ) ) + K * ( 1. - alpha_f );
    if ( nodes->giveConstraints()->isActive() ) {
        nodes->giveConstraints()->transformToConstraintSpace(Keff);
    }
}

//////////////////////////////////////////////////////////
//void TransientLinearTransportSolver :: updateFeff() {
//    nodes->giveReducedDoFArray(v,v_red);
//    feff =  f - ( C * v_red ) * ( 1. + ( alpha_m - 1. ) / gamma ); //not used in the code, available for checking
//}

//////////////////////////////////////////////////////////
void TransientLinearTransportSolver :: updateFieldVariables() {
    SteadyStateLinearSolver :: updateFieldVariables();
    v = ( trial_r - r ) * ( 1. / ( dt * gamma ) ) + v_old * ( 1. - 1. / gamma );
}

//////////////////////////////////////////////////////////
void TransientLinearTransportSolver :: solve() {
    SteadyStateLinearSolver :: solve();
}

//////////////////////////////////////////////////////////
void TransientLinearTransportSolver :: runBeforeEachStep() {
    SteadyStateLinearSolver :: runBeforeEachStep();
    v_old = v;
}

//////////////////////////////////////////////////////////
void TransientLinearTransportSolver :: runAfterEachStep() {
    SteadyStateLinearSolver :: runAfterEachStep();
};

//////////////////////////////////////////////////////////
Solver *TransientLinearTransportSolver :: readFromFile(const string filename) {
    SteadyStateNonLinearSolver :: readFromFile(filename);

    double num;
    string param, line;
    ifstream inputfile(filename.c_str() );
    if ( inputfile.is_open() ) {
        while ( getline(inputfile >> std :: ws, line) ) {
            if ( line.empty() ) {
                continue;
            }
            if ( line.at(0) == '#' ) {
                continue;
            }
            istringstream iss(line);
            iss >> param;
            if ( param.compare("spectral_radius") == 0 ) {
                iss >> num;
                applySpectralRadius(num);
            }
        }
        inputfile.close();
    }
    return this;
};

//////////////////////////////////////////////////////////
bool TransientLinearTransportSolver :: updateSystemMatrices(string matrixType, unsigned iteration) {
    bool updated0 = SteadyStateNonLinearSolver :: updateSystemMatrices(matrixType, iteration);
    bool updated1 = false;
    if ( dampingMatrixUpdate == 0 || ( dampingMatrixUpdate > 0 && iteration % dampingMatrixUpdate == 0 ) ) {
        elems->updateDampingMatrix(C);
        updated1 = true;
    }
    return ( updated0 || updated1 );
}

//////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////
// TRANSIENT NON-LINEAR TRANSPORT SOLVER
TransientNonLinearTransportSolver :: TransientNonLinearTransportSolver() {
    name = "TransientNonLinearTransportSolver";
}

//////////////////////////////////////////////////////////
TransientNonLinearTransportSolver :: ~TransientNonLinearTransportSolver() {}

//////////////////////////////////////////////////////////
void TransientNonLinearTransportSolver :: init(string init_r_file, string init_v_file, const bool initial) {
    TransientLinearTransportSolver :: init(init_r_file, init_v_file, initial);
}


//////////////////////////////////////////////////////////
void TransientNonLinearTransportSolver :: runBeforeEachStep() {
    SteadyStateNonLinearSolver :: runBeforeEachStep();
    v_old = v;
}

//////////////////////////////////////////////////////////
void TransientNonLinearTransportSolver :: runAfterEachStep() {
    SteadyStateNonLinearSolver :: runAfterEachStep();
};

//////////////////////////////////////////////////////////
void TransientNonLinearTransportSolver :: solve() {
    SteadyStateNonLinearSolver :: solve();
}

//////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////
// TRANSIENT LINEAR MECHANICAL SOLVER
TransientLinearMechanicalSolver :: TransientLinearMechanicalSolver() {
    name = "TransientLinearMechanicalSolver";
}

//////////////////////////////////////////////////////////
TransientLinearMechanicalSolver :: ~TransientLinearMechanicalSolver() {}

//////////////////////////////////////////////////////////
void TransientLinearMechanicalSolver :: solve() {
    SteadyStateLinearSolver :: solve();
}

//////////////////////////////////////////////////////////
void TransientLinearMechanicalSolver :: prepareSystemMatricesAndInitialField(string init_r_file, string init_v_file, const bool initial) {
    //initial conditions
    if ( init_v_file.compare("") != 0 ) {
        v = nodes->readInitialConditions(init_r_file);
    } else {
        v = Vector(totalDoFnum);
    }
    v_old = Vector(totalDoFnum);
    a_old = Vector(totalDoFnum);
    elems->prepareMassMatrix(M);
    elems->updateMassMatrix(M);

    TransientLinearTransportSolver :: prepareSystemMatricesAndInitialField(init_r_file, init_v_file, initial);
}

//////////////////////////////////////////////////////////
void TransientLinearMechanicalSolver :: init(string init_r_file, string init_v_file, const bool initial) {
    prepareSystemMatricesAndInitialField(init_r_file, init_v_file, initial);
    computeKeff();

    //compute initial acceleration
    nodes->giveReducedForceArray(residuals, f);
    CoordinateIndexedSparseMatrix Cred(C);
    CoordinateIndexedSparseMatrix Mred(M);
    if ( nodes->giveConstraints()->isActive() ) {
        nodes->giveConstraints()->transformToConstraintSpace(Cred);
        nodes->giveConstraints()->transformToConstraintSpace(Mred);
    }
    Vector v_red(freeDoFnum - nodes->giveNumConstrDoFs() );
    nodes->giveReducedDoFArray(v, v_red);
    terminated = !LinalgSymmetricSolver(Mred, ddr, f - Cred * v_red,  ddr, conj_grad_precision, conj_grad_relative_maxit, symsolver_type);
    a = Vector(totalDoFnum);
    nodes->giveFullDoFArray(ddr, a);
}

//////////////////////////////////////////////////////////
void TransientLinearMechanicalSolver :: applySpectralRadius(double rhoinfty) {
    //set up the generalized-alpha method according to Chung and Hulbert 1993
    if ( abs(rhoinfty - 0.5) > 0.5 ) {
        cerr << "Error in solver: spectral radius must be inside interval [0,1]" << endl;
        exit(1);
    }

    alpha_m = 0.5 * ( 3. - rhoinfty ) / ( 1. + rhoinfty );
    alpha_f = 1. / ( 1. + rhoinfty );
    gamma = 1. / 2. - alpha_m + alpha_f;
    beta = 0.25 * pow(1. + alpha_m - alpha_f, 2);
}

//////////////////////////////////////////////////////////
void TransientLinearMechanicalSolver :: computeKeff() {
    Keff = K * ( 1. - alpha_f ) + C * ( ( 1 - alpha_f ) * gamma / dt / beta ) + M * ( ( 1 - alpha_m ) / dt / dt / beta );
    if ( nodes->giveConstraints()->isActive() ) {
        nodes->giveConstraints()->transformToConstraintSpace(Keff);
    }
}

//////////////////////////////////////////////////////////
//void TransientLinearMechanicalSolver :: updateFeff() {
//feff = f - C * ( ( 1. + gamma * ( alpha_f - 1 ) / beta ) * v_red + dt * ( 1. - alpha_f ) * ( 1. - gamma / ( 2. * beta ) ) * a_red ) - M * ( ( alpha_m - 1. ) / ( dt * beta ) * v_red  + ( 1. + ( alpha_m - 1 ) / ( 2. * beta ) ) * a_red );
//}

//////////////////////////////////////////////////////////
void TransientLinearMechanicalSolver :: updateFieldVariables() {
    SteadyStateLinearSolver :: updateFieldVariables();
    v = ( trial_r - r ) * ( gamma / ( dt * beta ) ) + v_old * ( 1. - gamma / beta ) + a_old * ( dt * ( 1. - gamma / ( 2. * beta ) ) );
    a = ( ( trial_r - r ) - v_old * dt - a_old * ( pow(dt, 2) * ( 1. / 2. - beta ) ) ) / ( pow(dt, 2) * beta );
}


//////////////////////////////////////////////////////////
void TransientLinearMechanicalSolver :: computeForcesAtIntegrationTime(const bool frozen) {
    elems->integrateDampingForces(v * ( 1. - alpha_f ) +  v_old * alpha_f, f_dam);
    elems->integrateInertiaForces(a * ( 1. - alpha_m ) +  a_old * alpha_m, f_acc);
    computeInternalExternalForces(r * alpha_f + trial_r * ( 1. - alpha_f ), load_old * alpha_f + load * ( 1. - alpha_f ), frozen, dt * (1.-alpha_f) );
    residuals -= f_dam + f_acc;
}

//////////////////////////////////////////////////////////
void TransientLinearMechanicalSolver :: computeForcesAtStepEnd(const bool frozen) {
    elems->integrateDampingForces(v, f_dam);
    //elems->integrateInertiaForces( a, f_acc );
    computeInternalExternalForces(trial_r, load, frozen, dt);
    residuals -= f_dam + f_acc;
}

//////////////////////////////////////////////////////////
void TransientLinearMechanicalSolver :: runBeforeEachStep() {
    TransientLinearTransportSolver :: runBeforeEachStep();
    a_old = a;
}

//////////////////////////////////////////////////////////
void TransientLinearMechanicalSolver :: runAfterEachStep() {
    TransientLinearTransportSolver :: runAfterEachStep();
};


//////////////////////////////////////////////////////////
bool TransientLinearMechanicalSolver :: updateSystemMatrices(string matrixType, unsigned iteration) {
    bool updated0 = TransientLinearTransportSolver :: updateSystemMatrices(matrixType, iteration);
    bool updated1 = false;
    if ( massMatrixUpdate == 0 || ( massMatrixUpdate > 0 && iteration % massMatrixUpdate == 0 ) ) {
        elems->updateMassMatrix(M);
        updated1 = true;
    }
    return ( updated0 || updated1 );
}


//////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////
// TRANSIENT NON LINEAR MECHANICAL SOLVER
TransientNonLinearMechanicalSolver :: TransientNonLinearMechanicalSolver() {
    name = "TransientNonLinearMechanicalSolver";
}

//////////////////////////////////////////////////////////
TransientNonLinearMechanicalSolver :: ~TransientNonLinearMechanicalSolver() {}

//////////////////////////////////////////////////////////
void TransientNonLinearMechanicalSolver :: init(string init_r_file, string init_v_file, const bool initial) {
    TransientLinearMechanicalSolver :: init(init_r_file, init_v_file, initial);
}

//////////////////////////////////////////////////////////
void TransientNonLinearMechanicalSolver :: runBeforeEachStep() {
    TransientNonLinearTransportSolver :: runBeforeEachStep();
    a_old = a;
}

//////////////////////////////////////////////////////////
void TransientNonLinearMechanicalSolver :: runAfterEachStep() {
    TransientNonLinearTransportSolver :: runAfterEachStep();
};


//////////////////////////////////////////////////////////
void TransientNonLinearMechanicalSolver :: solve() {
    SteadyStateNonLinearSolver :: solve();
}
