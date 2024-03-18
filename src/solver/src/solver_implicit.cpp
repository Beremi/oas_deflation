#include "solver_implicit.h"
#include "adaptivity.h"
#define numPhysicalFields 4

using namespace std;

//////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////
// STEADY STATE LINEAR SOLVER
SteadyStateLinearSolver :: SteadyStateLinearSolver() {
    name = "SteadyStateLinearSolver";
    conj_grad_precision = 1e-14;
    conj_grad_relative_maxit = 0.85;
    isTimeReal = false;
}

//////////////////////////////////////////////////////////
SteadyStateLinearSolver :: ~SteadyStateLinearSolver() {
    //delete linalgsolver;
}

//////////////////////////////////////////////////////////
void SteadyStateLinearSolver :: prepareSystemMatricesAndInitialField(string init_r_file, string init_v_file, const bool initial) {
    ( void ) init_v_file;
    ( void ) initial;
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
    Solver :: init(init_r_file, init_v_file, initial);
    prepareSystemMatricesAndInitialField(init_r_file, init_v_file, initial);
    computeKeff();
}

//////////////////////////////////////////////////////////
void SteadyStateLinearSolver :: rebuild() {
    Solver :: rebuild();
    prepareSystemMatricesAndInitialField("", "", false);
    computeKeff();
}

//////////////////////////////////////////////////////////
void SteadyStateLinearSolver :: computeKeff() {
    Keff = K;
    if ( nodes->giveConstraints()->isActive() ) {
        nodes->giveConstraints()->transformToConstraintSpace(Keff);
    }
    factorizeLinearSystem();
}

//////////////////////////////////////////////////////////
void SteadyStateLinearSolver :: reset() {
    load.setZero();
    nodes->addRHS_nodalLoad(load, time); //add nodal load
    nodes->updateDirrichletBC(trial_r, time); //give prescribed DoFs
    computeForcesAtIntegrationTime(true);
    nodes->giveReducedForceArray(residuals, f);

    /*
     * if ( LinalgSymmetricSolver(Keff, ddr, f, ddr, conj_grad_precision, conj_grad_relative_maxit, symsolver_type) == false ) {
     *  std :: cerr << "Conjugate gradients did not converge, attempt to restart step" << endl;
     *  std :: cerr << "solver restart did not work" << endl;
     *  exit(1);
     * }
     */
    linalgsolver->solve(ddr, f);

    updateFieldVariables();
    computeForcesAtIntegrationTime(true);

    runAfterEachStep();
}

//////////////////////////////////////////////////////////
Solver *SteadyStateLinearSolver :: readFromFile(const string filename) {
    string param, line;
    bool bdt, bttime;
    bdt = bttime = false;
    ifstream inputfile( filename.c_str() );
    if ( inputfile.is_open() ) {
        while ( getline(inputfile >> std :: ws, line) ) {
            if ( line.empty() || ( line.at(0) == '#' ) ) {
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
            } else if ( param.compare("conj_grad_precision") == 0 ) {
                iss >> conj_grad_precision;
            } else if ( param.compare("conj_grad_relative_maxit") == 0 ) {
                iss >> conj_grad_relative_maxit;
            } else if ( param.compare("init_time") == 0 ) {
                iss >> this->init_time;
            } else if ( param.compare("init_step") == 0 ) {
                iss >> this->init_step;
            } else if ( param.compare("solver_type") == 0 ) {
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
    updateFieldVariables();      //with ddr=0
    computeForcesAtIntegrationTime(true);

    //solve linear system
    /*
     * nodes->giveReducedForceArray(residuals, f);
     * if ( LinalgSymmetricSolver(Keff, ddr, f, ddr, conj_grad_precision, conj_grad_relative_maxit, symsolver_type) == false ) {
     *  terminated = true;
     *  cerr << "Conjugate gradients did not converge during initialization of solver" << endl;
     *  return;
     * }
     */
    linalgsolver->solve(ddr, f);
   
    /*
     * cout << "----- K ----" << endl;
     * Keff.print();
     * cout << "----- d ----" << endl;
     * for(auto p:ddr ) cout << p<< endl;
     * cout << "----- f ----" << endl;
     * for(auto p:f ) cout << p<< endl;
     */

    updateFieldVariables(); //calculate master fields at the step end
    computeForcesAtStepEnd(false); //to obtain the actual stress, fluxes, ...
}


//////////////////////////////////////////////////////////
void SteadyStateLinearSolver :: runBeforeEachStep() {
    Solver :: runBeforeEachStep();
    trial_r = r;
    cout << "######### Solving step " << step << " at time " << time << "; time step " << dt << " #########" << endl;
}

//////////////////////////////////////////////////////////
void SteadyStateLinearSolver :: runAfterEachStep() {
    Solver :: runAfterEachStep();
}


//////////////////////////////////////////////////////////
void SteadyStateLinearSolver :: factorizeLinearSystem() {

    cout << "factorizing system matrix" << endl;

    if ( linalgsolver == nullptr ) {
        if ( symsolver_type == "EigenConj" ) {
            std::unique_ptr<ConjGradSolver> cgs = std::make_unique< ConjGradSolver >();
            cgs->setPrecisionAndRelMaxIters(conj_grad_precision, conj_grad_relative_maxit);
            linalgsolver = std::move(cgs);
        } else if  ( symsolver_type == "EigenLDLT" ) {
            linalgsolver = std::make_unique< LDLTSolver >();
        } else if  ( symsolver_type == "EigenLLT" ) {
            linalgsolver = std::make_unique< LLTSolver >();
        } else if  ( symsolver_type == "EigenSparseLU" ) {
            linalgsolver = std::make_unique< LUSolver >();
        } else {
            cerr << "Solver type " << symsolver_type << " is not implemented" << endl;
            exit(1);
        }
    }    
    linalgsolver->factorize(Keff);
}

//////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////
// STEADY STATE NONLINEAR SOLVER

SteadyStateNonLinearSolver :: SteadyStateNonLinearSolver() {
    name = "SteadyStateNonLinearSolver";
    idc = nullptr;

    EPS2.resize(4);
    EPS2 [ 0 ] = 1e-20; //mechanics
    EPS2 [ 1 ] = 1e-12; //transport
    EPS2 [ 2 ] = 1e-17; //temperature
    EPS2 [ 3 ] = 1e-18; //humidity

    maxIt = 30;
    minIt = 1;
    enlargeIt = shortenIt = 0;
    maxDisErr = maxResErr = maxEneErr = 1e-5;
    limitEneErr = limitResErr = limitDisErr = 0;
    step_increase = 1.25;
    step_decrease = 0.8;
    critical_step_decrease = 0.5;
    stiffnessMatrixUpdate = 1e3;
    dampingMatrixUpdate = -1;
    massMatrixUpdate = -1;

    it = 0;
    restarts = 0;
}

//////////////////////////////////////////////////////////
SteadyStateNonLinearSolver :: ~SteadyStateNonLinearSolver() {
    if ( idc ) {
        delete idc;
    }
}

//////////////////////////////////////////////////////////
void SteadyStateNonLinearSolver :: init(string init_r_file, string init_v_file, const bool initial) {
    elems->readMatStatsFromFile(this->init_time, this->init_step, this->initdt, this->init_idc_time);
    SteadyStateLinearSolver :: init(init_r_file, init_v_file, initial);
    this->fully_converged = false;
    if ( idc ) {
        idc->init(nodes, funcs, initial);   //indirect displacement control
        ddf = Vector :: Zero(freeDoFnum);
        full_ddf = Vector :: Zero(totalDoFnum);
        f_last_iter = Vector :: Zero(freeDoFnum);
        if ( initial ) {
            idc_time = this->init_idc_time;
            idc_dt = 1e-6;
            idc_time_converged = this->init_idc_time;
        } else {
            // JK: during solver initialization after geometry update, idc_time must be set to time of previously converged step
            idc_time = idc_time_converged;
        }
    }

    computeForcesAtIntegrationTime(true); //to initialize all fields in the model
}

//////////////////////////////////////////////////////////
Solver *SteadyStateNonLinearSolver :: readFromFile(const string filename) {
    SteadyStateLinearSolver :: readFromFile(filename);

    string param, line;
    dtmax = dtmin = dt;
    bool bdtmin = false;
    bool bdtmax = false;
    bool ben = false;
    bool bsh = false;
    unsigned helpuint;
    double valueIN;
    ifstream inputfile(filename.c_str(), ios :: in | ios :: binary);
    if ( inputfile.is_open() ) {
        while ( getline(inputfile >> std :: ws, line) ) {
            if ( line.empty() || ( line.at(0) == '#' ) ) {
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
                iss >> maxDisErr;
                maxResErr = maxEneErr = maxDisErr;
            } else if ( param.compare("stiffness_matrix_update") == 0 ) {
                iss >> valueIN;
                stiffnessMatrixUpdate = int( valueIN );
            } else if ( param.compare("mass_matrix_update") == 0 ) {
                iss >> valueIN;
                massMatrixUpdate = int( valueIN );
            } else if ( param.compare("damping_matrix_update") == 0 ) {
                iss >> valueIN;
                dampingMatrixUpdate = int( valueIN );
            } else if ( param.compare("limit_tolerance") == 0 ) {
                iss >> valueIN;
                limitEneErr = limitResErr = limitDisErr = valueIN;
            } else if ( param.compare("minIt") == 0 || param.compare("min_iterations") == 0 ) {
                iss >> minIt;
            } else if ( param.compare("maxIt") == 0 || param.compare("max_iterations") == 0 ) {
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
        //Jan E: NO, let's tolerate that first step can be out of bounds
        //if ( dt < dtmin ) {
        //    std :: cerr << "dt < dtmin, setting dtmin = dt" << '\n';
        //    dtmin = dt;
        //} else if ( dt > dtmax ) {
        //    std :: cerr << "dt > dtmax, setting dtmax = dt" << '\n';
        //    dtmax = dt;
        //}
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
void SteadyStateNonLinearSolver :: giveValues(string code, Vector &result) const {
    if ( code.compare("iterations") == 0 ) {
        result.resize(1);
        result [ 0 ] = it;
    } else if ( code.compare("restarts") == 0 ) {
        result.resize(1);
        result [ 0 ] = restarts;
    } else if ( code.compare("error_dofs") == 0 || code.compare("error_displacements") == 0 ) {
        result.resize(1);
        result [ 0 ] = disErr;
    } else if ( code.compare("error_residuals") == 0 ) {
        result.resize(1);
        result [ 0 ] = resErr;
    } else if ( code.compare("error_energy") == 0 ) {
        result.resize(1);
        result [ 0 ] = eneErr;
    } else if ( code.compare("idc_time") == 0 ) {
        result.resize(1);
        result [ 0 ] = idc_time;
    } else if ( code.compare("converged") == 0 ) {
        result.resize(1);
        result [ 0 ] = fully_converged;
    } else {
        SteadyStateLinearSolver :: giveValues(code, result);
    }
}

//////////////////////////////////////////////////////////
bool SteadyStateNonLinearSolver :: updateSystemMatrices(string matrixType, unsigned iteration) {
    if ( stiffnessMatrixUpdate == 0 || ( stiffnessMatrixUpdate > 0 && iteration % abs(stiffnessMatrixUpdate) == 0 ) ) {
        elems->updateStiffnessMatrix(K, matrixType);
        return true;
    } else {
        return false;
    }
}


//////////////////////////////////////////////////////////
void SteadyStateNonLinearSolver :: evaluateErrors() {
    computeTotalInternalAndExternalAndKineticEnergy();

    vector< unsigned >pf  = nodes->givePhysicalFieldsOfDoFs();
    Vector f_extPF = Vector :: Zero(numPhysicalFields);
    Vector f_intPF = Vector :: Zero(numPhysicalFields);
    Vector f_damPF = Vector :: Zero(numPhysicalFields);
    Vector f_accPF = Vector :: Zero(numPhysicalFields);
    Vector residualPF = Vector :: Zero(numPhysicalFields);
    Vector full_ddrPF = Vector :: Zero(numPhysicalFields);
    Vector trial_rPF = Vector :: Zero(numPhysicalFields);
    Vector energyPF = Vector :: Zero(numPhysicalFields);
    
    unsigned pff;
    for ( unsigned i = 0; i < totalDoFnum; i++ ) {
        pff = pf [ i ];
        residualPF [ pff ] += pow(residuals [ i ], 2);
        f_extPF [ pff ] += pow(f_ext [ i ], 2);
        f_intPF [ pff ] += pow(f_int [ i ], 2);
        f_damPF [ pff ] += pow(f_dam [ i ], 2);
        f_accPF [ pff ] += pow(f_acc [ i ], 2);
        full_ddrPF [ pff ] += pow(full_ddr [ i ], 2);
        trial_rPF [ pff ] += pow(trial_r [ i ], 2);
        energyPF [ pff ] += residuals [ i ] * full_ddr [ i ];
    }  
 
    resErr = disErr = eneErr = 0;
    for ( unsigned i = 0; i < numPhysicalFields; i++ ) {
        resErr += residualPF [ i ] / max(max( max(f_extPF [ i ], f_intPF [ i ]), max(f_damPF [ i ], f_accPF [ i ]) ), EPS2 [ i ]);
        disErr += full_ddrPF [ i ] / max(trial_rPF [ i ], EPS2 [ i ]);
        eneErr += abs(energyPF [ i ]) / max(max(max(abs(W_ext [ i ]), abs(W_int [ i ])), abs(W_kin [ i ])), EPS2 [ i ]);
    }
    resErr = sqrt(resErr);
    disErr = sqrt(disErr);
    eneErr = sqrt(eneErr);
}

//////////////////////////////////////////////////////////
void SteadyStateNonLinearSolver :: reset() {
    load.setZero();
    ddr.setZero();
    double reset_time = time;
    if ( idc ) {
        reset_time = idc_time;
    }

    bool converged = false;

    restarts = 0;
    while ( !converged ) {
        nodes->addRHS_nodalLoad(load, reset_time); //add nodal load
        nodes->updateDirrichletBC(trial_r, reset_time); //give prescribed DoFs
        updateFieldVariables();     //with ddr=0
        computeForcesAtIntegrationTime(true);

        it = 0;
        while ( !converged && it < maxIt ) {
            if ( updateSystemMatrices("secant", it) ) {
                computeKeff();                                    //only if required
            }
            nodes->giveReducedForceArray(residuals, f);   // NOTE JK when IDC applied and step reset, residuals from the last iteration are used here //JE: no, they are actually computed again here

            /*
             * if ( LinalgSymmetricSolver(Keff, ddr, f, ddr, conj_grad_precision, conj_grad_relative_maxit, symsolver_type) == false ) {
             *  std :: cerr << "Conjugate gradients did not converge, attempt to restart step" << endl;
             *  it = maxIt;
             *  resErr = 1e10;
             *  break;
             * }
             */

            linalgsolver->solve(ddr, f);

            //update DoFs
            updateFieldVariables();
            //compute residuals
            computeForcesAtIntegrationTime(true); //to obtain the actual stress, fluxes, ...

            //compute and print errors
            evaluateErrors();
            if ( it == 0 ) {
                disErr = 0;                        //error in displacement change, only from second iteration
            }
            cout << setw(6) << it << setw(15) << resErr;
            if ( it == 0 ) {
                cout << setw(15) << "---";
            } else {
                cout << setw(15) << disErr;
            }
            cout << setw(15) << eneErr << endl;

            if ( std :: isnan(resErr) || std :: isnan(disErr) || std :: isnan(eneErr) ) {
                std :: cerr << "calculating with NaN in ";
                if ( std :: isnan(resErr) ) {
                    std :: cerr << "\tresiduals ";
                }
                if ( std :: isnan(disErr) ) {
                    std :: cerr << "\tdisplacements ";
                }
                if ( std :: isnan(eneErr) ) {
                    std :: cerr << "\tenergies ";
                }
                std :: cerr << endl;
                it = maxIt;
                resErr = 1e10;
                break;
            }

            if ( disErr <= maxDisErr && resErr <= maxResErr && eneErr <= maxEneErr ) {
                converged = true;
            } else {
                converged = false;
            }
            it++;
        }

        if ( converged ) {
            this->fully_converged = true;
        } else if ( !converged ) {
            if ( disErr < limitDisErr && resErr < limitResErr && eneErr < limitEneErr ) {
                std :: cout << "tolerance increased in this step" << '\n';
                converged = true;
                this->fully_converged = false;
                computeForcesAtStepEnd(false); //to obtain the actual stress, fluxes, ...
            } else {
                std :: cerr << "Error: " << name << " did not converge to the solution" << endl;
                terminated = true;
                return;
            }
        }
    }
    runAfterEachStep();
}

//////////////////////////////////////////////////////////
void SteadyStateNonLinearSolver :: solve() {
    double load_mult;
    bool converged = false;
    bool restarted = false;
    bool restart_now = false;
    Vector help_idc_r, help_idc_f;

    //MyVector reset_residuals = residuals;   ///> if step restarted when IDC applied, residuals need to be reset to stage before the step start
    //JE: no, these are recomputed
    restarts = 0;
    while ( !converged ) {
        //setup loading

        if ( !idc ) {
            nodes->addRHS_nodalLoad(load, time); //add nodal load
            nodes->updateDirrichletBC(trial_r, time); //give prescribed DoFs
            updateFieldVariables();      //with ddr=0
            computeForcesAtIntegrationTime(true);
        } else {
            help_idc_r =  Vector :: Zero(totalDoFnum);
            help_idc_f =  Vector :: Zero(totalDoFnum);
        }

        it = 0;
        while ( !converged && it < maxIt ) {
            if ( updateSystemMatrices("secant", it) ) {
                computeKeff();                                    //only if required
            }
            nodes->giveReducedForceArray(residuals, f);   // NOTE JK when IDC applied and step reset, residuals from the last iteration are used here //JE: no, they are actually computed again here

            if ( idc ) {      //indirect displacement control
                f_last_iter = f;
                load.setZero();
                nodes->addRHS_nodalLoad(load, idc_time + idc_dt); //add nodal load
                nodes->updateDirrichletBC(trial_r, idc_time + idc_dt); //give prescribed DoFs
                //updateFieldVariables(); //if used ddr needs to be set to zero
                computeForcesAtIntegrationTime(true);
                nodes->giveReducedForceArray(residuals, f);

                /*
                 * if ( LinalgSymmetricSolver(Keff, ddr, f_last_iter, ddr, conj_grad_precision, conj_grad_relative_maxit, symsolver_type) == false ) {
                 *  std :: cerr << "Conjugate gradients did not converge, attempt to restart step" << endl;
                 *  it = maxIt;
                 *  resErr = 1e10;
                 *  break;
                 * }
                 * if ( LinalgSymmetricSolver(Keff, ddf, f - f_last_iter, ddf, conj_grad_precision, conj_grad_relative_maxit, symsolver_type) == false ) {
                 *  std :: cerr << "Conjugate gradients did not converge, attempt to restart step" << endl;
                 *  it = maxIt;
                 *  resErr = 1e10;
                 *  break;
                 * }
                 */
                linalgsolver->solve(ddr, f_last_iter);
                linalgsolver->solve(ddf, f - f_last_iter);

                if ( !restart_now ) {
                    nodes->giveFullDoFArray(ddr, full_ddr);
                    nodes->giveFullDoFArray(ddf, full_ddf);
                    nodes->updateDirrichletBC(full_ddf, idc_time + idc_dt); //give prescribed DoFs, BC are driven by time
                    nodes->updateDirrichletBC(help_idc_f, idc_time); //give prescribed DoFs, BC are driven by time
                    full_ddf = full_ddf - help_idc_f;
                    help_idc_r = trial_r + full_ddr;
                    load_mult = idc->giveMultiplierCorrection(help_idc_r, full_ddf, time);

                    ddr += load_mult * ddf;
                    idc_time += idc_dt * load_mult;

                    load.setZero();
                    nodes->addRHS_nodalLoad(load, idc_time); //add nodal load
                    nodes->updateDirrichletBC(trial_r, idc_time); //give prescribed DoFs
                }
            } else {         //direct controll
                /*
                 * if ( LinalgSymmetricSolver(Keff, ddr, f, ddr, conj_grad_precision, conj_grad_relative_maxit, symsolver_type) == false ) {
                 *  std :: cerr << "Conjugate gradients did not converge, attempt to restart step" << endl;
                 *  it = maxIt;
                 *  resErr = 1e10;
                 *  break;
                 * }
                 */
                linalgsolver->solve(ddr, f);
            }

            //update DoFs
            updateFieldVariables();
            //compute residuals
            computeForcesAtIntegrationTime(false); //to obtain the actual stress, fluxes, ...

            //compute and print errors
            evaluateErrors();
            if ( it == 0 ) {
                disErr = 0;                        //error in displacement change, only from second iteration
            }
            cout << setw(6) << it << setw(15) << resErr;
            if ( it == 0 ) {
                cout << setw(15) << "---";
            } else {
                cout << setw(15) << disErr;
            }
            cout << setw(15) << eneErr << endl;

            // This check works only if flag "-ffast-math" is removed from CMake
            if ( std :: isnan(resErr) || std :: isnan(disErr) || std :: isnan(eneErr) ) {
                std :: cerr << "calculating with NaN in ";
                if ( std :: isnan(resErr) ) {
                    std :: cerr << "\tresidua ";
                }
                if ( std :: isnan(disErr) ) {
                    std :: cerr << "\tdisplacements ";
                }
                if ( std :: isnan(eneErr) ) {
                    std :: cerr << "\tenergies ";
                }
                std :: cerr << endl;
                it = maxIt;
                resErr = 1e10;
                break;
            }

            it++;
            if ( disErr <= maxDisErr && resErr <= maxResErr && eneErr <= maxEneErr && it >= minIt ) {
                converged = true;
            } else {
                converged = false;
            }
        }

        if ( converged ) {
            this->fully_converged = true;
            computeForcesAtStepEnd(false); //to obtain the actual stress, fluxes, ...
        }

        if ( !converged && dt > dtmin * 1.00001 ) {
            time -= dt;
            dt = fmax(dt * critical_step_decrease, dtmin);
            trial_r = r;
            f_int = f_int_old;
            f_ext = f_ext_old;
            load.setZero(); // std :: fill(begin(load), end(load), 0);
            ddr.setZero(); // std :: fill(begin(ddr), end(ddr), 0);            //ddr *= 0;
            elems->resetMaterialStatuses();   ///> reset material internal vars to the last converged state
            if ( idc ) { //idc solver needs residuals from last converged step
                idc_time = idc_time_converged;
                nodes->addRHS_nodalLoad(load, idc_time); //add nodal load
                nodes->updateDirrichletBC(trial_r, idc_time); //give prescribed DoFs
                updateFieldVariables();      //with ddr=0
                computeForcesAtIntegrationTime(true);
            }

            time += dt;
            std :: cout << "Restarting step, timestep = " << dt << ", time = " << time << endl;
            restarts++;
            restarted = true;
            restart_now = false;
        } else if ( !converged ) {
            if ( disErr < limitDisErr && resErr < limitResErr && eneErr < limitEneErr ) {
                std :: cout << "tolerance increased in this step" << '\n';
                converged = true;
                this->fully_converged = false;
                computeForcesAtStepEnd(false); //to obtain the actual stress, fluxes, ...
            } else {
                std :: cerr << "Error: " << name << " did not converge to the solution" << endl;
                terminated = true;
                return;
            }
        } else if ( ( !restarted ) && converged && it < enlargeIt ) {
            dt = fmin(dt * step_increase, dtmax);
            std :: cout << "enlarging step, timestep = " << dt << '\n';
        } else if ( converged && it > shortenIt && dt > dtmin ) {
            dt = fmax(dt * step_decrease, dtmin);
            std :: cout << "shortening step, timestep = " << dt << '\n';
        }
        if  ( dt > dtmax ) {
            dt = dtmax;
            std :: cout << "shortening step to the maximum one: " << dt << '\n';
        }
        if  ( dt < dtmin ) {
            dt = dtmin;
            std :: cout << "enlarging step to the minimum one: " << dt << '\n';
        }
    }
}


//////////////////////////////////////////////////////////
void SteadyStateNonLinearSolver :: runBeforeEachStep() {
    if ( time == init_time && idc ) {
        nodes->addRHS_nodalLoad(load, time); //add nodal load
        nodes->updateDirrichletBC(trial_r, time); //give prescribed DoFs at time 0
        computeInternalExternalForces(trial_r, load, true, -1.);
    }


    SteadyStateLinearSolver :: runBeforeEachStep();

    cout <<  scientific; //cout << setprecision(8);
    cout << "----------------------------------------------------" << endl;
    cout << setw(6) << "iter." << setw(15) << "residual" << setw(15) << "displacement" << setw(15) << "energy error" << endl;
    cout << setw(6) << " " << setw(15) << maxResErr << setw(15) << maxDisErr << setw(15) << maxEneErr << endl;
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
void SteadyStateNonLinearSolver :: checkAllVectorsForNaNs() {
    // JK left for testing
    bool trial_r_nan, r_nan, full_ddr_nan, ddr_nan, f_int_nan, f_ext_nan, load_nan, f_nan;
    trial_r_nan = r_nan = full_ddr_nan = ddr_nan = f_int_nan = f_ext_nan = load_nan = f_nan = false;
    for ( unsigned i = 0; i < trial_r.size(); i++ ) {
        if ( std :: isnan(trial_r [ i ]) || std :: isinf(trial_r [ i ]) ) {
            trial_r_nan = true;
        }
        if ( std :: isnan(r [ i ]) || std :: isinf(r [ i ]) ) {
            r_nan = true;
        }
        if ( std :: isnan(full_ddr [ i ]) || std :: isinf(full_ddr [ i ]) ) {
            full_ddr_nan = true;
        }
        if ( i < ddr.size() ) {
            if ( std :: isnan(ddr [ i ]) || std :: isinf(ddr [ i ]) ) {
                ddr_nan = true;
            }
        }

        if ( std :: isnan(f_int [ i ]) || std :: isinf(f_int [ i ]) ) {
            f_int_nan = true;
        }
        if ( std :: isnan(f_ext [ i ]) || std :: isinf(f_ext [ i ]) ) {
            f_ext_nan = true;
        }
        if ( std :: isnan(load [ i ]) || std :: isinf(load [ i ]) ) {
            load_nan = true;
        }
        if ( i < f.size() ) {
            if ( std :: isnan(f [ i ]) || std :: isinf(f [ i ]) ) {
                f_nan = true;
            }
        }
    }
    std :: cout << "trial_r_nan " << trial_r_nan << " r_nan " << r_nan << " full_ddr_nan " << full_ddr_nan << " ddr_nan " << ddr_nan << " f_int_nan " << f_int_nan << " f_ext_nan " << f_ext_nan << " load_nan " << load_nan << " f_nan " << f_nan << '\n';
}

//////////////////////////////////////////////////////////
void SteadyStateNonLinearSolver :: runAfterEachStep() {
    if ( !terminated ) {
        SteadyStateLinearSolver :: runAfterEachStep();

        cout << "----------------------------------------------------" << endl;

        if ( idc ) {
            idc_time_converged = idc_time;
        }
    }
}

//////////////////////////////////////////////////////////
double SteadyStateNonLinearSolver :: giveIDCtime(const bool converged) {
    if ( converged ) {
        return this->idc_time_converged;
    }
    return this->idc_time;
}

//////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////
// TRANSIENT LINEAR TRANSPORT SOLVER
TransientLinearTransportSolver :: TransientLinearTransportSolver() {
    name = "TransientLinearTransportSolver";
    timeIntM = 0; //generalized-alpha
    isTimeReal = true;
    applySpectralRadius(0.8);
}

//////////////////////////////////////////////////////////
void TransientLinearTransportSolver :: applySpectralRadius(double rhoinfty) {
    //set up the generalized-alpha method according to Chung and Hulbert 1993
    if ( abs(rhoinfty - 0.5) > 0.5 ) {
       cerr << "Error in solver: spectral radius must be inside interval [0,1]" << endl;
       exit(1);
    }

    //according to Chung and Hulbert, 1993, JAM
    alpha_m = ( 2. * rhoinfty - 1. ) / ( 1. + rhoinfty );
    alpha_f = rhoinfty / ( 1. + rhoinfty );
    gamma = 1. / 2. - alpha_m + alpha_f;
}

//////////////////////////////////////////////////////////
void TransientLinearTransportSolver :: checkIntegrationParams(){
    if (timeIntM==0){   //generalized alpha
        //TODO
    } else if (timeIntM==1){    //HHT method
        //TODO
    } else if (timeIntM==2){    //Newmark method
        if (alpha_m!=0 || alpha_f != 0) {
            cerr << "Solver Error: Newmark method requires alpha_m=alpha_f = 0"<< endl;
            exit(1);
        }        
        if (abs(gamma-0.5)>0.5) {
            cerr << "Solver Error: Newmark method requires gamma withn interva 0-1"<< endl;
            exit(1);
        }    
        if (abs(beta-0.25)>0.25) {
            cerr << "Solver Error: Newmark method requires gamma withn interva 0-0.5"<< endl;
            exit(1);
        }   
    } else {
        cerr << "Solver Error: unknowns method for time integration" << endl;
        exit(1);
    }
}
//////////////////////////////////////////////////////////
void TransientLinearTransportSolver :: setDefaultIntegrationParams(){
    if (timeIntM==0){   //generalized alpha
        applySpectralRadius(0.8);
    } else if (timeIntM==1){    //HHT method
        alpha_m = 0;
        alpha_f = 0;        
        gamma = 0.5;
    } else if (timeIntM==2){    //Newmark method
        alpha_m = alpha_f = 0;        
        gamma = 0.5;
    } else {
        cerr << "Solver Error: unknowns method for time integration" << endl;
        exit(1);
    }
}

//////////////////////////////////////////////////////////
TransientLinearTransportSolver :: ~TransientLinearTransportSolver() {}

//////////////////////////////////////////////////////////
void TransientLinearTransportSolver :: prepareSystemMatricesAndInitialField(string init_r_file, string init_v_file, const bool initial) {
    SteadyStateLinearSolver :: prepareSystemMatricesAndInitialField(init_r_file, init_v_file, initial);

    v_old = Vector :: Zero(totalDoFnum);
    elems->prepareDampingMatrix(C);
    elems->updateDampingMatrix(C);
}

//////////////////////////////////////////////////////////
void TransientLinearTransportSolver :: init(string init_r_file, string init_v_file, const bool initial) {

    checkIntegrationParams();

    prepareSystemMatricesAndInitialField(init_r_file, init_v_file, initial);
    computeKeff();

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
    v = Vector :: Zero(totalDoFnum);
    nodes->giveFullDoFArray(ddr, v);
}

//////////////////////////////////////////////////////////
void TransientLinearTransportSolver :: rebuild() {
    Solver :: rebuild();
    prepareSystemMatricesAndInitialField("", "", false);
    computeKeff();
}

//////////////////////////////////////////////////////////
void TransientLinearTransportSolver :: computeForcesAtIntegrationTime(const bool frozen) {
    elems->integrateDampingForces(v * ( 1. - alpha_m ) +  v_old * alpha_m, f_dam);
    Vector ll = load_old * alpha_f + load * ( 1. - alpha_f ); 
    computeInternalExternalForces( r * alpha_f + trial_r * ( 1. - alpha_f ), ll, frozen, dt * ( 1. - alpha_f ) );
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
    factorizeLinearSystem();
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
    ifstream inputfile( filename.c_str() );
    if ( inputfile.is_open() ) {
        while ( getline(inputfile >> std :: ws, line) ) {
            if ( line.empty() || ( line.at(0) == '#' ) ) {
                continue;
            }
            istringstream iss(line);
            iss >> param;
            if ( param.compare("integration_method") == 0 ) {
                iss >> param;   
                if ( param.compare("newmark") == 0 ) {
                    timeIntM = 2;
                } else if ( param.compare("hht") == 0 ) {
                    timeIntM = 1;
                } else if ( param.compare("generalized_alpha") == 0 ) {
                    timeIntM = 0;
                } else {
                    cerr << "Error: unknown integration method " << param << endl;
                    exit(1);
                }
                setDefaultIntegrationParams();
            } else if ( param.compare("spectral_radius") == 0 ) {
                iss >> num;
                applySpectralRadius(num);
            } else if ( param.compare("gamma") == 0 ) {
                iss >> gamma;
            } else if ( param.compare("beta") == 0 ) {
                iss >> beta;
            } else if ( param.compare("alpha_m") == 0 ) {
                iss >> alpha_m;
            } else if ( param.compare("alpha_f") == 0 ) {
                iss >> alpha_f;
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
    if ( dampingMatrixUpdate == 0 || ( dampingMatrixUpdate > 0 && iteration % abs(dampingMatrixUpdate) == 0 ) ) {
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
        v = Vector :: Zero(totalDoFnum);
    }
    v_old = Vector :: Zero(totalDoFnum);
    a_old = Vector :: Zero(totalDoFnum);
    elems->prepareMassMatrix(M);
    elems->updateMassMatrix(M);

    TransientLinearTransportSolver :: prepareSystemMatricesAndInitialField(init_r_file, init_v_file, initial);
}

//////////////////////////////////////////////////////////
void TransientLinearMechanicalSolver :: init(string init_r_file, string init_v_file, const bool initial) {
    Solver :: init(init_r_file, init_v_file, initial);

    checkIntegrationParams();

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
    v = Vector :: Zero(totalDoFnum);
    Vector v_red = Vector :: Zero(freeDoFnum);
    nodes->giveReducedDoFArray(v, v_red);
    terminated = !LinalgSymmetricSolver(Mred, ddr, f - Cred * v_red,  ddr, conj_grad_precision, conj_grad_relative_maxit, symsolver_type);
    //terminated = !linalgsolver->solve(ddr, f_last_iter);
    a = Vector :: Zero(totalDoFnum);
    nodes->giveFullDoFArray(ddr, a);
}

//////////////////////////////////////////////////////////
void TransientLinearMechanicalSolver :: applySpectralRadius(double rhoinfty) {
    //set up the generalized-alpha method according to Chung and Hulbert 1993
    if ( abs(rhoinfty - 0.5) > 0.5 ) {
       cerr << "Error in solver: spectral radius must be inside interval [0,1]" << endl;
       exit(1);
    }
    //according to Chung and Hulbert, 1993, JAM
    alpha_m = ( 2. * rhoinfty - 1. ) / ( 1. + rhoinfty );
    alpha_f = rhoinfty / ( 1. + rhoinfty );
    gamma = 1. / 2. - alpha_m + alpha_f;
    beta = 1. / 4. * pow(1 - alpha_m + alpha_f, 2);
}

//////////////////////////////////////////////////////////
void TransientLinearMechanicalSolver :: checkIntegrationParams(){
    if (timeIntM==0){   //generalized alpha
        if (alpha_m>0.5 || alpha_f>0.5 || alpha_m<0 || alpha_f < 0 ) {
            cerr << "Solver Error: Generalized-alpha method requires alpha_m<=0.5 and alpha_f<=0.5"<< endl;
            exit(1);
        }    
        if (abs(gamma-0.5)>0.5) {
            cerr << "Solver Error: Generalized-alpha method requires gamma withn interva 0-1"<< endl;
            exit(1);
        }    
        if (beta>0.5 || beta<0.25+0.5*(alpha_f-alpha_m)) {
            cerr << "Solver Error: Generalized-alpha method requires beta withn interva 0.25+0.5(alpha_f-alpha_m) " << beta<< endl;
            exit(1);
        }             

    } else if (timeIntM==1){    //HHT method
        if (alpha_m!=0) {
            cerr << "Solver Error: HHT method requires alpha_m=0"<< endl;
            exit(1);
        }  
        if ( alpha_f>1./3. || alpha_f < 0 ) {
            cerr << "Solver Error: HHT method requires alpha_f inside interva 0-1./3"<< endl;
            exit(1);
        }    
        if (abs(gamma-0.5)>0.5) {
            cerr << "Solver Error: HHT method requires gamma withn interva 0-1"<< endl;
            exit(1);
        }    
        if (beta>0.5 || beta<0.25+0.5*(alpha_f)) {
            cerr << "Solver Error: HHT method requires beta withn interva 0.25+0.5*alpha_f"<< endl;
            exit(1);
        }     
    } else if (timeIntM==2){    //Newmark method
        if (alpha_m!=0 || alpha_f != 0) {
            cerr << "Solver Error: Newmark method requires alpha_m=alpha_f = 0 " << alpha_m << " " << alpha_f<< endl;
            exit(1);
        }        
        if (abs(gamma-0.5)>0.5) {
            cerr << "Solver Error: Newmark method requires gamma withn interva 0-1"<< endl;
            exit(1);
        }    
        if (abs(beta-0.25)>0.25) {
            cerr << "Solver Error: Newmark method requires gamma withn interva 0-0.5"<< endl;
            exit(1);
        }   
    } else {
        cerr << "Solver Error: unknowns method for time integration" << endl;
        exit(1);
    }
}
//////////////////////////////////////////////////////////
void TransientLinearMechanicalSolver :: setDefaultIntegrationParams(){
    if (timeIntM==0){   //generalized alpha
        applySpectralRadius(0.8);
    } else if (timeIntM==1){    //HHT method
        alpha_m = 0.;
        alpha_f = 0.2;        
        gamma = 0.5+alpha_f;
        beta = 0.25*pow(1+alpha_f,2);
    } else if (timeIntM==2){    //Newmark method
        alpha_m = 0.;
        alpha_f = 0.;        
        gamma = 0.5;
        beta = 0.25;
    } else {
        cerr << "Solver Error: unknowns method for time integration" << endl;
        exit(1);
    }
}

//////////////////////////////////////////////////////////
void TransientLinearMechanicalSolver :: computeKeff() {
    Keff = K * ( 1. - alpha_f ) + C * ( ( 1 - alpha_f ) * gamma / dt / beta ) + M * ( ( 1 - alpha_m ) / dt / dt / beta );
    if ( nodes->giveConstraints()->isActive() ) {
        nodes->giveConstraints()->transformToConstraintSpace(Keff);
    }
    factorizeLinearSystem();
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
    Vector ll = load_old * alpha_f + load * ( 1. - alpha_f ); 
    computeInternalExternalForces( r * alpha_f + trial_r * ( 1. - alpha_f ), ll, frozen, dt * ( 1. - alpha_f ) );
    residuals -= f_dam + f_acc;
}

//////////////////////////////////////////////////////////
void TransientLinearMechanicalSolver :: computeForcesAtStepEnd(const bool frozen) {
    elems->integrateDampingForces(v, f_dam);
    elems->integrateInertiaForces(a, f_acc);
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
    if ( massMatrixUpdate == 0 || ( massMatrixUpdate > 0 && iteration % abs(massMatrixUpdate) == 0 ) ) {
        elems->updateMassMatrix(M);
        updated1 = true;
    }
    return ( updated0 || updated1 );
}

//////////////////////////////////////////////////////////
void TransientLinearMechanicalSolver :: computeTotalKineticEnergy(){
    Vector vel = (trial_r - r)/dt;
    W_kin [ 0 ] = (M*vel).dot(vel)/2.;
    W_kin [ 1 ] = (M*v).dot(v)/2.; //TODO: not stable, shall be revised
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
