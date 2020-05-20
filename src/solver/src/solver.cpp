#include "solver.h"
#define EPS2 1e-30

//////////////////////////////////////////////////////////
Solver *Solver :: readFromLine(istringstream &iss) {
    string param;
    iss >> param;
    if ( param.compare("SteadyStateLinearSolver") == 0 ) {
        SteadyStateLinearSolver *newsolver = new SteadyStateLinearSolver();
        newsolver->readFromLine(iss);
        return newsolver;
    } else if ( param.compare("SteadyStateNonLinearSolver") == 0 ) {
        SteadyStateNonLinearSolver *newsolver = new SteadyStateNonLinearSolver();
        newsolver->readFromLine(iss);
        return newsolver;
    } else if ( param.compare("TransientLinearMechanicalSolver") == 0 ) {
        TransientLinearMechanicalSolver *newsolver = new TransientLinearMechanicalSolver();
        newsolver->readFromLine(iss);
        return newsolver;
    } else if ( param.compare("TransientLinearTransportSolver") == 0 ) {
        TransientLinearTransportSolver *newsolver = new TransientLinearTransportSolver();
        newsolver->readFromLine(iss);
        return newsolver;
    } else {
        cerr << "Error: Solver " << param << " is not implemented" << endl;
        exit(EXIT_FAILURE);
    };
}

//////////////////////////////////////////////////////////
void Solver :: setNextStepTime() {
    double nextExtremeTime = funcs->giveTimeOfNextExtreme(time);
    // NOTE 1/4 of time step added to prevent next step extremely short
    if ( nextExtremeTime < time + 1.25 * dt ) {
        time = nextExtremeTime;
    } else {
        time += dt;
    }
}

//////////////////////////////////////////////////////////
void Solver :: runBeforeEachStep() {
    step += 1;
    setNextStepTime();
}

//////////////////////////////////////////////////////////
void Solver :: runAfterEachStep() {
    if ( time + dt > termination_time ) {
        dt = termination_time - time;
    }
    if ( dt < 1e-15 ) {
        terminated = true;
    }
}

//////////////////////////////////////////////////////////
void Solver :: init() {
    step = 0;
    time = 0.;

    terminated = false;

    freeDoFnum = nodes->giveNumFreeDoFs();
    fixedDoFnum = ( nodes->giveTotalNumDoFs() - freeDoFnum );
    totalDoFnum = freeDoFnum + fixedDoFnum;

    elems->prepareSteadyStateMatrix(Kini);
    K = Kini;
    elems->updateSteadyStateMatrix(K, "elastic");
    f_ext = Vector(totalDoFnum);
    load = Vector(totalDoFnum);
    f_int = Vector(totalDoFnum);
    r = Vector(totalDoFnum);
    pbc = Vector(fixedDoFnum);
    f = Vector(freeDoFnum - nodes->giveNumConstrDoFs() );
    ddr = Vector(freeDoFnum - nodes->giveNumConstrDoFs() );
    full_ddr = Vector(totalDoFnum);
}

//////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////
// STEADY STATE LINEAR SOLVER
SteadyStateLinearSolver :: SteadyStateLinearSolver() {
    name = "SteadyStateLinearSolver";
    conj_grad_precission = 1e-16;
    conj_grad_relative_maxit = 0.85;
}

//////////////////////////////////////////////////////////
SteadyStateLinearSolver :: ~SteadyStateLinearSolver() {}

//////////////////////////////////////////////////////////
void SteadyStateLinearSolver :: init() {
    Solver :: init();
}

//////////////////////////////////////////////////////////
Solver *SteadyStateLinearSolver :: readFromLine(istringstream &iss) {
    string param;
    bool bdt, bttime;
    bdt = bttime = false;

    while ( !iss.eof() ) {
        iss >> param;
        if ( param.compare("time_step") == 0 ) {
            bdt = true;
            iss >> dt;
        } else if ( param.compare("total_time") == 0 ) {
            bttime = true;
            iss >> termination_time;
        } else if ( param.compare("conj_grad_precission") == 0 ) {
            iss >> conj_grad_precission;
        } else if ( param.compare("conj_grad_relative_maxit") == 0 ) {
            iss >> conj_grad_relative_maxit;
        }
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
    ;
    cout << name << " succesfully loaded" << endl;
    return this;
};

//////////////////////////////////////////////////////////
void SteadyStateLinearSolver :: solve() {
    nodes->addRHS_nodalLoad(load, time);  //add nodal load
    nodes->updateDirrichletBC(r, time); //give prescribed DoFs
    computeInternalExternalForces(r);

    //solve linear system
    nodes->giveReducedDoFArray(f_ext - f_int, f);
    terminated = !ConjGrad(K, ddr, f, ddr, conj_grad_precission, conj_grad_relative_maxit);
    // if ( ConjGrad(K, ddr, f, ddr, conj_grad_precission, conj_grad_relative_maxit) == false ) {
    //     terminated = true;
    //     cerr << "Conjugate gradients did not converge" << endl;
    // }
    nodes->giveFullDoFArray(ddr, full_ddr);

    for ( unsigned i = 0; i < totalDoFnum; i++ ) {
        r [ i ] += full_ddr [ i ];
    }

    /*
     * for ( unsigned i = 0; i < freeDoFnum; i++ ) {
     *  for ( unsigned j = 0; j < freeDoFnum; j++ ) cout << K [ i ][j] << "\t";
     *  cout << "|\t" << ddr[i] << "\t=\t" << f[i] << endl;
     * }
     * exit(1);
     */

    computeInternalExternalForces(r);
}

//////////////////////////////////////////////////////////
void SteadyStateLinearSolver :: computeInternalExternalForces(Vector &rr) {
    elems->giveInternalForces(rr, f_int);
    nodes->updateExteranlForcesByReactions(f_int, load, f_ext);     //give prescribed DoFs
    // #constr_old
    // if (nodes->giveConstraints()->isActive()){
    //   // also ddr must be updated, because displacement error is calculated from relative change (||ddr||/||r||)
    //   // matters especially when force control on master is applied
    //   nodes->giveConstraints()->calculateDependentDoFs(full_ddr);
    // }
}

//////////////////////////////////////////////////////////
void SteadyStateLinearSolver :: runBeforeEachStep() {
    Solver :: runBeforeEachStep();
    cout << "######### Solving step " << step << " at time " << time << "; time step " << dt << " #########" << endl;
    load *= 0;  //clear nodal load
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
}

//////////////////////////////////////////////////////////
SteadyStateNonLinearSolver :: ~SteadyStateNonLinearSolver() {}

//////////////////////////////////////////////////////////
void SteadyStateNonLinearSolver :: init() {
    Solver :: init();
    f_int_old = Vector(totalDoFnum);
    f_ext_old = Vector(totalDoFnum);
    trial_r = Vector(totalDoFnum);
    residual = Vector(totalDoFnum);
    W_ext_old = 0;
    W_int_old = 0;
}

//////////////////////////////////////////////////////////
Solver *SteadyStateNonLinearSolver :: readFromLine(istringstream &iss) {
    SteadyStateLinearSolver :: readFromLine(iss);
    iss.clear(); // clear string stream
    iss.seekg(0, iss.beg); //reset position in string stream

    maxIt = 30;
    disErr = resErr = eneErr = 1e-5;
    limitEneErr = limitResErr = limitDisErr = 0;
    step_increase = 1.25;
    step_decrease = 0.8;
    critical_step_decrease = 0.5;

    string param;
    dtmax = dtmin = dt;
    bool bdtmin = false;
    bool bdtmax = false;
    double valueIN;

    while ( !iss.eof() ) {
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
        } else if ( param.compare("limit_tolerance") == 0 ) {
            iss >> valueIN;
            limitEneErr = limitResErr = limitDisErr = valueIN;
        } else if ( param.compare("maxIt") == 0 ) {
            iss >> maxIt;
            if ( maxIt < 1 ) {
              std :: cerr << "number of itteration cannot be smaller than " << maxIt << "!!!,";
              maxIt = 30;
              std::cout << " setting to default value: maxIt = " << maxIt << '\n';
            } else if ( maxIt < 3 ) {
              std :: cout << "solver parameter maxIt set to " << maxIt << ", be carefull with such a small number" << '\n';
            }
        } else if ( param.compare("step_increase") == 0 ) {
            iss >> valueIN;
            if (valueIN < 1){
              std::cerr << "step_increase cannot be smaller than 1! leaving default value " << step_increase << '\n';
            } else {
              step_increase = valueIN;
            }
        } else if ( param.compare("step_decrease") == 0 ) {
            iss >> valueIN;
            if (valueIN > 1){
              std::cerr << "step_decrease cannot be greater than 1! leaving default value " << step_decrease << '\n';
            } else {
              step_decrease = valueIN;
            }
        } else if ( param.compare("critical_step_decrease") == 0 ) {
          iss >> valueIN;
          if (valueIN > 1){
            std::cerr << "critical_step_decrease cannot be greater than 1! leaving default value " << critical_step_decrease << '\n';
          } else {
            critical_step_decrease = valueIN;
          }
        }
    }
    cout << name << " succesfully loaded, ";
    ;
    if ( !bdtmin && !bdtmax ) {
        cout << "fixed time step used" << endl;
    } else if ( bdtmin || bdtmax ) {
        cout << "adaptive time step used" << endl;
    } else {
        cout << endl;
    }
    ;
    return this;
};

//////////////////////////////////////////////////////////
void SteadyStateNonLinearSolver :: solve() {
    bool converged = false;
    bool restarted = false;
    double displa_error = 0;
    double energy_error = 0;
    double residu_error = 0;
    while ( !converged ) {
        //setup loading
        nodes->addRHS_nodalLoad(load, time); //add nodal load
        nodes->updateDirrichletBC(trial_r, time); //give prescribed DoFs
        computeInternalExternalForces(trial_r);

        // std::cout << " ---------------------- initial " << '\n';
        // this->printAllVectors();

        unsigned it = 0;
        // unsigned maxIt = 30;
        while ( !converged && it < maxIt ) {
            K = Kini;
            elems->updateSteadyStateMatrix(K, "secant");

            //solve linear system
            nodes->giveReducedDoFArray(f_ext - f_int, f);
            // terminated = !ConjGrad(K, ddr, f, ddr, conj_grad_precission, conj_grad_relative_maxit);
            if ( ConjGrad(K, ddr, f, ddr, conj_grad_precission, conj_grad_relative_maxit) == false ) {
                terminated = true;
                cerr << "Conjugate gradients did not converge" << endl;
                return;
            }
            // sem přidat constraint z elem container
            nodes->giveFullDoFArray(ddr, full_ddr);

            //update DoFs
            for ( unsigned i = 0; i < totalDoFnum; i++ ) {
                trial_r [ i ] += full_ddr [ i ];
            }

            // std::cout << "before ----------------------" << '\n';
            // this->printAllVectors();
            //compute internal forces
            computeInternalExternalForces(trial_r);
            // std::cout << "after ----------------------" << '\n';
            // this->printAllVectors();

            //compute residuals
            W_int = W_int_old;
            W_ext = W_ext_old;
            for ( unsigned i = 0; i < totalDoFnum; i++ ) {
                residual [ i ] = f_int [ i ] - f_ext [ i ];
                W_int += abs(0.5 * ( f_int [ i ] + f_int_old [ i ] ) * ( r [ i ] - trial_r [ i ] ) );
                W_ext += abs(0.5 * ( f_ext [ i ] + f_ext_old [ i ] ) * ( r [ i ] - trial_r [ i ] ) );
            }

            //compute errors
            residu_error =  l2_norm(residual) / max(max(l2_norm(f_ext), l2_norm(f_int) ), EPS2);
            displa_error = ( it == 0 ) ? 0. : l2_norm(full_ddr) / max(l2_norm(trial_r), EPS2); //error in displacement change, only from second iteration
            energy_error =  abs(inner_product(& residual [ 0 ], & residual [ totalDoFnum ], & full_ddr [ 0 ], ( double ) ( 0 ) ) ) / max(max(W_ext, W_int), EPS2);

            cout << setw(6) << it << setw(15) << residu_error;
            if ( it == 0 ) {
                cout << setw(15) << "---";
            } else {
                cout << setw(15) << displa_error;
            }
            cout << setw(15) << energy_error;

            cout << endl;

            if ( std :: isnan(residu_error) || std :: isnan(displa_error) || std :: isnan(energy_error) ) {
                cerr << "calculating with NaN - exit" << '\n';
                exit(1);
            }

            if ( displa_error > disErr || residu_error > resErr || energy_error > eneErr ) {
                converged = false;
            } else {
                converged = true;
            }
            it++;
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
        } else if ( ( !restarted ) && converged && it < maxIt / 3 && dt < dtmax ) {
            dt = fmin(dt * step_increase, dtmax);
            std :: cout << "enlarging step, timestep = " << dt << '\n';
        } else if ( converged && it > maxIt / 2 && dt > dtmin ) {
            dt = fmin(dt * step_decrease, dtmax);
            std :: cout << "shortening step, timestep = " << dt << '\n';
        }
        // std::cerr << "number of iterations: " << it << '\n';
    }
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
        for ( unsigned i = 0; i < totalDoFnum; i++ ) {
            r [ i ] = trial_r [ i ];
            f_int_old [ i ] = f_int [ i ];
            f_ext_old [ i ] = f_ext [ i ];
        }
        W_int_old = W_int;
        W_ext_old = W_ext;
        elems->updateMaterialStatuses();
        cout << "----------------------------------------------------" << endl;
    }
}

//////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////
// TRANSIENT LINEAR MECHANICAL SOLVER
TransientLinearMechanicalSolver :: TransientLinearMechanicalSolver() {
    name = "TransientLinearMechanicalSolver";

    applySpectralRadius(0.8);
}

//////////////////////////////////////////////////////////
TransientLinearMechanicalSolver :: ~TransientLinearMechanicalSolver() {}

//////////////////////////////////////////////////////////
void TransientLinearMechanicalSolver :: init() {
    SteadyStateLinearSolver :: init();
    elems->prepareMassMatrix(M);
    elems->updateMassMatrix(M);
    C = M * 0 + K * 0; //no damping, possibly change to something if needed
    updateKeff();

    r = Vector(totalDoFnum); //initial conitions, assumed zero for now
    r_old = r;
    v = Vector(totalDoFnum); //initial conitions, assumed zero for now
    r_red = Vector(freeDoFnum - nodes->giveNumConstrDoFs() );
    nodes->giveReducedDoFArray(r, r_red);
    v_red = Vector(freeDoFnum - nodes->giveNumConstrDoFs() );
    nodes->giveReducedDoFArray(v, v_red);

    //compute initial acceleration
    a = Vector(totalDoFnum);
    nodes->addRHS_nodalLoad(load, 0);
    nodes->updateDirrichletBC(r, 0);
    computeInternalExternalForces(r); //at time 0
    nodes->giveReducedDoFArray(f_ext - f_int, f);
    terminated = !ConjGrad(M, a_red, f - C * v_red,  a_red, conj_grad_precission, conj_grad_relative_maxit);
}

//////////////////////////////////////////////////////////
Solver *TransientLinearMechanicalSolver :: readFromLine(istringstream &iss) {
    string param;

    SteadyStateLinearSolver :: readFromLine(iss);
    iss.clear(); // clear string stream
    iss.seekg(0, iss.beg); //reset position in string stream

    while ( !iss.eof() ) {
        iss >> param;
        if ( param.compare("alpha_f") == 0 ) {
            iss >> alpha_f;
        } else if ( param.compare("alpha_m") == 0 ) {
            iss >> alpha_m;
        } else if ( param.compare("gamma") == 0 ) {
            iss >> gamma;
        } else if ( param.compare("beta") == 0 ) {
            iss >> gamma;
        } else if ( param.compare("spectral_radius") == 0 ) {
            double rhoinfty;
            iss >> rhoinfty;
            applySpectralRadius(rhoinfty);
        }
    }
    if ( alpha_m > 0.5 ) {
        cerr << "Error in solver: alpha_m (" << alpha_m << ") cannot exceed 0.5" << endl;
        exit(1);
    }
    if ( alpha_f > 0.5 ) {
        cerr << "Error in solver: alpha_f (" << alpha_f << ") cannot exceed 0.5" << endl;
        exit(1);
    }
    if ( alpha_m > alpha_f ) {
        cerr << "Error in solver: alpha_f (" << alpha_f << ") must be larger than alpha_m (" << alpha_m << ")" << endl;
        exit(1);
    }
    if ( gamma < 0.5 ) {
        cerr << "Error in solver: gamma (" << gamma << ") must be larger than 0.5" << endl;
        exit(1);
    }
    if ( beta < 0.25 + 0.5 * ( alpha_f - alpha_m ) ) {
        cerr << "Error in solver: beta (" << beta << ") must be larger than 0.25 + 0.5*(alpha_f-alpha_m)" << endl;
        exit(1);
    }
    cout << name << " succesfully loaded, ";
    return this;
};

//////////////////////////////////////////////////////////
void TransientLinearMechanicalSolver :: solve() {
    nodes->addRHS_nodalLoad(load, time - alpha_f * dt);  //add nodal load, at time alpha_f
    nodes->updateDirrichletBC(r, time);   //give prescribed DoFs, at time t
    r_f = r_old * alpha_f + r * ( 1. - alpha_f );
    computeInternalExternalForces(r_f); //at time alpha_f

    //solve linear system
    nodes->giveReducedDoFArray(f_ext - f_int, f);
    updateFeff();
    terminated = !ConjGrad(Keff, ddr, feff, ddr, conj_grad_precission, conj_grad_relative_maxit);

    updateFieldVariables();
    r_old = r;
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
void TransientLinearMechanicalSolver :: updateKeff() {
    Keff = K * ( 1. - alpha_f ) + C * ( ( 1 - alpha_f ) * gamma / dt / beta ) + M * ( ( 1 - alpha_m ) / dt / dt / beta );
}

//////////////////////////////////////////////////////////
void TransientLinearMechanicalSolver :: updateFeff() {
    feff = f - C * ( ( 1. + gamma * ( alpha_f - 1 ) / beta ) * v_red + dt * ( 1. - alpha_f ) * ( 1. - gamma / ( 2. * beta ) ) * a_red ) - M * ( ( alpha_m - 1. ) / ( dt * beta ) * v_red  + ( 1. + ( alpha_m - 1 ) / ( 2. * beta ) ) * a_red );
}

//////////////////////////////////////////////////////////
void TransientLinearMechanicalSolver :: updateFieldVariables() {
    r_red += ddr;
    v_red = ddr * ( gamma / ( dt * beta ) ) + v_red * ( 1. - gamma / beta ) + a_red * ( dt * ( 1. - gamma / ( 2. * beta ) ) );
    a_red = ( ddr - v_red * dt - a_red * ( pow(dt, 2) * ( 1. / 2. - beta ) ) ) / ( pow(dt, 2) * beta );

    nodes->giveFullDoFArray(r_red, r);
    nodes->giveFullDoFArray(v_red, v);
    nodes->giveFullDoFArray(a_red, a);
}

//////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////
// TRANSIENT LINEAR TRANSPORT SOLVER
TransientLinearTransportSolver :: TransientLinearTransportSolver() {
    name = "TransientLinearTransportSolver";

    applySpectralRadius(0.8);
}

//////////////////////////////////////////////////////////
TransientLinearTransportSolver :: ~TransientLinearTransportSolver() {}

//////////////////////////////////////////////////////////
void TransientLinearTransportSolver :: init() {
    SteadyStateLinearSolver :: init();
    elems->prepareCapacityMatrix(C);
    elems->updateCapacityMatrix(C);
    updateKeff();

    r = Vector(totalDoFnum); //initial conitions, assumed zero for now
    r_old = r;
    r_red = Vector(freeDoFnum - nodes->giveNumConstrDoFs() );
    nodes->giveReducedDoFArray(r, r_red);

    //compute initial presure derivative
    v = Vector(totalDoFnum);
    nodes->addRHS_nodalLoad(load, 0);
    nodes->updateDirrichletBC(r, 0);
    computeInternalExternalForces(r); //at time 0
    nodes->giveReducedDoFArray(f_ext - f_int, f);
    terminated = !ConjGrad(C, v_red, f,  v_red, conj_grad_precission, conj_grad_relative_maxit);
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
void TransientLinearTransportSolver :: updateKeff() {
    Keff = C * ( ( 1. - alpha_m ) / ( dt * gamma ) ) + K * ( 1. - alpha_f );
}

//////////////////////////////////////////////////////////
void TransientLinearTransportSolver :: updateFeff() {
    feff = f - ( C * v_red ) * ( 1. + ( alpha_m - 1. ) / gamma );
}

//////////////////////////////////////////////////////////
void TransientLinearTransportSolver :: updateFieldVariables() {
    r_red += ddr;
    v_red = ddr * ( 1. / ( dt * gamma ) ) + v_red * ( 1. - 1. / gamma );

    nodes->giveFullDoFArray(r_red, r);
    nodes->giveFullDoFArray(v_red, v);
}
