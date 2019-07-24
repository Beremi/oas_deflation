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
    } else if ( param.compare("SteadyStateNonLinearSolver") == 0 )    {
        SteadyStateNonLinearSolver *newsolver = new SteadyStateNonLinearSolver();
        newsolver->readFromLine(iss);
        return newsolver;
    } else  {
        cerr << "Error: Solver " << param << " is not implemented" << endl;
        exit(0);
    };
}

//////////////////////////////////////////////////////////
void Solver :: runBeforeEachStep() {
    step += 1;
    time += dt;
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

    elems->prepareSteadyStateMatrices(K);
    elems->updateSteadyStateMatrices(K, "elastic");
    f_ext = Vector(totalDoFnum);
    load = Vector(totalDoFnum);
    f_int = Vector(totalDoFnum);
    r = Vector(totalDoFnum);
    pbc = Vector(fixedDoFnum);
    f = Vector(freeDoFnum);
    ddr = Vector(freeDoFnum);
    full_ddr = Vector(totalDoFnum);
}

//////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////
// STEADY STATE LINEAR SOLVER
SteadyStateLinearSolver :: SteadyStateLinearSolver() {
    name = "SteadyStateLinearSolver";
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
        } else if ( param.compare("total_time") == 0 )    {
            bttime = true;
            iss >> termination_time;
        }
    }
    if ( !bdt ) {
        cerr << name << ": solver parameter 'time_step' was not specified" << endl;
        exit(0);
    }
    ;
    if ( !bttime ) {
        cerr << name << ": solver parameter 'total_time' was not specified" << endl;
        exit(0);
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
    if ( ConjGrad(K, ddr, f, ddr) == false ) {
        cerr << "Conjugate gradients did not converge" << endl;
    }
    nodes->giveFullDoFArray(ddr, full_ddr);

    for ( unsigned i = 0; i < totalDoFnum; i++ ) {
        r [ i ] += full_ddr [ i ];
    }
    computeInternalExternalForces(r);
}

//////////////////////////////////////////////////////////
void SteadyStateLinearSolver :: computeInternalExternalForces(Vector &rr) {
    elems->giveInternalForces(rr, f_int);
    nodes->updateExteranlForcesByReactions(f_int, load, f_ext);     //give prescribed DoFs
}

//////////////////////////////////////////////////////////
void SteadyStateLinearSolver :: runBeforeEachStep() {
    Solver :: runBeforeEachStep();
    cout << "######### Solving step " << step << " at time " << time << "; time step " << dt << " #########" << endl;
    load = load * 0.;  //clear nodal load
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
    disErr = resErr = eneErr = 1e-2;
}

//////////////////////////////////////////////////////////
Solver *SteadyStateNonLinearSolver :: readFromLine(istringstream &iss) {
    string param;
    bool bdt, bdtmax, bdtmin, bttime;
    bdt = bttime = bdtmax = bdtmin = false;

    while ( !iss.eof() ) {
        iss >> param;
        if ( param.compare("time_step") == 0 ) {
            bdt = true;
            iss >> dt;
        } else if ( param.compare("max_time_step") == 0 ) {
            bdtmax = true;
            iss >> dtmax;
        } else if ( param.compare("min_time_step") == 0 ) {
            bdtmin = true;
            iss >> dtmin;
        } else if ( param.compare("total_time") == 0 )    {
            bttime = true;
            iss >> termination_time;
        }
    }
    if ( !bdt ) {
        cerr << name << ": solver parameter 'time_step' was not specified" << endl;
        exit(0);
    }
    ;
    if ( !bttime ) {
      cerr << name << ": solver parameter 'total_time' was not specified" << endl;
      exit(0);
    }
    ;
    if ( !bdtmax ) {
        // cout << name << ": solver parameter 'max_time_step' was not specified, setting to timestep" << endl;
        dtmax = dt;
    }
    ;
    if ( !bdtmin ) {
        // cout << name << ": solver parameter 'min_time_step' was not specified, setting to timestep" << endl;
        dtmin = dt;
    }
    ;
    cout << name << " succesfully loaded, ";
    ;
    if ( !bdtmin && !bdtmax){
      cout << "fixed time step used" << endl;
    } else if ( bdtmin || bdtmax){
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
    while ( !converged ){
      //setup loading
      nodes->addRHS_nodalLoad(load, time);  //add nodal load
      nodes->updateDirrichletBC(trial_r, time); //give prescribed DoFs
      computeInternalExternalForces(trial_r);

      unsigned it = 0;
      unsigned maxIt = 30;
      while ( !converged && it < maxIt ) {
          elems->updateSteadyStateMatrices(K, "secant");

          //solve linear system
          nodes->giveReducedDoFArray(f_ext - f_int, f);
          if ( ConjGrad(K, ddr, f, ddr) == false ) {
              cerr << "Conjugate gradients did not converge" << endl;
          }
          nodes->giveFullDoFArray(ddr, full_ddr);

          //update DoFs
          for ( unsigned i = 0; i < totalDoFnum; i++ ) {
              trial_r [ i ] += full_ddr [ i ];
          }

          //compute internal forces
          computeInternalExternalForces(trial_r);

          //compute residuals
          W_int = W_int_old;
          W_ext = W_ext_old;
          for ( unsigned i = 0; i < totalDoFnum; i++ ) {
              residual [ i ] = f_int [ i ] - f_ext [ i ];
              W_int += abs(0.5 * ( f_int [ i ] + f_int_old [ i ] ) * ( r [ i ] - trial_r [ i ] ) );
              W_ext += abs(0.5 * ( f_ext [ i ] + f_ext_old [ i ] ) * ( r [ i ] - trial_r [ i ] ) );
          }

          //compute errors
          double residu_error =  l2_norm(residual) / max(max(l2_norm(f_ext), l2_norm(f_int) ), EPS2);
          double displa_error = ( it == 0 ) ? 0. : l2_norm(full_ddr) / max(l2_norm(trial_r), EPS2);   //error in displacement change, only from second iteration
          double energy_error =  abs(inner_product(& residual [ 0 ], & residual [ totalDoFnum ], & full_ddr [ 0 ], ( double ) ( 0 ) ) ) / max(max(W_ext, W_int), EPS2);

          cout << setw(6) << it << setw(15) << residu_error;
          if ( it == 0 ) {
              cout << setw(15) << "---";
          } else {
              cout << setw(15) << displa_error;
          }
          cout << setw(15) << energy_error;
          cout << endl;

          if (std :: isnan(residu_error) || std :: isnan(displa_error) || std :: isnan(energy_error) ){
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
          dt = fmax(dt/2, dtmin);
          time += dt;
          cout << "Restarting step, timestep = " << dt << ", time = " << time << endl;
      } else if ( !converged ) {
          cerr << "Error: Nonlinear static solver did not converge to the solution" << endl;
          exit(1);
      } else if ( converged && it < maxIt/3){
          dt = fmin(dt*2, dtmax);
          std::cout << "enlarging step, timestep = " << dt << '\n';
      }
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
//////////////////////////////////////////////////////////
void SteadyStateNonLinearSolver :: runAfterEachStep() {
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
