#include "solver_explicit.h"
#include "adaptivity.h"
#define numPhysicalFields 4

using namespace std;


//////////////////////////////////////////////////////////
TransientCentralDifferenceMechanicalSolver :: TransientCentralDifferenceMechanicalSolver(){

}

//////////////////////////////////////////////////////////
TransientCentralDifferenceMechanicalSolver :: ~TransientCentralDifferenceMechanicalSolver(){
}

//////////////////////////////////////////////////////////
void TransientCentralDifferenceMechanicalSolver :: init(std :: string init_r_file, std :: string init_v_file, const bool initial){
    Solver :: init(init_r_file, init_v_file, initial);

    v_old = Vector :: Zero(totalDoFnum);
    v = Vector :: Zero(totalDoFnum);
    accel = Vector :: Zero(totalDoFnum);  

    nodes->updateDirrichletBC(trial_r, time); //give prescribed DoFs
    nodes->addRHS_nodalLoad(load, time); //add nodal load      
}

//////////////////////////////////////////////////////////
void TransientCentralDifferenceMechanicalSolver :: solve(){
    //according to Belytschko, Liu, Moran, page 71

    /*
    elems->integrateDampingForces(v, f_dam);
    elems->integrateInertiaForces(a, f_acc);
    computeInternalExternalForces(r);
    residuals -= f_dam + f_acc;
    
    //computeAcceleration();

    v_old = v;
    v = dt*accel + v_old;
    ddr = dt*v;

    nodes->updateDirrichletBC(trial_r, time); //give prescribed DoFs
    nodes->addRHS_nodalLoad(load, time); //add nodal load

    //update DoFs
    updateFieldVariables();
    //compute residuals
    computeInternalExternalForces(trial_r, load, frozen, dt);

    computeForcesAtStepEnd(false); //to obtain the actual stress, fluxes, ...   
    */
}

//////////////////////////////////////////////////////////
void TransientCentralDifferenceMechanicalSolver :: giveValues(std :: string code, Vector &result) const{
    if ( code.compare("???????") == 0 ) {
        result.resize(1);
        result [ 0 ] = 0;
    } else {
        Solver :: giveValues(code, result);
    }
}

