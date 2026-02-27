#include "solver_explicit.h"
#include "adaptivity.h"
#define numPhysicalFields 4

using namespace std;


//////////////////////////////////////////////////////////
TransientCentralDifferenceMechanicalSolver :: TransientCentralDifferenceMechanicalSolver() {
    name = "TransientCentralDifferenceMechanicalSolver";
    show_period = 1e3;
    showTime = false;
}

//////////////////////////////////////////////////////////
TransientCentralDifferenceMechanicalSolver :: ~TransientCentralDifferenceMechanicalSolver() {}

//////////////////////////////////////////////////////////
void TransientCentralDifferenceMechanicalSolver :: init(std :: string init_r_file, std :: string init_v_file, const bool initial) {
    Solver :: init(init_r_file, init_v_file, initial);

    v_red_old = Vector :: Zero(freeDoFnum);
    v_red = Vector :: Zero(freeDoFnum);
    a_red = Vector :: Zero(freeDoFnum);

    v = Vector :: Zero(totalDoFnum);
    a = Vector :: Zero(totalDoFnum);


    elems->prepareMassMatrix(M, 1);
    elems->updateMassMatrix(M, 1);
    if ( nodes->giveConstraints()->isActive() ) {
        nodes->giveConstraints()->transformToConstraintSpace(M);
    }
    lumpedM = M.diagonal();

    nodes->updateDirrichletBC(trial_r, time); //give prescribed DoFs
    nodes->addRHS_nodalLoad(load, time); //add nodal load

    //compute actions at the end of the last time step
    //elems->integrateDampingForces(v, f_dam);
    //elems->integrateInertiaForces(a, f_acc);
    computeInternalExternalForces(r, load, 0, time, dt);
    residuals -= f_dam + f_acc;
}

//////////////////////////////////////////////////////////
void TransientCentralDifferenceMechanicalSolver :: computeAcceleration() {
    nodes->giveReducedForceArray(residuals, f);
    for ( unsigned i = 0; i < freeDoFnum; i++ ) {
        a_red [ i ] = f [ i ] / lumpedM [ i ];
    }
}

//////////////////////////////////////////////////////////
void TransientCentralDifferenceMechanicalSolver :: runBeforeEachStep() {
    Solver :: runBeforeEachStep();
    trial_r = r;
    if ( step % show_period == 0 ) {
        cout << "######### Solving step " << step << " at time " << time << "; time step " << dt << " #########" << endl;
        showTime = true;
    } else {
        showTime = false;
    }
}

//////////////////////////////////////////////////////////
Solver *TransientCentralDifferenceMechanicalSolver :: readFromFile(const string filename) {
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
            } else if ( param.compare("show_every") == 0 ) {
                iss >> show_period;
            } else if ( param.compare("init_time") == 0 ) {
                iss >> this->init_time;
            } else if ( param.compare("init_step") == 0 ) {
                iss >> this->init_step;
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
void TransientCentralDifferenceMechanicalSolver :: solve() {
    //according to Belytschko, Liu, Moran, page 71

    computeAcceleration();

    v_red_old = v_red;
    if ( time == dt ) { //the first time step, only half velocity action
        v_red = 0.5 * dt * a_red + v_red_old;
    } else {
        v_red = dt * a_red + v_red_old;
    }
    ddr = dt * v_red;


    //update DoFs
    updateFieldVariables();
    nodes->updateDirrichletBC(trial_r, time); //give prescribed DoFs
    nodes->addRHS_nodalLoad(load, time); //add nodal load

    //update velocity and acceleration
    nodes->giveFullDoFArray(v_red, v);
    nodes->giveFullDoFArray(a_red, a);
    vector< double >blocked = bcs->giveBlockedDoFValues(time);
    double new_v;
    for ( auto &k: blocked ) {
        new_v = ( trial_r [ k ] - r [ k ] ) / dt;
        a [ k ] = ( new_v - v [ k ] ) / dt;
        v [ k ] = new_v;
    }

    //DOES NOT SUPPORT CONSTRAINT WITH CONJUGATE VARIABLES AND FUNCTIONS
    nodes->giveConstraints()->calculateDependentDoFs(v);
    nodes->giveConstraints()->calculateDependentDoFs(a);

    //compute residuals
    elems->integrateDampingForces(v, f_dam);
    computeInternalExternalForces(trial_r, load, 0, time, dt);
    residuals -= f_dam;
}

//////////////////////////////////////////////////////////
void TransientCentralDifferenceMechanicalSolver :: giveValues(std :: string code, Vector &result) const {
    if ( code.compare("???????") == 0 ) {
        result.resize(1);
        result [ 0 ] = 0;
    } else {
        Solver :: giveValues(code, result);
    }
}
