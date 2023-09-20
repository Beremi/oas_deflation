#include "solver.h"
#include "solver_implicit.h"
#include "solver_explicit.h"
#include "adaptivity.h"
#define numPhysicalFields 4

using namespace std;

//////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////
// BASIC SOLVER CLASS
Solver :: Solver() {
    name = "basic solver";
    time = init_time;
}


void Solver :: setContainers(ElementContainer *e, NodeContainer *n, FunctionContainer *functions, BCContainer *bc) {
    elems = e;
    nodes = n;
    funcs = functions;
    bcs = bc;
}

//////////////////////////////////////////////////////////
Solver *Solver :: readFromFile(const string filename) {
    string param, paramA, line;
    ifstream inputfile(filename.c_str() );
    if ( inputfile.is_open() ) {
        while ( getline(inputfile >> std :: ws, line) ) {
            if ( line.empty() || ( line.at(0) == '#' ) ) {
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
    double nextBCTime = bcs->giveTimeOfNextChange(time);
    double nextCritTime = min(nextExtremeTime, nextBCTime);

    if ( abs( time - bcs->giveTimeOfNextChange(time - 1e-12) ) < 1e-12 ) {
        masterModel->jumpToNextStage();
    }

    // NOTE 1/4 of time step added to prevent next step extremely short
    if ( nextCritTime < time + 1.25 * dt ) {
        dt = nextCritTime - time;
        time = nextCritTime;
    } else {
        time += dt;
    }
}

//////////////////////////////////////////////////////////
void Solver :: runBeforeEachStep() {
    step += 1;

    setNextStepTime();

    load_old = load; //copy old load to be used in generalized alpha method
    load.setZero();  //clear nodal load
    ddr.setZero(); //clear step contribution;
}

//////////////////////////////////////////////////////////
void Solver :: runAfterEachStep() {
    r = trial_r;
    f_int_old = f_int;
    f_ext_old = f_ext;

    elems->updateMaterialStatuses();

    if ( time + dt > termination_time ) {
        dt = termination_time - time;
    }
    if ( dt < 1e-12 ) {
        terminated = true;
    }
}

//////////////////////////////////////////////////////////
void Solver :: setTime(double t) {
    time = t;
    if ( time + 1e-12 > termination_time ) {
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

    if ( initial ) {
        step = init_step;
        dt = initdt;
    }

    terminated = false;
    fully_converged = true;  ///> only make sense for nonlin solvers, true otherwise

    freeDoFnum = nodes->giveNumFreeDoFs();
    totalDoFnum = nodes->giveTotalNumDoFs();

    load = Vector :: Zero(totalDoFnum);
    r = Vector :: Zero(totalDoFnum);
    f_ext = Vector :: Zero(totalDoFnum);
    f_int = Vector :: Zero(totalDoFnum);
    f_dam = Vector :: Zero(totalDoFnum);
    f_acc = Vector :: Zero(totalDoFnum);
    pbc = Vector :: Zero(totalDoFnum - freeDoFnum - nodes->giveNumConstrDoFs() );
    f = Vector :: Zero(freeDoFnum);

    trial_r = Vector :: Zero(totalDoFnum);
    residuals = Vector :: Zero(totalDoFnum);
    ddr = Vector :: Zero(freeDoFnum);
    full_ddr = Vector :: Zero(totalDoFnum);
    f_int_old = Vector :: Zero(totalDoFnum);
    f_ext_old = Vector :: Zero(totalDoFnum);
}

//////////////////////////////////////////////////////////
void Solver :: updateFieldVariables() {
    nodes->giveFullDoFArray(ddr, full_ddr);
    nodes->updateFullDoFsByDependenciesOnConjugates(full_ddr, trial_r, f_ext);
    trial_r += full_ddr;
}

//////////////////////////////////////////////////////////
void Solver :: giveValues(string code, Vector &result) const {
    if ( code.compare("simulation_time") == 0 ) {
        result.resize(1);
        result [ 0 ] = time;
    } else if ( code.compare("elapsed_time") == 0 ) {
        result.resize(1);
        result [ 0 ] = chrono :: duration_cast< std :: chrono :: duration< double > >(std :: chrono :: system_clock :: now() - masterModel->giveStartTime() ).count();
    } else if ( code.compare("number_of_dof") == 0 ) {
        result.resize(1);
        result [ 0 ] = freeDoFnum;
    } else {
        result.resize(0);
    }
}

//////////////////////////////////////////////////////////
void Solver :: rebuild() {
    freeDoFnum = nodes->giveNumFreeDoFs();
    pbc = Vector :: Zero(totalDoFnum - freeDoFnum - nodes->giveNumConstrDoFs() );
    f = Vector :: Zero(freeDoFnum);
    ddr = Vector :: Zero(freeDoFnum);
}
