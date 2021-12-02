#ifndef SOLVER_H
#define SOLVER_H

#include <vector>
#include <iostream>
#include <fstream>
#include <iomanip>      // std::setw
#include "linalg.h"
#include "element_container.h"
#include "node_container.h"
#include "indirect_displ_control.h"

//////////////////////////////////////////////////////////
class Solver
{
protected:
    string name;
    ElementContainer *elems;
    NodeContainer *nodes;
    FunctionContainer *funcs;
    double time, dt, initdt, termination_time;
    double init_time = 0.0;  ///> when starting from previously calculated results
    Vector f_ext, load, load_old, f_int, pbc, r, f, full_ddr, ddr, residuals;
    Vector f_int_old, f_ext_old, f_dam, f_acc, trial_r;
    unsigned freeDoFnum, totalDoFnum;
    int step;
    unsigned init_step = 0;  ///> when starting from previously calculated results
    bool terminated, fully_converged;
    string symsolver_type = "EigenConj";

    virtual void updateFieldVariables();

public:
    Solver();
    virtual ~Solver() {};
    virtual void init(string init_r_file, string init_v_file, const bool initial = true);
    virtual Solver *readFromFile(const string filename);
    virtual void solveStep() { runBeforeEachStep(); solve();  runAfterEachStep(); };
    void setContainers(ElementContainer *e, NodeContainer *n, FunctionContainer *functions);
    string giveName() const { return name; }
    bool isTerminated() { return terminated; }
    bool convergedToTolerance() const { return fully_converged; };
    Vector giveDoFValues() const { return r; }
    Vector giveTrialDoFValues() const { return trial_r; }
    Vector giveNodalForces() { return f_ext; }
    int giveStepNumber() const { return step; };
    double giveTime() const { return time; };
    int giveTerminationStatus() const { return ( termination_time - time > 1e-15 ); };
    void setTime(double t);
    void setStep(unsigned t) { step = t; };
    void setTimeStep(double timeStep) { dt = timeStep; };
    void setInitialTimeStep(double timeStep) { initdt = timeStep; };
    double giveTimeStep() const { return dt; };
    void setInitialTimeAndStep(double t, unsigned s) { this->init_time = t; this->init_step = s; }
    virtual void setNextStepTime();
    virtual void runBeforeEachStep();
    virtual void runAfterEachStep();
    virtual void solve() {};
};

//////////////////////////////////////////////////////////
class SteadyStateLinearSolver : public Solver
{
protected:
    double conj_grad_precision;
    double conj_grad_relative_maxit;
    CoordinateIndexedSparseMatrix Keff, K;

    virtual void computeForcesAtIntegrationTime(const bool frozen)  { computeInternalExternalForces(trial_r, load, frozen, -1); }; //do not use dt as this is quasistatic simulation
    virtual void computeForcesAtStepEnd(const bool frozen) { computeInternalExternalForces(trial_r, load, frozen, -1); };
    virtual void computeInternalExternalForces(const Vector &rr, const Vector &ll, const bool frozen, double timeStep);
    virtual void computeKeff();
    virtual void prepareSystemMatricesAndInitialField(string init_r_file, string init_v_file, const bool initial);
private:
public:
    SteadyStateLinearSolver();
    virtual ~SteadyStateLinearSolver();     //destructor
    virtual void init(string init_r_file, string init_v_file, const bool initial = true);
    virtual Solver *readFromFile(const string filename);
    virtual void runBeforeEachStep();
    virtual void runAfterEachStep();
    virtual void solve();
};

//////////////////////////////////////////////////////////
class SteadyStateNonLinearSolver : public SteadyStateLinearSolver
{
protected:
    double dtmax, dtmin;  // for adaptive step
    double W_ext_oldM, W_int_oldM, W_extM, W_intM, W_kinM; //mechanics
    double W_ext_oldT, W_int_oldT, W_extT, W_intT, W_kinT; //transport
    double disErr, resErr, eneErr;
    double limitDisErr, limitResErr, limitEneErr;
    unsigned maxIt; ///> maximum number of iteration
    unsigned enlargeIt, shortenIt; ///> if lower/higher numIt -> enlarge/shorten time-step
    ///> time step is changed according to actual fraction of number of iteration vs maxIt
    double step_increase;  ///> increased in case the number of iteration is in lower 1/3
    double step_decrease;  ///> decreased in case the number of iteration is in upper 1/2
    double critical_step_decrease;  ///> decreased in case of step restart (not satisfying tolerance)

    IndirectDC *idc;        //indirect displacement control
    Vector ddf, full_ddf, f_last_iter;
    double idc_time, idc_dt, idc_time_converged; //time in which load advancements are measured
    double init_idc_time = 0.0;  ///> when starting from previously calculated results
    int stiffnessMatrixUpdate, dampingMatrixUpdate, massMatrixUpdate; //update matrices every X iteration
    void printAllVectors();
    void evaluateErrors(double *displa_error, double *energy_error, double *residu_error);
    virtual bool updateSystemMatrices(string matrixType, unsigned iteration);

public:
    SteadyStateNonLinearSolver();
    virtual ~SteadyStateNonLinearSolver();     //destructor
    virtual void init(string init_r_file, string init_v_file, const bool initial = true);
    virtual Solver *readFromFile(const string filename);
    virtual Vector giveNodalForces() { return f_ext_old; };
    virtual void runBeforeEachStep();
    virtual void runAfterEachStep();
    virtual void solve();
    double giveIDCtime(const bool converged = true);
};

//////////////////////////////////////////////////////////
class TransientLinearTransportSolver : public SteadyStateNonLinearSolver /// solver of first order differential equation in time
{
protected:
    double alpha_f, alpha_m, gamma, beta;
    CoordinateIndexedSparseMatrix C;
    Vector v, v_old;

    virtual void applySpectralRadius(double rhoinfty);
    virtual void computeKeff();
    virtual void prepareSystemMatricesAndInitialField(string init_r_file, string init_v_file, const bool initial);
    virtual bool updateSystemMatrices(string matrixType, unsigned iteration);
    virtual void updateFieldVariables();
    virtual void computeForcesAtIntegrationTime(const bool frozen);
    virtual void computeForcesAtStepEnd(const bool frozen);

public:
    TransientLinearTransportSolver();
    virtual ~TransientLinearTransportSolver();     //destructor
    virtual void init(string init_r_file, string init_v_file, const bool initial = true);
    virtual void solve();
    virtual void runBeforeEachStep();
    virtual void runAfterEachStep();
    virtual Solver *readFromFile(const string filename);
};

//////////////////////////////////////////////////////////
class TransientNonLinearTransportSolver : public TransientLinearTransportSolver /// solver of first order differential equation in time
{
protected:

public:
    TransientNonLinearTransportSolver();
    virtual ~TransientNonLinearTransportSolver();     //destructor
    virtual void init(string init_r_file, string init_v_file, const bool initial = true);
    virtual void solve();
    virtual void runBeforeEachStep();
    virtual void runAfterEachStep();
};

//////////////////////////////////////////////////////////
class TransientLinearMechanicalSolver : public TransientNonLinearTransportSolver /// solver of second order differential equation in time
{
protected:
    CoordinateIndexedSparseMatrix M;
    Vector a, a_old;

    virtual void applySpectralRadius(double rhoinfty);
    virtual void computeKeff();
    virtual void prepareSystemMatricesAndInitialField(string init_r_file, string init_v_file, const bool initial);
    virtual bool updateSystemMatrices(string matrixType, unsigned iteration);
    virtual void updateFieldVariables();
    virtual void computeForcesAtIntegrationTime(const bool frozen);
    virtual void computeForcesAtStepEnd(const bool frozen);
public:
    TransientLinearMechanicalSolver();
    virtual ~TransientLinearMechanicalSolver();     //destructor
    virtual void init(string init_r_file, string init_v_file, const bool initial = true);
    virtual void solve();
    virtual void runBeforeEachStep();
    virtual void runAfterEachStep();
};

//////////////////////////////////////////////////////////
class TransientNonLinearMechanicalSolver : public TransientLinearMechanicalSolver /// solver of second order differential equation in time
{
protected:

public:
    TransientNonLinearMechanicalSolver();
    virtual ~TransientNonLinearMechanicalSolver();     //destructor
    virtual void init(string init_r_file, string init_v_file, const bool initial = true);
    virtual void solve();
    virtual void runBeforeEachStep();
    virtual void runAfterEachStep();
};


#endif
