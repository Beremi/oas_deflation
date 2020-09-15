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
    double time, dt, termination_time;
    CoordinateIndexedSparseMatrix K, Kini;
    Vector f_ext, load, f_int, pbc, r, f, full_ddr, ddr;
    unsigned freeDoFnum, fixedDoFnum, totalDoFnum;
    int step;
    bool terminated;
    virtual void setNextStepTime();
    virtual void runBeforeEachStep();
    virtual void runAfterEachStep();
    virtual void solve() {};
    virtual void computeInternalExternalForces(Vector &rr);
    virtual void computeInternalExternalForcesWithFrozenIntVariables(Vector &rr);

public:
    Solver() { name = "basic solver"; };
    virtual ~Solver() {};
    virtual void init();
    virtual Solver *readFromFile(const string filename);
    virtual void solveStep() { runBeforeEachStep(); solve(); runAfterEachStep(); };
    void setContainers(ElementContainer *e, NodeContainer *n, FunctionContainer *functions) { elems = e; nodes = n; funcs = functions;}
    string giveName() const { return name; }
    bool isTerminated() { return terminated; }
    Vector giveDoFValues() { return r; }
    Vector giveNodalForces() { return f_ext; }
    int giveStepNumber() const { return step; };
    double giveTime() const { return time; };
    int giveTerminationStatus() const { return ( termination_time - time > 1e-15 ); };
};

//////////////////////////////////////////////////////////
class SteadyStateLinearSolver : public Solver
{
protected:
    virtual void runBeforeEachStep();
    virtual void runAfterEachStep();
    virtual void solve();
    double conj_grad_precission;
    double conj_grad_relative_maxit;
private:

public:
    SteadyStateLinearSolver();
    virtual ~SteadyStateLinearSolver();     //destructor
    virtual void init();
    virtual Solver *readFromFile(const string filename);
};

//////////////////////////////////////////////////////////
class SteadyStateNonLinearSolver : public SteadyStateLinearSolver
{
protected:
    double dtmax, dtmin;  // for adaptive step
    Vector f_int_old, f_ext_old, residual;
    Vector trial_r;
    double W_ext_oldM, W_int_oldM, W_extM, W_intM; //mechanics
    double W_ext_oldT, W_int_oldT, W_extT, W_intT; //transport
    double disErr, resErr, eneErr;
    double limitDisErr, limitResErr, limitEneErr;
    unsigned maxIt;
    ///> time step is changed according to actual fraction of number of iteration vs maxIt
    double step_increase;  ///> increased in case the number of iteration is in lower 1/3
    double step_decrease;  ///> decreased in case the number of iteration is in upper 1/2
    double critical_step_decrease;  ///> decreased in case of step restart (not satisfying tolerance)

    IndirectDC *idc;        //indirect displacement control
    Vector ddf, full_ddf, f_last_iter;
    double idc_time, idc_dt, idc_time_converged; //time in which load advancements are masured

    virtual void runBeforeEachStep();
    virtual void runAfterEachStep();
    virtual void solve();
    void printAllVectors();
    void evaluateErrors(double *displa_error, double *energy_error, double *residu_error);
private:

public:
    SteadyStateNonLinearSolver();
    virtual ~SteadyStateNonLinearSolver();     //destructor
    virtual void init();
    virtual Solver *readFromFile(const string filename);
    virtual Vector giveNodalForces() { return f_ext_old; };
};

//////////////////////////////////////////////////////////
class TransientLinearMechanicalSolver : public SteadyStateLinearSolver /// solver of second order differential equation in time
{
protected:
    double alpha_f, alpha_m, gamma, beta, rhoinfty;
    CoordinateIndexedSparseMatrix M, C, Keff;
    Vector a, v, a_red, v_red, r_red, feff, r_old, r_f;
    virtual void solve();
    virtual void updateKeff();
    virtual void updateFeff();
    virtual void updateFieldVariables();
    virtual void applySpectralRadius(double rhoinfty);
private:

public:
    TransientLinearMechanicalSolver();
    virtual ~TransientLinearMechanicalSolver();     //destructor
    virtual void init();
    virtual Solver *readFromFile(const string filename);
};

//////////////////////////////////////////////////////////
class TransientLinearTransportSolver : public TransientLinearMechanicalSolver /// solver of first order differential equation in time
{
protected:
    virtual void updateKeff();
    virtual void updateFeff();
    virtual void updateFieldVariables();
    virtual void applySpectralRadius(double rhoinfty);
private:

public:
    TransientLinearTransportSolver();
    virtual ~TransientLinearTransportSolver();     //destructor
    virtual void init();
};

#endif
