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
    CoordinateIndexedSparseMatrix K, Kini;
    Vector f_ext, load, f_int, pbc, r, f, full_ddr, ddr, residuals;
    Vector f_int_old, f_ext_old, f_dam, f_acc, trial_r;
    unsigned freeDoFnum, fixedDoFnum, totalDoFnum;
    int step;
    bool terminated;

    virtual void computeInternalExternalForces(const Vector &rr, bool frozen);

public:
    Solver();
    virtual ~Solver() {};
    virtual void init(const bool &initial = true);
    virtual Solver *readFromFile(const string filename);
    virtual void solveStep() { runBeforeEachStep(); solve();  runAfterEachStep(); };
    void setContainers(ElementContainer *e, NodeContainer *n, FunctionContainer *functions) { elems = e; nodes = n; funcs = functions; }
    string giveName() const { return name; }
    bool isTerminated() { return terminated; }
    Vector giveDoFValues() const { return r; }
    Vector giveNodalForces() { return f_ext; }
    int giveStepNumber() const { return step; };
    double giveTime() const { return time; };
    int giveTerminationStatus() const { return ( termination_time - time > 1e-15 ); };
    void setTime(double t);
    void setStep(unsigned t) { step = t; };
    virtual void setNextStepTime();
    virtual void runBeforeEachStep();
    virtual void runAfterEachStep();
    virtual void solve() { cout << "SOLVING CRAP " << endl; cout.flush();};
};

//////////////////////////////////////////////////////////
class SteadyStateLinearSolver : public Solver
{
protected:
    double conj_grad_precission;
    double conj_grad_relative_maxit;
private:

public:
    SteadyStateLinearSolver();
    virtual ~SteadyStateLinearSolver();     //destructor
    virtual void init(const bool &initial = true);
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

    void printAllVectors();
    void evaluateErrors(double *displa_error, double *energy_error, double *residu_error);

public:
    SteadyStateNonLinearSolver();
    virtual ~SteadyStateNonLinearSolver();     //destructor
    virtual void init(const bool &initial = true);
    virtual Solver *readFromFile(const string filename);
    virtual Vector giveNodalForces() { return f_ext_old; };
    virtual void runBeforeEachStep();
    virtual void runAfterEachStep();
    virtual void solve();
};


//////////////////////////////////////////////////////////
class TransientLinearTransportSolver : public SteadyStateNonLinearSolver /// solver of first order differential equation in time
{
protected:
    double alpha_f, alpha_m, gamma, beta, rhoinfty;
    CoordinateIndexedSparseMatrix C, Keff;
    Vector v, v_old;

    virtual void applySpectralRadius(double rhoinfty);
    virtual void updateKeff();
    //virtual void updateFeff();
    virtual void updateFieldVariables();
    virtual void computeInternalExternalDampingForces(const Vector &rr, const Vector &vv, bool frozen);

public:
    TransientLinearTransportSolver();
    virtual ~TransientLinearTransportSolver();     //destructor
    virtual void init(const bool &initial = true);
    virtual void solve();
    virtual void runBeforeEachStep();
    virtual void runAfterEachStep();
};



//////////////////////////////////////////////////////////
class TransientNonLinearTransportSolver : public TransientLinearTransportSolver /// solver of first order differential equation in time
{
protected:

public:
    TransientNonLinearTransportSolver();
    virtual ~TransientNonLinearTransportSolver();     //destructor
    virtual void init(const bool &initial = true);
    virtual void solve();
    virtual void runBeforeEachStep();
    virtual void runAfterEachStep();
};


/*
//////////////////////////////////////////////////////////
class TransientNonLinearMechanicalSolver : public SteadyStateNonLinearSolver /// solver of second order differential equation in time
{
protected:
    double alpha_f, alpha_m, gamma, beta, rhoinfty;
    CoordinateIndexedSparseMatrix M, C, Keff;
    Vector a, v, a_old, v_old, feff, r_f, dr;

    virtual void computeInternalExternalForces(Vector &rr);
    virtual void computeInternalExternalForcesWithFrozenIntVariables(Vector &rr);

private:

public:
    TransientNonLinearMechanicalSolver();
    virtual ~TransientNonLinearMechanicalSolver();     //destructor
    virtual void init(const bool &initial = true);
    virtual Solver *readFromFile(const string filename);
    virtual void solve();
    virtual void runAfterEachStep(){SteadyStateNonLinearSolver :: runAfterEachStep(); computeInternalExternalForces(r);};
};

//////////////////////////////////////////////////////////
class TransientLinearMechanicalSolver : public TransientNonLinearMechanicalSolver /// solver of second order differential equation in time
{
protected:

public:
    TransientLinearMechanicalSolver();
    virtual ~TransientLinearMechanicalSolver();     //destructor
    virtual void solve();
    virtual void runBeforeEachStep(){SteadyStateLinearSolver :: runBeforeEachStep();};
    virtual void runAfterEachStep(){SteadyStateLinearSolver :: runAfterEachStep();};
};
*/


#endif
