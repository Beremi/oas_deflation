#ifndef SOLVER_H
#define SOLVER_H

#include <vector>
#include <iostream>
#include <fstream>
#include <iomanip>      // std::setw
#include "linear_algebra.h"
#include "element_container.h"
#include "node_container.h"

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

public:
    Solver() { name = "basic solver"; };
    virtual ~Solver() {};
    virtual void init();
    virtual Solver *readFromLine(istringstream &iss);
    virtual void solveStep() { runBeforeEachStep(); solve(); runAfterEachStep(); };
    void setElementContainer(ElementContainer *e) { elems = e; }
    void setNodeContainer(NodeContainer *n) { nodes = n; };
    void setFunctionContainer(FunctionContainer *functions) { this->funcs = functions; };
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
    virtual Solver *readFromLine(istringstream &iss);
    virtual void computeInternalExternalForces(Vector &rr);
};

//////////////////////////////////////////////////////////
class SteadyStateNonLinearSolver : public SteadyStateLinearSolver
{
protected:
    double dtmax, dtmin;  // for adaptive step
    Vector f_int_old, f_ext_old, residual;
    Vector trial_r;
    double W_ext_old, W_int_old, W_ext, W_int;
    double disErr, resErr, eneErr;
    double limitDisErr, limitResErr, limitEneErr;
    unsigned maxIt;

    virtual void runBeforeEachStep();
    virtual void runAfterEachStep();
    virtual void solve();
    void printAllVectors();
private:

public:
    SteadyStateNonLinearSolver();
    virtual ~SteadyStateNonLinearSolver();     //destructor
    virtual void init();
    virtual Solver *readFromLine(istringstream &iss);
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
    virtual Solver *readFromLine(istringstream &iss);
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
