#ifndef SOLVER_H
#define SOLVER_H

#include <vector>
#include <iostream>
#include <fstream>
#include <iomanip>      // std::setw
#include "linear_algebra.h"
#include "element_container.h"
#include "node_container.h"


class Solver
{
protected:
    string name;
    ElementContainer *elems;
    NodeContainer *nodes;
    double time, dt, termination_time;
    CoordinateIndexedSparseMatrix K;
    Vector f_ext, load, f_int, pbc, r, f, full_ddr, ddr;
    unsigned freeDoFnum, fixedDoFnum, totalDoFnum;
    int step;
    bool terminated;
    virtual void runBeforeEachStep();
    virtual void runAfterEachStep();

public:
    Solver() { name = "basic solver"; };
    virtual ~Solver() {};
    virtual void init();
    virtual Solver *readFromLine(istringstream &iss);
    virtual void solveStep() {};
    void setElementContainer(ElementContainer *e) { elems = e; }
    void setNodeContainer(NodeContainer *n) { nodes = n; };
    string giveName() const { return name; }
    bool isTerminated() { return terminated; }
    Vector giveDoFValues() { return r; }
    Vector giveNodalForces() { return f_ext; }
    int giveStepNumber() const { return step; };
};

class SteadyStateLinearSolver : public Solver
{
protected:
    virtual void runBeforeEachStep();
    virtual void runAfterEachStep();
    virtual void solve();
private:

public:
    SteadyStateLinearSolver();
    virtual ~SteadyStateLinearSolver();     //destructor
    virtual void init();
    virtual void solveStep() { runBeforeEachStep(); solve(); runAfterEachStep(); };
    virtual Solver *readFromLine(istringstream &iss);
    virtual void computeInternalExternalForces(Vector &rr);
};

class SteadyStateNonLinearSolver : public SteadyStateLinearSolver
{
protected:
    Vector f_int_old, f_ext_old, residual;
    Vector trial_r;
    double W_ext_old, W_int_old, W_ext, W_int;
    double disErr, resErr, eneErr;

    virtual void runBeforeEachStep();
    virtual void runAfterEachStep();
    virtual void solve();
private:

public:
    SteadyStateNonLinearSolver();
    virtual ~SteadyStateNonLinearSolver();     //destructor
    virtual void init();
    virtual void solveStep() { runBeforeEachStep(); solve(); runAfterEachStep(); };
    virtual Solver *readFromLine(istringstream &iss);
};

#endif
