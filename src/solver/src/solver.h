#ifndef _SOLVER_H
#define _SOLVER_H

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
    std :: string name;
    ElementContainer *elems;
    NodeContainer *nodes;
    FunctionContainer *funcs;
    BCContainer *bcs;
    double time, dt, initdt, termination_time;
    double init_time = 0.0;  ///> when starting from previously calculated results
    Vector f_ext, load, load_old, f_int, pbc, r, f, full_ddr, ddr, residuals;
    Vector f_int_old, f_ext_old, f_dam, f_acc, trial_r;
    unsigned freeDoFnum, totalDoFnum;
    int step;
    unsigned init_step = 0;  ///> when starting from previously calculated results
    bool terminated, fully_converged;
    std :: string symsolver_type = "EigenConj";

    bool showTime;

    virtual void updateFieldVariables();
    virtual void computeInternalExternalForces(const Vector &rr, const Vector &ll, const bool frozen, double timeStep);

public:
    Solver();
    virtual ~Solver() {};
    virtual void init(std :: string init_r_file, std :: string init_v_file, const bool initial = true);
    virtual Solver *readFromFile(const std :: string filename);
    virtual void solveStep() { runBeforeEachStep(); solve();  runAfterEachStep(); };
    void setContainers(ElementContainer *e, NodeContainer *n, FunctionContainer *functions, BCContainer *bc);
    std :: string giveName() const { return name; }
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
    virtual void giveValues(std :: string code, Vector &result) const;
    virtual void rebuild();
    Vector giveResiduals() const { return residuals; };
    bool showStepTime()const { return showTime; };
};


#endif /* _SOLVER_H */
