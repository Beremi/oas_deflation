#ifndef _SOLVER_H
#define _SOLVER_H

#include <vector>
#include <iostream>
#include <fstream>
#include <iomanip>      // std::setw
#include "linalg.h"
#include "linalg_profile.h"
#include "element_container.h"
#include "node_container.h"
#include "indirect_control.h"

class Pertrubation; //forward declaration
class ExporterContainer; //forward declaration

//////////////////////////////////////////////////////////
class Solver
{
protected:
    std :: string name;
    ElementContainer *elems;
    NodeContainer *nodes;
    FunctionContainer *funcs;
    ExporterContainer *exporters;
    BCContainer *bcs;
    double time = 0;
    double dt, initdt, termination_time;
    double init_time = 0.0;  ///> when starting from previously calculated results
    Vector f_ext, load, load_old, f_int, pbc, r, f, full_ddr, ddr, residuals;
    Vector f_int_old, f_ext_old, f_dam, f_acc, trial_r;
    Vector v;
    unsigned freeDoFnum, totalDoFnum;
    int step = 0;
    unsigned init_step = 0;  ///> when starting from previously calculated results
    bool isTimeReal;
    bool terminated, fully_converged;

    Vector W_ext, W_int, W_kin, W_int_old, W_ext_old;
    virtual void computeTotalKineticEnergy() {};
    virtual void computeTotalInternalAndExternalAndKineticEnergy();

    bool showTime;
    bool silent;
    RuntimePhaseProfiler runtimeProfiler;
    bool runtimeProfileEnabled = false;
    std :: string runtimeProfileFileBase = "runtime_profile";
    int runtimeProfileIteration = -1;
    unsigned runtimeProfileCumulIteration = 0;
    std :: string runtimeProfileSystemKind = "solver";

    virtual void updateFieldVariables();
    void initializeRuntimeProfiler();
    void setRuntimeProfileContext(std :: string systemKind, int iteration, unsigned cumulIteration);
    void recordRuntimePhase(std :: string phase, double durationSeconds, std :: string detail = "");
    
    std :: vector< Pertrubation * >pertrubations;

    //Vector lumpMatrix(CoordinateIndexedSparseMatrix &Q) const;

public:
    Solver();
    virtual ~Solver();
    virtual void init(std :: string init_r_file, std :: string init_v_file, const bool initial = true);
    virtual Solver *readFromFile(const std :: string filename);
    virtual void solveStep() { runBeforeEachStep(); solve();  runAfterEachStep(); };
    void setContainers(ElementContainer *e, NodeContainer *n, FunctionContainer *functions, BCContainer *bc, ExporterContainer *exp);
    std :: string giveName() const { return name; }
    bool isTerminated() { return terminated; }
    bool convergedToTolerance() const { return fully_converged; };
    Vector giveDoFValues() const { return r; }
    double giveDoFValue(unsigned i) const { return r [ i ]; }
    Vector giveTrialDoFValues() const { return trial_r; }
    Vector giveExternalForces() const {  return f_ext; }
    double giveTrialDoFValue(unsigned k) const { return trial_r [ k ]; }
    double giveDoFVelocity(unsigned k) const { return v [ k ]; }
    double giveExternalForce(unsigned k) const;
    Vector giveNodalForces() { return f_ext; }
    double giveDoFInertiaForce(unsigned i) const { return f_acc [ i ]; }
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
    bool showStepTime() const { return showTime; };
    bool isSilent() const { return silent; };
    unsigned giveTotalDoFNumber() const { return totalDoFnum; };
    void computeInternalExternalForces(const Vector &rr, Vector &ll, const bool frozen, double t, double timeStep);
    Vector giveInternalForces() const {return f_int;};
    virtual void updateSystemMatrices() {};    
    virtual CoordinateIndexedSparseMatrix giveSystemMatrix() const {return CoordinateIndexedSparseMatrix(freeDoFnum, freeDoFnum);};  
    void setAllDoFsExternally(Vector DoFs);    
};

//////////////////////////////////////////////////////////
class Pertrubation
{
protected:
    std :: string name;
    double finalized;
    double time;
    double magnitude;
    double seed;
public:
    Pertrubation() {};
    virtual ~Pertrubation() { finalized = false; name = "Pertrubation"; };
    bool shouldBeApplied(double solverTime) const;
    Vector pertrube(unsigned size);
    virtual void readFromLine(std :: istringstream &iss);
};

#endif /* _SOLVER_H */
