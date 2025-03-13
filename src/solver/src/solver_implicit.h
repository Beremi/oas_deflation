#ifndef _SOLVER_IMPLICIT_H
#define _SOLVER_IMPLICIT_H

#include "solver.h"

//////////////////////////////////////////////////////////
class SteadyStateLinearSolver : public Solver
{
protected:
    double conj_grad_precision;
    double conj_grad_relative_maxit;
    CoordinateIndexedSparseMatrix Keff, K;
    std :: unique_ptr< LinAlgSolver >linalgsolver;
    std :: string symsolver_type = "EigenConj";

    int stiffnessMatrixUpdate; //update matrices every X iteration
    std :: string stiffMatType, stiffMatTypeFirstIT;

    virtual bool updateSystemMatrices(unsigned iteration, bool enforce);

    virtual void computeForcesAtIntegrationTime(const bool frozen)  { computeInternalExternalForces(trial_r, load, frozen, -1); }; //do not use dt as this is quasistatic simulation
    virtual void computeForcesAtStepEnd(const bool frozen) { computeInternalExternalForces(trial_r, load, frozen, -1); };
    virtual void computeKeff();
    virtual void prepareSystemMatricesAndInitialField(std :: string init_r_file, std :: string init_v_file, const bool initial);
private:
public:
    SteadyStateLinearSolver();
    virtual ~SteadyStateLinearSolver();     //destructor
    virtual void init(std :: string init_r_file, std :: string init_v_file, const bool initial = true);
    virtual Solver *readFromFile(const std :: string filename);
    virtual void runBeforeEachStep();
    virtual void runAfterEachStep();
    virtual void solve();
    virtual void reset();
    virtual void rebuild();
    virtual void factorizeLinearSystem();
};

//////////////////////////////////////////////////////////
class SteadyStateNonLinearSolver : public SteadyStateLinearSolver
{
protected:
    unsigned it, restarts; //number of iterations, number of restarts
    double dtmax, dtmin;  // for adaptive step
    Vector EPS2;
    double disErr, resErr, eneErr;
    double maxDisErr, maxResErr, maxEneErr;
    double limitDisErr, limitResErr, limitEneErr;
    unsigned maxIt; ///> maximum number of iteration
    unsigned minIt; ///> minimum number of iteration
    unsigned enlargeIt, shortenIt; ///> if lower/higher numIt -> enlarge/shorten time-step
    ///> time step is changed according to actual fraction of number of iteration vs maxIt
    double step_increase;  ///> increased in case the number of iteration is in lower 1/3
    double step_decrease;  ///> decreased in case the number of iteration is in upper 1/2
    double critical_step_decrease;  ///> decreased in case of step restart (not satisfying tolerance)

    IndirectDC *idc;        //indirect displacement control
    Vector ddf, full_ddf, f_last_iter;
    double idc_time, idc_dt, idc_time_converged; //time in which load advancements are measured
    double init_idc_time = 0.0;  ///> when starting from previously calculated results
    void printAllVectors();
    void checkAllVectorsForNaNs();
    void evaluateErrors();
    virtual void reset();
    virtual void giveValues(std :: string code, Vector &result) const;

public:
    SteadyStateNonLinearSolver();
    virtual ~SteadyStateNonLinearSolver();     //destructor
    virtual void init(std :: string init_r_file, std :: string init_v_file, const bool initial = true);
    virtual Solver *readFromFile(const std :: string filename);
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
    Vector v_old;
    unsigned timeIntM; //0 - generalized alphal; 1 - HHT; 2 - Newmark
    bool check_time_integr_params;
    int dampingMatrixUpdate; //update matrices every X iteration

    virtual void applySpectralRadius(double rhoinfty);
    virtual void checkIntegrationParams();
    virtual void setDefaultIntegrationParams();
    virtual void computeKeff();
    virtual void prepareSystemMatricesAndInitialField(std :: string init_r_file, std :: string init_v_file, const bool initial);
    virtual bool updateSystemMatrices(unsigned iteration, bool enforce);
    virtual void updateFieldVariables();
    virtual void computeForcesAtIntegrationTime(const bool frozen);
    virtual void computeForcesAtStepEnd(const bool frozen);

public:
    TransientLinearTransportSolver();
    virtual ~TransientLinearTransportSolver();     //destructor
    virtual void init(std :: string init_r_file, std :: string init_v_file, const bool initial = true);
    virtual void solve();
    virtual void runBeforeEachStep();
    virtual void runAfterEachStep();
    virtual Solver *readFromFile(const std :: string filename);
    virtual void rebuild();
};

//////////////////////////////////////////////////////////
class TransientNonLinearTransportSolver : public TransientLinearTransportSolver /// solver of first order differential equation in time
{
protected:

public:
    TransientNonLinearTransportSolver();
    virtual ~TransientNonLinearTransportSolver();     //destructor
    virtual void init(std :: string init_r_file, std :: string init_v_file, const bool initial = true);
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
    virtual void computeTotalKineticEnergy();
    virtual void applySpectralRadius(double rhoinfty);
    virtual void checkIntegrationParams();
    virtual void setDefaultIntegrationParams();
    virtual void computeKeff();
    virtual void prepareSystemMatricesAndInitialField(std :: string init_r_file, std :: string init_v_file, const bool initial);
    virtual bool updateSystemMatrices(unsigned iteration, bool enforce);
    virtual void updateFieldVariables();
    virtual void computeForcesAtIntegrationTime(const bool frozen);
    virtual void computeForcesAtStepEnd(const bool frozen);
    int massMatrixUpdate; //update matrices every X iteration
    bool lumpMassM;
public:
    TransientLinearMechanicalSolver();
    virtual ~TransientLinearMechanicalSolver();     //destructor
    virtual void init(std :: string init_r_file, std :: string init_v_file, const bool initial = true);
    virtual void solve();
    virtual void runBeforeEachStep();
    virtual void runAfterEachStep();
    virtual Solver *readFromFile(const std :: string filename);
};

//////////////////////////////////////////////////////////
class TransientNonLinearMechanicalSolver : public TransientLinearMechanicalSolver /// solver of second order differential equation in time
{
protected:

public:
    TransientNonLinearMechanicalSolver();
    virtual ~TransientNonLinearMechanicalSolver();     //destructor
    virtual void init(std :: string init_r_file, std :: string init_v_file, const bool initial = true);
    virtual void solve();
    virtual void runBeforeEachStep();
    virtual void runAfterEachStep();
};


#endif /* _SOLVER_IMPLICIT_H */
