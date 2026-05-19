#ifndef _SOLVER_IMPLICIT_H
#define _SOLVER_IMPLICIT_H

#include "solver.h"
#include <memory>

//////////////////////////////////////////////////////////
class SteadyStateLinearSolver : public Solver
{
protected:
    double conj_grad_precision;
    double conj_grad_relative_maxit;
    DeflatedFGMRESOptions dfgmresOptions;
    HypreBoomerAMGOptions hypreOptions;
    double elasticNearNullspaceCoordinateScale = -1.;
    CoordinateIndexedSparseMatrix Keff, K;
    std :: unique_ptr< LinAlgSolver >linalgsolver;
    std :: string symsolver_type = "EigenConj";

    int stiffnessMatrixIterUpdate; //update matrices every X iteration
    int stiffnessMatrixCumulIterUpdate; //update matrices every X iteration cumulatively
    int stiffnessMatrixStepUpdate; //update matrices every X step
    std :: string stiffMatType, stiffMatTypeFirstIT;
    double stiffMatElasticBlendBeta = 0.;

    virtual bool updateSystemMatrices(unsigned iteration, unsigned cumul_iteration, bool enforce);

    virtual void computeForcesAtIntegrationTime(const bool frozen)  { computeInternalExternalForces(trial_r, load, frozen, time, -1); }; //do not use dt as this is quasistatic simulation
    virtual void computeForcesAtStepEnd(const bool frozen) { computeInternalExternalForces(trial_r, load, frozen, time, -1); };
    virtual void computeKeff();
    virtual void prepareSystemMatricesAndInitialField(std :: string init_r_file, std :: string init_v_file, const bool initial);
    void collectLinearDeflationVector(const Vector &x, bool success);
    ElasticDofMap buildElasticDofMap() const;
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
    virtual void updateSystemMatrices() {updateSystemMatrices(0, 0, 1); computeKeff();};    
    virtual CoordinateIndexedSparseMatrix giveSystemMatrix() const {return Keff;};    
};

//////////////////////////////////////////////////////////
class SteadyStateNonLinearSolver : public SteadyStateLinearSolver
{
protected:
    enum class NonlinearDampingType { Off, Fixed, Adaptive, RollbackAdaptive };
    enum class NonlinearLineSearchType { Off, Backtracking, Bisection };
    enum class NonlinearMeritType { Residual, Energy, Mixed };
    enum class NonlinearLineSearchEvaluationMode { Frozen, Actual, FrozenThenActual };
    enum class NonlinearTrustRegionType { Off, StepNorm };

    struct NonlinearStateSnapshot {
        Vector trial_r;
        Vector ddr;
        Vector full_ddr;
        Vector f;
        Vector f_int;
        Vector f_ext;
        Vector f_dam;
        Vector f_acc;
        Vector residuals;
        Vector load;
        Vector W_int;
        Vector W_ext;
        Vector W_kin;
        double disErr = 0.;
        double resErr = 0.;
        double eneErr = 0.;
        std :: shared_ptr< ElementContainer :: MaterialStatusSnapshot > materialStatuses;
        std :: uint64_t materialStatusHash = 0;
        bool materialStatusHashValid = false;
    };

    struct NonlinearTrialResult {
        bool accepted = false;
        double alpha = 1.;
        unsigned trials = 0;
        double meritBefore = 0.;
        double meritAfter = 0.;
    };

    unsigned it, restarts, cumul_it; //number of iterations, number of restarts
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

    IndirectControl *idc;        //indirect displacement control
    Vector ddf, full_ddf, f_last_iter;
    Vector eigen_trial_rPF, eigen_f_extPF, eigen_WextPF; //RVE fields to transfer macroscopic data to error evaluation
    double idc_time, idc_dt, idc_time_converged; //time in which load advancements are measured
    double init_idc_time = 0.0;  ///> when starting from previously calculated results
    NonlinearDampingType nonlinearDampingType = NonlinearDampingType :: Off;
    NonlinearLineSearchType nonlinearLineSearchType = NonlinearLineSearchType :: Off;
    NonlinearMeritType nonlinearMeritType = NonlinearMeritType :: Mixed;
    double nonlinearDampingFactor = 1.;
    double nonlinearDampingMin = 0.03125;
    double nonlinearDampingMax = 1.;
    double nonlinearDampingIncrease = 1.25;
    double nonlinearDampingDecrease = 0.5;
    double currentNonlinearDampingFactor = 1.;
    double nonlinearLineSearchReduction = 0.5;
    double nonlinearLineSearchMinAlpha = 0.03125;
    double nonlinearLineSearchMaxAlpha = 1.0;
    double nonlinearLineSearchBisectionTolerance = 1e-3;
    unsigned nonlinearLineSearchMaxTrials = 6;
    double nonlinearLineSearchArmijo = 0.;
    bool nonlinearLineSearchAcceptAnyDecrease = true;
    bool nonlinearLineSearchFreezeMaterial = true;
    bool nonlinearLineSearchCutbackOnFail = true;
    bool nonlinearLineSearchReportTrials = false;
    NonlinearLineSearchEvaluationMode nonlinearLineSearchEvaluationMode = NonlinearLineSearchEvaluationMode :: Frozen;
    bool nonlinearStagnationCutback = false;
    unsigned nonlinearStagnationIterations = 8;
    double nonlinearStagnationRatio = 0.95;
    double nonlinearGrowthCutback = 1.25;
    bool nonlinearAdaptiveMatrixUpdate = false;
    double nonlinearRebuildOnSmallAlpha = 0.;
    bool nonlinearRebuildOnMeritGrowth = false;
    bool nonlinearRebuildOnStagnation = false;
    bool nonlinearForceMatrixRebuild = false;
    bool lastNonlinearMatrixRebuild = false;
    NonlinearTrustRegionType nonlinearTrustRegionType = NonlinearTrustRegionType :: Off;
    double nonlinearTrustRadiusInitial = 0.;
    double nonlinearTrustRadiusMin = 1e-12;
    double nonlinearTrustRadiusMax = 1e100;
    double nonlinearTrustShrink = 0.5;
    double nonlinearTrustExpand = 1.5;
    unsigned nonlinearTrustMaxTrials = 8;
    double currentNonlinearTrustRadius = 0.;
    double lastNonlinearTrustRadius = 0.;
    double nonlinearRollbackGrowthRatio = 1.10;
    double nonlinearRollbackIncreaseRatio = 1.00;
    double nonlinearRollbackGoodReductionRatio = 0.50;
    unsigned nonlinearRollbackIgnoreInitialIterations = 0;
    bool nonlinearRollbackEnable = true;
    bool nonlinearRollbackDecreaseOnAnyIncrease = true;
    unsigned nonlinearRollbackMaxGrowthDecreases = 1000000000;
    unsigned nonlinearRollbackGrowthDecreaseCounter = 0;
    bool nonlinearMaterialSnapshotRollback = false;
    bool nonlinearMaterialSnapshotVerify = false;
    bool nonlinearTangentCheck = false;
    unsigned nonlinearTangentCheckStep = 0;
    unsigned nonlinearTangentCheckIteration = 0;
    double nonlinearTangentCheckEps = 1e-6;
    unsigned nonlinearTangentCheckRandomVectors = 0;
    bool nonlinearTangentCheckIncludeNewton = true;
    bool nonlinearTangentCheckStopAfter = false;
    std :: string nonlinearTangentCheckOutput = "tangent_check.tsv";
    bool nonlinearTangentCheckDone = false;
    double nonlinearBestMerit = 0.;
    unsigned nonlinearStagnationCounter = 0;
    double lastNonlinearAlpha = 1.;
    unsigned lastNonlinearLineSearchTrials = 0;
    double lastNonlinearMeritBefore = 0.;
    double lastNonlinearMeritAfter = 0.;
    double lastNonlinearConvergenceMerit = 0.;
    std :: string nonlinearCutbackReason;
    bool warnedGlobalizationWithIDC = false;
    void printAllVectors();
    void checkAllVectorsForNaNs();
    void evaluateErrors();
    bool nonlinearGlobalizationActive() const;
    bool nonlinearConvergenceCriteriaSatisfied() const;
    double currentNonlinearMerit() const;
    double currentNonlinearGlobalizationMerit() const;
    double currentNonlinearConvergenceMerit() const;
    bool nonlinearMeritAccepted(double trialMerit, double baseMerit, double alpha) const;
    NonlinearStateSnapshot saveNonlinearState() const;
    void restoreNonlinearState(const NonlinearStateSnapshot &snapshot, bool resetMaterialStatuses);
    void resetNonlinearGlobalizationAttempt();
    bool applyScaledIncrementAndEvaluate(const NonlinearStateSnapshot &baseState, const Vector &increment, double alpha, bool frozen, bool resetMaterialStatuses = true);
    NonlinearTrialResult performBacktrackingLineSearch(const NonlinearStateSnapshot &baseState, const Vector &increment, double meritBefore);
    NonlinearTrialResult performBisectionLineSearch(const NonlinearStateSnapshot &baseState, const Vector &increment, double meritBefore);
    NonlinearTrialResult performStepNormTrustRegion(const NonlinearStateSnapshot &baseState, const Vector &increment, double meritBefore);
    double currentNonlinearDampingAlpha() const;
    void updateAdaptiveNonlinearDamping(double meritBefore, double meritAfter);
    bool shouldCutbackForNonlinearProgress(double meritBefore, double meritAfter);
    void requestAdaptiveMatrixRebuildIfNeeded(double meritBefore, double meritAfter);
    bool maybeRunNonlinearTangentCheck(const Vector &newtonIncrement);
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
    void setEigenErrorValue_rPF(unsigned pf, double value);
    void setEigenErrorValue_fext(unsigned pf, double value);
    void setEigenErrorValue_Wext(unsigned pf, double value);
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
    int dampingMatrixIterUpdate; //update matrices every X iteration
    int dampingMatrixCumulIterUpdate; //update matrices every X iteration cumulatively
    int dampingMatrixStepUpdate; //update matrices every X step

    virtual void applySpectralRadius(double rhoinfty);
    virtual void checkIntegrationParams();
    virtual void setDefaultIntegrationParams();
    virtual void computeKeff();
    virtual void prepareSystemMatricesAndInitialField(std :: string init_r_file, std :: string init_v_file, const bool initial);
    virtual bool updateSystemMatrices(unsigned iteration, unsigned cumul_iteration, bool enforce);
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
    virtual bool updateSystemMatrices(unsigned iteration, unsigned cumul_iteration, bool enforce);
    virtual void updateFieldVariables();
    virtual void computeForcesAtIntegrationTime(const bool frozen);
    virtual void computeForcesAtStepEnd(const bool frozen);
    int massMatrixIterUpdate; //update matrices every X iteration
    int massMatrixCumulIterUpdate; //update matrices every X iteration
    int massMatrixStepUpdate; //update matrices every X step
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
