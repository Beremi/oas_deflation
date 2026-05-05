#ifndef _LINALG_H
#define _LINALG_H

#include <valarray>
#include <map>
#include <chrono>
#include <ctime>
#include <complex>
#include <algorithm>
#include <string>
#include <stdio.h>
#include <assert.h>
#include <numeric>
#include <iostream>
#include <iomanip>
#include <vector>
#include <memory>
#include <Eigen/Core>
#include <Eigen/Dense>
#include <Eigen/IterativeLinearSolvers>
#include <Eigen/SparseLU>
#include <Eigen/SparseCholesky>
#include <Eigen/Eigenvalues>
#include <Eigen/SparseCore>
#ifdef SUPERLU_FOUND
 #include <Eigen/SuperLUSupport>
#endif
#ifdef PARDISO_FOUND
 #include <Eigen/PardisoSupport>
#endif
#ifdef CHOLMOD_FOUND
 #include <Eigen/CholmodSupport>
#endif
#include <Spectra/SymEigsSolver.h>
#include <Spectra/SymGEigsSolver.h>
#include <Spectra/MatOp/SparseSymMatProd.h>
#include <Spectra/MatOp/SparseGenMatProd.h>
#include <Spectra/MatOp/SparseCholesky.h>

//  \[ ([\w\s\+\*]*) \] \[ ([\w\s\+\*]*) \]
#include "globals.h"

using Ttripletd = Eigen :: Triplet< double >;

using Vector = Eigen :: VectorXd;
using Point = Eigen :: Vector3d;
using Matrix = Eigen :: MatrixXd;
typedef typename Eigen :: SparseMatrix< double, Eigen :: ColMajor >CoordinateIndexedSparseMatrix;   // row-major-sparse * dense vector/matrix products - multi-threading

const static Eigen :: IOFormat VectorSemicolonFmt(Eigen :: FullPrecision, Eigen :: DontAlignCols, "; ", "; ", "", "", "", "");

//////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////
// SOLVER FOR LINEAR ALGEBRA, BASE CLASS
class LinAlgSolver
{
protected:
    std :: string name;
public:
    LinAlgSolver() {};
    virtual ~LinAlgSolver() {};
    virtual bool analyzePattern(const CoordinateIndexedSparseMatrix &A) { ( void ) A; return false; };
    virtual bool factorize(const CoordinateIndexedSparseMatrix &A) { ( void ) A; name = "null solver, base class"; return false; };
    virtual bool solve(Vector &x, const Vector &b) { ( void ) x; ( void ) b; return false; };
    std :: string giveName()const { return name; };
    virtual long long giveLastIterations() const { return 0; };
    virtual double giveLastError() const { return 0.; };
    virtual double giveLastTrueRelativeResidual() const { return giveLastError(); };
    virtual void setRuntimeTolerance(double tolerance, double trueTolerance) { ( void ) tolerance; ( void ) trueTolerance; };
    virtual double giveRuntimeTolerance() const { return 0.; };
    virtual double giveRuntimeTrueTolerance() const { return giveRuntimeTolerance(); };
    virtual void collectDeflationVector(const Vector &v) { ( void ) v; };
    virtual unsigned giveDeflationBasisSize() const { return 0; };
    virtual unsigned giveDeflationDiscardedCount() const { return 0; };
    virtual unsigned giveDeflationRawCandidateCount() const { return 0; };
    virtual unsigned giveDeflationCapacityEvictionCount() const { return 0; };
    virtual unsigned giveDeflationLowANormDiscardCount() const { return 0; };
    virtual double giveDeflationLastInitialANorm() const { return 0.; };
    virtual double giveDeflationLastFinalANorm() const { return 0.; };
    virtual double giveDeflationOrthogonalityMaxOffdiag() const { return 0.; };
    virtual double giveDeflationOrthogonalityMaxDiagError() const { return 0.; };
    virtual std :: string giveDeflationLastDiscardReason() const { return ""; };
    virtual double giveLastPreconditionerApplySeconds() const { return 0.; };
    virtual double giveLastOrthogonalizationSeconds() const { return 0.; };
    virtual double giveLastLeastSquaresSeconds() const { return 0.; };
    virtual double giveLastMatvecSeconds() const { return 0.; };
    virtual double giveLastDeflationSeconds() const { return 0.; };
};

//////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////
// Reduced-to-elastic-block metadata for vector AMG setup
struct ElasticDofMap
{
    unsigned fullRows = 0;
    unsigned blockSize = 0;
    unsigned dimension = 0;
    std :: vector< unsigned >reducedToFull;
    std :: vector< double >coordinates;
    std :: vector< double >nearNullspace;
    int nearNullspaceColumns = 0;

    bool isValid() const {
        return fullRows > 0 && blockSize > 0 && reducedToFull.size() > 0;
    }
};

//////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////
// AMGCL solver options
struct AmgclSolverOptions
{
    double tolerance = 1e-6;
    double trueTolerance = -1.;
    unsigned maxIterations = 500;
    double epsStrong = 0.0;
    double relax = 1.0;
    unsigned blockSize = 1;
    unsigned coarseEnough = 0;
    unsigned npre = 1;
    unsigned npost = 1;
    unsigned ncycle = 1;
    std :: string coarsening = "smoothed_aggregation";
    std :: string krylov = "cg";
    unsigned gmresRestart = 80;
    bool estimateSpectralRadius = false;
    unsigned powerIterations = 0;
    bool nearNullspace = true;
    bool elasticFullLift = true;
    bool diagonalScale = false;
    bool useBlockBackend = true;
    std :: string backend = "auto";
    std :: string blockRelaxation = "ilu0";
    bool reuseInitialGuess = false;
    bool checkMatrix = false;
    bool warnOnFailure = true;
    bool verbose = false;
};

//////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////
// Deflated FGMRES solver options
struct DeflatedFGMRESOptions
{
    double tolerance = 1e-6;
    double trueTolerance = 1e-6;
    unsigned maxIterations = 500;
    unsigned restart = 80;
    unsigned deflationVectors = 20;
    double deflationEps = 1e-3;
    bool collectNewtonSteps = true;
    bool reorthogonalizeOnMatrixChange = true;
    bool reorthogonalizeKrylov = true;
    std :: string preconditioner = "amgcl";
    bool verbose = false;
};

//////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////
// Hypre BoomerAMG solver options
struct HypreBoomerAMGOptions
{
    double tolerance = 1e-6;
    unsigned maxIterations = 500;
    int coarsenType = 8;       // HMIS
    int interpType = 6;        // extended+i
    double strongThreshold = 0.5;
    int nodal = 4;
    int numFunctions = 0;
    int relaxType = 6;         // symmetric hybrid Gauss-Seidel/Jacobi
    int relaxOrder = 0;
    int pMaxElmts = 4;
    int aggNumLevels = 0;
    int boomerMaxIterations = 1;
    int printLevel = 0;
    int skipBreak = 0;
    bool flex = false;
    bool recomputeResidual = false;
    int recomputeResidualPeriod = 0;
    double nonGalerkinTol = -1.;
    bool useDofFunctions = false;
    bool useInterpVectors = false;
    int interpVecVariant = 2;
    bool checkMatrix = false;
};

//////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////
// SOLVER FOR LINEAR ALGEBRA, CONJUGATE GRADIENTS
class ConjGradSolver : public LinAlgSolver
{
protected:
    Vector initialGuess;
    CoordinateIndexedSparseMatrix activeMatrix;
    double relMaxIT, maxIT, precision;
    long long lastIterations;
    double lastError;
    double lastTrueRelativeResidual;
    Eigen :: ConjugateGradient< Eigen :: SparseMatrix< double >, Eigen :: Lower | Eigen :: Upper, Eigen :: DiagonalPreconditioner< double > >cgK;
    //ConjugateGradient< SparseMatrix< double >, Lower | Upper, IncompleteCholesky< double > >cgK;
public:
    ConjGradSolver();
    virtual ~ConjGradSolver();
    virtual bool analyzePattern(const CoordinateIndexedSparseMatrix &A);
    virtual bool factorize(const CoordinateIndexedSparseMatrix &A);
    virtual bool solve(Vector &x, const Vector &b);
    void setPrecisionAndRelMaxIters(double p, double rmi);
    virtual long long giveLastIterations() const { return lastIterations; };
    virtual double giveLastError() const { return lastError; };
    virtual double giveLastTrueRelativeResidual() const { return lastTrueRelativeResidual; };
    virtual void setRuntimeTolerance(double tolerance, double trueTolerance);
    virtual double giveRuntimeTolerance() const { return precision; };
    virtual double giveRuntimeTrueTolerance() const { return precision; };
};

//////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////
// SOLVER FOR LINEAR ALGEBRA, AMGCL-PRECONDITIONED CONJUGATE GRADIENTS
#ifdef AMGCL_FOUND
class AmgclCGSolver : public LinAlgSolver
{
protected:
    struct Impl;
    std :: unique_ptr< Impl >impl;
    AmgclSolverOptions options;
    ElasticDofMap elasticMap;
    long long lastIterations;
    double lastError;
    double lastTrueRelativeResidual;
public:
    AmgclCGSolver();
    virtual ~AmgclCGSolver();
    virtual bool analyzePattern(const CoordinateIndexedSparseMatrix &A);
    virtual bool factorize(const CoordinateIndexedSparseMatrix &A);
    virtual bool solve(Vector &x, const Vector &b);
    void setOptions(const AmgclSolverOptions &opts);
    void setElasticDofMap(ElasticDofMap map);
    virtual long long giveLastIterations() const { return lastIterations; };
    virtual double giveLastError() const { return lastError; };
    virtual double giveLastTrueRelativeResidual() const { return lastTrueRelativeResidual; };
    virtual void setRuntimeTolerance(double tolerance, double trueTolerance);
    virtual double giveRuntimeTolerance() const { return options.tolerance; };
    virtual double giveRuntimeTrueTolerance() const { return options.trueTolerance > 0. ? options.trueTolerance : options.tolerance; };
};
#endif

//////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////
// SOLVER FOR LINEAR ALGEBRA, DEFLATED FLEXIBLE GMRES WITH AMGCL PRECONDITIONER
#ifdef AMGCL_FOUND
class DeflatedFGMRESSolver : public LinAlgSolver
{
protected:
    struct Impl;
    std :: unique_ptr< Impl >impl;
    DeflatedFGMRESOptions options;
    AmgclSolverOptions amgclOptions;
    HypreBoomerAMGOptions hypreOptions;
    ElasticDofMap elasticMap;
    long long lastIterations;
    double lastError;
    double lastTrueRelativeResidual;
    double lastPreconditionerApplySeconds;
    double lastOrthogonalizationSeconds;
    double lastLeastSquaresSeconds;
    double lastMatvecSeconds;
    double lastDeflationSeconds;
    Vector reducedVectorToActiveUnknown(const Vector &v) const;
    bool appendRawDeflationVector(const Vector &v);
    bool appendActiveDeflationVector(const Vector &rawReducedVector, bool storeRawVector);
    void updateDeflationOrthogonalityDiagnostics();
    void rebuildDeflationBasis();
    Vector projectDeflation(const Vector &v) const;
    Vector applyPreconditioner(const Vector &v);
public:
    DeflatedFGMRESSolver();
    virtual ~DeflatedFGMRESSolver();
    virtual bool analyzePattern(const CoordinateIndexedSparseMatrix &A);
    virtual bool factorize(const CoordinateIndexedSparseMatrix &A);
    virtual bool solve(Vector &x, const Vector &b);
    void setOptions(const DeflatedFGMRESOptions &opts);
    void setAmgclOptions(const AmgclSolverOptions &opts);
    void setHypreOptions(const HypreBoomerAMGOptions &opts);
    void setElasticDofMap(ElasticDofMap map);
    virtual void collectDeflationVector(const Vector &v);
    virtual long long giveLastIterations() const { return lastIterations; };
    virtual double giveLastError() const { return lastError; };
    virtual double giveLastTrueRelativeResidual() const { return lastTrueRelativeResidual; };
    virtual void setRuntimeTolerance(double tolerance, double trueTolerance);
    virtual double giveRuntimeTolerance() const { return options.tolerance; };
    virtual double giveRuntimeTrueTolerance() const { return options.trueTolerance > 0. ? options.trueTolerance : options.tolerance; };
    virtual unsigned giveDeflationBasisSize() const;
    virtual unsigned giveDeflationDiscardedCount() const;
    virtual unsigned giveDeflationRawCandidateCount() const;
    virtual unsigned giveDeflationCapacityEvictionCount() const;
    virtual unsigned giveDeflationLowANormDiscardCount() const;
    virtual double giveDeflationLastInitialANorm() const;
    virtual double giveDeflationLastFinalANorm() const;
    virtual double giveDeflationOrthogonalityMaxOffdiag() const;
    virtual double giveDeflationOrthogonalityMaxDiagError() const;
    virtual std :: string giveDeflationLastDiscardReason() const;
    virtual double giveLastPreconditionerApplySeconds() const { return lastPreconditionerApplySeconds; };
    virtual double giveLastOrthogonalizationSeconds() const { return lastOrthogonalizationSeconds; };
    virtual double giveLastLeastSquaresSeconds() const { return lastLeastSquaresSeconds; };
    virtual double giveLastMatvecSeconds() const { return lastMatvecSeconds; };
    virtual double giveLastDeflationSeconds() const { return lastDeflationSeconds; };
};
#endif

//////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////
// SOLVER FOR LINEAR ALGEBRA, HYPRE BOOMERAMG-PRECONDITIONED CG
#ifdef HYPRE_FOUND
class HypreBoomerAMGCGSolver : public LinAlgSolver
{
protected:
    struct Impl;
    std :: unique_ptr< Impl >impl;
    HypreBoomerAMGOptions options;
    ElasticDofMap elasticMap;
    long long lastIterations;
    double lastError;
    double lastTrueRelativeResidual;
public:
    HypreBoomerAMGCGSolver();
    virtual ~HypreBoomerAMGCGSolver();
    virtual bool analyzePattern(const CoordinateIndexedSparseMatrix &A);
    virtual bool factorize(const CoordinateIndexedSparseMatrix &A);
    virtual bool solve(Vector &x, const Vector &b);
    void setOptions(const HypreBoomerAMGOptions &opts);
    void setElasticDofMap(ElasticDofMap map);
    virtual long long giveLastIterations() const { return lastIterations; };
    virtual double giveLastError() const { return lastError; };
    virtual double giveLastTrueRelativeResidual() const { return lastTrueRelativeResidual; };
    virtual void setRuntimeTolerance(double tolerance, double trueTolerance);
    virtual double giveRuntimeTolerance() const { return options.tolerance; };
    virtual double giveRuntimeTrueTolerance() const { return options.tolerance; };
};
#endif

//////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////
// SOLVER FOR LINEAR ALGEBRA, LDLT SOLVER
class LDLTSolver : public LinAlgSolver
{
protected:
    Eigen :: SimplicialLDLT< Eigen :: SparseMatrix< double > >ldlt;
public:
    LDLTSolver();
    virtual ~LDLTSolver();
    virtual bool analyzePattern(const CoordinateIndexedSparseMatrix &A);
    virtual bool factorize(const CoordinateIndexedSparseMatrix &A);
    virtual bool solve(Vector &x, const Vector &b);
};

//////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////
// SOLVER FOR LINEAR ALGEBRA, LDLT SOLVER
class LLTSolver : public LinAlgSolver
{
protected:
    Eigen :: SimplicialLLT< Eigen :: SparseMatrix< double > >llt;
public:
    LLTSolver();
    virtual ~LLTSolver();
    virtual bool analyzePattern(const CoordinateIndexedSparseMatrix &A);
    virtual bool factorize(const CoordinateIndexedSparseMatrix &A);
    virtual bool solve(Vector &x, const Vector &b);
};
//////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////
// SOLVER FOR LINEAR ALGEBRA, LU SOLVER
class LUSolver : public LinAlgSolver
{
protected:
    Eigen :: SparseLU< Eigen :: SparseMatrix< double >, Eigen :: COLAMDOrdering< int > >lu;
public:
    LUSolver();
    virtual ~LUSolver();
    virtual bool analyzePattern(const CoordinateIndexedSparseMatrix &A);
    virtual bool factorize(const CoordinateIndexedSparseMatrix &A);
    virtual bool solve(Vector &x, const Vector &b);
};
//////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////
// SOLVER FOR LINEAR ALGEBRA, SUPERLU SOLVER
#ifdef SUPERLU_FOUND
class SuperLUSolver : public LinAlgSolver
{
protected:
    Eigen :: SuperLU< Eigen :: SparseMatrix< double > >superlu;
public:
    SuperLUSolver();
    virtual ~SuperLUSolver();
    virtual bool analyzePattern(const CoordinateIndexedSparseMatrix &A);
    virtual bool factorize(const CoordinateIndexedSparseMatrix &A);
    virtual bool solve(Vector &x, const Vector &b);
};
#endif

//////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////
// SOLVER FOR LINEAR ALGEBRA, PARDISO SOLVERS
#ifdef PARDISO_FOUND
class PardisoLUSolver : public LinAlgSolver
{
protected:
    Eigen :: PardisoLU< Eigen :: SparseMatrix< double > >pardiso;
public:
    PardisoLUSolver();
    virtual ~PardisoLUSolver();
    virtual bool analyzePattern(const CoordinateIndexedSparseMatrix &A);
    virtual bool factorize(const CoordinateIndexedSparseMatrix &A);
    virtual bool solve(Vector &x, const Vector &b);
};

class PardisoLDLTSolver : public LinAlgSolver
{
protected:
    Eigen :: PardisoLDLT< Eigen :: SparseMatrix< double > >pardiso;
public:
    PardisoLDLTSolver();
    virtual ~PardisoLDLTSolver();
    virtual bool analyzePattern(const CoordinateIndexedSparseMatrix &A);
    virtual bool factorize(const CoordinateIndexedSparseMatrix &A);
    virtual bool solve(Vector &x, const Vector &b);
};

class PardisoLLTSolver : public LinAlgSolver
{
protected:
    Eigen :: PardisoLLT< Eigen :: SparseMatrix< double > >pardiso;
public:
    PardisoLLTSolver();
    virtual ~PardisoLLTSolver();
    virtual bool analyzePattern(const CoordinateIndexedSparseMatrix &A);
    virtual bool factorize(const CoordinateIndexedSparseMatrix &A);
    virtual bool solve(Vector &x, const Vector &b);
};
#endif

//////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////
// SOLVER FOR LINEAR ALGEBRA, CHOLMOD SOLVERS
#ifdef CHOLMOD_FOUND
class CholmodLLTSolver : public LinAlgSolver
{
protected:
    Eigen :: CholmodSimplicialLLT< Eigen :: SparseMatrix< double > >cholmod;
public:
    CholmodLLTSolver();
    virtual ~CholmodLLTSolver();
    virtual bool analyzePattern(const CoordinateIndexedSparseMatrix &A);
    virtual bool factorize(const CoordinateIndexedSparseMatrix &A);
    virtual bool solve(Vector &x, const Vector &b);
};

class CholmodLDLTSolver : public LinAlgSolver
{
protected:
    Eigen :: CholmodSimplicialLDLT< Eigen :: SparseMatrix< double > >cholmod;
public:
    CholmodLDLTSolver();
    virtual ~CholmodLDLTSolver();
    virtual bool analyzePattern(const CoordinateIndexedSparseMatrix &A);
    virtual bool factorize(const CoordinateIndexedSparseMatrix &A);
    virtual bool solve(Vector &x, const Vector &b);
};

class CholmodSupernodalLLTSolver : public LinAlgSolver
{
protected:
    Eigen :: CholmodSupernodalLLT< Eigen :: SparseMatrix< double > >cholmod;
public:
    CholmodSupernodalLLTSolver();
    virtual ~CholmodSupernodalLLTSolver();
    virtual bool analyzePattern(const CoordinateIndexedSparseMatrix &A);
    virtual bool factorize(const CoordinateIndexedSparseMatrix &A);
    virtual bool solve(Vector &x, const Vector &b);
};
#endif

//////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////
bool StandAloneLinalgSolver(const CoordinateIndexedSparseMatrix &A, Vector &ddr, const Vector &f, double precision, double relmaxit, std::string solver_type);

bool LinalgNonSymmetricSolver(const CoordinateIndexedSparseMatrix &A, Vector &x, const Vector &b, const Vector x0, double precision, double relmaxit);

bool LinalgEigenSolver(const Vector &A, Vector &eigenvalues, std :: vector< Vector > &eigevectors);

bool LinalgEigenSolver(const Matrix &mat, Vector &eigenvalues, std :: vector< Vector > &eigevectors);

bool LinalgEigenSpectraSolver(const CoordinateIndexedSparseMatrix &mat, Vector &eigenvalues, Matrix &eigenvectors, int n_eigen_vals);

bool LinalgEigenSpectraGENSolver(const CoordinateIndexedSparseMatrix &mat, const CoordinateIndexedSparseMatrix &matB, Vector &eigenvalues, Matrix &eigenvectors, int n_eigen_vals);

bool LinalgLUSolver(const CoordinateIndexedSparseMatrix &A, Vector &x, const Vector &b);

double checkCoplanarity(const Point &ptA, const Point &ptB, const Point &ptC, const Point &ptD);

Matrix dyadicProduct(const Vector &a, const Vector &b);

double triArea2D(const Point *a, const Point *b, const Point *c);

double triArea3D(const Point *a, const Point *b, const Point *c);

double find_intesection_of_segment_and_triangle(const Point *A, const Point *B, const Point *a, const Point *b, const Point *c);

double triInertia2D(const Point *a, const Point *b, const Point *c);

double tetraVolumeSigned(const Point *a, const Point *b, const Point *c, const Point *d);

Matrix tetraInertia3D(const Point *a, const Point *b, const Point *c, const Point *d);

bool is_positive_integer(const std :: string &s);

void giveGaussIntegrationPointAndWeights(unsigned n, Vector &locs, Vector &weis);
#endif /* _LINALG_H */
