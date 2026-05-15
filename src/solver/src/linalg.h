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
    virtual void collectDeflationVector(const Vector &v) { ( void ) v; };
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
        return fullRows > 0 && blockSize > 0 && !reducedToFull.empty();
    }
};

//////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////
// Deflated flexible GMRES options
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
    bool elasticFullLift = true;
    int elasticReorder = 0;      // 0 = none, 1 = node-major, 2 = coordinate-sorted node-major
    std :: string preconditioner = "hypre";
    bool verbose = false;
};

//////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////
// hypre BoomerAMG options used as a DFGMRES preconditioner
struct HypreBoomerAMGOptions
{
    double tolerance = 1e-6;
    unsigned maxIterations = 500;
    int coarsenType = 8;
    int interpType = 6;
    double strongThreshold = 0.5;
    int nodal = 4;
    int nodalDiag = 0;
    int numFunctions = 0;
    int relaxType = 6;
    int relaxOrder = 0;
    int numSweeps = 0;
    int pMaxElmts = 4;
    int aggNumLevels = 0;
    int boomerMaxIterations = 1;
    int chebyOrder = 0;
    double chebyFraction = -1.;
    int elasticReorder = 0;      // 0 = none, 1 = node-major, 2 = coordinate-sorted node-major
    int printLevel = 0;
    double nonGalerkinTol = -1.;
    bool useDofFunctions = false;
    bool useInterpVectors = false;
    int interpVecVariant = 2;
    bool checkMatrix = false;
    int threads = 0;             // 0 = use the current OpenMP thread count
};

//////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////
// SOLVER FOR LINEAR ALGEBRA, DEFLATED FLEXIBLE GMRES
#ifdef HYPRE_FOUND
class DeflatedFGMRESSolver : public LinAlgSolver
{
protected:
    struct Impl;
    std :: unique_ptr< Impl >impl;
    DeflatedFGMRESOptions options;
    HypreBoomerAMGOptions hypreOptions;
    ElasticDofMap elasticMap;
    long long lastIterations;
    double lastError;
    double lastTrueRelativeResidual;

    Vector reducedVectorToActiveUnknown(const Vector &v) const;
    bool appendRawDeflationVector(const Vector &v);
    bool appendActiveDeflationVector(const Vector &rawReducedVector, bool storeRawVector);
    void updateDeflationOrthogonalityDiagnostics();
    void rebuildDeflationBasis();
    Vector projectDeflation(const Vector &v) const;
    Vector applyPreconditioner(const Vector &v);
    Vector activeMatvec(const Vector &v) const;

public:
    DeflatedFGMRESSolver();
    virtual ~DeflatedFGMRESSolver();
    virtual bool analyzePattern(const CoordinateIndexedSparseMatrix &A);
    virtual bool factorize(const CoordinateIndexedSparseMatrix &A);
    virtual bool solve(Vector &x, const Vector &b);
    void setOptions(const DeflatedFGMRESOptions &opts);
    void setHypreOptions(const HypreBoomerAMGOptions &opts);
    void setElasticDofMap(ElasticDofMap map);
    virtual void collectDeflationVector(const Vector &v);
    virtual long long giveLastIterations() const { return lastIterations; };
    virtual double giveLastError() const { return lastError; };
    virtual double giveLastTrueRelativeResidual() const { return lastTrueRelativeResidual; };
};
#endif

//////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////
// SOLVER FOR LINEAR ALGEBRA, CONJUGATE GRADIENTS
class ConjGradSolver : public LinAlgSolver
{
protected:
    Vector initialGuess;
    double relMaxIT, maxIT, precision;
    Eigen :: ConjugateGradient< Eigen :: SparseMatrix< double >, Eigen :: Lower | Eigen :: Upper, Eigen :: DiagonalPreconditioner< double > >cgK;
    //ConjugateGradient< SparseMatrix< double >, Lower | Upper, IncompleteCholesky< double > >cgK;
public:
    ConjGradSolver();
    virtual ~ConjGradSolver();
    virtual bool analyzePattern(const CoordinateIndexedSparseMatrix &A);
    virtual bool factorize(const CoordinateIndexedSparseMatrix &A);
    virtual bool solve(Vector &x, const Vector &b);
    void setPrecisionAndRelMaxIters(double p, double rmi);
};

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
