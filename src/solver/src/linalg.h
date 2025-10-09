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
    virtual bool analyzePattern(const CoordinateIndexedSparseMatrix &A){ return false; };
    virtual bool factorize(const CoordinateIndexedSparseMatrix &A) { ( void ) A; name = "null solver, base class"; return false; };
    virtual bool solve(Vector &x, const Vector &b) { ( void ) x; ( void ) b; return false; };
    std :: string giveName()const { return name; };
};

//////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////
// SOLVER FOR LINEAR ALGEBRA, CONJUGATE GRADIENTS
class ConjGradSolver : public LinAlgSolver
{
protected:
    Vector initialGuess;
    double relMaxIT, maxIT, precision;
    Eigen :: ConjugateGradient< Eigen :: SparseMatrix< double >, Eigen :: Lower | Eigen :: Upper >cgK;
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
bool LinalgSymmetricSolver(const CoordinateIndexedSparseMatrix &A, Vector &x, const Vector &b, const Vector &x0, double precision, double relmaxit, std :: string solver_type);

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
