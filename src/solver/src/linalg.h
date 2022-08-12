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

//  \[ ([\w\s\+\*]*) \] \[ ([\w\s\+\*]*) \]
#include "globals.h"

using Ttripletd = Eigen :: Triplet< double >;

using Vector = Eigen :: VectorXd;
using Point = Eigen :: Vector3d;
using Matrix = Eigen :: MatrixXd;
typedef typename Eigen :: SparseMatrix< double, Eigen :: RowMajor >CoordinateIndexedSparseMatrix;  // row-major-sparse * dense vector/matrix products - multi-threading

const static Eigen :: IOFormat VectorSemicolonFmt(Eigen :: FullPrecision, Eigen :: DontAlignCols, "; ", "; ", "", "", "", "");

bool LinalgSymmetricSolver(const CoordinateIndexedSparseMatrix &A, Vector &x, const Vector &b, const Vector &x0, double precision, double relmaxit, std :: string solver_type);

bool LinalgNonSymmetricSolver(const CoordinateIndexedSparseMatrix &A, Vector &x, const Vector &b, const Vector x0, double precision, double relmaxit);

bool LinalgEigenSolver(const Vector &A, Vector &eigenvalues, std :: vector< Vector > &eigevalues);

bool LinalgEigenSolver(const Matrix &mat, Vector &eigenvalues, std :: vector< Vector > &eigevalues);

double checkCoplanarity(const Point &ptA, const Point &ptB, const Point &ptC, const Point &ptD);

Matrix dyadicProduct(const Vector &a, const Vector &b);

double triArea2D(const Point *a, const Point *b, const Point *c);

double triArea3D(const Point *a, const Point *b, const Point *c);

double triInertia2D(const Point *a, const Point *b, const Point *c);

double tetraVolumeSigned(const Point *a, const Point *b, const Point *c, const Point *d);

Matrix tetraInertia3D(const Point *a, const Point *b, const Point *c, const Point *d);

#endif /* _LINALG_H */
