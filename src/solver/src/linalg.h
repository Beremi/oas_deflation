
#ifndef LINALG_H
#define LINALG_H

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

//  \[ ([\w\s\+\*]*) \] \[ ([\w\s\+\*]*) \]
#include "globals.h"

using Ttripletd = Eigen :: Triplet< double >;

//typedef typename Eigen :: Matrix<double, Eigen :: Dynamic, 1> MyVector;
//typedef Eigen :: MyMatrix<double, 3, 1> Point;
using Vector = Eigen :: VectorXd;
using Point = Eigen :: Vector3d;
using Matrix = Eigen :: MatrixXd;
//typedef typename Eigen :: MyMatrix<double, Eigen :: Dynamic, Eigen :: Dynamic> MyMatrix;
typedef typename Eigen :: SparseMatrix< double, Eigen :: RowMajor > CoordinateIndexedSparseMatrix;

//Eigen :: IOFormat CommaInitFmt(Eigen :: StreamPrecision, Eigen :: DontAlignCols, ", ", ", ", "", "", " << ", "");
const static Eigen :: IOFormat SemicolonInitFmt(Eigen :: FullPrecision, Eigen :: DontAlignCols, "; ", "; ", "", "", "", "");
//Eigen :: IOFormat CleanFmt(4, 0, ", ", "\n", "[", "]");
//Eigen :: IOFormat OctaveFmt(Eigen :: StreamPrecision, 0, ", ", ";\n", "", "", "[", "]");
//Eigen :: IOFormat HeavyFmt(Eigen :: FullPrecision, 0, ", ", ";\n", "[", "]", "[", "]");

//bool LinalgSymmetricSolver(const CoordinateIndexedSparseMatrix &A, Eigen :: Ref<MyVector> x, const Eigen::Ref<const MyVector> b, const Eigen::Ref<const MyVector> x0, double precision, double relmaxit, string solver_type);
bool LinalgSymmetricSolver(const CoordinateIndexedSparseMatrix &A, Vector &x, const Vector &b, const Vector &x0, double precision, double relmaxit, std :: string solver_type);

bool LinalgNonSymmetricSolver(const CoordinateIndexedSparseMatrix &A, Vector &x, const Vector &b, const Vector x0, double precision, double relmaxit);

bool LinalgEigenSolver(const Vector &A, Vector &eigenvalues, std :: vector< Vector > &eigevalues);

bool LinalgEigenSolver(const Matrix& mat, Vector& eigenvalues, std::vector< Vector > &eigevalues);

double checkCoplanarity(const Point &ptA, const Point &ptB, const Point &ptC, const Point &ptD);

Matrix dyadicProduct(const Vector &a, const Vector &b);

double triArea2D(const Point *a, const Point *b, const Point *c);

double triArea3D(const Point *a, const Point *b, const Point *c);



#endif
