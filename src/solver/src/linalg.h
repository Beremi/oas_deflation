
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

#include "linear_algebra.h"
#include "globals.h"


bool LinalgSymmetricSolver(const CoordinateIndexedSparseMatrix &A, Vector &x, const Vector &b, const Vector x0, double precision, double relmaxit, string solver_type);


bool LinalgNonSymmetricSolver(const CoordinateIndexedSparseMatrix &A, Vector &x, const Vector &b, const Vector x0, double precision, double relmaxit);

#endif
