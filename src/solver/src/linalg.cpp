#include "linalg.h"

#ifdef PARDISO_FOUND
 #include <mkl.h>     // For mkl_set_num_threads
#endif

#ifdef AMGCL_FOUND
 #include <algorithm>
 #include <cctype>
 #include <cmath>
 #include <deque>
 #include <limits>
 #include <tuple>
 #include <amgcl/backend/builtin.hpp>
 #include <amgcl/backend/builtin_hybrid.hpp>
 #include <amgcl/adapter/crs_tuple.hpp>
 #include <amgcl/adapter/block_matrix.hpp>
 #include <amgcl/make_solver.hpp>
 #include <amgcl/amg.hpp>
 #include <amgcl/coarsening/smoothed_aggregation.hpp>
 #include <amgcl/coarsening/smoothed_aggr_emin.hpp>
 #include <amgcl/coarsening/as_scalar.hpp>
 #include <amgcl/relaxation/ilu0.hpp>
 #include <amgcl/relaxation/spai0.hpp>
 #include <amgcl/solver/cg.hpp>
 #include <amgcl/solver/fgmres.hpp>
 #include <amgcl/value_type/static_matrix.hpp>
#endif

#ifdef HYPRE_FOUND
 #include <HYPRE.h>
 #include <HYPRE_IJ_mv.h>
 #include <HYPRE_krylov.h>
 #include <HYPRE_parcsr_ls.h>
 #include <HYPRE_utilities.h>
 #include <mpi.h>
#endif

using namespace std;

namespace {

double trueRelativeResidual(const CoordinateIndexedSparseMatrix &A, const Vector &x, const Vector &b) {
    if ( b.size() == 0 ) {
        return 0.;
    }
    const double bNorm = std :: max(b.norm(), 1e-300);
    return ( A * x - b ).norm() / bNorm;
}

void makeRowMajorCrs(
    const CoordinateIndexedSparseMatrix &A,
    ptrdiff_t &rows,
    std :: vector< ptrdiff_t > &ptr,
    std :: vector< ptrdiff_t > &col,
    std :: vector< double > &val
) {
    Eigen :: SparseMatrix< double, Eigen :: RowMajor, ptrdiff_t >rowMajor = A;
    rowMajor.makeCompressed();
    rows = rowMajor.rows();
    ptr.assign(rowMajor.outerIndexPtr(), rowMajor.outerIndexPtr() + rowMajor.outerSize() + 1);
    col.assign(rowMajor.innerIndexPtr(), rowMajor.innerIndexPtr() + rowMajor.nonZeros());
    val.assign(rowMajor.valuePtr(), rowMajor.valuePtr() + rowMajor.nonZeros());
}

CoordinateIndexedSparseMatrix liftMatrixToElasticSpace(const CoordinateIndexedSparseMatrix &A, const ElasticDofMap &map) {
    if ( !map.isValid() ) {
        return A;
    }
    std :: vector< Ttripletd >triplets;
    triplets.reserve(A.nonZeros() + map.fullRows);
    std :: vector< char >hasReducedRow(map.fullRows, 0);
    for ( unsigned row = 0; row < map.reducedToFull.size(); row++ ) {
        if ( map.reducedToFull [ row ] < map.fullRows ) {
            hasReducedRow [ map.reducedToFull [ row ] ] = 1;
        }
    }
    for ( int outer = 0; outer < A.outerSize(); outer++ ) {
        for ( CoordinateIndexedSparseMatrix :: InnerIterator it(A, outer); it; ++it ) {
            const unsigned fullRow = map.reducedToFull [ it.row() ];
            const unsigned fullCol = map.reducedToFull [ it.col() ];
            triplets.emplace_back(fullRow, fullCol, it.value() );
        }
    }
    for ( unsigned row = 0; row < map.fullRows; row++ ) {
        if ( !hasReducedRow [ row ] ) {
            triplets.emplace_back(row, row, 1.);
        }
    }
    CoordinateIndexedSparseMatrix lifted(map.fullRows, map.fullRows);
    lifted.setFromTriplets(triplets.begin(), triplets.end());
    lifted.makeCompressed();
    return lifted;
}

CoordinateIndexedSparseMatrix symmetricallyScaleMatrix(const CoordinateIndexedSparseMatrix &A, const std :: vector< double > &scale) {
    std :: vector< Ttripletd >triplets;
    triplets.reserve(A.nonZeros() );
    for ( int outer = 0; outer < A.outerSize(); outer++ ) {
        for ( CoordinateIndexedSparseMatrix :: InnerIterator it(A, outer); it; ++it ) {
            triplets.emplace_back(it.row(), it.col(), it.value() * scale [ it.row() ] * scale [ it.col() ]);
        }
    }
    CoordinateIndexedSparseMatrix scaled(A.rows(), A.cols() );
    scaled.setFromTriplets(triplets.begin(), triplets.end());
    scaled.makeCompressed();
    return scaled;
}

Vector liftVectorToElasticSpace(const Vector &v, const ElasticDofMap &map) {
    if ( !map.isValid() ) {
        return v;
    }
    Vector full = Vector :: Zero(map.fullRows);
    for ( unsigned row = 0; row < map.reducedToFull.size(); row++ ) {
        full [ map.reducedToFull [ row ] ] = v [ row ];
    }
    return full;
}

Vector restrictVectorFromElasticSpace(const Vector &v, const ElasticDofMap &map) {
    if ( !map.isValid() ) {
        return v;
    }
    Vector reduced = Vector :: Zero(map.reducedToFull.size());
    for ( unsigned row = 0; row < map.reducedToFull.size(); row++ ) {
        reduced [ row ] = v [ map.reducedToFull [ row ] ];
    }
    return reduced;
}

std :: string lowerCopy(std :: string value) {
    std :: transform(value.begin(), value.end(), value.begin(), [](unsigned char c) {
        return static_cast< char >(std :: tolower(c) );
    });
    return value;
}

}

#ifdef HYPRE_FOUND
namespace {
void ensureHypreRuntimeInitialized();
bool hypreCheck(HYPRE_Int error, const char *operation);
bool makeHypreVector(
    const std :: vector< HYPRE_BigInt > &indices,
    const std :: vector< double > &values,
    HYPRE_IJVector &ijVector,
    HYPRE_ParVector &parVector
);
}
#endif

#ifdef AMGCL_FOUND
struct AmgclCGSolver :: Impl
{
    using ScalarBackend = amgcl :: backend :: builtin< double >;
    using ScalarSolver = amgcl :: make_solver<
        amgcl :: amg<
            ScalarBackend,
            amgcl :: coarsening :: smoothed_aggregation,
            amgcl :: relaxation :: spai0
        >,
        amgcl :: solver :: cg< ScalarBackend >
    >;
    using Block3 = amgcl :: static_matrix< double, 3, 3 >;
    using Block6 = amgcl :: static_matrix< double, 6, 6 >;
    using Block3Backend = amgcl :: backend :: builtin< Block3 >;
    using Block6Backend = amgcl :: backend :: builtin< Block6 >;
    using Hybrid3Backend = amgcl :: backend :: builtin_hybrid< Block3 >;
    using Hybrid6Backend = amgcl :: backend :: builtin_hybrid< Block6 >;
    using Hybrid3Solver = amgcl :: make_solver<
        amgcl :: amg<
            Hybrid3Backend,
            amgcl :: coarsening :: smoothed_aggregation,
            amgcl :: relaxation :: spai0
        >,
        amgcl :: solver :: cg< Hybrid3Backend >
    >;
    using Hybrid6Solver = amgcl :: make_solver<
        amgcl :: amg<
            Hybrid6Backend,
            amgcl :: coarsening :: smoothed_aggregation,
            amgcl :: relaxation :: spai0
        >,
        amgcl :: solver :: cg< Hybrid6Backend >
    >;
    using Hybrid3FgmresSolver = amgcl :: make_solver<
        amgcl :: amg<
            Hybrid3Backend,
            amgcl :: coarsening :: smoothed_aggregation,
            amgcl :: relaxation :: spai0
        >,
        amgcl :: solver :: fgmres< Hybrid3Backend >
    >;
    using Hybrid6FgmresSolver = amgcl :: make_solver<
        amgcl :: amg<
            Hybrid6Backend,
            amgcl :: coarsening :: smoothed_aggregation,
            amgcl :: relaxation :: spai0
        >,
        amgcl :: solver :: fgmres< Hybrid6Backend >
    >;
    using Hybrid3EminSolver = amgcl :: make_solver<
        amgcl :: amg<
            Hybrid3Backend,
            amgcl :: coarsening :: smoothed_aggr_emin,
            amgcl :: relaxation :: spai0
        >,
        amgcl :: solver :: cg< Hybrid3Backend >
    >;
    using Hybrid6EminSolver = amgcl :: make_solver<
        amgcl :: amg<
            Hybrid6Backend,
            amgcl :: coarsening :: smoothed_aggr_emin,
            amgcl :: relaxation :: spai0
        >,
        amgcl :: solver :: cg< Hybrid6Backend >
    >;
    using Block3IluSolver = amgcl :: make_solver<
        amgcl :: amg<
            Block3Backend,
            amgcl :: coarsening :: as_scalar< amgcl :: coarsening :: smoothed_aggregation > :: type,
            amgcl :: relaxation :: ilu0
        >,
        amgcl :: solver :: cg< Block3Backend >
    >;
    using Block6IluSolver = amgcl :: make_solver<
        amgcl :: amg<
            Block6Backend,
            amgcl :: coarsening :: as_scalar< amgcl :: coarsening :: smoothed_aggregation > :: type,
            amgcl :: relaxation :: ilu0
        >,
        amgcl :: solver :: cg< Block6Backend >
    >;
    using Block3SpaiSolver = amgcl :: make_solver<
        amgcl :: amg<
            Block3Backend,
            amgcl :: coarsening :: as_scalar< amgcl :: coarsening :: smoothed_aggregation > :: type,
            amgcl :: relaxation :: spai0
        >,
        amgcl :: solver :: cg< Block3Backend >
    >;
    using Block6SpaiSolver = amgcl :: make_solver<
        amgcl :: amg<
            Block6Backend,
            amgcl :: coarsening :: as_scalar< amgcl :: coarsening :: smoothed_aggregation > :: type,
            amgcl :: relaxation :: spai0
        >,
        amgcl :: solver :: cg< Block6Backend >
    >;

    ptrdiff_t rows = 0;
    std :: vector< ptrdiff_t >ptr;
    std :: vector< ptrdiff_t >col;
    std :: vector< double >val;
    CoordinateIndexedSparseMatrix reducedMatrix;
    CoordinateIndexedSparseMatrix activeMatrix;
    std :: vector< double >scale;
    std :: unique_ptr< ScalarSolver >scalarSolver;
    std :: unique_ptr< Block3IluSolver >block3IluSolver;
    std :: unique_ptr< Block6IluSolver >block6IluSolver;
    std :: unique_ptr< Block3SpaiSolver >block3SpaiSolver;
    std :: unique_ptr< Block6SpaiSolver >block6SpaiSolver;
    std :: unique_ptr< Hybrid3Solver >hybrid3Solver;
    std :: unique_ptr< Hybrid6Solver >hybrid6Solver;
    std :: unique_ptr< Hybrid3FgmresSolver >hybrid3FgmresSolver;
    std :: unique_ptr< Hybrid6FgmresSolver >hybrid6FgmresSolver;
    std :: unique_ptr< Hybrid3EminSolver >hybrid3EminSolver;
    std :: unique_ptr< Hybrid6EminSolver >hybrid6EminSolver;
    std :: vector< double >initialGuess;
    unsigned activeBlockSize = 1;
    std :: string activeBackend = "scalar";

    auto matrixTuple() {
        return std :: tie(rows, ptr, col, val);
    }
};
#endif

//////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////
// CONJUGATE GRADIENT SOLVER
//////////////////////////////////////////////////////////

ConjGradSolver :: ConjGradSolver() {
    name = "EigenConj";
    relMaxIT = 1.;
    precision = 1e-12;
    lastIterations = -1;
    lastError = -1.;
    lastTrueRelativeResidual = -1.;
}

//////////////////////////////////////////////////////////
ConjGradSolver :: ~ConjGradSolver() {}

//////////////////////////////////////////////////////////
void ConjGradSolver :: setPrecisionAndRelMaxIters(double p, double rmi) {
    relMaxIT = rmi;
    precision = p;
}

//////////////////////////////////////////////////////////
void ConjGradSolver :: setRuntimeTolerance(double tolerance, double trueTolerance) {
    ( void ) trueTolerance;
    if ( tolerance > 0. ) {
        precision = tolerance;
        cgK.setTolerance(precision);
    }
}

//////////////////////////////////////////////////////////
bool ConjGradSolver :: factorize(const CoordinateIndexedSparseMatrix &A) {
#if PRINT_DEBUG_TIME
    auto start = std :: chrono :: system_clock :: now();
#endif
    //JacobiSVD<MatrixXd> svd(cgb);
    //double cond = svd.singularValues()(0) / svd.singularValues()(svd.singularValues().size()-1);
    //cout << "condition number is " << cond<< " " << svd.singularValues()(0) << " " << svd.singularValues()(svd.singularValues().size()-1) << endl;

    if ( A.rows() > 0 ) {
        activeMatrix = A;
        cgK.setMaxIterations( relMaxIT * A.cols() );
        cgK.setTolerance(precision);
        cgK.factorize(A);
        initialGuess = Vector :: Zero( A.cols() );
        maxIT = relMaxIT * A.cols();
    }

#if PRINT_DEBUG_TIME
    now = std :: chrono :: system_clock :: now();

    elapsed_seconds = now - start;
    std :: cout << "linalg solver decomposition duration: " << convertTimeToString_(elapsed_seconds) << std :: endl;
    cout.flush();
#endif
    return true;
}

//////////////////////////////////////////////////////////
bool ConjGradSolver :: analyzePattern(const CoordinateIndexedSparseMatrix &A) {
#if PRINT_DEBUG_TIME
    auto start = std :: chrono :: system_clock :: now();
#endif

    if ( A.rows() > 0 ) {
        cgK.analyzePattern(A);
    }

#if PRINT_DEBUG_TIME
    now = std :: chrono :: system_clock :: now();

    elapsed_seconds = now - start;
    std :: cout << "linalg solver decomposition duration: " << convertTimeToString_(elapsed_seconds) << std :: endl;
    cout.flush();
#endif
    return true;
}



//////////////////////////////////////////////////////////
bool ConjGradSolver :: solve(Vector &x, const Vector &b) {
#if PRINT_DEBUG_TIME
    auto start = std :: chrono :: system_clock :: now();
#endif

    bool result = false;
    if ( b.size() > 0 ) {
        x = cgK.solveWithGuess(b, initialGuess);
        lastIterations = cgK.iterations();
        lastError = cgK.error();
        lastTrueRelativeResidual = trueRelativeResidual(activeMatrix, x, b);
        result = size_t( cgK.iterations() ) < maxIT;
        if ( !result ) {
            cerr << "Eigen Conjugate Gradients performed " << cgK.iterations() << " iterations and reached error " << cgK.error() << ", required precision is " << precision << endl;
            exit(1);
        }
        initialGuess = x;
    } else {
        lastIterations = 0;
        lastError = 0.;
        lastTrueRelativeResidual = 0.;
        result = true;
    }

#if PRINT_DEBUG_TIME
    now = std :: chrono :: system_clock :: now();

    elapsed_seconds = now - start;
    std :: cout << "linalg solver duration: " << convertTimeToString_(elapsed_seconds) << std :: endl;
    cout.flush();
#endif
    return result;
}

#ifdef AMGCL_FOUND
//////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////
// AMGCL-PRECONDITIONED CONJUGATE GRADIENT SOLVER
//////////////////////////////////////////////////////////

AmgclCGSolver :: AmgclCGSolver() {
    name = "AmgclCGElasticSolver";
    impl = std :: make_unique< Impl >();
    lastIterations = -1;
    lastError = -1.;
    lastTrueRelativeResidual = -1.;
}

//////////////////////////////////////////////////////////
AmgclCGSolver :: ~AmgclCGSolver() {}

//////////////////////////////////////////////////////////
void AmgclCGSolver :: setOptions(const AmgclSolverOptions &opts) {
    options = opts;
}

//////////////////////////////////////////////////////////
void AmgclCGSolver :: setRuntimeTolerance(double tolerance, double trueTolerance) {
    if ( tolerance > 0. ) {
        options.tolerance = tolerance;
    }
    if ( trueTolerance > 0. ) {
        options.trueTolerance = trueTolerance;
    } else if ( tolerance > 0. ) {
        options.trueTolerance = tolerance;
    }
}

//////////////////////////////////////////////////////////
void AmgclCGSolver :: setElasticDofMap(ElasticDofMap map) {
    elasticMap = std :: move(map);
}

//////////////////////////////////////////////////////////
bool AmgclCGSolver :: analyzePattern(const CoordinateIndexedSparseMatrix &A) {
    ( void ) A;
    return true;
}

//////////////////////////////////////////////////////////
bool AmgclCGSolver :: factorize(const CoordinateIndexedSparseMatrix &A) {
    if ( A.rows() <= 0 ) {
        return true;
    }

    const bool useElasticLift = options.elasticFullLift && elasticMap.isValid();
    impl->reducedMatrix = A;
    impl->activeMatrix = useElasticLift ? liftMatrixToElasticSpace(A, elasticMap) : A;
    impl->scale.assign(impl->activeMatrix.rows(), 1.);
    if ( options.diagonalScale ) {
        double minScale = std :: numeric_limits< double > :: infinity();
        double maxScale = 0.;
        long long nonPositiveDiagonal = 0;
        for ( int i = 0; i < impl->activeMatrix.rows(); ++i ) {
            const double diag = impl->activeMatrix.coeff(i, i);
            if ( diag > 0. && std :: isfinite(diag) ) {
                impl->scale [ i ] = 1. / sqrt(diag);
                minScale = std :: min(minScale, impl->scale [ i ]);
                maxScale = std :: max(maxScale, impl->scale [ i ]);
            } else {
                ++nonPositiveDiagonal;
            }
        }
        impl->activeMatrix = symmetricallyScaleMatrix(impl->activeMatrix, impl->scale);
        std :: cout << "AMGCL-CG: symmetric diagonal scaling enabled, min_scale=" << minScale
                    << ", max_scale=" << maxScale
                    << ", non_positive_diagonal_count=" << nonPositiveDiagonal << std :: endl;
    }
    makeRowMajorCrs(impl->activeMatrix, impl->rows, impl->ptr, impl->col, impl->val);

    if ( options.checkMatrix ) {
        CoordinateIndexedSparseMatrix diff = impl->activeMatrix - CoordinateIndexedSparseMatrix( impl->activeMatrix.transpose() );
        const double matrixNorm = impl->activeMatrix.norm();
        const double symmetryRelativeNorm = matrixNorm > 0. ? diff.norm() / matrixNorm : 0.;
        long long nonPositiveDiagonal = 0;
        double minDiagonal = std :: numeric_limits< double > :: infinity();
        double maxDiagonal = -std :: numeric_limits< double > :: infinity();
        for ( int i = 0; i < impl->activeMatrix.rows(); ++i ) {
            const double value = impl->activeMatrix.coeff(i, i);
            minDiagonal = std :: min(minDiagonal, value);
            maxDiagonal = std :: max(maxDiagonal, value);
            if ( value <= 0. ) {
                ++nonPositiveDiagonal;
            }
        }
        std :: cout << "AMGCL-CG matrix check: symmetry_relative_norm=" << symmetryRelativeNorm
                    << ", min_diagonal=" << minDiagonal
                    << ", max_diagonal=" << maxDiagonal
                    << ", non_positive_diagonal_count=" << nonPositiveDiagonal << std :: endl;
    }

    Impl :: ScalarSolver :: params params;
    params.solver.tol = options.tolerance;
    params.solver.maxiter = options.maxIterations;
    params.solver.verbose = options.verbose;
    params.precond.coarsening.aggr.eps_strong = options.epsStrong;
    params.precond.coarsening.aggr.block_size = options.blockSize;
    params.precond.coarsening.relax = options.relax;
    params.precond.coarsening.estimate_spectral_radius = options.estimateSpectralRadius;
    params.precond.coarsening.power_iters = options.powerIterations;
    if ( options.coarseEnough > 0 ) {
        params.precond.coarse_enough = options.coarseEnough;
    }
    params.precond.npre = options.npre;
    params.precond.npost = options.npost;
    params.precond.ncycle = options.ncycle;

    if ( options.nearNullspace && elasticMap.nearNullspaceColumns > 0 && elasticMap.nearNullspace.size() == size_t( impl->rows * elasticMap.nearNullspaceColumns ) ) {
        params.precond.coarsening.nullspace.cols = elasticMap.nearNullspaceColumns;
        params.precond.coarsening.nullspace.B = elasticMap.nearNullspace;
        std :: cout << "AMGCL-CG: using " << elasticMap.nearNullspaceColumns << " near-nullspace modes" << std :: endl;
    } else if ( options.nearNullspace ) {
        std :: cout << "AMGCL-CG: near-nullspace requested but unavailable; continuing without it" << std :: endl;
    }

    auto matrix = impl->matrixTuple();
    impl->scalarSolver.reset();
    impl->block3IluSolver.reset();
    impl->block6IluSolver.reset();
    impl->block3SpaiSolver.reset();
    impl->block6SpaiSolver.reset();
    impl->hybrid3Solver.reset();
    impl->hybrid6Solver.reset();
    impl->hybrid3FgmresSolver.reset();
    impl->hybrid6FgmresSolver.reset();
    impl->hybrid3EminSolver.reset();
    impl->hybrid6EminSolver.reset();
    impl->activeBlockSize = 1;
    impl->activeBackend = "scalar";

    auto fillBlockParams = [&](auto &blockParams) {
        blockParams.solver.tol = options.tolerance;
        blockParams.solver.maxiter = options.maxIterations;
        blockParams.solver.verbose = options.verbose;
        blockParams.precond.coarsening.aggr.eps_strong = options.epsStrong;
        blockParams.precond.coarsening.relax = options.relax;
        blockParams.precond.coarsening.estimate_spectral_radius = options.estimateSpectralRadius;
        blockParams.precond.coarsening.power_iters = options.powerIterations;
        blockParams.precond.coarsening.nullspace.cols = params.precond.coarsening.nullspace.cols;
        blockParams.precond.coarsening.nullspace.B = params.precond.coarsening.nullspace.B;
        if ( options.coarseEnough > 0 ) {
            blockParams.precond.coarse_enough = options.coarseEnough;
        }
        blockParams.precond.npre = options.npre;
        blockParams.precond.npost = options.npost;
        blockParams.precond.ncycle = options.ncycle;
    };

    auto fillEminParams = [&](auto &eminParams) {
        eminParams.solver.tol = options.tolerance;
        eminParams.solver.maxiter = options.maxIterations;
        eminParams.solver.verbose = options.verbose;
        eminParams.precond.coarsening.aggr.eps_strong = options.epsStrong;
        eminParams.precond.coarsening.aggr.block_size = options.blockSize;
        eminParams.precond.coarsening.nullspace.cols = params.precond.coarsening.nullspace.cols;
        eminParams.precond.coarsening.nullspace.B = params.precond.coarsening.nullspace.B;
        if ( options.coarseEnough > 0 ) {
            eminParams.precond.coarse_enough = options.coarseEnough;
        }
        eminParams.precond.npre = options.npre;
        eminParams.precond.npost = options.npost;
        eminParams.precond.ncycle = options.ncycle;
    };

    const std :: string backendMode = lowerCopy(options.backend);
    const std :: string coarseningMode = lowerCopy(options.coarsening);
    const std :: string krylovMode = lowerCopy(options.krylov);
    const bool useEminCoarsening = coarseningMode == "emin" || coarseningMode == "smoothed_aggr_emin" || coarseningMode == "smoothed-aggregation-emin" || coarseningMode == "smoothed_aggregation_emin";
    const bool useFgmres = krylovMode == "fgmres" || krylovMode == "flexible_gmres" || krylovMode == "flexible-gmres";
    const bool preferHybrid = ( backendMode == "auto" || backendMode == "hybrid" ) && options.useBlockBackend;
    const bool preferPlainBlock = ( backendMode == "block" || backendMode == "plain_block" || backendMode == "plain-block" ) && options.useBlockBackend;

    if ( useFgmres && useEminCoarsening ) {
        std :: cout << "AMGCL-CG: fgmres requested with smoothed_aggr_emin; falling back to smoothed_aggregation for this test path" << std :: endl;
    }

    if ( useFgmres && preferHybrid && useElasticLift && elasticMap.blockSize == 3 && impl->rows % 3 == 0 ) {
        Impl :: Hybrid3FgmresSolver :: params fgmresParams;
        fillBlockParams(fgmresParams);
        fgmresParams.solver.M = options.gmresRestart;
        impl->hybrid3FgmresSolver = std :: make_unique< Impl :: Hybrid3FgmresSolver >(matrix, fgmresParams);
        impl->activeBlockSize = 3;
        impl->activeBackend = "hybrid_fgmres";
    } else if ( useFgmres && preferHybrid && useElasticLift && elasticMap.blockSize == 6 && impl->rows % 6 == 0 ) {
        Impl :: Hybrid6FgmresSolver :: params fgmresParams;
        fillBlockParams(fgmresParams);
        fgmresParams.solver.M = options.gmresRestart;
        impl->hybrid6FgmresSolver = std :: make_unique< Impl :: Hybrid6FgmresSolver >(matrix, fgmresParams);
        impl->activeBlockSize = 6;
        impl->activeBackend = "hybrid_fgmres";
    } else if ( useEminCoarsening && preferHybrid && useElasticLift && elasticMap.blockSize == 3 && impl->rows % 3 == 0 ) {
        Impl :: Hybrid3EminSolver :: params eminParams;
        fillEminParams(eminParams);
        impl->hybrid3EminSolver = std :: make_unique< Impl :: Hybrid3EminSolver >(matrix, eminParams);
        impl->activeBlockSize = 3;
        impl->activeBackend = "hybrid_emin";
    } else if ( useEminCoarsening && preferHybrid && useElasticLift && elasticMap.blockSize == 6 && impl->rows % 6 == 0 ) {
        Impl :: Hybrid6EminSolver :: params eminParams;
        fillEminParams(eminParams);
        impl->hybrid6EminSolver = std :: make_unique< Impl :: Hybrid6EminSolver >(matrix, eminParams);
        impl->activeBlockSize = 6;
        impl->activeBackend = "hybrid_emin";
    } else if ( preferHybrid && useElasticLift && elasticMap.blockSize == 3 && impl->rows % 3 == 0 ) {
        Impl :: Hybrid3Solver :: params hybridParams;
        fillBlockParams(hybridParams);
        impl->hybrid3Solver = std :: make_unique< Impl :: Hybrid3Solver >(matrix, hybridParams);
        impl->activeBlockSize = 3;
        impl->activeBackend = "hybrid";
    } else if ( preferHybrid && useElasticLift && elasticMap.blockSize == 6 && impl->rows % 6 == 0 ) {
        Impl :: Hybrid6Solver :: params hybridParams;
        fillBlockParams(hybridParams);
        impl->hybrid6Solver = std :: make_unique< Impl :: Hybrid6Solver >(matrix, hybridParams);
        impl->activeBlockSize = 6;
        impl->activeBackend = "hybrid";
    } else if ( preferPlainBlock && useElasticLift && elasticMap.blockSize == 3 && impl->rows % 3 == 0 ) {
        auto blockMatrix = amgcl :: adapter :: block_matrix< Impl :: Block3 >(matrix);
        if ( options.blockRelaxation == "spai0" ) {
            Impl :: Block3SpaiSolver :: params blockParams;
            fillBlockParams(blockParams);
            impl->block3SpaiSolver = std :: make_unique< Impl :: Block3SpaiSolver >(blockMatrix, blockParams);
        } else {
            Impl :: Block3IluSolver :: params blockParams;
            fillBlockParams(blockParams);
            impl->block3IluSolver = std :: make_unique< Impl :: Block3IluSolver >(blockMatrix, blockParams);
        }
        impl->activeBlockSize = 3;
        impl->activeBackend = "plain_block";
    } else if ( preferPlainBlock && useElasticLift && elasticMap.blockSize == 6 && impl->rows % 6 == 0 ) {
        auto blockMatrix = amgcl :: adapter :: block_matrix< Impl :: Block6 >(matrix);
        if ( options.blockRelaxation == "spai0" ) {
            Impl :: Block6SpaiSolver :: params blockParams;
            fillBlockParams(blockParams);
            impl->block6SpaiSolver = std :: make_unique< Impl :: Block6SpaiSolver >(blockMatrix, blockParams);
        } else {
            Impl :: Block6IluSolver :: params blockParams;
            fillBlockParams(blockParams);
            impl->block6IluSolver = std :: make_unique< Impl :: Block6IluSolver >(blockMatrix, blockParams);
        }
        impl->activeBlockSize = 6;
        impl->activeBackend = "plain_block";
    } else {
        impl->scalarSolver = std :: make_unique< Impl :: ScalarSolver >(matrix, params);
        impl->activeBackend = "scalar";
    }
    impl->initialGuess.assign(impl->rows, 0.);

    std :: cout << "AMGCL-CG: setup complete, tolerance=" << options.tolerance
                << ", true_tolerance=" << ( options.trueTolerance > 0. ? options.trueTolerance : options.tolerance )
                << ", max_iterations=" << options.maxIterations
                << ", eps_strong=" << options.epsStrong
                << ", relax=" << options.relax
                << ", requested_block_size=" << options.blockSize
                << ", active_block_size=" << impl->activeBlockSize
                << ", backend=" << impl->activeBackend
                << ", coarsening=" << options.coarsening
                << ", krylov=" << options.krylov
                << ", diagonal_scale=" << ( options.diagonalScale ? 1 : 0 )
                << ", block_relaxation=" << options.blockRelaxation
                << ", lifted_rows=" << impl->rows
                << ", npre=" << options.npre
                << ", npost=" << options.npost
                << ", ncycle=" << options.ncycle << std :: endl;
    return true;
}

//////////////////////////////////////////////////////////
bool AmgclCGSolver :: solve(Vector &x, const Vector &b) {
    if ( b.size() == 0 ) {
        lastIterations = 0;
        lastError = 0.;
        lastTrueRelativeResidual = 0.;
        return true;
    }
    if ( !impl->scalarSolver && !impl->block3IluSolver && !impl->block6IluSolver && !impl->block3SpaiSolver && !impl->block6SpaiSolver && !impl->hybrid3Solver && !impl->hybrid6Solver && !impl->hybrid3FgmresSolver && !impl->hybrid6FgmresSolver && !impl->hybrid3EminSolver && !impl->hybrid6EminSolver ) {
        std :: cerr << "AMGCL-CG Error: solve called before factorize" << std :: endl;
        return false;
    }

    Vector activeB = options.elasticFullLift && elasticMap.isValid() ? liftVectorToElasticSpace(b, elasticMap) : b;
    if ( options.diagonalScale ) {
        for ( Eigen :: Index i = 0; i < activeB.size(); ++i ) {
            activeB [ i ] *= impl->scale [ i ];
        }
    }
    std :: vector< double >rhs(activeB.data(), activeB.data() + activeB.size());
    if ( impl->initialGuess.size() != size_t( activeB.size() ) ) {
        impl->initialGuess.assign(activeB.size(), 0.);
    } else if ( !options.reuseInitialGuess ) {
        std :: fill(impl->initialGuess.begin(), impl->initialGuess.end(), 0.);
    }

    auto matrix = impl->matrixTuple();
    size_t iterations = 0;
    double error = 0.;
    if ( impl->hybrid3FgmresSolver ) {
        std :: tie(iterations, error) = ( * impl->hybrid3FgmresSolver )(matrix, rhs, impl->initialGuess);
    } else if ( impl->hybrid6FgmresSolver ) {
        std :: tie(iterations, error) = ( * impl->hybrid6FgmresSolver )(matrix, rhs, impl->initialGuess);
    } else if ( impl->hybrid3EminSolver ) {
        std :: tie(iterations, error) = ( * impl->hybrid3EminSolver )(matrix, rhs, impl->initialGuess);
    } else if ( impl->hybrid6EminSolver ) {
        std :: tie(iterations, error) = ( * impl->hybrid6EminSolver )(matrix, rhs, impl->initialGuess);
    } else if ( impl->hybrid3Solver ) {
        std :: tie(iterations, error) = ( * impl->hybrid3Solver )(matrix, rhs, impl->initialGuess);
    } else if ( impl->hybrid6Solver ) {
        std :: tie(iterations, error) = ( * impl->hybrid6Solver )(matrix, rhs, impl->initialGuess);
    } else if ( impl->block3IluSolver ) {
        auto blockMatrix = amgcl :: adapter :: block_matrix< Impl :: Block3 >(matrix);
        auto blockRhs = amgcl :: backend :: reinterpret_as_rhs< Impl :: Block3 >(rhs);
        auto blockX = amgcl :: backend :: reinterpret_as_rhs< Impl :: Block3 >(impl->initialGuess);
        std :: tie(iterations, error) = ( * impl->block3IluSolver )(blockMatrix, blockRhs, blockX);
    } else if ( impl->block3SpaiSolver ) {
        auto blockMatrix = amgcl :: adapter :: block_matrix< Impl :: Block3 >(matrix);
        auto blockRhs = amgcl :: backend :: reinterpret_as_rhs< Impl :: Block3 >(rhs);
        auto blockX = amgcl :: backend :: reinterpret_as_rhs< Impl :: Block3 >(impl->initialGuess);
        std :: tie(iterations, error) = ( * impl->block3SpaiSolver )(blockMatrix, blockRhs, blockX);
    } else if ( impl->block6IluSolver ) {
        auto blockMatrix = amgcl :: adapter :: block_matrix< Impl :: Block6 >(matrix);
        auto blockRhs = amgcl :: backend :: reinterpret_as_rhs< Impl :: Block6 >(rhs);
        auto blockX = amgcl :: backend :: reinterpret_as_rhs< Impl :: Block6 >(impl->initialGuess);
        std :: tie(iterations, error) = ( * impl->block6IluSolver )(blockMatrix, blockRhs, blockX);
    } else if ( impl->block6SpaiSolver ) {
        auto blockMatrix = amgcl :: adapter :: block_matrix< Impl :: Block6 >(matrix);
        auto blockRhs = amgcl :: backend :: reinterpret_as_rhs< Impl :: Block6 >(rhs);
        auto blockX = amgcl :: backend :: reinterpret_as_rhs< Impl :: Block6 >(impl->initialGuess);
        std :: tie(iterations, error) = ( * impl->block6SpaiSolver )(blockMatrix, blockRhs, blockX);
    } else {
        std :: tie(iterations, error) = ( * impl->scalarSolver )(matrix, rhs, impl->initialGuess);
    }

    Vector activeX = Eigen :: Map< Vector >(impl->initialGuess.data(), impl->initialGuess.size());
    if ( options.diagonalScale ) {
        for ( Eigen :: Index i = 0; i < activeX.size(); ++i ) {
            activeX [ i ] *= impl->scale [ i ];
        }
    }
    x = options.elasticFullLift && elasticMap.isValid() ? restrictVectorFromElasticSpace(activeX, elasticMap) : activeX;
    lastIterations = static_cast< long long >(iterations);
    lastError = trueRelativeResidual(impl->reducedMatrix, x, b);
    lastTrueRelativeResidual = lastError;

    const double trueTolerance = options.trueTolerance > 0. ? options.trueTolerance : options.tolerance;
    const bool converged = lastTrueRelativeResidual <= trueTolerance;
    if ( !converged && options.warnOnFailure ) {
        std :: cerr << "AMGCL-CG warning: performed " << iterations
                    << " iterations and reached AMGCL reported residual " << error
                    << ", true relative residual " << lastTrueRelativeResidual
                    << ", required true tolerance is " << trueTolerance << std :: endl;
    }
    return converged;
}
#endif

#ifdef AMGCL_FOUND
//////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////
// DEFLATED FLEXIBLE GMRES SOLVER
//////////////////////////////////////////////////////////

struct DeflatedFGMRESSolver :: Impl
{
    using ScalarBackend = amgcl :: backend :: builtin< double >;
    using ScalarPrecond = amgcl :: amg<
        ScalarBackend,
        amgcl :: coarsening :: smoothed_aggregation,
        amgcl :: relaxation :: spai0
    >;
    using Block3 = amgcl :: static_matrix< double, 3, 3 >;
    using Block6 = amgcl :: static_matrix< double, 6, 6 >;
    using Hybrid3Backend = amgcl :: backend :: builtin_hybrid< Block3 >;
    using Hybrid6Backend = amgcl :: backend :: builtin_hybrid< Block6 >;
    using Hybrid3Precond = amgcl :: amg<
        Hybrid3Backend,
        amgcl :: coarsening :: smoothed_aggregation,
        amgcl :: relaxation :: spai0
    >;
    using Hybrid6Precond = amgcl :: amg<
        Hybrid6Backend,
        amgcl :: coarsening :: smoothed_aggregation,
        amgcl :: relaxation :: spai0
    >;

    ptrdiff_t rows = 0;
    std :: vector< ptrdiff_t >ptr;
    std :: vector< ptrdiff_t >col;
    std :: vector< double >val;
    CoordinateIndexedSparseMatrix reducedMatrix;
    CoordinateIndexedSparseMatrix activeMatrix;
    std :: vector< double >scale;
    std :: unique_ptr< ScalarPrecond >scalarPrecond;
    std :: unique_ptr< Hybrid3Precond >hybrid3Precond;
    std :: unique_ptr< Hybrid6Precond >hybrid6Precond;
    std :: unique_ptr< AmgclCGSolver >innerCgSolver;
#ifdef HYPRE_FOUND
    HYPRE_IJMatrix hypreIjMatrix = nullptr;
    HYPRE_ParCSRMatrix hypreParMatrix = nullptr;
    HYPRE_Solver hyprePrecond = nullptr;
    HYPRE_IJVector hypreSetupIjB = nullptr;
    HYPRE_IJVector hypreSetupIjX = nullptr;
    HYPRE_ParVector hypreSetupParB = nullptr;
    HYPRE_ParVector hypreSetupParX = nullptr;
    std :: vector< HYPRE_BigInt >hypreIndices;
#endif
    std :: vector< Vector >reducedBasis;
    std :: vector< Vector >basisU;
    std :: vector< Vector >basisC;
    unsigned discardedBasisVectors = 0;
    unsigned capacityEvictions = 0;
    unsigned nonfiniteDiscardCount = 0;
    unsigned lowNormDiscardCount = 0;
    unsigned lowANormDiscardCount = 0;
    unsigned reorthogonalizationCount = 0;
    double lastCandidateInitialANorm = 0.;
    double lastCandidateFinalANorm = 0.;
    double basisOrthogonalityMaxOffdiag = 0.;
    double basisOrthogonalityMaxDiagError = 0.;
    std :: string lastDiscardReason = "";
    std :: string activePreconditioner = "none";
    bool matrixReady = false;

    auto matrixTuple() {
        return std :: tie(rows, ptr, col, val);
    }

#ifdef HYPRE_FOUND
    void clearHyprePreconditioner() {
        if ( hyprePrecond ) {
            HYPRE_BoomerAMGDestroy(hyprePrecond);
            hyprePrecond = nullptr;
        }
        if ( hypreSetupIjB ) {
            HYPRE_IJVectorDestroy(hypreSetupIjB);
            hypreSetupIjB = nullptr;
            hypreSetupParB = nullptr;
        }
        if ( hypreSetupIjX ) {
            HYPRE_IJVectorDestroy(hypreSetupIjX);
            hypreSetupIjX = nullptr;
            hypreSetupParX = nullptr;
        }
        if ( hypreIjMatrix ) {
            HYPRE_IJMatrixDestroy(hypreIjMatrix);
            hypreIjMatrix = nullptr;
            hypreParMatrix = nullptr;
        }
        hypreIndices.clear();
    }
#endif

    void clearPreconditioner() {
        scalarPrecond.reset();
        hybrid3Precond.reset();
        hybrid6Precond.reset();
        innerCgSolver.reset();
#ifdef HYPRE_FOUND
        clearHyprePreconditioner();
#endif
        activePreconditioner = "none";
        matrixReady = false;
    }

    ~Impl() {
        clearPreconditioner();
    }
};

//////////////////////////////////////////////////////////
DeflatedFGMRESSolver :: DeflatedFGMRESSolver() {
    name = "DeflatedFGMRES";
    impl = std :: make_unique< Impl >();
    lastIterations = -1;
    lastError = -1.;
    lastTrueRelativeResidual = -1.;
    lastPreconditionerApplySeconds = 0.;
    lastOrthogonalizationSeconds = 0.;
    lastLeastSquaresSeconds = 0.;
    lastMatvecSeconds = 0.;
    lastDeflationSeconds = 0.;
}

//////////////////////////////////////////////////////////
DeflatedFGMRESSolver :: ~DeflatedFGMRESSolver() {}

//////////////////////////////////////////////////////////
void DeflatedFGMRESSolver :: setOptions(const DeflatedFGMRESOptions &opts) {
    options = opts;
}

//////////////////////////////////////////////////////////
void DeflatedFGMRESSolver :: setRuntimeTolerance(double tolerance, double trueTolerance) {
    if ( tolerance > 0. ) {
        options.tolerance = tolerance;
    }
    if ( trueTolerance > 0. ) {
        options.trueTolerance = trueTolerance;
    } else if ( tolerance > 0. ) {
        options.trueTolerance = tolerance;
    }
}

//////////////////////////////////////////////////////////
void DeflatedFGMRESSolver :: setAmgclOptions(const AmgclSolverOptions &opts) {
    amgclOptions = opts;
}

//////////////////////////////////////////////////////////
void DeflatedFGMRESSolver :: setHypreOptions(const HypreBoomerAMGOptions &opts) {
    hypreOptions = opts;
}

//////////////////////////////////////////////////////////
void DeflatedFGMRESSolver :: setElasticDofMap(ElasticDofMap map) {
    elasticMap = std :: move(map);
}

//////////////////////////////////////////////////////////
bool DeflatedFGMRESSolver :: analyzePattern(const CoordinateIndexedSparseMatrix &A) {
    ( void ) A;
    return true;
}

//////////////////////////////////////////////////////////
Vector DeflatedFGMRESSolver :: reducedVectorToActiveUnknown(const Vector &v) const {
    const bool useElasticLift = amgclOptions.elasticFullLift && elasticMap.isValid();
    Vector active = useElasticLift ? liftVectorToElasticSpace(v, elasticMap) : v;
    if ( amgclOptions.diagonalScale && impl->scale.size() == size_t( active.size() ) ) {
        for ( Eigen :: Index i = 0; i < active.size(); ++i ) {
            if ( std :: abs(impl->scale [ i ]) > 0. ) {
                active [ i ] /= impl->scale [ i ];
            }
        }
    }
    return active;
}

//////////////////////////////////////////////////////////
bool DeflatedFGMRESSolver :: appendRawDeflationVector(const Vector &v) {
    if ( options.deflationVectors == 0 || !impl->matrixReady ) {
        return false;
    }
    if ( v.size() != impl->reducedMatrix.rows() || !v.allFinite() ) {
        impl->discardedBasisVectors++;
        impl->nonfiniteDiscardCount++;
        impl->lastDiscardReason = "nonfinite_or_size";
        return false;
    }
    if ( v.norm() <= 1e-300 ) {
        impl->discardedBasisVectors++;
        impl->lowNormDiscardCount++;
        impl->lastDiscardReason = "low_vector_norm";
        return false;
    }

    if ( impl->reducedBasis.size() >= options.deflationVectors ) {
        impl->reducedBasis.erase(impl->reducedBasis.begin() );
        if ( !impl->basisU.empty() ) {
            impl->basisU.erase(impl->basisU.begin() );
        }
        if ( !impl->basisC.empty() ) {
            impl->basisC.erase(impl->basisC.begin() );
        }
        impl->capacityEvictions++;
        impl->lastDiscardReason = "capacity_eviction";
    }
    const bool accepted = appendActiveDeflationVector(v, true);
    updateDeflationOrthogonalityDiagnostics();
    return accepted;
}

//////////////////////////////////////////////////////////
bool DeflatedFGMRESSolver :: appendActiveDeflationVector(const Vector &rawReducedVector, bool storeRawVector) {
    Vector u = reducedVectorToActiveUnknown(rawReducedVector);
    if ( !u.allFinite() || u.norm() <= 1e-300 ) {
        impl->discardedBasisVectors++;
        impl->lowNormDiscardCount++;
        impl->lastDiscardReason = "low_vector_norm";
        return false;
    }

    Vector au = impl->activeMatrix * u;
    const double initialANorm = u.dot(au);
    impl->lastCandidateInitialANorm = initialANorm;
    if ( !std :: isfinite(initialANorm) ) {
        impl->discardedBasisVectors++;
        impl->nonfiniteDiscardCount++;
        impl->lastDiscardReason = "nonfinite_a_norm";
        return false;
    }

    // Modified Gram-Schmidt in the A inner product. Since basisC stores A*basisU
    // and basisU^T*basisC is kept close to identity, no coarse Gram inverse is needed.
    for ( int pass = 0; pass < 2; ++pass ) {
        for ( size_t basisIndex = 0; basisIndex < impl->basisU.size(); ++basisIndex ) {
            const double coefficient = impl->basisU [ basisIndex ].dot(au);
            u.noalias() -= coefficient * impl->basisU [ basisIndex ];
            au.noalias() -= coefficient * impl->basisC [ basisIndex ];
        }
    }

    const double finalANorm = u.dot(au);
    impl->lastCandidateFinalANorm = finalANorm;
    if ( !std :: isfinite(finalANorm) || std :: abs(finalANorm) <= options.deflationEps ) {
        impl->discardedBasisVectors++;
        impl->lowANormDiscardCount++;
        impl->lastDiscardReason = "low_a_norm";
        return false;
    }

    const double invANorm = 1. / sqrt(std :: abs(finalANorm) );
    u *= invANorm;
    au *= invANorm;
    if ( !u.allFinite() || !au.allFinite() ) {
        impl->discardedBasisVectors++;
        impl->nonfiniteDiscardCount++;
        impl->lastDiscardReason = "nonfinite_orthogonalized";
        return false;
    }

    impl->basisU.push_back(u);
    impl->basisC.push_back(au);
    if ( storeRawVector ) {
        impl->reducedBasis.push_back(rawReducedVector);
    }
    impl->lastDiscardReason = "accepted";
    return true;
}

//////////////////////////////////////////////////////////
void DeflatedFGMRESSolver :: updateDeflationOrthogonalityDiagnostics() {
    impl->basisOrthogonalityMaxOffdiag = 0.;
    impl->basisOrthogonalityMaxDiagError = 0.;

    for ( size_t i = 0; i < impl->basisU.size(); ++i ) {
        for ( size_t j = 0; j < impl->basisU.size(); ++j ) {
            const double value = impl->basisU [ i ].dot(impl->basisC [ j ]);
            if ( i == j ) {
                impl->basisOrthogonalityMaxDiagError = std :: max(impl->basisOrthogonalityMaxDiagError, std :: abs(value - 1.) );
            } else {
                impl->basisOrthogonalityMaxOffdiag = std :: max(impl->basisOrthogonalityMaxOffdiag, std :: abs(value) );
            }
        }
    }
}

//////////////////////////////////////////////////////////
void DeflatedFGMRESSolver :: rebuildDeflationBasis() {
    std :: vector< Vector >rawBasis = impl->reducedBasis;
    impl->basisU.clear();
    impl->basisC.clear();
    impl->reducedBasis.clear();
    impl->basisOrthogonalityMaxOffdiag = 0.;
    impl->basisOrthogonalityMaxDiagError = 0.;
    impl->lastDiscardReason = "none";

    if ( !impl->matrixReady || rawBasis.empty() || options.deflationVectors == 0 ) {
        return;
    }

    for ( const Vector &rawVector : rawBasis ) {
        if ( impl->reducedBasis.size() >= options.deflationVectors ) {
            impl->reducedBasis.erase(impl->reducedBasis.begin() );
            if ( !impl->basisU.empty() ) {
                impl->basisU.erase(impl->basisU.begin() );
            }
            if ( !impl->basisC.empty() ) {
                impl->basisC.erase(impl->basisC.begin() );
            }
            impl->capacityEvictions++;
        }
        appendActiveDeflationVector(rawVector, true);
    }
    updateDeflationOrthogonalityDiagnostics();
    impl->reorthogonalizationCount++;
}

//////////////////////////////////////////////////////////
Vector DeflatedFGMRESSolver :: projectDeflation(const Vector &v) const {
    if ( impl->basisU.empty() ) {
        return v;
    }
    Vector projected = v;
    for ( size_t i = 0; i < impl->basisU.size(); ++i ) {
        const double coefficient = impl->basisC [ i ].dot(projected);
        projected.noalias() -= coefficient * impl->basisU [ i ];
    }
    return projected;
}

//////////////////////////////////////////////////////////
Vector DeflatedFGMRESSolver :: applyPreconditioner(const Vector &v) {
    const auto start = std :: chrono :: steady_clock :: now();
    auto finish = [&](const Vector &result) {
        lastPreconditionerApplySeconds += std :: chrono :: duration< double >(std :: chrono :: steady_clock :: now() - start).count();
        return result;
    };

    if ( lowerCopy(options.preconditioner) == "none" || lowerCopy(options.preconditioner) == "identity" ) {
        return finish(v);
    }
    if ( impl->innerCgSolver ) {
        Vector result = Vector :: Zero(v.size() );
        const bool ok = impl->innerCgSolver->solve(result, v);
        if ( !ok && options.verbose ) {
            std :: cerr << "DeflatedFGMRES inner AMGCL-CG preconditioner did not reach its requested tolerance; continuing with available correction" << std :: endl;
        }
        return finish(result);
    }
#ifdef HYPRE_FOUND
    if ( impl->hyprePrecond && impl->hypreParMatrix ) {
        if ( v.size() != Eigen :: Index(impl->hypreIndices.size() ) ) {
            if ( options.verbose ) {
                std :: cerr << "DeflatedFGMRES hypre preconditioner size mismatch; using identity fallback" << std :: endl;
            }
            return finish(v);
        }

        std :: vector< double >rhs(v.data(), v.data() + v.size() );
        std :: vector< double >guess(v.size(), 0.);
        std :: vector< double >answer(v.size(), 0.);
        HYPRE_IJVector ijB = nullptr;
        HYPRE_IJVector ijX = nullptr;
        HYPRE_ParVector parB = nullptr;
        HYPRE_ParVector parX = nullptr;

        auto cleanup = [&]() {
            if ( ijB ) {
                HYPRE_IJVectorDestroy(ijB);
            }
            if ( ijX ) {
                HYPRE_IJVectorDestroy(ijX);
            }
        };

        if ( !makeHypreVector(impl->hypreIndices, rhs, ijB, parB) || !makeHypreVector(impl->hypreIndices, guess, ijX, parX) ) {
            cleanup();
            if ( options.verbose ) {
                std :: cerr << "DeflatedFGMRES hypre vector setup failed; using identity fallback" << std :: endl;
            }
            return finish(v);
        }

        const HYPRE_Int solveError = HYPRE_BoomerAMGSolve(impl->hyprePrecond, impl->hypreParMatrix, parB, parX);
        if ( solveError != 0 ) {
            HYPRE_ClearAllErrors();
            cleanup();
            if ( options.verbose ) {
                std :: cerr << "DeflatedFGMRES hypre preconditioner apply failed with error " << solveError << "; using identity fallback" << std :: endl;
            }
            return finish(v);
        }
        if ( !hypreCheck(HYPRE_IJVectorGetValues(ijX, static_cast< HYPRE_Int >( impl->hypreIndices.size() ), impl->hypreIndices.data(), answer.data() ), "HYPRE_IJVectorGetValues") ) {
            cleanup();
            return finish(v);
        }
        cleanup();
        return finish(Eigen :: Map< Vector >(answer.data(), answer.size() ) );
    }
#endif
    if ( !impl->scalarPrecond && !impl->hybrid3Precond && !impl->hybrid6Precond ) {
        return finish(v);
    }

    std :: vector< double >rhs(v.data(), v.data() + v.size() );
    std :: vector< double >result(v.size(), 0.);
    if ( impl->hybrid3Precond ) {
        impl->hybrid3Precond->apply(rhs, result);
    } else if ( impl->hybrid6Precond ) {
        impl->hybrid6Precond->apply(rhs, result);
    } else {
        impl->scalarPrecond->apply(rhs, result);
    }
    return finish(Eigen :: Map< Vector >(result.data(), result.size() ) );
}

//////////////////////////////////////////////////////////
bool DeflatedFGMRESSolver :: factorize(const CoordinateIndexedSparseMatrix &A) {
    if ( A.rows() <= 0 ) {
        return true;
    }

    impl->clearPreconditioner();
    const bool useElasticLift = amgclOptions.elasticFullLift && elasticMap.isValid();
    impl->reducedMatrix = A;
    impl->activeMatrix = useElasticLift ? liftMatrixToElasticSpace(A, elasticMap) : A;
    impl->scale.assign(impl->activeMatrix.rows(), 1.);

    if ( amgclOptions.diagonalScale ) {
        double minScale = std :: numeric_limits< double > :: infinity();
        double maxScale = 0.;
        long long nonPositiveDiagonal = 0;
        for ( int i = 0; i < impl->activeMatrix.rows(); ++i ) {
            const double diag = impl->activeMatrix.coeff(i, i);
            if ( diag > 0. && std :: isfinite(diag) ) {
                impl->scale [ i ] = 1. / sqrt(diag);
                minScale = std :: min(minScale, impl->scale [ i ]);
                maxScale = std :: max(maxScale, impl->scale [ i ]);
            } else {
                nonPositiveDiagonal++;
            }
        }
        impl->activeMatrix = symmetricallyScaleMatrix(impl->activeMatrix, impl->scale);
        std :: cout << "DeflatedFGMRES: symmetric diagonal scaling enabled, min_scale=" << minScale
                    << ", max_scale=" << maxScale
                    << ", non_positive_diagonal_count=" << nonPositiveDiagonal << std :: endl;
    }

    makeRowMajorCrs(impl->activeMatrix, impl->rows, impl->ptr, impl->col, impl->val);
    auto matrix = impl->matrixTuple();

    auto fillPrecondParams = [&](auto &params) {
        params.coarsening.aggr.eps_strong = amgclOptions.epsStrong;
        params.coarsening.aggr.block_size = amgclOptions.blockSize;
        params.coarsening.relax = amgclOptions.relax;
        params.coarsening.estimate_spectral_radius = amgclOptions.estimateSpectralRadius;
        params.coarsening.power_iters = amgclOptions.powerIterations;
        if ( amgclOptions.coarseEnough > 0 ) {
            params.coarse_enough = amgclOptions.coarseEnough;
        }
        params.npre = amgclOptions.npre;
        params.npost = amgclOptions.npost;
        params.ncycle = amgclOptions.ncycle;
        if ( amgclOptions.nearNullspace && elasticMap.nearNullspaceColumns > 0 && elasticMap.nearNullspace.size() == size_t( impl->rows * elasticMap.nearNullspaceColumns ) ) {
            params.coarsening.nullspace.cols = elasticMap.nearNullspaceColumns;
            params.coarsening.nullspace.B = elasticMap.nearNullspace;
        }
    };

    const std :: string preconditionerMode = lowerCopy(options.preconditioner);
    const std :: string backendMode = lowerCopy(amgclOptions.backend);
    const bool preferHybrid = preconditionerMode != "none" && preconditionerMode != "identity"
        && ( backendMode == "auto" || backendMode == "hybrid" )
        && amgclOptions.useBlockBackend;

    const bool useInnerCg = preconditionerMode == "amgcl_cg" || preconditionerMode == "amgclcg"
        || preconditionerMode == "amgcl_solver" || preconditionerMode == "amgcl-solver";
    const bool useHypreApply = preconditionerMode == "hypre" || preconditionerMode == "boomeramg"
        || preconditionerMode == "hypre_boomeramg" || preconditionerMode == "hypre-boomeramg";

    try {
        if ( preconditionerMode == "none" || preconditionerMode == "identity" ) {
            impl->activePreconditioner = "identity";
        } else if ( useHypreApply ) {
#ifdef HYPRE_FOUND
            ensureHypreRuntimeInitialized();
            if ( impl->rows <= 0 || impl->rows > std :: numeric_limits< HYPRE_Int > :: max() ) {
                std :: cerr << "DeflatedFGMRES hypre preconditioner supports positive row counts fitting HYPRE_Int; rows=" << impl->rows << std :: endl;
                impl->clearPreconditioner();
                return false;
            }

            const HYPRE_BigInt firstRow = 0;
            const HYPRE_BigInt lastRow = impl->rows - 1;
            impl->hypreIndices.resize(impl->rows);
            for ( ptrdiff_t row = 0; row < impl->rows; row++ ) {
                impl->hypreIndices [ row ] = row;
            }

            if ( !hypreCheck(HYPRE_IJMatrixCreate(MPI_COMM_SELF, firstRow, lastRow, firstRow, lastRow, &impl->hypreIjMatrix), "HYPRE_IJMatrixCreate") ) {
                impl->clearPreconditioner();
                return false;
            }
            if ( !hypreCheck(HYPRE_IJMatrixSetObjectType(impl->hypreIjMatrix, HYPRE_PARCSR), "HYPRE_IJMatrixSetObjectType") ) {
                impl->clearPreconditioner();
                return false;
            }
            std :: vector< HYPRE_Int >rowSizes(impl->rows, 0);
            for ( ptrdiff_t row = 0; row < impl->rows; row++ ) {
                rowSizes [ row ] = static_cast< HYPRE_Int >( impl->ptr [ row + 1 ] - impl->ptr [ row ] );
            }
            HYPRE_IJMatrixSetRowSizes(impl->hypreIjMatrix, rowSizes.data() );
            HYPRE_IJMatrixSetOMPFlag(impl->hypreIjMatrix, 1);
            if ( !hypreCheck(HYPRE_IJMatrixInitialize(impl->hypreIjMatrix), "HYPRE_IJMatrixInitialize") ) {
                impl->clearPreconditioner();
                return false;
            }

            std :: vector< HYPRE_BigInt >rowColumns;
            for ( ptrdiff_t row = 0; row < impl->rows; row++ ) {
                HYPRE_Int rowNnz = static_cast< HYPRE_Int >( impl->ptr [ row + 1 ] - impl->ptr [ row ] );
                if ( rowNnz == 0 ) {
                    continue;
                }
                rowColumns.resize(rowNnz);
                for ( HYPRE_Int j = 0; j < rowNnz; j++ ) {
                    rowColumns [ j ] = impl->col [ impl->ptr [ row ] + j ];
                }
                const HYPRE_BigInt rowId = row;
                if ( !hypreCheck(HYPRE_IJMatrixSetValues(impl->hypreIjMatrix, 1, &rowNnz, &rowId, rowColumns.data(), impl->val.data() + impl->ptr [ row ] ), "HYPRE_IJMatrixSetValues") ) {
                    impl->clearPreconditioner();
                    return false;
                }
            }
            if ( !hypreCheck(HYPRE_IJMatrixAssemble(impl->hypreIjMatrix), "HYPRE_IJMatrixAssemble") ) {
                impl->clearPreconditioner();
                return false;
            }
            if ( !hypreCheck(HYPRE_IJMatrixGetObject(impl->hypreIjMatrix, reinterpret_cast< void ** >( &impl->hypreParMatrix ) ), "HYPRE_IJMatrixGetObject") ) {
                impl->clearPreconditioner();
                return false;
            }
            if ( !hypreCheck(HYPRE_BoomerAMGCreate(&impl->hyprePrecond), "HYPRE_BoomerAMGCreate") ) {
                impl->clearPreconditioner();
                return false;
            }
            HYPRE_BoomerAMGSetTol(impl->hyprePrecond, 0.);
            HYPRE_BoomerAMGSetMaxIter(impl->hyprePrecond, hypreOptions.boomerMaxIterations);
            HYPRE_BoomerAMGSetCoarsenType(impl->hyprePrecond, hypreOptions.coarsenType);
            HYPRE_BoomerAMGSetInterpType(impl->hyprePrecond, hypreOptions.interpType);
            HYPRE_BoomerAMGSetStrongThreshold(impl->hyprePrecond, hypreOptions.strongThreshold);
            HYPRE_BoomerAMGSetRelaxType(impl->hyprePrecond, hypreOptions.relaxType);
            HYPRE_BoomerAMGSetRelaxOrder(impl->hyprePrecond, hypreOptions.relaxOrder);
            HYPRE_BoomerAMGSetPMaxElmts(impl->hyprePrecond, hypreOptions.pMaxElmts);
            HYPRE_BoomerAMGSetAggNumLevels(impl->hyprePrecond, hypreOptions.aggNumLevels);
            HYPRE_BoomerAMGSetPrintLevel(impl->hyprePrecond, hypreOptions.printLevel);
            if ( hypreOptions.nonGalerkinTol >= 0. ) {
                HYPRE_Real tolerances [ 3 ] = { 0., static_cast< HYPRE_Real >( hypreOptions.nonGalerkinTol ), static_cast< HYPRE_Real >( hypreOptions.nonGalerkinTol ) };
                HYPRE_BoomerAMGSetNonGalerkTol(impl->hyprePrecond, 3, tolerances);
            }
            if ( useElasticLift && elasticMap.blockSize > 1 ) {
                const HYPRE_Int numFunctions = static_cast< HYPRE_Int >( hypreOptions.numFunctions > 0 ? hypreOptions.numFunctions : elasticMap.blockSize );
                HYPRE_BoomerAMGSetNumFunctions(impl->hyprePrecond, numFunctions);
                if ( hypreOptions.nodal > 0 ) {
                    HYPRE_BoomerAMGSetNodal(impl->hyprePrecond, hypreOptions.nodal);
                }
                if ( hypreOptions.useDofFunctions && options.verbose ) {
                    std :: cout << "DeflatedFGMRES hypre: explicit dof_func is disabled; using hypre's default cyclic node-major mapping." << std :: endl;
                }
                if ( hypreOptions.useInterpVectors && elasticMap.nearNullspaceColumns > 0 && elasticMap.nearNullspace.size() == size_t( impl->rows * elasticMap.nearNullspaceColumns ) && options.verbose ) {
                    std :: cout << "DeflatedFGMRES hypre: interpolation vectors requested but skipped pending an ownership-safe wrapper." << std :: endl;
                }
            }

            std :: vector< double >setupRhs(impl->rows, 1.);
            std :: vector< double >setupGuess(impl->rows, 0.);
            if ( !makeHypreVector(impl->hypreIndices, setupRhs, impl->hypreSetupIjB, impl->hypreSetupParB) ) {
                impl->clearPreconditioner();
                return false;
            }
            if ( !makeHypreVector(impl->hypreIndices, setupGuess, impl->hypreSetupIjX, impl->hypreSetupParX) ) {
                impl->clearPreconditioner();
                return false;
            }
            if ( !hypreCheck(HYPRE_BoomerAMGSetup(impl->hyprePrecond, impl->hypreParMatrix, impl->hypreSetupParB, impl->hypreSetupParX), "HYPRE_BoomerAMGSetup") ) {
                impl->clearPreconditioner();
                return false;
            }
            impl->activePreconditioner = "hypre_boomeramg_apply";
#else
            std :: cerr << "DeflatedFGMRES hypre preconditioner requested, but this OAS build was compiled without hypre support" << std :: endl;
            impl->clearPreconditioner();
            return false;
#endif
        } else if ( useInnerCg ) {
            AmgclSolverOptions innerOptions = amgclOptions;
            ElasticDofMap innerMap;
            if ( useElasticLift && elasticMap.isValid() ) {
                innerMap.fullRows = impl->activeMatrix.rows();
                innerMap.blockSize = elasticMap.blockSize;
                innerMap.dimension = elasticMap.dimension;
                innerMap.reducedToFull.resize(innerMap.fullRows);
                std :: iota(innerMap.reducedToFull.begin(), innerMap.reducedToFull.end(), 0u);
                innerMap.coordinates = elasticMap.coordinates;
                innerMap.nearNullspace = elasticMap.nearNullspace;
                innerMap.nearNullspaceColumns = elasticMap.nearNullspaceColumns;
                innerOptions.elasticFullLift = true;
            } else {
                innerOptions.elasticFullLift = false;
            }
            innerOptions.diagonalScale = false;
            innerOptions.reuseInitialGuess = false;
            innerOptions.warnOnFailure = false;
            innerOptions.krylov = "cg";
            impl->innerCgSolver = std :: make_unique< AmgclCGSolver >();
            impl->innerCgSolver->setOptions(innerOptions);
            impl->innerCgSolver->setElasticDofMap(std :: move(innerMap) );
            if ( !impl->innerCgSolver->factorize(impl->activeMatrix) ) {
                impl->clearPreconditioner();
                return false;
            }
            impl->activePreconditioner = "amgcl_cg_inner";
        } else if ( preferHybrid && useElasticLift && elasticMap.blockSize == 3 && impl->rows % 3 == 0 ) {
            Impl :: Hybrid3Precond :: params params;
            fillPrecondParams(params);
            impl->hybrid3Precond = std :: make_unique< Impl :: Hybrid3Precond >(matrix, params);
            impl->activePreconditioner = "amgcl_hybrid3_spai0";
        } else if ( preferHybrid && useElasticLift && elasticMap.blockSize == 6 && impl->rows % 6 == 0 ) {
            Impl :: Hybrid6Precond :: params params;
            fillPrecondParams(params);
            impl->hybrid6Precond = std :: make_unique< Impl :: Hybrid6Precond >(matrix, params);
            impl->activePreconditioner = "amgcl_hybrid6_spai0";
        } else {
            Impl :: ScalarPrecond :: params params;
            fillPrecondParams(params);
            impl->scalarPrecond = std :: make_unique< Impl :: ScalarPrecond >(matrix, params);
            impl->activePreconditioner = "amgcl_scalar_spai0";
        }
    } catch ( const std :: exception &e ) {
        std :: cerr << "DeflatedFGMRES preconditioner setup failed: " << e.what() << std :: endl;
        impl->clearPreconditioner();
        return false;
    }

    impl->matrixReady = true;
    if ( options.reorthogonalizeOnMatrixChange ) {
        rebuildDeflationBasis();
    } else if ( !impl->reducedBasis.empty() ) {
        rebuildDeflationBasis();
    }

    std :: cout << "DeflatedFGMRES: setup complete, tolerance=" << options.tolerance
                << ", true_tolerance=" << options.trueTolerance
                << ", max_iterations=" << options.maxIterations
                << ", restart=" << options.restart
                << ", deflation_vectors=" << options.deflationVectors
                << ", active_basis_size=" << impl->basisU.size()
                << ", discarded_basis_vectors=" << impl->discardedBasisVectors
                << ", deflation_eps=" << options.deflationEps
                << ", preconditioner=" << impl->activePreconditioner
                << ", lifted_rows=" << impl->rows
                << ", npre=" << amgclOptions.npre
                << ", npost=" << amgclOptions.npost
                << ", ncycle=" << amgclOptions.ncycle << std :: endl;
    return true;
}

//////////////////////////////////////////////////////////
bool DeflatedFGMRESSolver :: solve(Vector &x, const Vector &b) {
    lastPreconditionerApplySeconds = 0.;
    lastOrthogonalizationSeconds = 0.;
    lastLeastSquaresSeconds = 0.;
    lastMatvecSeconds = 0.;
    lastDeflationSeconds = 0.;

    if ( b.size() == 0 ) {
        lastIterations = 0;
        lastError = 0.;
        lastTrueRelativeResidual = 0.;
        return true;
    }
    if ( !impl->matrixReady ) {
        std :: cerr << "DeflatedFGMRES Error: solve called before factorize" << std :: endl;
        return false;
    }

    const bool useElasticLift = amgclOptions.elasticFullLift && elasticMap.isValid();
    Vector activeB = useElasticLift ? liftVectorToElasticSpace(b, elasticMap) : b;
    if ( amgclOptions.diagonalScale ) {
        for ( Eigen :: Index i = 0; i < activeB.size(); ++i ) {
            activeB [ i ] *= impl->scale [ i ];
        }
    }

    const Eigen :: Index n = activeB.size();
    const unsigned restart = std :: max(1u, std :: min(options.restart, options.maxIterations) );
    Vector activeX = Vector :: Zero(n);
    const double bNorm = std :: max(activeB.norm(), 1e-300);

    if ( !impl->basisU.empty() ) {
        const auto deflationStart = std :: chrono :: steady_clock :: now();
        for ( size_t i = 0; i < impl->basisU.size(); ++i ) {
            const double coefficient = impl->basisU [ i ].dot(activeB);
            activeX.noalias() += coefficient * impl->basisU [ i ];
        }
        lastDeflationSeconds += std :: chrono :: duration< double >(std :: chrono :: steady_clock :: now() - deflationStart).count();
    }

    auto matvecStart = std :: chrono :: steady_clock :: now();
    Vector residual = activeB - impl->activeMatrix * activeX;
    lastMatvecSeconds += std :: chrono :: duration< double >(std :: chrono :: steady_clock :: now() - matvecStart).count();
    double trueResidual = residual.norm() / bNorm;
    if ( trueResidual <= options.trueTolerance ) {
        lastIterations = 0;
        lastError = trueResidual;
        lastTrueRelativeResidual = trueResidual;
    } else {
        unsigned totalIterations = 0;
        double gmresResidual = trueResidual;
        while ( totalIterations < options.maxIterations && trueResidual > options.trueTolerance ) {
            const double beta = residual.norm();
            if ( beta / bNorm <= options.tolerance ) {
                break;
            }

            const unsigned innerLimit = std :: min(restart, options.maxIterations - totalIterations);
            Eigen :: MatrixXd V = Eigen :: MatrixXd :: Zero(n, innerLimit + 1);
            Eigen :: MatrixXd Z = Eigen :: MatrixXd :: Zero(n, innerLimit);
            Eigen :: MatrixXd H = Eigen :: MatrixXd :: Zero(innerLimit + 1, innerLimit);
            Eigen :: VectorXd g = Eigen :: VectorXd :: Zero(innerLimit + 1);
            V.col(0) = residual / beta;
            g [ 0 ] = beta;

            Eigen :: VectorXd y;
            unsigned usedInner = 0;
            bool happyBreakdown = false;
            for ( unsigned j = 0; j < innerLimit; ++j ) {
                Vector basisVector = V.col(j);
                Vector z = applyPreconditioner(basisVector);
                const auto deflationStart = std :: chrono :: steady_clock :: now();
                z = projectDeflation(z);
                lastDeflationSeconds += std :: chrono :: duration< double >(std :: chrono :: steady_clock :: now() - deflationStart).count();
                matvecStart = std :: chrono :: steady_clock :: now();
                Vector w = impl->activeMatrix * z;
                lastMatvecSeconds += std :: chrono :: duration< double >(std :: chrono :: steady_clock :: now() - matvecStart).count();

                const auto orthogonalizationStart = std :: chrono :: steady_clock :: now();
                for ( unsigned i = 0; i <= j; ++i ) {
                    H(i, j) = V.col(i).dot(w);
                    w.noalias() -= H(i, j) * V.col(i);
                }
                if ( options.reorthogonalizeKrylov ) {
                    for ( unsigned i = 0; i <= j; ++i ) {
                        const double correction = V.col(i).dot(w);
                        H(i, j) += correction;
                        w.noalias() -= correction * V.col(i);
                    }
                }
                lastOrthogonalizationSeconds += std :: chrono :: duration< double >(std :: chrono :: steady_clock :: now() - orthogonalizationStart).count();

                H(j + 1, j) = w.norm();
                Z.col(j) = z;
                usedInner = j + 1;
                if ( H(j + 1, j) > 1e-14 ) {
                    V.col(j + 1) = w / H(j + 1, j);
                }

                Eigen :: MatrixXd subH = H.block(0, 0, j + 2, j + 1);
                Eigen :: VectorXd subg = g.head(j + 2);
                const auto leastSquaresStart = std :: chrono :: steady_clock :: now();
                y = subH.colPivHouseholderQr().solve(subg);
                lastLeastSquaresSeconds += std :: chrono :: duration< double >(std :: chrono :: steady_clock :: now() - leastSquaresStart).count();
                gmresResidual = ( subg - subH * y ).norm() / bNorm;
                totalIterations++;

                if ( options.verbose ) {
                    std :: cout << "DeflatedFGMRES iter " << totalIterations
                                << " gmres_residual=" << gmresResidual
                                << " basis_size=" << impl->basisU.size() << std :: endl;
                }
                if ( H(j + 1, j) <= 1e-14 ) {
                    happyBreakdown = true;
                    break;
                }
                if ( gmresResidual <= options.tolerance || totalIterations >= options.maxIterations ) {
                    break;
                }
            }

            if ( usedInner > 0 ) {
                if ( y.size() != Eigen :: Index(usedInner) ) {
                    Eigen :: MatrixXd subH = H.block(0, 0, usedInner + 1, usedInner);
                    Eigen :: VectorXd subg = g.head(usedInner + 1);
                    const auto leastSquaresStart = std :: chrono :: steady_clock :: now();
                    y = subH.colPivHouseholderQr().solve(subg);
                    lastLeastSquaresSeconds += std :: chrono :: duration< double >(std :: chrono :: steady_clock :: now() - leastSquaresStart).count();
                }
                activeX.noalias() += Z.leftCols(usedInner) * y;
            }

            matvecStart = std :: chrono :: steady_clock :: now();
            residual = activeB - impl->activeMatrix * activeX;
            lastMatvecSeconds += std :: chrono :: duration< double >(std :: chrono :: steady_clock :: now() - matvecStart).count();
            trueResidual = residual.norm() / bNorm;
            if ( trueResidual <= options.trueTolerance || happyBreakdown || usedInner == 0 || totalIterations >= options.maxIterations ) {
                break;
            }
        }
        lastIterations = totalIterations;
        lastError = gmresResidual;
        lastTrueRelativeResidual = trueResidual;
    }

    if ( amgclOptions.diagonalScale ) {
        for ( Eigen :: Index i = 0; i < activeX.size(); ++i ) {
            activeX [ i ] *= impl->scale [ i ];
        }
    }
    x = useElasticLift ? restrictVectorFromElasticSpace(activeX, elasticMap) : activeX;

    lastTrueRelativeResidual = trueRelativeResidual(impl->reducedMatrix, x, b);
    lastError = lastTrueRelativeResidual;
    const bool converged = lastTrueRelativeResidual <= options.trueTolerance;
    if ( !converged ) {
        std :: cerr << "DeflatedFGMRES warning: performed " << lastIterations
                    << " iterations and reached true relative residual " << lastTrueRelativeResidual
                    << ", required true tolerance is " << options.trueTolerance
                    << ", active_basis_size=" << impl->basisU.size() << std :: endl;
    }
    return converged;
}

//////////////////////////////////////////////////////////
void DeflatedFGMRESSolver :: collectDeflationVector(const Vector &v) {
    if ( !options.collectNewtonSteps || options.deflationVectors == 0 || !impl->matrixReady ) {
        return;
    }
    appendRawDeflationVector(v);
    std :: cout << "DeflatedFGMRES: collected Newton increment, raw_candidate_count=" << impl->reducedBasis.size()
                << ", active_basis_size=" << impl->basisU.size()
                << ", discarded_basis_vectors=" << impl->discardedBasisVectors
                << ", low_a_norm_discards=" << impl->lowANormDiscardCount
                << ", capacity_evictions=" << impl->capacityEvictions
                << ", last_initial_a_norm=" << impl->lastCandidateInitialANorm
                << ", last_final_a_norm=" << impl->lastCandidateFinalANorm
                << ", basis_orthogonality_max_offdiag=" << impl->basisOrthogonalityMaxOffdiag
                << ", basis_orthogonality_max_diag_error=" << impl->basisOrthogonalityMaxDiagError
                << ", last_discard_reason=" << ( impl->lastDiscardReason.empty() ? "-" : impl->lastDiscardReason )
                << std :: endl;
}

//////////////////////////////////////////////////////////
unsigned DeflatedFGMRESSolver :: giveDeflationBasisSize() const {
    return static_cast< unsigned >( impl->basisU.size() );
}

//////////////////////////////////////////////////////////
unsigned DeflatedFGMRESSolver :: giveDeflationDiscardedCount() const {
    return impl->discardedBasisVectors;
}

//////////////////////////////////////////////////////////
unsigned DeflatedFGMRESSolver :: giveDeflationRawCandidateCount() const {
    return static_cast< unsigned >( impl->reducedBasis.size() );
}

//////////////////////////////////////////////////////////
unsigned DeflatedFGMRESSolver :: giveDeflationCapacityEvictionCount() const {
    return impl->capacityEvictions;
}

//////////////////////////////////////////////////////////
unsigned DeflatedFGMRESSolver :: giveDeflationLowANormDiscardCount() const {
    return impl->lowANormDiscardCount;
}

//////////////////////////////////////////////////////////
double DeflatedFGMRESSolver :: giveDeflationLastInitialANorm() const {
    return impl->lastCandidateInitialANorm;
}

//////////////////////////////////////////////////////////
double DeflatedFGMRESSolver :: giveDeflationLastFinalANorm() const {
    return impl->lastCandidateFinalANorm;
}

//////////////////////////////////////////////////////////
double DeflatedFGMRESSolver :: giveDeflationOrthogonalityMaxOffdiag() const {
    return impl->basisOrthogonalityMaxOffdiag;
}

//////////////////////////////////////////////////////////
double DeflatedFGMRESSolver :: giveDeflationOrthogonalityMaxDiagError() const {
    return impl->basisOrthogonalityMaxDiagError;
}

//////////////////////////////////////////////////////////
std :: string DeflatedFGMRESSolver :: giveDeflationLastDiscardReason() const {
    return impl->lastDiscardReason;
}
#endif

#ifdef HYPRE_FOUND
//////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////
// HYPRE BOOMERAMG-PRECONDITIONED CONJUGATE GRADIENT SOLVER
//////////////////////////////////////////////////////////

struct HypreBoomerAMGCGSolver :: Impl
{
    HYPRE_IJMatrix ijMatrix = nullptr;
    HYPRE_ParCSRMatrix parMatrix = nullptr;
    HYPRE_Solver solver = nullptr;
    HYPRE_Solver precond = nullptr;
    HYPRE_IJVector setupIjB = nullptr;
    HYPRE_IJVector setupIjX = nullptr;
    HYPRE_ParVector setupParB = nullptr;
    HYPRE_ParVector setupParX = nullptr;
    CoordinateIndexedSparseMatrix reducedMatrix;
    CoordinateIndexedSparseMatrix activeMatrix;
    HYPRE_BigInt firstRow = 0;
    HYPRE_BigInt lastRow = -1;
    std :: vector< HYPRE_Int >dofFunc;
    std :: vector< HYPRE_BigInt >indices;
    bool matrixReady = false;

    void clearSolverSetup() {
        if ( solver ) {
            HYPRE_ParCSRPCGDestroy(solver);
            solver = nullptr;
        }
        if ( precond ) {
            HYPRE_BoomerAMGDestroy(precond);
            precond = nullptr;
        }
        if ( setupIjB ) {
            HYPRE_IJVectorDestroy(setupIjB);
            setupIjB = nullptr;
            setupParB = nullptr;
        }
        if ( setupIjX ) {
            HYPRE_IJVectorDestroy(setupIjX);
            setupIjX = nullptr;
            setupParX = nullptr;
        }
    }

    void clearMatrix() {
        if ( ijMatrix ) {
            HYPRE_IJMatrixDestroy(ijMatrix);
            ijMatrix = nullptr;
            parMatrix = nullptr;
        }
        matrixReady = false;
    }

    ~Impl() {
        clearSolverSetup();
        clearMatrix();
    }
};

namespace {

void ensureHypreRuntimeInitialized() {
    static bool initializedHere = false;
    int mpiInitialized = 0;
    MPI_Initialized(&mpiInitialized);
    if ( !mpiInitialized ) {
        int provided = 0;
        MPI_Init_thread(nullptr, nullptr, MPI_THREAD_FUNNELED, &provided);
    }
    if ( !HYPRE_Initialized() ) {
        HYPRE_Initialize();
        HYPRE_SetMemoryLocation(HYPRE_MEMORY_HOST);
        HYPRE_SetExecutionPolicy(HYPRE_EXEC_HOST);
    }
    initializedHere = true;
    ( void ) initializedHere;
}

bool hypreCheck(HYPRE_Int error, const char *operation) {
    if ( error == 0 ) {
        return true;
    }
    std :: cerr << "Hypre error " << error << " during " << operation << std :: endl;
    HYPRE_ClearAllErrors();
    return false;
}

HYPRE_Int hyprePreconditionerAlreadySetup(HYPRE_Solver, HYPRE_ParCSRMatrix, HYPRE_ParVector, HYPRE_ParVector) {
    return 0;
}

bool makeHypreVector(
    const std :: vector< HYPRE_BigInt > &indices,
    const std :: vector< double > &values,
    HYPRE_IJVector &ijVector,
    HYPRE_ParVector &parVector
) {
    if ( indices.empty() ) {
        return false;
    }
    const HYPRE_BigInt first = indices.front();
    const HYPRE_BigInt last = indices.back();
    if ( !hypreCheck(HYPRE_IJVectorCreate(MPI_COMM_SELF, first, last, &ijVector), "HYPRE_IJVectorCreate") ) {
        return false;
    }
    if ( !hypreCheck(HYPRE_IJVectorSetObjectType(ijVector, HYPRE_PARCSR), "HYPRE_IJVectorSetObjectType") ) {
        return false;
    }
    if ( !hypreCheck(HYPRE_IJVectorInitialize(ijVector), "HYPRE_IJVectorInitialize") ) {
        return false;
    }
    if ( !hypreCheck(HYPRE_IJVectorSetValues(ijVector, static_cast< HYPRE_Int >( indices.size() ), indices.data(), values.data() ), "HYPRE_IJVectorSetValues") ) {
        return false;
    }
    if ( !hypreCheck(HYPRE_IJVectorAssemble(ijVector), "HYPRE_IJVectorAssemble") ) {
        return false;
    }
    if ( !hypreCheck(HYPRE_IJVectorGetObject(ijVector, reinterpret_cast< void ** >( &parVector ) ), "HYPRE_IJVectorGetObject") ) {
        return false;
    }
    return true;
}

}

//////////////////////////////////////////////////////////
HypreBoomerAMGCGSolver :: HypreBoomerAMGCGSolver() {
    name = "HypreBoomerAMGCGSolver";
    impl = std :: make_unique< Impl >();
    lastIterations = -1;
    lastError = -1.;
    lastTrueRelativeResidual = -1.;
}

//////////////////////////////////////////////////////////
HypreBoomerAMGCGSolver :: ~HypreBoomerAMGCGSolver() {}

//////////////////////////////////////////////////////////
void HypreBoomerAMGCGSolver :: setOptions(const HypreBoomerAMGOptions &opts) {
    options = opts;
}

//////////////////////////////////////////////////////////
void HypreBoomerAMGCGSolver :: setRuntimeTolerance(double tolerance, double trueTolerance) {
    const double selectedTolerance = trueTolerance > 0. ? trueTolerance : tolerance;
    if ( selectedTolerance > 0. ) {
        options.tolerance = selectedTolerance;
        if ( impl && impl->solver ) {
            HYPRE_ParCSRPCGSetTol(impl->solver, options.tolerance);
        }
    }
}

//////////////////////////////////////////////////////////
void HypreBoomerAMGCGSolver :: setElasticDofMap(ElasticDofMap map) {
    elasticMap = std :: move(map);
}

//////////////////////////////////////////////////////////
bool HypreBoomerAMGCGSolver :: analyzePattern(const CoordinateIndexedSparseMatrix &A) {
    ( void ) A;
    return true;
}

//////////////////////////////////////////////////////////
bool HypreBoomerAMGCGSolver :: factorize(const CoordinateIndexedSparseMatrix &A) {
    if ( A.rows() <= 0 ) {
        return true;
    }
    ensureHypreRuntimeInitialized();

    impl->clearSolverSetup();
    impl->clearMatrix();

    const bool useElasticLift = elasticMap.isValid();
    impl->reducedMatrix = A;
    impl->activeMatrix = useElasticLift ? liftMatrixToElasticSpace(A, elasticMap) : A;

    if ( options.checkMatrix ) {
        CoordinateIndexedSparseMatrix diff = impl->activeMatrix - CoordinateIndexedSparseMatrix( impl->activeMatrix.transpose() );
        const double matrixNorm = impl->activeMatrix.norm();
        const double symmetryRelativeNorm = matrixNorm > 0. ? diff.norm() / matrixNorm : 0.;
        long long nonPositiveDiagonal = 0;
        double minDiagonal = std :: numeric_limits< double > :: infinity();
        double maxDiagonal = -std :: numeric_limits< double > :: infinity();
        for ( int i = 0; i < impl->activeMatrix.rows(); ++i ) {
            const double value = impl->activeMatrix.coeff(i, i);
            minDiagonal = std :: min(minDiagonal, value);
            maxDiagonal = std :: max(maxDiagonal, value);
            if ( value <= 0. ) {
                ++nonPositiveDiagonal;
            }
        }
        std :: cout << "Hypre BoomerAMG-CG matrix check: symmetry_relative_norm=" << symmetryRelativeNorm
                    << ", min_diagonal=" << minDiagonal
                    << ", max_diagonal=" << maxDiagonal
                    << ", non_positive_diagonal_count=" << nonPositiveDiagonal << std :: endl;
    }

    ptrdiff_t rows = 0;
    std :: vector< ptrdiff_t >ptr;
    std :: vector< ptrdiff_t >col;
    std :: vector< double >val;
    makeRowMajorCrs(impl->activeMatrix, rows, ptr, col, val);
    if ( rows <= 0 || rows > std :: numeric_limits< HYPRE_Int > :: max() ) {
        std :: cerr << "Hypre BoomerAMG-CG supports positive row counts fitting HYPRE_Int; rows=" << rows << std :: endl;
        return false;
    }

    impl->firstRow = 0;
    impl->lastRow = rows - 1;
    impl->indices.resize(rows);
    for ( ptrdiff_t row = 0; row < rows; row++ ) {
        impl->indices [ row ] = row;
    }

    if ( !hypreCheck(HYPRE_IJMatrixCreate(MPI_COMM_SELF, impl->firstRow, impl->lastRow, impl->firstRow, impl->lastRow, &impl->ijMatrix), "HYPRE_IJMatrixCreate") ) {
        return false;
    }
    if ( !hypreCheck(HYPRE_IJMatrixSetObjectType(impl->ijMatrix, HYPRE_PARCSR), "HYPRE_IJMatrixSetObjectType") ) {
        return false;
    }
    std :: vector< HYPRE_Int >rowSizes(rows, 0);
    for ( ptrdiff_t row = 0; row < rows; row++ ) {
        rowSizes [ row ] = static_cast< HYPRE_Int >( ptr [ row + 1 ] - ptr [ row ] );
    }
    HYPRE_IJMatrixSetRowSizes(impl->ijMatrix, rowSizes.data() );
    HYPRE_IJMatrixSetOMPFlag(impl->ijMatrix, 1);
    if ( !hypreCheck(HYPRE_IJMatrixInitialize(impl->ijMatrix), "HYPRE_IJMatrixInitialize") ) {
        return false;
    }

    std :: vector< HYPRE_BigInt >rowColumns;
    for ( ptrdiff_t row = 0; row < rows; row++ ) {
        HYPRE_Int rowNnz = static_cast< HYPRE_Int >( ptr [ row + 1 ] - ptr [ row ] );
        if ( rowNnz == 0 ) {
            continue;
        }
        rowColumns.resize(rowNnz);
        for ( HYPRE_Int j = 0; j < rowNnz; j++ ) {
            rowColumns [ j ] = col [ ptr [ row ] + j ];
        }
        const HYPRE_BigInt rowId = row;
        if ( !hypreCheck(HYPRE_IJMatrixSetValues(impl->ijMatrix, 1, &rowNnz, &rowId, rowColumns.data(), val.data() + ptr [ row ] ), "HYPRE_IJMatrixSetValues") ) {
            return false;
        }
    }

    if ( !hypreCheck(HYPRE_IJMatrixAssemble(impl->ijMatrix), "HYPRE_IJMatrixAssemble") ) {
        return false;
    }
    if ( !hypreCheck(HYPRE_IJMatrixGetObject(impl->ijMatrix, reinterpret_cast< void ** >( &impl->parMatrix ) ), "HYPRE_IJMatrixGetObject") ) {
        return false;
    }

    impl->dofFunc.clear();
    if ( useElasticLift && elasticMap.blockSize > 1 ) {
        impl->dofFunc.resize(rows, 0);
        for ( ptrdiff_t row = 0; row < rows; row++ ) {
            impl->dofFunc [ row ] = static_cast< HYPRE_Int >( row % elasticMap.blockSize );
        }
    }
    if ( !hypreCheck(HYPRE_ParCSRPCGCreate(MPI_COMM_SELF, &impl->solver), "HYPRE_ParCSRPCGCreate") ) {
        impl->clearSolverSetup();
        return false;
    }
    HYPRE_ParCSRPCGSetTol(impl->solver, options.tolerance);
    HYPRE_ParCSRPCGSetAbsoluteTol(impl->solver, 0.);
    HYPRE_ParCSRPCGSetMaxIter(impl->solver, static_cast< HYPRE_Int >( options.maxIterations ) );
    HYPRE_ParCSRPCGSetTwoNorm(impl->solver, 1);
    HYPRE_ParCSRPCGSetRelChange(impl->solver, 0);
    HYPRE_PCGSetSkipBreak(impl->solver, options.skipBreak);
    HYPRE_PCGSetFlex(impl->solver, options.flex ? 1 : 0);
    HYPRE_PCGSetRecomputeResidual(impl->solver, options.recomputeResidual ? 1 : 0);
    if ( options.recomputeResidualPeriod > 0 ) {
        HYPRE_PCGSetRecomputeResidualP(impl->solver, options.recomputeResidualPeriod);
    }
    HYPRE_ParCSRPCGSetPrintLevel(impl->solver, options.printLevel);
    HYPRE_ParCSRPCGSetLogging(impl->solver, 1);

    if ( !hypreCheck(HYPRE_BoomerAMGCreate(&impl->precond), "HYPRE_BoomerAMGCreate") ) {
        impl->clearSolverSetup();
        return false;
    }
    HYPRE_BoomerAMGSetTol(impl->precond, 0.);
    HYPRE_BoomerAMGSetMaxIter(impl->precond, options.boomerMaxIterations);
    HYPRE_BoomerAMGSetCoarsenType(impl->precond, options.coarsenType);
    HYPRE_BoomerAMGSetInterpType(impl->precond, options.interpType);
    HYPRE_BoomerAMGSetStrongThreshold(impl->precond, options.strongThreshold);
    HYPRE_BoomerAMGSetRelaxType(impl->precond, options.relaxType);
    HYPRE_BoomerAMGSetRelaxOrder(impl->precond, options.relaxOrder);
    HYPRE_BoomerAMGSetPMaxElmts(impl->precond, options.pMaxElmts);
    HYPRE_BoomerAMGSetAggNumLevels(impl->precond, options.aggNumLevels);
    HYPRE_BoomerAMGSetPrintLevel(impl->precond, options.printLevel);
    if ( options.nonGalerkinTol >= 0. ) {
        HYPRE_Real tolerances [ 3 ] = { 0., static_cast< HYPRE_Real >( options.nonGalerkinTol ), static_cast< HYPRE_Real >( options.nonGalerkinTol ) };
        HYPRE_BoomerAMGSetNonGalerkTol(impl->precond, 3, tolerances);
    }
    if ( useElasticLift && elasticMap.blockSize > 1 ) {
        const HYPRE_Int numFunctions = static_cast< HYPRE_Int >( options.numFunctions > 0 ? options.numFunctions : elasticMap.blockSize );
        HYPRE_BoomerAMGSetNumFunctions(impl->precond, numFunctions);
        if ( options.useDofFunctions ) {
            std :: cout << "Hypre BoomerAMG-CG: explicit dof_func is disabled; using hypre's default cyclic node-major mapping." << std :: endl;
        }
        if ( options.nodal > 0 ) {
            HYPRE_BoomerAMGSetNodal(impl->precond, options.nodal);
        }
        if ( options.useInterpVectors && elasticMap.nearNullspaceColumns > 0 && elasticMap.nearNullspace.size() == size_t( rows * elasticMap.nearNullspaceColumns ) ) {
            std :: cout << "Hypre BoomerAMG-CG: interpolation vectors requested but skipped pending an ownership-safe wrapper." << std :: endl;
        }
    }

    HYPRE_ParCSRPCGSetPrecond(impl->solver, HYPRE_BoomerAMGSolve, hyprePreconditionerAlreadySetup, impl->precond);

    std :: vector< double >setupRhs(rows, 1.);
    std :: vector< double >setupGuess(rows, 0.);
    if ( !makeHypreVector(impl->indices, setupRhs, impl->setupIjB, impl->setupParB) ) {
        impl->clearSolverSetup();
        return false;
    }
    if ( !makeHypreVector(impl->indices, setupGuess, impl->setupIjX, impl->setupParX) ) {
        impl->clearSolverSetup();
        return false;
    }
    if ( !hypreCheck(HYPRE_BoomerAMGSetup(impl->precond, impl->parMatrix, impl->setupParB, impl->setupParX), "HYPRE_BoomerAMGSetup") ) {
        impl->clearSolverSetup();
        return false;
    }
    if ( !hypreCheck(HYPRE_ParCSRPCGSetup(impl->solver, impl->parMatrix, impl->setupParB, impl->setupParX), "HYPRE_ParCSRPCGSetup") ) {
        impl->clearSolverSetup();
        return false;
    }

    impl->matrixReady = true;

    std :: cout << "Hypre BoomerAMG-CG: matrix and preconditioner setup complete, rows=" << rows
                << ", nnz=" << impl->activeMatrix.nonZeros()
                << ", tolerance=" << options.tolerance
                << ", max_iterations=" << options.maxIterations
                << ", coarsen_type=" << options.coarsenType
                << ", interp_type=" << options.interpType
                << ", strong_threshold=" << options.strongThreshold
                << ", nodal=" << options.nodal
                << ", flex=" << ( options.flex ? 1 : 0 )
                << ", skip_break=" << options.skipBreak
                << ", block_size=" << ( useElasticLift ? elasticMap.blockSize : 1 ) << std :: endl;
    return true;
}

//////////////////////////////////////////////////////////
bool HypreBoomerAMGCGSolver :: solve(Vector &x, const Vector &b) {
    if ( b.size() == 0 ) {
        lastIterations = 0;
        lastError = 0.;
        lastTrueRelativeResidual = 0.;
        return true;
    }
    if ( !impl->matrixReady || !impl->parMatrix || !impl->solver || !impl->precond ) {
        std :: cerr << "Hypre BoomerAMG-CG Error: solve called before factorize" << std :: endl;
        return false;
    }

    Vector activeB = elasticMap.isValid() ? liftVectorToElasticSpace(b, elasticMap) : b;
    std :: vector< double >rhs(activeB.data(), activeB.data() + activeB.size() );
    std :: vector< double >guess(activeB.size(), 0.);
    std :: vector< double >answer(activeB.size(), 0.);

    HYPRE_IJVector ijB = nullptr;
    HYPRE_IJVector ijX = nullptr;
    HYPRE_ParVector parB = nullptr;
    HYPRE_ParVector parX = nullptr;

    auto cleanup = [&]() {
        if ( ijB ) {
            HYPRE_IJVectorDestroy(ijB);
        }
        if ( ijX ) {
            HYPRE_IJVectorDestroy(ijX);
        }
    };

    if ( !makeHypreVector(impl->indices, rhs, ijB, parB) ) {
        cleanup();
        return false;
    }
    if ( !makeHypreVector(impl->indices, guess, ijX, parX) ) {
        cleanup();
        return false;
    }

    if ( !hypreCheck(HYPRE_ParCSRPCGSetup(impl->solver, impl->parMatrix, parB, parX), "HYPRE_ParCSRPCGSetup") ) {
        cleanup();
        return false;
    }
    HYPRE_Int solveError = HYPRE_ParCSRPCGSolve(impl->solver, impl->parMatrix, parB, parX);
    if ( solveError != 0 ) {
        HYPRE_ClearAllErrors();
    }

    HYPRE_Int iterations = -1;
    HYPRE_Real hypreResidual = -1.;
    HYPRE_ParCSRPCGGetNumIterations(impl->solver, &iterations);
    HYPRE_ParCSRPCGGetFinalRelativeResidualNorm(impl->solver, &hypreResidual);
    HYPRE_IJVectorGetValues(ijX, static_cast< HYPRE_Int >( impl->indices.size() ), impl->indices.data(), answer.data() );

    Vector activeX = Eigen :: Map< Vector >(answer.data(), answer.size() );
    x = elasticMap.isValid() ? restrictVectorFromElasticSpace(activeX, elasticMap) : activeX;
    lastIterations = iterations;
    lastError = trueRelativeResidual(impl->reducedMatrix, x, b);
    lastTrueRelativeResidual = lastError;
    bool ok = ( solveError == 0 || solveError == HYPRE_ERROR_CONV ) && lastTrueRelativeResidual <= options.tolerance;
    if ( !ok ) {
        if ( solveError != 0 ) {
            std :: cerr << "Hypre error " << solveError << " during HYPRE_ParCSRPCGSolve" << std :: endl;
        }
        std :: cerr << "Hypre BoomerAMG-CG warning: performed " << iterations
                    << " iterations and reached hypre residual " << hypreResidual
                    << ", true relative residual " << lastTrueRelativeResidual
                    << ", required tolerance is " << options.tolerance << std :: endl;
    }

    cleanup();
    return ok;
}
#endif

//////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////
// LDLT SOLVER
//////////////////////////////////////////////////////////

LDLTSolver :: LDLTSolver() {
    name = "LDLTSolver";
}

//////////////////////////////////////////////////////////
LDLTSolver :: ~LDLTSolver() {}

//////////////////////////////////////////////////////////
bool LDLTSolver :: factorize(const CoordinateIndexedSparseMatrix &A) {
#if PRINT_DEBUG_TIME
    auto start = std :: chrono :: system_clock :: now();
#endif
    if ( A.rows() > 0 ) {
        ldlt.factorize(A);
    }
#if PRINT_DEBUG_TIME
    now = std :: chrono :: system_clock :: now();

    elapsed_seconds = now - start;
    std :: cout << "linalg solver decomposition duration: " << convertTimeToString_(elapsed_seconds) << std :: endl;
    cout.flush();
#endif
    return true;
}


//////////////////////////////////////////////////////////
bool LDLTSolver :: analyzePattern(const CoordinateIndexedSparseMatrix &A) {
#if PRINT_DEBUG_TIME
    auto start = std :: chrono :: system_clock :: now();
#endif

    if ( A.rows() > 0 ) {
        ldlt.analyzePattern(A);
    }

#if PRINT_DEBUG_TIME
    now = std :: chrono :: system_clock :: now();

    elapsed_seconds = now - start;
    std :: cout << "linalg solver decomposition duration: " << convertTimeToString_(elapsed_seconds) << std :: endl;
    cout.flush();
#endif
    return true;
}



//////////////////////////////////////////////////////////
bool LDLTSolver :: solve(Vector &x, const Vector &b) {
#if PRINT_DEBUG_TIME
    auto start = std :: chrono :: system_clock :: now();
#endif

    if ( b.size() > 0 ) {
        x = ldlt.solve(b);
    }

#if PRINT_DEBUG_TIME
    now = std :: chrono :: system_clock :: now();

    elapsed_seconds = now - start;
    std :: cout << "linalg solver duration: " << convertTimeToString_(elapsed_seconds) << std :: endl;
    cout.flush();
#endif
    return true;
}


//////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////
// LU SOLVER
//////////////////////////////////////////////////////////

LUSolver :: LUSolver() {
    name = "LUSolver";
}

//////////////////////////////////////////////////////////
LUSolver :: ~LUSolver() {}

//////////////////////////////////////////////////////////
bool LUSolver :: factorize(const CoordinateIndexedSparseMatrix &A) {
#if PRINT_DEBUG_TIME
    auto start = std :: chrono :: system_clock :: now();
#endif
    if ( A.rows() > 0 ) {
        lu.factorize(A);
    }

#if PRINT_DEBUG_TIME
    now = std :: chrono :: system_clock :: now();

    elapsed_seconds = now - start;
    std :: cout << "linalg solver decomposition duration: " << convertTimeToString_(elapsed_seconds) << std :: endl;
    cout.flush();
#endif
    return true;
}

//////////////////////////////////////////////////////////
bool LUSolver :: analyzePattern(const CoordinateIndexedSparseMatrix &A) {
#if PRINT_DEBUG_TIME
    auto start = std :: chrono :: system_clock :: now();
#endif

    if ( A.rows() > 0 ) {
        lu.analyzePattern(A);
    }

#if PRINT_DEBUG_TIME
    now = std :: chrono :: system_clock :: now();

    elapsed_seconds = now - start;
    std :: cout << "linalg solver decomposition duration: " << convertTimeToString_(elapsed_seconds) << std :: endl;
    cout.flush();
#endif
    return true;
}


//////////////////////////////////////////////////////////
bool LUSolver :: solve(Vector &x, const Vector &b) {
#if PRINT_DEBUG_TIME
    auto start = std :: chrono :: system_clock :: now();
#endif

    if ( b.size() > 0 ) {
        x = lu.solve(b);
    }

#if PRINT_DEBUG_TIME
    now = std :: chrono :: system_clock :: now();

    elapsed_seconds = now - start;
    std :: cout << "linalg solver duration: " << convertTimeToString_(elapsed_seconds) << std :: endl;
    cout.flush();
#endif
    return true;
}


//////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////
// SUPERLU SOLVER
//////////////////////////////////////////////////////////

#ifdef SUPERLU_FOUND
SuperLUSolver :: SuperLUSolver() {
    name = "SuperLUSolver";
}

//////////////////////////////////////////////////////////
SuperLUSolver :: ~SuperLUSolver() {}

//////////////////////////////////////////////////////////
bool SuperLUSolver :: factorize(const CoordinateIndexedSparseMatrix &A) {
 #if PRINT_DEBUG_TIME
    auto start = std :: chrono :: system_clock :: now();
 #endif
    if ( A.rows() > 0 ) {
        superlu.factorize(A);
    }

 #if PRINT_DEBUG_TIME
    now = std :: chrono :: system_clock :: now();

    elapsed_seconds = now - start;
    std :: cout << "linalg solver decomposition duration: " << convertTimeToString_(elapsed_seconds) << std :: endl;
    cout.flush();
 #endif
    return true;
}


//////////////////////////////////////////////////////////
bool SuperLUSolver :: analyzePattern(const CoordinateIndexedSparseMatrix &A) {
 #if PRINT_DEBUG_TIME
    auto start = std :: chrono :: system_clock :: now();
 #endif
    if ( A.rows() > 0 ) {
        superlu.analyzePattern(A);
    }

 #if PRINT_DEBUG_TIME
    now = std :: chrono :: system_clock :: now();

    elapsed_seconds = now - start;
    std :: cout << "linalg solver decomposition duration: " << convertTimeToString_(elapsed_seconds) << std :: endl;
    cout.flush();
 #endif
    return true;
}

//////////////////////////////////////////////////////////
bool SuperLUSolver :: solve(Vector &x, const Vector &b) {
 #if PRINT_DEBUG_TIME
    auto start = std :: chrono :: system_clock :: now();
 #endif

    if ( b.size() > 0 ) {
        x = superlu.solve(b);
    }

 #if PRINT_DEBUG_TIME
    now = std :: chrono :: system_clock :: now();

    elapsed_seconds = now - start;
    std :: cout << "linalg solver duration: " << convertTimeToString_(elapsed_seconds) << std :: endl;
    cout.flush();
 #endif
    return true;
}
#endif

//////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////
// LLT SOLVER
//////////////////////////////////////////////////////////

LLTSolver :: LLTSolver() {
    name = "LLTSolver";
}

//////////////////////////////////////////////////////////
LLTSolver :: ~LLTSolver() {}

//////////////////////////////////////////////////////////
bool LLTSolver :: factorize(const CoordinateIndexedSparseMatrix &A) {
#if PRINT_DEBUG_TIME
    auto start = std :: chrono :: system_clock :: now();
#endif

    if ( A.rows() > 0 ) {
        llt.factorize(A);
    }

#if PRINT_DEBUG_TIME
    now = std :: chrono :: system_clock :: now();

    elapsed_seconds = now - start;
    std :: cout << "linalg solver decomposition duration: " << convertTimeToString_(elapsed_seconds) << std :: endl;
    cout.flush();
#endif
    return true;
}


//////////////////////////////////////////////////////////
bool LLTSolver :: analyzePattern(const CoordinateIndexedSparseMatrix &A) {
#if PRINT_DEBUG_TIME
    auto start = std :: chrono :: system_clock :: now();
#endif

    if ( A.rows() > 0 ) {
        llt.analyzePattern(A);
    }

#if PRINT_DEBUG_TIME
    now = std :: chrono :: system_clock :: now();

    elapsed_seconds = now - start;
    std :: cout << "linalg solver decomposition duration: " << convertTimeToString_(elapsed_seconds) << std :: endl;
    cout.flush();
#endif
    return true;
}


//////////////////////////////////////////////////////////
bool LLTSolver :: solve(Vector &x, const Vector &b) {
#if PRINT_DEBUG_TIME
    auto start = std :: chrono :: system_clock :: now();
#endif

    if ( b.size() > 0 ) {
        x = llt.solve(b);
    }

#if PRINT_DEBUG_TIME
    now = std :: chrono :: system_clock :: now();

    elapsed_seconds = now - start;
    std :: cout << "linalg solver duration: " << convertTimeToString_(elapsed_seconds) << std :: endl;
    cout.flush();
#endif
    return true;
}

//////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////
// PARDISO LU SOLVER
//////////////////////////////////////////////////////////

#ifdef PARDISO_FOUND
PardisoLUSolver :: PardisoLUSolver() {
    name = "PardisoLUSolver";

    // Set MKL to use all available threads for Pardiso
    // This can be overridden by MKL_NUM_THREADS environment variable
    int max_threads = mkl_get_max_threads();

    // Check if MKL_NUM_THREADS or OMP_NUM_THREADS is set
    const char *mkl_threads_env = std :: getenv("MKL_NUM_THREADS");
    const char *omp_threads_env = std :: getenv("OMP_NUM_THREADS");

    if ( mkl_threads_env != nullptr ) {
        int requested_threads = std :: atoi(mkl_threads_env);
        if ( requested_threads > 0 && requested_threads <= max_threads ) {
            mkl_set_num_threads(requested_threads);
            std :: cout << "PardisoLU: Using " << requested_threads
                        << " threads (from MKL_NUM_THREADS)" << std :: endl;
        }
    } else if ( omp_threads_env != nullptr ) {
        int requested_threads = std :: atoi(omp_threads_env);
        if ( requested_threads > 0 && requested_threads <= max_threads ) {
            mkl_set_num_threads(requested_threads);
            std :: cout << "PardisoLU: Using " << requested_threads
                        << " threads (from OMP_NUM_THREADS)" << std :: endl;
        }
    } else {
        // Use all available threads by default
        mkl_set_num_threads(max_threads);
        std :: cout << "PardisoLU: Using " << max_threads
                    << " threads (default - all available)" << std :: endl;
    }
}

//////////////////////////////////////////////////////////
PardisoLUSolver :: ~PardisoLUSolver() {}

//////////////////////////////////////////////////////////
bool PardisoLUSolver :: factorize(const CoordinateIndexedSparseMatrix &A) {
 #if PRINT_DEBUG_TIME
    auto start = std :: chrono :: system_clock :: now();
 #endif
    if ( A.rows() > 0 ) {
        pardiso.factorize(A);
    }

 #if PRINT_DEBUG_TIME
    now = std :: chrono :: system_clock :: now();

    elapsed_seconds = now - start;
    std :: cout << "linalg solver decomposition duration: " << convertTimeToString_(elapsed_seconds) << std :: endl;
    cout.flush();
 #endif
    return true;
}


//////////////////////////////////////////////////////////
bool PardisoLUSolver :: analyzePattern(const CoordinateIndexedSparseMatrix &A) {
 #if PRINT_DEBUG_TIME
    auto start = std :: chrono :: system_clock :: now();
 #endif
    if ( A.rows() > 0 ) {
        pardiso.analyzePattern(A);
    }

 #if PRINT_DEBUG_TIME
    now = std :: chrono :: system_clock :: now();

    elapsed_seconds = now - start;
    std :: cout << "linalg solver decomposition duration: " << convertTimeToString_(elapsed_seconds) << std :: endl;
    cout.flush();
 #endif
    return true;
}

//////////////////////////////////////////////////////////
bool PardisoLUSolver :: solve(Vector &x, const Vector &b) {
 #if PRINT_DEBUG_TIME
    auto start = std :: chrono :: system_clock :: now();
 #endif

    if ( b.size() > 0 ) {
        x = pardiso.solve(b);
    }

 #if PRINT_DEBUG_TIME
    now = std :: chrono :: system_clock :: now();

    elapsed_seconds = now - start;
    std :: cout << "linalg solver duration: " << convertTimeToString_(elapsed_seconds) << std :: endl;
    cout.flush();
 #endif
    return true;
}

//////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////
// PARDISO LDLT SOLVER
//////////////////////////////////////////////////////////

PardisoLDLTSolver :: PardisoLDLTSolver() {
    name = "PardisoLDLTSolver";

    // Set MKL to use all available threads for Pardiso
    int max_threads = mkl_get_max_threads();

    const char *mkl_threads_env = std :: getenv("MKL_NUM_THREADS");
    const char *omp_threads_env = std :: getenv("OMP_NUM_THREADS");

    if ( mkl_threads_env != nullptr ) {
        int requested_threads = std :: atoi(mkl_threads_env);
        if ( requested_threads > 0 && requested_threads <= max_threads ) {
            mkl_set_num_threads(requested_threads);
            std :: cout << "PardisoLDLT: Using " << requested_threads
                        << " threads (from MKL_NUM_THREADS)" << std :: endl;
        }
    } else if ( omp_threads_env != nullptr ) {
        int requested_threads = std :: atoi(omp_threads_env);
        if ( requested_threads > 0 && requested_threads <= max_threads ) {
            mkl_set_num_threads(requested_threads);
            std :: cout << "PardisoLDLT: Using " << requested_threads
                        << " threads (from OMP_NUM_THREADS)" << std :: endl;
        }
    } else {
        mkl_set_num_threads(max_threads);
        std :: cout << "PardisoLDLT: Using " << max_threads
                    << " threads (default - all available)" << std :: endl;
    }
}

//////////////////////////////////////////////////////////
PardisoLDLTSolver :: ~PardisoLDLTSolver() {}

//////////////////////////////////////////////////////////
bool PardisoLDLTSolver :: factorize(const CoordinateIndexedSparseMatrix &A) {
 #if PRINT_DEBUG_TIME
    auto start = std :: chrono :: system_clock :: now();
 #endif
    if ( A.rows() > 0 ) {
        pardiso.factorize(A);
    }

 #if PRINT_DEBUG_TIME
    now = std :: chrono :: system_clock :: now();

    elapsed_seconds = now - start;
    std :: cout << "linalg solver decomposition duration: " << convertTimeToString_(elapsed_seconds) << std :: endl;
    cout.flush();
 #endif
    return true;
}


//////////////////////////////////////////////////////////
bool PardisoLDLTSolver :: analyzePattern(const CoordinateIndexedSparseMatrix &A) {
 #if PRINT_DEBUG_TIME
    auto start = std :: chrono :: system_clock :: now();
 #endif
    if ( A.rows() > 0 ) {
        pardiso.analyzePattern(A);
    }

 #if PRINT_DEBUG_TIME
    now = std :: chrono :: system_clock :: now();

    elapsed_seconds = now - start;
    std :: cout << "linalg solver decomposition duration: " << convertTimeToString_(elapsed_seconds) << std :: endl;
    cout.flush();
 #endif
    return true;
}

//////////////////////////////////////////////////////////
bool PardisoLDLTSolver :: solve(Vector &x, const Vector &b) {
 #if PRINT_DEBUG_TIME
    auto start = std :: chrono :: system_clock :: now();
 #endif

    if ( b.size() > 0 ) {
        x = pardiso.solve(b);
    }

 #if PRINT_DEBUG_TIME
    now = std :: chrono :: system_clock :: now();

    elapsed_seconds = now - start;
    std :: cout << "linalg solver duration: " << convertTimeToString_(elapsed_seconds) << std :: endl;
    cout.flush();
 #endif
    return true;
}

//////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////
// PARDISO LLT SOLVER
//////////////////////////////////////////////////////////

PardisoLLTSolver :: PardisoLLTSolver() {
    name = "PardisoLLTSolver";

    // Set MKL to use all available threads for Pardiso
    int max_threads = mkl_get_max_threads();

    const char *mkl_threads_env = std :: getenv("MKL_NUM_THREADS");
    const char *omp_threads_env = std :: getenv("OMP_NUM_THREADS");

    if ( mkl_threads_env != nullptr ) {
        int requested_threads = std :: atoi(mkl_threads_env);
        if ( requested_threads > 0 && requested_threads <= max_threads ) {
            mkl_set_num_threads(requested_threads);
            std :: cout << "PardisoLLT: Using " << requested_threads
                        << " threads (from MKL_NUM_THREADS)" << std :: endl;
        }
    } else if ( omp_threads_env != nullptr ) {
        int requested_threads = std :: atoi(omp_threads_env);
        if ( requested_threads > 0 && requested_threads <= max_threads ) {
            mkl_set_num_threads(requested_threads);
            std :: cout << "PardisoLLT: Using " << requested_threads
                        << " threads (from OMP_NUM_THREADS)" << std :: endl;
        }
    } else {
        mkl_set_num_threads(max_threads);
        std :: cout << "PardisoLLT: Using " << max_threads
                    << " threads (default - all available)" << std :: endl;
    }
}

//////////////////////////////////////////////////////////
PardisoLLTSolver :: ~PardisoLLTSolver() {}

//////////////////////////////////////////////////////////
bool PardisoLLTSolver :: factorize(const CoordinateIndexedSparseMatrix &A) {
 #if PRINT_DEBUG_TIME
    auto start = std :: chrono :: system_clock :: now();
 #endif
    if ( A.rows() > 0 ) {
        pardiso.factorize(A);
    }

 #if PRINT_DEBUG_TIME
    now = std :: chrono :: system_clock :: now();

    elapsed_seconds = now - start;
    std :: cout << "linalg solver decomposition duration: " << convertTimeToString_(elapsed_seconds) << std :: endl;
    cout.flush();
 #endif
    return true;
}


//////////////////////////////////////////////////////////
bool PardisoLLTSolver :: analyzePattern(const CoordinateIndexedSparseMatrix &A) {
 #if PRINT_DEBUG_TIME
    auto start = std :: chrono :: system_clock :: now();
 #endif
    if ( A.rows() > 0 ) {
        pardiso.analyzePattern(A);
    }

 #if PRINT_DEBUG_TIME
    now = std :: chrono :: system_clock :: now();

    elapsed_seconds = now - start;
    std :: cout << "linalg solver decomposition duration: " << convertTimeToString_(elapsed_seconds) << std :: endl;
    cout.flush();
 #endif
    return true;
}

//////////////////////////////////////////////////////////
bool PardisoLLTSolver :: solve(Vector &x, const Vector &b) {
 #if PRINT_DEBUG_TIME
    auto start = std :: chrono :: system_clock :: now();
 #endif

    if ( b.size() > 0 ) {
        x = pardiso.solve(b);
    }

 #if PRINT_DEBUG_TIME
    now = std :: chrono :: system_clock :: now();

    elapsed_seconds = now - start;
    std :: cout << "linalg solver duration: " << convertTimeToString_(elapsed_seconds) << std :: endl;
    cout.flush();
 #endif
    return true;
}
#endif

//////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////
// CHOLMOD SIMPLICIAL LLT SOLVER
//////////////////////////////////////////////////////////

#ifdef CHOLMOD_FOUND
CholmodLLTSolver :: CholmodLLTSolver() {
    name = "CholmodLLTSolver";
    std :: cout << "CholmodLLT: Using CHOLMOD simplicial LLT factorization" << std :: endl;
}

//////////////////////////////////////////////////////////
CholmodLLTSolver :: ~CholmodLLTSolver() {}

//////////////////////////////////////////////////////////
bool CholmodLLTSolver :: factorize(const CoordinateIndexedSparseMatrix &A) {
 #if PRINT_DEBUG_TIME
    auto start = std :: chrono :: system_clock :: now();
 #endif
    if ( A.rows() > 0 ) {
        cholmod.factorize(A);
    }

 #if PRINT_DEBUG_TIME
    now = std :: chrono :: system_clock :: now();

    elapsed_seconds = now - start;
    std :: cout << "linalg solver decomposition duration: " << convertTimeToString_(elapsed_seconds) << std :: endl;
    cout.flush();
 #endif
    return true;
}

//////////////////////////////////////////////////////////
bool CholmodLLTSolver :: analyzePattern(const CoordinateIndexedSparseMatrix &A) {
 #if PRINT_DEBUG_TIME
    auto start = std :: chrono :: system_clock :: now();
 #endif
    if ( A.rows() > 0 ) {
        cholmod.analyzePattern(A);
    }

 #if PRINT_DEBUG_TIME
    now = std :: chrono :: system_clock :: now();

    elapsed_seconds = now - start;
    std :: cout << "linalg solver decomposition duration: " << convertTimeToString_(elapsed_seconds) << std :: endl;
    cout.flush();
 #endif
    return true;
}

//////////////////////////////////////////////////////////
bool CholmodLLTSolver :: solve(Vector &x, const Vector &b) {
 #if PRINT_DEBUG_TIME
    auto start = std :: chrono :: system_clock :: now();
 #endif

    if ( b.size() > 0 ) {
        x = cholmod.solve(b);
    }

 #if PRINT_DEBUG_TIME
    now = std :: chrono :: system_clock :: now();

    elapsed_seconds = now - start;
    std :: cout << "linalg solver duration: " << convertTimeToString_(elapsed_seconds) << std :: endl;
    cout.flush();
 #endif
    return true;
}

//////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////
// CHOLMOD SIMPLICIAL LDLT SOLVER
//////////////////////////////////////////////////////////


CholmodLDLTSolver :: CholmodLDLTSolver() {
    name = "CholmodLDLT";
    cout << "Creating CholmodLDLTSolver (Simplicial LDLT)" << endl;

    // Access internal CHOLMOD settings using correct Eigen API
    cholmod_common &c = cholmod.cholmod();

    // Try different ordering methods (may help with some matrices)
    // AMD ordering (default)
    c.nmethods = 1;
    //c.method[0].ordering = CHOLMOD_AMD;

    // OR try METIS if available (sometimes better for large problems)
    c.method [ 0 ].ordering = CHOLMOD_METIS;

    // OR try nested dissection
    //c.method[0].ordering = CHOLMOD_NESDIS;

    // Ensure single-threaded
    c.useGPU = 0;

    //cout << "CHOLMOD LDLT using ordering: " << c.method[0].ordering << endl;
}

//////////////////////////////////////////////////////////
CholmodLDLTSolver :: ~CholmodLDLTSolver() {}

//////////////////////////////////////////////////////////
bool CholmodLDLTSolver :: factorize(const CoordinateIndexedSparseMatrix &A) {
 #if PRINT_DEBUG_TIME
    auto start = std :: chrono :: system_clock :: now();
 #endif

    if ( A.rows() > 0 ) {
 #if PRINT_DEBUG_TIME
        std :: cout << "CHOLMOD LDLT factorization:" << std :: endl;
        std :: cout << "  Matrix size: " << A.rows() << "x" << A.cols() << std :: endl;
        std :: cout << "  Non-zeros in A: " << A.nonZeros() << std :: endl;
 #endif
        cholmod.factorize(A);

 #if PRINT_DEBUG_TIME
        if ( cholmod.info() != Eigen :: Success ) {
            std :: cerr << "CHOLMOD factorization failed!" << std :: endl;
            return false;
        }
 #endif
    }

 #if PRINT_DEBUG_TIME
    auto now = std :: chrono :: system_clock :: now();
    auto elapsed = std :: chrono :: duration< double >(now - start);
    std :: cout << "  Factorization time: " << elapsed.count() << " seconds" << std :: endl;
 #endif

    return true;
}

//////////////////////////////////////////////////////////
bool CholmodLDLTSolver :: analyzePattern(const CoordinateIndexedSparseMatrix &A) {
 #if PRINT_DEBUG_TIME
    auto start = std :: chrono :: system_clock :: now();
 #endif

    if ( A.rows() > 0 ) {
 #if PRINT_DEBUG_TIME
        // Print matrix statistics
        std :: cout << "CHOLMOD LDLT analyzing matrix pattern:" << std :: endl;
        std :: cout << "  Size: " << A.rows() << " x " << A.cols() << std :: endl;
        std :: cout << "  Non-zeros: " << A.nonZeros() << std :: endl;
        std :: cout << "  Density: " << ( double ) A.nonZeros() / ( A.rows() * A.cols() ) * 100.0 << "%" << std :: endl;
 #endif
        cholmod.analyzePattern(A);
    }

 #if PRINT_DEBUG_TIME
    auto now = std :: chrono :: system_clock :: now();
    auto elapsed = std :: chrono :: duration< double >(now - start);
    std :: cout << "  Pattern analysis time: " << elapsed.count() << " seconds" << std :: endl;
 #endif

    return true;
}

//////////////////////////////////////////////////////////
bool CholmodLDLTSolver :: solve(Vector &x, const Vector &b) {
 #if PRINT_DEBUG_TIME
    auto start = std :: chrono :: system_clock :: now();
 #endif

    if ( b.size() > 0 ) {
        x = cholmod.solve(b);
    }

 #if PRINT_DEBUG_TIME
    now = std :: chrono :: system_clock :: now();

    elapsed_seconds = now - start;
    std :: cout << "linalg solver duration: " << convertTimeToString_(elapsed_seconds) << std :: endl;
    cout.flush();
 #endif
    return true;
}

//////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////
// CHOLMOD SUPERNODAL LLT SOLVER
//////////////////////////////////////////////////////////

CholmodSupernodalLLTSolver :: CholmodSupernodalLLTSolver() {
    name = "CholmodSupernodalLLTSolver";
    std :: cout << "CholmodSupernodalLLT: Using CHOLMOD supernodal LLT factorization" << std :: endl;
}

//////////////////////////////////////////////////////////
CholmodSupernodalLLTSolver :: ~CholmodSupernodalLLTSolver() {}

//////////////////////////////////////////////////////////
bool CholmodSupernodalLLTSolver :: factorize(const CoordinateIndexedSparseMatrix &A) {
 #if PRINT_DEBUG_TIME
    auto start = std :: chrono :: system_clock :: now();
 #endif
    if ( A.rows() > 0 ) {
        cholmod.factorize(A);
    }

 #if PRINT_DEBUG_TIME
    now = std :: chrono :: system_clock :: now();

    elapsed_seconds = now - start;
    std :: cout << "linalg solver decomposition duration: " << convertTimeToString_(elapsed_seconds) << std :: endl;
    cout.flush();
 #endif
    return true;
}

//////////////////////////////////////////////////////////
bool CholmodSupernodalLLTSolver :: analyzePattern(const CoordinateIndexedSparseMatrix &A) {
 #if PRINT_DEBUG_TIME
    auto start = std :: chrono :: system_clock :: now();
 #endif
    if ( A.rows() > 0 ) {
        cholmod.analyzePattern(A);
    }

 #if PRINT_DEBUG_TIME
    now = std :: chrono :: system_clock :: now();

    elapsed_seconds = now - start;
    std :: cout << "linalg solver decomposition duration: " << convertTimeToString_(elapsed_seconds) << std :: endl;
    cout.flush();
 #endif
    return true;
}

//////////////////////////////////////////////////////////
bool CholmodSupernodalLLTSolver :: solve(Vector &x, const Vector &b) {
 #if PRINT_DEBUG_TIME
    auto start = std :: chrono :: system_clock :: now();
 #endif

    if ( b.size() > 0 ) {
        x = cholmod.solve(b);
    }

 #if PRINT_DEBUG_TIME
    now = std :: chrono :: system_clock :: now();

    elapsed_seconds = now - start;
    std :: cout << "linalg solver duration: " << convertTimeToString_(elapsed_seconds) << std :: endl;
    cout.flush();
 #endif
    return true;
}
#endif

//////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////





/*
 *
 *
 *
 *  } else if ( solver_type == "EigenLDLT" ) {
 *      Eigen :: SimplicialLDLT< Eigen :: SparseMatrix< double > >simplicial_ldlt_solver;
 *      x = simplicial_ldlt_solver.compute(A).solve(b);
 *      //cout << "error " << ( A * x - b ).lpNorm< Eigen :: Infinity >() << endl;
 *      //result = ( A * x - b ).lpNorm< Eigen :: Infinity >() < precision;
 *      result = 1;
 *  } else if ( solver_type == "EigenLLT" ) {
 *      Eigen :: SimplicialLLT< Eigen :: SparseMatrix< double > >simplicial_llt_solver;
 *      x = simplicial_llt_solver.compute(A).solve(b);
 *      //result = ( A * x - b ).lpNorm< Eigen :: Infinity >() < precision;
 *      result = 1;
 *  } else if ( solver_type == "EigenSparseLU" ) {
 *      Eigen :: SparseLU< Eigen :: SparseMatrix< double >, Eigen :: COLAMDOrdering< int > >sparseLU_solver;
 *      sparseLU_solver.analyzePattern(A);
 *      sparseLU_solver.factorize(A);
 *
 *      x = sparseLU_solver.solve(b);
 *      //result = ( A * x - b ).lpNorm< Eigen :: Infinity >() < precision;
 *      result = 1;
 *
 *
 * }
 *
 *
 *
 */



bool StandAloneLinalgSolver(const CoordinateIndexedSparseMatrix &A, Vector &ddr, const Vector &f, double precision, double relmaxit, string solver_type) {
    std :: cout << "Using solver type: " << solver_type << std :: endl;
#if PRINT_DEBUG_TIME
    auto start = std :: chrono :: system_clock :: now();
#endif
    if ( f.size() == 0 ) {
        return true;                   // when problem is completely constrained (e.g. single facet)
    }
    std :: unique_ptr< LinAlgSolver >linalgsolver;
    if ( linalgsolver == nullptr ) {
        if ( solver_type == "EigenConj" ) {
            std :: unique_ptr< ConjGradSolver >cgs = std :: make_unique< ConjGradSolver >();
            cgs->setPrecisionAndRelMaxIters(precision, relmaxit);
            linalgsolver = std :: move(cgs);
            linalgsolver->analyzePattern(A);
        } else if  ( solver_type == "EigenLDLT" ) {
            linalgsolver = std :: make_unique< LDLTSolver >();
            linalgsolver->analyzePattern(A);
        } else if  ( solver_type == "EigenLLT" ) {
            linalgsolver = std :: make_unique< LLTSolver >();
            linalgsolver->analyzePattern(A);
        } else if  ( solver_type == "EigenLU" || solver_type == "EigenSparseLU" ) {
            linalgsolver = std :: make_unique< LUSolver >();
            linalgsolver->analyzePattern(A);
#ifdef SUPERLU_FOUND
        } else if  ( solver_type == "SuperLU" ) {
            linalgsolver = std :: make_unique< SuperLUSolver >();
            linalgsolver->analyzePattern(A);
#endif
#ifdef PARDISO_FOUND
        } else if  ( solver_type == "PardisoLU" ) {
            linalgsolver = std :: make_unique< PardisoLUSolver >();
            linalgsolver->analyzePattern(A);
        } else if  ( solver_type == "PardisoLDLT" ) {
            linalgsolver = std :: make_unique< PardisoLDLTSolver >();
            linalgsolver->analyzePattern(A);
        } else if  ( solver_type == "PardisoLLT" ) {
            linalgsolver = std :: make_unique< PardisoLLTSolver >();
            linalgsolver->analyzePattern(A);
#endif
#ifdef CHOLMOD_FOUND
        } else if  ( solver_type == "CholmodLLT" ) {
            linalgsolver = std :: make_unique< CholmodLLTSolver >();
            linalgsolver->analyzePattern(A);
        } else if  ( solver_type == "CholmodLDLT" ) {
            linalgsolver = std :: make_unique< CholmodLDLTSolver >();
            linalgsolver->analyzePattern(A);
        } else if  ( solver_type == "CholmodSupernodalLLT" ) {
            linalgsolver = std :: make_unique< CholmodSupernodalLLTSolver >();
            linalgsolver->analyzePattern(A);
#endif
        } else {
            cerr << "Solver type " << solver_type << " is not implemented" << endl;
            exit(1);
        }
    }
    linalgsolver->factorize(A);
    linalgsolver->solve(ddr, f);
    

#if PRINT_DEBUG_TIME
    now = std :: chrono :: system_clock :: now();

    elapsed_seconds = now - start;
    std :: cout << "linalg solver duration: " << convertTimeToString_(elapsed_seconds) << std :: endl;
    cout.flush();
#endif

    return true;
}


bool LinalgNonSymmetricSolver(const CoordinateIndexedSparseMatrix &A, Vector &x, const Vector &b, const Vector x0, double precision, double relmaxit) {
#if PRINT_DEBUG_TIME
    auto start = std :: chrono :: system_clock :: now();
#endif

    size_t Maxit = b.size() * relmaxit;

    //BiCGSTAB<SparseMatrix<double>, Eigen::IncompleteLUT<double>> bicg;
    Eigen :: BiCGSTAB< Eigen :: SparseMatrix< double > >bicg;
    bicg.setMaxIterations(Maxit);
    bicg.setTolerance(precision);

    bicg.compute(A);

    x = bicg.solveWithGuess(b, x0);
    //VectorXd cgx = bicg.solve(cgb);

#if PRINT_DEBUG_TIME
    now = std :: chrono :: system_clock :: now();

    elapsed_seconds = now - start;
    std :: cout << "linalg solver duration: " << convertTimeToString(elapsed_seconds) << " nit: " << bicg.iterations() << std :: endl;
    cout.flush();
#endif

    bool result = size_t( bicg.iterations() ) < Maxit;
    return result;
}

//sorterd eigenvalues and eigenvectors
bool LinalgEigenSolver(const Vector &A, Vector &eigenvalues, vector< Vector > &eigenvectors) {
    size_t ndim;
    bool sym;
    if ( A.size() == 3 ) { //2D
        ndim = 2;
        sym = true;
    } else if ( A.size() == 6 ) {      //3D
        ndim = 3;
        sym = true;
    } else if ( A.size() == 4 ) {      //2D
        ndim = 2;
        sym = false;
    } else if ( A.size() == 9 ) {      //3D
        ndim = 3;
        sym = false;
    } else {
        cerr << "Error: LinalgEigenSolver implemented only for vectorized matrices of size 2 or 3, submitted size " << A.size() << endl;
        exit(1);
    }
    Eigen :: MatrixXd mat = Eigen :: MatrixXd :: Zero(ndim, ndim);

    if ( ndim == 2 && sym ) {
        mat(0, 0) = A [ 0 ];
        mat(1, 1) = A [ 1 ];
        mat(1, 0) = mat(0, 1) = A [ 2 ];
    } else if ( ndim == 2 && !sym ) {
        mat(0, 0) = A [ 0 ];
        mat(1, 1) = A [ 1 ];
        mat(1, 0) = mat(0, 1) = ( A [ 2 ] + A [ 3 ] ) / 2.;
    } else if ( ndim == 3 && sym ) {
        mat(0, 0) = A [ 0 ];
        mat(1, 1) = A [ 1 ];
        mat(2, 2) = A [ 2 ];
        mat(2, 1) = mat(1, 2) = A [ 3 ];
        mat(2, 0) = mat(0, 2) = A [ 4 ];
        mat(1, 0) = mat(0, 1) = A [ 5 ];
    } else if ( ndim == 3 && !sym ) {
        mat(0, 0) = A [ 0 ];
        mat(1, 1) = A [ 1 ];
        mat(2, 2) = A [ 2 ];
        mat(2, 1) = mat(1, 2) = ( A [ 4 ] + A [ 3 ] ) / 2.;
        mat(2, 0) = mat(0, 2) = ( A [ 6 ] + A [ 5 ] ) / 2.;
        mat(1, 0) = mat(0, 1) = ( A [ 8 ] + A [ 7 ] ) / 2.;
    }

    return LinalgEigenSolver(mat, eigenvalues, eigenvectors);
}

bool LinalgEigenSolver(const Matrix &mat, Vector &eigenvalues, vector< Vector > &eigenvectors) {
    Eigen :: EigenSolver< Matrix >es;
    es.compute(mat, /* computeEigenvectors = */ true);

    unsigned ndim = mat.rows();
    vector< double >eigenvalsvector(ndim);
    eigenvalues.resize(ndim);
    eigenvectors.resize(ndim);

    for ( unsigned i = 0; i < ndim; i++ ) {
        eigenvalsvector [ i ] = ( es.eigenvalues() [ i ] ).real();
    }

    // initialize original index locations
    vector< size_t >idx(ndim);
    iota(idx.begin(), idx.end(), 0);
    stable_sort(idx.begin(), idx.end(), [ & eigenvalsvector ](size_t i1, size_t i2) {
        return eigenvalsvector [ i1 ] > eigenvalsvector [ i2 ];
    });

    for ( unsigned i = 0; i < ndim; i++ ) {
        eigenvalues [ i ] = eigenvalsvector [ idx [ i ] ];
        eigenvectors [ i ].resize(ndim);
        Eigen :: VectorXcd v = es.eigenvectors().col(idx [ i ]);
        for ( unsigned j = 0; j < ndim; j++ ) {
            eigenvectors [ i ] [ j ] = ( v [ j ] ).real();
        }
    }

    return true;
}

bool LinalgEigenSpectraSolver(const CoordinateIndexedSparseMatrix &mat, Vector &eigenvalues, Matrix &eigenvectors, int n_eigen_vals) {
    // Define matrix operation for Spectra
    Spectra :: SparseSymMatProd< double >op(mat);

    Spectra :: SymEigsSolver< Spectra :: SparseSymMatProd< double > >eigs(op, n_eigen_vals, 3 * n_eigen_vals);
    eigs.init();
    eigs.compute(Spectra :: SortRule :: SmallestAlge);  // Largest algebraic eigenvalues

    if ( eigs.info() == Spectra :: CompInfo :: Successful ) {
        eigenvalues = eigs.eigenvalues();
        eigenvectors = eigs.eigenvectors();
        std :: cout << "Eigenvalues:\n" << eigenvalues << std :: endl;
    } else {
        std :: cerr << "Eigenvalue computation failed!" << std :: endl;
    }

    return true;
}

bool LinalgEigenSpectraGENSolver(const CoordinateIndexedSparseMatrix &mat, const CoordinateIndexedSparseMatrix &matB, Vector &eigenvalues, Matrix &eigenvectors, int n_eigen_vals) {
    // Define matrix operation for Spectra
    Spectra :: SparseGenMatProd< double >op(mat);
    Spectra :: SparseCholesky< double >op_B(matB);

    // Define the solver: Solve Ax = λBx
    //Spectra :: SymGEigsSolver<Spectra :: SparseGenMatProd<double>, Spectra :: SparseSymMatProd<double>, Spectra :: GEigsMode::Cholesky> eigs(op, op_B, n_eigen_vals, 2 * n_eigen_vals);
    Spectra :: SymGEigsSolver< Spectra :: SparseGenMatProd< double >, Spectra :: SparseCholesky< double >, Spectra :: GEigsMode :: Cholesky >eigs(op, op_B, n_eigen_vals, 2 * n_eigen_vals);
    eigs.init();
    eigs.compute(Spectra :: SortRule :: SmallestAlge);  // Largest algebraic eigenvalues

    if ( eigs.info() == Spectra :: CompInfo :: Successful ) {
        eigenvalues = eigs.eigenvalues();
        eigenvectors = eigs.eigenvectors();
        std :: cout << "Eigenvalues:\n" << eigenvalues << std :: endl;
    } else {
        std :: cerr << "Eigenvalue computation failed!" << std :: endl;
    }

    return true;
}

bool LinalgLUSolver(const CoordinateIndexedSparseMatrix &A, Vector &x, const Vector &b) {
    ( void ) A;
    ( void ) x;
    ( void ) b;
    //return LinalgSymmetricSolver(A, x, b, x, 1e-12, 2, "EigenConj");


    //Eigen::SparseLU<Eigen::SparseMatrix<double>, Eigen::COLAMDOrdering<int> > solver;

    // fill A and b;
    // Compute the ordering permutation vector from the structural pattern of A
    //solver.analyzePattern(A);

    // Compute the numerical factorization
    //solver.factorize(A);

    //Use the factors to solve the linear system
    //x = solver.solve(b);

    return 0;
}




//JM: Coplanarity check of 4 points (for 3d faces)
double checkCoplanarity(const Point &ptA, const Point &ptB, const Point &ptC, const Point &ptD) {
    Point AB = ptB - ptA;
    Point AC = ptC - ptA;
    Point AD = ptD - ptA;
    //triple scalar product AB*(ACxAD) =>0
    double coplanarityError = AB.dot( AC.cross(AD) );
    return coplanarityError;
}


Matrix dyadicProduct(const Vector &a, const Vector &b) {
    return a * b.transpose();
}

double triArea2D(const Point *a, const Point *b, const Point *c) { //points in counter clockwise direction
    return 0.5 * ( a->x() * ( b->y() - c->y() ) + b->x() * ( c->y() - a->y() ) + c->x() * ( a->y() - b->y() ) );
}

double triArea3D(const Point *a, const Point *b, const Point *c) { //points
    Point AB = ( * b ) - ( * a );
    Point AC = ( * c ) - ( * a );
    return AB.cross(AC).norm() * 0.5;

    //return abs(0.5 * pow(pow( ( b->y() - a->y() ) * ( c->z() - a->z() ) - ( b->z() - a->z() ) * ( c->y() - a->y() ), 2 ) + pow( ( b->z() - a->z() ) * ( c->x() - a->x() ) - ( b->x() - a->x() ) * ( c->z() - a->z() ), 2 ) + pow( ( b->x() - a->x() ) * ( c->y() - a->y() ) - ( b->y() - a->y() ) * ( c->x() - a->x() ), 2 ), 0.5) );
}

double triInertia2D(const Point *a, const Point *b, const Point *c) { // Inertia of triangle relative to global [0, 0]
    double Ix;
    double Iy;
    if ( a->y() == b->y() && a->y() == c->y() ) {
        Ix = 0;
    } else if ( a->x() == b->x() && a->x() == c->x() ) {
        Ix = 0;
    } else if ( b->y() == c->y() ) {
        Ix = ( ( b->x() - a->x() ) / ( b->y() - a->y() ) * ( pow(b->y(), 4) - pow(a->y(), 4) ) / 4. ) -   ( a->y() * ( b->x() - a->x() ) / ( b->y() - a->y() ) * ( pow(b->y(), 3) - pow(a->y(), 3) ) / 3. ) - ( ( c->x() - a->x() ) / ( c->y() - a->y() ) * ( pow(b->y(), 4) - pow(a->y(), 4) ) / 4. ) + ( a->y() * ( c->x() - a->x() ) / ( c->y() - a->y() ) * ( pow(b->y(), 3) - pow(a->y(), 3) ) / 3. );
    } else if ( a->y() == b->y() ) {
        Ix = ( ( b->x() - c->x() ) / ( b->y() - c->y() ) * ( pow(c->y(), 4) - pow(a->y(), 4) ) / 4. ) - ( c->y() * ( b->x() - c->x() ) / ( b->y() - c->y() ) * ( pow(c->y(), 3) - pow(a->y(), 3) ) / 3. ) + ( c->x() * ( pow(c->y(), 3) - pow(a->y(), 3) ) / 3. ) - ( ( c->x() - a->x() ) / ( c->y() - a->y() ) * ( pow(c->y(), 4) - pow(a->y(), 4) ) / 4. ) + ( a->y() * ( c->x() - a->x() ) / ( c->y() - a->y() ) * ( pow(c->y(), 3) - pow(a->y(), 3) ) / 3. ) - ( a->x() * ( pow(c->y(), 3) - pow(a->y(), 3) ) / 3. );
    } else if ( a->y() == c->y() ) {
        Ix = ( ( c->x() - b->x() ) / ( c->y() - b->y() ) * ( pow(b->y(), 4) - pow(a->y(), 4) ) / 4. ) - ( b->y() * ( c->x() - b->x() ) / ( c->y() - b->y() ) * ( pow(b->y(), 3) - pow(a->y(), 3) ) / 3. ) + ( b->x() * ( pow(b->y(), 3) - pow(a->y(), 3) ) / 3. ) - ( ( b->x() - a->x() ) / ( b->y() - a->y() ) * ( pow(b->y(), 4) - pow(a->y(), 4) ) / 4. ) + ( a->y() * ( b->x() - a->x() ) / ( b->y() - a->y() ) * ( pow(b->y(), 3) - pow(a->y(), 3) ) / 3. ) - ( a->x() * ( pow(b->y(), 3) - pow(a->y(), 3) ) / 3. );
    } else {
        Ix = ( ( b->x() - c->x() ) / ( b->y() - c->y() ) * ( pow(c->y(), 4) - pow(a->y(), 4) ) / 4. ) - ( c->y() * ( b->x() - c->x() ) / ( b->y() - c->y() ) * ( pow(c->y(), 3) - pow(a->y(), 3) ) / 3. ) + ( c->x() * ( pow(c->y(), 3) - pow(a->y(), 3) ) / 3. ) - ( ( c->x() - a->x() ) / ( c->y() - a->y() ) * ( pow(c->y(), 4) - pow(a->y(), 4) ) / 4. ) + ( a->y() * ( c->x() - a->x() ) / ( c->y() - a->y() ) * ( pow(c->y(), 3) - pow(a->y(), 3) ) / 3. ) - ( a->x() * ( pow(c->y(), 3) - pow(a->y(), 3) ) / 3. ) - ( ( b->x() - c->x() ) / ( b->y() - c->y() ) * ( pow(b->y(), 4) - pow(a->y(), 4) ) / 4. ) + ( c->y() * ( b->x() - c->x() ) / ( b->y() - c->y() ) * ( pow(b->y(), 3) - pow(a->y(), 3) ) / 3. ) - ( c->x() * ( pow(b->y(), 3) - pow(a->y(), 3) ) / 3. ) + ( ( b->x() - a->x() ) / ( b->y() - a->y() ) * ( pow(b->y(), 4) - pow(a->y(), 4) ) / 4. ) - ( a->y() * ( b->x() - a->x() ) / ( b->y() - a->y() ) * ( pow(b->y(), 3) - pow(a->y(), 3) ) / 3. ) + ( a->x() * ( pow(b->y(), 3) - pow(a->y(), 3) ) / 3. );
    }

    if ( a->y() == b->y() && a->y() == c->y() ) {
        Iy = 0;
    } else if ( a->x() == b->x() && a->x() == c->x() ) {
        Iy = 0;
    } else if ( b->x() == c->x() ) {
        Iy = ( ( b->y() - a->y() ) / ( b->x() - a->x() ) * ( pow(b->x(), 4) - pow(a->x(), 4) ) / 4. ) - ( a->x() * ( b->y() - a->y() ) / ( b->x() - a->x() ) * ( pow(b->x(), 3) - pow(a->x(), 3) ) / 3. ) - ( ( c->y() - a->y() ) / ( c->x() - a->x() ) * ( pow(b->x(), 4) - pow(a->x(), 4) ) / 4. ) + ( a->x() * ( c->y() - a->y() ) / ( c->x() - a->x() ) * ( pow(b->x(), 3) - pow(a->x(), 3) ) / 3. );
    } else if ( a->x() == b->x() ) {
        Iy = ( ( b->y() - c->y() ) / ( b->x() - c->x() ) * ( pow(c->x(), 4) - pow(a->x(), 4) ) / 4. ) - ( c->x() * ( b->y() - c->y() ) / ( b->x() - c->x() ) * ( pow(c->x(), 3) - pow(a->x(), 3) ) / 3. ) + ( c->y() * ( pow(c->x(), 3) - pow(a->x(), 3) ) / 3. ) - ( ( c->y() - a->y() ) / ( c->x() - a->x() ) * ( pow(c->x(), 4) - pow(a->x(), 4) ) / 4. ) + ( a->x() * ( c->y() - a->y() ) / ( c->x() - a->x() ) * ( pow(c->x(), 3) - pow(a->x(), 3) ) / 3. ) - ( a->y() * ( pow(c->x(), 3) - pow(a->x(), 3) ) / 3. );
    } else if ( a->x() == c->x() ) {
        Iy = ( ( c->y() - b->y() ) / ( c->x() - b->x() ) * ( pow(b->x(), 4) - pow(a->x(), 4) ) / 4. ) - ( b->x() * ( c->y() - b->y() ) / ( c->x() - b->x() ) * ( pow(b->x(), 3) - pow(a->x(), 3) ) / 3. ) + ( b->y() * ( pow(b->x(), 3) - pow(a->x(), 3) ) / 3. ) - ( ( b->y() - a->y() ) / ( b->x() - a->x() ) * ( pow(b->x(), 4) - pow(a->x(), 4) ) / 4. ) + ( a->x() * ( b->y() - a->y() ) / ( b->x() - a->x() ) * ( pow(b->x(), 3) - pow(a->x(), 3) ) / 3. ) - ( a->y() * ( pow(b->x(), 3) - pow(a->x(), 3) ) / 3. );
    } else {
        Iy = ( ( b->y() - c->y() ) / ( b->x() - c->x() ) * ( pow(c->x(), 4) - pow(a->x(), 4) ) / 4. ) - ( c->x() * ( b->y() - c->y() ) / ( b->x() - c->x() ) * ( pow(c->x(), 3) - pow(a->x(), 3) ) / 3. ) + ( c->y() * ( pow(c->x(), 3) - pow(a->x(), 3) ) / 3. ) - ( ( c->y() - a->y() ) / ( c->x() - a->x() ) * ( pow(c->x(), 4) - pow(a->x(), 4) ) / 4. ) + ( a->x() * ( c->y() - a->y() ) / ( c->x() - a->x() ) * ( pow(c->x(), 3) - pow(a->x(), 3) ) / 3. ) - ( a->y() * ( pow(c->x(), 3) - pow(a->x(), 3) ) / 3. ) - ( ( b->y() - c->y() ) / ( b->x() - c->x() ) * ( pow(b->x(), 4) - pow(a->x(), 4) ) / 4. ) + ( c->x() * ( b->y() - c->y() ) / ( b->x() - c->x() ) * ( pow(b->x(), 3) - pow(a->x(), 3) ) / 3. ) - ( c->y() * ( pow(b->x(), 3) - pow(a->x(), 3) ) / 3. ) + ( ( b->y() - a->y() ) / ( b->x() - a->x() ) * ( pow(b->x(), 4) - pow(a->x(), 4) ) / 4. ) - ( a->x() * ( b->y() - a->y() ) / ( b->x() - a->x() ) * ( pow(b->x(), 3) - pow(a->x(), 3) ) / 3. ) + ( a->y() * ( pow(b->x(), 3) - pow(a->x(), 3) ) / 3. );
    }

    double Iz = abs(Ix) + abs(Iy);
    return Iz;
}

Matrix tetraInertia3D(const Point *a, const Point *b, const Point *c, const Point *d) { // Inertia relative to the centroid of the tetrahedron, where the centroid must be [0,0,0]
    // Jacobian of transformation to the reference tetrahedron
    double detJ = ( b->x() - a->x() ) * ( c->y() - a->y() ) * ( d->z() - a->z() ) +  ( c->x() - a->x() ) * ( d->y() - a->y() ) * ( b->z() - a->z() ) + ( d->x() - a->x() ) * ( b->y() - a->y() ) * ( c->z() - a->z() ) - ( b->x() - a->x() ) * ( d->y() - a->y() ) * ( c->z() - a->z() ) - ( c->x() - a->x() ) * ( b->y() - a->y() ) * ( d->z() - a->z() ) - ( d->x() - a->x() ) * ( c->y() - a->y() ) * ( b->z() - a->z() );
    // Inertia
    double Ixx = abs(detJ) * ( pow(a->y(), 2) + a->y() * b->y() + pow(b->y(), 2) + a->y() * c->y() + b->y() * c->y() + pow(c->y(), 2) + a->y() * d->y() + b->y() * d->y() + c->y() * d->y() + pow(d->y(), 2) + pow(a->z(), 2) + a->z() * b->z() + pow(b->z(), 2) + a->z() * c->z() + b->z() * c->z() + pow(c->z(), 2) + a->z() * d->z() + b->z() * d->z() + c->z() * d->z() + pow(d->z(), 2) ) / 60.;

    double Iyy = abs(detJ) * ( pow(a->x(), 2) + a->x() * b->x() + pow(b->x(), 2) + a->x() * c->x() + b->x() * c->x() + pow(c->x(), 2) + a->x() * d->x() + b->x() * d->x() + c->x() * d->x() + pow(d->x(), 2) + pow(a->z(), 2) + a->z() * b->z() + pow(b->z(), 2) + a->z() * c->z() + b->z() * c->z() + pow(c->z(), 2) + a->z() * d->z() + b->z() * d->z() + c->z() * d->z() + pow(d->z(), 2) ) / 60.;

    double Izz = abs(detJ) * ( pow(a->x(), 2) + a->x() * b->x() + pow(b->x(), 2) + a->x() * c->x() + b->x() * c->x() + pow(c->x(), 2) + a->x() * d->x() + b->x() * d->x() + c->x() * d->x() + pow(d->x(), 2) + pow(a->y(), 2) + a->y() * b->y() + pow(b->y(), 2) + a->y() * c->y() + b->y() * c->y() + pow(c->y(), 2) + a->y() * d->y() + b->y() * d->y() + c->y() * d->y() + pow(d->y(), 2) ) / 60.;

    double Ixz = abs(detJ) * ( 2. * a->x() * a->z() + b->x() * a->z() + c->x() * a->z() + d->x() * a->z() + a->x() * b->z() + 2. * b->x() * b->z() + c->x() * b->z() + d->x() * b->z() + a->x() * c->z() * b->x() * c->z() + 2. * c->x() * c->z() + d->x() * c->z() + a->x() * d->z() + b->x() * d->z() + c->x() * d->z() + 2. * d->x() * d->z() ) / 120.;

    double Iyz = abs(detJ) * ( 2. * a->y() * a->z() + b->y() * a->z() + c->y() * a->z() + d->y() * a->z() + a->y() * b->z() + 2. * b->y() * b->z() + c->y() * b->z() + d->y() * b->z() + a->y() * c->z() * b->y() * c->z() + 2. * c->y() * c->z() + d->y() * c->z() + a->y() * d->z() + b->y() * d->z() + c->y() * d->z() + 2. * d->y() * d->z() ) / 120.;

    double Ixy = abs(detJ) * ( 2. * a->x() * a->y() + b->x() * a->y() + c->x() * a->y() + d->x() * a->y() + a->x() * b->y() + 2. * b->x() * b->y() + c->x() * b->y() + d->x() * b->y() + a->x() * c->y() * b->x() * c->y() + 2. * c->x() * c->y() + d->x() * c->y() + a->x() * d->y() + b->x() * d->y() + c->x() * d->y() + 2. * d->x() * d->y() ) / 120.;

    // Inertia matrix
    Matrix I = Matrix :: Zero(3, 3);

    I(0, 0) = Ixx;
    I(1, 1) = Iyy;
    I(2, 2) = Izz;
    I(0, 1) = I(1, 0) = -Ixy;
    I(0, 2) = I(2, 0) = -Ixz;
    I(1, 2) = I(2, 1) = -Iyz;

    return I;
}

double tetraVolumeSigned(const Point *a, const Point *b, const Point *c, const Point *d) {
    return ( ( * d ) - ( * a ) ).dot( ( ( * b ) - ( * d ) ).cross( ( * c ) - ( * d ) ) ) / 6.;
}

bool is_positive_integer(const std :: string &s)
{
    return !s.empty() && std :: find_if(s.begin(),
                                        s.end(), [](unsigned char c) {
        return !std :: isdigit(c);
    }) == s.end();
}

void giveGaussIntegrationPointAndWeights(unsigned n, Vector &locs, Vector &weis) {
    locs.resize(n);
    weis.resize(n);
    if ( n == 1 ) {
        locs [ 0 ] = 0;
        weis [ 0 ] = 2;
    } else if ( n == 2 ) {
        locs [ 0 ] = -1. / sqrt(3);
        locs [ 1 ] = -locs [ 0 ];
        weis [ 0 ] = weis [ 1 ] = 1;
    } else if ( n == 3 ) {
        locs [ 0 ] = -sqrt(3. / 5.);
        locs [ 1 ] = 0;
        locs [ 2 ] = -locs [ 0 ];
        weis [ 0 ] = weis [ 2 ] = 5. / 9.;
        weis [ 1 ] = 8. / 9.;
    } else if ( n == 4 ) {
        locs [ 0 ] = -sqrt(3. / 7. + 2. / 7. * sqrt(6. / 5.) );
        locs [ 1 ] = -sqrt(3. / 7. - 2. / 7. * sqrt(6. / 5.) );
        locs [ 2 ] = -locs [ 1 ];
        locs [ 3 ] = -locs [ 0 ];
        weis [ 0 ] = weis [ 3 ] = ( 18. - sqrt(30.) ) / 36;
        weis [ 1 ] = weis [ 2 ] = ( 18. + sqrt(30.) ) / 36;
    } else if ( n == 5 ) {
        locs [ 0 ] = -sqrt(5. + 2. * sqrt(10. / 7.) ) / 3;
        locs [ 1 ] = -sqrt(5. - 2. * sqrt(10. / 7.) ) / 3;
        locs [ 2 ] = 0.;
        locs [ 3 ] = -locs [ 1 ];
        locs [ 4 ] = -locs [ 0 ];
        weis [ 0 ] = weis [ 4 ] = ( 322. - 13. * sqrt(70.) ) / 900;
        weis [ 1 ] = weis [ 3 ] = ( 322. + 13. * sqrt(70.) ) / 900;
        weis [ 2 ] = 128. / 225.;
    } else {
        cerr << "Gauss integration for n=" << n << " not implemented" << endl;
        exit(1);
    }
}


double find_intesection_of_segment_and_triangle(const Point *A, const Point *B, const Point *a, const Point *b, const Point *c) {
    //test bounding box
    double maxp, minp;
    for (unsigned i = 0; i < 3; i++) {
        maxp = ( * a ) [ i ];
        maxp = max(maxp, ( * b ) [ i ]);
        maxp = max(maxp, ( * c ) [ i ]);
        minp = ( * a ) [ i ];
        minp = min(minp, ( * b ) [ i ]);
        minp = min(minp, ( * c ) [ i ]);
        if ( min( ( * A ) [ i ], ( * B ) [ i ]) > maxp || max( ( * A ) [ i ], ( * B ) [ i ]) < minp ) {
            return -1.;
        }
    }

    //test intersection exists
    if ( tetraVolumeSigned(A, a, b, c) * tetraVolumeSigned(B, a, b, c) > 0 ) {
        return -1;
    }
    double v1 = tetraVolumeSigned(A, B, a, b);
    double v2 = tetraVolumeSigned(A, B, b, c);
    if ( v1 * v2 < 0 ) {
        return -1;
    }
    double v3 = tetraVolumeSigned(A, B, c, a);
    if ( v1 * v3 < 0 || v2 * v3 < 0 ) {
        return -1;
    }

    //compute intersection with facet plane
    Point dirvec = ( * B ) - ( * A );
    double length = dirvec.norm();
    dirvec.normalize();
    Point normal = Point( ( ( * b ) [ 1 ] - ( * a ) [ 1 ] ) * ( ( * c ) [ 2 ] - ( * a ) [ 2 ] ) - ( ( * c ) [ 1 ] - ( * a ) [ 1 ] ) * ( ( * b ) [ 2 ] - ( * a ) [ 2 ] ), ( ( * b ) [ 2 ] - ( * a ) [ 2 ] ) * ( ( * c ) [ 0 ] - ( * a ) [ 0 ] ) - ( ( * c ) [ 2 ] - ( * a ) [ 2 ] ) * ( ( * b ) [ 0 ] - ( * a ) [ 0 ] ), ( ( * b ) [ 0 ] - ( * a ) [ 0 ] ) * ( ( * c ) [ 1 ] - ( * a ) [ 1 ] ) - ( ( * c ) [ 0 ] - ( * a ) [ 0 ] ) * ( ( * b ) [ 1 ] - ( * a ) [ 1 ] ) );
    normal.normalize();
    double d = -( ( * a ) [ 0 ] * normal [ 0 ] + ( * a ) [ 1 ] * normal [ 1 ] + ( * a ) [ 2 ] * normal [ 2 ] );
    double t = abs( ( normal [ 0 ] * ( * A ) [ 0 ] + normal [ 1 ] * ( * A ) [ 1 ] + normal [ 2 ] * ( * A ) [ 2 ] + d ) / ( normal [ 0 ] * dirvec [ 0 ] + normal [ 1 ] * dirvec [ 1 ] + normal [ 2 ] * dirvec [ 2 ] ) );
    if ( t >= 0 && t <= length ) {
        return t;
    } else {
        return -1.;
    }
}
