#include "linalg.h"

#ifdef PARDISO_FOUND
 #include <mkl.h>     // For mkl_set_num_threads
#endif

#ifdef HYPRE_FOUND
 #include <HYPRE.h>
 #include <HYPRE_IJ_mv.h>
 #include <HYPRE_parcsr_ls.h>
 #include <HYPRE_utilities.h>
 #include <mpi.h>
#endif

#ifdef _OPENMP
 #include <omp.h>
#endif

#include <cctype>
#include <cstddef>
#include <cmath>
#include <cstdlib>
#include <limits>

using namespace std;

namespace {

#ifdef HYPRE_FOUND

double trueRelativeResidual(const CoordinateIndexedSparseMatrix &A, const Vector &x, const Vector &b) {
    if ( b.size() == 0 ) {
        return 0.;
    }
    const double bNorm = std :: max(b.norm(), 1e-300);
    return ( A * x - b ).norm() / bNorm;
}

std :: string lowerCopy(std :: string value) {
    std :: transform(value.begin(), value.end(), value.begin(), [](unsigned char c) {
        return static_cast< char >(std :: tolower(c) );
    });
    return value;
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

std :: vector< unsigned > makeElasticOldToNewPermutation(const ElasticDofMap &map, int reorderMode) {
    if ( reorderMode <= 0 || !map.isValid() || map.blockSize == 0 || map.fullRows == 0 || map.fullRows % map.blockSize != 0 ) {
        return {};
    }
    const unsigned nodeCount = map.fullRows / map.blockSize;
    std :: vector< unsigned >nodeOrder(nodeCount);
    std :: iota(nodeOrder.begin(), nodeOrder.end(), 0u);

    if ( reorderMode >= 2 && map.dimension > 0 && map.coordinates.size() >= size_t(nodeCount) * map.dimension ) {
        const unsigned dim = map.dimension;
        std :: stable_sort(nodeOrder.begin(), nodeOrder.end(), [&](unsigned a, unsigned b) {
            for ( unsigned d = 0; d < dim; d++ ) {
                const double av = map.coordinates [ size_t(a) * dim + d ];
                const double bv = map.coordinates [ size_t(b) * dim + d ];
                if ( av < bv ) {
                    return true;
                }
                if ( av > bv ) {
                    return false;
                }
            }
            return a < b;
        });
    }

    std :: vector< unsigned >oldToNew(map.fullRows);
    for ( unsigned newNode = 0; newNode < nodeCount; newNode++ ) {
        const unsigned oldNode = nodeOrder [ newNode ];
        for ( unsigned dof = 0; dof < map.blockSize; dof++ ) {
            oldToNew [ oldNode * map.blockSize + dof ] = newNode * map.blockSize + dof;
        }
    }
    return oldToNew;
}

std :: vector< unsigned > invertPermutation(const std :: vector< unsigned > &oldToNew) {
    std :: vector< unsigned >newToOld(oldToNew.size() );
    for ( unsigned oldIndex = 0; oldIndex < oldToNew.size(); oldIndex++ ) {
        newToOld [ oldToNew [ oldIndex ] ] = oldIndex;
    }
    return newToOld;
}

bool permutationIsIdentity(const std :: vector< unsigned > &oldToNew) {
    for ( unsigned i = 0; i < oldToNew.size(); i++ ) {
        if ( oldToNew [ i ] != i ) {
            return false;
        }
    }
    return true;
}

CoordinateIndexedSparseMatrix symmetricallyPermuteMatrix(const CoordinateIndexedSparseMatrix &A, const std :: vector< unsigned > &oldToNew) {
    if ( oldToNew.empty() || A.rows() != Eigen :: Index(oldToNew.size() ) || A.cols() != Eigen :: Index(oldToNew.size() ) ) {
        return A;
    }
    std :: vector< Ttripletd >triplets;
    triplets.reserve(A.nonZeros() );
    for ( int outer = 0; outer < A.outerSize(); outer++ ) {
        for ( CoordinateIndexedSparseMatrix :: InnerIterator it(A, outer); it; ++it ) {
            triplets.emplace_back(oldToNew [ it.row() ], oldToNew [ it.col() ], it.value() );
        }
    }
    CoordinateIndexedSparseMatrix permuted(A.rows(), A.cols() );
    permuted.setFromTriplets(triplets.begin(), triplets.end());
    permuted.makeCompressed();
    return permuted;
}

Vector permuteVectorToNewOrdering(const Vector &oldVector, const std :: vector< unsigned > &oldToNew) {
    if ( oldToNew.empty() || oldVector.size() != Eigen :: Index(oldToNew.size() ) ) {
        return oldVector;
    }
    Vector newVector(oldVector.size() );
    for ( Eigen :: Index oldIndex = 0; oldIndex < oldVector.size(); oldIndex++ ) {
        newVector [ oldToNew [ oldIndex ] ] = oldVector [ oldIndex ];
    }
    return newVector;
}

Vector permuteVectorToOldOrdering(const Vector &newVector, const std :: vector< unsigned > &newToOld) {
    if ( newToOld.empty() || newVector.size() != Eigen :: Index(newToOld.size() ) ) {
        return newVector;
    }
    Vector oldVector(newVector.size() );
    for ( Eigen :: Index newIndex = 0; newIndex < newVector.size(); newIndex++ ) {
        oldVector [ newToOld [ newIndex ] ] = newVector [ newIndex ];
    }
    return oldVector;
}

const char *elasticReorderName(int mode) {
    if ( mode <= 0 ) {
        return "none";
    }
    if ( mode == 1 ) {
        return "node_major";
    }
    if ( mode == 2 ) {
        return "coordinate_node_major";
    }
    return "coordinate_node_major";
}

#ifdef HYPRE_FOUND
void ensureHypreRuntimeInitialized() {
    static bool mpiInitializedHere = false;
    int mpiInitialized = 0;
    MPI_Initialized(&mpiInitialized);
    if ( !mpiInitialized ) {
        int provided = MPI_THREAD_SINGLE;
        const int mpiError = MPI_Init_thread(nullptr, nullptr, MPI_THREAD_FUNNELED, &provided);
        if ( mpiError != MPI_SUCCESS ) {
            std :: cerr << "MPI_Init_thread failed before hypre initialization" << std :: endl;
        } else {
            mpiInitializedHere = true;
        }
    }

    static bool initializedHere = false;
    if ( !HYPRE_Initialized() ) {
        HYPRE_Initialize();
        initializedHere = true;
    }
    if ( initializedHere ) {
        HYPRE_SetMemoryLocation(HYPRE_MEMORY_HOST);
        HYPRE_SetExecutionPolicy(HYPRE_EXEC_HOST);
    }

    static bool registeredFinalize = false;
    if ( mpiInitializedHere && !registeredFinalize ) {
        std :: atexit([]() {
            if ( HYPRE_Initialized() ) {
                HYPRE_Finalize();
            }
            int finalized = 0;
            MPI_Finalized(&finalized);
            if ( !finalized ) {
                MPI_Finalize();
            }
        });
        registeredFinalize = true;
    }
}

bool hypreCheck(HYPRE_Int error, const char *operation) {
    if ( error == 0 ) {
        return true;
    }
    std :: cerr << "hypre error " << error << " during " << operation << std :: endl;
    HYPRE_ClearAllErrors();
    return false;
}

int maxOpenMPThreadCount() {
#ifdef _OPENMP
    return std :: max(1, omp_get_max_threads() );
#else
    return 1;
#endif
}

bool hypreRelaxTypeUsesThreadStableSmoother(int relaxType) {
    switch ( relaxType ) {
    case 0:
    case 7:
    case 16:
    case 18:
        return true;
    default:
        return false;
    }
}

int desiredHypreThreadCount(const HypreBoomerAMGOptions &options) {
    const int maxThreads = maxOpenMPThreadCount();
    if ( options.threads > 0 ) {
        return std :: max(1, std :: min(options.threads, maxThreads) );
    }
    if ( !hypreRelaxTypeUsesThreadStableSmoother(options.relaxType) ) {
        return 1;
    }
    return maxThreads;
}

class ScopedOpenMPThreadLimit
{
public:
    explicit ScopedOpenMPThreadLimit(int threads) {
#ifdef _OPENMP
        previousThreads = omp_get_max_threads();
        if ( threads > 0 && previousThreads != threads ) {
            omp_set_num_threads(threads);
            changed = true;
        }
#else
        ( void ) threads;
#endif
    }

    ~ScopedOpenMPThreadLimit() {
#ifdef _OPENMP
        if ( changed ) {
            omp_set_num_threads(previousThreads);
        }
#endif
    }

private:
    int previousThreads = 1;
    bool changed = false;
};

bool setHypreMatrixValuesFromRowMajor(
    HYPRE_IJMatrix matrix,
    const Eigen :: SparseMatrix< double, Eigen :: RowMajor, ptrdiff_t > &rowMajor
) {
    const ptrdiff_t rows = rowMajor.rows();
    if ( rows <= 0 || rows > static_cast< ptrdiff_t >( std :: numeric_limits< HYPRE_Int > :: max() ) ) {
        return false;
    }

    std :: vector< HYPRE_Int >rowSizes(rows, 0);
    std :: vector< HYPRE_BigInt >rowIds(rows, 0);
    std :: vector< HYPRE_BigInt >rowColumns(rowMajor.nonZeros(), 0);
    for ( ptrdiff_t row = 0; row < rows; row++ ) {
        rowSizes [ row ] = static_cast< HYPRE_Int >( rowMajor.outerIndexPtr() [ row + 1 ] - rowMajor.outerIndexPtr() [ row ] );
        rowIds [ row ] = static_cast< HYPRE_BigInt >( row );
    }
    for ( ptrdiff_t i = 0; i < rowMajor.nonZeros(); i++ ) {
        rowColumns [ i ] = static_cast< HYPRE_BigInt >( rowMajor.innerIndexPtr() [ i ] );
    }

    return hypreCheck(
        HYPRE_IJMatrixSetValues(
            matrix,
            static_cast< HYPRE_Int >( rows ),
            rowSizes.data(),
            rowIds.data(),
            rowColumns.data(),
            rowMajor.valuePtr()
        ),
        "HYPRE_IJMatrixSetValues"
    );
}

const HYPRE_BigInt *hypreIndexPointer(const std :: vector< HYPRE_BigInt > &indices) {
    return indices.empty() ? nullptr : indices.data();
}

bool makeHypreVector(
    const std :: vector< HYPRE_BigInt > &indices,
    const std :: vector< double > &values,
    HYPRE_IJVector &ijVector,
    HYPRE_ParVector &parVector
) {
    if ( indices.empty() || indices.size() != values.size() ) {
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
    if ( !hypreCheck(HYPRE_IJVectorSetValues(ijVector, static_cast< HYPRE_Int >( indices.size() ), hypreIndexPointer(indices), values.data() ), "HYPRE_IJVectorSetValues") ) {
        return false;
    }
    if ( !hypreCheck(HYPRE_IJVectorAssemble(ijVector), "HYPRE_IJVectorAssemble") ) {
        return false;
    }
    return hypreCheck(HYPRE_IJVectorGetObject(ijVector, reinterpret_cast< void ** >( &parVector ) ), "HYPRE_IJVectorGetObject");
}
#endif

#endif

} // namespace

//////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////
// CONJUGATE GRADIENT SOLVER
//////////////////////////////////////////////////////////

ConjGradSolver :: ConjGradSolver() {
    name = "EigenConj";
    relMaxIT = 1.;
    precision = 1e-12;
}

//////////////////////////////////////////////////////////
ConjGradSolver :: ~ConjGradSolver() {}

//////////////////////////////////////////////////////////
void ConjGradSolver :: setPrecisionAndRelMaxIters(double p, double rmi) {
    relMaxIT = rmi;
    precision = p;
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
        result = size_t( cgK.iterations() ) < maxIT;
        if ( !result ) {
            cerr << "Eigen Conjugate Gradients performed " << cgK.iterations() << " iterations and reached error " << cgK.error() << ", required precision is " << precision << endl;
            exit(1);
        }
        initialGuess = x;
    } else {
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

#ifdef HYPRE_FOUND
//////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////
// DEFLATED FLEXIBLE GMRES SOLVER
//////////////////////////////////////////////////////////

struct DeflatedFGMRESSolver :: Impl
{
    Eigen :: Index rows = 0;
    CoordinateIndexedSparseMatrix reducedMatrix;
    CoordinateIndexedSparseMatrix activeMatrix;
    std :: vector< unsigned >activeOldToNew;
    std :: vector< unsigned >activeNewToOld;
    HYPRE_IJMatrix hypreIjMatrix = nullptr;
    HYPRE_ParCSRMatrix hypreParMatrix = nullptr;
    HYPRE_Solver hyprePrecond = nullptr;
    HYPRE_IJVector hypreSetupIjB = nullptr;
    HYPRE_IJVector hypreSetupIjX = nullptr;
    HYPRE_ParVector hypreSetupParB = nullptr;
    HYPRE_ParVector hypreSetupParX = nullptr;
    HYPRE_IJVector hypreApplyIjB = nullptr;
    HYPRE_IJVector hypreApplyIjX = nullptr;
    HYPRE_ParVector hypreApplyParB = nullptr;
    HYPRE_ParVector hypreApplyParX = nullptr;
    std :: vector< HYPRE_BigInt >hypreIndices;
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
        if ( hypreApplyIjB ) {
            HYPRE_IJVectorDestroy(hypreApplyIjB);
            hypreApplyIjB = nullptr;
            hypreApplyParB = nullptr;
        }
        if ( hypreApplyIjX ) {
            HYPRE_IJVectorDestroy(hypreApplyIjX);
            hypreApplyIjX = nullptr;
            hypreApplyParX = nullptr;
        }
        if ( hypreIjMatrix ) {
            HYPRE_IJMatrixDestroy(hypreIjMatrix);
            hypreIjMatrix = nullptr;
            hypreParMatrix = nullptr;
        }
        hypreIndices.clear();
    }

    void clearPreconditioner() {
        clearHyprePreconditioner();
        activeOldToNew.clear();
        activeNewToOld.clear();
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
}

//////////////////////////////////////////////////////////
DeflatedFGMRESSolver :: ~DeflatedFGMRESSolver() {}

//////////////////////////////////////////////////////////
void DeflatedFGMRESSolver :: setOptions(const DeflatedFGMRESOptions &opts) {
    options = opts;
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
    Vector active = ( options.elasticFullLift && elasticMap.isValid() ) ? liftVectorToElasticSpace(v, elasticMap) : v;
    if ( !impl->activeOldToNew.empty() ) {
        active = permuteVectorToNewOrdering(active, impl->activeOldToNew);
    }
    return active;
}

//////////////////////////////////////////////////////////
Vector DeflatedFGMRESSolver :: activeMatvec(const Vector &v) const {
    return impl->activeMatrix * v;
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

    Vector au = activeMatvec(u);
    const double initialANorm = u.dot(au);
    impl->lastCandidateInitialANorm = initialANorm;
    if ( !std :: isfinite(initialANorm) ) {
        impl->discardedBasisVectors++;
        impl->nonfiniteDiscardCount++;
        impl->lastDiscardReason = "nonfinite_a_norm";
        return false;
    }

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
    if ( lowerCopy(options.preconditioner) == "none" || lowerCopy(options.preconditioner) == "identity" ) {
        return v;
    }
    if ( impl->hyprePrecond && impl->hypreParMatrix ) {
        if ( v.size() != Eigen :: Index(impl->hypreIndices.size() )
             || !impl->hypreApplyIjB || !impl->hypreApplyIjX
             || !impl->hypreApplyParB || !impl->hypreApplyParX ) {
            return v;
        }

        const HYPRE_BigInt *indexPointer = hypreIndexPointer(impl->hypreIndices);
        const HYPRE_Int vectorSize = static_cast< HYPRE_Int >( impl->hypreIndices.size() );
        if ( !hypreCheck(HYPRE_IJVectorSetValues(impl->hypreApplyIjB, vectorSize, indexPointer, v.data() ), "HYPRE_IJVectorSetValues") ) {
            return v;
        }
        if ( !hypreCheck(HYPRE_IJVectorSetConstantValues(impl->hypreApplyIjX, 0.), "HYPRE_IJVectorSetConstantValues") ) {
            return v;
        }

        HYPRE_Int solveError = 0;
        {
            ScopedOpenMPThreadLimit hypreThreadLimit(desiredHypreThreadCount(hypreOptions) );
            solveError = HYPRE_BoomerAMGSolve(impl->hyprePrecond, impl->hypreParMatrix, impl->hypreApplyParB, impl->hypreApplyParX);
        }
        if ( solveError != 0 ) {
            HYPRE_ClearAllErrors();
            if ( options.verbose ) {
                std :: cerr << "DeflatedFGMRES hypre preconditioner apply failed with error " << solveError << "; using identity fallback" << std :: endl;
            }
            return v;
        }
        Vector answer(v.size() );
        if ( !hypreCheck(HYPRE_IJVectorGetValues(impl->hypreApplyIjX, vectorSize, indexPointer, answer.data() ), "HYPRE_IJVectorGetValues") ) {
            return v;
        }
        return answer;
    }
    return v;
}

//////////////////////////////////////////////////////////
bool DeflatedFGMRESSolver :: factorize(const CoordinateIndexedSparseMatrix &A) {
    if ( A.rows() <= 0 ) {
        return true;
    }

    impl->clearPreconditioner();
    const bool useElasticLift = options.elasticFullLift && elasticMap.isValid();
    const std :: string preconditionerMode = lowerCopy(options.preconditioner);
    const bool useHypreApply = preconditionerMode == "hypre" || preconditionerMode == "boomeramg"
        || preconditionerMode == "hypre_boomeramg" || preconditionerMode == "hypre-boomeramg";
    const int activeElasticReorderMode = options.elasticReorder > 0 ? options.elasticReorder : hypreOptions.elasticReorder;

    impl->reducedMatrix = A;
    impl->activeMatrix = useElasticLift ? liftMatrixToElasticSpace(A, elasticMap) : A;
    if ( activeElasticReorderMode > 0 && useElasticLift && elasticMap.blockSize > 1 ) {
        impl->activeOldToNew = makeElasticOldToNewPermutation(elasticMap, activeElasticReorderMode);
        if ( impl->activeOldToNew.size() == size_t( impl->activeMatrix.rows() ) && !permutationIsIdentity(impl->activeOldToNew) ) {
            impl->activeNewToOld = invertPermutation(impl->activeOldToNew);
            impl->activeMatrix = symmetricallyPermuteMatrix(impl->activeMatrix, impl->activeOldToNew);
        } else {
            impl->activeOldToNew.clear();
            impl->activeNewToOld.clear();
        }
    }
    impl->rows = impl->activeMatrix.rows();

    if ( hypreOptions.checkMatrix ) {
        Vector probe(impl->activeMatrix.cols() );
        for ( Eigen :: Index i = 0; i < probe.size(); i++ ) {
            probe [ i ] = sin(0.001 * static_cast< double >( i + 1 ) );
        }
        const Vector product = activeMatvec(probe);
        const Vector reference = impl->activeMatrix * probe;
        const double relativeDifference = ( product - reference ).norm() / std :: max(reference.norm(), 1e-300);
        std :: cout << "DeflatedFGMRES: active matvec check relative_difference=" << relativeDifference << std :: endl;
        if ( !std :: isfinite(relativeDifference) || relativeDifference > 1e-10 ) {
            std :: cerr << "DeflatedFGMRES Error: active matvec does not match Eigen sparse matvec" << std :: endl;
            impl->clearPreconditioner();
            return false;
        }
    }

    if ( useHypreApply ) {
        ensureHypreRuntimeInitialized();
        if ( impl->rows <= 0 || impl->rows > std :: numeric_limits< HYPRE_Int > :: max() ) {
            std :: cerr << "DeflatedFGMRES hypre preconditioner supports positive row counts fitting HYPRE_Int; rows=" << impl->rows << std :: endl;
            impl->clearPreconditioner();
            return false;
        }

        const HYPRE_BigInt firstRow = 0;
        const HYPRE_BigInt lastRow = impl->rows - 1;
        impl->hypreIndices.resize(impl->rows);
        for ( Eigen :: Index row = 0; row < impl->rows; row++ ) {
            impl->hypreIndices [ row ] = row;
        }

        Eigen :: SparseMatrix< double, Eigen :: RowMajor, ptrdiff_t >rowMajor = impl->activeMatrix;
        rowMajor.makeCompressed();
        if ( !hypreCheck(HYPRE_IJMatrixCreate(MPI_COMM_SELF, firstRow, lastRow, firstRow, lastRow, &impl->hypreIjMatrix), "HYPRE_IJMatrixCreate") ) {
            impl->clearPreconditioner();
            return false;
        }
        if ( !hypreCheck(HYPRE_IJMatrixSetObjectType(impl->hypreIjMatrix, HYPRE_PARCSR), "HYPRE_IJMatrixSetObjectType") ) {
            impl->clearPreconditioner();
            return false;
        }
        std :: vector< HYPRE_Int >rowSizes(impl->rows, 0);
        for ( Eigen :: Index row = 0; row < impl->rows; row++ ) {
            rowSizes [ row ] = static_cast< HYPRE_Int >( rowMajor.outerIndexPtr() [ row + 1 ] - rowMajor.outerIndexPtr() [ row ] );
        }
        HYPRE_IJMatrixSetRowSizes(impl->hypreIjMatrix, rowSizes.data() );
        HYPRE_IJMatrixSetOMPFlag(impl->hypreIjMatrix, 1);
        if ( !hypreCheck(HYPRE_IJMatrixInitialize(impl->hypreIjMatrix), "HYPRE_IJMatrixInitialize") ) {
            impl->clearPreconditioner();
            return false;
        }
        if ( !setHypreMatrixValuesFromRowMajor(impl->hypreIjMatrix, rowMajor) ) {
            impl->clearPreconditioner();
            return false;
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
        if ( hypreOptions.numSweeps > 0 ) {
            HYPRE_BoomerAMGSetNumSweeps(impl->hyprePrecond, hypreOptions.numSweeps);
        }
        HYPRE_BoomerAMGSetPMaxElmts(impl->hyprePrecond, hypreOptions.pMaxElmts);
        HYPRE_BoomerAMGSetAggNumLevels(impl->hyprePrecond, hypreOptions.aggNumLevels);
        if ( hypreOptions.nodalDiag != 0 ) {
            HYPRE_BoomerAMGSetNodalDiag(impl->hyprePrecond, hypreOptions.nodalDiag);
        }
        if ( hypreOptions.chebyOrder > 0 ) {
            HYPRE_BoomerAMGSetChebyOrder(impl->hyprePrecond, hypreOptions.chebyOrder);
        }
        if ( hypreOptions.chebyFraction > 0. ) {
            HYPRE_BoomerAMGSetChebyFraction(impl->hyprePrecond, hypreOptions.chebyFraction);
        }
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
            if ( hypreOptions.useInterpVectors && options.verbose ) {
                std :: cout << "DeflatedFGMRES hypre: interpolation vectors requested but skipped in this minimal port." << std :: endl;
            }
        }

        std :: vector< double >setupRhs(impl->rows, 1.);
        std :: vector< double >setupGuess(impl->rows, 0.);
        if ( !makeHypreVector(impl->hypreIndices, setupRhs, impl->hypreSetupIjB, impl->hypreSetupParB)
             || !makeHypreVector(impl->hypreIndices, setupGuess, impl->hypreSetupIjX, impl->hypreSetupParX)
             || !makeHypreVector(impl->hypreIndices, setupGuess, impl->hypreApplyIjB, impl->hypreApplyParB)
             || !makeHypreVector(impl->hypreIndices, setupGuess, impl->hypreApplyIjX, impl->hypreApplyParX) ) {
            impl->clearPreconditioner();
            return false;
        }
        {
            ScopedOpenMPThreadLimit hypreThreadLimit(desiredHypreThreadCount(hypreOptions) );
            if ( !hypreCheck(HYPRE_BoomerAMGSetup(impl->hyprePrecond, impl->hypreParMatrix, impl->hypreSetupParB, impl->hypreSetupParX), "HYPRE_BoomerAMGSetup") ) {
                impl->clearPreconditioner();
                return false;
            }
        }
        impl->activePreconditioner = "hypre_boomeramg_apply";
    } else {
        impl->activePreconditioner = "identity";
    }

    impl->matrixReady = true;
    if ( options.reorthogonalizeOnMatrixChange ) {
        rebuildDeflationBasis();
    }
    if ( options.verbose || useHypreApply ) {
        std :: cout << "DeflatedFGMRES setup complete: rows=" << impl->rows
                    << ", nnz=" << impl->activeMatrix.nonZeros()
                    << ", reduced_rows=" << A.rows()
                    << ", preconditioner=" << impl->activePreconditioner
                    << ", tolerance=" << options.tolerance
                    << ", true_tolerance=" << options.trueTolerance
                    << ", restart=" << options.restart
                    << ", max_iterations=" << options.maxIterations
                    << ", basis_size=" << impl->basisU.size()
                    << ", elastic_lift=" << ( useElasticLift ? 1 : 0 )
                    << ", elastic_reorder=" << elasticReorderName(activeElasticReorderMode)
                    << ", hypre_threads=" << ( useHypreApply ? desiredHypreThreadCount(hypreOptions) : 0 )
                    << std :: endl;
    }
    return true;
}

//////////////////////////////////////////////////////////
bool DeflatedFGMRESSolver :: solve(Vector &x, const Vector &b) {
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

    const bool useElasticLift = options.elasticFullLift && elasticMap.isValid();
    Vector activeB = useElasticLift ? liftVectorToElasticSpace(b, elasticMap) : b;
    if ( !impl->activeOldToNew.empty() ) {
        activeB = permuteVectorToNewOrdering(activeB, impl->activeOldToNew);
    }

    const Eigen :: Index n = activeB.size();
    const unsigned restart = std :: max(1u, std :: min(options.restart, options.maxIterations) );
    Vector activeX = Vector :: Zero(n);
    const double bNorm = std :: max(activeB.norm(), 1e-300);

    if ( !impl->basisU.empty() ) {
        for ( size_t i = 0; i < impl->basisU.size(); ++i ) {
            const double coefficient = impl->basisU [ i ].dot(activeB);
            activeX.noalias() += coefficient * impl->basisU [ i ];
        }
    }

    Vector residual = activeB - activeMatvec(activeX);
    double trueResidual = residual.norm() / bNorm;
    double gmresResidual = trueResidual;
    unsigned totalIterations = 0;

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
            z = projectDeflation(z);
            Vector w = activeMatvec(z);

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

            H(j + 1, j) = w.norm();
            Z.col(j) = z;
            usedInner = j + 1;
            if ( H(j + 1, j) > 1e-14 ) {
                V.col(j + 1) = w / H(j + 1, j);
            }

            Eigen :: MatrixXd subH = H.block(0, 0, j + 2, j + 1);
            Eigen :: VectorXd subg = g.head(j + 2);
            y = subH.colPivHouseholderQr().solve(subg);
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
                y = subH.colPivHouseholderQr().solve(subg);
            }
            activeX.noalias() += Z.leftCols(usedInner) * y.head(usedInner);
        }

        residual = activeB - activeMatvec(activeX);
        trueResidual = residual.norm() / bNorm;
        if ( trueResidual <= options.trueTolerance || happyBreakdown || usedInner == 0 || totalIterations >= options.maxIterations ) {
            break;
        }
    }

    if ( !impl->activeNewToOld.empty() ) {
        activeX = permuteVectorToOldOrdering(activeX, impl->activeNewToOld);
    }
    x = useElasticLift ? restrictVectorFromElasticSpace(activeX, elasticMap) : activeX;

    lastIterations = totalIterations;
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
