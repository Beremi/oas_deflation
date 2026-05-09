#include "linalg_profile.h"

#include <cmath>
#include <cstdint>
#include <cstring>
#include <sstream>

using namespace std;

namespace {
constexpr uint64_t FNV_OFFSET = 1469598103934665603ULL;
constexpr uint64_t FNV_PRIME = 1099511628211ULL;
}

//////////////////////////////////////////////////////////
LinearSolverProfiler :: ~LinearSolverProfiler() {
    flushIterationSummary();
    if ( eventsFile.is_open() ) {
        eventsFile.close();
    }
}

//////////////////////////////////////////////////////////
void LinearSolverProfiler :: configure(bool enableProfile, bool matrixDelta, fs :: path resultDir, const string &fileBase) {
    enabled = enableProfile;
    matrixDeltaEnabled = matrixDelta;
    if ( !enabled ) {
        return;
    }

    fs :: create_directories(resultDir);
    string base = fileBase.empty() ? "linear_profile" : fileBase;
    eventsPath = resultDir / ( base + "_events.tsv" );
    iterationsPath = resultDir / ( base + "_iterations.tsv" );
    eventsFile.open(eventsPath.string(), ios :: out);
    headerWritten = false;
    eventIndex = 0;
    lastMatrixValues.clear();
    lastRhsByKind.clear();
    iterationStats.clear();
    writeEventHeader();
}

//////////////////////////////////////////////////////////
string LinearSolverProfiler :: hashToString(uint64_t value) {
    ostringstream oss;
    oss << hex << value;
    return oss.str();
}

//////////////////////////////////////////////////////////
uint64_t LinearSolverProfiler :: hashBytes(uint64_t hash, const void *data, size_t size) {
    const unsigned char *bytes = static_cast< const unsigned char * >(data);
    for ( size_t i = 0; i < size; i++ ) {
        hash ^= static_cast< uint64_t >(bytes [ i ]);
        hash *= FNV_PRIME;
    }
    return hash;
}

//////////////////////////////////////////////////////////
string LinearSolverProfiler :: makeIterationKey(int step, int iteration, const string &systemKind) {
    return to_string(step) + "\t" + to_string(iteration) + "\t" + systemKind;
}

//////////////////////////////////////////////////////////
double LinearSolverProfiler :: vectorNorm(const Vector &v) {
    if ( v.size() == 0 ) {
        return 0.;
    }
    return v.norm();
}

//////////////////////////////////////////////////////////
string LinearSolverProfiler :: joinKinds(const set< string > &values) {
    string result;
    for ( const auto &value : values ) {
        if ( !result.empty() ) {
            result += ",";
        }
        result += value;
    }
    return result;
}

//////////////////////////////////////////////////////////
void LinearSolverProfiler :: writeEventHeader() {
    if ( !enabled || headerWritten || !eventsFile.good() ) {
        return;
    }
    eventsFile
        << "event_index\tphase\tstep\titeration\tcumulative_iteration\tsystem_kind\trhs_kind"
        << "\tsolver_type\tsolver_name\trows\tcols\tnnz\tpattern_hash\tvalue_hash"
        << "\tmatrix_norm\tmatrix_relative_delta\trhs_norm\trhs_relative_delta"
        << "\tsolution_norm\tsolver_iterations\tsolver_error\trequested_tolerance\tdeflation_basis_size"
        << "\tdeflation_discarded_count\tdeflation_raw_candidate_count"
        << "\tdeflation_capacity_eviction_count\tdeflation_low_anorm_discard_count"
        << "\tdeflation_last_initial_anorm\tdeflation_last_final_anorm"
        << "\tdeflation_orthogonality_max_offdiag\tdeflation_orthogonality_max_diag_error"
        << "\tdeflation_last_discard_reason\tpreconditioner_apply_seconds"
        << "\torthogonalization_seconds\tleast_squares_seconds\tmatvec_seconds"
        << "\tdeflation_seconds\tsuccess\tduration_seconds\n";
    headerWritten = true;
}

//////////////////////////////////////////////////////////
LinearProfileMatrixStats LinearSolverProfiler :: inspectMatrix(const CoordinateIndexedSparseMatrix &A, bool updateMatrixDeltaState) {
    LinearProfileMatrixStats stats;
    stats.rows = A.rows();
    stats.cols = A.cols();
    stats.nnz = A.nonZeros();

    uint64_t patternHash = FNV_OFFSET;
    patternHash = hashBytes(patternHash, &stats.rows, sizeof(stats.rows) );
    patternHash = hashBytes(patternHash, &stats.cols, sizeof(stats.cols) );
    patternHash = hashBytes(patternHash, &stats.nnz, sizeof(stats.nnz) );
    const long long outerSize = A.outerSize() + 1;
    if ( A.outerIndexPtr() != nullptr && outerSize > 0 ) {
        patternHash = hashBytes(patternHash, A.outerIndexPtr(), sizeof(CoordinateIndexedSparseMatrix :: StorageIndex) * outerSize);
    }
    if ( A.innerIndexPtr() != nullptr && stats.nnz > 0 ) {
        patternHash = hashBytes(patternHash, A.innerIndexPtr(), sizeof(CoordinateIndexedSparseMatrix :: StorageIndex) * stats.nnz);
    }
    stats.patternHash = hashToString(patternHash);

    uint64_t valueHash = FNV_OFFSET;
    double norm2 = 0.;
    vector< double >values;
    if ( A.valuePtr() != nullptr && stats.nnz > 0 ) {
        values.assign(A.valuePtr(), A.valuePtr() + stats.nnz);
        valueHash = hashBytes(valueHash, A.valuePtr(), sizeof(double) * stats.nnz);
        for ( double value : values ) {
            norm2 += value * value;
        }
    }
    stats.valueHash = hashToString(valueHash);
    stats.matrixNorm = sqrt(norm2);
    stats.relativeMatrixDelta = -1.;

    if ( updateMatrixDeltaState && matrixDeltaEnabled ) {
        if ( !lastMatrixValues.empty() && lastMatrixValues.size() == values.size() && lastMatrixStats.patternHash == stats.patternHash ) {
            double delta2 = 0.;
            for ( size_t i = 0; i < values.size(); i++ ) {
                const double diff = values [ i ] - lastMatrixValues [ i ];
                delta2 += diff * diff;
            }
            stats.relativeMatrixDelta = sqrt(delta2) / max(lastMatrixStats.matrixNorm, 1e-300);
        }
        lastMatrixValues = std :: move(values);
    }

    if ( updateMatrixDeltaState ) {
        lastMatrixStats = stats;
    }
    return stats;
}

//////////////////////////////////////////////////////////
void LinearSolverProfiler :: recordEvent(
    const string &phase,
    int step,
    int iteration,
    unsigned cumulativeIteration,
    const string &systemKind,
    const string &rhsKind,
    const string &solverType,
    const string &solverName,
    const LinearProfileMatrixStats &matrixStats,
    double rhsNorm,
    double rhsRelativeDelta,
    double solutionNorm,
    long long solverIterations,
    double solverError,
    double requestedTolerance,
    unsigned deflationBasisSize,
    unsigned deflationDiscardedCount,
    unsigned deflationRawCandidateCount,
    unsigned deflationCapacityEvictionCount,
    unsigned deflationLowANormDiscardCount,
    double deflationLastInitialANorm,
    double deflationLastFinalANorm,
    double deflationOrthogonalityMaxOffdiag,
    double deflationOrthogonalityMaxDiagError,
    const string &deflationLastDiscardReason,
    double preconditionerApplySeconds,
    double orthogonalizationSeconds,
    double leastSquaresSeconds,
    double matvecSeconds,
    double deflationSeconds,
    bool success,
    double durationSeconds
) {
    if ( !enabled || !eventsFile.good() ) {
        return;
    }
    writeEventHeader();
    eventsFile << scientific;
    eventsFile.precision(16);
    eventsFile
        << eventIndex++ << '\t'
        << phase << '\t'
        << step << '\t'
        << iteration << '\t'
        << cumulativeIteration << '\t'
        << systemKind << '\t'
        << ( rhsKind.empty() ? "-" : rhsKind ) << '\t'
        << solverType << '\t'
        << solverName << '\t'
        << matrixStats.rows << '\t'
        << matrixStats.cols << '\t'
        << matrixStats.nnz << '\t'
        << matrixStats.patternHash << '\t'
        << matrixStats.valueHash << '\t'
        << matrixStats.matrixNorm << '\t'
        << matrixStats.relativeMatrixDelta << '\t'
        << rhsNorm << '\t'
        << rhsRelativeDelta << '\t'
        << solutionNorm << '\t'
        << solverIterations << '\t'
        << solverError << '\t'
        << requestedTolerance << '\t'
        << deflationBasisSize << '\t'
        << deflationDiscardedCount << '\t'
        << deflationRawCandidateCount << '\t'
        << deflationCapacityEvictionCount << '\t'
        << deflationLowANormDiscardCount << '\t'
        << deflationLastInitialANorm << '\t'
        << deflationLastFinalANorm << '\t'
        << deflationOrthogonalityMaxOffdiag << '\t'
        << deflationOrthogonalityMaxDiagError << '\t'
        << ( deflationLastDiscardReason.empty() ? "-" : deflationLastDiscardReason ) << '\t'
        << preconditionerApplySeconds << '\t'
        << orthogonalizationSeconds << '\t'
        << leastSquaresSeconds << '\t'
        << matvecSeconds << '\t'
        << deflationSeconds << '\t'
        << ( success ? 1 : 0 ) << '\t'
        << durationSeconds << '\n';
    eventsFile.flush();

    updateIterationStats(
        phase,
        step,
        iteration,
        cumulativeIteration,
        systemKind,
        rhsKind,
        solverType,
        solverName,
        matrixStats,
        durationSeconds
    );
}

//////////////////////////////////////////////////////////
void LinearSolverProfiler :: updateIterationStats(
    const string &phase,
    int step,
    int iteration,
    unsigned cumulativeIteration,
    const string &systemKind,
    const string &rhsKind,
    const string &solverType,
    const string &solverName,
    const LinearProfileMatrixStats &matrixStats,
    double durationSeconds
) {
    if ( step < 0 || iteration < 0 ) {
        return;
    }
    const string key = makeIterationKey(step, iteration, systemKind);
    LinearProfileIterationStats &stats = iterationStats [ key ];
    stats.step = step;
    stats.iteration = iteration;
    stats.cumulativeIteration = cumulativeIteration;
    stats.systemKind = systemKind;
    stats.solverType = solverType;
    stats.solverName = solverName;
    stats.rows = matrixStats.rows;
    stats.cols = matrixStats.cols;
    stats.nnz = matrixStats.nnz;
    stats.patternHash = matrixStats.patternHash;
    if ( !rhsKind.empty() && rhsKind != "-" ) {
        stats.rhsKinds.insert(rhsKind);
    }
    if ( phase == "analyze" ) {
        stats.analyzeCount++;
        stats.analyzeSeconds += durationSeconds;
    } else if ( phase == "factorize" ) {
        stats.factorizeCount++;
        stats.factorizeSeconds += durationSeconds;
    } else if ( phase == "solve" ) {
        stats.solveCount++;
        stats.solveSeconds += durationSeconds;
    }
}

//////////////////////////////////////////////////////////
void LinearSolverProfiler :: recordMatrixEvent(
    const string &phase,
    int step,
    int iteration,
    unsigned cumulativeIteration,
    const string &systemKind,
    const string &solverType,
    const string &solverName,
    const CoordinateIndexedSparseMatrix &A,
    double durationSeconds,
    bool success
) {
    if ( !enabled ) {
        return;
    }
    const bool updateState = ( phase == "factorize" );
    LinearProfileMatrixStats matrixStats = inspectMatrix(A, updateState);
    recordEvent(
        phase,
        step,
        iteration,
        cumulativeIteration,
        systemKind,
        "-",
        solverType,
        solverName,
        matrixStats,
        -1.,
        -1.,
        -1.,
        -1,
        -1.,
        -1.,
        0,
        0,
        0,
        0,
        0,
        0.,
        0.,
        0.,
        0.,
        "-",
        0.,
        0.,
        0.,
        0.,
        0.,
        success,
        durationSeconds
    );
}

//////////////////////////////////////////////////////////
void LinearSolverProfiler :: recordSolveEvent(
    int step,
    int iteration,
    unsigned cumulativeIteration,
    const string &systemKind,
    const string &rhsKind,
    const string &solverType,
    const string &solverName,
    const Vector &rhs,
    const Vector &solution,
    long long solverIterations,
    double solverError,
    double requestedTolerance,
    unsigned deflationBasisSize,
    unsigned deflationDiscardedCount,
    unsigned deflationRawCandidateCount,
    unsigned deflationCapacityEvictionCount,
    unsigned deflationLowANormDiscardCount,
    double deflationLastInitialANorm,
    double deflationLastFinalANorm,
    double deflationOrthogonalityMaxOffdiag,
    double deflationOrthogonalityMaxDiagError,
    const string &deflationLastDiscardReason,
    double preconditionerApplySeconds,
    double orthogonalizationSeconds,
    double leastSquaresSeconds,
    double matvecSeconds,
    double deflationSeconds,
    double durationSeconds,
    bool success
) {
    if ( !enabled ) {
        return;
    }
    const double rhsNorm = vectorNorm(rhs);
    double rhsRelativeDelta = -1.;
    auto found = lastRhsByKind.find(rhsKind);
    if ( found != lastRhsByKind.end() && found->second.size() == rhs.size() ) {
        rhsRelativeDelta = ( rhs - found->second ).norm() / max(found->second.norm(), 1e-300);
    }
    lastRhsByKind [ rhsKind ] = rhs;

    recordEvent(
        "solve",
        step,
        iteration,
        cumulativeIteration,
        systemKind,
        rhsKind,
        solverType,
        solverName,
        lastMatrixStats,
        rhsNorm,
        rhsRelativeDelta,
        vectorNorm(solution),
        solverIterations,
        solverError,
        requestedTolerance,
        deflationBasisSize,
        deflationDiscardedCount,
        deflationRawCandidateCount,
        deflationCapacityEvictionCount,
        deflationLowANormDiscardCount,
        deflationLastInitialANorm,
        deflationLastFinalANorm,
        deflationOrthogonalityMaxOffdiag,
        deflationOrthogonalityMaxDiagError,
        deflationLastDiscardReason,
        preconditionerApplySeconds,
        orthogonalizationSeconds,
        leastSquaresSeconds,
        matvecSeconds,
        deflationSeconds,
        success,
        durationSeconds
    );
}

//////////////////////////////////////////////////////////
void LinearSolverProfiler :: flushIterationSummary() {
    if ( !enabled || iterationsPath.empty() ) {
        return;
    }
    ofstream output(iterationsPath.string(), ios :: out);
    if ( !output.good() ) {
        return;
    }
    output
        << "step\titeration\tcumulative_iteration\tsystem_kind\trhs_kinds\tsolver_type\tsolver_name"
        << "\trows\tcols\tnnz\tpattern_hash\tanalyze_count\tfactorize_count\tsolve_count"
        << "\tanalyze_seconds\tfactorize_seconds\tsolve_seconds\ttotal_linear_seconds\n";
    output << scientific;
    output.precision(16);
    for ( const auto &item : iterationStats ) {
        const LinearProfileIterationStats &stats = item.second;
        output
            << stats.step << '\t'
            << stats.iteration << '\t'
            << stats.cumulativeIteration << '\t'
            << stats.systemKind << '\t'
            << joinKinds(stats.rhsKinds) << '\t'
            << stats.solverType << '\t'
            << stats.solverName << '\t'
            << stats.rows << '\t'
            << stats.cols << '\t'
            << stats.nnz << '\t'
            << stats.patternHash << '\t'
            << stats.analyzeCount << '\t'
            << stats.factorizeCount << '\t'
            << stats.solveCount << '\t'
            << stats.analyzeSeconds << '\t'
            << stats.factorizeSeconds << '\t'
            << stats.solveSeconds << '\t'
            << stats.analyzeSeconds + stats.factorizeSeconds + stats.solveSeconds << '\n';
    }
}

//////////////////////////////////////////////////////////
RuntimePhaseProfiler :: ~RuntimePhaseProfiler() {
    flushSummary();
    if ( eventsFile.is_open() ) {
        eventsFile.close();
    }
}

//////////////////////////////////////////////////////////
string RuntimePhaseProfiler :: makeSummaryKey(const string &phase, const string &detail) {
    return phase + "\t" + detail;
}

//////////////////////////////////////////////////////////
void RuntimePhaseProfiler :: configure(bool enableProfile, fs :: path resultDir, const string &fileBase) {
    const bool requestedEnabled = enableProfile;
    if ( !requestedEnabled ) {
        enabled = false;
        return;
    }

    fs :: create_directories(resultDir);
    string base = fileBase.empty() ? "runtime_profile" : fileBase;
    fs :: path newEventsPath = resultDir / ( base + "_events.tsv" );
    fs :: path newSummaryPath = resultDir / ( base + "_summary.tsv" );
    if ( enabled && eventsFile.is_open() && eventsPath == newEventsPath && summaryPath == newSummaryPath ) {
        return;
    }

    enabled = true;
    if ( eventsFile.is_open() ) {
        eventsFile.close();
    }
    eventsPath = newEventsPath;
    summaryPath = newSummaryPath;
    eventsFile.open(eventsPath.string(), ios :: out);
    headerWritten = false;
    eventIndex = 0;
    summary.clear();
    writeEventHeader();
}

//////////////////////////////////////////////////////////
void RuntimePhaseProfiler :: writeEventHeader() {
    if ( !enabled || headerWritten || !eventsFile.good() ) {
        return;
    }
    eventsFile
        << "event_index\tstep\titeration\tcumulative_iteration\tsystem_kind"
        << "\tphase\tdetail\tduration_seconds\tsample_count\n";
    headerWritten = true;
}

//////////////////////////////////////////////////////////
void RuntimePhaseProfiler :: recordPhase(
    int step,
    int iteration,
    unsigned cumulativeIteration,
    const string &systemKind,
    const string &phase,
    const string &detail,
    double durationSeconds
) {
    recordPhaseSamples(step, iteration, cumulativeIteration, systemKind, phase, detail, durationSeconds, 1);
}

//////////////////////////////////////////////////////////
void RuntimePhaseProfiler :: recordPhaseSamples(
    int step,
    int iteration,
    unsigned cumulativeIteration,
    const string &systemKind,
    const string &phase,
    const string &detail,
    double durationSeconds,
    long long sampleCount
) {
    if ( !enabled ) {
        return;
    }
    if ( sampleCount <= 0 && durationSeconds <= 0. ) {
        return;
    }
    sampleCount = max(1ll, sampleCount);
    writeEventHeader();
    const string normalizedDetail = detail.empty() ? "-" : detail;
    if ( eventsFile.good() ) {
        eventsFile << scientific;
        eventsFile.precision(16);
        eventsFile
            << eventIndex++ << '\t'
            << step << '\t'
            << iteration << '\t'
            << cumulativeIteration << '\t'
            << systemKind << '\t'
            << phase << '\t'
            << normalizedDetail << '\t'
            << durationSeconds << '\t'
            << sampleCount << '\n';
        eventsFile.flush();
    }

    RuntimePhaseSummary &stats = summary [ makeSummaryKey(phase, normalizedDetail) ];
    stats.count += sampleCount;
    stats.seconds += durationSeconds;
    stats.maxSeconds = max(stats.maxSeconds, durationSeconds / sampleCount);
}

//////////////////////////////////////////////////////////
void RuntimePhaseProfiler :: flushSummary() {
    if ( !enabled || summaryPath.empty() ) {
        return;
    }
    ofstream output(summaryPath.string(), ios :: out);
    if ( !output.good() ) {
        return;
    }
    output << "phase\tdetail\tcount\ttotal_seconds\tmean_seconds\tmax_seconds\n";
    output << scientific;
    output.precision(16);
    for ( const auto &item : summary ) {
        const string &key = item.first;
        const RuntimePhaseSummary &stats = item.second;
        const size_t sep = key.find('\t');
        const string phase = ( sep == string :: npos ) ? key : key.substr(0, sep);
        const string detail = ( sep == string :: npos ) ? "-" : key.substr(sep + 1);
        output
            << phase << '\t'
            << detail << '\t'
            << stats.count << '\t'
            << stats.seconds << '\t'
            << ( stats.count > 0 ? stats.seconds / stats.count : 0. ) << '\t'
            << stats.maxSeconds << '\n';
    }
}
