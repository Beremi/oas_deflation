#ifndef _LINALG_PROFILE_H
#define _LINALG_PROFILE_H

#include "linalg.h"
#include <fstream>
#include <map>
#include <set>
#include <vector>

struct LinearProfileMatrixStats
{
    long long rows = 0;
    long long cols = 0;
    long long nnz = 0;
    std :: string patternHash = "";
    std :: string valueHash = "";
    double matrixNorm = 0.;
    double relativeMatrixDelta = -1.;
};

struct LinearProfileIterationStats
{
    int step = -1;
    int iteration = -1;
    unsigned cumulativeIteration = 0;
    std :: string systemKind;
    long long analyzeCount = 0;
    long long factorizeCount = 0;
    long long solveCount = 0;
    double analyzeSeconds = 0.;
    double factorizeSeconds = 0.;
    double solveSeconds = 0.;
    std :: string solverType;
    std :: string solverName;
    long long rows = 0;
    long long cols = 0;
    long long nnz = 0;
    std :: string patternHash;
    std :: set< std :: string >rhsKinds;
};

class LinearSolverProfiler
{
protected:
    bool enabled = false;
    bool matrixDeltaEnabled = false;
    bool headerWritten = false;
    long long eventIndex = 0;
    fs :: path eventsPath;
    fs :: path iterationsPath;
    std :: ofstream eventsFile;
    LinearProfileMatrixStats lastMatrixStats;
    std :: vector< double >lastMatrixValues;
    std :: map< std :: string, Vector >lastRhsByKind;
    std :: map< std :: string, LinearProfileIterationStats >iterationStats;

    static std :: string hashToString(uint64_t value);
    static uint64_t hashBytes(uint64_t hash, const void *data, size_t size);
    static std :: string makeIterationKey(int step, int iteration, const std :: string &systemKind);
    static double vectorNorm(const Vector &v);
    static std :: string joinKinds(const std :: set< std :: string > &values);
    void writeEventHeader();
    void recordEvent(
        const std :: string &phase,
        int step,
        int iteration,
        unsigned cumulativeIteration,
        const std :: string &systemKind,
        const std :: string &rhsKind,
        const std :: string &solverType,
        const std :: string &solverName,
        const LinearProfileMatrixStats &matrixStats,
        double rhsNorm,
        double rhsRelativeDelta,
        double solutionNorm,
        long long solverIterations,
        double solverError,
        bool success,
        double durationSeconds
    );
    void updateIterationStats(
        const std :: string &phase,
        int step,
        int iteration,
        unsigned cumulativeIteration,
        const std :: string &systemKind,
        const std :: string &rhsKind,
        const std :: string &solverType,
        const std :: string &solverName,
        const LinearProfileMatrixStats &matrixStats,
        double durationSeconds
    );
public:
    LinearSolverProfiler() {};
    ~LinearSolverProfiler();
    void configure(bool enabled, bool matrixDelta, fs :: path resultDir, const std :: string &fileBase);
    bool isEnabled() const { return enabled; };
    LinearProfileMatrixStats inspectMatrix(const CoordinateIndexedSparseMatrix &A, bool updateMatrixDeltaState);
    void recordMatrixEvent(
        const std :: string &phase,
        int step,
        int iteration,
        unsigned cumulativeIteration,
        const std :: string &systemKind,
        const std :: string &solverType,
        const std :: string &solverName,
        const CoordinateIndexedSparseMatrix &A,
        double durationSeconds,
        bool success
    );
    void recordSolveEvent(
        int step,
        int iteration,
        unsigned cumulativeIteration,
        const std :: string &systemKind,
        const std :: string &rhsKind,
        const std :: string &solverType,
        const std :: string &solverName,
        const Vector &rhs,
        const Vector &solution,
        long long solverIterations,
        double solverError,
        double durationSeconds,
        bool success
    );
    void flushIterationSummary();
};

#endif /* _LINALG_PROFILE_H */
