#include <Eigen/Core>

#include <HYPRE.h>
#include <HYPRE_IJ_mv.h>
#include <HYPRE_parcsr_ls.h>
#include <mpi.h>

#ifdef _OPENMP
#include <omp.h>
#endif

#include <algorithm>
#include <array>
#include <cerrno>
#include <chrono>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <deque>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>
#include <map>
#include <numeric>
#include <sstream>
#include <stdexcept>
#include <string>
#include <unordered_map>
#include <vector>

namespace {

struct Options {
    std::string replayDir;
    std::string representation = "reduced-original";
    std::string nullspace = "off";
    std::string solver = "fgmres";
    int threads = 1;
    int nodal = 4;
    int interpType = 6;
    double strongThreshold = 0.25;
    int relaxType = 16;
    int numSweeps = 3;
    int pMaxElmts = 4;
    int aggNumLevels = 0;
    int boomerMaxIterations = 1;
    int chebyOrder = 3;
    double chebyFraction = -1.0;
    double nonGalerkinTol = -1.0;
    int interpVecVariant = 2;
    int interpVecQMax = 0;
    double interpVecAbsQTrunc = 0.0;
    int nodalDiag = 0;
    int keepSameSign = 0;
    double nullspaceCoordinateScale = std::numeric_limits<double>::quiet_NaN();
    double tolerance = 1.0e-1;
    int maxIterations = 500;
    int applyWarmup = 1;
    int applyRepeat = 5;
    bool csv = false;
};

struct CsrMatrix {
    int rows = 0;
    int cols = 0;
    std::vector<long long> ptr;
    std::vector<int> col;
    std::vector<double> val;

    long long nnz() const {
        return static_cast<long long>(val.size());
    }
};

struct MetadataRow {
    int reducedRow = -1;
    int fullDof = -1;
    int nodeId = -1;
    int relativeDof = -1;
    int physicalField = -1;
    double x = 0.0;
    double y = 0.0;
    double z = 0.0;
    int fullElasticRow = -1;
};

struct ReplayData {
    CsrMatrix reducedMatrix;
    Eigen::VectorXd rhs;
    Eigen::VectorXd referenceSolution;
    std::vector<MetadataRow> metadata;
    std::vector<int> reducedToFull;
    std::vector<int> fullToReduced;
    std::vector<double> nearNullspace;
    int nearNullspaceColumns = 0;
    int fullRows = 0;
    int blockSize = 0;
    int dimension = 0;
    std::vector<std::array<double, 3>> nodeCoordinates;
    std::vector<unsigned char> nodeHasCoordinate;
};

struct VariantData {
    CsrMatrix matrix;
    Eigen::VectorXd rhs;
    Eigen::VectorXd referenceInActiveOrder;
    std::vector<int> activeToReduced;
    std::vector<int> dofFunc;
    std::vector<double> nearNullspace;
    int nearNullspaceColumns = 0;
    int blockSize = 1;
    int dimension = 0;
    bool systemAmgCompatible = false;
    bool interpolationVectorsCompatible = false;
    std::string dofMode = "none";
};

struct HypreObjects {
    HYPRE_IJMatrix ijMatrix = nullptr;
    HYPRE_ParCSRMatrix parMatrix = nullptr;
    HYPRE_IJVector ijB = nullptr;
    HYPRE_IJVector ijX = nullptr;
    HYPRE_ParVector parB = nullptr;
    HYPRE_ParVector parX = nullptr;
    HYPRE_Solver precond = nullptr;
    HYPRE_Solver solver = nullptr;
    std::string solverKind = "fgmres";
    std::vector<HYPRE_Int> dofFuncStorage;
    std::vector<HYPRE_IJVector> interpIjVectors;
    std::vector<HYPRE_ParVector> interpParVectors;

    ~HypreObjects() {
        if (solver) {
            if (solverKind == "pcg") {
                HYPRE_ParCSRPCGDestroy(solver);
            } else {
                HYPRE_ParCSRFlexGMRESDestroy(solver);
            }
        }
        if (precond) HYPRE_BoomerAMGDestroy(precond);
        for (HYPRE_IJVector vector : interpIjVectors) {
            if (vector) HYPRE_IJVectorDestroy(vector);
        }
        if (ijB) HYPRE_IJVectorDestroy(ijB);
        if (ijX) HYPRE_IJVectorDestroy(ijX);
        if (ijMatrix) HYPRE_IJMatrixDestroy(ijMatrix);
    }
};

struct Timings {
    double matrixSeconds = 0.0;
    double vectorSeconds = 0.0;
    double preconditionerSetupSeconds = 0.0;
    double solverSetupSeconds = 0.0;
    double applyMeanSeconds = 0.0;
    double solveSeconds = 0.0;
};

struct ResultRow {
    bool success = true;
    std::string message = "ok";
    Options options;
    int ompMaxThreads = 1;
    std::string ompNumThreads;
    std::string ompDynamic;
    std::string ompProcBind;
    std::string ompPlaces;
    std::string dofMode;
    int rows = 0;
    int cols = 0;
    long long nnz = 0;
    int fullRows = 0;
    int blockSize = 0;
    int nearNullspaceColumns = 0;
    int usedInterpVectors = 0;
    double symmetryProbeRelative = 0.0;
    Timings timings;
    int iterations = -1;
    double hypreResidual = -1.0;
    double trueResidual = -1.0;
    double solutionError = -1.0;
};

double secondsSince(const std::chrono::steady_clock::time_point &start) {
    return std::chrono::duration<double>(std::chrono::steady_clock::now() - start).count();
}

std::string envOrEmpty(const char *name) {
    const char *value = std::getenv(name);
    return value ? value : "";
}

std::string pathJoin(const std::string &dir, const std::string &name) {
    return dir.empty() || dir.back() == '/' ? dir + name : dir + "/" + name;
}

std::vector<std::string> splitTab(const std::string &line) {
    std::vector<std::string> values;
    std::stringstream stream(line);
    std::string item;
    while (std::getline(stream, item, '\t')) {
        values.push_back(item);
    }
    return values;
}

void printUsage(const char *argv0) {
    std::cerr
        << "Usage: " << argv0 << " --replay-dir DIR [options]\n"
        << "  --representation reduced-original|full-node-major|full-coordinate-node-major|full-rcm-node-major|full-component-major|compressed-free-node-major\n"
        << "  --nullspace off|rbm6|rbm6-normalized|rbm6-physical|rbm6-scale|rbm6-rot-neg|rbm6-rot-zero|transrot6\n"
        << "  --solver fgmres|pcg\n"
        << "  --threads N\n"
        << "  --nodal N --interp-type N --strong-threshold X --relax-type N --num-sweeps N\n"
        << "  --non-galerkin-tol VALUE, negative disables\n"
        << "  --interp-vec-qmax N --interp-vec-abs-qtrunc X\n"
        << "  --nodal-diag N --keep-same-sign N\n"
        << "  --csv\n";
}

Options parseOptions(int argc, char **argv) {
    Options options;
    for (int i = 1; i < argc; i++) {
        const std::string arg = argv[i];
        auto value = [&](const char *name) -> std::string {
            if (i + 1 >= argc) throw std::runtime_error(std::string("missing value for ") + name);
            return argv[++i];
        };
        if (arg == "--replay-dir") options.replayDir = value("--replay-dir");
        else if (arg == "--representation") options.representation = value("--representation");
        else if (arg == "--nullspace") options.nullspace = value("--nullspace");
        else if (arg == "--solver") options.solver = value("--solver");
        else if (arg == "--threads") options.threads = std::stoi(value("--threads"));
        else if (arg == "--nodal") options.nodal = std::stoi(value("--nodal"));
        else if (arg == "--interp-type") options.interpType = std::stoi(value("--interp-type"));
        else if (arg == "--strong-threshold") options.strongThreshold = std::stod(value("--strong-threshold"));
        else if (arg == "--relax-type") options.relaxType = std::stoi(value("--relax-type"));
        else if (arg == "--num-sweeps") options.numSweeps = std::stoi(value("--num-sweeps"));
        else if (arg == "--p-max") options.pMaxElmts = std::stoi(value("--p-max"));
        else if (arg == "--agg-levels") options.aggNumLevels = std::stoi(value("--agg-levels"));
        else if (arg == "--boomer-max-iterations") options.boomerMaxIterations = std::stoi(value("--boomer-max-iterations"));
        else if (arg == "--cheby-order") options.chebyOrder = std::stoi(value("--cheby-order"));
        else if (arg == "--cheby-fraction") options.chebyFraction = std::stod(value("--cheby-fraction"));
        else if (arg == "--non-galerkin-tol") options.nonGalerkinTol = std::stod(value("--non-galerkin-tol"));
        else if (arg == "--interp-vec-variant") options.interpVecVariant = std::stoi(value("--interp-vec-variant"));
        else if (arg == "--interp-vec-qmax") options.interpVecQMax = std::stoi(value("--interp-vec-qmax"));
        else if (arg == "--interp-vec-abs-qtrunc") options.interpVecAbsQTrunc = std::stod(value("--interp-vec-abs-qtrunc"));
        else if (arg == "--nodal-diag") options.nodalDiag = std::stoi(value("--nodal-diag"));
        else if (arg == "--keep-same-sign") options.keepSameSign = std::stoi(value("--keep-same-sign"));
        else if (arg == "--nullspace-coordinate-scale") options.nullspaceCoordinateScale = std::stod(value("--nullspace-coordinate-scale"));
        else if (arg == "--tolerance") options.tolerance = std::stod(value("--tolerance"));
        else if (arg == "--max-iterations") options.maxIterations = std::stoi(value("--max-iterations"));
        else if (arg == "--apply-warmup") options.applyWarmup = std::stoi(value("--apply-warmup"));
        else if (arg == "--apply-repeat") options.applyRepeat = std::stoi(value("--apply-repeat"));
        else if (arg == "--csv") options.csv = true;
        else if (arg == "--help" || arg == "-h") {
            printUsage(argv[0]);
            std::exit(0);
        } else {
            throw std::runtime_error("unknown argument: " + arg);
        }
    }
    if (options.replayDir.empty()) throw std::runtime_error("--replay-dir is required");
    options.threads = std::max(1, options.threads);
    options.applyRepeat = std::max(1, options.applyRepeat);
    return options;
}

void hypreCheck(HYPRE_Int error, const char *operation) {
    if (error != 0) {
        std::ostringstream out;
        out << "hypre error " << error << " during " << operation;
        HYPRE_ClearAllErrors();
        throw std::runtime_error(out.str());
    }
}

void initializeHypre(int threads) {
    int mpiInitialized = 0;
    MPI_Initialized(&mpiInitialized);
    if (!mpiInitialized) {
        int provided = 0;
        MPI_Init_thread(nullptr, nullptr, MPI_THREAD_FUNNELED, &provided);
    }
    if (!HYPRE_Initialized()) {
        HYPRE_Initialize();
        HYPRE_SetMemoryLocation(HYPRE_MEMORY_HOST);
        HYPRE_SetExecutionPolicy(HYPRE_EXEC_HOST);
    }
#ifdef _OPENMP
    omp_set_num_threads(std::max(1, threads));
#else
    (void)threads;
#endif
}

CsrMatrix loadMatrixMarket(const std::string &path) {
    auto open = [](const std::string &filePath) -> FILE * {
        FILE *file = std::fopen(filePath.c_str(), "r");
        if (!file) throw std::runtime_error("cannot open " + filePath);
        return file;
    };

    int rows = 0;
    int cols = 0;
    long long inputNnz = 0;
    char line[4096];
    {
        FILE *file = open(path);
        while (std::fgets(line, sizeof(line), file)) {
            if (line[0] == '%') continue;
            if (std::sscanf(line, "%d %d %lld", &rows, &cols, &inputNnz) == 3) break;
        }
        if (rows <= 0 || cols <= 0 || inputNnz <= 0) {
            std::fclose(file);
            throw std::runtime_error("bad MatrixMarket header in " + path);
        }
        CsrMatrix matrix;
        matrix.rows = rows;
        matrix.cols = cols;
        matrix.ptr.assign(static_cast<size_t>(rows) + 1, 0);
        while (std::fgets(line, sizeof(line), file)) {
            if (line[0] == '%' || line[0] == '\n') continue;
            char *cursor = line;
            const long row = std::strtol(cursor, &cursor, 10);
            std::strtol(cursor, &cursor, 10);
            std::strtod(cursor, &cursor);
            if (row <= 0 || row > rows) {
                std::fclose(file);
                throw std::runtime_error("bad MatrixMarket row in " + path);
            }
            matrix.ptr[static_cast<size_t>(row)]++;
        }
        std::fclose(file);

        for (int row = 0; row < rows; row++) {
            matrix.ptr[static_cast<size_t>(row + 1)] += matrix.ptr[static_cast<size_t>(row)];
        }
        matrix.col.assign(static_cast<size_t>(inputNnz), 0);
        matrix.val.assign(static_cast<size_t>(inputNnz), 0.0);
        std::vector<long long> cursor = matrix.ptr;

        file = open(path);
        while (std::fgets(line, sizeof(line), file)) {
            if (line[0] == '%') continue;
            int dummyRows = 0;
            int dummyCols = 0;
            long long dummyNnz = 0;
            if (std::sscanf(line, "%d %d %lld", &dummyRows, &dummyCols, &dummyNnz) == 3) break;
        }
        while (std::fgets(line, sizeof(line), file)) {
            if (line[0] == '%' || line[0] == '\n') continue;
            char *entry = line;
            errno = 0;
            const long row = std::strtol(entry, &entry, 10);
            const long col = std::strtol(entry, &entry, 10);
            const double value = std::strtod(entry, &entry);
            if (errno != 0 || row <= 0 || col <= 0 || row > rows || col > cols) {
                std::fclose(file);
                throw std::runtime_error("bad MatrixMarket entry in " + path);
            }
            const long long index = cursor[static_cast<size_t>(row - 1)]++;
            matrix.col[static_cast<size_t>(index)] = static_cast<int>(col - 1);
            matrix.val[static_cast<size_t>(index)] = value;
        }
        std::fclose(file);
        return matrix;
    }
}

Eigen::VectorXd loadVector(const std::string &path, int size) {
    std::ifstream in(path);
    if (!in) throw std::runtime_error("cannot open " + path);
    Eigen::VectorXd vector = Eigen::VectorXd::Zero(size);
    std::string line;
    std::getline(in, line);
    while (std::getline(in, line)) {
        if (line.empty()) continue;
        const std::vector<std::string> values = splitTab(line);
        if (values.size() < 2) continue;
        const int index = std::stoi(values[0]);
        if (index < 0 || index >= size) throw std::runtime_error("bad vector index in " + path);
        vector[index] = std::stod(values[1]);
    }
    return vector;
}

std::map<std::string, std::string> loadSummary(const std::string &path) {
    std::map<std::string, std::string> summary;
    std::ifstream in(path);
    if (!in) return summary;
    std::string line;
    std::getline(in, line);
    while (std::getline(in, line)) {
        const std::vector<std::string> values = splitTab(line);
        if (values.size() >= 2) summary[values[0]] = values[1];
    }
    return summary;
}

std::vector<MetadataRow> loadMetadata(const std::string &path, int rows) {
    std::ifstream in(path);
    if (!in) throw std::runtime_error("cannot open " + path);
    std::vector<MetadataRow> metadata(static_cast<size_t>(rows));
    std::string line;
    std::getline(in, line);
    while (std::getline(in, line)) {
        if (line.empty()) continue;
        const std::vector<std::string> values = splitTab(line);
        if (values.size() < 9) throw std::runtime_error("bad metadata row in " + path);
        MetadataRow row;
        row.reducedRow = std::stoi(values[0]);
        row.fullDof = std::stoi(values[1]);
        row.nodeId = std::stoi(values[2]);
        row.relativeDof = std::stoi(values[3]);
        row.physicalField = std::stoi(values[4]);
        row.x = std::stod(values[5]);
        row.y = std::stod(values[6]);
        row.z = std::stod(values[7]);
        row.fullElasticRow = std::stoi(values[8]);
        if (row.reducedRow < 0 || row.reducedRow >= rows) throw std::runtime_error("metadata reduced row out of range");
        metadata[static_cast<size_t>(row.reducedRow)] = row;
    }
    return metadata;
}

void loadNearNullspace(const std::string &path, ReplayData &data) {
    std::ifstream in(path);
    if (!in) throw std::runtime_error("cannot open " + path);
    std::string line;
    std::getline(in, line);
    const std::vector<std::string> header = splitTab(line);
    if (header.size() < 2) throw std::runtime_error("near_nullspace.tsv has no vector columns");
    data.nearNullspaceColumns = static_cast<int>(header.size()) - 1;
    data.nearNullspace.assign(static_cast<size_t>(data.fullRows) * data.nearNullspaceColumns, 0.0);
    while (std::getline(in, line)) {
        if (line.empty()) continue;
        const std::vector<std::string> values = splitTab(line);
        if (values.size() < header.size()) throw std::runtime_error("bad near-nullspace row");
        const int row = std::stoi(values[0]);
        if (row < 0 || row >= data.fullRows) throw std::runtime_error("near-nullspace row out of range");
        for (int col = 0; col < data.nearNullspaceColumns; col++) {
            data.nearNullspace[static_cast<size_t>(row) * data.nearNullspaceColumns + col] = std::stod(values[static_cast<size_t>(col + 1)]);
        }
    }
}

ReplayData loadReplay(const std::string &replayDir) {
    ReplayData data;
    data.reducedMatrix = loadMatrixMarket(pathJoin(replayDir, "matrix.mtx"));
    data.rhs = loadVector(pathJoin(replayDir, "rhs.tsv"), data.reducedMatrix.rows);
    data.referenceSolution = loadVector(pathJoin(replayDir, "solution.tsv"), data.reducedMatrix.rows);
    const std::map<std::string, std::string> summary = loadSummary(pathJoin(replayDir, "summary.tsv"));
    data.fullRows = summary.count("elastic_full_rows") ? std::stoi(summary.at("elastic_full_rows")) : data.reducedMatrix.rows;
    data.blockSize = summary.count("elastic_block_size") ? std::stoi(summary.at("elastic_block_size")) : 1;
    data.dimension = summary.count("elastic_dimension") ? std::stoi(summary.at("elastic_dimension")) : 0;
    data.metadata = loadMetadata(pathJoin(replayDir, "metadata.tsv"), data.reducedMatrix.rows);
    data.reducedToFull.assign(static_cast<size_t>(data.reducedMatrix.rows), -1);
    data.fullToReduced.assign(static_cast<size_t>(data.fullRows), -1);
    const int nodeCount = data.blockSize > 0 ? data.fullRows / data.blockSize : 0;
    data.nodeCoordinates.assign(static_cast<size_t>(std::max(0, nodeCount)), {0.0, 0.0, 0.0});
    data.nodeHasCoordinate.assign(static_cast<size_t>(std::max(0, nodeCount)), 0);
    for (const MetadataRow &row : data.metadata) {
        data.reducedToFull[static_cast<size_t>(row.reducedRow)] = row.fullElasticRow;
        if (row.fullElasticRow >= 0 && row.fullElasticRow < data.fullRows) {
            data.fullToReduced[static_cast<size_t>(row.fullElasticRow)] = row.reducedRow;
            const int node = data.blockSize > 0 ? row.fullElasticRow / data.blockSize : -1;
            if (node >= 0 && node < nodeCount) {
                data.nodeCoordinates[static_cast<size_t>(node)] = {row.x, row.y, row.z};
                data.nodeHasCoordinate[static_cast<size_t>(node)] = 1;
            }
        }
    }
    loadNearNullspace(pathJoin(replayDir, "near_nullspace.tsv"), data);
    return data;
}

Eigen::VectorXd spmv(const CsrMatrix &matrix, const Eigen::VectorXd &x) {
    Eigen::VectorXd y = Eigen::VectorXd::Zero(matrix.rows);
#ifdef _OPENMP
#pragma omp parallel for schedule(static)
#endif
    for (int row = 0; row < matrix.rows; row++) {
        double sum = 0.0;
        for (long long k = matrix.ptr[static_cast<size_t>(row)]; k < matrix.ptr[static_cast<size_t>(row + 1)]; k++) {
            sum += matrix.val[static_cast<size_t>(k)] * x[matrix.col[static_cast<size_t>(k)]];
        }
        y[row] = sum;
    }
    return y;
}

Eigen::VectorXd spmvTranspose(const CsrMatrix &matrix, const Eigen::VectorXd &x) {
    Eigen::VectorXd y = Eigen::VectorXd::Zero(matrix.cols);
    for (int row = 0; row < matrix.rows; row++) {
        const double xr = x[row];
        for (long long k = matrix.ptr[static_cast<size_t>(row)]; k < matrix.ptr[static_cast<size_t>(row + 1)]; k++) {
            y[matrix.col[static_cast<size_t>(k)]] += matrix.val[static_cast<size_t>(k)] * xr;
        }
    }
    return y;
}

double symmetryProbeRelative(const CsrMatrix &matrix) {
    if (matrix.rows != matrix.cols || matrix.rows <= 0) return std::numeric_limits<double>::infinity();
    Eigen::VectorXd x(matrix.rows);
    for (int i = 0; i < matrix.rows; i++) x[i] = std::sin(0.000137 * double(i + 1));
    const Eigen::VectorXd ax = spmv(matrix, x);
    const Eigen::VectorXd atx = spmvTranspose(matrix, x);
    return (ax - atx).norm() / std::max(ax.norm(), 1.0e-300);
}

CsrMatrix permuteCsr(const CsrMatrix &matrix, const std::vector<int> &oldToNew) {
    if (oldToNew.empty()) return matrix;
    CsrMatrix out;
    out.rows = matrix.rows;
    out.cols = matrix.cols;
    out.ptr.assign(static_cast<size_t>(out.rows) + 1, 0);
    for (int oldRow = 0; oldRow < matrix.rows; oldRow++) {
        const int newRow = oldToNew[static_cast<size_t>(oldRow)];
        out.ptr[static_cast<size_t>(newRow + 1)] += matrix.ptr[static_cast<size_t>(oldRow + 1)] - matrix.ptr[static_cast<size_t>(oldRow)];
    }
    for (int row = 0; row < out.rows; row++) out.ptr[static_cast<size_t>(row + 1)] += out.ptr[static_cast<size_t>(row)];
    out.col.assign(matrix.col.size(), 0);
    out.val.assign(matrix.val.size(), 0.0);
    std::vector<long long> cursor = out.ptr;
    for (int oldRow = 0; oldRow < matrix.rows; oldRow++) {
        const int newRow = oldToNew[static_cast<size_t>(oldRow)];
        for (long long k = matrix.ptr[static_cast<size_t>(oldRow)]; k < matrix.ptr[static_cast<size_t>(oldRow + 1)]; k++) {
            const long long index = cursor[static_cast<size_t>(newRow)]++;
            out.col[static_cast<size_t>(index)] = oldToNew[static_cast<size_t>(matrix.col[static_cast<size_t>(k)])];
            out.val[static_cast<size_t>(index)] = matrix.val[static_cast<size_t>(k)];
        }
    }
    return out;
}

CsrMatrix liftToFull(const CsrMatrix &matrix, const ReplayData &data) {
    CsrMatrix out;
    out.rows = data.fullRows;
    out.cols = data.fullRows;
    out.ptr.assign(static_cast<size_t>(out.rows) + 1, 0);
    std::vector<unsigned char> hasRow(static_cast<size_t>(data.fullRows), 0);
    for (int row = 0; row < matrix.rows; row++) {
        const int fullRow = data.reducedToFull[static_cast<size_t>(row)];
        if (fullRow < 0 || fullRow >= data.fullRows) throw std::runtime_error("invalid reduced-to-full row map");
        out.ptr[static_cast<size_t>(fullRow + 1)] += matrix.ptr[static_cast<size_t>(row + 1)] - matrix.ptr[static_cast<size_t>(row)];
        hasRow[static_cast<size_t>(fullRow)] = 1;
    }
    for (int row = 0; row < data.fullRows; row++) {
        if (!hasRow[static_cast<size_t>(row)]) out.ptr[static_cast<size_t>(row + 1)]++;
    }
    for (int row = 0; row < out.rows; row++) out.ptr[static_cast<size_t>(row + 1)] += out.ptr[static_cast<size_t>(row)];
    out.col.assign(static_cast<size_t>(out.ptr.back()), 0);
    out.val.assign(static_cast<size_t>(out.ptr.back()), 0.0);
    std::vector<long long> cursor = out.ptr;
    for (int row = 0; row < matrix.rows; row++) {
        const int fullRow = data.reducedToFull[static_cast<size_t>(row)];
        for (long long k = matrix.ptr[static_cast<size_t>(row)]; k < matrix.ptr[static_cast<size_t>(row + 1)]; k++) {
            const long long index = cursor[static_cast<size_t>(fullRow)]++;
            const int reducedCol = matrix.col[static_cast<size_t>(k)];
            out.col[static_cast<size_t>(index)] = data.reducedToFull[static_cast<size_t>(reducedCol)];
            out.val[static_cast<size_t>(index)] = matrix.val[static_cast<size_t>(k)];
        }
    }
    for (int row = 0; row < data.fullRows; row++) {
        if (!hasRow[static_cast<size_t>(row)]) {
            const long long index = cursor[static_cast<size_t>(row)]++;
            out.col[static_cast<size_t>(index)] = row;
            out.val[static_cast<size_t>(index)] = 1.0;
        }
    }
    return out;
}

std::vector<int> invertPermutation(const std::vector<int> &oldToNew) {
    std::vector<int> newToOld(oldToNew.size(), -1);
    for (int oldIndex = 0; oldIndex < static_cast<int>(oldToNew.size()); oldIndex++) {
        newToOld[static_cast<size_t>(oldToNew[static_cast<size_t>(oldIndex)])] = oldIndex;
    }
    return newToOld;
}

std::vector<int> coordinateNodePermutation(const ReplayData &data) {
    const int nodeCount = data.blockSize > 0 ? data.fullRows / data.blockSize : 0;
    std::vector<int> nodes(static_cast<size_t>(nodeCount));
    std::iota(nodes.begin(), nodes.end(), 0);
    std::stable_sort(nodes.begin(), nodes.end(), [&](int a, int b) {
        const bool ah = data.nodeHasCoordinate[static_cast<size_t>(a)] != 0;
        const bool bh = data.nodeHasCoordinate[static_cast<size_t>(b)] != 0;
        if (ah != bh) return ah > bh;
        if (!ah) return a < b;
        const auto &ca = data.nodeCoordinates[static_cast<size_t>(a)];
        const auto &cb = data.nodeCoordinates[static_cast<size_t>(b)];
        if (ca[0] != cb[0]) return ca[0] < cb[0];
        if (ca[1] != cb[1]) return ca[1] < cb[1];
        if (ca[2] != cb[2]) return ca[2] < cb[2];
        return a < b;
    });
    std::vector<int> oldToNew(static_cast<size_t>(data.fullRows), 0);
    for (int newNode = 0; newNode < nodeCount; newNode++) {
        const int oldNode = nodes[static_cast<size_t>(newNode)];
        for (int dof = 0; dof < data.blockSize; dof++) {
            oldToNew[static_cast<size_t>(oldNode * data.blockSize + dof)] = newNode * data.blockSize + dof;
        }
    }
    return oldToNew;
}

std::vector<int> componentMajorPermutation(const ReplayData &data) {
    const int nodeCount = data.blockSize > 0 ? data.fullRows / data.blockSize : 0;
    std::vector<int> oldToNew(static_cast<size_t>(data.fullRows), 0);
    for (int node = 0; node < nodeCount; node++) {
        for (int dof = 0; dof < data.blockSize; dof++) {
            oldToNew[static_cast<size_t>(node * data.blockSize + dof)] = dof * nodeCount + node;
        }
    }
    return oldToNew;
}

std::vector<int> rcmPermutationFromMatrix(const CsrMatrix &matrix, int blockSize) {
    if (blockSize <= 0 || matrix.rows % blockSize != 0) return {};
    const int nodeCount = matrix.rows / blockSize;
    std::vector<std::vector<int>> adjacency(static_cast<size_t>(nodeCount));
    for (int row = 0; row < matrix.rows; row++) {
        const int rowNode = row / blockSize;
        for (long long k = matrix.ptr[static_cast<size_t>(row)]; k < matrix.ptr[static_cast<size_t>(row + 1)]; k++) {
            const int colNode = matrix.col[static_cast<size_t>(k)] / blockSize;
            if (rowNode == colNode || colNode < 0 || colNode >= nodeCount) continue;
            adjacency[static_cast<size_t>(rowNode)].push_back(colNode);
            adjacency[static_cast<size_t>(colNode)].push_back(rowNode);
        }
    }
    std::vector<int> degree(static_cast<size_t>(nodeCount), 0);
    for (int node = 0; node < nodeCount; node++) {
        auto &neighbors = adjacency[static_cast<size_t>(node)];
        std::sort(neighbors.begin(), neighbors.end());
        neighbors.erase(std::unique(neighbors.begin(), neighbors.end()), neighbors.end());
        degree[static_cast<size_t>(node)] = static_cast<int>(neighbors.size());
    }
    std::vector<int> candidates(static_cast<size_t>(nodeCount));
    std::iota(candidates.begin(), candidates.end(), 0);
    std::stable_sort(candidates.begin(), candidates.end(), [&](int a, int b) {
        if (degree[static_cast<size_t>(a)] != degree[static_cast<size_t>(b)]) {
            return degree[static_cast<size_t>(a)] < degree[static_cast<size_t>(b)];
        }
        return a < b;
    });
    std::vector<unsigned char> visited(static_cast<size_t>(nodeCount), 0);
    std::vector<int> cmOrder;
    cmOrder.reserve(static_cast<size_t>(nodeCount));
    std::deque<int> queue;
    std::vector<int> next;
    for (int start : candidates) {
        if (visited[static_cast<size_t>(start)]) continue;
        visited[static_cast<size_t>(start)] = 1;
        queue.push_back(start);
        while (!queue.empty()) {
            const int node = queue.front();
            queue.pop_front();
            cmOrder.push_back(node);
            next.clear();
            for (int neighbor : adjacency[static_cast<size_t>(node)]) {
                if (!visited[static_cast<size_t>(neighbor)]) next.push_back(neighbor);
            }
            std::stable_sort(next.begin(), next.end(), [&](int a, int b) {
                if (degree[static_cast<size_t>(a)] != degree[static_cast<size_t>(b)]) {
                    return degree[static_cast<size_t>(a)] < degree[static_cast<size_t>(b)];
                }
                return a < b;
            });
            for (int neighbor : next) {
                if (!visited[static_cast<size_t>(neighbor)]) {
                    visited[static_cast<size_t>(neighbor)] = 1;
                    queue.push_back(neighbor);
                }
            }
        }
    }
    std::reverse(cmOrder.begin(), cmOrder.end());
    std::vector<int> oldToNew(static_cast<size_t>(matrix.rows), 0);
    for (int newNode = 0; newNode < nodeCount; newNode++) {
        const int oldNode = cmOrder[static_cast<size_t>(newNode)];
        for (int dof = 0; dof < blockSize; dof++) {
            oldToNew[static_cast<size_t>(oldNode * blockSize + dof)] = newNode * blockSize + dof;
        }
    }
    return oldToNew;
}

template <class F>
Eigen::VectorXd mapVector(int size, F oldIndexForNew, const Eigen::VectorXd &oldVector) {
    Eigen::VectorXd out(size);
    for (int newIndex = 0; newIndex < size; newIndex++) out[newIndex] = oldVector[oldIndexForNew(newIndex)];
    return out;
}

std::vector<double> permuteNullspace(const std::vector<double> &oldValues, int columns, const std::vector<int> &oldToNew) {
    if (oldToNew.empty()) return oldValues;
    std::vector<double> out(oldValues.size(), 0.0);
    for (int oldRow = 0; oldRow < static_cast<int>(oldToNew.size()); oldRow++) {
        const int newRow = oldToNew[static_cast<size_t>(oldRow)];
        for (int col = 0; col < columns; col++) {
            out[static_cast<size_t>(newRow) * columns + col] = oldValues[static_cast<size_t>(oldRow) * columns + col];
        }
    }
    return out;
}

struct NullspaceBuildOptions {
    bool usePhysicalCoordinates = false;
    bool negativeRotationalRows = false;
    bool zeroRotationalRows = false;
    bool translationRotationOnly = false;
    double coordinateScale = std::numeric_limits<double>::quiet_NaN();
};

std::pair<std::array<double, 3>, double> coordinateCenterAndScale(const ReplayData &data) {
    std::array<double, 3> center = {0.0, 0.0, 0.0};
    int count = 0;
    for (int node = 0; node < static_cast<int>(data.nodeCoordinates.size()); node++) {
        if (!data.nodeHasCoordinate[static_cast<size_t>(node)]) continue;
        const auto &point = data.nodeCoordinates[static_cast<size_t>(node)];
        center[0] += point[0];
        center[1] += point[1];
        center[2] += point[2];
        count++;
    }
    if (count > 0) {
        center[0] /= count;
        center[1] /= count;
        center[2] /= count;
    }

    double radius2 = 0.0;
    for (int node = 0; node < static_cast<int>(data.nodeCoordinates.size()); node++) {
        if (!data.nodeHasCoordinate[static_cast<size_t>(node)]) continue;
        const auto &point = data.nodeCoordinates[static_cast<size_t>(node)];
        const double x = point[0] - center[0];
        const double y = point[1] - center[1];
        const double z = point[2] - center[2];
        radius2 += x * x + y * y + z * z;
    }
    const double scale = count > 0 ? std::sqrt(radius2 / count) : 1.0;
    return {center, scale > 0.0 ? scale : 1.0};
}

void orthonormalizeColumns(std::vector<double> &values, int rows, int columns) {
    std::vector<double> orthonormal(values.size(), 0.0);
    int accepted = 0;
    for (int col = 0; col < columns; col++) {
        std::vector<double> candidate(static_cast<size_t>(rows), 0.0);
        for (int row = 0; row < rows; row++) {
            candidate[static_cast<size_t>(row)] = values[static_cast<size_t>(row) * columns + col];
        }
        for (int previous = 0; previous < accepted; previous++) {
            double dot = 0.0;
            for (int row = 0; row < rows; row++) {
                dot += candidate[static_cast<size_t>(row)]
                    * orthonormal[static_cast<size_t>(row) * columns + previous];
            }
            for (int row = 0; row < rows; row++) {
                candidate[static_cast<size_t>(row)] -= dot * orthonormal[static_cast<size_t>(row) * columns + previous];
            }
        }
        double norm2 = 0.0;
        for (double value : candidate) norm2 += value * value;
        const double norm = std::sqrt(norm2);
        if (norm <= 1.0e-12) continue;
        for (int row = 0; row < rows; row++) {
            orthonormal[static_cast<size_t>(row) * columns + accepted] =
                candidate[static_cast<size_t>(row)] / norm;
        }
        accepted++;
    }
    if (accepted != columns) {
        throw std::runtime_error("rigid-body near-nullspace lost rank during orthonormalization");
    }
    values.swap(orthonormal);
}

std::vector<double> buildRigidParticleNullspace(
    const ReplayData &data,
    const VariantData &variant,
    const NullspaceBuildOptions &options
) {
    if (data.dimension != 3 || data.blockSize < 6) {
        throw std::runtime_error("rigid-body near-nullspace variants currently require 3D six-DOF particles");
    }
    const int columns = 6;
    std::vector<double> values(static_cast<size_t>(variant.matrix.rows) * columns, 0.0);
    const auto centerScale = coordinateCenterAndScale(data);
    const auto &center = centerScale.first;
    double scale = options.usePhysicalCoordinates ? 1.0 : centerScale.second;
    if (std::isfinite(options.coordinateScale) && options.coordinateScale > 0.0) {
        scale = options.coordinateScale;
    }
    const double rotationalSign = options.negativeRotationalRows ? -1.0 : 1.0;

    for (int active = 0; active < static_cast<int>(variant.activeToReduced.size()); active++) {
        const int reduced = variant.activeToReduced[static_cast<size_t>(active)];
        if (reduced < 0) continue;
        const MetadataRow &row = data.metadata[static_cast<size_t>(reduced)];
        const int relative = row.relativeDof;
        const double x = (row.x - center[0]) / scale;
        const double y = (row.y - center[1]) / scale;
        const double z = (row.z - center[2]) / scale;
        double *target = &values[static_cast<size_t>(active) * columns];

        if (options.translationRotationOnly) {
            if (relative < 3) {
                target[relative] = 1.0;
            } else if (relative < 6) {
                target[relative] = rotationalSign;
            }
        } else {
            if (relative < 3) {
                target[relative] = 1.0;
            }
            if (relative == 0) {
                target[4] = z;
                target[5] = -y;
            } else if (relative == 1) {
                target[3] = -z;
                target[5] = x;
            } else if (relative == 2) {
                target[3] = y;
                target[4] = -x;
            } else if (!options.zeroRotationalRows && relative >= 3 && relative < 6) {
                target[relative] = rotationalSign;
            }
        }
    }

    orthonormalizeColumns(values, variant.matrix.rows, columns);
    return values;
}

bool nullspaceRequested(const std::string &mode) {
    return mode != "off";
}

void applyNullspaceMode(const ReplayData &data, VariantData &variant, const Options &runOptions) {
    const std::string &mode = runOptions.nullspace;
    if (mode == "off" || mode == "rbm6") {
        return;
    }

    NullspaceBuildOptions options;
    if (mode == "rbm6-normalized") {
    } else if (mode == "rbm6-physical") {
        options.usePhysicalCoordinates = true;
    } else if (mode == "rbm6-scale") {
        options.usePhysicalCoordinates = true;
    } else if (mode == "rbm6-rot-neg") {
        options.negativeRotationalRows = true;
    } else if (mode == "rbm6-rot-zero") {
        options.zeroRotationalRows = true;
    } else if (mode == "transrot6") {
        options.translationRotationOnly = true;
    } else {
        throw std::runtime_error("unknown nullspace mode: " + mode);
    }
    options.coordinateScale = runOptions.nullspaceCoordinateScale;
    variant.nearNullspaceColumns = 6;
    variant.nearNullspace = buildRigidParticleNullspace(data, variant, options);
}

std::vector<double> restrictNullspaceToReduced(const ReplayData &data, const std::vector<int> &activeToReduced) {
    std::vector<double> out(static_cast<size_t>(activeToReduced.size()) * data.nearNullspaceColumns, 0.0);
    for (int active = 0; active < static_cast<int>(activeToReduced.size()); active++) {
        const int reduced = activeToReduced[static_cast<size_t>(active)];
        const int full = data.reducedToFull[static_cast<size_t>(reduced)];
        for (int col = 0; col < data.nearNullspaceColumns; col++) {
            out[static_cast<size_t>(active) * data.nearNullspaceColumns + col] =
                data.nearNullspace[static_cast<size_t>(full) * data.nearNullspaceColumns + col];
        }
    }
    return out;
}

VariantData makeVariant(const ReplayData &data, const std::string &representation) {
    VariantData out;
    out.blockSize = data.blockSize;
    out.dimension = data.dimension;
    out.nearNullspaceColumns = data.nearNullspaceColumns;

    if (representation == "reduced-original") {
        out.matrix = data.reducedMatrix;
        out.rhs = data.rhs;
        out.referenceInActiveOrder = data.referenceSolution;
        out.activeToReduced.resize(static_cast<size_t>(data.reducedMatrix.rows));
        std::iota(out.activeToReduced.begin(), out.activeToReduced.end(), 0);
        out.nearNullspace = restrictNullspaceToReduced(data, out.activeToReduced);
        out.systemAmgCompatible = false;
        out.interpolationVectorsCompatible = false;
        out.dofMode = "scalar";
        return out;
    }

    if (representation == "compressed-free-node-major") {
        std::vector<int> reducedOrder(static_cast<size_t>(data.reducedMatrix.rows));
        std::iota(reducedOrder.begin(), reducedOrder.end(), 0);
        std::stable_sort(reducedOrder.begin(), reducedOrder.end(), [&](int a, int b) {
            return data.reducedToFull[static_cast<size_t>(a)] < data.reducedToFull[static_cast<size_t>(b)];
        });
        std::vector<int> oldToNew(static_cast<size_t>(data.reducedMatrix.rows), 0);
        out.activeToReduced.assign(static_cast<size_t>(data.reducedMatrix.rows), 0);
        for (int newIndex = 0; newIndex < data.reducedMatrix.rows; newIndex++) {
            const int oldIndex = reducedOrder[static_cast<size_t>(newIndex)];
            oldToNew[static_cast<size_t>(oldIndex)] = newIndex;
            out.activeToReduced[static_cast<size_t>(newIndex)] = oldIndex;
        }
        out.matrix = permuteCsr(data.reducedMatrix, oldToNew);
        out.rhs = mapVector(out.matrix.rows, [&](int newIndex) { return out.activeToReduced[static_cast<size_t>(newIndex)]; }, data.rhs);
        out.referenceInActiveOrder = mapVector(out.matrix.rows, [&](int newIndex) { return out.activeToReduced[static_cast<size_t>(newIndex)]; }, data.referenceSolution);
        out.nearNullspace = restrictNullspaceToReduced(data, out.activeToReduced);
        std::vector<int> activeToFull(out.activeToReduced.size(), 0);
        for (int i = 0; i < static_cast<int>(activeToFull.size()); i++) {
            activeToFull[static_cast<size_t>(i)] = data.reducedToFull[static_cast<size_t>(out.activeToReduced[static_cast<size_t>(i)])];
        }
        (void)activeToFull;
        out.systemAmgCompatible = false;
        out.interpolationVectorsCompatible = false;
        out.dofMode = "scalar_free_node_major";
        return out;
    }

    out.matrix = liftToFull(data.reducedMatrix, data);
    out.rhs = Eigen::VectorXd::Zero(data.fullRows);
    out.referenceInActiveOrder = Eigen::VectorXd::Zero(data.fullRows);
    out.activeToReduced.assign(static_cast<size_t>(data.fullRows), -1);
    for (int reduced = 0; reduced < data.reducedMatrix.rows; reduced++) {
        const int full = data.reducedToFull[static_cast<size_t>(reduced)];
        out.rhs[full] = data.rhs[reduced];
        out.referenceInActiveOrder[full] = data.referenceSolution[reduced];
        out.activeToReduced[static_cast<size_t>(full)] = reduced;
    }
    out.nearNullspace = data.nearNullspace;
    out.systemAmgCompatible = true;
    out.interpolationVectorsCompatible = true;
    out.dofMode = "cyclic_node_major";

    std::vector<int> oldToNew;
    if (representation == "full-coordinate-node-major") {
        oldToNew = coordinateNodePermutation(data);
    } else if (representation == "full-rcm-node-major") {
        oldToNew = rcmPermutationFromMatrix(out.matrix, data.blockSize);
    } else if (representation == "full-component-major") {
        oldToNew = componentMajorPermutation(data);
        out.dofMode = "explicit_component_major";
    } else if (representation != "full-node-major") {
        throw std::runtime_error("unknown representation: " + representation);
    }

    if (!oldToNew.empty()) {
        const std::vector<int> newToOld = invertPermutation(oldToNew);
        out.matrix = permuteCsr(out.matrix, oldToNew);
        out.rhs = mapVector(out.matrix.rows, [&](int newIndex) { return newToOld[static_cast<size_t>(newIndex)]; }, out.rhs);
        out.referenceInActiveOrder = mapVector(out.matrix.rows, [&](int newIndex) { return newToOld[static_cast<size_t>(newIndex)]; }, out.referenceInActiveOrder);
        out.nearNullspace = permuteNullspace(out.nearNullspace, out.nearNullspaceColumns, oldToNew);
        std::vector<int> remapped(out.activeToReduced.size(), -1);
        for (int oldIndex = 0; oldIndex < static_cast<int>(out.activeToReduced.size()); oldIndex++) {
            remapped[static_cast<size_t>(oldToNew[static_cast<size_t>(oldIndex)])] = out.activeToReduced[static_cast<size_t>(oldIndex)];
        }
        out.activeToReduced.swap(remapped);
    }

    if (representation == "full-component-major") {
        const int nodeCount = data.blockSize > 0 ? data.fullRows / data.blockSize : 0;
        out.dofFunc.assign(static_cast<size_t>(data.fullRows), 0);
        for (int dof = 0; dof < data.blockSize; dof++) {
            for (int node = 0; node < nodeCount; node++) {
                out.dofFunc[static_cast<size_t>(dof * nodeCount + node)] = dof;
            }
        }
        // hypre's nodal BoomerAMG assumes a node-compatible layout even when
        // dof_func is supplied. Treat component-major as a scalar reordered
        // diagnostic in this replay harness to avoid unsafe hypre ownership
        // paths observed with explicit dof_func on this matrix.
        out.dofFunc.clear();
        out.systemAmgCompatible = false;
        out.interpolationVectorsCompatible = false;
        out.dofMode = "scalar_component_major";
    }
    return out;
}

void createHypreVector(const Eigen::VectorXd &values, HYPRE_IJVector &ijVector, HYPRE_ParVector &parVector) {
    const HYPRE_BigInt first = 0;
    const HYPRE_BigInt last = values.size() - 1;
    hypreCheck(HYPRE_IJVectorCreate(MPI_COMM_SELF, first, last, &ijVector), "HYPRE_IJVectorCreate");
    hypreCheck(HYPRE_IJVectorSetObjectType(ijVector, HYPRE_PARCSR), "HYPRE_IJVectorSetObjectType");
    hypreCheck(HYPRE_IJVectorInitialize(ijVector), "HYPRE_IJVectorInitialize");
    std::vector<HYPRE_BigInt> indices(static_cast<size_t>(values.size()));
    for (int i = 0; i < values.size(); i++) indices[static_cast<size_t>(i)] = i;
    hypreCheck(HYPRE_IJVectorSetValues(ijVector, static_cast<HYPRE_Int>(values.size()), indices.data(), values.data()), "HYPRE_IJVectorSetValues");
    hypreCheck(HYPRE_IJVectorAssemble(ijVector), "HYPRE_IJVectorAssemble");
    hypreCheck(HYPRE_IJVectorGetObject(ijVector, reinterpret_cast<void **>(&parVector)), "HYPRE_IJVectorGetObject");
}

void buildHypreMatrix(const CsrMatrix &matrix, HypreObjects &hypre) {
    const HYPRE_BigInt first = 0;
    const HYPRE_BigInt last = matrix.rows - 1;
    hypreCheck(HYPRE_IJMatrixCreate(MPI_COMM_SELF, first, last, first, last, &hypre.ijMatrix), "HYPRE_IJMatrixCreate");
    hypreCheck(HYPRE_IJMatrixSetObjectType(hypre.ijMatrix, HYPRE_PARCSR), "HYPRE_IJMatrixSetObjectType");
    std::vector<HYPRE_Int> rowSizes(static_cast<size_t>(matrix.rows));
    for (int row = 0; row < matrix.rows; row++) {
        rowSizes[static_cast<size_t>(row)] = static_cast<HYPRE_Int>(matrix.ptr[static_cast<size_t>(row + 1)] - matrix.ptr[static_cast<size_t>(row)]);
    }
    HYPRE_IJMatrixSetRowSizes(hypre.ijMatrix, rowSizes.data());
    HYPRE_IJMatrixSetOMPFlag(hypre.ijMatrix, 1);
    hypreCheck(HYPRE_IJMatrixInitialize(hypre.ijMatrix), "HYPRE_IJMatrixInitialize");
    std::vector<HYPRE_BigInt> rowIds(static_cast<size_t>(matrix.rows));
    std::vector<HYPRE_BigInt> cols(matrix.col.size());
    for (int row = 0; row < matrix.rows; row++) rowIds[static_cast<size_t>(row)] = row;
#ifdef _OPENMP
#pragma omp parallel for schedule(static)
#endif
    for (long long k = 0; k < static_cast<long long>(matrix.col.size()); k++) {
        cols[static_cast<size_t>(k)] = matrix.col[static_cast<size_t>(k)];
    }
    hypreCheck(
        HYPRE_IJMatrixSetValues(
            hypre.ijMatrix,
            static_cast<HYPRE_Int>(matrix.rows),
            rowSizes.data(),
            rowIds.data(),
            cols.data(),
            matrix.val.data()
        ),
        "HYPRE_IJMatrixSetValues"
    );
    hypreCheck(HYPRE_IJMatrixAssemble(hypre.ijMatrix), "HYPRE_IJMatrixAssemble");
    hypreCheck(HYPRE_IJMatrixGetObject(hypre.ijMatrix, reinterpret_cast<void **>(&hypre.parMatrix)), "HYPRE_IJMatrixGetObject");
}

void setBoomerOptions(HYPRE_Solver precond, const Options &options, const VariantData &variant, int &usedInterpVectors, HypreObjects &hypre) {
    HYPRE_BoomerAMGSetTol(precond, 0.0);
    HYPRE_BoomerAMGSetMaxIter(precond, options.boomerMaxIterations);
    HYPRE_BoomerAMGSetCoarsenType(precond, 8);
    HYPRE_BoomerAMGSetInterpType(precond, options.interpType);
    HYPRE_BoomerAMGSetStrongThreshold(precond, options.strongThreshold);
    HYPRE_BoomerAMGSetRelaxType(precond, options.relaxType);
    HYPRE_BoomerAMGSetRelaxOrder(precond, 0);
    if (options.numSweeps > 0) HYPRE_BoomerAMGSetNumSweeps(precond, options.numSweeps);
    HYPRE_BoomerAMGSetPMaxElmts(precond, options.pMaxElmts);
    HYPRE_BoomerAMGSetAggNumLevels(precond, options.aggNumLevels);
    if (options.nodalDiag != 0) HYPRE_BoomerAMGSetNodalDiag(precond, options.nodalDiag);
    if (options.keepSameSign != 0) HYPRE_BoomerAMGSetKeepSameSign(precond, options.keepSameSign);
    if (options.chebyOrder > 0) HYPRE_BoomerAMGSetChebyOrder(precond, options.chebyOrder);
    if (options.chebyFraction > 0.0) HYPRE_BoomerAMGSetChebyFraction(precond, options.chebyFraction);
    if (options.nonGalerkinTol >= 0.0) {
        HYPRE_Real tolerances[3] = {0.0, options.nonGalerkinTol, options.nonGalerkinTol};
        HYPRE_BoomerAMGSetNonGalerkTol(precond, 3, tolerances);
    }
    if (variant.systemAmgCompatible && variant.blockSize > 1) {
        HYPRE_BoomerAMGSetNumFunctions(precond, variant.blockSize);
        if (!variant.dofFunc.empty()) {
            hypre.dofFuncStorage.assign(variant.dofFunc.begin(), variant.dofFunc.end());
            HYPRE_BoomerAMGSetDofFunc(precond, hypre.dofFuncStorage.data());
        }
        if (options.nodal > 0) HYPRE_BoomerAMGSetNodal(precond, options.nodal);
    }
    if (nullspaceRequested(options.nullspace)) {
        if (!variant.interpolationVectorsCompatible || options.nodal <= 0 || variant.nearNullspaceColumns <= 0) {
            throw std::runtime_error("near-nullspace requires a system-compatible nodal representation");
        }
        HYPRE_BoomerAMGSetInterpVecVariant(precond, options.interpVecVariant);
        if (options.interpVecQMax > 0) HYPRE_BoomerAMGSetInterpVecQMax(precond, options.interpVecQMax);
        if (options.interpVecAbsQTrunc > 0.0) HYPRE_BoomerAMGSetInterpVecAbsQTrunc(precond, options.interpVecAbsQTrunc);
        hypre.interpIjVectors.resize(static_cast<size_t>(variant.nearNullspaceColumns), nullptr);
        hypre.interpParVectors.resize(static_cast<size_t>(variant.nearNullspaceColumns), nullptr);
        for (int col = 0; col < variant.nearNullspaceColumns; col++) {
            Eigen::VectorXd mode(variant.matrix.rows);
            for (int row = 0; row < variant.matrix.rows; row++) {
                mode[row] = variant.nearNullspace[static_cast<size_t>(row) * variant.nearNullspaceColumns + col];
            }
            createHypreVector(mode, hypre.interpIjVectors[static_cast<size_t>(col)], hypre.interpParVectors[static_cast<size_t>(col)]);
        }
        HYPRE_BoomerAMGSetInterpVectors(precond, variant.nearNullspaceColumns, hypre.interpParVectors.data());
        usedInterpVectors = variant.nearNullspaceColumns;
    }
}

Eigen::VectorXd readHypreVector(HYPRE_IJVector vector, int rows) {
    Eigen::VectorXd values(rows);
    std::vector<HYPRE_BigInt> indices(static_cast<size_t>(rows));
    for (int i = 0; i < rows; i++) indices[static_cast<size_t>(i)] = i;
    hypreCheck(HYPRE_IJVectorGetValues(vector, rows, indices.data(), values.data()), "HYPRE_IJVectorGetValues");
    return values;
}

void configureSolver(HypreObjects &hypre, const Options &options) {
    hypre.solverKind = options.solver;
    if (options.solver == "pcg") {
        hypreCheck(HYPRE_ParCSRPCGCreate(MPI_COMM_SELF, &hypre.solver), "HYPRE_ParCSRPCGCreate");
        HYPRE_ParCSRPCGSetTol(hypre.solver, options.tolerance);
        HYPRE_ParCSRPCGSetAbsoluteTol(hypre.solver, 0.0);
        HYPRE_ParCSRPCGSetMaxIter(hypre.solver, options.maxIterations);
        HYPRE_ParCSRPCGSetTwoNorm(hypre.solver, 1);
        HYPRE_ParCSRPCGSetRelChange(hypre.solver, 0);
        HYPRE_ParCSRPCGSetPrintLevel(hypre.solver, 0);
        HYPRE_ParCSRPCGSetLogging(hypre.solver, 1);
        HYPRE_ParCSRPCGSetPrecond(hypre.solver, HYPRE_BoomerAMGSolve, HYPRE_BoomerAMGSetup, hypre.precond);
    } else {
        hypreCheck(HYPRE_ParCSRFlexGMRESCreate(MPI_COMM_SELF, &hypre.solver), "HYPRE_ParCSRFlexGMRESCreate");
        HYPRE_ParCSRFlexGMRESSetTol(hypre.solver, options.tolerance);
        HYPRE_ParCSRFlexGMRESSetAbsoluteTol(hypre.solver, 0.0);
        HYPRE_ParCSRFlexGMRESSetMaxIter(hypre.solver, options.maxIterations);
        HYPRE_ParCSRFlexGMRESSetKDim(hypre.solver, 80);
        HYPRE_ParCSRFlexGMRESSetPrintLevel(hypre.solver, 0);
        HYPRE_ParCSRFlexGMRESSetLogging(hypre.solver, 1);
        HYPRE_ParCSRFlexGMRESSetPrecond(hypre.solver, HYPRE_BoomerAMGSolve, HYPRE_BoomerAMGSetup, hypre.precond);
    }
}

void setupKrylov(HypreObjects &hypre, const Options &options) {
    if (options.solver == "pcg") {
        hypreCheck(HYPRE_ParCSRPCGSetup(hypre.solver, hypre.parMatrix, hypre.parB, hypre.parX), "HYPRE_ParCSRPCGSetup");
    } else {
        hypreCheck(HYPRE_ParCSRFlexGMRESSetup(hypre.solver, hypre.parMatrix, hypre.parB, hypre.parX), "HYPRE_ParCSRFlexGMRESSetup");
    }
}

void solveKrylov(HypreObjects &hypre, const Options &options) {
    if (options.solver == "pcg") {
        hypreCheck(HYPRE_ParCSRPCGSolve(hypre.solver, hypre.parMatrix, hypre.parB, hypre.parX), "HYPRE_ParCSRPCGSolve");
    } else {
        hypreCheck(HYPRE_ParCSRFlexGMRESSolve(hypre.solver, hypre.parMatrix, hypre.parB, hypre.parX), "HYPRE_ParCSRFlexGMRESSolve");
    }
}

void getKrylovStatus(HypreObjects &hypre, const Options &options, int &iterations, double &residual) {
    HYPRE_Int hypreIterations = -1;
    HYPRE_Real hypreResidual = -1.0;
    if (options.solver == "pcg") {
        HYPRE_ParCSRPCGGetNumIterations(hypre.solver, &hypreIterations);
        HYPRE_ParCSRPCGGetFinalRelativeResidualNorm(hypre.solver, &hypreResidual);
    } else {
        HYPRE_ParCSRFlexGMRESGetNumIterations(hypre.solver, &hypreIterations);
        HYPRE_ParCSRFlexGMRESGetFinalRelativeResidualNorm(hypre.solver, &hypreResidual);
    }
    iterations = hypreIterations;
    residual = hypreResidual;
}

ResultRow runBenchmark(const Options &options, const ReplayData &replay) {
    ResultRow row;
    row.options = options;
    row.fullRows = replay.fullRows;
    row.blockSize = replay.blockSize;
    row.nearNullspaceColumns = replay.nearNullspaceColumns;
    try {
        initializeHypre(options.threads);
#ifdef _OPENMP
        row.ompMaxThreads = omp_get_max_threads();
#else
        row.ompMaxThreads = 1;
#endif
        row.ompNumThreads = envOrEmpty("OMP_NUM_THREADS");
        row.ompDynamic = envOrEmpty("OMP_DYNAMIC");
        row.ompProcBind = envOrEmpty("OMP_PROC_BIND");
        row.ompPlaces = envOrEmpty("OMP_PLACES");
        const auto variantStart = std::chrono::steady_clock::now();
        VariantData variant = makeVariant(replay, options.representation);
        applyNullspaceMode(replay, variant, options);
        row.timings.matrixSeconds = secondsSince(variantStart);
        row.rows = variant.matrix.rows;
        row.cols = variant.matrix.cols;
        row.nnz = variant.matrix.nnz();
        row.dofMode = variant.dofMode;
        row.symmetryProbeRelative = symmetryProbeRelative(variant.matrix);
        if (options.solver == "pcg" && row.symmetryProbeRelative > 1.0e-8) {
            throw std::runtime_error("pcg skipped because matrix symmetry probe is too large");
        }
        HypreObjects hypre;
        const auto hypreMatrixStart = std::chrono::steady_clock::now();
        buildHypreMatrix(variant.matrix, hypre);
        row.timings.matrixSeconds += secondsSince(hypreMatrixStart);
        const auto vectorStart = std::chrono::steady_clock::now();
        createHypreVector(variant.rhs, hypre.ijB, hypre.parB);
        createHypreVector(Eigen::VectorXd::Zero(variant.matrix.rows), hypre.ijX, hypre.parX);
        row.timings.vectorSeconds = secondsSince(vectorStart);
        hypreCheck(HYPRE_BoomerAMGCreate(&hypre.precond), "HYPRE_BoomerAMGCreate");
        setBoomerOptions(hypre.precond, options, variant, row.usedInterpVectors, hypre);

        const auto precondStart = std::chrono::steady_clock::now();
        hypreCheck(HYPRE_BoomerAMGSetup(hypre.precond, hypre.parMatrix, hypre.parB, hypre.parX), "HYPRE_BoomerAMGSetup");
        row.timings.preconditionerSetupSeconds = secondsSince(precondStart);

        for (int i = 0; i < options.applyWarmup; i++) {
            HYPRE_IJVectorSetConstantValues(hypre.ijX, 0.0);
            hypreCheck(HYPRE_BoomerAMGSolve(hypre.precond, hypre.parMatrix, hypre.parB, hypre.parX), "HYPRE_BoomerAMGSolve warmup");
        }
        const auto applyStart = std::chrono::steady_clock::now();
        for (int i = 0; i < options.applyRepeat; i++) {
            HYPRE_IJVectorSetConstantValues(hypre.ijX, 0.0);
            hypreCheck(HYPRE_BoomerAMGSolve(hypre.precond, hypre.parMatrix, hypre.parB, hypre.parX), "HYPRE_BoomerAMGSolve repeat");
        }
        row.timings.applyMeanSeconds = secondsSince(applyStart) / double(options.applyRepeat);

        HYPRE_IJVectorSetConstantValues(hypre.ijX, 0.0);
        configureSolver(hypre, options);
        const auto solverSetupStart = std::chrono::steady_clock::now();
        setupKrylov(hypre, options);
        row.timings.solverSetupSeconds = secondsSince(solverSetupStart);
        const auto solveStart = std::chrono::steady_clock::now();
        solveKrylov(hypre, options);
        row.timings.solveSeconds = secondsSince(solveStart);
        getKrylovStatus(hypre, options, row.iterations, row.hypreResidual);

        const Eigen::VectorXd activeX = readHypreVector(hypre.ijX, variant.matrix.rows);
        row.trueResidual = (variant.rhs - spmv(variant.matrix, activeX)).norm() / std::max(variant.rhs.norm(), 1.0e-300);
        row.solutionError = (activeX - variant.referenceInActiveOrder).norm() / std::max(variant.referenceInActiveOrder.norm(), 1.0e-300);
    } catch (const std::exception &e) {
        row.success = false;
        row.message = e.what();
    }
    return row;
}

void printCsvHeader() {
    std::cout
        << "success,message,representation,nullspace,solver,threads,omp_max_threads,omp_num_threads,omp_dynamic,omp_proc_bind,omp_places,nodal,interp_type,strong_threshold,relax_type,num_sweeps,"
        << "non_galerkin_tol,interp_vec_variant,interp_vec_qmax,interp_vec_abs_qtrunc,nodal_diag,keep_same_sign,dof_mode,rows,cols,nnz,full_rows,block_size,near_nullspace_columns,"
        << "used_interp_vectors,symmetry_probe_relative,matrix_seconds,vector_seconds,preconditioner_setup_seconds,"
        << "solver_setup_seconds,apply_mean_seconds,solve_seconds,iterations,hypre_residual,true_residual,solution_error\n";
}

void printCsvRow(const ResultRow &row) {
    auto quote = [](std::string value) {
        std::replace(value.begin(), value.end(), ',', ';');
        return value;
    };
    const Options &o = row.options;
    std::cout
        << (row.success ? 1 : 0) << ","
        << quote(row.message) << ","
        << o.representation << ","
        << o.nullspace << ","
        << o.solver << ","
        << o.threads << ","
        << row.ompMaxThreads << ","
        << quote(row.ompNumThreads) << ","
        << quote(row.ompDynamic) << ","
        << quote(row.ompProcBind) << ","
        << quote(row.ompPlaces) << ","
        << o.nodal << ","
        << o.interpType << ","
        << std::setprecision(12) << o.strongThreshold << ","
        << o.relaxType << ","
        << o.numSweeps << ","
        << o.nonGalerkinTol << ","
        << o.interpVecVariant << ","
        << o.interpVecQMax << ","
        << o.interpVecAbsQTrunc << ","
        << o.nodalDiag << ","
        << o.keepSameSign << ","
        << row.dofMode << ","
        << row.rows << ","
        << row.cols << ","
        << row.nnz << ","
        << row.fullRows << ","
        << row.blockSize << ","
        << row.nearNullspaceColumns << ","
        << row.usedInterpVectors << ","
        << row.symmetryProbeRelative << ","
        << row.timings.matrixSeconds << ","
        << row.timings.vectorSeconds << ","
        << row.timings.preconditionerSetupSeconds << ","
        << row.timings.solverSetupSeconds << ","
        << row.timings.applyMeanSeconds << ","
        << row.timings.solveSeconds << ","
        << row.iterations << ","
        << row.hypreResidual << ","
        << row.trueResidual << ","
        << row.solutionError
        << "\n";
}

} // namespace

int main(int argc, char **argv) {
    try {
        Options options = parseOptions(argc, argv);
        ReplayData replay = loadReplay(options.replayDir);
        ResultRow row = runBenchmark(options, replay);
        if (options.csv) {
            printCsvHeader();
            printCsvRow(row);
        } else {
            std::cout << "success=" << row.success << "\n"
                      << "message=" << row.message << "\n"
                      << "representation=" << options.representation << "\n"
                      << "nullspace=" << options.nullspace << "\n"
                      << "threads=" << options.threads << "\n"
                      << "rows=" << row.rows << "\n"
                      << "nnz=" << row.nnz << "\n"
                      << "setup_seconds=" << row.timings.preconditionerSetupSeconds << "\n"
                      << "apply_mean_seconds=" << row.timings.applyMeanSeconds << "\n"
                      << "solve_seconds=" << row.timings.solveSeconds << "\n"
                      << "iterations=" << row.iterations << "\n"
                      << "true_residual=" << row.trueResidual << "\n"
                      << "solution_error=" << row.solutionError << "\n";
        }
        HYPRE_Finalize();
        int finalized = 0;
        MPI_Finalized(&finalized);
        if (!finalized) MPI_Finalize();
        return row.success ? 0 : 2;
    } catch (const std::exception &e) {
        std::cerr << "benchmark-hypre-replay error: " << e.what() << "\n";
        return 1;
    }
}
