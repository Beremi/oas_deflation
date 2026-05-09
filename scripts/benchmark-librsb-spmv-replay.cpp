#include <Eigen/Core>
#include <Eigen/SparseCore>

#include <rsb.h>

#ifdef _OPENMP
#include <omp.h>
#endif

#include <algorithm>
#include <cerrno>
#include <chrono>
#include <cctype>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <iomanip>
#include <iostream>
#include <limits>
#include <memory>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

using SparseRowMatrix = Eigen::SparseMatrix<double, Eigen::RowMajor, int>;
using Triplet = Eigen::Triplet<double>;

namespace {

struct Options {
    std::string replayDir;
    std::string backend = "all";
    std::string vectorName = "solution";
    std::vector<int> threads = {1};
    std::vector<std::string> rsbFlags = {"noflags"};
    std::vector<double> subdivisions = {1.0};
    int warmup = 3;
    int repeat = 20;
    double validateTol = 1.0e-10;
    bool csv = false;
    bool tune = false;
    int tuneRounds = 2;
    double tuneSeconds = 2.0;
};

struct MatrixInfo {
    int rows = 0;
    int cols = 0;
    long long inputNnz = 0;
};

struct RsbInfo {
    long requestedThreads = 0;
    long actualThreads = -1;
    long rows = 0;
    long cols = 0;
    long nnz = 0;
    double bytes = 0.0;
    long leaves = 0;
    unsigned long flags = 0;
    double indexBytesPerNnz = 0.0;
};

struct BenchmarkResult {
    std::string backend;
    int threads = 1;
    int ompMaxThreads = 1;
    std::string rsbFlagsName;
    std::string vectorSource;
    double subdivision = 1.0;
    bool tune = false;
    double tuneSpeedup = 0.0;
    long tuneBestThreads = 0;
    int warmup = 0;
    int repeat = 0;
    int rows = 0;
    int cols = 0;
    long long nnz = 0;
    RsbInfo rsb;
    double secondsTotal = 0.0;
    double secondsPerSpmv = 0.0;
    double gflops = 0.0;
    double gibPerSecond = 0.0;
    double relativeError = 0.0;
    double checksum = 0.0;
    bool success = true;
    std::string message;
};

struct RsbDeleter {
    void operator()(rsb_mtx_t *matrix) const {
        if (matrix) {
            rsb_mtx_free(matrix);
        }
    }
};

using RsbMatrixPtr = std::unique_ptr<rsb_mtx_t, RsbDeleter>;

std::string getEnv(const char *name) {
    const char *value = std::getenv(name);
    return value ? std::string(value) : std::string();
}

std::vector<std::string> split(const std::string &text, char delimiter) {
    std::vector<std::string> values;
    std::stringstream stream(text);
    std::string item;
    while (std::getline(stream, item, delimiter)) {
        item.erase(item.begin(), std::find_if(item.begin(), item.end(), [](unsigned char c) {
            return !std::isspace(c);
        }));
        item.erase(std::find_if(item.rbegin(), item.rend(), [](unsigned char c) {
            return !std::isspace(c);
        }).base(), item.end());
        if (!item.empty()) {
            values.push_back(item);
        }
    }
    return values;
}

std::vector<int> parseIntList(const std::string &text) {
    std::vector<int> values;
    for (const std::string &item : split(text, ',')) {
        values.push_back(std::stoi(item));
    }
    return values;
}

std::vector<double> parseDoubleList(const std::string &text) {
    std::vector<double> values;
    for (const std::string &item : split(text, ',')) {
        values.push_back(std::stod(item));
    }
    return values;
}

void printUsage(const char *argv0) {
    std::cerr
        << "Usage: " << argv0 << " --replay-dir DIR [options]\n"
        << "Options:\n"
        << "  --backend eigen|manual|librsb|all\n"
        << "  --vector solution|rhs|ones\n"
        << "  --threads 1,2,4,8,16\n"
        << "  --rsb-flags noflags,default\n"
        << "  --subdivisions 0.5,1,2,4\n"
        << "  --warmup N\n"
        << "  --repeat N\n"
        << "  --validate-tol VALUE\n"
        << "  --tune\n"
        << "  --tune-rounds N\n"
        << "  --tune-seconds VALUE\n"
        << "  --csv\n";
}

Options parseOptions(int argc, char **argv) {
    Options options;
    for (int i = 1; i < argc; i++) {
        const std::string arg = argv[i];
        auto requireValue = [&](const char *name) -> std::string {
            if (i + 1 >= argc) {
                throw std::runtime_error(std::string("missing value for ") + name);
            }
            return argv[++i];
        };

        if (arg == "--replay-dir") {
            options.replayDir = requireValue("--replay-dir");
        } else if (arg == "--backend") {
            options.backend = requireValue("--backend");
        } else if (arg == "--vector") {
            options.vectorName = requireValue("--vector");
        } else if (arg == "--threads") {
            options.threads = parseIntList(requireValue("--threads"));
        } else if (arg == "--rsb-flags") {
            options.rsbFlags = split(requireValue("--rsb-flags"), ',');
        } else if (arg == "--subdivisions") {
            options.subdivisions = parseDoubleList(requireValue("--subdivisions"));
        } else if (arg == "--warmup") {
            options.warmup = std::stoi(requireValue("--warmup"));
        } else if (arg == "--repeat") {
            options.repeat = std::stoi(requireValue("--repeat"));
        } else if (arg == "--validate-tol") {
            options.validateTol = std::stod(requireValue("--validate-tol"));
        } else if (arg == "--tune") {
            options.tune = true;
        } else if (arg == "--tune-rounds") {
            options.tuneRounds = std::stoi(requireValue("--tune-rounds"));
        } else if (arg == "--tune-seconds") {
            options.tuneSeconds = std::stod(requireValue("--tune-seconds"));
        } else if (arg == "--csv") {
            options.csv = true;
        } else if (arg == "--help" || arg == "-h") {
            printUsage(argv[0]);
            std::exit(0);
        } else {
            throw std::runtime_error("unknown argument: " + arg);
        }
    }

    if (options.replayDir.empty()) {
        throw std::runtime_error("--replay-dir is required");
    }
    if (options.warmup < 0 || options.repeat <= 0) {
        throw std::runtime_error("--warmup must be nonnegative and --repeat must be positive");
    }
    return options;
}

std::string pathJoin(const std::string &dir, const std::string &name) {
    if (dir.empty() || dir.back() == '/') {
        return dir + name;
    }
    return dir + "/" + name;
}

bool parseMatrixMarketSize(FILE *file, int &rows, int &cols, long long &nnz) {
    char line[4096];
    while (std::fgets(line, sizeof(line), file)) {
        if (line[0] == '%') {
            continue;
        }
        char *cursor = line;
        errno = 0;
        const long parsedRows = std::strtol(cursor, &cursor, 10);
        const long parsedCols = std::strtol(cursor, &cursor, 10);
        const long long parsedNnz = std::strtoll(cursor, &cursor, 10);
        if (errno != 0 || parsedRows <= 0 || parsedCols <= 0 || parsedNnz < 0) {
            return false;
        }
        if (parsedRows > std::numeric_limits<int>::max()
                || parsedCols > std::numeric_limits<int>::max()) {
            throw std::runtime_error("matrix dimensions exceed Eigen int index range");
        }
        rows = static_cast<int>(parsedRows);
        cols = static_cast<int>(parsedCols);
        nnz = parsedNnz;
        return true;
    }
    return false;
}

SparseRowMatrix loadMatrixMarket(const std::string &path, MatrixInfo &info) {
    FILE *file = std::fopen(path.c_str(), "r");
    if (!file) {
        throw std::runtime_error("cannot open matrix: " + path);
    }

    int rows = 0;
    int cols = 0;
    long long inputNnz = 0;
    if (!parseMatrixMarketSize(file, rows, cols, inputNnz)) {
        std::fclose(file);
        throw std::runtime_error("cannot parse MatrixMarket size from " + path);
    }

    std::vector<Triplet> triplets;
    const size_t reserveNnz = static_cast<size_t>(inputNnz);
    if (static_cast<long long>(reserveNnz) != inputNnz) {
        std::fclose(file);
        throw std::runtime_error("matrix nnz exceeds size_t range");
    }
    triplets.reserve(reserveNnz);

    char line[4096];
    long long parsed = 0;
    while (std::fgets(line, sizeof(line), file)) {
        if (line[0] == '%' || line[0] == '\n' || line[0] == '\0') {
            continue;
        }
        char *cursor = line;
        errno = 0;
        const long row = std::strtol(cursor, &cursor, 10);
        const long col = std::strtol(cursor, &cursor, 10);
        const double value = std::strtod(cursor, &cursor);
        if (errno != 0 || row <= 0 || col <= 0 || row > rows || col > cols) {
            std::fclose(file);
            throw std::runtime_error("bad MatrixMarket entry in " + path);
        }
        triplets.emplace_back(static_cast<int>(row - 1), static_cast<int>(col - 1), value);
        parsed++;
    }
    std::fclose(file);

    if (parsed != inputNnz) {
        std::ostringstream message;
        message << "MatrixMarket nnz mismatch: header=" << inputNnz << " parsed=" << parsed;
        throw std::runtime_error(message.str());
    }

    SparseRowMatrix matrix(rows, cols);
    matrix.setFromTriplets(triplets.begin(), triplets.end());
    matrix.makeCompressed();

    info.rows = rows;
    info.cols = cols;
    info.inputNnz = inputNnz;
    return matrix;
}

bool fileExists(const std::string &path) {
    FILE *file = std::fopen(path.c_str(), "r");
    if (!file) {
        return false;
    }
    std::fclose(file);
    return true;
}

Eigen::VectorXd loadVectorFile(const std::string &path, int size) {
    FILE *file = std::fopen(path.c_str(), "r");
    if (!file) {
        throw std::runtime_error("cannot open vector: " + path);
    }

    Eigen::VectorXd vector = Eigen::VectorXd::Zero(size);
    char line[4096];
    while (std::fgets(line, sizeof(line), file)) {
        if (std::strncmp(line, "index", 5) == 0 || line[0] == '\n' || line[0] == '\0') {
            continue;
        }
        char *cursor = line;
        errno = 0;
        const long index = std::strtol(cursor, &cursor, 10);
        const double value = std::strtod(cursor, &cursor);
        if (errno != 0 || index < 0 || index >= size) {
            std::fclose(file);
            throw std::runtime_error("bad vector entry in " + path);
        }
        vector[static_cast<int>(index)] = value;
    }
    std::fclose(file);
    return vector;
}

Eigen::VectorXd loadInputVector(const Options &options, int size, std::string &source) {
    if (options.vectorName == "ones") {
        source = "ones";
        return Eigen::VectorXd::Ones(size);
    }

    const std::string requested = pathJoin(options.replayDir, options.vectorName + ".tsv");
    if (fileExists(requested)) {
        source = options.vectorName;
        return loadVectorFile(requested, size);
    }

    const std::string solution = pathJoin(options.replayDir, "solution.tsv");
    if (fileExists(solution)) {
        source = "solution";
        return loadVectorFile(solution, size);
    }

    const std::string rhs = pathJoin(options.replayDir, "rhs.tsv");
    if (fileExists(rhs)) {
        source = "rhs";
        return loadVectorFile(rhs, size);
    }

    source = "ones";
    return Eigen::VectorXd::Ones(size);
}

void setThreadCount(int threads) {
#ifdef _OPENMP
    omp_set_num_threads(threads);
#else
    (void)threads;
#endif
    Eigen::setNbThreads(threads);
}

int ompMaxThreads() {
#ifdef _OPENMP
    return omp_get_max_threads();
#else
    return 1;
#endif
}

double estimatedSpmvBytes(long long nnz, int rows) {
    const double perNnz = static_cast<double>(sizeof(double) + sizeof(int) + sizeof(double));
    const double perRow = static_cast<double>(sizeof(double) + sizeof(int));
    return static_cast<double>(nnz) * perNnz + static_cast<double>(rows) * perRow;
}

double relativeError(const Eigen::VectorXd &reference, const Eigen::VectorXd &candidate) {
    const double denominator = std::max(reference.norm(), 1.0e-300);
    return (candidate - reference).norm() / denominator;
}

template <typename Function>
double timeRepeated(const Options &options, Function function) {
    for (int i = 0; i < options.warmup; i++) {
        function();
    }

    const auto start = std::chrono::steady_clock::now();
    for (int i = 0; i < options.repeat; i++) {
        function();
    }
    const auto end = std::chrono::steady_clock::now();
    return std::chrono::duration<double>(end - start).count();
}

void fillCommonResult(
    BenchmarkResult &result,
    const Options &options,
    const SparseRowMatrix &matrix,
    int threads,
    const std::string &vectorSource,
    double totalSeconds,
    const Eigen::VectorXd &reference,
    const Eigen::VectorXd &candidate
) {
    result.threads = threads;
    result.ompMaxThreads = ompMaxThreads();
    result.vectorSource = vectorSource;
    result.warmup = options.warmup;
    result.repeat = options.repeat;
    result.rows = static_cast<int>(matrix.rows());
    result.cols = static_cast<int>(matrix.cols());
    result.nnz = static_cast<long long>(matrix.nonZeros());
    result.secondsTotal = totalSeconds;
    result.secondsPerSpmv = totalSeconds / static_cast<double>(options.repeat);
    result.gflops = 2.0 * static_cast<double>(matrix.nonZeros()) / result.secondsPerSpmv / 1.0e9;
    result.gibPerSecond = estimatedSpmvBytes(result.nnz, result.rows)
        / result.secondsPerSpmv / 1024.0 / 1024.0 / 1024.0;
    result.relativeError = relativeError(reference, candidate);
    result.checksum = candidate.sum();
    result.success = result.relativeError <= options.validateTol;
    if (!result.success) {
        std::ostringstream message;
        message << "relative error " << result.relativeError
                << " exceeds tolerance " << options.validateTol;
        result.message = message.str();
    }
}

BenchmarkResult benchmarkEigen(
    const Options &options,
    const SparseRowMatrix &matrix,
    const Eigen::VectorXd &x,
    const Eigen::VectorXd &reference,
    int threads,
    const std::string &vectorSource
) {
    Eigen::VectorXd y(matrix.rows());
    const double totalSeconds = timeRepeated(options, [&]() {
        y.noalias() = matrix * x;
    });

    BenchmarkResult result;
    result.backend = "eigen-rowmajor";
    fillCommonResult(result, options, matrix, threads, vectorSource, totalSeconds, reference, y);
    return result;
}

BenchmarkResult benchmarkManual(
    const Options &options,
    const SparseRowMatrix &matrix,
    const Eigen::VectorXd &x,
    const Eigen::VectorXd &reference,
    int threads,
    const std::string &vectorSource
) {
    Eigen::VectorXd y(matrix.rows());
    const int *rowPtr = matrix.outerIndexPtr();
    const int *colInd = matrix.innerIndexPtr();
    const double *values = matrix.valuePtr();
    const int rows = static_cast<int>(matrix.rows());

    const double totalSeconds = timeRepeated(options, [&]() {
#pragma omp parallel for schedule(static)
        for (int row = 0; row < rows; row++) {
            double sum = 0.0;
            for (int idx = rowPtr[row]; idx < rowPtr[row + 1]; idx++) {
                sum += values[idx] * x[colInd[idx]];
            }
            y[row] = sum;
        }
    });

    BenchmarkResult result;
    result.backend = "manual-crs-openmp";
    fillCommonResult(result, options, matrix, threads, vectorSource, totalSeconds, reference, y);
    return result;
}

bool ensureRsbInitialized() {
    static bool initialized = false;
    if (initialized) {
        return true;
    }
    const rsb_err_t error = rsb_lib_init(RSB_NULL_INIT_OPTIONS);
    initialized = (error == RSB_ERR_NO_ERROR);
    return initialized;
}

rsb_flags_t parseRsbFlags(const std::string &name) {
    if (name == "noflags") {
        return RSB_FLAG_NOFLAGS;
    }
    if (name == "default") {
        return RSB_FLAG_DEFAULT_RSB_MATRIX_FLAGS;
    }
    std::ostringstream message;
    message << "unsupported --rsb-flags value: " << name;
    throw std::runtime_error(message.str());
}

bool queryRsbInfo(const rsb_mtx_t *matrix, RsbInfo &info) {
    if (!matrix) {
        return false;
    }

    rsb_int_t actualThreads = -1;
    rsb_err_t error = rsb_lib_get_opt(RSB_IO_WANT_EXECUTING_THREADS, &actualThreads);
    if (error == RSB_ERR_NO_ERROR) {
        info.actualThreads = static_cast<long>(actualThreads);
    }

    rsb_coo_idx_t rows = 0;
    rsb_coo_idx_t cols = 0;
    rsb_nnz_idx_t nnz = 0;
    size_t bytes = 0;
    rsb_blk_idx_t leaves = 0;
    rsb_flags_t flags = RSB_FLAG_NOFLAGS;
    rsb_real_t indexBytesPerNnz = 0.0;

    error = RSB_ERR_NO_ERROR;
    error |= rsb_mtx_get_info(matrix, RSB_MIF_MATRIX_ROWS__TO__RSB_COO_INDEX_T, &rows);
    error |= rsb_mtx_get_info(matrix, RSB_MIF_MATRIX_COLS__TO__RSB_COO_INDEX_T, &cols);
    error |= rsb_mtx_get_info(matrix, RSB_MIF_MATRIX_NNZ__TO__RSB_NNZ_INDEX_T, &nnz);
    error |= rsb_mtx_get_info(matrix, RSB_MIF_TOTAL_SIZE__TO__SIZE_T, &bytes);
    error |= rsb_mtx_get_info(matrix, RSB_MIF_LEAVES_COUNT__TO__RSB_BLK_INDEX_T, &leaves);
    error |= rsb_mtx_get_info(matrix, RSB_MIF_MATRIX_FLAGS__TO__RSB_FLAGS_T, &flags);
    error |= rsb_mtx_get_info(
        matrix,
        RSB_MIF_INDEX_STORAGE_IN_BYTES_PER_NNZ__TO__RSB_REAL_T,
        &indexBytesPerNnz
    );

    if (error != RSB_ERR_NO_ERROR) {
        return false;
    }

    info.rows = static_cast<long>(rows);
    info.cols = static_cast<long>(cols);
    info.nnz = static_cast<long>(nnz);
    info.bytes = static_cast<double>(bytes);
    info.leaves = static_cast<long>(leaves);
    info.flags = static_cast<unsigned long>(flags);
    info.indexBytesPerNnz = static_cast<double>(indexBytesPerNnz);
    return true;
}

RsbMatrixPtr buildRsbMatrix(
    const SparseRowMatrix &matrix,
    int threads,
    const std::string &flagsName,
    double subdivision,
    RsbInfo &info
) {
    if (!ensureRsbInitialized()) {
        throw std::runtime_error("rsb_lib_init failed");
    }

    rsb_int_t requestedThreads = static_cast<rsb_int_t>(threads);
    rsb_err_t error = rsb_lib_set_opt(RSB_IO_WANT_EXECUTING_THREADS, &requestedThreads);
    if (error != RSB_ERR_NO_ERROR) {
        throw std::runtime_error("rsb_lib_set_opt(RSB_IO_WANT_EXECUTING_THREADS) failed");
    }

    rsb_real_t requestedSubdivision = static_cast<rsb_real_t>(subdivision);
    error = rsb_lib_set_opt(RSB_IO_WANT_SUBDIVISION_MULTIPLIER, &requestedSubdivision);
    if (error != RSB_ERR_NO_ERROR) {
        throw std::runtime_error("rsb_lib_set_opt(RSB_IO_WANT_SUBDIVISION_MULTIPLIER) failed");
    }

    if (matrix.rows() > std::numeric_limits<rsb_coo_idx_t>::max()
            || matrix.cols() > std::numeric_limits<rsb_coo_idx_t>::max()
            || matrix.nonZeros() > std::numeric_limits<rsb_nnz_idx_t>::max()) {
        throw std::runtime_error("matrix exceeds librsb index range");
    }

    std::vector<rsb_coo_idx_t> rowPtr(static_cast<size_t>(matrix.rows()) + 1);
    std::vector<rsb_coo_idx_t> colInd(static_cast<size_t>(matrix.nonZeros()));
    for (int i = 0; i <= matrix.rows(); i++) {
        rowPtr[static_cast<size_t>(i)] = static_cast<rsb_coo_idx_t>(matrix.outerIndexPtr()[i]);
    }
    for (int i = 0; i < matrix.nonZeros(); i++) {
        colInd[static_cast<size_t>(i)] = static_cast<rsb_coo_idx_t>(matrix.innerIndexPtr()[i]);
    }

    error = RSB_ERR_NO_ERROR;
    rsb_mtx_t *raw = rsb_mtx_alloc_from_csr_const(
        matrix.valuePtr(),
        rowPtr.data(),
        colInd.data(),
        static_cast<rsb_nnz_idx_t>(matrix.nonZeros()),
        RSB_NUMERICAL_TYPE_DOUBLE,
        static_cast<rsb_coo_idx_t>(matrix.rows()),
        static_cast<rsb_coo_idx_t>(matrix.cols()),
        1,
        1,
        parseRsbFlags(flagsName),
        &error
    );
    if (!raw || error != RSB_ERR_NO_ERROR) {
        throw std::runtime_error("rsb_mtx_alloc_from_csr_const failed");
    }

    RsbMatrixPtr rsb(raw);
    info.requestedThreads = threads;
    queryRsbInfo(rsb.get(), info);
    return rsb;
}

void tuneRsbMatrix(
    RsbMatrixPtr &matrix,
    const Options &options,
    const Eigen::VectorXd &x,
    Eigen::VectorXd &y,
    BenchmarkResult &result
) {
    if (!options.tune) {
        return;
    }

    double alpha = 1.0;
    double beta = 0.0;
    rsb_real_t speedup = 0.0;
    rsb_int_t bestThreads = static_cast<rsb_int_t>(result.threads);
    rsb_mtx_t *raw = matrix.release();

    const rsb_err_t error = rsb_tune_spmm(
        &raw,
        &speedup,
        &bestThreads,
        static_cast<rsb_int_t>(options.tuneRounds),
        static_cast<rsb_time_t>(options.tuneSeconds),
        RSB_TRANSPOSITION_N,
        &alpha,
        nullptr,
        1,
        RSB_FLAG_WANT_COLUMN_MAJOR_ORDER,
        x.data(),
        0,
        &beta,
        y.data(),
        0
    );
    matrix.reset(raw);

    if (error != RSB_ERR_NO_ERROR) {
        result.success = false;
        result.message = "rsb_tune_spmm failed";
        return;
    }

    result.tune = true;
    result.tuneSpeedup = static_cast<double>(speedup);
    result.tuneBestThreads = static_cast<long>(bestThreads);
    queryRsbInfo(matrix.get(), result.rsb);
}

BenchmarkResult benchmarkRsb(
    const Options &options,
    const SparseRowMatrix &matrix,
    const Eigen::VectorXd &x,
    const Eigen::VectorXd &reference,
    int threads,
    const std::string &vectorSource,
    const std::string &flagsName,
    double subdivision
) {
    BenchmarkResult result;
    result.backend = "librsb";
    result.rsbFlagsName = flagsName;
    result.subdivision = subdivision;
    result.threads = threads;
    result.ompMaxThreads = ompMaxThreads();
    result.vectorSource = vectorSource;
    result.warmup = options.warmup;
    result.repeat = options.repeat;
    result.rows = static_cast<int>(matrix.rows());
    result.cols = static_cast<int>(matrix.cols());
    result.nnz = static_cast<long long>(matrix.nonZeros());

    Eigen::VectorXd y(matrix.rows());
    try {
        RsbMatrixPtr rsb = buildRsbMatrix(matrix, threads, flagsName, subdivision, result.rsb);
        tuneRsbMatrix(rsb, options, x, y, result);
        if (!result.success) {
            return result;
        }

        const double alpha = 1.0;
        const double beta = 0.0;
        const double totalSeconds = timeRepeated(options, [&]() {
            const rsb_err_t error = rsb_spmv(
                RSB_TRANSPOSITION_N,
                &alpha,
                rsb.get(),
                x.data(),
                1,
                &beta,
                y.data(),
                1
            );
            if (error != RSB_ERR_NO_ERROR) {
                throw std::runtime_error("rsb_spmv failed");
            }
        });
        fillCommonResult(result, options, matrix, threads, vectorSource, totalSeconds, reference, y);
        result.backend = "librsb";
        result.rsbFlagsName = flagsName;
        result.subdivision = subdivision;
        result.tune = options.tune;
        result.rsb.requestedThreads = threads;
        queryRsbInfo(rsb.get(), result.rsb);
    } catch (const std::exception &error) {
        result.success = false;
        result.message = error.what();
    }
    return result;
}

std::string csvEscape(const std::string &value) {
    if (value.find_first_of(",\"\n") == std::string::npos) {
        return value;
    }
    std::string escaped = "\"";
    for (char c : value) {
        if (c == '"') {
            escaped += '"';
        }
        escaped += c;
    }
    escaped += '"';
    return escaped;
}

void printCsvHeader() {
    std::cout
        << "backend,threads_requested,omp_max_threads,rsb_num_threads_env,omp_proc_bind,omp_places,omp_dynamic,"
        << "rsb_actual_threads,rsb_flags_name,rsb_flags_hex,subdivision,tuned,tune_speedup,tune_best_threads,warmup,repeat,"
        << "vector_source,rows,cols,nnz,rsb_rows,rsb_cols,rsb_nnz,rsb_bytes,rsb_leaves,"
        << "rsb_matrix_flags_hex,rsb_index_bytes_per_nnz,seconds_total,seconds_per_spmv,gflops,gib_per_s,"
        << "relative_error,checksum,success,message\n";
}

void printCsvResult(const BenchmarkResult &result) {
    std::cout
        << result.backend << ','
        << result.threads << ','
        << result.ompMaxThreads << ','
        << csvEscape(getEnv("RSB_NUM_THREADS")) << ','
        << csvEscape(getEnv("OMP_PROC_BIND")) << ','
        << csvEscape(getEnv("OMP_PLACES")) << ','
        << csvEscape(getEnv("OMP_DYNAMIC")) << ','
        << result.rsb.actualThreads << ','
        << csvEscape(result.rsbFlagsName) << ','
        << "0x" << std::hex << result.rsb.flags << std::dec << ','
        << result.subdivision << ','
        << (result.tune ? 1 : 0) << ','
        << result.tuneSpeedup << ','
        << result.tuneBestThreads << ','
        << result.warmup << ','
        << result.repeat << ','
        << csvEscape(result.vectorSource) << ','
        << result.rows << ','
        << result.cols << ','
        << result.nnz << ','
        << result.rsb.rows << ','
        << result.rsb.cols << ','
        << result.rsb.nnz << ','
        << result.rsb.bytes << ','
        << result.rsb.leaves << ','
        << "0x" << std::hex << result.rsb.flags << std::dec << ','
        << result.rsb.indexBytesPerNnz << ','
        << std::setprecision(12) << result.secondsTotal << ','
        << result.secondsPerSpmv << ','
        << result.gflops << ','
        << result.gibPerSecond << ','
        << result.relativeError << ','
        << result.checksum << ','
        << (result.success ? 1 : 0) << ','
        << csvEscape(result.message)
        << '\n';
}

void printHumanHeader(const Options &options, const MatrixInfo &info, const std::string &vectorSource) {
    std::cerr << "replay_dir=" << options.replayDir << '\n';
    std::cerr << "matrix_rows=" << info.rows
              << " matrix_cols=" << info.cols
              << " input_nnz=" << info.inputNnz << '\n';
    std::cerr << "vector_source=" << vectorSource << '\n';
    std::cerr << "RSB_NUM_THREADS=" << getEnv("RSB_NUM_THREADS")
              << " OMP_PROC_BIND=" << getEnv("OMP_PROC_BIND")
              << " OMP_PLACES=" << getEnv("OMP_PLACES")
              << " OMP_DYNAMIC=" << getEnv("OMP_DYNAMIC") << '\n';
}

bool wantsBackend(const Options &options, const std::string &backend) {
    return options.backend == "all" || options.backend == backend;
}

} // namespace

int main(int argc, char **argv) {
    try {
        const Options options = parseOptions(argc, argv);

        MatrixInfo matrixInfo;
        SparseRowMatrix matrix = loadMatrixMarket(pathJoin(options.replayDir, "matrix.mtx"), matrixInfo);
        std::string vectorSource;
        Eigen::VectorXd x = loadInputVector(options, static_cast<int>(matrix.cols()), vectorSource);
        Eigen::VectorXd reference(matrix.rows());
        reference.noalias() = matrix * x;

        printHumanHeader(options, matrixInfo, vectorSource);
        if (options.csv) {
            printCsvHeader();
        }

        bool allSucceeded = true;
        for (int threads : options.threads) {
            setThreadCount(threads);

            if (wantsBackend(options, "eigen")) {
                BenchmarkResult result = benchmarkEigen(options, matrix, x, reference, threads, vectorSource);
                allSucceeded = allSucceeded && result.success;
                printCsvResult(result);
            }

            if (wantsBackend(options, "manual")) {
                BenchmarkResult result = benchmarkManual(options, matrix, x, reference, threads, vectorSource);
                allSucceeded = allSucceeded && result.success;
                printCsvResult(result);
            }

            if (wantsBackend(options, "librsb")) {
                for (const std::string &flagsName : options.rsbFlags) {
                    for (double subdivision : options.subdivisions) {
                        BenchmarkResult result = benchmarkRsb(
                            options,
                            matrix,
                            x,
                            reference,
                            threads,
                            vectorSource,
                            flagsName,
                            subdivision
                        );
                        allSucceeded = allSucceeded && result.success;
                        printCsvResult(result);
                    }
                }
            }
        }

        return allSucceeded ? 0 : 2;
    } catch (const std::exception &error) {
        std::cerr << "error: " << error.what() << '\n';
        printUsage(argv[0]);
        return 1;
    }
}
