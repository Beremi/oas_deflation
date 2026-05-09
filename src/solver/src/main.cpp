#ifdef PYBIND_FOUND
 #include <pybind11/embed.h>
 #include <pybind11/pybind11.h>
#endif

#include "model.h"
#include "version.h"
#include <iostream>

using namespace std;
fs :: path GlobPaths :: BASEDIR;

Model *masterModel;

#ifdef PYBIND_FOUND
namespace py = pybind11;
int test(fs :: path python_script_path);
#endif

int main(int argc, char **argv) {
#ifdef PYBIND_FOUND
    // Initialize the Python interpreter with pybind11
    // This should be done only once, at the start of main()
    // and before any threads are spawned
    try {
        py :: scoped_interpreter guard{};    // init ONCE, here
        // optional: print something to verify we got here
        std :: cout << "Python interpreter initialized (pybind11)\n";

        // Print Python version and executable path
        py :: module sys = py :: module :: import("sys");
        std :: cout << "Python version: " << py :: str( sys.attr("version") ).cast< std :: string >() << std :: endl;
        std :: cout << "Python executable: " << py :: str( sys.attr("executable") ).cast< std :: string >() << std :: endl;

        // Additional diagnostics to identify the mismatch
        std :: cout << "Python implementation: " << py :: str( sys.attr("implementation").attr("name") ).cast< std :: string >() << std :: endl;

        // Check if running in virtual/conda environment
        auto venv = sys.attr("prefix");
        auto base = sys.attr("base_prefix");
        bool in_venv = !venv.equal(base);
        std :: cout << "Running in virtual/conda environment: " << ( in_venv ? "Yes" : "No" ) << std :: endl;
        std :: cout << "Python prefix: " << py :: str( sys.attr("prefix") ).cast< std :: string >() << std :: endl;
        std :: cout << "Python base_prefix: " << py :: str( sys.attr("base_prefix") ).cast< std :: string >() << std :: endl;
        std :: cout << "Python sys.path: " << py :: str( sys.attr("path") ).cast< std :: string >() << std :: endl;

        // Check which libpython is actually loaded
        py :: module sysconfig = py :: module :: import("sysconfig");
        std :: cout << "Python config prefix: " << py :: str(sysconfig.attr("get_config_var")("prefix") ).cast< std :: string > ( ) << std :: endl;
                                                             std :: cout << "Python LIBDIR: " << py :: str(sysconfig.attr("get_config_var")("LIBDIR") ).cast< std :: string > () << std :: endl;

                                                                                                           // If you use threads anywhere: ensure the main thread owns the GIL initially
                                                                                                           py :: gil_scoped_acquire gil;
    } catch ( const py :: error_already_set &e ) {
        std :: cerr << "Python error:\n" << e.what() << std :: endl;
        return 1;
    } catch ( const std :: exception &e ) {
        std :: cerr << "C++ exception:\n" << e.what() << std :: endl;
        return 1;
    }
#endif

    vector< string >args(argv + 1, argv + argc);
    int num = 0;
    fs :: path master_filename;

    if ( argc == 1 ) {
        fprintf(stdout, "Expected argument after options\n");
        fprintf(stdout, "Usage: %s [parameters: -x value] path/to/master.inp\n",
                argv [ 0 ]);
        fprintf(stdout, "     : [-j num] has effect for Eigen conjugate gradients\n");
#ifdef PYBIND_FOUND
        fprintf(stdout, "     : [-f script.py] path to Python script\n");
#endif
        cerr << string(80, '=') << endl;
        cerr << version_info() << endl;
        exit(EXIT_SUCCESS);
    }

#ifdef PYBIND_FOUND
    fs :: path python_script_path;
#endif

    for ( size_t i = 0; i != args.size(); ++i ) {
        if ( args [ i ] == "-h" || args [ i ] == "--help" ) {
            fprintf(stdout, "Usage: %s [-j num] path/to/master.inp\n", argv [ 0 ]);
            exit(EXIT_SUCCESS);
        } else if ( args [ i ] == "-j" ) {
            if ( i + 1 < args.size() ) {
                num = atoi( ( args [ i + 1 ] ).c_str() );
            } else {
                fprintf(stderr, "Usage: %s [-j num] path/to/master.inp\n", argv [ 0 ]);
                exit(EXIT_FAILURE);
            };
        }
#ifdef PYBIND_FOUND
        else if ( args [ i ] == "-f" ) {
            if ( i + 1 < args.size() ) {
                python_script_path = fs :: path(args [ i + 1 ]);
            } else {
                fprintf(stderr, "Usage: %s [-f script.py] path/to/master.inp\n", argv [ 0 ]);
                exit(EXIT_FAILURE);
            };
        }
#endif
    }
    master_filename = args.back();

#ifdef PYBIND_FOUND
    if ( !python_script_path.empty() ) {
        test(python_script_path);
    }
#endif

    // OMP set 1 thread by default
    omp_set_dynamic(0);
    char *val = getenv("OMP_NUM_THREADS");
    if ( num != 0 ) {
        omp_set_num_threads(num);
        cout << "Number of threads = " << num << endl;
    } else if ( val != nullptr ) {
        cout << "Number of threads = " << val << endl;
    } else {
        omp_set_num_threads(1);
        cout << "Number of threads = 1" << endl;
    }

    fs :: path input = fs :: absolute(master_filename);
    if ( !fs :: exists(input) ) {
        fprintf(stderr, "The problem with input file: %s. (Does not exist, wrong path) \n", input.string().c_str() );
        exit(EXIT_FAILURE);
    }
    GlobPaths :: BASEDIR = input.parent_path();

    //initial time
    std :: chrono :: time_point< std :: chrono :: system_clock >now;
    std :: chrono :: duration< double >elapsed_seconds;
    masterModel = new Model(PRINT_TIME);
    if ( PRINT_TIME ) {
        std :: time_t time_now = std :: chrono :: system_clock :: to_time_t(masterModel->giveStartTime() );
        string nowstring = ctime(& time_now);
        std :: cout << "######### start of calculation on: " << nowstring.substr(0, nowstring.length() - 1) << " #########" << endl;
    }

    auto phaseStart = std :: chrono :: steady_clock :: now();
    masterModel->readFromFile( input.string() );
    double readDuration = std :: chrono :: duration< double >(std :: chrono :: steady_clock :: now() - phaseStart).count();

    // check if exists or create directory for results
    phaseStart = std :: chrono :: steady_clock :: now();
    if ( !fs :: exists(masterModel->resultDir) ) {
        fs :: create_directories(masterModel->resultDir);
    }
    // save version information to result dir
    ofstream outputfile( ( masterModel->resultDir / "version.txt" ).string() );
    if ( outputfile.is_open() ) {
        outputfile << std :: scientific;
        outputfile << version_info(true);
        outputfile << endl;
    }
    outputfile.close();
    double resultSetupDuration = std :: chrono :: duration< double >(std :: chrono :: steady_clock :: now() - phaseStart).count();

    if ( masterModel->giveSolver() ) {
        masterModel->giveSolver()->initializeRuntimeProfilerForExternalTiming();
        masterModel->giveSolver()->recordExternalRuntimePhase("model.read_input", readDuration);
        masterModel->giveSolver()->recordExternalRuntimePhase("model.result_setup", resultSetupDuration);
    }

    phaseStart = std :: chrono :: steady_clock :: now();
    masterModel->init();
    if ( masterModel->giveSolver() ) {
        masterModel->giveSolver()->recordExternalRuntimePhase("model.init_total", std :: chrono :: duration< double >(std :: chrono :: steady_clock :: now() - phaseStart).count() );
    }
    phaseStart = std :: chrono :: steady_clock :: now();
    masterModel->solve();
    if ( masterModel->giveSolver() ) {
        masterModel->giveSolver()->recordExternalRuntimePhase("model.solve_total", std :: chrono :: duration< double >(std :: chrono :: steady_clock :: now() - phaseStart).count() );
    }


    if ( PRINT_TIME ) {
        now = std :: chrono :: system_clock :: now();
        elapsed_seconds = now - masterModel->giveStartTime();
        std :: time_t time_now = std :: chrono :: system_clock :: to_time_t(now);
        string nowstring = ctime(& time_now);
        std :: cout << "######### end of calculation on: " << nowstring.substr(0, nowstring.length() - 1) << " #########" << endl;
        std :: cout << "######### total duration: " << convertTimeToString(elapsed_seconds) << " #########" << endl;
    }

    // when run automatically, termination status needs to be distinguished
    int terminationStatus = masterModel->giveSolver()->giveTerminationStatus();
    //std :: cout << "termination status = " << terminationStatus << '\n';
    //return terminationStatus;

    delete masterModel; //delete mater model class

    return terminationStatus;
}
