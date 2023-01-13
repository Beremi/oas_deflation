#include "model.h"
#include "version.h"

using namespace std;
fs :: path GlobPaths :: BASEDIR;

Model *masterModel;

int main(int argc, char **argv) {
    vector< string >args(argv + 1, argv + argc);
    int num = 0;
    string master_filename;

    if ( argc == 1 ) {
        fprintf(stdout, "Expected argument after options\n");
        fprintf(stdout, "Usage: %s [-j num] path/to/master.inp\n",
                argv [ 0 ]);
        fprintf(stdout, "     : [-j num] has effect for Eigen conjugate gradients\n");
        cerr << string(80, '=') << endl;
        cerr << version_info() << endl;
        exit(EXIT_SUCCESS);
    }

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
    }
    master_filename = args.back();

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
        fprintf(stderr, "The problem with input file: %s. (Does not exist, wrong path) \n", input.c_str() );
        exit(EXIT_FAILURE);
    }
    GlobPaths :: BASEDIR = input.parent_path();

    //initial time
    std :: chrono :: time_point< std :: chrono :: system_clock >now;
    std :: chrono :: duration< double >elapsed_seconds;
    masterModel = new Model(PRINT_TIME);
    if ( PRINT_TIME ) {
        std :: time_t time_now = std :: chrono :: system_clock :: to_time_t( masterModel->giveStartTime() );
        string nowstring = ctime(& time_now);
        std :: cout << "######### start of calculation on: " << nowstring.substr(0, nowstring.length() - 1) << " #########" << endl;
    }



    masterModel->readFromFile(input.string() );

    // check if exists or create directory for results
    if ( !fs :: exists(masterModel->resultDir) ) {
        fs :: create_directories(masterModel->resultDir);
    }
    // save version information to result dir
    ofstream outputfile( ( masterModel->resultDir / "version.txt" ).string() );
    if ( outputfile.is_open() ) {
        outputfile << std :: scientific;
        outputfile << version_info();
        outputfile << endl;
    }
    outputfile.close();


    masterModel->init();
    masterModel->solve();


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
