#include "model.h"
#include "version.h"
#include <getopt.h>

using namespace std;
fs :: path GlobPaths :: BASEDIR;

Model *masterModel;

int main(int argc, char **argv) {
    int numfnd, opt;
    int num;

    num = 0;
    numfnd = 0;
    while ( ( opt = getopt(argc, argv, "j:") ) != -1 ) {
        switch ( opt ) {
        case 'j':
            num = atoi(optarg);
            numfnd = 1;
            break;
        default: /* '?' */
            fprintf(stderr, "Usage: %s [-j num] path/to/master.inp\n",
                    argv [ 0 ]);
            exit(EXIT_FAILURE);
        }
    }

    if ( optind >= argc ) {
        fprintf(stderr, "Expected argument after options\n");
        fprintf(stderr, "Usage: %s [-j num] path/to/master.inp\n",
                argv [ 0 ]);
        fprintf(stderr, "     : [-j num] has effect for Eigen conjugate gradients\n");
        exit(EXIT_FAILURE);
    }

    // OMP set 1 thread by default
    omp_set_dynamic(0);
    char *val = getenv("OMP_NUM_THREADS");
    if ( num != 0 ) {
        omp_set_num_threads(num);
        cout << "Number of threads = " << num << endl;
    } else if ( val != nullptr )  {
        cout << "Number of threads = " << val << endl;
    } else {
        omp_set_num_threads(1);
        cout << "Number of threads = 1" << endl;
    }

    fs :: path input = fs :: absolute(argv [ optind ]);
    GlobPaths :: BASEDIR = input.parent_path();

    //initial time
    auto start = std :: chrono :: system_clock :: now();
    auto now = start;
    std :: chrono :: duration< double >elapsed_seconds;
    if ( PRINT_TIME ) {
        std :: time_t time_now = std :: chrono :: system_clock :: to_time_t(start);
        string nowstring = ctime(& time_now);
        std :: cout << "######### start of calculation on: " << nowstring.substr(0, nowstring.length() - 1) << " #########" << endl;
    }


    masterModel = new Model(PRINT_TIME);
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
        elapsed_seconds = now - start;
        std :: time_t time_now = std :: chrono :: system_clock :: to_time_t(now);
        string nowstring = ctime(& time_now);
        std :: cout << "######### end of calculation on: " << nowstring.substr(0, nowstring.length() - 1) << " #########" << endl;
        std :: cout << "######### total duration: " << convertTimeToString(elapsed_seconds) << " #########" << endl;
    }

    //int terminationStatus = solver->giveTerminationStatus();
    //std :: cout << "termination status = " << terminationStatus << '\n';
    //return terminationStatus;

    delete masterModel; //delete mater model class

    return 0;
}
