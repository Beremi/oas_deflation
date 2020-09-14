#include "model.h"

using namespace std;

int main(int argc, char **argv) {
    if ( argc == 1 ) {
        cerr << "Error: no input file specified, please try again with input file" << endl;
        cerr << "Usage: DiscreteModel [path/to/master.inp]" << endl;
        //cerr << std::setfill('=') << setw(50) << "" << endl;
        cerr << string(80, '=') << endl;
        //cerr << version_info() << endl;
        exit(EXIT_FAILURE);
    }
    
    // OMP set 1 thread by default
    omp_set_dynamic(0);
    char *val = getenv("OMP_NUM_THREADS");
    if (val != nullptr)
        cout << "OMP_NUM_THREADS = " << val << endl;
    else
        omp_set_num_threads(1);

    //initial time
    auto start = std :: chrono :: system_clock :: now();
    auto now = start;
    std :: chrono :: duration< double >elapsed_seconds;
    if ( PRINT_TIME ) {
        std :: time_t time_now = std :: chrono :: system_clock :: to_time_t(start);
        string nowstring = ctime(& time_now);
        std :: cout << "######### start of calculation on: " << nowstring.substr(0, nowstring.length() - 1) << " #########" << endl;
    }

    Model model(PRINT_TIME);
    model.readFromFile(argv[1]);
    model.init();
    model.solve();


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

    return 0;
}
