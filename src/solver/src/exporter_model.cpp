#include "exporter_model.h"
#include "model.h"
#include "solver_implicit.h"
#include <cstdio>

using namespace std;

//////////////////////////////////////////////////////////
std :: string ElementStatsExporter :: giveFileName(unsigned step, int iteration) const {
    std :: ostringstream oss;
    if ( iteration < 0 ) {
        oss << filename << "_" << std :: setfill('0') << std :: setw(5) << step << ".dat";
    } else {
        oss << filename << "_" << std :: setfill('0') << std :: setw(5) << step
            << "_iter_" << std :: setfill('0') << std :: setw(5) << iteration << ".dat";
    }
    return oss.str();
}

//////////////////////////////////////////////////////////
void ElementStatsExporter :: readFromLine(istringstream &iss) {
    string param;
    this->remove_previous = false;
    last_saved_file = "none";

    for ( unsigned i = 0; i < elems->giveSize(); i++ ) {
        elems_to_save.push_back(i);
    }


    iss >> filename;
    while ( iss >> param ) {
        if ( param.compare("remove_previous") == 0 ) {
            iss >> this->remove_previous;
        }
    }
    DataExporter :: readFromLine(iss);
}


//////////////////////////////////////////////////////////
void ElementStatsExporter :: exportData(unsigned step, int iteration, fs :: path resultDir) const {
    std :: string fname = giveFileName(step, iteration);
    std :: string this_file_path = ( resultDir / fname ).string();
    double idc_time_converged = 0;
    double time_step;
    time_step = masterModel->giveSolver()->giveTimeStep();
    if ( masterModel->giveSolver()->giveName().compare("SteadyStateNonLinearSolver") == 0 ) {
        idc_time_converged = static_cast< SteadyStateNonLinearSolver * >( masterModel->giveSolver() )->giveIDCtime();
    }
    elems->saveElemStatsToFile(this_file_path, elems_to_save, 0, step, false, idc_time_converged, time_step);

    if ( remove_previous ) {
        if ( !masterModel->giveSolver()->isTerminated() && masterModel->giveSolver()->convergedToTolerance() ) {  // remove previous only if not terminated and reached tolerance in this step
            if ( this->last_saved_file.compare("none") != 0 ) {
                if ( std :: remove( last_saved_file.c_str() ) == 0 ) {
                    std :: cout << "element statuses saved to file " << this_file_path << '\n';
                } else {
                    std :: cerr << "previous file with elem statuses not removed" << '\n';
                }
            }
        }
        // here I need to modify private var in scope of const method
        // https://stackoverflow.com/questions/32686220/is-it-possible-change-value-of-member-variable-within-const-function
        ElementStatsExporter *ese = const_cast< ElementStatsExporter * >( this );
        ese->last_saved_file = this_file_path;
    }
}
