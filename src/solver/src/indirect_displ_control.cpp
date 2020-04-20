#include "indirect_displ_control.h"

//////////////////////////////////////////////////////////
void IndirectDC :: readFromLine(istringstream &iss) {
    string param;
    iss >> param;
    /*
    if ( param.compare("SteadyStateLinearSolver") == 0 ) {
        SteadyStateLinearSolver *newsolver = new SteadyStateLinearSolver();
        newsolver->readFromLine(iss);
        return newsolver;
    } else if ( param.compare("SteadyStateNonLinearSolver") == 0 ) {
        SteadyStateNonLinearSolver *newsolver = new SteadyStateNonLinearSolver();
        newsolver->readFromLine(iss);
        return newsolver;
    } else if ( param.compare("TransientLinearMechanicalSolver") == 0 ) {
        TransientLinearMechanicalSolver *newsolver = new TransientLinearMechanicalSolver();
        newsolver->readFromLine(iss);
        return newsolver;
    } else if ( param.compare("TransientLinearTransportSolver") == 0 ) {
        TransientLinearTransportSolver *newsolver = new TransientLinearTransportSolver();
        newsolver->readFromLine(iss);
        return newsolver;
    } else {
        cerr << "Error: Solver " << param << " is not implemented" << endl;
        exit(EXIT_FAILURE);
    };
    */
}

//////////////////////////////////////////////////////////
void IndirectDC :: init() {

}

