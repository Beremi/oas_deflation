#include "solver_eigenvalue.h"
#include "adaptivity.h"
#define numPhysicalFields 4

using namespace std;


//////////////////////////////////////////////////////////
EigenvalueMechanicalSolver :: EigenvalueMechanicalSolver() {
    name = "EigenvalueMechanicalSolver";
    //use_lumped_M = true;
    num_eigs = 0;
    time = 0;
    dt = 1;
    termination_time = 1;
}

//////////////////////////////////////////////////////////
EigenvalueMechanicalSolver :: ~EigenvalueMechanicalSolver() {}

//////////////////////////////////////////////////////////
void EigenvalueMechanicalSolver :: init(std :: string init_r_file, std :: string init_v_file, const bool initial) {
    Solver :: init(init_r_file, init_v_file, initial);

    elems->prepareMassMatrix(M, 1);
    elems->updateMassMatrix(M, 1);
    if ( nodes->giveConstraints()->isActive() ) {
        nodes->giveConstraints()->transformToConstraintSpace(M);
    }
    M.makeCompressed();
    elems->prepareStiffnessMatrix(K);
    K.makeCompressed();
    lumpedM = M.diagonal();
}

//////////////////////////////////////////////////////////
void EigenvalueMechanicalSolver :: runBeforeEachStep() {
    Solver :: runBeforeEachStep();
    cout << "######### Solving eigenvalue problem #########" << endl;
}

//////////////////////////////////////////////////////////
Solver *EigenvalueMechanicalSolver :: readFromFile(const string filename) {
    string param, line;
    bool bnum = false;
    ifstream inputfile( filename.c_str() );
    if ( inputfile.is_open() ) {
        while ( getline(inputfile >> std :: ws, line) ) {
            if ( line.empty() || ( line.at(0) == '#' ) ) {
                continue;
            }
            istringstream iss(line);
            iss >> param;
            if ( param.compare("num_of_eigenvalues") == 0 ) {
                bnum = true;
                iss >> num_eigs;
            }
        }
        inputfile.close();
    }
    if ( !bnum ) {
        cerr << name << ": number of eigenvalues 'eigen_num' not specified" << endl;
        exit(EXIT_FAILURE);
    }
    return this;
};

//////////////////////////////////////////////////////////
void EigenvalueMechanicalSolver :: solve() {
    Vector eigenvalues;
    Matrix eigenvectors;
    //LinalgEigenSpectraSolver(M, eigenvalues, eigenvectors, num_eigs);
    //LinalgEigenSpectraGENSolver(M, K, eigenvalues, eigenvectors, num_eigs);

    //cout << "Eigenvalues: " << eigenvalues << endl;

    //update velocity and acceleration
    //nodes->giveFullDoFArray(v_red, v);
    //nodes->giveFullDoFArray(a_red, a);
    //vector< double >blocked = bcs->giveBlockedDoFValues(time);

    //DOES NOT SUPPORT CONSTRAINT WITH CONJUGATE VARIABLES AND FUNCTIONS
    //nodes->giveConstraints()->calculateDependentDoFs(v);
    //nodes->giveConstraints()->calculateDependentDoFs(a);
}

//////////////////////////////////////////////////////////
void EigenvalueMechanicalSolver :: giveValues(std :: string code, Vector &result) const {
    if ( code.compare("???????") == 0 ) {
        result.resize(1);
        result [ 0 ] = 0;
    } else {
        Solver :: giveValues(code, result);
    }
}
