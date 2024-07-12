#include <torch/script.h>
#include "element_superelem.h"
#include "boundary_condition.h"

using namespace std;

//////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////
// ML ELEMENT
MLMechElement :: MLMechElement(unsigned dim) : Element(dim) {
    numOfNodes = 0;
    name = "MLElement";
    vtk_cell_type = 0;
    shafunc = new NullShapeF(ndim);
    inttype = new EmptyIntegration();
    physicalFields [ 0 ] = true; //mechanics
}

//////////////////////////////////////////////////////////
Matrix MLMechElement :: readDataNormalizationMatrix(int size) const {
    vector<double> matrixEntries;
    // ifstream matrixDataFile("C:/Users/209050/OAS_data/Plasticity/ML/ML_Torch_first_RT/data_normalization.txt");
    ifstream matrixDataFile(nm_path);
    string matrixRowString;
    string matrixEntry;
    int matrixRowNumber = 0;
    int matrixColumnNumber = size;
    while (getline(matrixDataFile, matrixRowString))
    {
        stringstream matrixRowStringStream(matrixRowString); 
        while (matrixRowStringStream >> matrixEntry)
        {
            matrixEntries.push_back(stod(matrixEntry));   //here we convert the string to double and fill in the row vector storing all the matrix entries
        }
        matrixRowNumber++; //update the row numbers
    } 

    // std::cout << "\nmatrixColumnNumber\n" << matrixColumnNumber << "\n" << std::flush;
    // std::cout << "\nmatrixRowNumber\n" << matrixRowNumber << "\n" << std::flush;
    // std::cout << "\nmatrixEntries.size()\n" << matrixEntries.size() << "\n" << std::flush;

    if (matrixEntries.size()!=matrixRowNumber*matrixColumnNumber) {
        cerr << "Loaded matrix size problem " << endl;
        exit(1);
    }
    Matrix A(matrixRowNumber,matrixColumnNumber);
    for(int i=0; i<matrixRowNumber; i++){
        for(int j=0; j<matrixColumnNumber; j++){
            A(i,j) = matrixEntries[i*matrixColumnNumber+j];
        }   
    }

    return A;
}

//////////////////////////////////////////////////////////
Matrix MLMechElement :: readStiffMatrixFromFile() const {
    vector< double >matrixEntries;
    ifstream matrixDataFile(sm_path);
    string matrixRowString;
    string matrixEntry;
    int matrixRowNumber = 0;
    while ( getline(matrixDataFile, matrixRowString) ) {
        stringstream matrixRowStringStream(matrixRowString);
        while ( matrixRowStringStream >> matrixEntry ) {
            matrixEntries.push_back( stod(matrixEntry) );   //here we convert the string to double and fill in the row vector storing all the matrix entries
        }
        matrixRowNumber++; //update the column numbers
    }

    if ( matrixEntries.size() != matrixRowNumber * matrixRowNumber ) {
        cerr << "Loaded matrix is not a square matrix, found " << matrixEntries.size() << ", expected " << matrixRowNumber * matrixRowNumber << endl;
        exit(1);
    }
    Matrix A(matrixRowNumber, matrixRowNumber);
    for ( int i = 0; i < matrixRowNumber; i++ ) {
        for ( int j = 0; j < matrixRowNumber; j++ ) {
            A(i, j) = matrixEntries [ i * matrixRowNumber + j ];
        }
    }

    return A;
}

//////////////////////////////////////////////////////////
void MLMechElement :: readFromLine(std :: istringstream &iss, NodeContainer *fullnodes, MaterialContainer *fullmatrs) {
    ( void ) fullmatrs;

    unsigned num;
    iss >> numOfNodes;
    nodes.resize(numOfNodes);
    for ( unsigned k = 0; k < numOfNodes; k++ ) {
        iss >> num;
        nodes [ k ] = fullnodes->giveNode(num);
    }

    string param;
    while (  iss >> param ) {
        if ( param.compare("poly_degree") == 0 ) {
            iss >> poly_degree;
        } else if ( param.compare("stiff_mat") == 0 ) {
            string filepath;
            iss >> filepath;
            sm_path = GlobPaths :: BASEDIR  / filepath;
        } else if ( param.compare("ml_model") == 0 ) {        
            string filepath;    
            iss >> filepath;
            ml_path = GlobPaths :: BASEDIR  / filepath;
            // std::cout << "\nnorm Matrix\n" << norm << "\n" << std::flush;
        } else if ( param.compare("norm_mat") == 0 ) {        
            string filepath;    
            iss >> filepath;
            nm_path = GlobPaths :: BASEDIR  / filepath;
        }
    }
}

//////////////////////////////////////////////////////////
void MLMechElement :: init() {
    Element :: init();

    //reduce stiffness matrix
    vector< unsigned >keep_ind;
    //nodal displacement
    for ( unsigned i = 0; i < 8; i++ ) {
        keep_ind.push_back(i);
    }
    // i runs over sides
    for ( unsigned i = 0; i < 8; i++ ) {
        // j runs over polynomial constants, linear basis means no additional dofs
        for ( unsigned j = 0; j < poly_degree - 1; j++ ) {
            keep_ind.push_back(8 + i * 4 + j);
        }
    }


    Matrix A = readStiffMatrixFromFile();
    stiffmat = A(keep_ind, keep_ind);

    if ( outDoFs != keep_ind.size() ) {
        cerr << "Error in MLMechElement: there are " << outDoFs << " DoFs on input of the element, but " << keep_ind.size() << " DoFs are required for polynom of degree " << poly_degree << endl;
        exit(1);
    }

    // load pytorch libtorch model
    module = torch::jit::load(ml_path.string()); 

    // load normalization matrix
    int size = keep_ind.size();
    norm = readDataNormalizationMatrix(size);

}

//////////////////////////////////////////////////////////
Matrix MLMechElement :: giveStiffnessMatrix(std :: string matrixType) const {
    ( void ) matrixType;
    return stiffmat;
}

//////////////////////////////////////////////////////////
Vector MLMechElement :: giveInternalForces(const Vector &DoFs, bool frozen, double timeStep){
    if (frozen){
    (void) timeStep;
    return stiffmat*DoFs;

    } else {
    (void) timeStep;
    
    int size = DoFs.size();
    int size_x = size - 2;
    // std::cout << "\nCheckpoint 00\n" << std::flush;
    // std::cout << "\nnorm Matrix\n" << norm << "\n" << std::flush;

    // Relative DoFs - with corner0 displacements = 0
    Vector DoFs_rel = DoFs.cast <double> ();
    double c0u = DoFs_rel[0];
    double c0v = DoFs_rel[1];

    for (int i = 0; i < 8; i+=2){
        DoFs_rel[i] -= c0u;
    }
    for (int i = 1; i < 8; i+=2){
        DoFs_rel[i] -= c0v;
    }

    // Normalization of input DoFs.
    Vector x_std = norm.row(1);
    Vector DoFs_norm = DoFs_rel.cast <double> ();
    DoFs_norm -= norm.row(0);
    // DoFs_norm = DoFs_norm.array() / norm.row(1).array(); // doesnt work for some reason
    DoFs_norm = DoFs_norm.array() / x_std.array();
    DoFs_norm.head(size_x) = DoFs_norm.segment(size - size_x, size_x);

    // bisymmetric log transformation
    // double Cx = 1 / log(10);
    // Eigen::VectorXd vec1 (size);
    // vec1.fill(1);
    // DoFs_norm = DoFs_norm.array().sign() * (vec1 + (DoFs_norm / Cx).cwiseAbs()).array().log10();
    // DoFs_norm.matrix();


    // Create a vector of inputs.
    std::vector<torch::jit::IValue> inputs;
    // inputs.push_back(torch::ones({1,24}));
    // Change DoFs type to float - Might be possible to avoid if a double type libtorch model is possible
    // Eigen::VectorXf DoFsFloat = DoFs.cast <float> (); 
    Eigen::VectorXf DoFsFloat = DoFs_norm.cast <float> (); 
    // std::cout << "\nDofs\n" << DoFs << "\n";
    // std::cout << "\nDofsFloat\n" << DoFsFloat << "\n";
  
    // // Create a torch Tensor populated by the DoFs values
    torch::Tensor inputs_torch = torch::from_blob(DoFsFloat.data(), {1,size_x}).clone(); // Populates torch Tensor with Eigen Vector/Matrix
    inputs.push_back(inputs_torch);

    // std::cout << "input tensor_torch\n" << inputs_torch << "\n";
    // std::cout << "input tensor\n" << inputs << "\n";

    // std::cout << "\nChceckpoint 1\n" << std::flush;
    // // Execute the model and turn its output into a tensor.
    auto output = module.forward(inputs).toTensor(); 

    // std::cout << "\nChceckpoint 2\n" << std::flush;
    // // Convert output Tensor to vector and to Eigen Vector and to double type
    std::vector<float> outVec(output.data_ptr<float>(), output.data_ptr<float>() + output.numel());
    Eigen::VectorXf outFloat = Eigen::VectorXf::Map(&outVec[0], size);
    Eigen::VectorXd forces_norm = outFloat.cast <double> ();
    // std::cout << "\nForces Norm\n" << forces_norm << "\n" << std::flush;

    // Denormalize output forces
    Vector y_std = norm.row(3);
    Vector forces = forces_norm.array() * y_std.array();
    forces += norm.row(2);

    // bisymmetric log inverse transfomration
    // double Cy = 1 / log(10);
    // Eigen::VectorXd vecm1 (size);
    // vecm1.fill(-1);
    // Eigen::VectorXd pow10 (size);
    // for (int i = 0; i < size; i++){
    //     pow10[i] = pow(10, fabs(forces[i]));
    // }
    // forces = (forces.array().sign().matrix() * Cy).array() * (vecm1 + pow10).array();


    // std::cout << "\nForces ML\n" << forces << "\n" << std::flush;
    // std::cout << "\nForces Ku\n" << stiffmat*DoFs << "\n";

    // return forces;

    // std::cout << "\nDofs\n" << DoFs << "\n";
    // std::cout << "\nStiffmat\n" << stiffmat << "\n";
    // std::cout << "\nForces Correct\n" << stiffmat*DoFs << "\n";


    // std::ofstream outFile("C:/Users/209050/OAS_data/Plasticity/Constrained/RegMesh/Rectangle_D/OutDofForce.txt");
    // outFile << "" ;
    // outFile.close();


    // std::ofstream outfile("C:/Users/209050/OAS_data/Plasticity/Constrained/RegMesh/Rectangle_D/OutDofForce.txt", std::ios::app);
    // if (!outfile.is_open()) {
    //     std::cerr << "Error opening file!" << std::endl;
    //     // return 1;
    // }

    // std::ofstream outFL("C:/Users/209050/OAS_data/Plasticity/Constrained/RegMesh/Rectangle_D/OutDoF.txt", std::ios::app);

    // // Write the array to the file

    // string outstirng;
    // outfile << "DoFs" << "\n";
    // for (double value : DoFs) {
    //     outstirng.append(std::to_string(value) + ", ");
    // }
    // outstirng.pop_back();
    // outstirng.pop_back();
    // outfile << outstirng << "\n";
    // outFL << outstirng << "\n";

    // outstirng = "";

    // outfile << "DoFs_REL" << "\n";
    // for (double value : DoFs_rel) {
    //     outstirng.append(std::to_string(value) + ", ");
    // }
    // outstirng.pop_back();
    // outstirng.pop_back();
    // outfile << outstirng << "\n";
    // outFL << outstirng << "\n";

    // outstirng = "";
    
    // outfile << "forces" << "\n";
    // for (double value : forces) {
    //     outstirng.append(std::to_string(value) + ", ");
    // }
    // outstirng.pop_back();
    // outstirng.pop_back();
    // outfile << outstirng << "\n";
    // outstirng = "";
    
    // Vector sdf = stiffmat*DoFs - forces;
    // outfile << "stiffmat*DoFs - forces" << "\n";
    // for (double value : sdf) {
    //     outstirng.append(std::to_string(value) + ", ");
    // }
    // outstirng.pop_back();
    // outstirng.pop_back();
    // outfile << outstirng << "\n";
    // outstirng = "";

    // Vector sd = stiffmat*DoFs;
    // outfile << "stiffmat*DoFs" << "\n";
    // for (double value : sd) {
    //     outstirng.append(std::to_string(value) + ", ");
    // }
    // outstirng.pop_back();
    // outstirng.pop_back();
    // outfile << outstirng << "\n";
    // outstirng = "";
    // outfile << "\n";


    // // outfile << forces << "\n";
    // // outfile << "stiffmat*DoFs - forces" << "\n";
    // // outfile << stiffmat*DoFs - forces << "\n";
    // // outfile << "stiffmat*DoFs" << "\n";
    // // outfile << stiffmat*DoFs << "\n";

    // // Close the file
    // outfile.close();
    // outFL.close();



    return stiffmat*DoFs - forces;
    // return stiffmat*DoFs;


    }
    

    // return stiffmat*DoFs;



    //  if (frozen){
        //    std::cout << "\nForces Correct\n" << stiffmat*DoFs << "\n";
    //     return stiffmat*DoFs;

    // } else {
    //        std::cout << "\nForces Correct\n" << stiffmat*DoFs << "\n";

    //     return stiffmat*DoFs;
    // }

}



//////////////////////////////////////////////////////////
Vector MLMechElement :: integrateInternalSources() {
    return Vector :: Zero(0);
}

//////////////////////////////////////////////////////////
void MLMechElement :: giveValues(std :: string code, Vector &result) const {
    ( void ) code;
    ( void ) result;
}


/*
 * //////////////////////////////////////////////////////////
 * Matrix MLMechElement :: giveBMatrix(const Point *x) const {
 *  cout << "MLMechElement :: giveBMatrix should not be called" << endl;
 *  exit(1);
 *  return Matrix::Zeros(0,0);
 * }
 *
 * //////////////////////////////////////////////////////////
 * Matrix MLMechElement :: giveHMatrix(const Point *x) const {
 *  cout << "MLMechElement :: giveHMatrix should not be called" << endl;
 *  exit(1);
 *  return Matrix::Zeros(0,0);
 * }
 */
//////////////////////////////////////////////////////////
Vector MLMechElement :: giveStrain(unsigned i, const Vector &DoFs) {
    ( void ) i;
    ( void ) DoFs;
    cout << "MLMechElement :: giveStrain should not be called" << endl;
    exit(1);
    return Vector :: Zero(0);
}


//////////////////////////////////////////////////////////
Matrix MLMechElement :: giveMassMatrix() {
    cout << "MLMechElement :: giveBMatrix should not be called" << endl;
    exit(1);
    return Matrix :: Zero(0, 0);
}

//////////////////////////////////////////////////////////
Matrix MLMechElement :: giveDampingMatrix() {
    cout << "MLMechElement :: giveBMatrix should not be called" << endl;
    exit(1);
    return Matrix :: Zero(0, 0);
}

//////////////////////////////////////////////////////////
Matrix MLMechElement :: giveLumpedMassMatrix() {
    cout << "MLMechElement :: giveBMatrix should not be called" << endl;
    exit(1);
    return Matrix :: Zero(0, 0);
}
