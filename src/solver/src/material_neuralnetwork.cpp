#include "material.h"
#include "element.h"
#include "material_neuralnetwork.h"
#include <torch/script.h>


using namespace std;


//////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////
// RECURRENT NETWORK TENSORIAL MECHANICAL MATERIAL

Vector NeuralNetworkMaterialStatus :: giveStress(const Vector &strain, double timeStep) {
    temp_strain = addEigenStrain(strain);

    NeuralNetworkMaterial *m = static_cast< NeuralNetworkMaterial * >( mat );

    // Normalization of input DoFs.
    Matrix norm = m->giveNormMatrix();
    Vector x_mean = norm.row(0).transpose();
    Vector x_std = norm.row(1).transpose();
    Vector y_mean = norm.row(2).head(6).transpose();
    Vector y_std = norm.row(3).head(6).transpose();

    // std::cout << "\nCheckpoint 01\n" << std::flush;

    // Vector strain_norm = ( temp_strain.array() - x_mean.array() ) / x_std.array();

    // std::cout << "Strain norm" << strain_norm << "\n" << std::flush;


    Vector input_norm = Eigen::VectorXd::Zero(10);
    
    input_norm.head(6) = temp_strain;
    input_norm(6) = m -> giveElasticModulus() / 1e6;
    input_norm(7) = m -> giveft() / 1e6;
    input_norm(8) = m -> giveGt();
    input_norm(9) = m -> giveRVEsize();

    Vector input_norm_N = ( input_norm.array() - x_mean.array() ) / x_std.array();

    // std::cout << "Input norm" << input_norm_N << "\n" << std::flush;


    Eigen :: VectorXf input_float = input_norm_N.cast< float > ();
    torch :: Tensor inputs_torch = torch :: from_blob(input_float.data(), { 1, 10 }).clone(); // Populates torch Tensor with Eigen Vector/Matrix
    // inputs_torch = inputs_torch.unsqueeze(0); // so that sequence len = 1

    inputs_torch.unsqueeze_(0);

    std :: vector< torch :: jit :: IValue >inputs;

    auto hc0_in = std::make_tuple(hc0[0], hc0[1]);
    auto hc1_in = std::make_tuple(hc1[0], hc1[1]);

    inputs.push_back(inputs_torch);
    // inputs.push_back(hc0);
    // inputs.push_back(hc1);
    inputs.push_back(hc0_in);
    inputs.push_back(hc1_in);

    // std::cout << "inputs build \n"  << "\n";

    auto outputs = m->giveNetwork().forward(inputs);
    auto tuple_output = outputs.toTuple();
    // std::cout << "prediction succesful\n"  << "\n";

    torch :: Tensor stress_norm_tensor = tuple_output->elements() [ 0 ].toTensor();
    // torch :: Tensor hidden_new = tuple_output->elements() [ 1 ].toTensor();

    auto hc0_tuple = tuple_output->elements()[1].toTuple();
    auto hc1_tuple = tuple_output->elements()[2].toTuple();

    // std::cout << "pre hidden update\n"  << "\n";


    // Update the hidden state
    torch::Tensor hidden0 = hc0_tuple->elements()[0].toTensor();
    // std::cout << hidden0  << "\n";
    // std::cout << hidden0.sizes()  << "\n";


    temp_hc0[0] = hc0_tuple->elements()[0].toTensor();
    temp_hc0[1] = hc0_tuple->elements()[1].toTensor();
    temp_hc1[0] = hc1_tuple->elements()[0].toTensor();
    temp_hc1[1] = hc1_tuple->elements()[1].toTensor();
   
    // std::cout << "post hidden update\n"  << "\n";


    // temp_hidden = hidden_new;

    stress_norm_tensor = stress_norm_tensor.squeeze(0); // gets rid of the "sequence" dimension
    stress_norm_tensor.squeeze_(0); 

    // std::cout << "\nCheckpoint 03\n" << std::flush;

    // Convert forces_norm Tensor to vector and to Eigen Vector and to double type
    std :: vector< float >outVec(stress_norm_tensor.data_ptr< float >(), stress_norm_tensor.data_ptr< float >() + stress_norm_tensor.numel() );
    Eigen :: VectorXf outFloat = Eigen :: VectorXf :: Map( & outVec [ 0 ], strain.size() );
    Eigen :: VectorXd stress_norm = outFloat.cast< double > ();
    // std::cout << "Stress Norm\n" << stress_norm << "\n" << std::flush;

    // // Denormalize output stress
    Vector stress = ( stress_norm.array() * y_std.array() ) + y_mean.array();

    // std::cout << "returning stress\n"  << "\n";

    return stress;
};

//////////////////////////////////////////////////////////
Matrix NeuralNetworkMaterialStatus :: giveDampingTensor() const {
    //TensMechMaterial *m = static_cast< TensMechMaterial * >( mat );
    Matrix M = Matrix :: Zero(1, 1);
    //M(0, 0) = m->giveDampingConstant();
    M(0, 0) = 0;
    return M;
}

//////////////////////////////////////////////////////////
double NeuralNetworkMaterialStatus :: giveMassConstant() const {
    NeuralNetworkMaterial *material = static_cast< NeuralNetworkMaterial * >( mat );
    return material->giveDensity();
}

//////////////////////////////////////////////////////////
Vector NeuralNetworkMaterialStatus :: giveStressWithFrozenIntVars(const Vector &strain, double timeStep) {
    ( void ) timeStep;
    temp_strain = addEigenStrain(strain);
    // temp_stress = giveStiffnessTensor("elastic") * temp_strain;
    temp_stress = giveStress(temp_strain, timeStep);
    return temp_stress;
};

//////////////////////////////////////////////////////////
Matrix NeuralNetworkMaterialStatus :: giveStiffnessTensor(string type) const {
    ( void ) type;

    NeuralNetworkMaterial *m = static_cast< NeuralNetworkMaterial * >( mat );

    if ( m->isStiffnessFromMatrix() ) {
        Matrix D = m->giveStiffnessMatrix();
        return D;
    } else {
        unsigned size = 1;
        unsigned dimension = mat->giveDimension();
        if ( dimension == 1 ) {
            size = 1;
        } else if ( dimension == 2 ) {
            size = 3;
        } else if ( dimension == 3 ) {
            return giveElasticStiffnessTensor3D();
        } else {
            cerr << name << ": unsupported dimension " << dimension << endl;
            exit(1);
        }
        Matrix D = Matrix :: Zero(size, size);

        if ( dimension == 1 ) {
            D(0, 0) = m->giveElasticModulus();
        } else if ( dimension == 2 ) {
            if ( m->isPlaneStress() ) {
                double factor = m->giveElasticModulus() / ( 1. - pow(m->givePoissonsRatio(), 2) );
                D(0, 0) = D(1, 1) = factor;
                D(0, 1) = D(1, 0) = m->givePoissonsRatio() * factor;
                D(2, 2) = ( 1. - m->givePoissonsRatio() ) / 2. * factor;
            } else {  //plane strain
                double factor = m->giveElasticModulus() / ( ( 1. - 2. * m->givePoissonsRatio() ) * ( 1. + m->givePoissonsRatio() ) );
                D(0, 0) = D(1, 1) = factor * ( 1. - m->givePoissonsRatio() );
                D(0, 1) = D(1, 0) = m->givePoissonsRatio() * factor;
                D(2, 2) = ( 1. - 2. * m->givePoissonsRatio() ) / 2. * factor;
            }
        } else {
            cerr << name << " error: dimension " << dimension << " not implemented" << endl;
            exit(1);
        }
        return D;
    }
};

//////////////////////////////////////////////////////////

Matrix NeuralNetworkMaterialStatus :: giveElasticStiffnessTensor3D() const {
    Matrix D = Matrix :: Zero(6, 6);
    NeuralNetworkMaterial *m = static_cast< NeuralNetworkMaterial * >( mat );
    double factor = m->giveElasticModulus() / ( 1. - 2. * m->givePoissonsRatio() ) / ( 1. + m->givePoissonsRatio() );
    D(0, 0) = D(1, 1) = D(2, 2) = factor * ( 1. - m->givePoissonsRatio() );
    D(0, 1) = D(1, 0) = D(0, 2) = D(2, 0) = D(2, 1) = D(1, 2) = m->givePoissonsRatio() * factor;
    D(3, 3) = D(4, 4) = D(5, 5) = ( 1. - 2. * m->givePoissonsRatio() ) / 2. * factor;
    return D;
};

//////////////////////////////////////////////////////////
void NeuralNetworkMaterialStatus :: update() {
    MaterialStatus :: update();
    // hidden = temp_hidden;
    hc0 = temp_hc0;
    hc1 = temp_hc1;

}

//////////////////////////////////////////////////////////
Matrix NeuralNetworkMaterialStatus :: giveMassTensor() const {
    NeuralNetworkMaterial *m = static_cast< NeuralNetworkMaterial * >( mat );
    unsigned dimension = mat->giveDimension();
    return Eigen :: MatrixXd :: Identity(dimension, dimension) * m->giveDensity();
}

//////////////////////////////////////////////////////////
// torch :: Tensor NeuralNetworkMaterialStatus :: giveHiddenState() const {
    // return hidden;
// }

std :: vector<std :: vector <torch::Tensor> > NeuralNetworkMaterialStatus :: giveHiddenState() const {
    std :: vector<std :: vector <torch::Tensor> > hs;
    hs.push_back(hc0);
    hs.push_back(hc1);
    return hs;

}


//////////////////////////////////////////////////////////
NeuralNetworkMaterialStatus :: NeuralNetworkMaterialStatus(NeuralNetworkMaterial *m, Element *e, unsigned ipnum) : MaterialStatus(m, e, ipnum) {
    name = "neural network tensorial mechanical mat. status";
    // std :: tuple< int, int >props = m->giveNetworkProps();
    // int num_layers = std :: get< 0 >(props);
    // int hidden_size = std :: get< 1 >(props);
    // hidden = torch :: zeros({ num_layers, hidden_size }); // initial hidden state as zero
    // temp_hidden = torch :: zeros({ num_layers, hidden_size }); // initial hidden state as zero

    std :: tuple< int, int >props = m->giveNetworkProps();
    int num_layers = std :: get< 0 >(props);
    int hidden_size = std :: get< 1 >(props);

    torch::Tensor h0 = torch :: zeros({ num_layers, hidden_size }); // initial hidden state as zero
    torch::Tensor c0 = torch :: zeros({ num_layers, hidden_size }); // initial hidden state as zero
    h0.unsqueeze_(0);
    c0.unsqueeze_(0);

    hc0.push_back(h0);
    hc0.push_back(c0);
    temp_hc0.push_back(h0);
    temp_hc0.push_back(c0);

    std :: tuple< int, int >props1 = m->giveNetworkProps1();
    int num_layers1 = std :: get< 0 >(props1);
    int hidden_size1 = std :: get< 1 >(props1);

    torch::Tensor h1 = torch :: zeros({ num_layers1, hidden_size1 }); // initial hidden state as zero
    torch::Tensor c1 = torch :: zeros({ num_layers1, hidden_size1 }); // initial hidden state as zero

    h1.unsqueeze_(0);
    c1.unsqueeze_(0);

    hc1.push_back(h1);
    hc1.push_back(c1);
    temp_hc1.push_back(h1);
    temp_hc1.push_back(c1);


}

//////////////////////////////////////////////////////////
bool NeuralNetworkMaterialStatus :: giveValues(string code, Vector &result) const {
    if ( code.compare("stress") == 0 || code.compare("stresses") == 0 ) {
        unsigned size = temp_stress.size();
        result.resize(size);
        for ( unsigned p = 0; p < size; p++ ) {
            result [ p ] = temp_stress [ p ];
        }
        return true;
    } else if ( code.compare("strain") == 0 || code.compare("strains") == 0 ) {
        unsigned size = temp_strain.size();
        result.resize(size);
        for ( unsigned p = 0; p < size; p++ ) {
            result [ p ] = temp_strain [ p ];
        }
        return true;
    } else {
        return MaterialStatus :: giveValues(code, result);
    }
}


Matrix NeuralNetworkMaterial :: readDataNormalizationMatrix(int size, fs :: path matrix_path) const {
    vector< double >matrixEntries;
    // ifstream matrixDataFile("C:/Users/209050/OAS_data/Plasticity/ML/ML_Torch_first_RT/data_normalization.txt");
    ifstream matrixDataFile(matrix_path);
    string matrixRowString;
    string matrixEntry;
    int matrixRowNumber = 0;
    int matrixColumnNumber = size;
    while ( getline(matrixDataFile, matrixRowString) ) {
        stringstream matrixRowStringStream(matrixRowString);
        while ( matrixRowStringStream >> matrixEntry ) {
            matrixEntries.push_back(stod(matrixEntry) );    //here we convert the string to double and fill in the row vector storing all the matrix entries
        }
        matrixRowNumber++; //update the row numbers
    }

    if ( matrixEntries.size() != matrixRowNumber * matrixColumnNumber ) {
        cerr << "Loaded matrix size problem " << endl;
        exit(1);
    }
    Matrix A(matrixRowNumber, matrixColumnNumber);
    for ( int i = 0; i < matrixRowNumber; i++ ) {
        for ( int j = 0; j < matrixColumnNumber; j++ ) {
            A(i, j) = matrixEntries [ i * matrixColumnNumber + j ];
        }
    }

    return A;
}


Matrix NeuralNetworkMaterial :: readStiffMatrixFromFile() const {
    vector< double >matrixEntries;
    ifstream matrixDataFile(sm_path);
    string matrixRowString;
    string matrixEntry;
    int matrixRowNumber = 0;
    while ( getline(matrixDataFile, matrixRowString) ) {
        stringstream matrixRowStringStream(matrixRowString);
        while ( matrixRowStringStream >> matrixEntry ) {
            matrixEntries.push_back( stod(matrixEntry) );    //here we convert the string to double and fill in the row vector storing all the matrix entries
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
void NeuralNetworkMaterial :: readFromLine(istringstream &iss) {
    string param;
    bool bstiffmat, bmlmodel, blayers, bhidden, blayers1, bhidden1, bnormmat, bE, bnu, bdensity, bft, bGt, bRVEsize;
    bstiffmat = bmlmodel = bnormmat = bE = bnu = blayers = bhidden = blayers1 = bhidden1 = bdensity = bft = bGt = bRVEsize = false;

    while (  iss >> param ) {
        if ( param.compare("E") == 0 ) {
            bE = true;
            iss >> E;
        } else if ( param.compare("nu") == 0 ) {
            bnu = true;
            iss >> nu;
        } else if ( param.compare("density") == 0 ) {
            bdensity = true;
            iss >> density;
        } else if ( param.compare("ft") == 0 ) {
            bft = true;
            iss >> ft;
        } else if ( param.compare("Gt") == 0 ) {
            bGt = true;
            iss >> Gt;
        } else if ( param.compare("RVEsize") == 0 ) {
            bRVEsize = true;
            iss >> RVEsize;
        } else if ( param.compare("planeStrain") == 0 ) {
            planeStress = false;
        } else if ( param.compare("stiff_mat") == 0 ) {
            bstiffmat = true;
            string filepath;
            iss >> filepath;
            sm_path = GlobPaths :: BASEDIR  / filepath;
            stiffmat_elastic = readStiffMatrixFromFile();
        } else if ( param.compare("ml_model") == 0 ) {
            bmlmodel = true;
            string filepath;
            iss >> filepath;
            ml_path = GlobPaths :: BASEDIR  / filepath;
            // std::cout << "\nnorm Matrix\n" << norm << "\n" << std::flush;
            network = torch :: jit :: load(ml_path.string() );
        } else if ( param.compare("norm_mat") == 0 ) {
            bnormmat = true;
            string filepath;
            iss >> filepath;
            nm_path = GlobPaths :: BASEDIR  / filepath;
            norm = readDataNormalizationMatrix(10, nm_path);
        } else if ( param.compare("num_layers") == 0 ) {
            blayers = true;
            iss >> num_layers;
        } else if ( param.compare("hidden_size") == 0 ) {
            bhidden = true;
            iss >> hidden_size;
        } else if ( param.compare("num_layers1") == 0 ) {
            blayers1 = true;
            iss >> num_layers1;
        } else if ( param.compare("hidden_size1") == 0 ) {
            bhidden1 = true;
            iss >> hidden_size1;
        }  else {
            cerr << "MLElement ERROR: " << param << " input parameter not defined \n";
        }
    }

    if ( !bstiffmat && !bE ) {
        cerr << name << ": no elastic stiffness matrix neither E, nu was specified" << endl;
        exit(EXIT_FAILURE);
    }
    if ( !bstiffmat && ( bE || bnu ) ) {
        cout << name << ": elasticity defined by E,nu" << endl;
        stiffness_from_matrix = false;
        if ( !bE ) {
            cerr << name << ": E was not specified" << endl;
            exit(EXIT_FAILURE);
        }
        if ( !bnu ) {
            cerr << name << ": nu was not specified" << endl;
            exit(EXIT_FAILURE);
        }
    }
    if ( !bE && !bnu ) {
        cout << name << ": elasticity defined by elastic stiffness matrix" << endl;
        stiffness_from_matrix = true;

        if ( !bstiffmat ) {
            cerr << name << ":  elastic stiffness matrix was not specified" << endl;
            exit(EXIT_FAILURE);
        }
    }
    if ( !bft ) {
        cerr << name << ": ft model was not specified" << endl;
        exit(EXIT_FAILURE);
    }
    if ( !bGt ) {
        cerr << name << ": Gt model was not specified" << endl;
        exit(EXIT_FAILURE);
    }
    if ( !bRVEsize ) {
        cerr << name << ": RVEsize model was not specified" << endl;
        exit(EXIT_FAILURE);
    }

    if ( !bmlmodel ) {
        cerr << name << ": ML model was not specified" << endl;
        exit(EXIT_FAILURE);
    }
    if ( !blayers ) {
        cerr << name << ": NN number of layers was not specified" << endl;
        exit(EXIT_FAILURE);
    }
    if ( !bhidden ) {
        cerr << name << ": NN hiden size was not specified" << endl;
        exit(EXIT_FAILURE);
    }
    if ( !bnormmat ) {
        cerr << name << ": normalization matrix was not specified" << endl;
        exit(EXIT_FAILURE);
    }
};

//////////////////////////////////////////////////////////
MaterialStatus *NeuralNetworkMaterial :: giveNewMaterialStatus(Element *e, unsigned ipnum) {
    NeuralNetworkMaterialStatus *newStatus = new NeuralNetworkMaterialStatus(this, e, ipnum); //needs to be deleted manually
    return newStatus;
};
