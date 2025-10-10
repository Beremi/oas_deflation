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


    Vector input_norm = Eigen :: VectorXd :: Zero(10);

    input_norm.head(6) = temp_strain;
    input_norm(6) = m->giveE0() / 1e6;
    input_norm(7) = m->giveft() / 1e6;
    input_norm(8) = m->giveGt();
    input_norm(9) = m->giveRVEsize();

    Vector input_norm_N = ( input_norm.array() - x_mean.array() ) / x_std.array();

    // std::cout << "Input norm" << input_norm_N << "\n" << std::flush;

    Eigen :: VectorXf input_float = input_norm_N.cast< float > ();
    torch :: Tensor inputs_torch = torch :: from_blob(input_float.data(), { 1, 10 }).clone(); // Populates torch Tensor with Eigen Vector/Matrix
    // inputs_torch = inputs_torch.unsqueeze(0); // so that sequence len = 1

    inputs_torch = inputs_torch.unsqueeze(0);

    // std::cout << "prev_strain" << prev_strain << "\n" << std::flush;
    // std::cout << "prev_strain size: " << prev_strain.size() << "\n";

    ///// strain increment - calculation, normalization, casting to torch tensor
    Vector strain_increment = temp_strain - prev_strain;
    Vector strain_increment_norm = ( strain_increment.array() - x_mean.head(6).array() ) / x_std.head(6).array();
    Eigen :: VectorXf increment_float = strain_increment_norm.cast< float > ();
    torch :: Tensor increment_torch = torch :: from_blob(increment_float.data(), { 1, 6 }).clone(); // Populates torch Tensor with Eigen Vector/Matrix
    increment_torch = increment_torch.unsqueeze(0);
    // std::cout << "increment_torch norm" << increment_torch << "\n" << std::flush;


    // inputs into form that libtorch wants
    std :: vector< torch :: jit :: IValue >inputs;
    inputs.push_back(inputs_torch);

    if (m -> give_use_strain_increment()) {
        inputs.push_back(increment_torch);
    }

    std::vector<Layer> layers = m->giveNetworkProps();
    for (size_t i = 0; i < layers.size(); ++i) {
        Layer& layer = layers[i]; // <- note the &
        if (layer.name == "LSTM") {
            auto hc_in = std::make_tuple(hc_vectors[i][0].clone(), hc_vectors[i][1].clone());
            inputs.push_back(hc_in);
        } else if (layer.name == "GRU") {
            auto h_in = hc_vectors[i][0].clone();
            inputs.push_back(h_in);
        } else if (layer.name == "LMSC") {
            // auto h_in = hc_vectors[i][0];
            auto h_in = hc_vectors[i][0].clone();

            inputs.push_back(h_in);
        }
    }

    // std::cout << "inputs \n"  << inputs << "\n";

    // network prediction
    auto outputs = m->giveNetwork().forward(inputs);
    // std::cout << "prediction succesful00\n"  << "\n";

    auto tuple_output = outputs.toTuple();
    // std::cout << "prediction succesful\n"  << "\n";

    // std::cout << "output \n"  << tuple_output << "\n";

    torch :: Tensor stress_norm_tensor = tuple_output->elements() [ 0 ].toTensor();
    // std::cout << "stress norm tensor \n"  << stress_norm_tensor << "\n";
    // torch :: Tensor hidden_new = tuple_output->elements() [ 1 ].toTensor();
    // std::cout << "hidden_new \n"  << hidden_new << "\n";

    // std::cout << "hc_vectors_prev \n"  << hc_vectors << "\n";
    
    // Update hidden states
    for (size_t i = 0; i < layers.size(); ++i) {
        Layer& layer = layers[i]; // <- note the &
        if (layer.name == "LSTM") {
            auto hc_tuple = tuple_output->elements()[i+1].toTuple();
            temp_hc_vectors[i][0] = hc_tuple->elements()[0].toTensor(); // update h tensor
            temp_hc_vectors[i][1] = hc_tuple->elements()[1].toTensor(); // update c tensor
        } else if ((layer.name == "GRU") || (layer.name == "LMSC")) {
            auto h_tensor = tuple_output->elements()[i+1].toTensor();
            temp_hc_vectors[i][0] = h_tensor;
        }
    }
;
    // std::cout << "hc_vectors_post \n"  << hc_vectors << "\n";

    // std::cout << "post hidden update\n"  << "\n";

    // temp_hidden = hidden_new;

    stress_norm_tensor = stress_norm_tensor.squeeze(0);
    stress_norm_tensor = stress_norm_tensor.squeeze(0); // !! do not use _squeeze() - it is IN-PLACE and messes up the tensor memory address
    
    // std::cout << "\nCheckpoint 03\n" << std::flush;

    // Convert forces_norm Tensor to vector and to Eigen Vector and to double type
    std :: vector< float >outVec( stress_norm_tensor.data_ptr< float >(), stress_norm_tensor.data_ptr< float >() + stress_norm_tensor.numel() );
    Eigen :: VectorXf outFloat = Eigen :: VectorXf :: Map(& outVec [ 0 ], strain.size() );
    Eigen :: VectorXd stress_norm = outFloat.cast< double > ();
    // std::cout << "Stress Norm\n" << stress_norm << "\n" << std::flush;

    // // Denormalize output stress
    temp_stress = ( stress_norm.array() * y_std.array() ) + y_mean.array();
    // std::cout << "Stress\n" << temp_stress << "\n" << std::flush;    

    // std::cout << "returning stress\n"  << "\n";
    return temp_stress;
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
    return giveStress(strain, timeStep);
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
    for (size_t i = 0; i < hc_vectors.size(); ++i) {
        hc_vectors [ i ] = temp_hc_vectors [ i ];
    }
    prev_strain = temp_strain;
}

//////////////////////////////////////////////////////////
Matrix NeuralNetworkMaterialStatus :: giveMassTensor() const {
    NeuralNetworkMaterial *m = static_cast< NeuralNetworkMaterial * >( mat );
    unsigned dimension = mat->giveDimension();
    return Eigen :: MatrixXd :: Identity(dimension, dimension) * m->giveDensity();
}

//////////////////////////////////////////////////////////
std :: vector< std :: vector< torch :: Tensor > >NeuralNetworkMaterialStatus :: giveHiddenState() const {
    return temp_hc_vectors;
}

//////////////////////////////////////////////////////////
NeuralNetworkMaterialStatus :: NeuralNetworkMaterialStatus(NeuralNetworkMaterial *m, Element *e, unsigned ipnum) : MaterialStatus(m, e, ipnum)  {
    name = "neural network tensorial mechanical mat. status";

    // prev_strain = Vector::Zero(temp_strain.size()); // doesn't work, needs to be done differently (with m -> giveDimension())
    prev_strain = Vector::Zero(6);

    std::vector<Layer> layers = m->giveNetworkProps();

    for (size_t i = 0; i < layers.size(); ++i) {
        Layer &layer = layers [ i ]; // <- note the &
        int num_layers = layer.num_layers;
        int hidden_size = layer.hidden_size;

        if ( layer.name == "LSTM" ) {
            std :: vector< torch :: Tensor >hc;
            torch :: Tensor h = torch :: zeros({ num_layers, hidden_size }); // initial hidden state as zero
            torch :: Tensor c = torch :: zeros({ num_layers, hidden_size }); // initial hidden state as zero
            h.unsqueeze_(0); // add the "batch_size" dimension (in this case batch_size = 1)
            c.unsqueeze_(0); // add the "batch_size" dimension (in this case batch_size = 1)
            hc.push_back(h);
            hc.push_back(c);
            hc_vectors.push_back(hc);
            temp_hc_vectors.push_back(hc);

        } else if ((layer.name == "GRU") || (layer.name == "LMSC") ) {
            std::vector<torch::Tensor> hc;
            torch::Tensor h = torch :: zeros({ num_layers, hidden_size }); // initial hidden state as zero
            h.unsqueeze_(0); // add the "batch_size" dimension (in this case batch_size = 1)
            hc.push_back(h);
            hc_vectors.push_back(hc);
            temp_hc_vectors.push_back(hc);
        }
        // std::cout << layer.name << "  " << layer.num_layers << "  " << layer.hidden_size << "\n" << std::flush;

        
    }
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
            matrixEntries.push_back( stod(matrixEntry) );    //here we convert the string to double and fill in the row vector storing all the matrix entries
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
            matrixEntries.push_back(stod(matrixEntry) );     //here we convert the string to double and fill in the row vector storing all the matrix entries
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
    bool bstiffmat, bmlmodel, bnormmat, bE, bnu, bdensity, bE0, bft, bGt, bRVEsize;
    bstiffmat = bmlmodel = bnormmat = bE = bnu = bdensity = bE0 = bft = bGt = bRVEsize = false;
    const int MAX_LAYERS = 5;

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
        } else if ( param.compare("E0") == 0 ) {
            bE0 = true;
            iss >> E0;
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
            try {
                network = torch :: jit :: load( ml_path.string() );
            } catch ( const c10 :: Error &e ) {
                std :: cerr << "Error loading the model: " << e.what() << std :: endl;
                exit(EXIT_FAILURE);
            }
        } else if ( param.compare("norm_mat") == 0 ) {
            bnormmat = true;
            string filepath;
            iss >> filepath;
            nm_path = GlobPaths :: BASEDIR  / filepath;
            norm = readDataNormalizationMatrix(10, nm_path);
        } else if ( param.find("layer_type") != std :: string :: npos ) {
            for (int i = 0; i < MAX_LAYERS; ++i) {
                std :: string layer_type = "layer_type" + std :: to_string(i);
                if ( param == layer_type ) {
                    layers.push_back({ "x", 0, 0 });
                    iss >> layers [ i ].name;
                }
            }
        } else if ( param.find("num_layers") != std :: string :: npos ) {
            for (int i = 0; i < MAX_LAYERS; ++i) {
                std :: string num_lrs = "num_layers" + std :: to_string(i);
                if ( param == num_lrs ) {
                    iss >> layers [ i ].num_layers;
                }
            }
        } else if ( param.find("hidden_size") != std :: string :: npos ) {
            for (int i = 0; i < MAX_LAYERS; ++i) {
                std :: string hdns = "hidden_size" + std :: to_string(i);
                if ( param == hdns ) {
                    iss >> layers [ i ].hidden_size;
                }
            }
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
    if ( !bE0 ) {
        cerr << name << ": E0 was not specified" << endl;
        exit(EXIT_FAILURE);
    }
    if ( !bft ) {
        cerr << name << ": ft was not specified" << endl;
        exit(EXIT_FAILURE);
    }
    if ( !bGt ) {
        cerr << name << ": Gt was not specified" << endl;
        exit(EXIT_FAILURE);
    }
    if ( !bRVEsize ) {
        cerr << name << ": RVEsize was not specified" << endl;
        exit(EXIT_FAILURE);
    }

    if ( !bmlmodel ) {
        cerr << name << ": ML model was not specified" << endl;
        exit(EXIT_FAILURE);
    }
    if ( !bnormmat ) {
        cerr << name << ": normalization matrix was not specified" << endl;
        exit(EXIT_FAILURE);
    }

    for (size_t i = 0; i < layers.size(); ++i) {
        Layer &layer = layers [ i ];
        if ( layer.name == "LSTM" || layer.name == "GRU" ) {
            if ( layer.hidden_size == 0 || layer.num_layers == 0 ) {
                cerr << "layer" << i << " has incorrectly specified num_layers of hidden_size" << endl;
                exit(EXIT_FAILURE);
            }
        }
        std::cout << layer.name << "  " << layer.num_layers << "  " << layer.hidden_size << "\n" << std::flush;

        if (layer.name == "LMSC") {
            use_strain_increment = true;
        }
    }

    
};

//////////////////////////////////////////////////////////
MaterialStatus *NeuralNetworkMaterial :: giveNewMaterialStatus(Element *e, unsigned ipnum) {
    NeuralNetworkMaterialStatus *newStatus = new NeuralNetworkMaterialStatus(this, e, ipnum); //needs to be deleted manually
    return newStatus;
};
