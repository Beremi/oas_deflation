#include "material.h"
#include "element.h"
#include "material_neuralnetwork.h"
#include <torch/script.h>


using namespace std;


//////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////
// RECURRENT NETWORK TENSORIAL MECHANICAL MATERIAL

void NeuralNetworkMaterialStatus :: computeStress(double timeStep) {
    computeConstitutiveStrain();

    NeuralNetworkMaterial *m = static_cast< NeuralNetworkMaterial * >( mat );

    Vector input_norm = Eigen :: VectorXd :: Zero(10);

    input_norm.head(6) = temp_strain;
    input_norm(6) = m->giveE0();
    input_norm(7) = m->giveAlpha() ;
    input_norm(8) = m->giveft();
    input_norm(9) = m->giveRVEsize();
    
    if (m->give_normalize_inputs()) { // Normalization of input DoFs.
        input_norm.array() -= m->x_mean.array();
        input_norm.array() /= m->x_std.array();
    }

    Eigen :: VectorXf input_float = input_norm.cast< float > (); // Cast to float because libtorch works with float tensors by default
    torch :: Tensor inputs_torch = torch :: from_blob(input_float.data(), { 1, 10 }).clone(); // Populates torch Tensor with Eigen Vector/Matrix
    // inputs_torch = inputs_torch.unsqueeze(0); // so that sequence len = 1

    inputs_torch = inputs_torch.unsqueeze(0);

    // inputs into form that libtorch wants
    std :: vector< torch :: jit :: IValue >inputs;
    inputs.push_back(inputs_torch);

    // std::cout << "prev_strain" << prev_strain << "\n" << std::flush;
    // std::cout << "prev_strain size: " << prev_strain.size() << "\n";

    ///// strain increment - calculation, normalization, casting to torch tensor
    if (m->give_use_strain_increment()) {
        Vector strain_increment = temp_strain - prev_strain;
        if (m->give_normalize_inputs()) { // Normalization of increment.
            strain_increment.array() /= m->x_std.head(strain_increment.size()).array(); // !! only divide by std, mean is not subtracted because it cancels out in the increment. 
        }
        Eigen :: VectorXf increment_float = strain_increment.cast< float > ();
        torch :: Tensor increment_torch = torch :: from_blob(increment_float.data(), { 1, 6 }).clone(); // Populates torch Tensor with Eigen Vector/Matrix, 
        // !! .clone() is safer (and here necessary because strain_increment is destroyed after end of if)
        increment_torch = increment_torch.unsqueeze(0);
        inputs.push_back(increment_torch);
    }

    inputs.push_back(hidden_state);

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

    // std::cout << "hidden_state_prev \n"  << hidden_state << "\n";

    // Update hidden states
    temp_hidden_state = tuple_output->elements() [ 1 ];
    // std::cout << "hidden_state_post \n"  << temp_hidden_state << "\n";

    // std::cout << "post hidden update\n"  << "\n";

    // temp_hidden = hidden_new;

    stress_norm_tensor = stress_norm_tensor.squeeze(0); // - first squeeze to remove batch dimension, second to remove sequence dimension
    stress_norm_tensor = stress_norm_tensor.squeeze(0); // !! do not use _squeeze() - it is IN-PLACE and messes up the tensor memory address

    // std::cout << "\nCheckpoint 03\n" << std::flush;

    // Convert forces_norm Tensor to vector and to Eigen Vector and to double type
    auto stress_ptr = stress_norm_tensor.data_ptr<float>();
    Eigen::Map<Eigen::VectorXf> outFloat(stress_ptr, temp_strain.size());
    Eigen::VectorXd stress_norm = outFloat.cast<double>();
    // std::cout << "Stress Norm\n" << stress_norm << "\n" << std::flush;

    if (m->give_normalize_inputs()) { // Denormalization of output stress.
        stress_norm.array() *= m->y_std.array();
        stress_norm.array() += m->y_mean.array();
    }

    temp_stress = stress_norm; // denormalization happens in network

    // std::cout << "returning stress\n"  << "\n";
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
void NeuralNetworkMaterialStatus :: computeStressWithFrozenIntVars(double timeStep) {
    ( void ) timeStep;
    computeStress(timeStep);
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
    // for (size_t i = 0; i < hc_vectors.size(); ++i) {
    //     hc_vectors [ i ] = temp_hc_vectors [ i ];
    // }
    hidden_state = temp_hidden_state;
    prev_strain = temp_strain;
}

//////////////////////////////////////////////////////////
Matrix NeuralNetworkMaterialStatus :: giveMassTensor() const {
    NeuralNetworkMaterial *m = static_cast< NeuralNetworkMaterial * >( mat );
    unsigned dimension = mat->giveDimension();
    return Eigen :: MatrixXd :: Identity(dimension, dimension) * m->giveDensity();
}

//////////////////////////////////////////////////////////
// std :: vector< std :: vector< torch :: Tensor > >NeuralNetworkMaterialStatus :: giveHiddenState() const {
//     return temp_hc_vectors;
// }

//////////////////////////////////////////////////////////
NeuralNetworkMaterialStatus :: NeuralNetworkMaterialStatus(NeuralNetworkMaterial *m, Element *e, unsigned ipnum) : MaterialStatus(m, e, ipnum)  {
    name = "neural network tensorial mechanical mat. status";

    // prev_strain = Vector::Zero(temp_strain.size()); // doesn't work, needs to be done differently (with m -> giveDimension())
    prev_strain = Vector :: Zero(6);

    torch::jit::Module inner = m->giveNetwork().attr("model").toModule(); // actual network stored in wrapper
     if (!inner.find_method("init_hidden")) {
         cerr << name << ": 'init_hidden' function not found in the model, cannot initialize hidden states" << endl;
         exit(EXIT_FAILURE);
     }

    hidden_state = inner.run_method("init_hidden", 1);
    
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
            matrixEntries.push_back(stod(matrixEntry) );     //here we convert the string to double and fill in the row vector storing all the matrix entries
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
            matrixEntries.push_back( stod(matrixEntry) );     //here we convert the string to double and fill in the row vector storing all the matrix entries
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
    bool bstiffmat, bmlmodel, bnormmat, bE, bnu, bdensity, bE0, balpha, bft, bGt, bRVEsize;
    bstiffmat = bmlmodel = bnormmat = bE = bnu = bdensity = bE0 = balpha = bft = bGt = bRVEsize = false;

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
        } else if ( param.compare("alpha") == 0 ) {
            balpha = true,
            iss >> alpha;
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
                network = torch :: jit :: load(ml_path.string() );
            } catch ( const c10 :: Error &e ) {
                std :: cerr << "Error loading the model: " << e.what() << std :: endl;
                exit(EXIT_FAILURE);
            }

            if (filepath.find("SUB") != std::string::npos) {
                use_strain_increment = true;
            }

            torch::jit::Module inner = network.attr("model").toModule(); // actual network stored in wrapper
            
            if (!inner.find_method("get_use_strain_increment")) {
                cout << name << ": 'get_use_strain_increment' function not found in the model, looking for LMSC and SUB in name" << endl;
                std::string last_part = ml_path.filename().string(); // 
                if (last_part.find("LMSC") != std::string::npos) { // LMSC in model name - strain increment needed
                    cout << name << ": LMSC found in model name, using strain increment" << endl;
                    use_strain_increment = true;
                } else if (last_part.find("SUB") != std::string::npos) { // substepping in model - strain increment needed even for non LMSC networks
                    cout << name << ": SUB found in model name, using strain increment" << endl;
                    use_strain_increment = true;
                } else {
                    cout << name << ": LMSC or SUB not found in model name, assuming strain increment is not used" << endl;
                    use_strain_increment = false;
                }
            } else {
                use_strain_increment = inner.run_method("get_use_strain_increment").toBool();
            }

        } else if ( param.compare("norm_mat") == 0 ) {
            bnormmat = true;
            normalize_inputs = true;
            string filepath;
            iss >> filepath;
            nm_path = GlobPaths :: BASEDIR  / filepath;
            norm = readDataNormalizationMatrix(10, nm_path);
            x_mean = norm.row(0).transpose();
            x_std  = norm.row(1).transpose();
            y_mean = norm.row(2).head(6).transpose();
            y_std  = norm.row(3).head(6).transpose();
        }  else {
            cerr << "NeuralNetworkMaterial ERROR: " << param << " input parameter not defined \n";
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
    if ( !balpha ) {
        cerr << name << ": alpha was not specified" << endl;
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
        cerr << name << ": normalization matrix was not specified, assuming normalization is a part of the network" << endl; 
    }

};

//////////////////////////////////////////////////////////
MaterialStatus *NeuralNetworkMaterial :: giveNewMaterialStatus(Element *e, unsigned ipnum) {
    NeuralNetworkMaterialStatus *newStatus = new NeuralNetworkMaterialStatus(this, e, ipnum); //needs to be deleted manually
    return newStatus;
};
