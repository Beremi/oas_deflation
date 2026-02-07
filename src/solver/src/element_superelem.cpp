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
Matrix MLMechElement :: readDataNormalizationMatrix(int size, fs :: path matrix_path) const {
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

    // std::cout << "\nmatrixColumnNumber\n" << matrixColumnNumber << "\n" << std::flush;
    // std::cout << "\nmatrixRowNumber\n" << matrixRowNumber << "\n" << std::flush;
    // std::cout << "\nmatrixEntries.size()\n" << matrixEntries.size() << "\n" << std::flush;

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
        } else if ( param.compare("norm_matL") == 0 ) {
            string filepath;
            iss >> filepath;
            nmL_path = GlobPaths :: BASEDIR  / filepath;
        }  else if ( param.compare("Ftype") == 0 ) {
            iss >> Ftype;
            if ( Ftype != "F" && Ftype != "Fp" ) {
                cerr << "MLElement ERROR:  incorrect value: " << Ftype << " for Ftype \n";
            }
        }  else if ( param.compare("Ktype") == 0 ) {
            iss >> Ktype;
            if ( Ktype != "L_tangent" && Ktype != "Lp_tangent" && Ktype != "K" && Ktype != "Kp" ) { // L_tangent, Lp_tangent, K, Kp
                cerr << "MLElement ERROR:  incorrect value: " << Ktype << " for Ktype \n";
            }
        }  else if ( param.compare("normalizationType") == 0 ) {
            iss >> normalizationType;
            if ( normalizationType != "normal" && normalizationType != "bisymLog" ) {
                cerr << "MLElement ERROR:  incorrect value: " << normalizationType << " for normalizationType \n";
            }
        }  else if ( param.compare("max_elastic_strain_energy") == 0 ) {
            iss >> max_elastic_strain_energy;
        }  else {
            cerr << "MLElement ERROR: " << param << " input parameter not defined \n";
        }
    }
}

//////////////////////////////////////////////////////////
void MLMechElement :: init() {
    Element :: init();

    // max_elastic_strain_energy = 1e4; //1e5

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
    stiffmat_elastic = A(keep_ind, keep_ind);
    stiffmat = stiffmat_elastic;

    if ( outDoFs != keep_ind.size() ) {
        cerr << "Error in MLMechElement: there are " << outDoFs << " DoFs on input of the element, but " << keep_ind.size() << " DoFs are required for polynom of degree " << poly_degree << endl;
        exit(1);
    }

    // load pytorch libtorch model
    module = torch :: jit :: load(ml_path.string() );

    // load normalization matrix
    int size  = keep_ind.size() - 3;
    int sizeL = size * ( size + 1 ) / 2;
    norm = readDataNormalizationMatrix(size, nm_path);
    if ( nmL_path.string().size() > 1 ) {
        normL = readDataNormalizationMatrix(sizeL, nmL_path);
    } else {
        normL = norm;
    }
}

//////////////////////////////////////////////////////////
Matrix MLMechElement :: giveStiffnessMatrix(std :: string matrixType) const {
    if ( matrixType == "elastic" ) {
        return stiffmat_elastic;
    } else {
        return stiffmat;
    }
    // cout << stiffmat << "\n";
}

//////////////////////////////////////////////////////////
Vector MLMechElement :: giveInternalForces() {
    if ( temp_frozen ) {
        return stiffmat_elastic * temp_DoFs;
    } else {
        Vector DoFs_rel = temp_DoFs.cast< double > ();
        double c0u = DoFs_rel [ 0 ];
        double c0v = DoFs_rel [ 1 ];

        // RIGID BODY MOTION
        for ( int i = 0; i < 8; i += 2 ) {
            DoFs_rel [ i ] -= c0u;
        }
        for ( int i = 1; i < 8; i += 2 ) {
            DoFs_rel [ i ] -= c0v;
        }

        double c1v = DoFs_rel [ 3 ];
        // RIGID BODY ROTATION
        DoFs_rel [ 3 ] += -c1v;
        DoFs_rel [ 4 ] += c1v;
        DoFs_rel [ 5 ] += -c1v;
        DoFs_rel [ 6 ] += c1v;


        // std::cout << "\nenergy:  " << 0.5 * DoFs_rel.transpose() * stiffmat_elastic * DoFs_rel << "\n" << std::flush;

        // always elastic if strain energy 1/2 u Ke uT < max_elastic_strain energy
        if ( 0.5 * DoFs_rel.transpose() * stiffmat_elastic * DoFs_rel <= max_elastic_strain_energy ) {
            Vector forces = stiffmat_elastic * DoFs_rel;
            stiffmat = stiffmat_elastic;
            // std::cout << "\nCheckpoint 00\n" << std::flush;
            return forces;
        } else {
            // std::cout << "\nChceckpoint 0\n" << std::flush;

            int size = DoFs.size();
            int size_x = size - 3; // size of network input

            // std::cout << "\nCheckpoint 00\n" << std::flush;
            // std::cout << "\nnorm Matrix\n" << norm << "\n" << std::flush;

            // WORKING WITH COUNTERCLOCK NUMBERING (order)

            // Normalization of input DoFs.
            Vector x_std = norm.row(1);
            Vector DoFs_norm(13);
            DoFs_norm(0) = DoFs_rel(2);
            DoFs_norm.tail(size_x - 1) = DoFs_rel.tail(size_x - 1);

            // std::cout << "\nCheckpoint 01\n" << std::flush;

            // bisymmetric log transformation
            if ( normalizationType == "bisymLog" ) {
                Vector Cx = norm.row(4);
                Eigen :: VectorXd vec1(size_x);
                vec1.fill(1);
                DoFs_norm = DoFs_norm.array().sign() * ( vec1.array() + ( DoFs_norm.array() / Cx.array() ).cwiseAbs() ).array().log10();
                DoFs_norm.matrix();
            }

            // std::cout << "\nCheckpoint 02\n" << std::flush;

            DoFs_norm -= norm.row(0);
            // DoFs_norm = DoFs_norm.array() / norm.row(1).array(); // doesnt work for some reason
            DoFs_norm = DoFs_norm.array() / x_std.array();


            // Create a vector of inputs.
            std :: vector< torch :: jit :: IValue >inputs;
            // inputs.push_back(torch::ones({1,24}));
            // Change DoFs type to float - Might be possible to avoid if a double type libtorch model is possible
            // Eigen::VectorXf DoFsFloat = DoFs.cast <float> ();
            Eigen :: VectorXf DoFsFloat = DoFs_norm.cast< float > ();
            // std::cout << "\nDofs\n" << DoFs << "\n";

            // std::cout << "\nDofsFloat\n" << DoFsFloat << "\n";

            // // // Create a torch Tensor populated by the DoFs values
            torch :: Tensor inputs_torch = torch :: from_blob(DoFsFloat.data(), { 1, size_x }).clone(); // Populates torch Tensor with Eigen Vector/Matrix
            // torch::Tensor inputs_torch = torch::ones({1, size_x}, torch::kFloat32); // for debugging

            inputs.push_back(inputs_torch);


            // //// NETWORK with (normalized L matrix, normalized forces) output
            // std::cout << "\nCheckpoint 01\n" << std::flush;

            auto outputs = module.forward(inputs);
            // std::cout << "\nCheckpoint 02\n" << std::flush;

            auto tuple_output = outputs.toTuple();
            torch :: Tensor forces_norm_tensor = tuple_output->elements() [ 1 ].toTensor();
            // std::cout << "\nCheckpoint 03\n" << std::flush;

            // Convert forces_norm Tensor to vector and to Eigen Vector and to double type
            std :: vector< float >outVec(forces_norm_tensor.data_ptr< float >(), forces_norm_tensor.data_ptr< float >() + forces_norm_tensor.numel() );
            Eigen :: VectorXf outFloat = Eigen :: VectorXf :: Map(& outVec [ 0 ], size);
            Eigen :: VectorXd forces_norm = outFloat.cast< double > ();
            // std::cout << "\nForces Norm\n" << forces_norm << "\n" << std::flush;

            // // Denormalize output forces
            Vector y_std = norm.row(3);
            Vector forces = forces_norm.array() * y_std.array();
            forces += norm.row(2);
            // std::cout << "\nForces Reduced CCorder\n" << forces << "\n" << std::flush;

            // bisymmetric log inverse transfomration
            if ( normalizationType == "bisymLog" ) {
                Vector Cy = norm.row(5);
                Eigen :: VectorXd vecm1(size_x);
                vecm1.fill(-1);
                Eigen :: VectorXd pow10(size_x);
                for (int i = 0; i < size_x; i++) {
                    pow10 [ i ] = pow( 10, fabs(forces [ i ]) );
                }
                forces = ( forces.array().sign().matrix().array() * Cy.array() ).array() * ( vecm1 + pow10 ).array();
            }


            // Compute back F0x, F0y, F1y (CCorder) from force and momentum equilibirium
            Eigen :: VectorXd forcesFull(size);
            forcesFull(2) = forces(0);       // Copy F1x
            forcesFull.tail(size_x - 1) = forces.tail(size_x - 1);    // Copy the rest of the forces
            forcesFull(3) = forces(1) - forces(2) + forces(3);   // Add F1y from M equlibirium
            forcesFull(0) = -( forcesFull(2) + forcesFull(4) + forcesFull(6) );  // Add F0x from Fx equlibirium
            forcesFull(1) = -( forcesFull(3) + forcesFull(5) + forcesFull(7) );  // Add F0y from Fy equlibirium
            // std::cout << "\nForces Full CCorder\n" << forcesFull << "\n" << std::flush;


            Matrix K;
            if ( Ktype.find('L') != std :: string :: npos ) {
                torch :: Tensor L_flat_norm_tensor = tuple_output->elements() [ 0 ].toTensor();

                // Convert L_norm Tensor to vector and to Eigen Vector and to double type
                std :: vector< float >outVec1(L_flat_norm_tensor.data_ptr< float >(), L_flat_norm_tensor.data_ptr< float >() + L_flat_norm_tensor.numel() );
                Eigen :: VectorXf outFloat1 = Eigen :: VectorXf :: Map(& outVec1 [ 0 ], size_x * ( size_x + 1 ) / 2); // n(n+1)/2 len of L_flat
                Eigen :: VectorXd L_flat_norm = outFloat1.cast< double > ();

                // Denormalize L flat
                Vector L_std = normL.row(1);
                Vector L_flat = L_flat_norm.array() * L_std.array();
                L_flat += normL.row(0);

                // std::cout << "\nCheckpoint 03\n" << std::flush;

                // Unflatten L
                Eigen :: MatrixXd L = Eigen :: MatrixXd :: Zero(size_x, size_x);
                int index = 0;
                for (int i = 0; i < size_x; ++i) {
                    for (int j = 0; j <= i; ++j) {
                        L(i, j) = L_flat(index++);
                    }
                }

                K = L * L.transpose();
            } else if ( Ktype.find('K') != std :: string :: npos ) {
                torch :: Tensor K_flat_norm_tensor = tuple_output->elements() [ 0 ].toTensor();

                // Convert K_norm Tensor to vector and to Eigen Vector and to double type
                std :: vector< float >outVec2(K_flat_norm_tensor.data_ptr< float >(), K_flat_norm_tensor.data_ptr< float >() + K_flat_norm_tensor.numel() );
                Eigen :: VectorXf outFloat2 = Eigen :: VectorXf :: Map(& outVec2 [ 0 ], size_x * size_x); // n(n+1)/2 len of L_flat
                Eigen :: VectorXd K_flat_norm = outFloat2.cast< double > ();

                // Unflatten K
                Matrix K_norm = Eigen :: Map< Eigen :: MatrixXd >(K_flat_norm.data(), size_x, size_x);

                // Denormalize K
                K = y_std.asDiagonal() * K_norm * x_std.cwiseInverse().asDiagonal();
            }

            std :: string filePathsL = "C:/Users/209050/OAS_data/Plasticity/ML/L_vals.txt";
            std :: ofstream outFilesL(filePathsL, std :: ios :: app);
            outFilesL << K << "\n";
            outFilesL.close();


            // calculate missing parts of matrix based on equilibrium equations
            Eigen :: MatrixXd Kfull(size, size);

            Kfull.block(3, 3, size_x, size_x) = K.block(0, 0, size_x, size_x); // block (start_x,start_y,len_x,len_y)
            Kfull.row(2).segment(3, size_x) =  Kfull.row(4).segment(3, size_x) - Kfull.row(5).segment(3, size_x) + Kfull.row(6).segment(3, size_x); //F1y
            Kfull.row(0).segment(3, size_x) =  -( Kfull.row(3).segment(3, size_x) + Kfull.row(4).segment(3, size_x) + Kfull.row(6).segment(3, size_x) ); //F0x
            Kfull.row(1).segment(3, size_x) =  -( Kfull.row(2).segment(3, size_x) + Kfull.row(5).segment(3, size_x) + Kfull.row(7).segment(3, size_x) ); //F0x
            Kfull.col(2) = Kfull.col(4) - Kfull.col(5) + Kfull.col(6);
            Kfull.col(0) = -( Kfull.col(3) + Kfull.col(4) + Kfull.col(6) );
            Kfull.col(1) = -( Kfull.col(2) + Kfull.col(5) + Kfull.col(7) );

            //swap the (2) and (3) columns and rows so that F0x and F0y are in the correct place
            Eigen :: VectorXd temp = Kfull.col(2);  // Store the 3rd column in a temporary vector
            Kfull.col(2) = Kfull.col(3);         // Set the 3rd column to the 4th column
            Kfull.col(3) = temp;                  // Set the 4th column to the original 3rd column (from temp)

            Eigen :: VectorXd temp_row = Kfull.row(2);  // Store the 3rd column in a temporary vector
            Kfull.row(2) = Kfull.row(3);         // Set the 3rd column to the 4th column
            Kfull.row(3) = temp_row;

            if ( Ktype.find('p') != std :: string :: npos ) {
                stiffmat = stiffmat_elastic - Kfull;
            } else {
                stiffmat = Kfull;
            }

            // std::cout << stiffmat << "\n" << std::flush;


            // // OLD

            if ( Ftype == "Fp" ) {
                return stiffmat_elastic * DoFs - forcesFull;
            } else {
                // std::cout << "direct Forces prediction" << "\n" << std::flush;

                return forcesFull;
            }

            // // NEW
            // return forcesFull;
            // return stiffmat * DoFs;
        }
        // return stiffmat*DoFs;
    }
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
void MLMechElement :: evaluateStrains(const Vector &DoFs) {
    temp_DoFs = DoFs;
}

//////////////////////////////////////////////////////////
void MLMechElement :: evaluateStresses(bool frozen, double timeStep) {
    temp_frozen = frozen;
    temp_timeStep = timeStep;
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
