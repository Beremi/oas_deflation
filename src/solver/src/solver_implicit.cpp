#include "solver_implicit.h"
#include "adaptivity.h"
#include <algorithm>
#include <array>
#include <cmath>
#include <fstream>
#include <limits>
#define numPhysicalFields 4

using namespace std;

//////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////
// STEADY STATE LINEAR SOLVER
SteadyStateLinearSolver :: SteadyStateLinearSolver() {
    name = "SteadyStateLinearSolver";
    conj_grad_precision = 1e-14;
    conj_grad_relative_maxit = 0.85;
    isTimeReal = false;
    stiffMatType = "elastic";
    stiffMatTypeFirstIT = "void";
    lastStiffMatType = stiffMatType;
}

//////////////////////////////////////////////////////////
SteadyStateLinearSolver :: ~SteadyStateLinearSolver() {
    //delete linalgsolver;
}

//////////////////////////////////////////////////////////
void SteadyStateLinearSolver :: prepareSystemMatricesAndInitialField(string init_r_file, string init_v_file, const bool initial) {
    ( void ) init_v_file;
    ( void ) init_r_file;    
    ( void ) initial;
    //nodes->addRHS_nodalLoad(load, 0); //to correctly account for abrupt initial change of BC
    //nodes->updateDirrichletBC(r, 0); //to correctly account for abrupt initial change of BC

    elems->prepareStiffnessMatrix(K);

    updateSystemMatrices(0, 0, 1);
}


//////////////////////////////////////////////////////////
void SteadyStateLinearSolver :: init(string init_r_file, string init_v_file, const bool initial) {
    Solver :: init(init_r_file, init_v_file, initial);
    prepareSystemMatricesAndInitialField(init_r_file, init_v_file, initial);
    computeKeff();
}

//////////////////////////////////////////////////////////
void SteadyStateLinearSolver :: rebuild() {
    Solver :: rebuild();
    prepareSystemMatricesAndInitialField("", "", false);
    computeKeff();
}

//////////////////////////////////////////////////////////
void SteadyStateLinearSolver :: collectLinearDeflationVector(const Vector &x, bool success) {
    if ( success && linalgsolver ) {
        linalgsolver->collectDeflationVector(x);
    }
}

//////////////////////////////////////////////////////////
ElasticDofMap SteadyStateLinearSolver :: buildElasticDofMap() const {
    ElasticDofMap map;
    if ( !nodes || freeDoFnum == 0 ) {
        return map;
    }

    unsigned dim = 0;
    std :: array< unsigned, numPhysicalFields >freeRowsByField {};
    unsigned unclassifiedRows = 0;
    for ( unsigned i = 0; i < freeDoFnum; i++ ) {
        const unsigned fullDoF = nodes->giveInvDoFid(i);
        const unsigned physicalField = nodes->givePhysicalFieldOfDoF(fullDoF);
        if ( physicalField < freeRowsByField.size() ) {
            freeRowsByField [ physicalField ]++;
        } else {
            unclassifiedRows++;
        }
        const Node *node = nodes->giveNodePointerOfDoFID(fullDoF);
        if ( node->doesMechanics() ) {
            dim = node->giveDimension();
        }
    }
    if ( dim != 2 && dim != 3 ) {
        return map;
    }

    unsigned blockSize = dim;
    std :: vector< const Node * >mechanicalNodes;
    std :: map< unsigned, unsigned >nodeIdToElasticBlock;
    for ( auto nodeIt = nodes->begin(); nodeIt != nodes->end(); ++nodeIt ) {
        const Node *node = *nodeIt;
        if ( !node->doesMechanics() ) {
            continue;
        }
        const unsigned nodeDofs = node->giveNumberOfDoFs();
        if ( dim == 2 && nodeDofs >= 3 ) {
            blockSize = 3;
        } else if ( dim == 3 && nodeDofs >= 6 ) {
            blockSize = 6;
        }
        nodeIdToElasticBlock [ node->giveID() ] = mechanicalNodes.size();
        mechanicalNodes.push_back(node);
    }

    if ( mechanicalNodes.empty() ) {
        return map;
    }

    map.dimension = dim;
    map.blockSize = blockSize;
    map.fullRows = mechanicalNodes.size() * blockSize;
    map.reducedToFull.assign(freeDoFnum, std :: numeric_limits< unsigned > :: max() );
    std :: vector< char >fullRowIsFree(map.fullRows, 0);

    unsigned translationalMechanicalRows = 0;
    unsigned rotationalMechanicalRows = 0;
    unsigned mappedRows = 0;
    unsigned unsupportedRows = 0;

    for ( unsigned i = 0; i < freeDoFnum; i++ ) {
        const unsigned fullDoF = nodes->giveInvDoFid(i);
        const Node *node = nodes->giveNodePointerOfDoFID(fullDoF);
        const unsigned relativeDoF = fullDoF - node->giveStartingDoF();
        if ( !node->doesMechanics() || relativeDoF >= blockSize ) {
            unsupportedRows++;
            continue;
        }
        const auto found = nodeIdToElasticBlock.find(node->giveID() );
        if ( found == nodeIdToElasticBlock.end() ) {
            unsupportedRows++;
            continue;
        }
        map.reducedToFull [ i ] = found->second * blockSize + relativeDoF;
        fullRowIsFree [ map.reducedToFull [ i ] ] = 1;
        mappedRows++;
        if ( relativeDoF < dim ) {
            translationalMechanicalRows++;
        } else {
            rotationalMechanicalRows++;
        }
    }

    if ( unsupportedRows > 0 || mappedRows != freeDoFnum ) {
        std :: cout << "AMG elasticity map unavailable: mapped_rows=" << mappedRows
                    << ", free_rows=" << freeDoFnum
                    << ", unsupported_rows=" << unsupportedRows << std :: endl;
        map.reducedToFull.clear();
        map.fullRows = 0;
        return map;
    }

    Point center = Point :: Zero();
    for ( const Node *node : mechanicalNodes ) {
        center += node->givePoint();
    }
    center /= double(mechanicalNodes.size() );

    double radius2 = 0.;
    map.coordinates.assign(size_t( mechanicalNodes.size() ) * dim, 0.);
    for ( unsigned nodeIndex = 0; nodeIndex < mechanicalNodes.size(); nodeIndex++ ) {
        Point point = mechanicalNodes [ nodeIndex ]->givePoint() - center;
        radius2 += point.squaredNorm();
        map.coordinates [ size_t(nodeIndex) * dim + 0 ] = point.x();
        map.coordinates [ size_t(nodeIndex) * dim + 1 ] = point.y();
        if ( dim == 3 ) {
            map.coordinates [ size_t(nodeIndex) * dim + 2 ] = point.z();
        }
    }

    const double scale = sqrt(radius2 / std :: max<size_t>(mechanicalNodes.size(), 1) );
    const bool hasRotationalDoFs = blockSize > dim;
    const double nearNullspaceScale = elasticNearNullspaceCoordinateScale > 0. ? elasticNearNullspaceCoordinateScale : scale;
    if ( nearNullspaceScale > 0. && ( !hasRotationalDoFs || elasticNearNullspaceCoordinateScale > 0. ) ) {
        for ( double &value : map.coordinates ) {
            value /= nearNullspaceScale;
        }
    }

    const int rawColumns = ( dim == 2 ) ? 3 : 6;
    std :: vector< double >raw(size_t( map.fullRows ) * rawColumns, 0.);
    for ( unsigned nodeIndex = 0; nodeIndex < mechanicalNodes.size(); nodeIndex++ ) {
        const double x = map.coordinates [ size_t(nodeIndex) * dim + 0 ];
        const double y = map.coordinates [ size_t(nodeIndex) * dim + 1 ];
        const double z = ( dim == 3 ) ? map.coordinates [ size_t(nodeIndex) * dim + 2 ] : 0.;
        for ( unsigned relativeDoF = 0; relativeDoF < blockSize; relativeDoF++ ) {
            const size_t row = size_t(nodeIndex) * blockSize + relativeDoF;
            if ( !fullRowIsFree [ row ] ) {
                continue;
            }
            if ( relativeDoF < dim ) {
                raw [ row * rawColumns + relativeDoF ] = 1.;
            }
            if ( dim == 2 ) {
                if ( relativeDoF == 0 ) {
                    raw [ row * rawColumns + 2 ] = -y;
                } else if ( relativeDoF == 1 ) {
                    raw [ row * rawColumns + 2 ] = x;
                } else if ( relativeDoF == 2 ) {
                    raw [ row * rawColumns + 2 ] = 1.;
                }
            } else {
                if ( relativeDoF == 0 ) {
                    raw [ row * rawColumns + 4 ] = z;
                    raw [ row * rawColumns + 5 ] = -y;
                } else if ( relativeDoF == 1 ) {
                    raw [ row * rawColumns + 3 ] = -z;
                    raw [ row * rawColumns + 5 ] = x;
                } else if ( relativeDoF == 2 ) {
                    raw [ row * rawColumns + 3 ] = y;
                    raw [ row * rawColumns + 4 ] = -x;
                } else if ( relativeDoF >= 3 && relativeDoF < 6 ) {
                    raw [ row * rawColumns + relativeDoF ] = 1.;
                }
            }
        }
    }

    int columns = 0;
    std :: vector< double >orthonormal(size_t( map.fullRows ) * rawColumns, 0.);
    for ( int col = 0; col < rawColumns; col++ ) {
        std :: vector< double >candidate(map.fullRows, 0.);
        for ( unsigned row = 0; row < map.fullRows; row++ ) {
            candidate [ row ] = raw [ size_t(row) * rawColumns + col ];
        }
        for ( int previous = 0; previous < columns; previous++ ) {
            double dot = 0.;
            for ( unsigned row = 0; row < map.fullRows; row++ ) {
                dot += candidate [ row ] * orthonormal [ size_t(row) * rawColumns + previous ];
            }
            for ( unsigned row = 0; row < map.fullRows; row++ ) {
                candidate [ row ] -= dot * orthonormal [ size_t(row) * rawColumns + previous ];
            }
        }

        double norm2 = 0.;
        for ( double value : candidate ) {
            norm2 += value * value;
        }
        const double norm = sqrt(norm2);
        if ( norm <= 1e-12 ) {
            continue;
        }
        for ( unsigned row = 0; row < map.fullRows; row++ ) {
            orthonormal [ size_t(row) * rawColumns + columns ] = candidate [ row ] / norm;
        }
        columns++;
    }

    if ( columns == 0 ) {
        map.reducedToFull.clear();
        map.fullRows = 0;
        return map;
    }

    map.nearNullspaceColumns = columns;
    map.nearNullspace.assign(size_t( map.fullRows ) * columns, 0.);
    for ( unsigned row = 0; row < map.fullRows; row++ ) {
        for ( int col = 0; col < columns; col++ ) {
            map.nearNullspace [ size_t(row) * columns + col ] = orthonormal [ size_t(row) * rawColumns + col ];
        }
    }

    std :: cout << "AMG elasticity map: full_rows=" << map.fullRows
                << ", reduced_rows=" << freeDoFnum
                << ", lifted_identity_rows=" << ( map.fullRows - mappedRows )
                << ", block_size=" << map.blockSize
                << ", dimension=" << map.dimension
                << ", centered_coordinate_scale=" << scale
                << ", near_nullspace_coordinate_scale="
                << ( elasticNearNullspaceCoordinateScale > 0. ? std :: to_string(elasticNearNullspaceCoordinateScale) : ( hasRotationalDoFs ? "physical" : "normalized" ) ) << std :: endl;
    std :: cout << "AMG elasticity map: prepared " << columns << " near-nullspace modes from "
                << translationalMechanicalRows << " translational and " << rotationalMechanicalRows
                << " rotational mechanical free DOF rows out of " << freeDoFnum << " free DOFs" << std :: endl;
    std :: cout << "AMG elasticity map: free reduced rows by physical field: mechanics=" << freeRowsByField [ 0 ]
                << ", transport=" << freeRowsByField [ 1 ]
                << ", thermal=" << freeRowsByField [ 2 ]
                << ", humidity=" << freeRowsByField [ 3 ]
                << ", unclassified=" << unclassifiedRows << std :: endl;
    return map;
}

//////////////////////////////////////////////////////////
void SteadyStateLinearSolver :: computeKeff() {
    Keff = K;
    if ( nodes->giveConstraints()->isActive() ) {
        nodes->giveConstraints()->transformToConstraintSpace(Keff);
    }
    factorizeLinearSystem();
}

//////////////////////////////////////////////////////////
void SteadyStateLinearSolver :: reset() {
    load.setZero();
    nodes->addRHS_nodalLoad(load, time); //add nodal load
    nodes->updateDirrichletBC(trial_r, time); //give prescribed DoFs
    computeForcesAtIntegrationTime(true);
    nodes->giveReducedForceArray(residuals, f);

    /*
     * if ( LinalgSymmetricSolver(Keff, ddr, f, ddr, conj_grad_precision, conj_grad_relative_maxit, symsolver_type) == false ) {
     *  std :: cerr << "Conjugate gradients did not converge, attempt to restart step" << endl;
     *  std :: cerr << "solver restart did not work" << endl;
     *  exit(1);
     * }
     */
    linalgsolver->solve(ddr, f);

    updateFieldVariables();
    computeForcesAtIntegrationTime(true);

    runAfterEachStep();
}

//////////////////////////////////////////////////////////
Solver *SteadyStateLinearSolver :: readFromFile(const string filename) {
    string param, line;
    bool bdt, bttime;
    bdt = bttime = false;
    ifstream inputfile( filename.c_str() );
    if ( inputfile.is_open() ) {
        while ( getline(inputfile >> std :: ws, line) ) {
            if ( line.empty() || ( line.at(0) == '#' ) ) {
                continue;
            }
            istringstream iss(line);
            iss >> param;
            if ( param.compare("time_step") == 0 ) {
                bdt = true;
                iss >> initdt;
                dt = initdt;
            } else if ( param.compare("total_time") == 0 ) {
                bttime = true;
                iss >> termination_time;
            } else if ( param.compare("conj_grad_precision") == 0 ) {
                iss >> conj_grad_precision;
            } else if ( param.compare("conj_grad_relative_maxit") == 0 ) {
                iss >> conj_grad_relative_maxit;
            } else if ( param.compare("init_time") == 0 ) {
                iss >> this->init_time;
            } else if ( param.compare("init_step") == 0 ) {
                iss >> this->init_step;
            } else if ( param.compare("solver_type") == 0 ) {
                iss >> symsolver_type;
            } else if ( param.compare("elastic_near_nullspace_coordinate_scale") == 0
                        || param.compare("amg_near_nullspace_coordinate_scale") == 0 ) {
                iss >> elasticNearNullspaceCoordinateScale;
            } else if ( param.compare("dfgmres_tolerance") == 0 ) {
                iss >> dfgmresOptions.tolerance;
            } else if ( param.compare("dfgmres_true_tolerance") == 0 ) {
                iss >> dfgmresOptions.trueTolerance;
            } else if ( param.compare("dfgmres_max_iterations") == 0 ) {
                iss >> dfgmresOptions.maxIterations;
            } else if ( param.compare("dfgmres_restart") == 0 ) {
                iss >> dfgmresOptions.restart;
            } else if ( param.compare("dfgmres_deflation_vectors") == 0 ) {
                iss >> dfgmresOptions.deflationVectors;
            } else if ( param.compare("dfgmres_deflation_eps") == 0 ) {
                iss >> dfgmresOptions.deflationEps;
            } else if ( param.compare("dfgmres_collect_newton_steps") == 0 ) {
                int value = 0;
                iss >> value;
                dfgmresOptions.collectNewtonSteps = ( value != 0 );
            } else if ( param.compare("dfgmres_preconditioner") == 0 ) {
                iss >> dfgmresOptions.preconditioner;
            } else if ( param.compare("dfgmres_reorthogonalize_on_matrix_change") == 0 ) {
                int value = 0;
                iss >> value;
                dfgmresOptions.reorthogonalizeOnMatrixChange = ( value != 0 );
            } else if ( param.compare("dfgmres_reorthogonalize") == 0 || param.compare("dfgmres_reorthogonalize_krylov") == 0 ) {
                int value = 0;
                iss >> value;
                dfgmresOptions.reorthogonalizeKrylov = ( value != 0 );
            } else if ( param.compare("dfgmres_elastic_reorder") == 0 ) {
                iss >> dfgmresOptions.elasticReorder;
            } else if ( param.compare("dfgmres_elastic_full_lift") == 0 ) {
                int value = 0;
                iss >> value;
                dfgmresOptions.elasticFullLift = ( value != 0 );
            } else if ( param.compare("dfgmres_verbose") == 0 ) {
                int value = 0;
                iss >> value;
                dfgmresOptions.verbose = ( value != 0 );
            } else if ( param.compare("hypre_tolerance") == 0 ) {
                iss >> hypreOptions.tolerance;
            } else if ( param.compare("hypre_max_iterations") == 0 ) {
                iss >> hypreOptions.maxIterations;
            } else if ( param.compare("hypre_coarsen_type") == 0 ) {
                iss >> hypreOptions.coarsenType;
            } else if ( param.compare("hypre_interp_type") == 0 ) {
                iss >> hypreOptions.interpType;
            } else if ( param.compare("hypre_strong_threshold") == 0 ) {
                iss >> hypreOptions.strongThreshold;
            } else if ( param.compare("hypre_nodal") == 0 ) {
                iss >> hypreOptions.nodal;
            } else if ( param.compare("hypre_nodal_diag") == 0 ) {
                iss >> hypreOptions.nodalDiag;
            } else if ( param.compare("hypre_num_functions") == 0 ) {
                iss >> hypreOptions.numFunctions;
            } else if ( param.compare("hypre_relax_type") == 0 ) {
                iss >> hypreOptions.relaxType;
            } else if ( param.compare("hypre_relax_order") == 0 ) {
                iss >> hypreOptions.relaxOrder;
            } else if ( param.compare("hypre_num_sweeps") == 0 || param.compare("hypre_sweeps") == 0 ) {
                iss >> hypreOptions.numSweeps;
            } else if ( param.compare("hypre_p_max") == 0 || param.compare("hypre_pmax") == 0 ) {
                iss >> hypreOptions.pMaxElmts;
            } else if ( param.compare("hypre_agg_levels") == 0 || param.compare("hypre_agg_num_levels") == 0 ) {
                iss >> hypreOptions.aggNumLevels;
            } else if ( param.compare("hypre_boomer_max_iterations") == 0 ) {
                iss >> hypreOptions.boomerMaxIterations;
            } else if ( param.compare("hypre_cheby_order") == 0 ) {
                iss >> hypreOptions.chebyOrder;
            } else if ( param.compare("hypre_cheby_fraction") == 0 ) {
                iss >> hypreOptions.chebyFraction;
            } else if ( param.compare("hypre_elastic_reorder") == 0 ) {
                iss >> hypreOptions.elasticReorder;
            } else if ( param.compare("hypre_print_level") == 0 ) {
                iss >> hypreOptions.printLevel;
            } else if ( param.compare("hypre_non_galerkin_tol") == 0 ) {
                iss >> hypreOptions.nonGalerkinTol;
            } else if ( param.compare("hypre_use_dof_functions") == 0 ) {
                int value = 0;
                iss >> value;
                hypreOptions.useDofFunctions = ( value != 0 );
            } else if ( param.compare("hypre_use_interp_vectors") == 0 ) {
                int value = 0;
                iss >> value;
                hypreOptions.useInterpVectors = ( value != 0 );
            } else if ( param.compare("hypre_interp_vec_variant") == 0 ) {
                iss >> hypreOptions.interpVecVariant;
            } else if ( param.compare("hypre_check_matrix") == 0 ) {
                int value = 0;
                iss >> value;
                hypreOptions.checkMatrix = ( value != 0 );
            } else if ( param.compare("hypre_threads") == 0 || param.compare("hypre_thread_count") == 0 ) {
                iss >> hypreOptions.threads;
            } else if ( param.compare("silent") == 0 ) {
                silent = true;
            } else if ( param.compare("pertrubation") == 0 ) {
                Pertrubation *p = new Pertrubation();
                p->readFromLine(iss);
                pertrubations.push_back(p);
            }
        }
        inputfile.close();
    }
    if ( !bdt ) {
        cerr << name << ": solver parameter 'time_step' was not specified" << endl;
        exit(EXIT_FAILURE);
    }
    ;
    if ( !bttime ) {
        cerr << name << ": solver parameter 'total_time' was not specified" << endl;
        exit(EXIT_FAILURE);
    }

    return this;
};

//////////////////////////////////////////////////////////
void SteadyStateLinearSolver :: solve() {
    nodes->addRHS_nodalLoad(load, time);  //add nodal load
    nodes->updateDirrichletBC(trial_r, time); //give prescribed DoFs
    updateFieldVariables();      //with ddr=0
    computeForcesAtIntegrationTime(true);

    nodes->giveReducedForceArray(residuals, f);
    //solve linear system
    /*
     *
     * if ( LinalgSymmetricSolver(Keff, ddr, f, ddr, conj_grad_precision, conj_grad_relative_maxit, symsolver_type) == false ) {
     *  terminated = true;
     *  cerr << "Conjugate gradients did not converge during initialization of solver" << endl;
     *  return;
     * }
     */

    linalgsolver->solve(ddr, f);

    /*
     * cout << "----- K ----" << endl;
     * Keff.print();
     * cout << "----- d ----" << endl;
     * for(auto p:ddr ) cout << p<< endl;
     * cout << "----- f ----" << endl;
     * for(auto p:f ) cout << p<< endl;
     */

    updateFieldVariables(); //calculate master fields at the step end
    computeForcesAtStepEnd(false); //to obtain the actual stress, fluxes, ...
}


//////////////////////////////////////////////////////////
bool SteadyStateLinearSolver :: updateSystemMatrices(unsigned iteration, unsigned cumul_iteration, bool enforce) {
    if ( enforce || ( iteration == 0 && ( stiffnessMatrixStepUpdate == 0 || ( stiffnessMatrixStepUpdate > 0 && step % stiffnessMatrixStepUpdate == 0 ) ) ) || stiffnessMatrixIterUpdate == 0 || ( stiffnessMatrixIterUpdate > 0 && iteration % stiffnessMatrixIterUpdate == 0 ) || ( stiffnessMatrixCumulIterUpdate > 0 && cumul_iteration % stiffnessMatrixCumulIterUpdate == 0 ) ) {
        std :: string matrixType = stiffMatType;
        if ( iteration == 0 && stiffMatTypeFirstIT.compare("void") != 0 ) {
            matrixType = stiffMatTypeFirstIT;
        }
        lastStiffMatType = matrixType;
        elems->updateStiffnessMatrix(K, matrixType);
        if ( stiffMatElasticBlendBeta > 0. && matrixType.compare("elastic") != 0 ) {
            CoordinateIndexedSparseMatrix elasticK(K);
            elems->updateStiffnessMatrix(elasticK, "elastic");
            K = ( 1. - stiffMatElasticBlendBeta ) * K + stiffMatElasticBlendBeta * elasticK;
            K.makeCompressed();
        }
        return true;
    } else {
        return false;
    }
}


//////////////////////////////////////////////////////////
void SteadyStateLinearSolver :: runBeforeEachStep() {
    Solver :: runBeforeEachStep();
    trial_r = r;
    if ( not silent ) {
        cout << "######### Solving step " << step << " at time " << time << "; time step " << dt << " #########" << endl;
    }
}

//////////////////////////////////////////////////////////
void SteadyStateLinearSolver :: runAfterEachStep() {
    Solver :: runAfterEachStep();
}


//////////////////////////////////////////////////////////
void SteadyStateLinearSolver :: factorizeLinearSystem() {
    if ( not silent ) {
        cout << "factorizing system matrix" << endl;
    }

    if ( linalgsolver == nullptr ) {
        if ( symsolver_type == "EigenConj" ) {
            std :: unique_ptr< ConjGradSolver >cgs = std :: make_unique< ConjGradSolver >();
            cgs->setPrecisionAndRelMaxIters(conj_grad_precision, conj_grad_relative_maxit);
            linalgsolver = std :: move(cgs);
            linalgsolver->analyzePattern(Keff);
        } else if  ( symsolver_type == "EigenLDLT" ) {
            linalgsolver = std :: make_unique< LDLTSolver >();
            linalgsolver->analyzePattern(Keff);
        } else if  ( symsolver_type == "EigenLLT" ) {
            linalgsolver = std :: make_unique< LLTSolver >();
            linalgsolver->analyzePattern(Keff);
        } else if  ( symsolver_type == "EigenLU" || symsolver_type == "EigenSparseLU" ) {
            linalgsolver = std :: make_unique< LUSolver >();
            linalgsolver->analyzePattern(Keff);
#ifdef SUPERLU_FOUND
        } else if  ( symsolver_type == "SuperLU" ) {
            linalgsolver = std :: make_unique< SuperLUSolver >();
            linalgsolver->analyzePattern(Keff);
#endif
#ifdef PARDISO_FOUND
        } else if  ( symsolver_type == "PardisoLU" ) {
            linalgsolver = std :: make_unique< PardisoLUSolver >();
            linalgsolver->analyzePattern(Keff);
        } else if  ( symsolver_type == "PardisoLDLT" ) {
            linalgsolver = std :: make_unique< PardisoLDLTSolver >();
            linalgsolver->analyzePattern(Keff);
        } else if  ( symsolver_type == "PardisoLLT" ) {
            linalgsolver = std :: make_unique< PardisoLLTSolver >();
            linalgsolver->analyzePattern(Keff);
#endif
#ifdef CHOLMOD_FOUND
        } else if  ( symsolver_type == "CholmodLLT" ) {
            linalgsolver = std :: make_unique< CholmodLLTSolver >();
            linalgsolver->analyzePattern(Keff);
        } else if  ( symsolver_type == "CholmodLDLT" ) {
            linalgsolver = std :: make_unique< CholmodLDLTSolver >();
            linalgsolver->analyzePattern(Keff);
        } else if  ( symsolver_type == "CholmodSupernodalLLT" ) {
            linalgsolver = std :: make_unique< CholmodSupernodalLLTSolver >();
            linalgsolver->analyzePattern(Keff);
#endif
#ifdef HYPRE_FOUND
        } else if  ( symsolver_type == "DeflatedFGMRES" || symsolver_type == "DFGMRES" || symsolver_type == "DeflatedFlexibleGMRES" ) {
            std :: unique_ptr< DeflatedFGMRESSolver >dfgmres = std :: make_unique< DeflatedFGMRESSolver >();
            dfgmres->setOptions(dfgmresOptions);
            dfgmres->setHypreOptions(hypreOptions);
            dfgmres->setElasticDofMap(buildElasticDofMap() );
            linalgsolver = std :: move(dfgmres);
            linalgsolver->analyzePattern(Keff);
#endif
        } else {
            cerr << "Solver type " << symsolver_type << " is not implemented" << endl;
            exit(1);
        }
    }
    linalgsolver->factorize(Keff);
}

//////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////
// STEADY STATE NONLINEAR SOLVER

SteadyStateNonLinearSolver :: SteadyStateNonLinearSolver() {
    name = "SteadyStateNonLinearSolver";
    idc = nullptr;

    EPS2.resize(4);
    EPS2 [ 0 ] = 1e-20; //mechanics
    EPS2 [ 1 ] = 1e-10; //transport
    EPS2 [ 2 ] = 1e-17; //temperature
    EPS2 [ 3 ] = 1e-18; //humidity

    maxIt = 30;
    minIt = 1;
    enlargeIt = shortenIt = 0;
    maxDisErr = maxResErr = maxEneErr = 1e-5;
    limitEneErr = limitResErr = limitDisErr = 0;
    step_increase = 1.25;
    step_decrease = 0.8;
    critical_step_decrease = 0.5;
    stiffnessMatrixIterUpdate = -1;
    stiffnessMatrixCumulIterUpdate = -1;
    stiffnessMatrixStepUpdate = -1;

    it = 0;
    restarts = 0;
    cumul_it = 0;
    stiffMatType = "secant";
    lastStiffMatType = stiffMatType;

    //eigen error fields
    eigen_trial_rPF = Vector :: Zero(numPhysicalFields);
    eigen_f_extPF = Vector :: Zero(numPhysicalFields);
    eigen_WextPF = Vector :: Zero(numPhysicalFields);
}

//////////////////////////////////////////////////////////
SteadyStateNonLinearSolver :: ~SteadyStateNonLinearSolver() {
    if ( idc ) {
        delete idc;
    }
}

//////////////////////////////////////////////////////////
void SteadyStateNonLinearSolver :: init(string init_r_file, string init_v_file, const bool initial) {
    elems->readMatStatsFromFile(this->init_time, this->init_step, this->initdt, this->init_idc_time);
    SteadyStateLinearSolver :: init(init_r_file, init_v_file, initial);
    this->fully_converged = false;
    if ( idc ) {
        idc->init(nodes, funcs, initial);   //indirect displacement control
        ddf = Vector :: Zero(freeDoFnum);
        full_ddf = Vector :: Zero(totalDoFnum);
        f_last_iter = Vector :: Zero(freeDoFnum);
        if ( initial ) {
            idc_time = this->init_idc_time;
            idc_dt = 1e-6;
            idc_time_converged = this->init_idc_time;
        } else {
            // JK: during solver initialization after geometry update, idc_time must be set to time of previously converged step
            idc_time = idc_time_converged;
        }
    }

    computeForcesAtIntegrationTime(true); //to initialize all fields in the model
}

//////////////////////////////////////////////////////////
Solver *SteadyStateNonLinearSolver :: readFromFile(const string filename) {
    SteadyStateLinearSolver :: readFromFile(filename);

    string param, line;
    dtmax = dtmin = dt;
    bool bdtmin = false;
    bool bdtmax = false;
    bool ben = false;
    bool bsh = false;
    unsigned helpuint;
    double valueIN;
    ifstream inputfile(filename.c_str(), ios :: in | ios :: binary);
    if ( inputfile.is_open() ) {
        while ( getline(inputfile >> std :: ws, line) ) {
            if ( line.empty() || ( line.at(0) == '#' ) ) {
                continue;
            }
            istringstream iss(line);
            iss >> param;
            if ( param.compare("max_time_step") == 0 ) {
                bdtmax = true;
                iss >> dtmax;
            } else if ( param.compare("min_time_step") == 0 ) {
                bdtmin = true;
                iss >> dtmin;
            } else if ( param.compare("tolerance") == 0 ) {
                iss >> maxDisErr;
                maxResErr = maxEneErr = maxDisErr;
            } else if ( param.compare("tolerance_residuals") == 0 ) {
                iss >> maxResErr;
            } else if ( param.compare("tolerance_energies") == 0 ) {
                iss >> maxEneErr;
            } else if ( param.compare("tolerance_increments") == 0 ) {
                iss >> maxDisErr;
            } else if ( param.compare("stiffness_matrix_update") == 0 || param.compare("stiffness_matrix_iter_update") == 0  ) {
                iss >> valueIN;
                stiffnessMatrixIterUpdate = int( valueIN );
            } else if ( param.compare("stiffness_matrix_cumul_iter_update") == 0  ) {
                iss >> valueIN;
                stiffnessMatrixCumulIterUpdate = int( valueIN );
            } else if ( param.compare("stiffness_matrix_step_update") == 0  ) {
                iss >> valueIN;
                stiffnessMatrixStepUpdate = int( valueIN );
            } else if ( param.compare("limit_tolerance") == 0 ) {
                iss >> valueIN;
                limitEneErr = limitResErr = limitDisErr = valueIN;
            } else if ( param.compare("minIt") == 0 || param.compare("min_iterations") == 0 ) {
                iss >> minIt;
            } else if ( param.compare("maxIt") == 0 || param.compare("max_iterations") == 0 ) {
                iss >> maxIt;
                if ( maxIt < 1 ) {
                    std :: cerr << "number of itteration cannot be smaller than " << maxIt << "!!!,";
                    maxIt = 30;
                    std :: cout << " setting to default value: maxIt = " << maxIt << '\n';
                } else if ( maxIt < 3 ) {
                    std :: cout << "solver parameter maxIt set to " << maxIt << ", be carefull with such a small number" << '\n';
                }
            } else if ( param.compare("enlargeIt") == 0 ) {
                iss >> enlargeIt;
                ben = true;
            } else if ( param.compare("shortenIt") == 0 ) {
                iss >> shortenIt;
                bsh = true;
            } else if ( param.compare("step_increase") == 0 ) {
                iss >> valueIN;
                if ( valueIN < 1 ) {
                    std :: cerr << "step_increase cannot be smaller than 1! leaving default value " << step_increase << '\n';
                } else {
                    step_increase = valueIN;
                }
            } else if ( param.compare("step_decrease") == 0 ) {
                iss >> valueIN;
                if ( valueIN > 1 ) {
                    std :: cerr << "step_decrease cannot be greater than 1! leaving default value " << step_decrease << '\n';
                } else {
                    step_decrease = valueIN;
                }
            } else if ( param.compare("critical_step_decrease") == 0 ) {
                iss >> valueIN;
                if ( valueIN > 1 ) {
                    std :: cerr << "critical_step_decrease cannot be greater than 1! leaving default value " << critical_step_decrease << '\n';
                } else {
                    critical_step_decrease = valueIN;
                }
            } else if ( param.compare("nonlinear_damping_type") == 0 ) {
                iss >> param;
                if ( param.compare("off") == 0 ) {
                    nonlinearDampingType = NonlinearDampingType :: Off;
                } else if ( param.compare("fixed") == 0 ) {
                    nonlinearDampingType = NonlinearDampingType :: Fixed;
                } else if ( param.compare("adaptive") == 0 ) {
                    nonlinearDampingType = NonlinearDampingType :: Adaptive;
                } else if ( param.compare("rollback_adaptive") == 0 || param.compare("adaptive_rollback") == 0 ) {
                    nonlinearDampingType = NonlinearDampingType :: RollbackAdaptive;
                } else {
                    std :: cerr << "unknown nonlinear_damping_type '" << param << "', using off" << '\n';
                    nonlinearDampingType = NonlinearDampingType :: Off;
                }
            } else if ( param.compare("nonlinear_damping_factor") == 0 ) {
                iss >> nonlinearDampingFactor;
            } else if ( param.compare("nonlinear_damping_min") == 0 ) {
                iss >> nonlinearDampingMin;
            } else if ( param.compare("nonlinear_damping_max") == 0 ) {
                iss >> nonlinearDampingMax;
            } else if ( param.compare("nonlinear_damping_increase") == 0 ) {
                iss >> nonlinearDampingIncrease;
            } else if ( param.compare("nonlinear_damping_decrease") == 0 ) {
                iss >> nonlinearDampingDecrease;
            } else if ( param.compare("nonlinear_line_search") == 0 ) {
                iss >> param;
                if ( param.compare("off") == 0 ) {
                    nonlinearLineSearchType = NonlinearLineSearchType :: Off;
                } else if ( param.compare("backtracking") == 0 ) {
                    nonlinearLineSearchType = NonlinearLineSearchType :: Backtracking;
                } else if ( param.compare("bisection") == 0 || param.compare("full_bisection") == 0 ) {
                    nonlinearLineSearchType = NonlinearLineSearchType :: Bisection;
                } else {
                    std :: cerr << "unknown nonlinear_line_search '" << param << "', using off" << '\n';
                    nonlinearLineSearchType = NonlinearLineSearchType :: Off;
                }
            } else if ( param.compare("nonlinear_line_search_merit") == 0 ) {
                iss >> param;
                if ( param.compare("residual") == 0 ) {
                    nonlinearMeritType = NonlinearMeritType :: Residual;
                } else if ( param.compare("energy") == 0 ) {
                    nonlinearMeritType = NonlinearMeritType :: Energy;
                } else if ( param.compare("mixed") == 0 ) {
                    nonlinearMeritType = NonlinearMeritType :: Mixed;
                } else {
                    std :: cerr << "unknown nonlinear_line_search_merit '" << param << "', using mixed" << '\n';
                    nonlinearMeritType = NonlinearMeritType :: Mixed;
                }
            } else if ( param.compare("nonlinear_line_search_reduction") == 0 ) {
                iss >> nonlinearLineSearchReduction;
            } else if ( param.compare("nonlinear_line_search_min_alpha") == 0 ) {
                iss >> nonlinearLineSearchMinAlpha;
            } else if ( param.compare("nonlinear_line_search_max_alpha") == 0 ) {
                iss >> nonlinearLineSearchMaxAlpha;
            } else if ( param.compare("nonlinear_line_search_bisection_tolerance") == 0 ) {
                iss >> nonlinearLineSearchBisectionTolerance;
            } else if ( param.compare("nonlinear_line_search_max_trials") == 0 ) {
                iss >> nonlinearLineSearchMaxTrials;
            } else if ( param.compare("nonlinear_line_search_armijo") == 0 ) {
                iss >> nonlinearLineSearchArmijo;
            } else if ( param.compare("nonlinear_line_search_accept_any_decrease") == 0 ) {
                int value = 0;
                iss >> value;
                nonlinearLineSearchAcceptAnyDecrease = ( value != 0 );
            } else if ( param.compare("nonlinear_line_search_freeze_material") == 0 ) {
                int value = 0;
                iss >> value;
                nonlinearLineSearchFreezeMaterial = ( value != 0 );
                nonlinearLineSearchEvaluationMode = nonlinearLineSearchFreezeMaterial ? NonlinearLineSearchEvaluationMode :: Frozen : NonlinearLineSearchEvaluationMode :: Actual;
            } else if ( param.compare("nonlinear_line_search_evaluation") == 0 || param.compare("nonlinear_line_search_trial_evaluation") == 0 ) {
                iss >> param;
                if ( param.compare("frozen") == 0 ) {
                    nonlinearLineSearchEvaluationMode = NonlinearLineSearchEvaluationMode :: Frozen;
                    nonlinearLineSearchFreezeMaterial = true;
                } else if ( param.compare("actual") == 0 ) {
                    nonlinearLineSearchEvaluationMode = NonlinearLineSearchEvaluationMode :: Actual;
                    nonlinearLineSearchFreezeMaterial = false;
                } else if ( param.compare("frozen_then_actual") == 0 || param.compare("frozen-then-actual") == 0 ) {
                    nonlinearLineSearchEvaluationMode = NonlinearLineSearchEvaluationMode :: FrozenThenActual;
                    nonlinearLineSearchFreezeMaterial = true;
                } else {
                    std :: cerr << "unknown nonlinear_line_search_evaluation '" << param << "', using frozen" << '\n';
                    nonlinearLineSearchEvaluationMode = NonlinearLineSearchEvaluationMode :: Frozen;
                    nonlinearLineSearchFreezeMaterial = true;
                }
            } else if ( param.compare("nonlinear_line_search_cutback_on_fail") == 0 ) {
                int value = 0;
                iss >> value;
                nonlinearLineSearchCutbackOnFail = ( value != 0 );
            } else if ( param.compare("nonlinear_line_search_report_trials") == 0 ) {
                int value = 0;
                iss >> value;
                nonlinearLineSearchReportTrials = ( value != 0 );
            } else if ( param.compare("nonlinear_stagnation_cutback") == 0 ) {
                int value = 0;
                iss >> value;
                nonlinearStagnationCutback = ( value != 0 );
            } else if ( param.compare("nonlinear_stagnation_iterations") == 0 ) {
                iss >> nonlinearStagnationIterations;
            } else if ( param.compare("nonlinear_stagnation_ratio") == 0 ) {
                iss >> nonlinearStagnationRatio;
            } else if ( param.compare("nonlinear_growth_cutback") == 0 ) {
                iss >> nonlinearGrowthCutback;
            } else if ( param.compare("nonlinear_adaptive_matrix_update") == 0 ) {
                int value = 0;
                iss >> value;
                nonlinearAdaptiveMatrixUpdate = ( value != 0 );
            } else if ( param.compare("nonlinear_rebuild_on_small_alpha") == 0 ) {
                iss >> nonlinearRebuildOnSmallAlpha;
            } else if ( param.compare("nonlinear_rebuild_on_merit_growth") == 0 ) {
                int value = 0;
                iss >> value;
                nonlinearRebuildOnMeritGrowth = ( value != 0 );
            } else if ( param.compare("nonlinear_rebuild_on_stagnation") == 0 ) {
                int value = 0;
                iss >> value;
                nonlinearRebuildOnStagnation = ( value != 0 );
            } else if ( param.compare("nonlinear_trust_region") == 0 ) {
                iss >> param;
                if ( param.compare("off") == 0 ) {
                    nonlinearTrustRegionType = NonlinearTrustRegionType :: Off;
                } else if ( param.compare("step_norm") == 0 || param.compare("step-norm") == 0 ) {
                    nonlinearTrustRegionType = NonlinearTrustRegionType :: StepNorm;
                } else {
                    std :: cerr << "unknown nonlinear_trust_region '" << param << "', using off" << '\n';
                    nonlinearTrustRegionType = NonlinearTrustRegionType :: Off;
                }
            } else if ( param.compare("nonlinear_trust_radius_initial") == 0 ) {
                iss >> nonlinearTrustRadiusInitial;
            } else if ( param.compare("nonlinear_trust_radius_min") == 0 ) {
                iss >> nonlinearTrustRadiusMin;
            } else if ( param.compare("nonlinear_trust_radius_max") == 0 ) {
                iss >> nonlinearTrustRadiusMax;
            } else if ( param.compare("nonlinear_trust_shrink") == 0 ) {
                iss >> nonlinearTrustShrink;
            } else if ( param.compare("nonlinear_trust_expand") == 0 ) {
                iss >> nonlinearTrustExpand;
            } else if ( param.compare("nonlinear_trust_max_trials") == 0 ) {
                iss >> nonlinearTrustMaxTrials;
            } else if ( param.compare("nonlinear_rollback_ignore_initial_iterations") == 0 ) {
                iss >> nonlinearRollbackIgnoreInitialIterations;
            } else if ( param.compare("nonlinear_rollback_growth_ratio") == 0 ) {
                iss >> nonlinearRollbackGrowthRatio;
            } else if ( param.compare("nonlinear_rollback_increase_ratio") == 0 ) {
                iss >> nonlinearRollbackIncreaseRatio;
            } else if ( param.compare("nonlinear_rollback_enable") == 0 ) {
                int value = 0;
                iss >> value;
                nonlinearRollbackEnable = ( value != 0 );
            } else if ( param.compare("nonlinear_rollback_decrease_on_any_increase") == 0 ) {
                int value = 0;
                iss >> value;
                nonlinearRollbackDecreaseOnAnyIncrease = ( value != 0 );
            } else if ( param.compare("nonlinear_rollback_max_growth_decreases") == 0 ) {
                iss >> nonlinearRollbackMaxGrowthDecreases;
            } else if ( param.compare("nonlinear_material_snapshot_rollback") == 0 ) {
                int value = 0;
                iss >> value;
                nonlinearMaterialSnapshotRollback = ( value != 0 );
            } else if ( param.compare("nonlinear_material_snapshot_verify") == 0 ) {
                int value = 0;
                iss >> value;
                nonlinearMaterialSnapshotVerify = ( value != 0 );
            } else if ( param.compare("nonlinear_tangent_check") == 0 ) {
                int value = 0;
                iss >> value;
                nonlinearTangentCheck = ( value != 0 );
            } else if ( param.compare("nonlinear_tangent_check_step") == 0 ) {
                iss >> nonlinearTangentCheckStep;
            } else if ( param.compare("nonlinear_tangent_check_iteration") == 0 ) {
                iss >> nonlinearTangentCheckIteration;
            } else if ( param.compare("nonlinear_tangent_check_eps") == 0 ) {
                iss >> nonlinearTangentCheckEps;
            } else if ( param.compare("nonlinear_tangent_check_random_vectors") == 0 ) {
                iss >> nonlinearTangentCheckRandomVectors;
            } else if ( param.compare("nonlinear_tangent_check_include_newton") == 0 ) {
                int value = 0;
                iss >> value;
                nonlinearTangentCheckIncludeNewton = ( value != 0 );
            } else if ( param.compare("nonlinear_tangent_check_stop_after") == 0 ) {
                int value = 0;
                iss >> value;
                nonlinearTangentCheckStopAfter = ( value != 0 );
            } else if ( param.compare("nonlinear_tangent_check_output") == 0 ) {
                iss >> nonlinearTangentCheckOutput;
            } else if ( param.compare("nonlinear_tangent_check_scope") == 0 ) {
                std :: string value;
                iss >> value;
                if ( value.compare("global") == 0 ) {
                    nonlinearTangentCheckScope = NonlinearTangentCheckScope :: Global;
                } else if ( value.compare("element_top") == 0 ) {
                    nonlinearTangentCheckScope = NonlinearTangentCheckScope :: ElementTop;
                } else {
                    std :: cerr << "unknown nonlinear_tangent_check_scope '" << value << "', using global" << '\n';
                    nonlinearTangentCheckScope = NonlinearTangentCheckScope :: Global;
                }
            } else if ( param.compare("nonlinear_tangent_check_top_elements") == 0 ) {
                iss >> nonlinearTangentCheckTopElements;
            } else if ( param.compare("nonlinear_tangent_check_element_output") == 0 ) {
                iss >> nonlinearTangentCheckElementOutput;
            } else if ( param.compare("nonlinear_tangent_check_matrix_type") == 0 ) {
                iss >> nonlinearTangentCheckMatrixType;
            } else if ( param.compare("indirect_displacement_control") == 0 || param.compare("indirect_control") == 0 ) {
                iss >> helpuint;
                if ( !idc ) {
                    idc = new IndirectControl();
                }
                idc->readFromStream(helpuint, inputfile);
            } else if ( param.compare("indirect_control_square_sum") == 0 ) {
                iss >> helpuint;
                if ( !idc ) {
                    idc = new IndirectControlSumOfSquares();
                }
                idc->readFromStream(helpuint, inputfile);
            } else if ( param.compare("stiff_matrix_type") == 0 ) {
                iss >> stiffMatType;
                if ( stiffMatType.compare("elastic") != 0 && stiffMatType.compare("secant") != 0  && stiffMatType.compare("tangent") != 0  && stiffMatType.compare("consistent") != 0 ) {
                    cerr << "Error: stiff_matrix_type must be 'elastic', 'secant', 'tangent', or 'consistent', entered value is " << stiffMatType << endl;
                    exit(1);
                }
            } else if ( param.compare("stiff_matrix_elastic_blend_beta") == 0 ) {
                iss >> stiffMatElasticBlendBeta;
            } else if ( param.compare("first_iteration_stiff_matrix_type") == 0 ) {
                iss >> stiffMatTypeFirstIT;
                if ( stiffMatTypeFirstIT.compare("elastic") != 0 && stiffMatTypeFirstIT.compare("secant") != 0  && stiffMatTypeFirstIT.compare("tangent") != 0  && stiffMatTypeFirstIT.compare("consistent") != 0 ) {
                    cerr << "Error: first_iteration_stiff_matrix_type must be 'elastic', 'secant', 'tangent', or 'consistent', entered value is " << stiffMatTypeFirstIT << endl;
                    exit(1);
                }
            }
        }
        inputfile.close();
    }
    if ( stiffMatElasticBlendBeta < 0. || stiffMatElasticBlendBeta > 1. ) {
        std :: cerr << "stiff_matrix_elastic_blend_beta must be in [0, 1], clamping" << std :: endl;
        stiffMatElasticBlendBeta = std :: min(1., std :: max(0., stiffMatElasticBlendBeta) );
    }
    if ( !bdtmin && !bdtmax ) {
        cout << "fixed time step used" << endl;
    } else if ( bdtmin || bdtmax ) {
        cout << "adaptive time step used" << endl;
        if ( dtmin > dtmax ) {
            std :: cerr << "dtmin > dtmax, swapping values" << '\n';
            double dm = dtmax;
            dtmax = dtmin;
            dtmin = dm;
        }
        //Jan E: NO, let's tolerate that first step can be out of bounds
        //if ( dt < dtmin ) {
        //    std :: cerr << "dt < dtmin, setting dtmin = dt" << '\n';
        //    dtmin = dt;
        //} else if ( dt > dtmax ) {
        //    std :: cerr << "dt > dtmax, setting dtmax = dt" << '\n';
        //    dtmax = dt;
        //}
    } else {
        cout << endl;
    }
    if ( enlargeIt > shortenIt ) {
        std :: cerr << "cannot set number of iterations for step enlargement higher than number of iterations for step shortening, setting back to default values" << '\n';
        ben = false;
        bsh = false;
    }
    if ( nonlinearDampingMin <= 0. ) {
        std :: cerr << "nonlinear_damping_min must be positive, setting to 0.03125" << '\n';
        nonlinearDampingMin = 0.03125;
    }
    if ( nonlinearDampingMax < nonlinearDampingMin ) {
        std :: cerr << "nonlinear_damping_max smaller than nonlinear_damping_min, setting max to min" << '\n';
        nonlinearDampingMax = nonlinearDampingMin;
    }
    if ( nonlinearDampingFactor <= 0. ) {
        std :: cerr << "nonlinear_damping_factor must be positive, setting to 1" << '\n';
        nonlinearDampingFactor = 1.;
    }
    if ( nonlinearDampingIncrease < 1. ) {
        std :: cerr << "nonlinear_damping_increase must be at least 1, setting to 1" << '\n';
        nonlinearDampingIncrease = 1.;
    }
    if ( nonlinearDampingDecrease <= 0. || nonlinearDampingDecrease > 1. ) {
        std :: cerr << "nonlinear_damping_decrease must be in (0, 1], setting to 0.5" << '\n';
        nonlinearDampingDecrease = 0.5;
    }
    if ( nonlinearLineSearchReduction <= 0. || nonlinearLineSearchReduction >= 1. ) {
        std :: cerr << "nonlinear_line_search_reduction must be in (0, 1), setting to 0.5" << '\n';
        nonlinearLineSearchReduction = 0.5;
    }
    if ( nonlinearLineSearchMinAlpha <= 0. ) {
        std :: cerr << "nonlinear_line_search_min_alpha must be positive, setting to 0.03125" << '\n';
        nonlinearLineSearchMinAlpha = 0.03125;
    }
    if ( nonlinearLineSearchMaxAlpha <= 0. ) {
        std :: cerr << "nonlinear_line_search_max_alpha must be positive, setting to 1" << '\n';
        nonlinearLineSearchMaxAlpha = 1.;
    }
    if ( nonlinearLineSearchBisectionTolerance <= 0. ) {
        std :: cerr << "nonlinear_line_search_bisection_tolerance must be positive, setting to 1e-3" << '\n';
        nonlinearLineSearchBisectionTolerance = 1e-3;
    }
    if ( nonlinearLineSearchMaxTrials < 1 ) {
        std :: cerr << "nonlinear_line_search_max_trials must be at least 1, setting to 1" << '\n';
        nonlinearLineSearchMaxTrials = 1;
    }
    if ( nonlinearLineSearchArmijo < 0. ) {
        std :: cerr << "nonlinear_line_search_armijo cannot be negative, setting to 0" << '\n';
        nonlinearLineSearchArmijo = 0.;
    }
    if ( nonlinearStagnationIterations < 1 ) {
        std :: cerr << "nonlinear_stagnation_iterations must be at least 1, setting to 1" << '\n';
        nonlinearStagnationIterations = 1;
    }
    if ( nonlinearStagnationRatio <= 0. || nonlinearStagnationRatio >= 1. ) {
        std :: cerr << "nonlinear_stagnation_ratio must be in (0, 1), setting to 0.95" << '\n';
        nonlinearStagnationRatio = 0.95;
    }
    if ( nonlinearRebuildOnSmallAlpha < 0. ) {
        std :: cerr << "nonlinear_rebuild_on_small_alpha cannot be negative, setting to 0" << '\n';
        nonlinearRebuildOnSmallAlpha = 0.;
    }
    if ( nonlinearTrustRadiusMin <= 0. ) {
        std :: cerr << "nonlinear_trust_radius_min must be positive, setting to 1e-12" << '\n';
        nonlinearTrustRadiusMin = 1e-12;
    }
    if ( nonlinearTrustRadiusMax < nonlinearTrustRadiusMin ) {
        std :: cerr << "nonlinear_trust_radius_max smaller than nonlinear_trust_radius_min, setting max to min" << '\n';
        nonlinearTrustRadiusMax = nonlinearTrustRadiusMin;
    }
    if ( nonlinearTrustShrink <= 0. || nonlinearTrustShrink >= 1. ) {
        std :: cerr << "nonlinear_trust_shrink must be in (0, 1), setting to 0.5" << '\n';
        nonlinearTrustShrink = 0.5;
    }
    if ( nonlinearTrustExpand < 1. ) {
        std :: cerr << "nonlinear_trust_expand must be at least 1, setting to 1" << '\n';
        nonlinearTrustExpand = 1.;
    }
    if ( nonlinearTrustMaxTrials < 1 ) {
        std :: cerr << "nonlinear_trust_max_trials must be at least 1, setting to 1" << '\n';
        nonlinearTrustMaxTrials = 1;
    }
    if ( nonlinearTangentCheckEps <= 0. ) {
        std :: cerr << "nonlinear_tangent_check_eps must be positive, setting to 1e-6" << '\n';
        nonlinearTangentCheckEps = 1e-6;
    }
    if ( nonlinearTangentCheckOutput.empty() ) {
        nonlinearTangentCheckOutput = "tangent_check.tsv";
    }
    if ( nonlinearTangentCheckElementOutput.empty() ) {
        nonlinearTangentCheckElementOutput = "tangent_check_elements.tsv";
    }
    if ( nonlinearTangentCheckTopElements < 1 ) {
        std :: cerr << "nonlinear_tangent_check_top_elements must be at least 1, setting to 1" << '\n';
        nonlinearTangentCheckTopElements = 1;
    }
    if ( nonlinearTangentCheckMatrixType.compare("current") != 0 &&
         nonlinearTangentCheckMatrixType.compare("elastic") != 0 &&
         nonlinearTangentCheckMatrixType.compare("secant") != 0 &&
         nonlinearTangentCheckMatrixType.compare("tangent") != 0 &&
         nonlinearTangentCheckMatrixType.compare("consistent") != 0 ) {
        std :: cerr << "nonlinear_tangent_check_matrix_type must be current, elastic, secant, tangent, or consistent; using current" << '\n';
        nonlinearTangentCheckMatrixType = "current";
    }
    if ( !ben ) {
        enlargeIt = maxIt / 3;
    }
    if ( !bsh ) {
        shortenIt = maxIt / 2;
    }

    return this;
};


//////////////////////////////////////////////////////////
void SteadyStateNonLinearSolver :: giveValues(string code, Vector &result) const {
    if ( code.compare("iterations") == 0 ) {
        result.resize(1);
        result [ 0 ] = it;
    } else if ( code.compare("restarts") == 0 ) {
        result.resize(1);
        result [ 0 ] = restarts;
    } else if ( code.compare("error_dofs") == 0 || code.compare("error_displacements") == 0 ) {
        result.resize(1);
        result [ 0 ] = disErr;
    } else if ( code.compare("error_residuals") == 0 ) {
        result.resize(1);
        result [ 0 ] = resErr;
    } else if ( code.compare("error_energy") == 0 ) {
        result.resize(1);
        result [ 0 ] = eneErr;
    } else if ( code.compare("idc_time") == 0 ) {
        result.resize(1);
        result [ 0 ] = idc_time;
    } else if ( code.compare("converged") == 0 ) {
        result.resize(1);
        result [ 0 ] = fully_converged;
    } else {
        SteadyStateLinearSolver :: giveValues(code, result);
    }
}

//////////////////////////////////////////////////////////
void SteadyStateNonLinearSolver :: evaluateErrors() {
    computeTotalInternalAndExternalAndKineticEnergy();

    vector< unsigned >pf  = nodes->givePhysicalFieldsOfDoFs();
    Vector f_extPF = Vector :: Zero(numPhysicalFields);
    Vector f_intPF = Vector :: Zero(numPhysicalFields);
    Vector f_damPF = Vector :: Zero(numPhysicalFields);
    Vector f_accPF = Vector :: Zero(numPhysicalFields);
    Vector residualPF = Vector :: Zero(numPhysicalFields);
    Vector full_ddrPF = Vector :: Zero(numPhysicalFields);
    Vector trial_rPF = Vector :: Zero(numPhysicalFields);
    Vector energyPF = Vector :: Zero(numPhysicalFields);

    unsigned pff;
    for ( unsigned i = 0; i < totalDoFnum; i++ ) {
        pff = pf [ i ];
        residualPF [ pff ] += pow(residuals [ i ], 2);
        f_extPF [ pff ] += pow(f_ext [ i ], 2);
        f_intPF [ pff ] += pow(f_int [ i ], 2);
        f_damPF [ pff ] += pow(f_dam [ i ], 2);
        f_accPF [ pff ] += pow(f_acc [ i ], 2);
        full_ddrPF [ pff ] += pow(full_ddr [ i ], 2);
        trial_rPF [ pff ] += pow(trial_r [ i ], 2);
        energyPF [ pff ] += residuals [ i ] * full_ddr [ i ];
    }

    resErr = disErr = eneErr = 0;
    for ( unsigned i = 0; i < numPhysicalFields; i++ ) {
        trial_rPF [ i ] += eigen_trial_rPF [ i ] + 2 * sqrt(eigen_trial_rPF [ i ] * trial_rPF [ i ]);
        f_extPF [ i ] += eigen_f_extPF [ i ] + 2 * sqrt(eigen_f_extPF [ i ] * f_extPF [ i ]);

        resErr += residualPF [ i ] / max(max( max(f_extPF [ i ], f_intPF [ i ]), max(f_damPF [ i ], f_accPF [ i ]) ), EPS2 [ i ]);
        disErr += full_ddrPF [ i ] / max(trial_rPF [ i ], EPS2 [ i ]);
        eneErr += abs(energyPF [ i ]) / max(max(max(abs(W_ext [ i ]) + eigen_WextPF [ i ], abs(W_int [ i ]) ), abs(W_kin [ i ]) ), EPS2 [ i ]);
        //cout << energyPF [ i ] << " "  << W_ext [ i ] << " "  << W_int [ i ] << " "  << EPS2 << endl;
    }
    resErr = sqrt(resErr);
    disErr = sqrt(disErr);
    eneErr = sqrt(eneErr);
}

//////////////////////////////////////////////////////////
bool SteadyStateNonLinearSolver :: nonlinearGlobalizationActive() const {
    return nonlinearDampingType != NonlinearDampingType :: Off
           || nonlinearLineSearchType != NonlinearLineSearchType :: Off
           || nonlinearTrustRegionType != NonlinearTrustRegionType :: Off
           || nonlinearAdaptiveMatrixUpdate
           || nonlinearStagnationCutback;
}

//////////////////////////////////////////////////////////
bool SteadyStateNonLinearSolver :: nonlinearConvergenceCriteriaSatisfied() const {
    return disErr <= maxDisErr && resErr <= maxResErr && eneErr <= maxEneErr;
}

//////////////////////////////////////////////////////////
double SteadyStateNonLinearSolver :: currentNonlinearMerit() const {
    return currentNonlinearGlobalizationMerit();
}

//////////////////////////////////////////////////////////
double SteadyStateNonLinearSolver :: currentNonlinearGlobalizationMerit() const {
    const double residualMerit = resErr / std :: max(maxResErr, 1e-300);
    const double energyMerit = eneErr / std :: max(maxEneErr, 1e-300);
    double merit = residualMerit;

    // The existing energy error is increment-based: before the first nonlinear
    // correction of a load step, full_ddr is zero and the energy term is
    // therefore not a meaningful line-search reference value. Use residual
    // merit on iteration 0, but keep the normal convergence check below.
    if ( it == 0 ) {
        if ( !std :: isfinite(residualMerit) ) {
            return std :: numeric_limits< double > :: infinity();
        }
        return residualMerit;
    }

    if ( nonlinearMeritType == NonlinearMeritType :: Energy ) {
        merit = energyMerit;
    } else if ( nonlinearMeritType == NonlinearMeritType :: Mixed ) {
        merit = std :: max(residualMerit, energyMerit);
    }

    if ( !std :: isfinite(merit) ) {
        return std :: numeric_limits< double > :: infinity();
    }
    return merit;
}

//////////////////////////////////////////////////////////
double SteadyStateNonLinearSolver :: currentNonlinearConvergenceMerit() const {
    const double residualMerit = resErr / std :: max(maxResErr, 1e-300);
    const double displacementMerit = disErr / std :: max(maxDisErr, 1e-300);
    const double energyMerit = eneErr / std :: max(maxEneErr, 1e-300);
    const double merit = std :: max(residualMerit, std :: max(displacementMerit, energyMerit) );
    if ( !std :: isfinite(merit) ) {
        return std :: numeric_limits< double > :: infinity();
    }
    return merit;
}

//////////////////////////////////////////////////////////
bool SteadyStateNonLinearSolver :: nonlinearMeritAccepted(double trialMerit, double baseMerit, double alpha) const {
    if ( !std :: isfinite(trialMerit) ) {
        return false;
    }
    if ( !std :: isfinite(baseMerit) || baseMerit <= 0. ) {
        return true;
    }
    if ( nonlinearLineSearchArmijo > 0. ) {
        const double armijoFactor = std :: max(0., 1. - nonlinearLineSearchArmijo * alpha);
        return trialMerit <= armijoFactor * baseMerit;
    }
    if ( nonlinearLineSearchAcceptAnyDecrease ) {
        return trialMerit < baseMerit;
    }
    return trialMerit <= baseMerit;
}

//////////////////////////////////////////////////////////
SteadyStateNonLinearSolver :: NonlinearStateSnapshot SteadyStateNonLinearSolver :: saveNonlinearState() const {
    NonlinearStateSnapshot snapshot;
    snapshot.trial_r = trial_r;
    snapshot.ddr = ddr;
    snapshot.full_ddr = full_ddr;
    snapshot.f = f;
    snapshot.f_int = f_int;
    snapshot.f_ext = f_ext;
    snapshot.f_dam = f_dam;
    snapshot.f_acc = f_acc;
    snapshot.residuals = residuals;
    snapshot.load = load;
    snapshot.W_int = W_int;
    snapshot.W_ext = W_ext;
    snapshot.W_kin = W_kin;
    snapshot.disErr = disErr;
    snapshot.resErr = resErr;
    snapshot.eneErr = eneErr;
    if ( nonlinearMaterialSnapshotRollback ) {
        snapshot.materialStatuses = std :: make_shared< ElementContainer :: MaterialStatusSnapshot >( elems->createMaterialStatusSnapshot() );
        if ( nonlinearMaterialSnapshotVerify ) {
            snapshot.materialStatusHash = ElementContainer :: materialStatusSnapshotHash( *snapshot.materialStatuses );
            snapshot.materialStatusHashValid = true;
        }
    }
    return snapshot;
}

//////////////////////////////////////////////////////////
void SteadyStateNonLinearSolver :: restoreNonlinearState(const NonlinearStateSnapshot &snapshot, bool resetMaterialStatuses) {
    trial_r = snapshot.trial_r;
    ddr = snapshot.ddr;
    full_ddr = snapshot.full_ddr;
    f = snapshot.f;
    f_int = snapshot.f_int;
    f_ext = snapshot.f_ext;
    f_dam = snapshot.f_dam;
    f_acc = snapshot.f_acc;
    residuals = snapshot.residuals;
    load = snapshot.load;
    W_int = snapshot.W_int;
    W_ext = snapshot.W_ext;
    W_kin = snapshot.W_kin;
    disErr = snapshot.disErr;
    resErr = snapshot.resErr;
    eneErr = snapshot.eneErr;
    if ( resetMaterialStatuses ) {
        if ( snapshot.materialStatuses ) {
            elems->restoreMaterialStatusSnapshot( *snapshot.materialStatuses );
            if ( nonlinearMaterialSnapshotVerify && snapshot.materialStatusHashValid ) {
                const std :: uint64_t restoredHash = elems->materialStatusStateHash();
                if ( restoredHash != snapshot.materialStatusHash ) {
                    std :: cerr << "Material status snapshot verification failed: expected hash "
                              << snapshot.materialStatusHash << ", restored hash " << restoredHash << std :: endl;
                    exit(1);
                }
            }
        } else {
            elems->resetMaterialStatuses();
        }
    }
}

//////////////////////////////////////////////////////////
void SteadyStateNonLinearSolver :: resetNonlinearGlobalizationAttempt() {
    currentNonlinearDampingFactor = std :: min(nonlinearDampingMax, std :: max(nonlinearDampingMin, nonlinearDampingFactor) );
    nonlinearBestMerit = std :: numeric_limits< double > :: infinity();
    nonlinearStagnationCounter = 0;
    lastNonlinearAlpha = ( nonlinearLineSearchType == NonlinearLineSearchType :: Bisection ) ? 0.5 * nonlinearLineSearchMaxAlpha : 1.;
    lastNonlinearLineSearchTrials = 0;
    lastNonlinearMeritBefore = 0.;
    lastNonlinearMeritAfter = 0.;
    lastNonlinearConvergenceMerit = 0.;
    nonlinearForceMatrixRebuild = false;
    lastNonlinearMatrixRebuild = false;
    currentNonlinearTrustRadius = nonlinearTrustRadiusInitial;
    lastNonlinearTrustRadius = currentNonlinearTrustRadius;
    nonlinearCutbackReason.clear();
    nonlinearRollbackGrowthDecreaseCounter = 0;
    nonlinearTangentCheckDone = false;
}

//////////////////////////////////////////////////////////
bool SteadyStateNonLinearSolver :: maybeRunNonlinearTangentCheck(const Vector &newtonIncrement) {
    if ( !nonlinearTangentCheck || nonlinearTangentCheckDone ) {
        return false;
    }
    if ( nonlinearTangentCheckStep > 0 && static_cast< unsigned >( step ) != nonlinearTangentCheckStep ) {
        return false;
    }
    if ( it != nonlinearTangentCheckIteration ) {
        return false;
    }

    nonlinearTangentCheckDone = true;

    const NonlinearStateSnapshot baseState = saveNonlinearState();
    ElementContainer :: MaterialStatusSnapshot materialSnapshot = elems->createMaterialStatusSnapshot();
    const Vector baseResidual = f;
    std :: vector< std :: pair< std :: string, Vector > > directions;

    if ( nonlinearTangentCheckIncludeNewton && newtonIncrement.size() == static_cast< int >( freeDoFnum ) ) {
        const double nrm = newtonIncrement.norm();
        if ( nrm > 0. && std :: isfinite(nrm) ) {
            directions.push_back( { "newton", newtonIncrement / nrm } );
        }
    }

    for ( unsigned r = 0; r < nonlinearTangentCheckRandomVectors; r++ ) {
        Vector p = Vector :: Zero(freeDoFnum);
        unsigned long long state = 1469598103934665603ULL + static_cast< unsigned long long >( r + 1 ) * 1099511628211ULL;
        for ( unsigned i = 0; i < freeDoFnum; i++ ) {
            state = state * 2862933555777941757ULL + 3037000493ULL;
            const double unit = static_cast< double >( ( state >> 11 ) & 0x1fffffULL ) / static_cast< double >( 0x1fffffULL );
            p [ i ] = 2. * unit - 1.;
        }
        const double nrm = p.norm();
        if ( nrm > 0. && std :: isfinite(nrm) ) {
            directions.push_back( { "random_" + std :: to_string(r + 1), p / nrm } );
        }
    }

    std :: ifstream input(nonlinearTangentCheckOutput);
    const bool writeHeader = !input.good() || input.peek() == std :: ifstream :: traits_type :: eof();
    input.close();
    std :: ofstream output(nonlinearTangentCheckOutput, std :: ios :: app);
    if ( writeHeader ) {
        output << "step\titeration\teps\tdirection\trelative_error\tcosine\tkp_norm\tfd_norm\tresidual_norm\n";
    }

    for ( const auto &direction : directions ) {
        restoreNonlinearState(baseState, false);
        elems->restoreMaterialStatusSnapshot(materialSnapshot);
        ddr = nonlinearTangentCheckEps * direction.second;
        updateFieldVariables();
        computeForcesAtIntegrationTime(false);
        Vector perturbedResidual = Vector :: Zero(freeDoFnum);
        nodes->giveReducedForceArray(residuals, perturbedResidual);

        const Vector kp = Keff * direction.second;
        const Vector fd = ( baseResidual - perturbedResidual ) / nonlinearTangentCheckEps;
        const double kpNorm = kp.norm();
        const double fdNorm = fd.norm();
        const double relativeError = ( kp - fd ).norm() / std :: max(fdNorm, 1e-300);
        const double cosine = kp.dot(fd) / std :: max(kpNorm * fdNorm, 1e-300);

        output << step << '\t'
               << it << '\t'
               << nonlinearTangentCheckEps << '\t'
               << direction.first << '\t'
               << relativeError << '\t'
               << cosine << '\t'
               << kpNorm << '\t'
               << fdNorm << '\t'
               << baseResidual.norm() << '\n';

        if ( !silent ) {
            std :: cout << "TANGENT_CHECK"
                        << " step " << step
                        << " it " << it
                        << " eps " << nonlinearTangentCheckEps
                        << " direction " << direction.first
                        << " relative_error " << relativeError
                        << " cosine " << cosine
                        << " kp_norm " << kpNorm
                        << " fd_norm " << fdNorm
                        << std :: endl;
        }

        if ( nonlinearTangentCheckScope == NonlinearTangentCheckScope :: ElementTop ) {
            writeElementTangentAttribution(direction.first, direction.second, baseState, materialSnapshot);
        }
    }
    output.close();

    restoreNonlinearState(baseState, false);
    elems->restoreMaterialStatusSnapshot(materialSnapshot);
    return true;
}

//////////////////////////////////////////////////////////
void SteadyStateNonLinearSolver :: writeElementTangentAttribution(const std :: string &directionName, const Vector &direction, const NonlinearStateSnapshot &baseState, const ElementContainer :: MaterialStatusSnapshot &materialSnapshot) {
    struct ElementTangentAttribution {
        unsigned elementId = 0;
        std :: string elementName;
        double relativeError = 0.;
        double cosine = 0.;
        double kpNorm = 0.;
        double fdNorm = 0.;
        double mismatchNorm = 0.;
    };

    std :: string matrixType = nonlinearTangentCheckMatrixType;
    if ( matrixType.compare("current") == 0 ) {
        matrixType = lastStiffMatType.empty() ? stiffMatType : lastStiffMatType;
    }

    std :: string matrixTypeLabel = matrixType;
    const bool useElasticBlend = stiffMatElasticBlendBeta > 0. && matrixType.compare("elastic") != 0;
    if ( useElasticBlend ) {
        matrixTypeLabel += "+elastic_blend_" + std :: to_string(stiffMatElasticBlendBeta);
    }

    Vector fullDirection = Vector :: Zero( trial_r.size() );
    nodes->giveFullDoFArray(direction, fullDirection);

    restoreNonlinearState(baseState, false);
    elems->restoreMaterialStatusSnapshot(materialSnapshot);

    std :: vector< ElementTangentAttribution > rows;
    rows.reserve(elems->giveSize() );

    for ( auto elemIt = elems->begin(); elemIt != elems->end(); ++elemIt ) {
        Element *elem = *elemIt;
        Element :: MaterialStatusSnapshot elementSnapshot = elem->createMaterialStatusSnapshot();

        std :: vector< unsigned >elDoFs = elem->giveDoFs();
        Vector baseElementDoFs = Vector :: Zero(elDoFs.size() );
        Vector elementDirection = Vector :: Zero(elDoFs.size() );
        for ( unsigned i = 0; i < elDoFs.size(); i++ ) {
            baseElementDoFs [ i ] = baseState.trial_r [ elDoFs [ i ] ];
            elementDirection [ i ] = fullDirection [ elDoFs [ i ] ];
        }

        Matrix elementK = elem->giveStiffnessMatrix(matrixType);
        if ( useElasticBlend ) {
            Matrix elasticK = elem->giveStiffnessMatrix("elastic");
            elementK = ( 1. - stiffMatElasticBlendBeta ) * elementK + stiffMatElasticBlendBeta * elasticK;
        }
        Vector kp = elementK * elementDirection;

        elem->restoreMaterialStatusSnapshot(elementSnapshot);
        elem->evaluateStrains(baseElementDoFs);
        elem->evaluateStresses(false, -1.);
        Vector baseForces = elem->giveInternalForces();

        elem->restoreMaterialStatusSnapshot(elementSnapshot);
        Vector perturbedElementDoFs = baseElementDoFs + nonlinearTangentCheckEps * elementDirection;
        elem->evaluateStrains(perturbedElementDoFs);
        elem->evaluateStresses(false, -1.);
        Vector perturbedForces = elem->giveInternalForces();
        elem->restoreMaterialStatusSnapshot(elementSnapshot);

        Vector fd = ( perturbedForces - baseForces ) / nonlinearTangentCheckEps;
        const double kpNorm = kp.norm();
        const double fdNorm = fd.norm();
        const double mismatchNorm = ( kp - fd ).norm();

        ElementTangentAttribution row;
        row.elementId = elem->giveID();
        row.elementName = elem->giveName();
        row.relativeError = mismatchNorm / std :: max(fdNorm, 1e-300);
        row.cosine = kp.dot(fd) / std :: max(kpNorm * fdNorm, 1e-300);
        row.kpNorm = kpNorm;
        row.fdNorm = fdNorm;
        row.mismatchNorm = mismatchNorm;
        rows.push_back(row);
    }

    std :: sort(rows.begin(), rows.end(), [](const ElementTangentAttribution &a, const ElementTangentAttribution &b) {
        return a.mismatchNorm > b.mismatchNorm;
    } );

    std :: ifstream input(nonlinearTangentCheckElementOutput);
    const bool writeHeader = !input.good() || input.peek() == std :: ifstream :: traits_type :: eof();
    input.close();
    std :: ofstream output(nonlinearTangentCheckElementOutput, std :: ios :: app);
    if ( writeHeader ) {
        output << "step\titeration\teps\tdirection\telement_id\telement_name\tmatrix_type\trelative_error\tcosine\tkp_norm\tfd_norm\tmismatch_norm\n";
    }

    const unsigned numRows = std :: min<unsigned>(nonlinearTangentCheckTopElements, rows.size() );
    for ( unsigned i = 0; i < numRows; i++ ) {
        const ElementTangentAttribution &row = rows [ i ];
        output << step << '\t'
               << it << '\t'
               << nonlinearTangentCheckEps << '\t'
               << directionName << '\t'
               << row.elementId << '\t'
               << row.elementName << '\t'
               << matrixTypeLabel << '\t'
               << row.relativeError << '\t'
               << row.cosine << '\t'
               << row.kpNorm << '\t'
               << row.fdNorm << '\t'
               << row.mismatchNorm << '\n';
    }
    output.close();

    restoreNonlinearState(baseState, false);
    elems->restoreMaterialStatusSnapshot(materialSnapshot);
}

//////////////////////////////////////////////////////////
bool SteadyStateNonLinearSolver :: applyScaledIncrementAndEvaluate(const NonlinearStateSnapshot &baseState, const Vector &increment, double alpha, bool frozen, bool resetMaterialStatuses) {
    restoreNonlinearState(baseState, resetMaterialStatuses);
    ddr = alpha * increment;
    updateFieldVariables();
    computeForcesAtIntegrationTime(frozen);
    evaluateErrors();
    if ( it == 0 ) {
        disErr = 0;                        //error in displacement change, only from second iteration
    }
    return std :: isfinite(resErr) && std :: isfinite(disErr) && std :: isfinite(eneErr);
}

//////////////////////////////////////////////////////////
SteadyStateNonLinearSolver :: NonlinearTrialResult SteadyStateNonLinearSolver :: performBacktrackingLineSearch(const NonlinearStateSnapshot &baseState, const Vector &increment, double meritBefore) {
    NonlinearTrialResult result;
    result.meritBefore = meritBefore;
    double alpha = nonlinearLineSearchMaxAlpha;

    for ( unsigned trial = 0; trial < nonlinearLineSearchMaxTrials && alpha >= nonlinearLineSearchMinAlpha * ( 1. - 1e-12 ); trial++ ) {
        const bool useFrozenTrial = nonlinearLineSearchEvaluationMode != NonlinearLineSearchEvaluationMode :: Actual;
        const bool requireActualAcceptance = nonlinearLineSearchEvaluationMode == NonlinearLineSearchEvaluationMode :: FrozenThenActual;
        const bool finiteFirstTrial = applyScaledIncrementAndEvaluate(baseState, increment, alpha, useFrozenTrial);
        double trialMerit = currentNonlinearGlobalizationMerit();
        const bool acceptedFirstTrial = finiteFirstTrial && nonlinearMeritAccepted(trialMerit, meritBefore, alpha);

        if ( !silent && nonlinearLineSearchReportTrials ) {
            std :: cout << "LS_TRIAL"
                        << " it " << it
                        << " trial " << trial + 1
                        << " alpha " << alpha
                        << " mode " << ( useFrozenTrial ? "frozen" : "actual" )
                        << " finite " << finiteFirstTrial
                        << " accepted " << acceptedFirstTrial
                        << " merit_before " << meritBefore
                        << " merit_trial " << trialMerit
                        << " res " << resErr
                        << " disp " << disErr
                        << " energy " << eneErr
                        << std :: endl;
        }

        if ( acceptedFirstTrial ) {
            if ( useFrozenTrial ) {
                const bool finiteActualTrial = applyScaledIncrementAndEvaluate(baseState, increment, alpha, false);
                trialMerit = currentNonlinearGlobalizationMerit();
                const bool acceptedActualTrial = finiteActualTrial && ( !requireActualAcceptance || nonlinearMeritAccepted(trialMerit, meritBefore, alpha) );
                if ( !silent && nonlinearLineSearchReportTrials ) {
                    std :: cout << "LS_ACTUAL"
                                << " it " << it
                                << " alpha " << alpha
                                << " finite " << finiteActualTrial
                                << " accepted " << acceptedActualTrial
                                << " merit_actual " << trialMerit
                                << " res " << resErr
                                << " disp " << disErr
                                << " energy " << eneErr
                                << std :: endl;
                }
                if ( !acceptedActualTrial ) {
                    restoreNonlinearState(baseState, true);
                    result.alpha = alpha;
                    result.trials = trial + 1;
                    result.meritAfter = trialMerit;
                    alpha *= nonlinearLineSearchReduction;
                    continue;
                }
            }
            result.accepted = true;
            result.alpha = alpha;
            result.trials = trial + 1;
            result.meritAfter = trialMerit;
            return result;
        }

        restoreNonlinearState(baseState, true);
        alpha *= nonlinearLineSearchReduction;
    }

    restoreNonlinearState(baseState, true);
    if ( !silent && nonlinearLineSearchReportTrials ) {
        std :: cout << "LS_FAIL"
                    << " it " << it
                    << " trials " << nonlinearLineSearchMaxTrials
                    << " merit_before " << meritBefore
                    << " cutback_on_fail " << nonlinearLineSearchCutbackOnFail
                    << std :: endl;
    }
    result.alpha = alpha;
    result.trials = nonlinearLineSearchMaxTrials;
    result.meritAfter = currentNonlinearMerit();
    return result;
}

//////////////////////////////////////////////////////////
SteadyStateNonLinearSolver :: NonlinearTrialResult SteadyStateNonLinearSolver :: performBisectionLineSearch(const NonlinearStateSnapshot &baseState, const Vector &increment, double meritBefore) {
    NonlinearTrialResult result;
    result.meritBefore = meritBefore;

    double lower = 0.;
    double upper = std :: max(nonlinearLineSearchBisectionTolerance, 2. * std :: max(lastNonlinearAlpha, nonlinearLineSearchBisectionTolerance) );
    upper = std :: min(upper, nonlinearLineSearchMaxAlpha);
    double acceptedAlpha = 0.;
    double acceptedMerit = std :: numeric_limits< double > :: infinity();
    const unsigned maxBisectionTrials = 64;
    unsigned trials = 0;

    while ( upper - lower > nonlinearLineSearchBisectionTolerance && trials < maxBisectionTrials ) {
        const double alpha = 0.5 * ( lower + upper );
        const bool useFrozenTrial = nonlinearLineSearchEvaluationMode != NonlinearLineSearchEvaluationMode :: Actual;
        const bool finiteTrial = applyScaledIncrementAndEvaluate(baseState, increment, alpha, useFrozenTrial);
        const double trialMerit = currentNonlinearGlobalizationMerit();
        trials++;

        if ( finiteTrial && nonlinearMeritAccepted(trialMerit, meritBefore, alpha) ) {
            lower = alpha;
            acceptedAlpha = alpha;
            acceptedMerit = trialMerit;
        } else {
            upper = alpha;
        }
        restoreNonlinearState(baseState, true);
    }

    result.trials = trials;
    if ( acceptedAlpha <= 0. ) {
        result.alpha = 0.;
        result.meritAfter = currentNonlinearMerit();
        return result;
    }

    const bool finiteActualTrial = applyScaledIncrementAndEvaluate(baseState, increment, acceptedAlpha, false);
    if ( !finiteActualTrial || ( nonlinearLineSearchEvaluationMode == NonlinearLineSearchEvaluationMode :: FrozenThenActual && !nonlinearMeritAccepted(currentNonlinearGlobalizationMerit(), meritBefore, acceptedAlpha) ) ) {
        restoreNonlinearState(baseState, true);
        result.alpha = acceptedAlpha;
        result.meritAfter = acceptedMerit;
        return result;
    }

    result.accepted = true;
    result.alpha = acceptedAlpha;
    result.meritAfter = currentNonlinearGlobalizationMerit();
    return result;
}

//////////////////////////////////////////////////////////
SteadyStateNonLinearSolver :: NonlinearTrialResult SteadyStateNonLinearSolver :: performStepNormTrustRegion(const NonlinearStateSnapshot &baseState, const Vector &increment, double meritBefore) {
    NonlinearTrialResult result;
    result.meritBefore = meritBefore;

    const double incrementNorm = increment.norm();
    if ( !std :: isfinite(incrementNorm) || incrementNorm <= 0. ) {
        result.alpha = 1.;
        result.trials = 1;
        const bool finiteTrial = applyScaledIncrementAndEvaluate(baseState, increment, 1., false);
        result.meritAfter = currentNonlinearGlobalizationMerit();
        result.accepted = finiteTrial;
        return result;
    }

    if ( currentNonlinearTrustRadius <= 0. || !std :: isfinite(currentNonlinearTrustRadius) ) {
        currentNonlinearTrustRadius = nonlinearTrustRadiusInitial > 0. ? nonlinearTrustRadiusInitial : incrementNorm;
    }
    currentNonlinearTrustRadius = std :: min(nonlinearTrustRadiusMax, std :: max(nonlinearTrustRadiusMin, currentNonlinearTrustRadius) );

    for ( unsigned trial = 0; trial < nonlinearTrustMaxTrials; trial++ ) {
        const double alpha = std :: min(1., currentNonlinearTrustRadius / incrementNorm);
        const bool finiteTrial = applyScaledIncrementAndEvaluate(baseState, increment, alpha, false);
        const double trialMerit = currentNonlinearGlobalizationMerit();
        result.trials = trial + 1;

        if ( finiteTrial && ( nonlinearConvergenceCriteriaSatisfied() || nonlinearMeritAccepted(trialMerit, meritBefore, alpha) ) ) {
            result.accepted = true;
            result.alpha = alpha;
            result.meritAfter = trialMerit;
            lastNonlinearTrustRadius = currentNonlinearTrustRadius;
            if ( std :: isfinite(meritBefore) && meritBefore > 0. && trialMerit < meritBefore ) {
                currentNonlinearTrustRadius = std :: min(nonlinearTrustRadiusMax, currentNonlinearTrustRadius * nonlinearTrustExpand);
            }
            return result;
        }

        restoreNonlinearState(baseState, true);
        currentNonlinearTrustRadius = std :: max(nonlinearTrustRadiusMin, currentNonlinearTrustRadius * nonlinearTrustShrink);
    }

    restoreNonlinearState(baseState, true);
    result.alpha = std :: min(1., currentNonlinearTrustRadius / incrementNorm);
    result.meritAfter = currentNonlinearGlobalizationMerit();
    lastNonlinearTrustRadius = currentNonlinearTrustRadius;
    return result;
}

//////////////////////////////////////////////////////////
double SteadyStateNonLinearSolver :: currentNonlinearDampingAlpha() const {
    if ( nonlinearDampingType == NonlinearDampingType :: Off ) {
        return 1.;
    }
    if ( nonlinearDampingType == NonlinearDampingType :: Adaptive || nonlinearDampingType == NonlinearDampingType :: RollbackAdaptive ) {
        return std :: min(nonlinearDampingMax, std :: max(nonlinearDampingMin, currentNonlinearDampingFactor) );
    }
    return std :: min(nonlinearDampingMax, std :: max(nonlinearDampingMin, nonlinearDampingFactor) );
}

//////////////////////////////////////////////////////////
void SteadyStateNonLinearSolver :: updateAdaptiveNonlinearDamping(double meritBefore, double meritAfter) {
    if ( nonlinearDampingType != NonlinearDampingType :: Adaptive ) {
        return;
    }
    if ( std :: isfinite(meritAfter) && ( !std :: isfinite(meritBefore) || meritAfter < meritBefore ) ) {
        currentNonlinearDampingFactor = std :: min(nonlinearDampingMax, currentNonlinearDampingFactor * nonlinearDampingIncrease);
    } else {
        currentNonlinearDampingFactor = std :: max(nonlinearDampingMin, currentNonlinearDampingFactor * nonlinearDampingDecrease);
    }
}

//////////////////////////////////////////////////////////
bool SteadyStateNonLinearSolver :: shouldCutbackForNonlinearProgress(double meritBefore, double meritAfter) {
    if ( !nonlinearStagnationCutback || !std :: isfinite(meritAfter) ) {
        return false;
    }
    if ( nonlinearGrowthCutback > 0. && std :: isfinite(meritBefore) && meritBefore > 0. && meritAfter > nonlinearGrowthCutback * meritBefore ) {
        nonlinearCutbackReason = "nonlinear_growth";
        if ( nonlinearAdaptiveMatrixUpdate && nonlinearRebuildOnMeritGrowth ) {
            nonlinearForceMatrixRebuild = true;
        }
        return true;
    }
    if ( !std :: isfinite(nonlinearBestMerit) ) {
        nonlinearBestMerit = meritAfter;
        nonlinearStagnationCounter = 0;
        return false;
    }
    if ( meritAfter <= nonlinearBestMerit * nonlinearStagnationRatio ) {
        nonlinearBestMerit = meritAfter;
        nonlinearStagnationCounter = 0;
        return false;
    }
    nonlinearStagnationCounter++;
    if ( nonlinearStagnationCounter >= nonlinearStagnationIterations ) {
        nonlinearCutbackReason = "nonlinear_stagnation";
        if ( nonlinearAdaptiveMatrixUpdate && nonlinearRebuildOnStagnation ) {
            nonlinearForceMatrixRebuild = true;
        }
        return true;
    }
    return false;
}

//////////////////////////////////////////////////////////
void SteadyStateNonLinearSolver :: requestAdaptiveMatrixRebuildIfNeeded(double meritBefore, double meritAfter) {
    if ( !nonlinearAdaptiveMatrixUpdate ) {
        return;
    }
    if ( nonlinearRebuildOnSmallAlpha > 0. && lastNonlinearAlpha > 0. && lastNonlinearAlpha < nonlinearRebuildOnSmallAlpha ) {
        nonlinearForceMatrixRebuild = true;
        if ( nonlinearCutbackReason.empty() ) {
            nonlinearCutbackReason = "matrix_rebuild_small_alpha";
        }
    }
    if ( nonlinearRebuildOnMeritGrowth && std :: isfinite(meritBefore) && meritBefore > 0. && std :: isfinite(meritAfter) && meritAfter > meritBefore ) {
        nonlinearForceMatrixRebuild = true;
        if ( nonlinearCutbackReason.empty() ) {
            nonlinearCutbackReason = "matrix_rebuild_merit_growth";
        }
    }
}

//////////////////////////////////////////////////////////
void SteadyStateNonLinearSolver :: setEigenErrorValue_rPF(unsigned pf, double value) {
    if ( pf >= numPhysicalFields ) {
        cout << "SteadyStateNonLinearSolver Error: value " << pf << " larger than number of physical fields " << numPhysicalFields << endl;
        exit(1);
    }
    eigen_trial_rPF [ pf ] = value;
}

//////////////////////////////////////////////////////////
void SteadyStateNonLinearSolver :: setEigenErrorValue_fext(unsigned pf, double value) {
    if ( pf >= numPhysicalFields ) {
        cout << "SteadyStateNonLinearSolver Error: value " << pf << " larger than number of physical fields " << numPhysicalFields << endl;
        exit(1);
    }
    eigen_f_extPF [ pf ] = value;
}

//////////////////////////////////////////////////////////
void SteadyStateNonLinearSolver :: setEigenErrorValue_Wext(unsigned pf, double value) {
    if ( pf >= numPhysicalFields ) {
        cout << "SteadyStateNonLinearSolver Error: value " << pf << " larger than number of physical fields " << numPhysicalFields << endl;
        exit(1);
    }
    eigen_WextPF [ pf ] = value;
}

//////////////////////////////////////////////////////////
void SteadyStateNonLinearSolver :: reset() {
    load.setZero();
    ddr.setZero();
    double reset_time = time;
    if ( idc ) {
        reset_time = idc_time;
    }

    bool converged = false;

    restarts = 0;
    while ( !converged ) {
        nodes->addRHS_nodalLoad(load, reset_time); //add nodal load
        nodes->updateDirrichletBC(trial_r, reset_time); //give prescribed DoFs
        updateFieldVariables();     //with ddr=0
        computeForcesAtIntegrationTime(true);

        it = 0;
        while ( !converged && it < maxIt ) {
            if ( updateSystemMatrices(it, cumul_it, false) ) {
                computeKeff();                                    //only if required
            }
            nodes->giveReducedForceArray(residuals, f);   // NOTE JK when IDC applied and step reset, residuals from the last iteration are used here //JE: no, they are actually computed again here

            /*
             * if ( LinalgSymmetricSolver(Keff, ddr, f, ddr, conj_grad_precision, conj_grad_relative_maxit, symsolver_type) == false ) {
             *  std :: cerr << "Conjugate gradients did not converge, attempt to restart step" << endl;
             *  it = maxIt;
             *  resErr = 1e10;
             *  break;
             * }
             */

            const bool directSolveSuccess = linalgsolver->solve(ddr, f);
            collectLinearDeflationVector(ddr, directSolveSuccess);

            //update DoFs
            updateFieldVariables();
            //compute residuals
            computeForcesAtIntegrationTime(true); //to obtain the actual stress, fluxes, ...

            //compute and print errors
            evaluateErrors();
            if ( it == 0 ) {
                disErr = 0;                        //error in displacement change, only from second iteration
            }
            if ( not silent ) {
                cout << setw(6) << it << setw(15) << resErr;
            }
            if ( it == 0 ) {
                if ( not silent ) {
                    cout << setw(15) << "---";
                }
            } else {
                if ( not silent ) {
                    cout << setw(15) << disErr;
                }
            }
            if ( not silent ) {
                cout << setw(15) << eneErr << endl;
            }

            if ( std :: isnan(resErr) || std :: isnan(disErr) || std :: isnan(eneErr) ) {
                std :: cerr << "calculating with NaN in ";
                if ( std :: isnan(resErr) ) {
                    std :: cerr << "\tresiduals ";
                }
                if ( std :: isnan(disErr) ) {
                    std :: cerr << "\tdisplacements ";
                }
                if ( std :: isnan(eneErr) ) {
                    std :: cerr << "\tenergies ";
                }
                std :: cerr << endl;
                it = maxIt;
                resErr = 1e10;
                break;
            }

            if ( disErr <= maxDisErr && resErr <= maxResErr && eneErr <= maxEneErr ) {
                converged = true;
            } else {
                converged = false;
            }
            it++;
            cumul_it++;
        }

        if ( converged ) {
            this->fully_converged = true;
        } else if ( !converged ) {
            if ( disErr < limitDisErr && resErr < limitResErr && eneErr < limitEneErr ) {
                if ( not silent ) {
                    std :: cout << "tolerance increased in this step" << '\n';
                }
                converged = true;
                this->fully_converged = false;
                computeForcesAtStepEnd(false); //to obtain the actual stress, fluxes, ...
            } else {
                std :: cerr << "Error: " << name << " did not converge to the solution" << endl;
                terminated = true;
                return;
            }
        }
    }
    runAfterEachStep();
}

//////////////////////////////////////////////////////////
void SteadyStateNonLinearSolver :: solve() {
    double load_mult;
    bool converged = false;
    bool restarted = false;
    bool restart_now = false;
    Vector help_idc_r, help_idc_f;

    //MyVector reset_residuals = residuals;   ///> if step restarted when IDC applied, residuals need to be reset to stage before the step start
    //JE: no, these are recomputed
    restarts = 0;
    while ( !converged ) {
        resetNonlinearGlobalizationAttempt();
        //setup loading

        if ( !idc ) {
            nodes->addRHS_nodalLoad(load, time); //add nodal load
            nodes->updateDirrichletBC(trial_r, time); //give prescribed DoFs
            updateFieldVariables();      //with ddr=0
            computeForcesAtIntegrationTime(true);
        } else {
            help_idc_r =  Vector :: Zero(totalDoFnum);
            help_idc_f =  Vector :: Zero(totalDoFnum);
        }

        it = 0;
        while ( !converged && it < maxIt ) {
            lastNonlinearMatrixRebuild = false;
            const bool forceMatrixRebuild = nonlinearForceMatrixRebuild;
            if ( ( step > 0 || it > 0 ) && updateSystemMatrices(it, cumul_it, forceMatrixRebuild) ) {
                computeKeff();                                    //only if required
                lastNonlinearMatrixRebuild = true;
            }
            nonlinearForceMatrixRebuild = false;
            nodes->giveReducedForceArray(residuals, f);

            const bool useNonlinearGlobalization = ( !idc ) && nonlinearGlobalizationActive();
            NonlinearStateSnapshot incrementBaseState;
            double meritBefore = lastNonlinearMeritBefore;
            if ( useNonlinearGlobalization ) {
                evaluateErrors();
                if ( it == 0 ) {
                    disErr = 0;
                }
                meritBefore = currentNonlinearGlobalizationMerit();
                incrementBaseState = saveNonlinearState();
                lastNonlinearMeritBefore = meritBefore;
                lastNonlinearConvergenceMerit = currentNonlinearConvergenceMerit();
                nonlinearCutbackReason.clear();
            } else if ( idc && nonlinearGlobalizationActive() && !warnedGlobalizationWithIDC ) {
                std :: cerr << "nonlinear globalization controls are not applied with indirect control; using legacy indirect-control update" << endl;
                warnedGlobalizationWithIDC = true;
            }

            if ( idc ) {      //indirect displacement control
                Vector trial_r_last_iter = trial_r;
                f_last_iter = f;

                load.setZero();
                nodes->addRHS_nodalLoad(load, idc_time + idc_dt); //add nodal load
                nodes->updateDirrichletBC(trial_r, idc_time + idc_dt); //give prescribed DoFs
                computeForcesAtIntegrationTime(true);
                nodes->giveReducedForceArray(residuals, f);

                linalgsolver->solve(ddr, f_last_iter); //A: original, no additiona load, time idc_time
                linalgsolver->solve(ddf, f - f_last_iter); //B: added load, time idc_time + idc_dt

                if ( !restart_now ) {
                    nodes->giveFullDoFArray(ddf, full_ddf);
                    nodes->giveFullDoFArray(ddr, full_ddr);

                    //compute B
                    Vector trial_r_B = trial_r_last_iter + full_ddf + full_ddr;
                    nodes->updateDirrichletBC(trial_r_B, idc_time + idc_dt); //give prescribed DoFs
                    Vector fext_B;
                    if ( idc->requireForces() ) {
                        load.setZero();
                        nodes->addRHS_nodalLoad(load, idc_time + idc_dt); //add nodal load
                        computeInternalExternalForces(trial_r_B, load, false, idc_time + idc_dt, idc_dt);
                        fext_B = f_ext;
                    }

                    //compute A
                    Vector trial_r_A = trial_r_last_iter + full_ddr;
                    Vector fext_A;
                    if ( idc->requireForces() ) {
                        load.setZero();
                        nodes->addRHS_nodalLoad(load, idc_time); //add nodal load
                        computeInternalExternalForces(trial_r_A, load, false, idc_time, idc_dt);
                        fext_A = f_ext;
                    }

                    //differences and load multiplier
                    trial_r_B -= trial_r_A;
                    fext_B -= fext_A;
                    load_mult = idc->giveMultiplierCorrection(trial_r_A, fext_A, trial_r_B, fext_B, time);

                    ddr += load_mult * ddf;
                    idc_time += idc_dt * load_mult;
                    trial_r = trial_r_last_iter;

                    load.setZero();
                    nodes->addRHS_nodalLoad(load, idc_time); //add nodal load
                    nodes->updateDirrichletBC(trial_r, idc_time); //give prescribed DoFs
                }
            } else {         //direct controll
                const bool directSolveSuccess = linalgsolver->solve(ddr, f);
                Vector ddrNewton = ddr;
                collectLinearDeflationVector(ddrNewton, directSolveSuccess);
                if ( maybeRunNonlinearTangentCheck(ddrNewton) && nonlinearTangentCheckStopAfter ) {
                    terminated = true;
                    return;
                }

                if ( useNonlinearGlobalization ) {
                    if ( nonlinearLineSearchType == NonlinearLineSearchType :: Backtracking || nonlinearLineSearchType == NonlinearLineSearchType :: Bisection ) {
                        NonlinearTrialResult lineSearchResult;
                        if ( nonlinearLineSearchType == NonlinearLineSearchType :: Bisection ) {
                            lineSearchResult = performBisectionLineSearch(incrementBaseState, ddrNewton, meritBefore);
                        } else {
                            lineSearchResult = performBacktrackingLineSearch(incrementBaseState, ddrNewton, meritBefore);
                        }
                        lastNonlinearAlpha = lineSearchResult.alpha;
                        lastNonlinearLineSearchTrials = lineSearchResult.trials;
                        lastNonlinearMeritAfter = lineSearchResult.meritAfter;
                        if ( !lineSearchResult.accepted ) {
                            nonlinearCutbackReason = "line_search_failed";
                            if ( nonlinearLineSearchCutbackOnFail ) {
                                restart_now = true;
                                break;
                            }
                            applyScaledIncrementAndEvaluate(incrementBaseState, ddrNewton, 1., false);
                            lastNonlinearAlpha = 1.;
                            lastNonlinearMeritAfter = currentNonlinearGlobalizationMerit();
                            nonlinearCutbackReason.clear();
                        }
                    } else if ( nonlinearTrustRegionType == NonlinearTrustRegionType :: StepNorm ) {
                        const NonlinearTrialResult trustRegionResult = performStepNormTrustRegion(incrementBaseState, ddrNewton, meritBefore);
                        lastNonlinearAlpha = trustRegionResult.alpha;
                        lastNonlinearLineSearchTrials = trustRegionResult.trials;
                        lastNonlinearMeritAfter = trustRegionResult.meritAfter;
                        if ( !trustRegionResult.accepted ) {
                            nonlinearCutbackReason = "trust_region_failed";
                            restart_now = true;
                            break;
                        }
                    } else {
                        const double alpha = currentNonlinearDampingAlpha();
                        const double previousAcceptedMerit = lastNonlinearMeritAfter;
                        applyScaledIncrementAndEvaluate(incrementBaseState, ddrNewton, alpha, false, false);
                        lastNonlinearAlpha = alpha;
                        lastNonlinearLineSearchTrials = 0;
                        lastNonlinearMeritAfter = currentNonlinearGlobalizationMerit();
                        if ( nonlinearDampingType == NonlinearDampingType :: RollbackAdaptive && it > nonlinearRollbackIgnoreInitialIterations && std :: isfinite(previousAcceptedMerit) && previousAcceptedMerit > 0. ) {
                            if ( lastNonlinearMeritAfter > nonlinearRollbackGrowthRatio * previousAcceptedMerit ) {
                                const double newAlpha = std :: max(nonlinearDampingMin, alpha * nonlinearDampingDecrease);
                                const bool canDecrease = nonlinearRollbackGrowthDecreaseCounter < nonlinearRollbackMaxGrowthDecreases
                                                         && newAlpha < alpha * ( 1. - 1e-12 );
                                if ( nonlinearRollbackEnable && canDecrease ) {
                                    if ( !silent ) {
                                        std :: cout << "NONLINEAR_DAMPING_ROLLBACK"
                                                    << " it " << it
                                                    << " alpha " << alpha
                                                    << " new_alpha " << newAlpha
                                                    << " previous_merit " << previousAcceptedMerit
                                                    << " trial_merit " << lastNonlinearMeritAfter
                                                    << std :: endl;
                                    }
                                    restoreNonlinearState(incrementBaseState, true);
                                    currentNonlinearDampingFactor = newAlpha;
                                    nonlinearRollbackGrowthDecreaseCounter++;
                                    lastNonlinearMeritAfter = previousAcceptedMerit;
                                    nonlinearCutbackReason = "damping_rollback";
                                    continue;
                                }
                                if ( canDecrease ) {
                                    currentNonlinearDampingFactor = newAlpha;
                                    nonlinearRollbackGrowthDecreaseCounter++;
                                    nonlinearCutbackReason = "damping_growth_decrease";
                                    if ( !silent ) {
                                        std :: cout << "NONLINEAR_DAMPING_GROWTH_DECREASE"
                                                    << " it " << it
                                                    << " alpha " << alpha
                                                    << " new_alpha " << newAlpha
                                                    << " previous_merit " << previousAcceptedMerit
                                                    << " trial_merit " << lastNonlinearMeritAfter
                                                    << std :: endl;
                                    }
                                }
                                if ( !canDecrease ) {
                                    if ( !silent ) {
                                        std :: cout << "NONLINEAR_DAMPING_GROWTH_LIMIT"
                                                    << " it " << it
                                                    << " alpha " << alpha
                                                    << " previous_merit " << previousAcceptedMerit
                                                    << " trial_merit " << lastNonlinearMeritAfter
                                                    << std :: endl;
                                    }
                                    currentNonlinearDampingFactor = alpha;
                                    nonlinearCutbackReason = "damping_growth_limit";
                                }
                            } else if ( nonlinearRollbackDecreaseOnAnyIncrease && lastNonlinearMeritAfter > nonlinearRollbackIncreaseRatio * previousAcceptedMerit ) {
                                if ( nonlinearRollbackGrowthDecreaseCounter < nonlinearRollbackMaxGrowthDecreases ) {
                                    const double newAlpha = std :: max(nonlinearDampingMin, alpha * nonlinearDampingDecrease);
                                    currentNonlinearDampingFactor = newAlpha;
                                    nonlinearRollbackGrowthDecreaseCounter++;
                                    nonlinearCutbackReason = "damping_increase_decrease";
                                    if ( !silent ) {
                                        std :: cout << "NONLINEAR_DAMPING_ALPHA_DECREASE"
                                                    << " it " << it
                                                    << " alpha " << alpha
                                                    << " new_alpha " << newAlpha
                                                    << " previous_merit " << previousAcceptedMerit
                                                    << " trial_merit " << lastNonlinearMeritAfter
                                                    << std :: endl;
                                    }
                                }
                            } else if ( lastNonlinearMeritAfter < previousAcceptedMerit && lastNonlinearMeritAfter > 0. ) {
                                const double progressFactor = std :: min(1.5, previousAcceptedMerit / lastNonlinearMeritAfter);
                                const double newAlpha = std :: min(nonlinearDampingMax, alpha * progressFactor);
                                currentNonlinearDampingFactor = newAlpha;
                                if ( newAlpha > alpha * ( 1. + 1e-12 ) && !silent ) {
                                    std :: cout << "NONLINEAR_DAMPING_ALPHA_INCREASE"
                                                << " it " << it
                                                << " alpha " << alpha
                                                << " new_alpha " << newAlpha
                                                << " previous_merit " << previousAcceptedMerit
                                                << " trial_merit " << lastNonlinearMeritAfter
                                                << std :: endl;
                                }
                            }
                        }
                    }
                    updateAdaptiveNonlinearDamping(meritBefore, lastNonlinearMeritAfter);
                    lastNonlinearConvergenceMerit = currentNonlinearConvergenceMerit();
                }
            }

            if ( !useNonlinearGlobalization ) {
                //update DoFs
                updateFieldVariables();
                //compute residuals
                computeForcesAtIntegrationTime(false); //to obtain the actual stress, fluxes, ...

                //compute errors
                evaluateErrors();
                if ( it == 0 ) {
                    disErr = 0;                        //error in displacement change, only from second iteration
                }
            }

            if ( useNonlinearGlobalization ) {
                requestAdaptiveMatrixRebuildIfNeeded(meritBefore, lastNonlinearMeritAfter);
            }

            if ( not silent ) {
                cout << setw(6) << it << setw(15) << resErr;
            }
            if ( it == 0 ) {
                if ( not silent ) {
                    cout << setw(15) << "---";
                }
            } else {
                if ( not silent ) {
                    cout << setw(15) << disErr;
                }
            }
            if ( not silent ) {
                cout << setw(15) << eneErr;
                if ( nonlinearGlobalizationActive() ) {
                    cout << setw(15) << lastNonlinearAlpha
                         << setw(15) << lastNonlinearLineSearchTrials
                         << setw(15) << lastNonlinearMeritAfter
                         << setw(15) << lastNonlinearConvergenceMerit
                         << setw(15) << ( nonlinearCutbackReason.empty() ? "-" : nonlinearCutbackReason )
                         << setw(15) << ( lastNonlinearMatrixRebuild ? 1 : 0 )
                         << setw(15) << lastNonlinearTrustRadius;
                }
                cout << endl;
            }

            //if (not silent) checkAllVectorsForNaNs();

            // This check works only if flag "-ffast-math" is removed from CMake
            if ( std :: isnan(resErr) || std :: isnan(disErr) || std :: isnan(eneErr) ) {
                std :: cerr << "calculating with NaN in ";
                if ( std :: isnan(resErr) ) {
                    std :: cerr << "\tresiduals ";
                }
                if ( std :: isnan(disErr) ) {
                    std :: cerr << "\tdisplacements ";
                }
                if ( std :: isnan(eneErr) ) {
                    std :: cerr << "\tenergies ";
                }
                std :: cerr << endl;
                it = maxIt;
                resErr = 1e10;
                break;
            }

            const double meritAfter = useNonlinearGlobalization ? currentNonlinearGlobalizationMerit() : 0.;
            it++;
            cumul_it++;
            if ( nonlinearConvergenceCriteriaSatisfied() && it >= minIt ) {
                converged = true;
            } else {
                converged = false;
            }

            if ( useNonlinearGlobalization && !converged && shouldCutbackForNonlinearProgress(meritBefore, meritAfter) ) {
                restart_now = true;
                break;
            }

            exporters->exportData(step, it, time, 0); //to export data during iterations
        }

        if ( converged ) {
            this->fully_converged = true;
            computeForcesAtStepEnd(false); //to obtain the actual stress, fluxes, ...
        }

        if ( !converged && dt > dtmin * 1.00001 ) {
            time -= dt;
            dt = fmax(dt * critical_step_decrease, dtmin);
            trial_r = r;
            f_int = f_int_old;
            f_ext = f_ext_old;
            load.setZero(); // std :: fill(begin(load), end(load), 0);
            ddr.setZero(); // std :: fill(begin(ddr), end(ddr), 0);            //ddr *= 0;
            elems->resetMaterialStatuses();   ///> reset material internal vars to the last converged state
            if ( idc ) { //idc solver needs residuals from last converged step
                idc_time = idc_time_converged;
                nodes->addRHS_nodalLoad(load, idc_time); //add nodal load
                nodes->updateDirrichletBC(trial_r, idc_time); //give prescribed DoFs
                updateFieldVariables();      //with ddr=0
                computeForcesAtIntegrationTime(true);
            }

            time += dt;
            if ( not silent ) {
                if ( !nonlinearCutbackReason.empty() ) {
                    std :: cout << "Nonlinear cutback reason: " << nonlinearCutbackReason << endl;
                }
                std :: cout << "Restarting step, timestep = " << dt << ", time = " << time << endl;
            }
            restarts++;
            restarted = true;
            restart_now = false;
            nonlinearCutbackReason.clear();
        } else if ( !converged ) {
            if ( disErr < limitDisErr && resErr < limitResErr && eneErr < limitEneErr ) {
                if ( not silent ) {
                    std :: cout << "tolerance increased in this step (fallback acceptance, not full convergence)" << '\n';
                }
                converged = true;
                this->fully_converged = false;
                computeForcesAtStepEnd(false); //to obtain the actual stress, fluxes, ...
            } else {
                std :: cerr << "Error: " << name << " did not converge to the solution" << endl;
                terminated = true;
                return;
            }
        } else if ( ( !restarted ) && converged && it < enlargeIt ) {
            dt = fmin(dt * step_increase, dtmax);
            if ( not silent ) {
                std :: cout << "enlarging step, timestep = " << dt << '\n';
            }
        } else if ( converged && it > shortenIt && dt > dtmin ) {
            dt = fmax(dt * step_decrease, dtmin);
            if ( not silent ) {
                std :: cout << "shortening step, timestep = " << dt << '\n';
            }
        }
        if  ( dt > dtmax ) {
            dt = dtmax;
            if ( not silent ) {
                std :: cout << "shortening step to the maximum one: " << dt << '\n';
            }
        }
        if  ( dt < dtmin ) {
            dt = dtmin;
            if ( not silent ) {
                std :: cout << "enlarging step to the minimum one: " << dt << '\n';
            }
        }
    }
}


//////////////////////////////////////////////////////////
void SteadyStateNonLinearSolver :: runBeforeEachStep() {
    if ( time == init_time && idc ) {
        nodes->addRHS_nodalLoad(load, time); //add nodal load
        nodes->updateDirrichletBC(trial_r, time); //give prescribed DoFs at time 0
        computeInternalExternalForces(trial_r, load, true, time, -1.);
    }


    SteadyStateLinearSolver :: runBeforeEachStep();

    if ( not silent ) {
        cout <<  scientific; //cout << setprecision(8);
        cout << "----------------------------------------------------" << endl;
        cout << setw(6) << "iter." << setw(15) << "residual" << setw(15) << "displacement" << setw(15) << "energy error" << endl;
        if ( nonlinearGlobalizationActive() ) {
            cout << setw(6) << " " << setw(15) << "alpha" << setw(15) << "ls_trials" << setw(15) << "glob_merit" << setw(15) << "conv_merit" << setw(15) << "cutback" << setw(15) << "K_rebuild" << setw(15) << "trust_radius" << endl;
        }
        cout << setw(6) << " " << setw(15) << maxResErr << setw(15) << maxDisErr << setw(15) << maxEneErr << endl;
        cout << "----------------------------------------------------" << endl;
    }
}

//////////////////////////////////////////////////////////
void SteadyStateNonLinearSolver :: printAllVectors() {
    std :: cout << "trial_r\t r\t full_ddr\t ddr\t f_int\t f_ext\t load\t f_reduced\t" << '\n';
    for ( unsigned i = 0; i < trial_r.size(); i++ ) {
        std :: cout << trial_r [ i ] << '\t';
        std :: cout << r [ i ] << '\t';
        std :: cout << full_ddr [ i ] << '\t';
        if ( i < ddr.size() ) {
            std :: cout << ddr [ i ] << '\t';
        } else {
            std :: cout << 'H' << '\t';
        }
        std :: cout << f_int [ i ] << '\t';
        std :: cout << f_ext [ i ] << '\t';
        std :: cout << load [ i ] << '\t';
        if ( i < f.size() ) {
            std :: cout << f [ i ] << '\n';
        } else {
            std :: cout << 'H' << '\n';
        }
    }
}

//////////////////////////////////////////////////////////
void SteadyStateNonLinearSolver :: checkAllVectorsForNaNs() {
    // JK left for testing
    bool trial_r_nan, r_nan, full_ddr_nan, ddr_nan, f_int_nan, f_ext_nan, load_nan, f_nan, residuals_nan;
    trial_r_nan = r_nan = full_ddr_nan = ddr_nan = f_int_nan = f_ext_nan = load_nan = f_nan = residuals_nan = false;
    for ( unsigned i = 0; i < trial_r.size(); i++ ) {
        if ( std :: isnan(trial_r [ i ]) || std :: isinf(trial_r [ i ]) ) {
            trial_r_nan = true;
        }
        if ( std :: isnan(r [ i ]) || std :: isinf(r [ i ]) ) {
            r_nan = true;
        }
        if ( std :: isnan(full_ddr [ i ]) || std :: isinf(full_ddr [ i ]) ) {
            full_ddr_nan = true;
        }
        if ( i < ddr.size() ) {
            if ( std :: isnan(ddr [ i ]) || std :: isinf(ddr [ i ]) ) {
                ddr_nan = true;
            }
        }

        if ( std :: isnan(f_int [ i ]) || std :: isinf(f_int [ i ]) ) {
            f_int_nan = true;
        }
        if ( std :: isnan(f_ext [ i ]) || std :: isinf(f_ext [ i ]) ) {
            f_ext_nan = true;
        }
        if ( std :: isnan(load [ i ]) || std :: isinf(load [ i ]) ) {
            load_nan = true;
        }
        if ( std :: isnan(residuals [ i ]) || std :: isinf(residuals [ i ]) ) {
            residuals_nan = true;
        }
        if ( i < f.size() ) {
            if ( std :: isnan(f [ i ]) || std :: isinf(f [ i ]) ) {
                f_nan = true;
            }
        }
    }
    std :: cout << "trial_r_nan " << trial_r_nan << " r_nan " << r_nan << " full_ddr_nan " << full_ddr_nan << " ddr_nan " << ddr_nan << " f_int_nan " << f_int_nan << " f_ext_nan " << f_ext_nan << " load_nan " << load_nan << " f_nan " << f_nan << " residuals_nan " <<  residuals_nan << '\n';
}

//////////////////////////////////////////////////////////
void SteadyStateNonLinearSolver :: runAfterEachStep() {
    if ( !terminated ) {
        SteadyStateLinearSolver :: runAfterEachStep();

        if ( not silent ) {
            cout << "----------------------------------------------------" << endl;
        }

        if ( idc ) {
            idc_time_converged = idc_time;
        }
    }
}

//////////////////////////////////////////////////////////
double SteadyStateNonLinearSolver :: giveIDCtime(const bool converged) {
    if ( converged ) {
        return this->idc_time_converged;
    }
    return this->idc_time;
}

//////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////
// TRANSIENT LINEAR TRANSPORT SOLVER
TransientLinearTransportSolver :: TransientLinearTransportSolver() {
    name = "TransientLinearTransportSolver";
    timeIntM = 0; //generalized-alpha
    isTimeReal = true;
    setDefaultIntegrationParams();
    check_time_integr_params = true;
    dampingMatrixIterUpdate = -1;
    dampingMatrixCumulIterUpdate = -1;
    dampingMatrixStepUpdate = -1;
    stiffMatType = "elastic";
}

//////////////////////////////////////////////////////////
void TransientLinearTransportSolver :: applySpectralRadius(double rhoinfty) {
    //set up the generalized-alpha method according to Chung and Hulbert 1993
    if ( abs(rhoinfty - 0.5) > 0.5 ) {
        cerr << "Error in solver: spectral radius must be inside interval [0,1]" << endl;
        exit(1);
    }

    //according to Chung and Hulbert, 1993, JAM
    alpha_m = ( 2. * rhoinfty - 1. ) / ( 1. + rhoinfty );
    alpha_f = rhoinfty / ( 1. + rhoinfty );
    gamma = 1. / 2. - alpha_m + alpha_f;
}

//////////////////////////////////////////////////////////
void TransientLinearTransportSolver :: checkIntegrationParams() {
    if ( !check_time_integr_params ) {
        return;
    }

    if ( timeIntM == 0 ) {  //generalized alpha
        //TODO
    } else if ( timeIntM == 2 ) {   //HHT method
        //TODO
    } else if ( timeIntM == 3 ) {   //Newmark method
        if ( alpha_m != 0 || alpha_f != 0 ) {
            cerr << "Solver Error: Newmark method requires alpha_m=alpha_f = 0" << endl;
            exit(1);
        }
        if ( abs(gamma - 0.5) > 0.5 ) {
            cerr << "Solver Error: Newmark method requires gamma withn interva 0-1" << endl;
            exit(1);
        }
        if ( abs(beta - 0.25) > 0.25 ) {
            cerr << "Solver Error: Newmark method requires gamma withn interva 0-0.5" << endl;
            exit(1);
        }
    } else {
        cerr << "Solver Error: unknowns method for time integration" << endl;
        exit(1);
    }
}

//////////////////////////////////////////////////////////
void TransientLinearTransportSolver :: setDefaultIntegrationParams() {
    if ( timeIntM == 0 ) {  //generalized alpha
        applySpectralRadius(0.8);
    } else if ( timeIntM == 2 ) {   //HHT method
        alpha_m = 0;
        alpha_f = 0;
        gamma = 0.5;
    } else if ( timeIntM == 3 ) {   //Newmark method
        alpha_m = alpha_f = 0;
        gamma = 0.5;
    } else {
        cerr << "Solver Error: unknowns method for time integration" << endl;
        exit(1);
    }
}

//////////////////////////////////////////////////////////
TransientLinearTransportSolver :: ~TransientLinearTransportSolver() {}

//////////////////////////////////////////////////////////
void TransientLinearTransportSolver :: prepareSystemMatricesAndInitialField(string init_r_file, string init_v_file, const bool initial) {
    v_old = Vector :: Zero(totalDoFnum);
    elems->prepareDampingMatrix(C);

    SteadyStateLinearSolver :: prepareSystemMatricesAndInitialField(init_r_file, init_v_file, initial);
}

//////////////////////////////////////////////////////////
void TransientLinearTransportSolver :: init(string init_r_file, string init_v_file, const bool initial) {
    Solver :: init(init_r_file, init_v_file, initial);

    checkIntegrationParams();

    prepareSystemMatricesAndInitialField(init_r_file, init_v_file, initial);
    computeKeff();

    nodes->giveReducedForceArray(residuals, f);
    CoordinateIndexedSparseMatrix Cred(C);
    if ( nodes->giveConstraints()->isActive() ) {
        nodes->giveConstraints()->transformToConstraintSpace(Cred);
    }
    //not good, might have zeros at the diagonal
    //StandAloneLinalgSolver(Cred, ddr, f, conj_grad_precision, conj_grad_relative_maxit, symsolver_type);
    
    //not ideal but working
    double k;
    for(unsigned i=0; i<Cred.rows(); i++){
        k = Cred.coeff(i,i);
        if (k == 0) ddr[i] = 0;
        else ddr[i] = f[i]/k;
    }
    nodes->giveFullDoFArray(ddr, v);
}

//////////////////////////////////////////////////////////
void TransientLinearTransportSolver :: rebuild() {
    Solver :: rebuild();
    prepareSystemMatricesAndInitialField("", "", false);
    computeKeff();
}

//////////////////////////////////////////////////////////
void TransientLinearTransportSolver :: computeForcesAtIntegrationTime(const bool frozen) {
    elems->integrateDampingForces(v * ( 1. - alpha_m ) +  v_old * alpha_m, f_dam);
    Vector ll = load_old * alpha_f + load * ( 1. - alpha_f );
    computeInternalExternalForces( r * alpha_f + trial_r * ( 1. - alpha_f ), ll, frozen, time - dt*alpha_f, dt * ( 1. - alpha_f ) );
    residuals -= f_dam;
}

//////////////////////////////////////////////////////////
void TransientLinearTransportSolver :: computeForcesAtStepEnd(const bool frozen) {
    elems->integrateDampingForces(v, f_dam);
    computeInternalExternalForces(trial_r, load, frozen, time, dt);
    residuals -= f_dam;
}

//////////////////////////////////////////////////////////
void TransientLinearTransportSolver :: computeKeff() {
    Keff = C * ( ( 1. - alpha_m ) / ( dt * gamma ) ) + K * ( 1. - alpha_f );
    if ( nodes->giveConstraints()->isActive() ) {
        nodes->giveConstraints()->transformToConstraintSpace(Keff);
    }
    factorizeLinearSystem();
}

//////////////////////////////////////////////////////////
//void TransientLinearTransportSolver :: updateFeff() {
//    nodes->giveReducedDoFArray(v,v_red);
//    feff =  f - ( C * v_red ) * ( 1. + ( alpha_m - 1. ) / gamma ); //not used in the code, available for checking
//}

//////////////////////////////////////////////////////////
void TransientLinearTransportSolver :: updateFieldVariables() {
    SteadyStateLinearSolver :: updateFieldVariables();
    v = ( trial_r - r ) * ( 1. / ( dt * gamma ) ) + v_old * ( 1. - 1. / gamma );
}

//////////////////////////////////////////////////////////
void TransientLinearTransportSolver :: solve() {
    SteadyStateLinearSolver :: solve();
}

//////////////////////////////////////////////////////////
void TransientLinearTransportSolver :: runBeforeEachStep() {
    SteadyStateLinearSolver :: runBeforeEachStep();
    v_old = v;
}

//////////////////////////////////////////////////////////
void TransientLinearTransportSolver :: runAfterEachStep() {
    SteadyStateLinearSolver :: runAfterEachStep();
};

//////////////////////////////////////////////////////////
Solver *TransientLinearTransportSolver :: readFromFile(const string filename) {
    SteadyStateNonLinearSolver :: readFromFile(filename);



    double num;
    int valueIN;
    string param, line;
    ifstream inputfile( filename.c_str() );
    if ( inputfile.is_open() ) {
        while ( getline(inputfile >> std :: ws, line) ) {
            if ( line.empty() || ( line.at(0) == '#' ) ) {
                continue;
            }
            istringstream iss(line);
            iss >> param;
            if ( param.compare("integration_method") == 0 ) {
                iss >> param;
                if ( param.compare("newmark") == 0 ) {
                    timeIntM = 3;
                } else if ( param.compare("hht") == 0 ) {
                    timeIntM = 2;
                } else if ( param.compare("wbz") == 0 ) {
                    timeIntM = 1;
                } else if ( param.compare("generalized_alpha") == 0 ) {
                    timeIntM = 0;
                } else {
                    cerr << "Error: unknown integration method " << param << endl;
                    exit(1);
                }
                setDefaultIntegrationParams();
            } else if ( param.compare("spectral_radius") == 0 ) {
                iss >> num;
                applySpectralRadius(num);
            } else if ( param.compare("gamma") == 0 ) {
                iss >> gamma;
            } else if ( param.compare("beta") == 0 ) {
                iss >> beta;
            } else if ( param.compare("alpha_m") == 0 ) {
                iss >> alpha_m;
            } else if ( param.compare("alpha_f") == 0 ) {
                iss >> alpha_f;
            } else if ( param.compare("do_not_check_time_integration_params") == 0 ) {
                check_time_integr_params = false;
            } else if ( param.compare("damping_matrix_update") == 0 || param.compare("damping_matrix_iter_update") == 0  ) {
                iss >> valueIN;
                dampingMatrixIterUpdate = int( valueIN );
            } else if ( param.compare("damping_matrix_cumul_iter_update") == 0  ) {
                iss >> valueIN;
                dampingMatrixCumulIterUpdate = int( valueIN );
            } else if ( param.compare("damping_matrix_step_update") == 0  ) {
                iss >> valueIN;
                dampingMatrixStepUpdate = int( valueIN );
            }
        }
        inputfile.close();
    }
    return this;
};

//////////////////////////////////////////////////////////
bool TransientLinearTransportSolver :: updateSystemMatrices(unsigned iteration, unsigned cumul_iteration, bool enforce) {
    bool updated0 = SteadyStateNonLinearSolver :: updateSystemMatrices(iteration, cumul_iteration, enforce);
    bool updated1 = false;
    if ( enforce || ( iteration == 0 && ( dampingMatrixStepUpdate == 0 || ( dampingMatrixStepUpdate > 0 && step % dampingMatrixStepUpdate == 0 ) ) ) || dampingMatrixIterUpdate == 0 || ( dampingMatrixIterUpdate > 0 && iteration % dampingMatrixIterUpdate == 0 ) || ( dampingMatrixCumulIterUpdate > 0 && cumul_iteration % dampingMatrixCumulIterUpdate == 0 ) ) {
        elems->updateDampingMatrix(C);
        updated1 = true;
    }
    return ( updated0 || updated1 );
}

//////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////
// TRANSIENT NON-LINEAR TRANSPORT SOLVER
TransientNonLinearTransportSolver :: TransientNonLinearTransportSolver() {
    name = "TransientNonLinearTransportSolver";
    stiffMatType = "secant";
}

//////////////////////////////////////////////////////////
TransientNonLinearTransportSolver :: ~TransientNonLinearTransportSolver() {}

//////////////////////////////////////////////////////////
void TransientNonLinearTransportSolver :: init(string init_r_file, string init_v_file, const bool initial) {
    TransientLinearTransportSolver :: init(init_r_file, init_v_file, initial);
}


//////////////////////////////////////////////////////////
void TransientNonLinearTransportSolver :: runBeforeEachStep() {
    SteadyStateNonLinearSolver :: runBeforeEachStep();
    v_old = v;
}

//////////////////////////////////////////////////////////
void TransientNonLinearTransportSolver :: runAfterEachStep() {
    SteadyStateNonLinearSolver :: runAfterEachStep();
};

//////////////////////////////////////////////////////////
void TransientNonLinearTransportSolver :: solve() {
    SteadyStateNonLinearSolver :: solve();
}

//////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////
// TRANSIENT LINEAR MECHANICAL SOLVER
TransientLinearMechanicalSolver :: TransientLinearMechanicalSolver() {
    name = "TransientLinearMechanicalSolver";
    setDefaultIntegrationParams(); //this always call method from TransientLinearMechanicalSolver
    lumpMassM = false;
    massMatrixIterUpdate = -1;
    massMatrixCumulIterUpdate = -1;
    massMatrixStepUpdate = -1;
    stiffMatType = "elastic";
}

//////////////////////////////////////////////////////////
TransientLinearMechanicalSolver :: ~TransientLinearMechanicalSolver() {}

//////////////////////////////////////////////////////////
void TransientLinearMechanicalSolver :: solve() {
    SteadyStateLinearSolver :: solve();
}

//////////////////////////////////////////////////////////
Solver *TransientLinearMechanicalSolver :: readFromFile(const string filename) {
    TransientLinearTransportSolver :: readFromFile(filename);

    int valueIN;
    string param, line;
    ifstream inputfile( filename.c_str() );
    if ( inputfile.is_open() ) {
        while ( getline(inputfile >> std :: ws, line) ) {
            if ( line.empty() || ( line.at(0) == '#' ) ) {
                continue;
            }
            istringstream iss(line);
            iss >> param;
            if ( param.compare("use_lumped_mass_matrix") == 0 ) {
                lumpMassM = true;
            } else if ( param.compare("mass_matrix_update") == 0 || param.compare("mass_matrix_iter_update") == 0  ) {
                iss >> valueIN;
                massMatrixIterUpdate = int( valueIN );
            } else if ( param.compare("mass_matrix_cumul_iter_update") == 0  ) {
                iss >> valueIN;
                massMatrixCumulIterUpdate = int( valueIN );
            } else if ( param.compare("mass_matrix_step_update") == 0  ) {
                iss >> valueIN;
                massMatrixStepUpdate = int( valueIN );
            }
        }
        inputfile.close();
    }
    return this;
};

//////////////////////////////////////////////////////////
void TransientLinearMechanicalSolver :: prepareSystemMatricesAndInitialField(string init_r_file, string init_v_file, const bool initial) {
    v_old = Vector :: Zero(totalDoFnum);
    a_old = Vector :: Zero(totalDoFnum);
    elems->prepareMassMatrix(M, lumpMassM);

    TransientLinearTransportSolver :: prepareSystemMatricesAndInitialField(init_r_file, init_v_file, initial);
}

//////////////////////////////////////////////////////////
void TransientLinearMechanicalSolver :: init(string init_r_file, string init_v_file, const bool initial) {
    Solver :: init(init_r_file, init_v_file, initial);

    checkIntegrationParams();

    prepareSystemMatricesAndInitialField(init_r_file, init_v_file, initial);
    computeKeff();

    //compute initial acceleration
    nodes->giveReducedForceArray(residuals, f);
    CoordinateIndexedSparseMatrix Cred(C);
    CoordinateIndexedSparseMatrix Mred(M);
    if ( nodes->giveConstraints()->isActive() ) {
        nodes->giveConstraints()->transformToConstraintSpace(Cred);
        nodes->giveConstraints()->transformToConstraintSpace(Mred);
    }
    Vector v_red = Vector :: Zero(freeDoFnum);
    nodes->giveReducedDoFArray(v, v_red);
    
    //not good, might have zeros at the diagonal
    //StandAloneLinalgSolver(Mred, ddr, f - Cred * v_red, conj_grad_precision, conj_grad_relative_maxit, symsolver_type)
    
    //not ideal but working
    double k;
    Vector q =  f - Cred * v_red;
    for(unsigned i=0; i<Mred.rows(); i++){
        k = Mred.coeff(i,i);
        if (k == 0) ddr[i] = 0;
        else ddr[i] = q[i]/k;
    }
    a = Vector :: Zero(totalDoFnum);
    nodes->giveFullDoFArray(ddr, a);
}

//////////////////////////////////////////////////////////
void TransientLinearMechanicalSolver :: applySpectralRadius(double rhoinfty) {
    //set up the generalized-alpha method according to Chung and Hulbert 1993
    if ( rhoinfty < 0 || rhoinfty > 1 ) {
        cerr << "Error in solver: spectral radius must be inside interval 0-1" << endl;
        exit(1);
    }
    if ( timeIntM == 0 ) { //generalized alpha
        //according to Chung and Hulbert, 1993, JAM
        alpha_m = ( 2. * rhoinfty - 1. ) / ( 1. + rhoinfty );
        alpha_f = rhoinfty / ( 1. + rhoinfty );
        gamma = 1. / 2. - alpha_m + alpha_f;
        beta = 1. / 4. * pow(1 - alpha_m + alpha_f, 2);
    } else if ( timeIntM == 1 ) { //WBZ method
        //Stability Analysis of Ubiquitous Direct Time Integration Methods, Mohamed Naguib and AF Ghaleb and Faraji Mollaie Amin,  2019
        alpha_m = ( rhoinfty - 1. ) / ( 1. + rhoinfty );
        alpha_f = 0;
        gamma = 1. / 2. - alpha_m;
        beta = 1. / 4. * pow(1 - alpha_m, 2);
    } else if ( timeIntM == 2 ) {  //HHT
        //Stability Analysis of Ubiquitous Direct Time Integration Methods, Mohamed Naguib and AF Ghaleb and Faraji Mollaie Amin,  2019
        alpha_m = 0.;
        alpha_f = ( 1. - rhoinfty ) / ( 1. + rhoinfty );
        gamma = 1. / 2. + alpha_f;
        beta = 1. / 4. * pow(1. + alpha_f, 2);
    } else if ( timeIntM == 3 ) {  //Newmark method
        //Stability Analysis of Ubiquitous Direct Time Integration Methods, Mohamed Naguib and AF Ghaleb and Faraji Mollaie Amin,  2019
        //check also Klaus-Jürgen Bathe and Gunwoo Noh 2012
        alpha_m = 0;
        alpha_f = 0;
        gamma = ( 3. - rhoinfty ) / ( 2. * rhoinfty + 2. );
        beta = 1. / pow(rhoinfty + 1., 2);
    } else {
        cerr << "Solver Error: unknown method for time integration" << endl;
        exit(1);
    }
}

//////////////////////////////////////////////////////////
void TransientLinearMechanicalSolver :: checkIntegrationParams() {
    if ( !check_time_integr_params ) {
        return;
    }
    if ( timeIntM == 0 ) {  //generalized alpha
        if ( alpha_m > 0.5 || alpha_f > 0.5 || alpha_m < 0 || alpha_f < 0 ) {
            cerr << "Solver Error: Generalized-alpha method requires alpha_m and alpha_f within interval 0-0.5, instead these parameters are " << alpha_m << " and " << alpha_f << endl;
            exit(1);
        }
        if ( abs(gamma - 0.5) > 0.5 ) {
            cerr << "Solver Error: Generalized-alpha method requires gamma withn interva 0-1, but gamma is set to " << gamma << endl;
            exit(1);
        }
        if ( beta > 0.5 || beta < 0.25 + 0.5 * ( alpha_f - alpha_m ) ) {
            cerr << "Solver Error: Generalized-alpha method requires beta withn interva 0.25+0.5(alpha_f-alpha_m), but gamma is set to " << beta << endl;
            exit(1);
        }
    } else if ( timeIntM == 1 ) {   //WBZ method
        ////
    } else if ( timeIntM == 2 ) {   //HHT method
        if ( alpha_m != 0 ) {
            cerr << "Solver Error: HHT method requires alpha_m=0" << endl;
            exit(1);
        }
        if ( alpha_f > 1. / 3. || alpha_f < 0 ) {
            cerr << "Solver Error: HHT method requires alpha_f inside interva 0-1./3" << endl;
            exit(1);
        }
        if ( abs(gamma - 0.5) > 0.5 ) {
            cerr << "Solver Error: HHT method requires gamma withn interva 0-1" << endl;
            exit(1);
        }
        if ( beta > 0.5 || beta < 0.25 + 0.5 * ( alpha_f ) ) {
            cerr << "Solver Error: HHT method requires beta withn interva 0.25+0.5*alpha_f" << endl;
            exit(1);
        }
    } else if ( timeIntM == 3 ) {   //Newmark method
        if ( alpha_m != 0 || alpha_f != 0 ) {
            cerr << "Solver Error: Newmark method requires alpha_m=alpha_f = 0, instead these parameters are " << alpha_m << " and " << alpha_f << endl;
            exit(1);
        }
        if ( abs(gamma - 0.5) > 0.5 ) {
            cerr << "Solver Error: Newmark method requires gamma withn interva 0-1, but gamma is set to " << gamma << endl;
            exit(1);
        }
        if ( abs(beta - 0.25) > 0.25 ) {
            cerr << "Solver Error: Newmark method requires beta withn interva 0-0.5, but beta is set to " << beta << endl;
            exit(1);
        }
    } else {
        cerr << "Solver Error: unknowns method for time integration" << endl;
        exit(1);
    }
}
//////////////////////////////////////////////////////////
void TransientLinearMechanicalSolver :: setDefaultIntegrationParams() {
    applySpectralRadius(0.8);
}

//////////////////////////////////////////////////////////
void TransientLinearMechanicalSolver :: computeKeff() {
    Keff = K * ( 1. - alpha_f ) + C * ( ( 1 - alpha_f ) * gamma / dt / beta ) + M * ( ( 1 - alpha_m ) / dt / dt / beta );
    if ( nodes->giveConstraints()->isActive() ) {
        nodes->giveConstraints()->transformToConstraintSpace(Keff);
    }
    factorizeLinearSystem();
}

//////////////////////////////////////////////////////////
//void TransientLinearMechanicalSolver :: updateFeff() {
//feff = f - C * ( ( 1. + gamma * ( alpha_f - 1 ) / beta ) * v_red + dt * ( 1. - alpha_f ) * ( 1. - gamma / ( 2. * beta ) ) * a_red ) - M * ( ( alpha_m - 1. ) / ( dt * beta ) * v_red  + ( 1. + ( alpha_m - 1 ) / ( 2. * beta ) ) * a_red );
//}

//////////////////////////////////////////////////////////
void TransientLinearMechanicalSolver :: updateFieldVariables() {
    SteadyStateLinearSolver :: updateFieldVariables();
    v = ( trial_r - r ) * ( gamma / ( dt * beta ) ) + v_old * ( 1. - gamma / beta ) + a_old * ( dt * ( 1. - gamma / ( 2. * beta ) ) );
    a = ( ( trial_r - r ) - v_old * dt - a_old * ( pow(dt, 2) * ( 1. / 2. - beta ) ) ) / ( pow(dt, 2) * beta );
}


//////////////////////////////////////////////////////////
void TransientLinearMechanicalSolver :: computeForcesAtIntegrationTime(const bool frozen) {
    elems->integrateDampingForces(v * ( 1. - alpha_f ) +  v_old * alpha_f, f_dam);
    elems->integrateInertiaForces(a * ( 1. - alpha_m ) +  a_old * alpha_m, f_acc);
    Vector ll = load_old * alpha_f + load * ( 1. - alpha_f );
    computeInternalExternalForces( r * alpha_f + trial_r * ( 1. - alpha_f ), ll, frozen, time-dt*alpha_f, dt * ( 1. - alpha_f ) );
    residuals -= f_dam + f_acc;
}

//////////////////////////////////////////////////////////
void TransientLinearMechanicalSolver :: computeForcesAtStepEnd(const bool frozen) {
    elems->integrateDampingForces(v, f_dam);
    elems->integrateInertiaForces(a, f_acc);
    computeInternalExternalForces(trial_r, load, frozen, time, dt);
    residuals -= f_dam + f_acc;
}

//////////////////////////////////////////////////////////
void TransientLinearMechanicalSolver :: runBeforeEachStep() {
    TransientLinearTransportSolver :: runBeforeEachStep();
    a_old = a;
}

//////////////////////////////////////////////////////////
void TransientLinearMechanicalSolver :: runAfterEachStep() {
    TransientLinearTransportSolver :: runAfterEachStep();
};


//////////////////////////////////////////////////////////
bool TransientLinearMechanicalSolver :: updateSystemMatrices(unsigned iteration, unsigned cumul_iteration, bool enforce) {
    bool updated0 = TransientLinearTransportSolver :: updateSystemMatrices(iteration, cumul_iteration, enforce);
    bool updated1 = false;
    if ( enforce || ( iteration == 0 && ( massMatrixStepUpdate == 0 || ( massMatrixStepUpdate > 0 && step % massMatrixStepUpdate == 0 ) ) )  || massMatrixIterUpdate == 0 || ( massMatrixIterUpdate > 0 && iteration % massMatrixIterUpdate == 0 ) || ( massMatrixCumulIterUpdate > 0 && cumul_iteration % massMatrixCumulIterUpdate == 0 ) ) {
        elems->updateMassMatrix(M, lumpMassM);
        if ( lumpMassM ) {
            elems->replaceTrueMassMatricesByLumpedOnes();
        }
        updated1 = true;
    }
    return ( updated0 || updated1 );
}

//////////////////////////////////////////////////////////
void TransientLinearMechanicalSolver :: computeTotalKineticEnergy() {
    W_kin [ 0 ] = elems->integrateKineticEnergy(v);
}



//////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////
// TRANSIENT NON LINEAR MECHANICAL SOLVER
TransientNonLinearMechanicalSolver :: TransientNonLinearMechanicalSolver() {
    name = "TransientNonLinearMechanicalSolver";
    stiffMatType = "secant";
}

//////////////////////////////////////////////////////////
TransientNonLinearMechanicalSolver :: ~TransientNonLinearMechanicalSolver() {}

//////////////////////////////////////////////////////////
void TransientNonLinearMechanicalSolver :: init(string init_r_file, string init_v_file, const bool initial) {
    TransientLinearMechanicalSolver :: init(init_r_file, init_v_file, initial);
}

//////////////////////////////////////////////////////////
void TransientNonLinearMechanicalSolver :: runBeforeEachStep() {
    TransientNonLinearTransportSolver :: runBeforeEachStep();
    a_old = a;
}

//////////////////////////////////////////////////////////
void TransientNonLinearMechanicalSolver :: runAfterEachStep() {
    TransientNonLinearTransportSolver :: runAfterEachStep();
};


//////////////////////////////////////////////////////////
void TransientNonLinearMechanicalSolver :: solve() {
    SteadyStateNonLinearSolver :: solve();
}
