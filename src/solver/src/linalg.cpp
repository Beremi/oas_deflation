#include "linalg.h"

using namespace std;
using namespace Eigen;

bool LinalgSymmetricSolver(const CoordinateIndexedSparseMatrix &A, Vector &x, const Vector &b, const Vector x0, double precision, double relmaxit, string solver_type) {
#if PRINT_DEBUG_TIME
    auto start = std :: chrono :: system_clock :: now();
#endif
    if ( b.size() == 0 ) {
        return true;                   // when problem is completely constrained (e.g. single facet)
    }

    if ( solver_type == "ConjGrad" ) {
        bool result = ConjGrad(A, x, b, x0, precision, relmaxit);

#if PRINT_DEBUG_TIME
        auto now = std :: chrono :: system_clock :: now();

        auto elapsed_seconds = now - start;
        std :: cout << "linalg solver preprocessing: " << convertTimeToString_(elapsed_seconds) << endl;
        cout.flush();
#endif
        return result;
    }

    size_t Maxit = fmax(b.size() * relmaxit, 1);

    auto rowsize = A.row_size.size();
    Eigen :: SparseMatrix< double, RowMajor >mat(rowsize, rowsize);

    typedef Eigen :: Triplet< double >T;
    std :: vector< T >tripletList;
    //tripletList.reserve(data.size());
    size_t idx = 0;
    for ( size_t i = 0; i < rowsize; i++ ) {
        for ( size_t c = 0; c < A.row_size [ i ]; c++ ) {
            //cout<<i << "--" <<  col[c]<<endl;
            //mat.insert(i, A.xx[idx]) = data[idx];
            tripletList.push_back( T(i, A.column_index [ idx ], A.array [ idx ]) );
            idx++;
        }
        ;
    }
    ;
    mat.setFromTriplets( tripletList.begin(), tripletList.end() );
    mat.makeCompressed();

    VectorXd cgb(rowsize);
    for ( size_t i = 0; i < rowsize; i++ ) {
        cgb [ i ] = b [ i ];
    }

#if PRINT_DEBUG_TIME
    auto now = std :: chrono :: system_clock :: now();
    auto elapsed_seconds = now - start;
    std :: cout << "linalg solver preprocessing: " << convertTimeToString_(elapsed_seconds) << endl;
    cout.flush();
#endif

    bool result = false;
    VectorXd cgx;
    if ( solver_type == "EigenConj" ) {
        VectorXd cgx0(rowsize);
        for ( size_t i = 0; i < rowsize; i++ ) {
            cgx0 [ i ] = x0 [ i ];
        }


        //JacobiSVD<MatrixXd> svd(cgb);
        //double cond = svd.singularValues()(0) / svd.singularValues()(svd.singularValues().size()-1);
        //cout << "condition number is " << cond<< " " << svd.singularValues()(0) << " " << svd.singularValues()(svd.singularValues().size()-1) << endl;

        ConjugateGradient< SparseMatrix< double >, Lower | Upper >cgK;
        //ConjugateGradient< SparseMatrix< double >, Lower | Upper, IncompleteCholesky< double > >cgK;
        cgK.setMaxIterations(Maxit);
        cgK.setTolerance(precision);

        cgK.compute(mat);

        cgx = cgK.solveWithGuess(cgb, cgx0);
        //VectorXd cgx = cgK.solve(cgb);
        result = size_t( cgK.iterations() ) < Maxit;
    } else if ( solver_type == "EigenLDLT" ) {
        SimplicialLDLT< SparseMatrix< double > >simplicial_ldlt_solver;
        cgx = simplicial_ldlt_solver.compute(mat).solve(cgb);
        cout << "error " << ( mat * cgx - cgb ).lpNorm< Infinity >() << endl;
        result = ( mat * cgx - cgb ).lpNorm< Infinity >() < precision;
    } else if ( solver_type == "EigenLLT" ) {
        SimplicialLLT< SparseMatrix< double > >simplicial_llt_solver;
        cgx = simplicial_llt_solver.compute(mat).solve(cgb);
        result = ( mat * cgx - cgb ).lpNorm< Infinity >() < precision;
    } else if ( solver_type == "EigenSparseLU" ) {
        SparseLU< SparseMatrix< double >, COLAMDOrdering< int > >sparseLU_solver;
        sparseLU_solver.analyzePattern(mat);
        sparseLU_solver.factorize(mat);

        cgx = sparseLU_solver.solve(cgb);
        result = ( mat * cgx - cgb ).lpNorm< Infinity >() < precision;
    }

    for ( size_t i = 0; i < rowsize; i++ ) {
        x [ i ] = cgx [ i ];
    }

#if PRINT_DEBUG_TIME
    now = std :: chrono :: system_clock :: now();

    elapsed_seconds = now - start;
    std :: cout << "linalg solver duration: " << convertTimeToString_(elapsed_seconds) << std :: endl;
    cout.flush();
#endif


    return result;
}


bool LinalgNonSymmetricSolver(const CoordinateIndexedSparseMatrix &A, Vector &x, const Vector &b, const Vector x0, double precision, double relmaxit) {
#if PRINT_DEBUG_TIME
    auto start = std :: chrono :: system_clock :: now();
#endif

    //return ConjGrad(A, x, b, x0, precision, relmaxit);

    size_t Maxit = b.size() * relmaxit;

    auto rowsize = A.row_size.size();
    Eigen :: SparseMatrix< double, RowMajor >mat(rowsize, rowsize);

    typedef Eigen :: Triplet< double >T;
    std :: vector< T >tripletList;
    //tripletList.reserve(data.size());
    size_t idx = 0;
    for ( size_t i = 0; i < rowsize; i++ ) {
        for ( size_t c = 0; c < A.row_size [ i ]; c++ ) {
            //cout<<i << "--" <<  col[c]<<endl;
            //mat.insert(i, A.xx[idx]) = data[idx];
            tripletList.push_back( T(i, A.column_index [ idx ], A.array [ idx ]) );
            idx++;
        }
        ;
    }
    ;
    mat.setFromTriplets( tripletList.begin(), tripletList.end() );
    mat.makeCompressed();

    VectorXd cgb(rowsize);
    for ( size_t i = 0; i < rowsize; i++ ) {
        cgb [ i ] = b [ i ];
    }

    VectorXd cgx0(rowsize);
    for ( size_t i = 0; i < rowsize; i++ ) {
        cgx0 [ i ] = x0 [ i ];
    }

#if PRINT_DEBUG_TIME
    auto now = std :: chrono :: system_clock :: now();

    auto elapsed_seconds = now - start;
    std :: cout << "linalg solver preprocessing: " << convertTimeToString(elapsed_seconds) << endl;
    cout.flush();
#endif

    //BiCGSTAB<SparseMatrix<double>, Eigen::IncompleteLUT<double>> bicg;
    BiCGSTAB< SparseMatrix< double > >bicg;
    bicg.setMaxIterations(Maxit);
    bicg.setTolerance(precision);

    bicg.compute(mat);

    VectorXd cgx = bicg.solveWithGuess(cgb, cgx0);
    //VectorXd cgx = bicg.solve(cgb);

#if PRINT_DEBUG_TIME
    now = std :: chrono :: system_clock :: now();

    elapsed_seconds = now - start;
    std :: cout << "linalg solver duration: " << convertTimeToString(elapsed_seconds) << " nit: " << bicg.iterations() << std :: endl;
    cout.flush();
#endif

    for ( size_t i = 0; i < rowsize; i++ ) {
        x [ i ] = cgx [ i ];
    }
    bool result = size_t( bicg.iterations() ) < Maxit;
    return result;
}

//sorterd eigenvalues and eigenvectors
bool LinalgEigenSolver(const Vector &A, Vector &eigenvalues, vector< Vector > &eigenvectors) {
    unsigned size = 3;

    MatrixXd mat = MatrixXd :: Random(size, size);
    mat(0, 0) = A [ 0 ];
    mat(1, 1) = A [ 1 ];
    mat(2, 2) = A [ 2 ];
    mat(2, 1) = mat(1, 2) = ( A [ 4 ] + A [ 3 ] ) / 2.;
    mat(2, 0) = mat(0, 2) = ( A [ 6 ] + A [ 5 ] ) / 2.;
    mat(1, 0) = mat(0, 1) = ( A [ 8 ] + A [ 7 ] ) / 2.;



    EigenSolver< MatrixXd >es(mat);
    vector< double >eigenvalsvector(size);
    eigenvalues.resize(size);
    eigenvectors.resize(size);

    for ( unsigned i = 0; i < size; i++ ) {
        eigenvalsvector [ i ] = real(es.eigenvalues() [ i ]);
    }

    // initialize original index locations
    vector< size_t >idx( eigenvalsvector.size() );
    iota(idx.begin(), idx.end(), 0);
    stable_sort(idx.begin(), idx.end(), [ & eigenvalsvector ](size_t i1, size_t i2) {
        return eigenvalsvector [ i1 ] > eigenvalsvector [ i2 ];
    });

    for ( unsigned i = 0; i < size; i++ ) {
        eigenvalues [ i ] = eigenvalsvector [ idx [ i ] ];
        eigenvectors [ i ].resize(size);
        VectorXcd v = es.eigenvectors().col(idx [ i ]);
        for ( unsigned j = 0; j < size; j++ ) {
            eigenvectors [ i ] [ j ] = real(v [ j ]);
        }
    }

    return true;
}
