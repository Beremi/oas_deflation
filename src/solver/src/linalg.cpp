#include "linalg.h"

using namespace std;
using namespace Eigen;

bool LinalgSymmetricSolver(const CoordinateIndexedSparseMatrix &A, Vector &x, const Vector &b, const Vector x0, double precision, double relmaxit) {
#if PRINT_DEBUG_TIME
    auto start = std :: chrono :: system_clock :: now();
#endif
    if ( b.size() == 0 ) {
        return true;                   // when problem is completely constrained (e.g. single facet)
    }
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
            tripletList.push_back(T(i, A.column_index [ idx ], A.array [ idx ]) );
            idx++;
        }
        ;
    }
    ;
    mat.setFromTriplets(tripletList.begin(), tripletList.end() );
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

    ConjugateGradient< SparseMatrix< double >, Lower | Upper >cgK;
    cgK.setMaxIterations(Maxit);
    cgK.setTolerance(precision);

    cgK.compute(mat);

    VectorXd cgx = cgK.solveWithGuess(cgb, cgx0);
    //VectorXd cgx = cgK.solve(cgb);

#if PRINT_DEBUG_TIME
    now = std :: chrono :: system_clock :: now();

    elapsed_seconds = now - start;
    std :: cout << "linalg solver duration: " << convertTimeToString(elapsed_seconds) << " nit: " << cgK.iterations() << std :: endl;
    cout.flush();
#endif

    for ( size_t i = 0; i < rowsize; i++ ) {
        x [ i ] = cgx [ i ];
    }

    bool result = size_t(cgK.iterations() ) < Maxit;
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
            tripletList.push_back(T(i, A.column_index [ idx ], A.array [ idx ]) );
            idx++;
        }
        ;
    }
    ;
    mat.setFromTriplets(tripletList.begin(), tripletList.end() );
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
    bool result = size_t(bicg.iterations() ) < Maxit;
    return result;
}
