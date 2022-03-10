#include "linalg.h"

using namespace std;
using namespace Eigen;

bool LinalgSymmetricSolver(const CoordinateIndexedSparseMatrix &A, MyVector &x, const MyVector &b, const MyVector &x0, double precision, double relmaxit, string solver_type) {
    cout << "A: " << A << endl;
    cout << "b: " << b << endl;
    cout << "x0: " << x0 << endl;

#if PRINT_DEBUG_TIME
    auto start = std :: chrono :: system_clock :: now();
#endif
    if ( b.size() == 0 ) {
        return true;                   // when problem is completely constrained (e.g. single facet)
    }

    size_t Maxit = fmax(b.size() * relmaxit, 1);

    bool result = false;
    if ( solver_type == "EigenConj" ) {
        //JacobiSVD<MatrixXd> svd(cgb);
        //double cond = svd.singularValues()(0) / svd.singularValues()(svd.singularValues().size()-1);
        //cout << "condition number is " << cond<< " " << svd.singularValues()(0) << " " << svd.singularValues()(svd.singularValues().size()-1) << endl;

        ConjugateGradient< SparseMatrix< double >, Lower | Upper >cgK;
        //ConjugateGradient< SparseMatrix< double >, Lower | Upper, IncompleteCholesky< double > >cgK;
        cgK.setMaxIterations(Maxit);
        cgK.setTolerance(precision);

        cgK.compute(A);

        x = cgK.solveWithGuess(b, x0);
        //VectorXd cgx = cgK.solve(cgb);
        result = size_t(cgK.iterations() ) < Maxit;
        if ( !result ) {
            cout << "Eigen Conjugate Gradients performed " << cgK.iterations() << " iterations and reached error " << cgK.error() << ", required precision is " << precision << endl;
        }
    } else if ( solver_type == "EigenLDLT" ) {
        SimplicialLDLT< SparseMatrix< double > >simplicial_ldlt_solver;
        x = simplicial_ldlt_solver.compute(A).solve(b);
        cout << "error " << ( A * x - b ).lpNorm< Infinity >() << endl;
        result = ( A * x - b ).lpNorm< Infinity >() < precision;
    } else if ( solver_type == "EigenLLT" ) {
        SimplicialLLT< SparseMatrix< double > >simplicial_llt_solver;
        x = simplicial_llt_solver.compute(A).solve(b);
        //result = ( A * x - b ).lpNorm< Infinity >() < precision;
    } else if ( solver_type == "EigenSparseLU" ) {
        SparseLU< SparseMatrix< double >, COLAMDOrdering< int > >sparseLU_solver;
        sparseLU_solver.analyzePattern(A);
        sparseLU_solver.factorize(A);

        x = sparseLU_solver.solve(b);
        //result = ( A * x - b ).lpNorm< Infinity >() < precision;
    }
    cout << "x: " << x << endl;

#if PRINT_DEBUG_TIME
    now = std :: chrono :: system_clock :: now();

    elapsed_seconds = now - start;
    std :: cout << "linalg solver duration: " << convertTimeToString_(elapsed_seconds) << std :: endl;
    cout.flush();
#endif


    return result;
}


bool LinalgNonSymmetricSolver(const CoordinateIndexedSparseMatrix &A, VectorXd &x, const VectorXd &b, const VectorXd x0, double precision, double relmaxit) {
#if PRINT_DEBUG_TIME
    auto start = std :: chrono :: system_clock :: now();
#endif

    //return ConjGrad(A, x, b, x0, precision, relmaxit);

    size_t Maxit = b.size() * relmaxit;

    //BiCGSTAB<SparseMatrix<double>, Eigen::IncompleteLUT<double>> bicg;
    BiCGSTAB< SparseMatrix< double > >bicg;
    bicg.setMaxIterations(Maxit);
    bicg.setTolerance(precision);

    bicg.compute(A);

    x = bicg.solveWithGuess(b, x0);
    //VectorXd cgx = bicg.solve(cgb);

#if PRINT_DEBUG_TIME
    now = std :: chrono :: system_clock :: now();

    elapsed_seconds = now - start;
    std :: cout << "linalg solver duration: " << convertTimeToString(elapsed_seconds) << " nit: " << bicg.iterations() << std :: endl;
    cout.flush();
#endif

    bool result = size_t(bicg.iterations() ) < Maxit;
    return result;
}

//sorterd eigenvalues and eigenvectors
bool LinalgEigenSolver(const VectorXd &A, VectorXd &eigenvalues, vector< VectorXd > &eigenvectors) {

    size_t ndim;
    bool sym;
    if (A.size()==3){ //2D
        ndim = 2; sym = true;
    }else if (A.size()==6){ //3D
        ndim =3; sym = true;
    }else if (A.size()==4){ //2D
        ndim = 2; sym = false;
    }else if (A.size()==9){ //3D
        ndim =3; sym = false;
    }else{
        cerr << "Error: LinalgEigenSolver implemented only for vectorized matrices of size 2 or 3, submitted size " << A.size() << endl;
        exit(1);
    }
    MatrixXd mat = MatrixXd :: Random(ndim, ndim);

    if (ndim==2 && sym){
        mat(0, 0) = A [ 0 ];
        mat(1, 1) = A [ 1 ];
        mat(1, 0) = mat(0, 1) = A [ 2 ];
    }else if (ndim==2 && !sym){
        mat(0, 0) = A [ 0 ];
        mat(1, 1) = A [ 1 ];
        mat(1, 0) = mat(0, 1) = (A [ 2 ] + A [ 3 ])/2.;
    }else if (ndim==3 && sym){
        mat(0, 0) = A [ 0 ];
        mat(1, 1) = A [ 1 ];
        mat(2, 2) = A [ 2 ];
        mat(2, 1) = mat(1, 2) = A [ 3 ];
        mat(2, 0) = mat(0, 2) = A [ 4 ];
        mat(1, 0) = mat(0, 1) = A [ 5 ];
    }else if (ndim==3 && !sym){
        mat(0, 0) = A [ 0 ];
        mat(1, 1) = A [ 1 ];
        mat(2, 2) = A [ 2 ];
        mat(2, 1) = mat(1, 2) = ( A [ 4 ] + A [ 3 ] ) / 2.;
        mat(2, 0) = mat(0, 2) = ( A [ 6 ] + A [ 5 ] ) / 2.;
        mat(1, 0) = mat(0, 1) = ( A [ 8 ] + A [ 7 ] ) / 2.;
    }

   return LinalgEigenSolver(mat, eigenvalues, eigenvectors);
}

bool LinalgEigenSolver(const MatrixXd &mat, VectorXd &eigenvalues, vector< VectorXd > &eigenvectors) {
    EigenSolver< MatrixXd >es(mat);

    unsigned ndim = mat.rows();
    vector< double > eigenvalsvector(ndim);
    eigenvalues.resize(ndim);
    eigenvectors.resize(ndim);

    for ( unsigned i = 0; i < ndim; i++ ) {
        eigenvalsvector [ i ] = real(es.eigenvalues() [ i ]);
    }

    // initialize original index locations
    vector< size_t >idx( ndim );
    iota(idx.begin(), idx.end(), 0);
    stable_sort(idx.begin(), idx.end(), [ & eigenvalsvector ](size_t i1, size_t i2) {
        return eigenvalsvector [ i1 ] > eigenvalsvector [ i2 ];
    });

    for ( unsigned i = 0; i < ndim; i++ ) {
        eigenvalues [ i ] = eigenvalsvector [ idx [ i ] ];
        eigenvectors [ i ].resize(ndim);
        VectorXcd v = es.eigenvectors().col(idx [ i ]);
        for ( unsigned j = 0; j < ndim; j++ ) {
            eigenvectors[i] [j] = real(v [ j ]);
        }
    }

    return true;
}



//JM: Coplanarity check of 4 points (for 3d faces)
double checkCoplanarity(const Point &ptA, const Point &ptB, const Point &ptC, const Point &ptD) {
    Point AB = ptB - ptA;
    Point AC = ptC - ptA;
    Point AD = ptD - ptA;
    //triple scalar product AB*(ACxAD) =>0
    double coplanarityError = AB.dot(AC.cross(AD) );
    return coplanarityError;
}


MatrixXd dyadicProduct(const VectorXd &a, const VectorXd &b) {
    return a * b.transpose();
}

double triArea2D(const Point *a, const Point *b, const Point *c) { //points in counter clockwise direction
    return 0.5 * ( a->x() * ( b->y() - c->y() ) + b->x() * ( c->y() - a->y() ) + c->x() * ( a->y() - b->y() ) );
}

double triArea3D(const Point *a, const Point *b, const Point *c) { //points
    return abs(0.5 * pow(pow( ( b->y() - a->y() ) * ( c->z() - a->z() ) - ( b->z() - a->z() ) * ( c->y() - a->y() ), 2 ) + pow( ( b->z() - a->z() ) * ( c->x() - a->x() ) - ( b->x() - a->x() ) * ( c->z() - a->z() ), 2 ) + pow( ( b->x() - a->x() ) * ( c->y() - a->y() ) - ( b->y() - a->y() ) * ( c->x() - a->x() ), 2 ), 0.5) );
}
