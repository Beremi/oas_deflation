#include "linalg.h"

using namespace std;

bool LinalgSymmetricSolver(const CoordinateIndexedSparseMatrix &A, Vector &x, const Vector &b, const Vector &x0, double precision, double relmaxit, string solver_type) {
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

        Eigen :: ConjugateGradient< Eigen :: SparseMatrix< double >, Eigen :: Lower | Eigen :: Upper >cgK;
        //ConjugateGradient< SparseMatrix< double >, Lower | Upper, IncompleteCholesky< double > >cgK;
        cgK.setMaxIterations(Maxit);
        cgK.setTolerance(precision);

        cgK.compute(A);

        x = cgK.solveWithGuess(b, x0);

        result = size_t( cgK.iterations() ) < Maxit;
        if ( !result ) {
            cout << "Eigen Conjugate Gradients performed " << cgK.iterations() << " iterations and reached error " << cgK.error() << ", required precision is " << precision << endl;
        }
    } else if ( solver_type == "EigenLDLT" ) {
        Eigen :: SimplicialLDLT< Eigen :: SparseMatrix< double > >simplicial_ldlt_solver;
        x = simplicial_ldlt_solver.compute(A).solve(b);
        //cout << "error " << ( A * x - b ).lpNorm< Eigen :: Infinity >() << endl;
        result = ( A * x - b ).lpNorm< Eigen :: Infinity >() < precision;
    } else if ( solver_type == "EigenLLT" ) {
        Eigen :: SimplicialLLT< Eigen :: SparseMatrix< double > >simplicial_llt_solver;
        x = simplicial_llt_solver.compute(A).solve(b);
        result = ( A * x - b ).lpNorm< Eigen :: Infinity >() < precision;
    } else if ( solver_type == "EigenSparseLU" ) {
        Eigen :: SparseLU< Eigen :: SparseMatrix< double >, Eigen :: COLAMDOrdering< int > >sparseLU_solver;
        sparseLU_solver.analyzePattern(A);
        sparseLU_solver.factorize(A);

        x = sparseLU_solver.solve(b);
        result = ( A * x - b ).lpNorm< Eigen :: Infinity >() < precision;
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

    size_t Maxit = b.size() * relmaxit;

    //BiCGSTAB<SparseMatrix<double>, Eigen::IncompleteLUT<double>> bicg;
    Eigen :: BiCGSTAB< Eigen :: SparseMatrix< double > >bicg;
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

    bool result = size_t( bicg.iterations() ) < Maxit;
    return result;
}

//sorterd eigenvalues and eigenvectors
bool LinalgEigenSolver(const Vector &A, Vector &eigenvalues, vector< Vector > &eigenvectors) {
    size_t ndim;
    bool sym;
    if ( A.size() == 3 ) { //2D
        ndim = 2;
        sym = true;
    } else if ( A.size() == 6 ) {      //3D
        ndim = 3;
        sym = true;
    } else if ( A.size() == 4 ) {      //2D
        ndim = 2;
        sym = false;
    } else if ( A.size() == 9 ) {      //3D
        ndim = 3;
        sym = false;
    } else {
        cerr << "Error: LinalgEigenSolver implemented only for vectorized matrices of size 2 or 3, submitted size " << A.size() << endl;
        exit(1);
    }
    Eigen :: MatrixXd mat = Eigen :: MatrixXd :: Zero(ndim, ndim);

    if ( ndim == 2 && sym ) {
        mat(0, 0) = A [ 0 ];
        mat(1, 1) = A [ 1 ];
        mat(1, 0) = mat(0, 1) = A [ 2 ];
    } else if ( ndim == 2 && !sym ) {
        mat(0, 0) = A [ 0 ];
        mat(1, 1) = A [ 1 ];
        mat(1, 0) = mat(0, 1) = ( A [ 2 ] + A [ 3 ] ) / 2.;
    } else if ( ndim == 3 && sym ) {
        mat(0, 0) = A [ 0 ];
        mat(1, 1) = A [ 1 ];
        mat(2, 2) = A [ 2 ];
        mat(2, 1) = mat(1, 2) = A [ 3 ];
        mat(2, 0) = mat(0, 2) = A [ 4 ];
        mat(1, 0) = mat(0, 1) = A [ 5 ];
    } else if ( ndim == 3 && !sym ) {
        mat(0, 0) = A [ 0 ];
        mat(1, 1) = A [ 1 ];
        mat(2, 2) = A [ 2 ];
        mat(2, 1) = mat(1, 2) = ( A [ 4 ] + A [ 3 ] ) / 2.;
        mat(2, 0) = mat(0, 2) = ( A [ 6 ] + A [ 5 ] ) / 2.;
        mat(1, 0) = mat(0, 1) = ( A [ 8 ] + A [ 7 ] ) / 2.;
    }

    return LinalgEigenSolver(mat, eigenvalues, eigenvectors);
}

bool LinalgEigenSolver(const Matrix &mat, Vector &eigenvalues, vector< Vector > &eigenvectors) {
    Eigen :: EigenSolver< Matrix >es;
    es.compute(mat, /* computeEigenvectors = */ true);

    unsigned ndim = mat.rows();
    vector< double >eigenvalsvector(ndim);
    eigenvalues.resize(ndim);
    eigenvectors.resize(ndim);

    for ( unsigned i = 0; i < ndim; i++ ) {
        eigenvalsvector [ i ] = ( es.eigenvalues() [ i ] ).real();
    }

    // initialize original index locations
    vector< size_t >idx(ndim);
    iota(idx.begin(), idx.end(), 0);
    stable_sort(idx.begin(), idx.end(), [ & eigenvalsvector ](size_t i1, size_t i2) {
        return eigenvalsvector [ i1 ] > eigenvalsvector [ i2 ];
    });

    for ( unsigned i = 0; i < ndim; i++ ) {
        eigenvalues [ i ] = eigenvalsvector [ idx [ i ] ];
        eigenvectors [ i ].resize(ndim);
        Eigen :: VectorXcd v = es.eigenvectors().col(idx [ i ]);
        for ( unsigned j = 0; j < ndim; j++ ) {
            eigenvectors [ i ] [ j ] = ( v [ j ] ).real();
        }
    }

    return true;
}

bool LinalgLUSolver(const CoordinateIndexedSparseMatrix &A, Vector &x, const Vector &b) {
    (void) A; (void) x; (void) b;
    //return LinalgSymmetricSolver(A, x, b, x, 1e-12, 2, "EigenConj");


    //Eigen::SparseLU<Eigen::SparseMatrix<double>, Eigen::COLAMDOrdering<int> > solver;

    // fill A and b;
    // Compute the ordering permutation vector from the structural pattern of A
    //solver.analyzePattern(A);

    // Compute the numerical factorization
    //solver.factorize(A);

    //Use the factors to solve the linear system
    //x = solver.solve(b);

    return 0;
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


Matrix dyadicProduct(const Vector &a, const Vector &b) {
    return a * b.transpose();
}

double triArea2D(const Point *a, const Point *b, const Point *c) { //points in counter clockwise direction
    return 0.5 * ( a->x() * ( b->y() - c->y() ) + b->x() * ( c->y() - a->y() ) + c->x() * ( a->y() - b->y() ) );
}

double triArea3D(const Point *a, const Point *b, const Point *c) { //points
    Point AB = ( * b ) - ( * a );
    Point AC = ( * c ) - ( * a );
    return AB.cross(AC).norm() * 0.5;

    //return abs(0.5 * pow(pow( ( b->y() - a->y() ) * ( c->z() - a->z() ) - ( b->z() - a->z() ) * ( c->y() - a->y() ), 2 ) + pow( ( b->z() - a->z() ) * ( c->x() - a->x() ) - ( b->x() - a->x() ) * ( c->z() - a->z() ), 2 ) + pow( ( b->x() - a->x() ) * ( c->y() - a->y() ) - ( b->y() - a->y() ) * ( c->x() - a->x() ), 2 ), 0.5) );
}

double triInertia2D(const Point *a, const Point *b, const Point *c) { // Inertia of triangle relative to global [0, 0]
    double Ix;
    double Iy;
    if ( a->y() == b->y() && a->y() == c->y() ) {
        Ix = 0;
    } else if ( a->x() == b->x() && a->x() == c->x() ) {
        Ix = 0;
    } else if ( b->y() == c->y() ) {
        Ix = ( ( b->x() - a->x() ) / ( b->y() - a->y() ) * ( pow(b->y(), 4) - pow(a->y(), 4) ) / 4. ) -   ( a->y() * ( b->x() - a->x() ) / ( b->y() - a->y() ) * ( pow(b->y(), 3) - pow(a->y(), 3) ) / 3. ) - ( ( c->x() - a->x() ) / ( c->y() - a->y() ) * ( pow(b->y(), 4) - pow(a->y(), 4) ) / 4. ) + ( a->y() * ( c->x() - a->x() ) / ( c->y() - a->y() ) * ( pow(b->y(), 3) - pow(a->y(), 3) ) / 3. );
    } else if ( a->y() == b->y() ) {
        Ix = ( ( b->x() - c->x() ) / ( b->y() - c->y() ) * ( pow(c->y(), 4) - pow(a->y(), 4) ) / 4. ) - ( c->y() * ( b->x() - c->x() ) / ( b->y() - c->y() ) * ( pow(c->y(), 3) - pow(a->y(), 3) ) / 3. ) + ( c->x() * ( pow(c->y(), 3) - pow(a->y(), 3) ) / 3. ) - ( ( c->x() - a->x() ) / ( c->y() - a->y() ) * ( pow(c->y(), 4) - pow(a->y(), 4) ) / 4. ) + ( a->y() * ( c->x() - a->x() ) / ( c->y() - a->y() ) * ( pow(c->y(), 3) - pow(a->y(), 3) ) / 3. ) - ( a->x() * ( pow(c->y(), 3) - pow(a->y(), 3) ) / 3. );
    } else if ( a->y() == c->y() ) {
        Ix = ( ( c->x() - b->x() ) / ( c->y() - b->y() ) * ( pow(b->y(), 4) - pow(a->y(), 4) ) / 4. ) - ( b->y() * ( c->x() - b->x() ) / ( c->y() - b->y() ) * ( pow(b->y(), 3) - pow(a->y(), 3) ) / 3. ) + ( b->x() * ( pow(b->y(), 3) - pow(a->y(), 3) ) / 3. ) - ( ( b->x() - a->x() ) / ( b->y() - a->y() ) * ( pow(b->y(), 4) - pow(a->y(), 4) ) / 4. ) + ( a->y() * ( b->x() - a->x() ) / ( b->y() - a->y() ) * ( pow(b->y(), 3) - pow(a->y(), 3) ) / 3. ) - ( a->x() * ( pow(b->y(), 3) - pow(a->y(), 3) ) / 3. );
    } else {
        Ix = ( ( b->x() - c->x() ) / ( b->y() - c->y() ) * ( pow(c->y(), 4) - pow(a->y(), 4) ) / 4. ) - ( c->y() * ( b->x() - c->x() ) / ( b->y() - c->y() ) * ( pow(c->y(), 3) - pow(a->y(), 3) ) / 3. ) + ( c->x() * ( pow(c->y(), 3) - pow(a->y(), 3) ) / 3. ) - ( ( c->x() - a->x() ) / ( c->y() - a->y() ) * ( pow(c->y(), 4) - pow(a->y(), 4) ) / 4. ) + ( a->y() * ( c->x() - a->x() ) / ( c->y() - a->y() ) * ( pow(c->y(), 3) - pow(a->y(), 3) ) / 3. ) - ( a->x() * ( pow(c->y(), 3) - pow(a->y(), 3) ) / 3. ) - ( ( b->x() - c->x() ) / ( b->y() - c->y() ) * ( pow(b->y(), 4) - pow(a->y(), 4) ) / 4. ) + ( c->y() * ( b->x() - c->x() ) / ( b->y() - c->y() ) * ( pow(b->y(), 3) - pow(a->y(), 3) ) / 3. ) - ( c->x() * ( pow(b->y(), 3) - pow(a->y(), 3) ) / 3. ) + ( ( b->x() - a->x() ) / ( b->y() - a->y() ) * ( pow(b->y(), 4) - pow(a->y(), 4) ) / 4. ) - ( a->y() * ( b->x() - a->x() ) / ( b->y() - a->y() ) * ( pow(b->y(), 3) - pow(a->y(), 3) ) / 3. ) + ( a->x() * ( pow(b->y(), 3) - pow(a->y(), 3) ) / 3. );
    }

    if ( a->y() == b->y() && a->y() == c->y() ) {
        Iy = 0;
    } else if ( a->x() == b->x() && a->x() == c->x() ) {
        Iy = 0;
    } else if ( b->x() == c->x() ) {
        Iy = ( ( b->y() - a->y() ) / ( b->x() - a->x() ) * ( pow(b->x(), 4) - pow(a->x(), 4) ) / 4. ) - ( a->x() * ( b->y() - a->y() ) / ( b->x() - a->x() ) * ( pow(b->x(), 3) - pow(a->x(), 3) ) / 3. ) - ( ( c->y() - a->y() ) / ( c->x() - a->x() ) * ( pow(b->x(), 4) - pow(a->x(), 4) ) / 4. ) + ( a->x() * ( c->y() - a->y() ) / ( c->x() - a->x() ) * ( pow(b->x(), 3) - pow(a->x(), 3) ) / 3. );
    } else if ( a->x() == b->x() ) {
        Iy = ( ( b->y() - c->y() ) / ( b->x() - c->x() ) * ( pow(c->x(), 4) - pow(a->x(), 4) ) / 4. ) - ( c->x() * ( b->y() - c->y() ) / ( b->x() - c->x() ) * ( pow(c->x(), 3) - pow(a->x(), 3) ) / 3. ) + ( c->y() * ( pow(c->x(), 3) - pow(a->x(), 3) ) / 3. ) - ( ( c->y() - a->y() ) / ( c->x() - a->x() ) * ( pow(c->x(), 4) - pow(a->x(), 4) ) / 4. ) + ( a->x() * ( c->y() - a->y() ) / ( c->x() - a->x() ) * ( pow(c->x(), 3) - pow(a->x(), 3) ) / 3. ) - ( a->y() * ( pow(c->x(), 3) - pow(a->x(), 3) ) / 3. );
    } else if ( a->x() == c->x() ) {
        Iy = ( ( c->y() - b->y() ) / ( c->x() - b->x() ) * ( pow(b->x(), 4) - pow(a->x(), 4) ) / 4. ) - ( b->x() * ( c->y() - b->y() ) / ( c->x() - b->x() ) * ( pow(b->x(), 3) - pow(a->x(), 3) ) / 3. ) + ( b->y() * ( pow(b->x(), 3) - pow(a->x(), 3) ) / 3. ) - ( ( b->y() - a->y() ) / ( b->x() - a->x() ) * ( pow(b->x(), 4) - pow(a->x(), 4) ) / 4. ) + ( a->x() * ( b->y() - a->y() ) / ( b->x() - a->x() ) * ( pow(b->x(), 3) - pow(a->x(), 3) ) / 3. ) - ( a->y() * ( pow(b->x(), 3) - pow(a->x(), 3) ) / 3. );
    } else {
        Iy = ( ( b->y() - c->y() ) / ( b->x() - c->x() ) * ( pow(c->x(), 4) - pow(a->x(), 4) ) / 4. ) - ( c->x() * ( b->y() - c->y() ) / ( b->x() - c->x() ) * ( pow(c->x(), 3) - pow(a->x(), 3) ) / 3. ) + ( c->y() * ( pow(c->x(), 3) - pow(a->x(), 3) ) / 3. ) - ( ( c->y() - a->y() ) / ( c->x() - a->x() ) * ( pow(c->x(), 4) - pow(a->x(), 4) ) / 4. ) + ( a->x() * ( c->y() - a->y() ) / ( c->x() - a->x() ) * ( pow(c->x(), 3) - pow(a->x(), 3) ) / 3. ) - ( a->y() * ( pow(c->x(), 3) - pow(a->x(), 3) ) / 3. ) - ( ( b->y() - c->y() ) / ( b->x() - c->x() ) * ( pow(b->x(), 4) - pow(a->x(), 4) ) / 4. ) + ( c->x() * ( b->y() - c->y() ) / ( b->x() - c->x() ) * ( pow(b->x(), 3) - pow(a->x(), 3) ) / 3. ) - ( c->y() * ( pow(b->x(), 3) - pow(a->x(), 3) ) / 3. ) + ( ( b->y() - a->y() ) / ( b->x() - a->x() ) * ( pow(b->x(), 4) - pow(a->x(), 4) ) / 4. ) - ( a->x() * ( b->y() - a->y() ) / ( b->x() - a->x() ) * ( pow(b->x(), 3) - pow(a->x(), 3) ) / 3. ) + ( a->y() * ( pow(b->x(), 3) - pow(a->x(), 3) ) / 3. );
    }

    double Iz = abs(Ix) + abs(Iy);
    return Iz;
}

Matrix tetraInertia3D(const Point *a, const Point *b, const Point *c, const Point *d) { // Inertia relative to the centroid of the tetrahedron, where the centroid must be [0,0,0]
    // Jacobian of transformation to the reference tetrahedron
    double detJ = ( b->x() - a->x() ) * ( c->y() - a->y() ) * ( d->z() - a->z() ) +  ( c->x() - a->x() ) * ( d->y() - a->y() ) * ( b->z() - a->z() ) + ( d->x() - a->x() ) * ( b->y() - a->y() ) * ( c->z() - a->z() ) - ( b->x() - a->x() ) * ( d->y() - a->y() ) * ( c->z() - a->z() ) - ( c->x() - a->x() ) * ( b->y() - a->y() ) * ( d->z() - a->z() ) - ( d->x() - a->x() ) * ( c->y() - a->y() ) * ( b->z() - a->z() );
    // Inertia
    double Ixx = abs(detJ) * ( pow(a->y(), 2) + a->y() * b->y() + pow(b->y(), 2) + a->y() * c->y() + b->y() * c->y() + pow(c->y(), 2) + a->y() * d->y() + b->y() * d->y() + c->y() * d->y() + pow(d->y(), 2) + pow(a->z(), 2) + a->z() * b->z() + pow(b->z(), 2) + a->z() * c->z() + b->z() * c->z() + pow(c->z(), 2) + a->z() * d->z() + b->z() * d->z() + c->z() * d->z() + pow(d->z(), 2) ) / 60.;

    double Iyy = abs(detJ) * ( pow(a->x(), 2) + a->x() * b->x() + pow(b->x(), 2) + a->x() * c->x() + b->x() * c->x() + pow(c->x(), 2) + a->x() * d->x() + b->x() * d->x() + c->x() * d->x() + pow(d->x(), 2) + pow(a->z(), 2) + a->z() * b->z() + pow(b->z(), 2) + a->z() * c->z() + b->z() * c->z() + pow(c->z(), 2) + a->z() * d->z() + b->z() * d->z() + c->z() * d->z() + pow(d->z(), 2) ) / 60.;

    double Izz = abs(detJ) * ( pow(a->x(), 2) + a->x() * b->x() + pow(b->x(), 2) + a->x() * c->x() + b->x() * c->x() + pow(c->x(), 2) + a->x() * d->x() + b->x() * d->x() + c->x() * d->x() + pow(d->x(), 2) + pow(a->y(), 2) + a->y() * b->y() + pow(b->y(), 2) + a->y() * c->y() + b->y() * c->y() + pow(c->y(), 2) + a->y() * d->y() + b->y() * d->y() + c->y() * d->y() + pow(d->y(), 2) ) / 60.;

    double Ixz = abs(detJ) * ( 2. * a->x() * a->z() + b->x() * a->z() + c->x() * a->z() + d->x() * a->z() + a->x() * b->z() + 2. * b->x() * b->z() + c->x() * b->z() + d->x() * b->z() + a->x() * c->z() * b->x() * c->z() + 2. * c->x() * c->z() + d->x() * c->z() + a->x() * d->z() + b->x() * d->z() + c->x() * d->z() + 2. * d->x() * d->z() ) / 120.;

    double Iyz = abs(detJ) * ( 2. * a->y() * a->z() + b->y() * a->z() + c->y() * a->z() + d->y() * a->z() + a->y() * b->z() + 2. * b->y() * b->z() + c->y() * b->z() + d->y() * b->z() + a->y() * c->z() * b->y() * c->z() + 2. * c->y() * c->z() + d->y() * c->z() + a->y() * d->z() + b->y() * d->z() + c->y() * d->z() + 2. * d->y() * d->z() ) / 120.;

    double Ixy = abs(detJ) * ( 2. * a->x() * a->y() + b->x() * a->y() + c->x() * a->y() + d->x() * a->y() + a->x() * b->y() + 2. * b->x() * b->y() + c->x() * b->y() + d->x() * b->y() + a->x() * c->y() * b->x() * c->y() + 2. * c->x() * c->y() + d->x() * c->y() + a->x() * d->y() + b->x() * d->y() + c->x() * d->y() + 2. * d->x() * d->y() ) / 120.;

    // Inertia matrix
    Matrix I = Matrix :: Zero(3, 3);

    I(0, 0) = Ixx;
    I(1, 1) = Iyy;
    I(2, 2) = Izz;
    I(0, 1) = I(1, 0) = -Ixy;
    I(0, 2) = I(2, 0) = -Ixz;
    I(1, 2) = I(2, 1) = -Iyz;

    return I;
}

double tetraVolumeSigned(const Point *a, const Point *b, const Point *c, const Point *d) {
    return ( ( * d ) - ( * a ) ).dot( ( ( * b ) - ( * d ) ).cross( ( * c ) - ( * d ) ) ) / 6.;
}

bool is_positive_integer(const std :: string &s)
{
    return !s.empty() && std :: find_if(s.begin(),
                                        s.end(), [](unsigned char c) {
        return !std :: isdigit(c);
    }) == s.end();
}

void giveGaussIntegrationPointAndWeights(unsigned n, Vector &locs, Vector &weis) {
    locs.resize(n);
    weis.resize(n);
    if ( n == 1 ) {
        locs [ 0 ] = 0;
        weis [ 0 ] = 2;
    } else if ( n == 2 ) {
        locs [ 0 ] = -1. / sqrt(3);
        locs [ 1 ] = -locs [ 0 ];
        weis [ 0 ] = weis [ 1 ] = 1;
    } else if ( n == 3 ) {
        locs [ 0 ] = -sqrt(3. / 5.);
        locs [ 1 ] = 0;
        locs [ 2 ] = -locs [ 0 ];
        weis [ 0 ] = weis [ 2 ] = 5. / 9.;
        weis [ 1 ] = 8. / 9.;
    } else if ( n == 4 ) {
        locs [ 0 ] = -sqrt( 3. / 7. + 2. / 7. * sqrt(6. / 5.) );
        locs [ 1 ] = -sqrt( 3. / 7. - 2. / 7. * sqrt(6. / 5.) );
        locs [ 2 ] = -locs [ 1 ];
        locs [ 3 ] = -locs [ 0 ];
        weis [ 0 ] = weis [ 3 ] = ( 18. - sqrt(30.) ) / 36;
        weis [ 1 ] = weis [ 2 ] = ( 18. + sqrt(30.) ) / 36;
    } else if ( n == 5 ) {
        locs [ 0 ] = -sqrt( 5. + 2. * sqrt(10. / 7.) ) / 3;
        locs [ 1 ] = -sqrt( 5. - 2. * sqrt(10. / 7.) ) / 3;
        locs [ 2 ] = 0.;
        locs [ 3 ] = -locs [ 1 ];
        locs [ 4 ] = -locs [ 0 ];
        weis [ 0 ] = weis [ 4 ] = ( 322. - 13. * sqrt(70.) ) / 900;
        weis [ 1 ] = weis [ 3 ] = ( 322. + 13. * sqrt(70.) ) / 900;
        weis [ 2 ] = 128. / 225.;
    } else {
        cerr << "Gauss integration for n=" << n << " not implemented" << endl;
        exit(1);
    }
}
