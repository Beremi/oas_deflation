
#ifndef LINEAR_ALGEBRA_H
#define LINEAR_ALGEBRA_H

#include <valarray>
#include <map>
#include <complex>
#include <algorithm>
#include <assert.h>
#include <numeric>
#include <iostream>
#include <vector>

using namespace std;

#define POINT_TOLERANCE 1e-8

typedef valarray< double >Vector;


class Point
{
public:
    double x;
    double y;
    double z;
    Point();
    Point(double x);
    Point(double x, double y);
    Point(double x, double y, double z);
    Point(const Point &p);
    ~Point(void);

    void setX(double x);
    void setY(double y);
    void setZ(double z);
    double getX() const;
    double getY() const;
    double getZ() const;
    double giveCoord(unsigned i) const { if ( i == 0 ) { return x; } else if ( i == 1 ) { return y; } else if ( i == 2 ) { return z; } else { cerr << "Point Error: attempt to read coordinate " <<  i << endl; exit(0); } };
    void setCoord(unsigned i, double val) { if ( i == 0 ) { x = val; } else if ( i == 1 ) { y = val; } else if ( i == 2 ) { z = val; } else { cerr << "Point Error: attempt to read coordinate " <<  i << endl; exit(0); } };
    void set(double x, double y, double z);
    void set(const Point &p);
    void set(const Point *p);
    Point &operator=(const Point &A);

    bool operator==(const Point &p) const;
    bool operator!=(const Point &p) const;

    Point operator-(const Point &p) const;
    Point operator+(const Point &p) const;
    Point operator/(const double p) const;
    Point operator*(const double p) const;
    double operator*(const Point &p) const;
    Point normalized() const;
    void normalize();

    void operator+=(const Point &p);
    void operator-=(const Point &p);

    void operator*=(const double d) {
        x *= d;
        y *= d;
        z *= d;
    }

    void operator/=(const double d) {
        x /= d;
        y /= d;
        z /= d;
    }
    double norm() const;
    double sqNorm() const;
    double sum() const;
    void print() const;
};

//Point cross(const Point &p, const Point &q);
Point cross(const Point p, const Point q);
double dot(const Point &p, const Point &q);
double squareDist(const Point &v1, const Point &v2);
double squareDist(const Point *v1, const Point *v2);
double determinant(const Point &u, const Point &v, const Point &w);

//JM coplanarity
double checkCoplanarity(const Point &ptA, const Point &ptB, const Point &ptC, const Point &ptD);

template< class T >
class Slice_iter
{
    valarray< T > *v;
    slice s;
    size_t curr;

    T &ref(size_t i) const { return ( * v ) [ s.start() + i * s.stride() ]; }

public:

    Slice_iter(valarray< T > *vv, slice ss) : v(vv), s(ss), curr(0) { }

    Slice_iter end() const
    {
        Slice_iter t = * this;
        t.curr = s.size();
        return t;
    }

    Slice_iter &operator++() { curr++; return * this; }
    Slice_iter operator++(int) { Slice_iter t = * this; curr++; return t; }

    T &operator[](size_t i) { return ref(i); }
    T &operator[](size_t i) const { return ref(i); }
    T &operator()(size_t i) { return ref(i); }
    T &operator()(size_t i) const { return ref(i); }
    T &operator*() { return ref(curr); }

    template< class C >
    friend bool operator==(const Slice_iter &p, const Slice_iter &q);
    template< class C >
    friend bool operator!=(const Slice_iter &p, const Slice_iter &q);
    template< class C >
    friend bool operator<(const Slice_iter &p, const Slice_iter &q);
};


template< class T >bool operator==(const Slice_iter< T > &p, const Slice_iter< T > &q)
{
    return p.curr == q.curr && p.s.stride() == q.s.stride() && p.s.start() == q.s.start();
};

template< class T >bool operator!=(const Slice_iter< T > &p, const Slice_iter< T > &q)
{
    return !( p == q );
};

template< class T >bool operator<(const Slice_iter< T > &p, const Slice_iter< T > &q)
{
    return p.curr < q.curr && p.s.stride() == q.s.stride() && p.s.start() == q.s.start();
};

template< class T >
class Cslice_iter
{
public:
    valarray< T > *v;
    slice s;
    size_t curr;

    T &ref(size_t i) const { return ( * v ) [ s.start() + i * s.stride() ]; }

public:

    Cslice_iter(valarray< T > *vv, slice ss) : v(vv), s(ss), curr(0) { }

    Cslice_iter end() const
    {
        Cslice_iter t = * this;
        t.curr = s.size();
        return t;
    }

    Cslice_iter &operator++() { curr++; return * this; }
    Cslice_iter operator++(int) { Cslice_iter t = * this; curr++; return t; }

    T &operator[](size_t i) const { return ref(i); }
    T &operator()(size_t i) const { return ref(i); }
    T &operator*() const { return ref(curr); }

    template< class C >
    friend bool operator==(const Cslice_iter &p, const Cslice_iter &q);
    template< class C >
    friend bool operator!=(const Cslice_iter &p, const Cslice_iter &q);
    template< class C >
    friend bool operator<(const Cslice_iter &p, const Cslice_iter &q);
};

template< class T >bool operator==(const Cslice_iter< T > &p, const Cslice_iter< T > &q)
{
    return p.curr == q.curr && p.s.stride() == q.s.stride() && p.s.start() == q.s.start();
}

template< class T >bool operator!=(const Cslice_iter< T > &p, const Cslice_iter< T > &q)
{
    return !( p == q );
}

template< class T >bool operator<(const Cslice_iter< T > &p, const Cslice_iter< T > &q)
{
    return p.curr < q.curr && p.s.stride() == q.s.stride() && p.s.start() == q.s.start();
}







struct MtM;

class Matrix
{
    size_t r, c;
    Vector *v;
public:



    /** Construct an M x N matrix.
     *
     * All Elements are initialised to 0.
     *
     * @param x number of rows
     * @param y number of columns
     */
    Matrix(size_t x, size_t y);

    /** Construct a 2x2 matrix
     */
    Matrix() {
        r = 2;
        c = 2;
        v = new Vector(0., 4);
    }

    /** Construct a square submatrix from a given source matrix.
     *
     * @param rl size of the submatrix
     * @param k staring row in the source matrix
     * @param l staring column in the source matrix
     * @param m source matrix
     */
    Matrix(size_t rl, size_t k, size_t l, const Matrix &m);


    /** construct a matrix from a given source matrix without specified row and column.
     *
     * @param i
     * @param j
     * @param m
     */
    Matrix(size_t i, size_t j, const Matrix &m);


    Matrix(const Matrix &);

    Matrix(const Point &p);

    virtual ~Matrix() {
        delete v;
    }

    Matrix &operator=(const Matrix &m);

    size_t size() const {
        return r * c;
    }

    size_t numCols() const {
        return c;
    }

    size_t numRows() const {
        return r;
    }

    Slice_iter< double >column(size_t i) {
        return Slice_iter< double >(v, slice(i, r, c) );
    }

    Cslice_iter< double >column(size_t i) const {
        return Cslice_iter< double >(v, slice(i, r, c) );
    }

    Slice_iter< double >row(size_t i) {
        return Slice_iter< double >(v, slice(i * c, c, 1) );
    }

    Cslice_iter< double >row(size_t i) const {
        return Cslice_iter< double >(v, slice(i * c, c, 1) );
    }

    Matrix transpose() const;

    double &operator()(size_t x, size_t y);
    double operator()(size_t x, size_t y) const;

    Slice_iter< double >operator()(size_t i) {
        return row(i);
    }

    Cslice_iter< double >operator()(size_t i) const {
        return row(i);
    }

    Slice_iter< double >operator[](size_t i) {
        return row(i);
    }

    Cslice_iter< double >operator[](size_t i) const {
        return row(i);
    }

    Matrix &operator*=(double);
    Matrix operator*(double) const;
    Matrix operator/(double) const;
    Matrix &operator/=(double);
    Matrix &operator*=(const Matrix &m);
    Vector &operator*=(const Vector &m);
    //  const Matrix& operator *(const Matrix &m) const ;
    Matrix &operator+=(const Matrix &m);
    Matrix &operator+=(const double x);
    Matrix operator+(const Matrix &m) const;
    Matrix &operator-=(const Matrix &m);
    Matrix operator-(const Matrix &m) const;
    Matrix &operator=(const MtM &m);

    bool operator==(const Matrix &m);
    bool operator!=(const Matrix &m);

    Vector &array() {
        return * v;
    }

    Vector array() const {
        return * v;
    }

    void print() const;

    double *GetArray() const;
};

struct MtV {
    const Matrix &m;
    const Vector &v;

    MtV(const Matrix &mm, const Vector &vv) : m(mm), v(vv) {}

    operator const Vector();
};

struct MtM {
    const Matrix &first;
    const Matrix &second;

    MtM(const Matrix &mm, const Matrix &mmm) : first(mm), second(mmm) {}

    operator const Matrix() const;
};

Matrix dyadicProduct(const Vector &a, const Vector &b);

Vector solveSystem(const Matrix &A, const Vector &b, Vector &x);

inline MtV operator*(const Matrix &mm, const Vector &v) {
    return MtV(mm, v);
};

inline MtM operator*(const Matrix &mm, const Matrix &mmm) {
    return MtM(mm, mmm);
};

inline const Matrix matrix_multiply(const Matrix &m0, const Matrix &m1) {
    assert(m0.numCols() == m1.numRows() );

    Matrix ret(m0.numRows(), m1.numCols() );

    for ( size_t i = 0; i < m0.numRows(); i++ ) {
        for ( size_t j = 0; j < m1.numCols(); j++ ) {
            const Cslice_iter< double > &ri = m0.row(i);
            const Cslice_iter< double > &cj = m1.column(j);
            ret [ i ] [ j ] = inner_product(& ri [ 0 ], & ri [ m0.numCols() ], cj, ( double ) ( 0 ) );
        }
    }
    return ret;
}

inline const Vector matrix_vector_multiply(const Matrix &m, const Vector &v) {
    assert(m.numCols() == v.size() );

    Vector ret(m.numRows() );

    for ( size_t i = 0; i < m.numRows(); i++ ) {
        const Cslice_iter< double > &ri = m.row(i);
        ret [ i ] = inner_product(& ri [ 0 ], & ri [ m.numCols() ], & v [ 0 ], ( double ) ( 0 ) );
    }
    return ret;
}

inline const Vector operator*(const Vector &v, const Matrix &m) {
    assert(m.numRows() == v.size() );

    Vector ret(v.size() );

    for ( size_t i = 0; i < m.numCols(); i++ ) {
        const Cslice_iter< double > &ri = m.column(i);
        ret [ i ] = inner_product(& ri [ 0 ], & ri [ m.numCols() ], & v [ 0 ], ( double ) ( 0 ) );
    }
    return ret;
}

//////////////////////////////////////////////////////////

struct SparseVector
{
public:
    Vector &val;
    valarray< unsigned > &idx;
    const size_t length;
    const size_t start;

    double zero;

public:
    SparseVector(Vector &v, valarray< unsigned > &idx, const size_t l, const size_t s);

    double operator[](const size_t) const;
    double &operator[](const size_t);
    double operator*(const Vector &) const;
    SparseVector operator*(const double d) const;
    double operator*(const SparseVector &) const;
    Vector operator+(const Vector &) const;

    void print() const;
};

//////////////////////////////////////////////////////////

struct ConstSparseVector
{
public:
    const Vector &val;
    const valarray< unsigned > &idx;
    const size_t length;
    const size_t start;

public:
    ConstSparseVector(const Vector &v,  const valarray< unsigned > &idx, const size_t l, const size_t s);

    inline double operator[](const size_t i) const
    {
        const unsigned *i_index_pointer = find(& idx [ start ], & idx [ min(start + length, idx.size() ) ], i);
        unsigned offset = i_index_pointer - & idx [ start ]; //todo: warning C4244: 'initializing': conversion from '__int64' to 'unsigned int', possible loss of data
        if ( i_index_pointer != & idx [ min(start + length, idx.size() ) ] ) {
            return ( val [ start + offset ] );
        }

        return 0;
    }
    double vectorMultiply(const double *vecA) const;

    double operator*(const Vector &) const;
    double operator*(const SparseVector &) const;
    double operator*(const ConstSparseVector &) const;
    Vector operator+(const Vector &) const;

    void print() const;
};

//////////////////////////////////////////////////////////

double innerProduct(const SparseVector &v0, const SparseVector &v1, const size_t end);
double innerProduct(const SparseVector &v0, const Vector &v1, const size_t end);
double innerProduct(const ConstSparseVector &v0, const Vector &v1, const size_t end);
inline void innerProductAssignAndAdd(const ConstSparseVector &v0, Vector &v1, double &t,  double toAdd, size_t end)
{
    for ( size_t j =  v0.start; v0.idx [ j ] < end; ++j ) {
        t += v1 [ v0.idx [ j ] ] * v0.val [ j ];
    }
    t += toAdd;
};

double reverseInnerProduct(const SparseVector &v0, const Vector &v1, const size_t start);
double reverseInnerProduct(const ConstSparseVector &v0, const Vector &v1, const size_t start);
inline void reverseInnerProductAssignAndAdd(const ConstSparseVector &v0, Vector &v1, double &t,  double toAdd, size_t start)
{
    for ( size_t j = v0.length + v0.start - 1; v0.idx [ j ] > start; --j ) {
        t += v1 [ v0.idx [ j ] ] * v0.val [ j ];
    }
    t += toAdd;
};

//////////////////////////////////////////////////////////

class CoordinateIndexedSparseMatrix
{
public:
    Vector array;
    unsigned RowCount;
    unsigned ColumnCount;
    valarray< unsigned >column_index;
    valarray< unsigned >row_size;
    valarray< unsigned >accumulated_row_size;
public:
    CoordinateIndexedSparseMatrix(map< pair< size_t, size_t >, double > &source, unsigned RowCountI, unsigned ColumnCountI);
    CoordinateIndexedSparseMatrix(map< pair< size_t, size_t >, Matrix > &source, unsigned RowCountI, unsigned ColumnCountI);
    CoordinateIndexedSparseMatrix(const valarray< unsigned > &, const valarray< unsigned > &, unsigned RowCountI, unsigned ColumnCountI);
    CoordinateIndexedSparseMatrix(const CoordinateIndexedSparseMatrix &source);
    CoordinateIndexedSparseMatrix();

    ~CoordinateIndexedSparseMatrix();

    SparseVector operator[](const size_t i);
    const ConstSparseVector operator[](const size_t i) const;
    double operator()(const size_t i, const size_t j) const;


    void vectorMultiply(const double *vecA, double *vecB) const;

    Vector operator*(const Vector &v) const;
    CoordinateIndexedSparseMatrix operator*(const double d) const;
    CoordinateIndexedSparseMatrix operator*(const CoordinateIndexedSparseMatrix &b) const;
    CoordinateIndexedSparseMatrix operator+(const CoordinateIndexedSparseMatrix &b) const;
    //CoordinateIndexedSparseMatrix operator -(const CoordinateIndexedSparseMatrix & b) const ;
    CoordinateIndexedSparseMatrix &operator=(const CoordinateIndexedSparseMatrix &);

    CoordinateIndexedSparseMatrix ExtendColumn(const Vector &c) const;
    CoordinateIndexedSparseMatrix ExtendRow(const Vector &c) const;

    Vector inverseDiagonal() const;
    Vector inverseDiagonalSquared() const;
    Vector diagonal() const;
    CoordinateIndexedSparseMatrix transpose() const;

    double froebeniusNorm() const;
    double infinityNorm() const;

    bool isThereNaN() const;
    void print(unsigned size1, unsigned size2);
    void print();
};

double l2_norm(Vector x);
bool ConjGrad(const CoordinateIndexedSparseMatrix &A, Vector &x, const Vector &b, const Vector x0);
bool ConjGrad(const CoordinateIndexedSparseMatrix &A, Vector &x, const Vector &b, const Vector x0, double precision, double relmaxit);
bool isMatrixSingular(const CoordinateIndexedSparseMatrix &A);


double triArea2D(const Point *a, const Point *b, const Point *c);
Point triNormal3D(const Point *a, const Point *b, const Point *c);
double triArea3D(const Point *a, const Point *b, const Point *c);
double tetVolume3D(const Point *a, const Point *b, const Point *c, const Point *d);


#endif
