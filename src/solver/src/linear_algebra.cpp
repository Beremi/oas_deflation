//////////////////////////////////////////////////////////
//////////////////J A N  E L I A S////////////////////////
///////////////////// 2 0 1 9 ////////////////////////////
//////////////////////////////////////////////////////////

#include "linear_algebra.h"

Point :: Point() {
    x = 0;
    y = 0;
    z = 0;
}

Point :: Point(double ox) {
    x = ox;
    y = 0;
    z = 0;
}

Point :: Point(double ox, double oy) {
    x = ox;
    y = oy;
    z = 0.0;
}

Point :: Point(double ox, double oy, double oz) {
    x = ox;
    y = oy;
    z = oz;
}

Point :: Point(const Point &p) {
    x = p.x;
    y = p.y;
    z = p.z;
}

Point :: ~Point(void) {}

void Point :: setX(double ox) {
    x = ox;
}

void Point :: setY(double oy) {
    y = oy;
}

void Point :: setZ(double oz) {
    z = oz;
}

double Point :: getX() const {
    return x;
}

double Point :: getY() const {
    return y;
}

double Point :: getZ() const {
    return z;
}

void Point :: set(double ox, double oy, double oz) {
    x = ox;
    y = oy;
    z = oz;
}

void Point :: set(const Point &p) {
    x = p.x;
    y = p.y;
    z = p.z;
}

void Point :: set(const Point *p) {
    x = p->x;
    y = p->y;
    z = p->z;
}

bool Point :: operator==(const Point &p) const {
    double delta = POINT_TOLERANCE;
    Point a(p);
    a.x += delta;
    a.y += delta;
    a.z += delta;
    Point b(p);
    b.x += delta;
    b.y += delta;
    b.z -= delta;
    Point c(p);
    c.x += delta;
    c.y -= delta;
    c.z += delta;
    Point d(p);
    d.x += delta;
    d.y -= delta;
    d.z -= delta;
    Point e(p);
    e.x -= delta;
    e.y += delta;
    e.z += delta;
    Point f(p);
    f.x -= delta;
    f.y += delta;
    f.z -= delta;
    Point g(p);
    g.x -= delta;
    g.y -= delta;
    g.z += delta;
    Point h(p);
    h.x -= delta;
    h.y -= delta;
    h.z -= delta;

    return squareDist(& p, this) < POINT_TOLERANCE * POINT_TOLERANCE ||
           squareDist(& a, this) < POINT_TOLERANCE * POINT_TOLERANCE ||
           squareDist(& b, this) < POINT_TOLERANCE * POINT_TOLERANCE ||
           squareDist(& c, this) < POINT_TOLERANCE * POINT_TOLERANCE ||
           squareDist(& d, this) < POINT_TOLERANCE * POINT_TOLERANCE ||
           squareDist(& e, this) < POINT_TOLERANCE * POINT_TOLERANCE ||
           squareDist(& f, this) < POINT_TOLERANCE * POINT_TOLERANCE ||
           squareDist(& g, this) < POINT_TOLERANCE * POINT_TOLERANCE ||
           squareDist(& h, this) < POINT_TOLERANCE * POINT_TOLERANCE;

    return abs(x - p.x) < .5 * POINT_TOLERANCE &&
           abs(y - p.y) < .5 * POINT_TOLERANCE &&
           abs(z - p.z) < .5 * POINT_TOLERANCE
    ;

    return squareDist(& p, this) < POINT_TOLERANCE * POINT_TOLERANCE;
}

bool Point :: operator!=(const Point &p) const {
    return squareDist(this, & p) > POINT_TOLERANCE * POINT_TOLERANCE;
}

Point Point :: operator-(const Point &p) const {
    Point ret( ( * this ) );
    ret.x -= p.x;
    ret.y -= p.y;
    ret.z -= p.z;
    return ret;
}

Point Point :: operator+(const Point &p) const {
    Point ret( ( * this ) );
    ret.x += p.x;
    ret.y += p.y;
    ret.z += p.z;
    return ret;
}

void Point :: operator+=(const Point &p) {
    x += p.x;
    y += p.y;
    z += p.z;
}

void Point :: operator-=(const Point &p) {
    x -= p.x;
    y -= p.y;
    z -= p.z;
}

Point Point :: operator/(const double p) const {
    Point ret( ( * this ) );
    double inv = 1. / p;
    ret.x *= inv;
    ret.y *= inv;
    ret.z *= inv;
    return ret;
}

Point Point :: operator*(const double p) const {
    Point ret( ( * this ) );
    ret.x *= p;
    ret.y *= p;
    ret.z *= p;
    return ret;
}

double Point :: operator*(const Point &p) const {
    return x * p.x + y * p.y + z * p.z;
}

Point Point :: normalized() const {
    Point ret( ( * this ) );
    ret.normalize();
    return ret;
}

void Point :: normalize() {
    double norm = this->norm();
    x = x / norm;
    y = y / norm;
    z = z / norm;
}

double Point :: norm() const {
    return sqrt(x * x + y * y + z * z);
}

double Point :: sqNorm() const {
    return x * x + y * y + z * z;
}

double Point :: sum() const {
    return x + y + z;
}

void Point :: print() const {
    cout << getX() << flush;
    cout << "; " << getY() << flush;
    cout << "; " << getZ() << flush << endl;
}

Point cross(const Point p, const Point q) {
    Point r( ( p.y *q.z - p.z *q.y ), ( p.z *q.x - p.x *q.z ), ( p.x *q.y - p.y *q.x ) );
    return r;
}

double dot(const Point &p, const Point &q) {
    double r(p.x *q.x + p.y *q.y + p.z *q.z);
    return r;
}

double squareDist(const Point &v1, const Point &v2) {
    double x = v2.getX() - v1.getX();
    double y = v2.getY() - v1.getY();
    double z = v2.getZ() - v1.getZ();

    return x * x + y * y + z * z;
}

double squareDist(const Point *v1, const Point *v2) {
    double x = v2->getX() - v1->getX();
    double y = v2->getY() - v1->getY();
    double z = v2->getZ() - v1->getZ();
    return x * x + y * y + z * z;
}

double determinant(const Point &u, const Point &v, const Point &w) {
    double D = 0.0;
    D = u.x * v.y * w.z - u.x * v.z * w.y + u.y * v.z * w.x - u.y * v.x * w.z + u.z * v.x * w.y - u.z * v.y * w.x;
    return D;
}

//JM: Coplanarity check of 4 points (for 3d faces)
double checkCoplanarity(const Point &ptA, const Point &ptB, const Point &ptC, const Point &ptD) {
    Point AB = ptB - ptA;
    Point AC = ptC - ptA;
    Point AD = ptD - ptA;
    //triple scalar product AB*(ACxAD) =>0
    double coplanarityError = dot(AB, cross(AC, AD) );
    return coplanarityError;
}



Matrix :: Matrix(size_t x, size_t y) {
    r = x;
    c = y;
    v = new Vector(0., x * y);
}

Matrix :: Matrix(size_t rl, size_t k, size_t l, const Matrix &m) {
    r = rl;
    c = rl;
    v = new Vector(rl * rl);
    for ( size_t i = 0; i < rl; i++ ) {
        for ( size_t j = 0; j < rl; j++ ) {
            ( * v ) [ i * rl + j ] = m [ k + i ] [ l + j ];
        }
    }
}

Matrix :: Matrix(size_t i, size_t j, const Matrix &m) {
    r = m.numRows() - 1;
    c = m.numCols() - 1;
    v = new Vector(r * c);
    for ( size_t l = 0; l < m.numRows(); l++ ) {
        if ( l != i ) {
            for ( size_t k = 0; k < m.numCols(); k++ ) {
                if ( k != j ) {
                    size_t r_index = l;
                    if ( l > i ) {
                        r_index--;
                    }
                    size_t c_index = k;
                    if ( k > j ) {
                        c_index--;
                    }

                    ( * this ) [ r_index ] [ c_index ] = m [ l ] [ k ];
                }
            }
        }
    }
}

Matrix :: Matrix(const Point &p) {
    r = 3;
    c = 1;
    v = new Vector(3);
    ( * v ) [ 0 ] = p.getX();
    ( * v ) [ 1 ] = p.getY();
    ( * v ) [ 2 ] = p.getZ();
}

double &Matrix :: operator()(size_t x, size_t y) {
    return row(x) [ y ];
}

double Matrix :: operator()(size_t x, size_t y) const {
    return row(x) [ y ];
}

Matrix Matrix :: transpose() const {
    Matrix ret(c, r);
    for ( size_t i = 0; i < r; i++ ) {
        for ( size_t j = 0; j < c; j++ ) {
            ret [ j ] [ i ] = ( * this ) [ i ] [ j ];
        }
    }

    return ret;
}

Matrix &Matrix :: operator*=(double d) {
    ( * v ) *= d;
    return * this;
}

Matrix Matrix :: operator*(double d) const {
    Matrix ret(* this);
    ret *= d;

    return ret;
}

Matrix Matrix :: operator/(double d) const {
    Matrix ret(* this);
    ret /= d;
    return ret;
}

Matrix &Matrix :: operator/=(double d) {
    ( * v ) /= d;
    return * this;
}

Matrix &Matrix :: operator=(const MtM &m) {
    ( * this ) = ( Matrix ) m;

    return * this;
}

Matrix &Matrix :: operator*=(const Matrix &m) {
    assert(m.numRows() == this->numCols() );

    Matrix ret = matrix_multiply( ( * this ), m );

    ( * v ) = ret.array();
    r = ret.numRows();
    c = ret.numCols();

    return * this;
}

Matrix &Matrix :: operator+=(const Matrix &m) {
    ( * v ) += ( m.array() );
    return * this;
}

Matrix Matrix :: operator+(const Matrix &m) const {
    Matrix ret(* this);
    ret += m;
    ;
    return ret;
}

Matrix &Matrix :: operator+=(const double x) {
    ( * v ) += x;
    return *this;
}

Matrix &Matrix :: operator-=(const Matrix &m) {
    ( * v ) -= ( m.array() );
    return * this;
}

Matrix Matrix :: operator-(const Matrix &m) const {
    Matrix ret(* this);
    ret -= m;
    ;
    return ret;
}

bool Matrix :: operator==(const Matrix &m) {
    if ( v->size() != m.array().size() ) {
        return false;
    }

    for ( size_t i = 0; i < v->size(); i++ ) {
        if ( ( * v ) [ i ] != m.array() [ i ] ) {
            return false;
        }
    }

    return true;
}

bool Matrix :: operator!=(const Matrix &m) {
    if ( v->size() != m.array().size() ) {
        return true;
    } else {
        for ( size_t i = 0; i < v->size(); i++ ) {
            if ( ( * v ) [ i ] != m.array() [ i ] ) {
                return true;
            }
        }
    }

    return false;
}

Matrix :: Matrix(const Matrix &m) : r(m.numRows() ), c(m.numCols() ) {
    v = new Vector(m.array() );
}

Matrix &Matrix :: operator=(const Matrix &m) {
    delete v;
    v = new Vector(m.array() );
    r = m.numRows();
    c = m.numCols();
    return * this;
}

void Matrix :: print() const {
    for ( size_t i = 0; i < numRows(); i++ ) {
        for ( size_t j = 0; j < numCols(); j++ ) {
            cout << ( * this ) [ i ] [ j ] << " " << flush;
        }

        cout << endl;
    }
}

double *Matrix :: GetArray() const {
    return & ( ( * v ) [ 0 ] );
}


MtM :: operator const Matrix() const {
    return matrix_multiply(first, second);
}

MtV :: operator const Vector() {
    return matrix_vector_multiply(m, v);
}


SparseVector :: SparseVector(Vector &v, valarray< unsigned > &i, const size_t l, const size_t s) : val(v), idx(i), length(l), start(s)
{
    zero = 0;
}

double SparseVector :: operator[](size_t i) const
{
    unsigned *__start__       = & idx [ start ];
    unsigned *__end__         = & idx [ start + length ];
    unsigned *i_index_pointer = lower_bound(__start__, __end__, i);
    unsigned offset            = i_index_pointer - __start__;

    if ( binary_search(__start__, __end__, i) ) {
        return val [ start + offset ];
    }

    return 0;
}

double &SparseVector :: operator[](const size_t i)
{
    zero = 0;
    unsigned *__start__       = & idx [ start ];
    unsigned *__end__         = & idx [ start + length ];
    unsigned *i_index_pointer = lower_bound(__start__, __end__, i);
    unsigned offset            = i_index_pointer - __start__;

    if ( binary_search(__start__, __end__, i) ) {
        return val [ start + offset ];
    }

    return zero;
}


double SparseVector :: operator*(const Vector &v) const
{
    if ( length != v.size() ) {
        cerr << "Sparse vector size does not match for multiplication with Vector" << endl;
    }
    double ret = 0;
    for ( size_t j = start; j < length + start; j++ ) {
        size_t index = idx [ j ];
        ret += v [ index ] * val [ j ];
    }
    return ret;
}

SparseVector SparseVector :: operator*(const double d) const
{
    SparseVector ret(* this);
    ret.val = ret.val * d;
    return ret;
}

double SparseVector :: operator*(const SparseVector &v) const
{
    if ( length != v.length ) {
        cerr << "Sparse vector sizes does not match for multiplication" << endl;
    }
    double ret = 0;
    size_t i = 0;
    size_t j = 0;
    while ( i < length && j < v.length ) {
        if ( idx [ start + i ] > v.idx [ v.start + j ] ) {
            j++;
        } else if ( idx [ start + i ] < v.idx [ v.start + j ] ) {
            i++;
        } else {
            ret += val [ start + i ] * v.val [ v.start + j ];
            i++;
            j++;
        }
    }
    return ret;

    //  for(size_t j = 0 ; j < length ; j++)
    //  {
    //      ret += v[(*(idx + j))]*(*(val + j)) ;
    //  }
    //  return ret ;
}

void SparseVector :: print() const
{
    for ( size_t j = 0; j < length; j++ ) {
        cout << idx [ start + j ] << " -> " << val [ start + j ] << endl;
    }
}

double innerProduct(const SparseVector &v0, const SparseVector &v1, const size_t end)
{
    double ret = 0;

    unsigned i = 0;
    unsigned j = 0;
    unsigned *i_index_pointer = lower_bound(& v0.idx [ v0.start ], & v0.idx [ v0.start + v0.length ], end);
    unsigned *j_index_pointer = lower_bound(& v1.idx [ v1.start ], & v1.idx [ v1.start + v1.length ], end);
    unsigned i_end = i_index_pointer - & v0.idx [ v0.start ];
    unsigned j_end = j_index_pointer - & v1.idx [ v1.start ];

    if ( v0.idx [ v0.start ] > v1.idx [ v1.start ] ) {
        j_index_pointer = lower_bound(& v1.idx [ v1.start + j ], & v1.idx [ v1.start + v1.length ], v0.idx [ v0.start + i ]);
        j = j_index_pointer - & v1.idx [ v1.start ];
    } else if ( v0.idx [ v0.start ] < v1.idx [ v1.start ] ) {
        i_index_pointer = lower_bound(& v0.idx [ v0.start + i ], & v0.idx [ v0.start + v0.length ], v1.idx [ v1.start + j ]);
        i = i_index_pointer - & v0.idx [ v0.start ];
    }

    while ( i < i_end &&  j < j_end ) {
        ret += v0.val [ v0.start + i ] * v1.val [ v1.start + j ];
        i++;
        j++;

        if ( v0.idx [ v0.start + i ] > v1.idx [ v1.start + j ] ) {
            j_index_pointer = lower_bound(& v1.idx [ v1.start + j ], & v1.idx [ v1.start + v1.length ], v0.idx [ v0.start + i ]);
            j = j_index_pointer - & v1.idx [ v1.start ];
        } else if ( v0.idx [ v0.start + i ] < v1.idx [ v1.start + j ] ) {
            i_index_pointer = lower_bound(& v0.idx [ v0.start + i ], & v0.idx [ v0.start + v0.length ], v1.idx [ v1.start + j ]);
            i = i_index_pointer - & v0.idx [ v0.start ];
        }
    }
    return ret;
}

double innerProduct(const SparseVector &v0, const Vector &v1, const size_t end)
{
    double ret = 0;
    for ( size_t j = v0.start; v0.idx [ j ] < end /*&& j < v0.length*/; ++j ) {
        ret += v1 [ v0.idx [ j ] ] * v0.val [ j ];
    }
    return ret;
}

double innerProduct(const ConstSparseVector &v0, const Vector &v1, const size_t end)
{
    double ret = 0;
    for ( size_t j =  v0.start; v0.idx [ j ] < end /*&& j < v0.length*/; ++j ) {
        ret += v1 [ v0.idx [ j ] ] * v0.val [ j ];
        //      ret += v1[v0.idx][j]*v0.val[j] ;
    }
    return ret;
}

// inline void innerProductAssignAndAdd(const Mu::ConstSparseVector & v0, Vector & v1, double &t,  double toAdd, size_t end)
// {
//  for(size_t j =  v0.start; v0.idx[j] < end ; ++j)
//  {
//      t += v1[v0.idx[j]]*v0.val[j] ;
//  }
//  t+=toAdd ;
// }

double reverseInnerProduct(const SparseVector &v0, const Vector &v1, const size_t s)
{
    double ret = 0;
    for ( size_t j = v0.length + v0.start - 1; /*j >0 && */ v0.idx [ j ] > s; --j ) {
        ret += v1 [ v0.idx [ j ] ] * v0.val [ j ];
    }
    return ret;
}

// inline void reverseInnerProductAssignAndAdd(const Mu::ConstSparseVector & v0, Vector & v1, double &t,  double toAdd, size_t start)
// {
//  for(size_t j = v0.length+v0.start-1 ; v0.idx[j] > start   ; --j)
//  {
//      t += v1[v0.idx[j]]*v0.val[j] ;
//  }
//  t+=toAdd ;
// }

double reverseInnerProduct(const ConstSparseVector &v0, const Vector &v1, const size_t s)
{
    double ret = 0;

    for ( size_t j = v0.length + v0.start - 1; v0.idx [ j ] > s; --j ) {
        ret += v1 [ v0.idx [ j ] ] * v0.val [ j ];
        //      ret += v1[v0.idx][j]*v0.val[j] ;
    }
    return ret;
}

Vector SparseVector :: operator+(const Vector &v) const
{
    Vector ret(v);

    for ( size_t j = 0; j < length; j++ ) {
        ret [ idx [ start + j ] ] += val [ start + j ];
    }

    return ret;
}


ConstSparseVector :: ConstSparseVector(const Vector &v, const valarray< unsigned > &i, const size_t l, const size_t s) : val(v), idx(i), length(l), start(s)
{}

// double ConstSparseVector::operator [](const size_t i) const
// {
//  const size_t *i_index_pointer = lower_bound(&idx[start], &idx[min(start+length,idx.size())], i) ;
//  size_t offset = i_index_pointer - &idx[start] ;
//  if(i_index_pointer != &idx[min(start+length,idx.size())])
//      return (val[start+offset]) ;
//
//  return 0 ;
// //   for(size_t j = 0 ; j < length ; j++)
// //   {
// //       if((*(idx + j)) == i)
// //           return (*(val + j)) ;
// //       if((*(idx + j)) > i)
// //           return 0 ;
// //   }
// //   return 0 ;
// }

double ConstSparseVector :: vectorMultiply(const double *vecA) const
{
    double ret = 0;

    for ( size_t j = start; j < length + start; j++ ) {
        size_t index = idx [ j ];
        ret += vecA [ index ] * val [ j ];
    }
    return ret;
}


double ConstSparseVector :: operator*(const Vector &v) const
{
    double ret = 0;

    for ( size_t j = start; j < length + start; j++ ) {
        size_t index = idx [ j ];
        ret += v [ index ] * val [ j ];
    }
    return ret;
    //  double ret = 0 ;
    //  for(size_t j = 0 ; j < length ; j++)
    //  {
    //      ret += v[(*(idx + j))]*(*(val + j)) ;
    //  }
    //  return ret ;
}

double ConstSparseVector :: operator*(const SparseVector &v) const
{
    double ret = 0;
    size_t i = start;
    size_t j = v.start;
    while ( i < length + start && j < v.length + v.start ) {
        while ( idx [ i ] > v.idx [ j ] ) {
            j++;
        }

        while ( idx [ i ] < v.idx [ j ] ) {
            i++;
        }

        ret += val [ i ] * v.val [ j ];
        i++;
        j++;
    }
    return ret;
}




double ConstSparseVector :: operator*(const ConstSparseVector &v) const
{
    double ret = 0;
    size_t i = start;
    size_t j = v.start;
    while ( i < length + start && j < v.length + v.start ) {
        while ( idx [ i ] > v.idx [ j ] ) {
            j++;
        }

        while ( idx [ i ] < v.idx [ j ] ) {
            i++;
        }

        ret += val [ i ] * v.val [ j ];
        i++;
        j++;
    }
    return ret;
}

void ConstSparseVector :: print() const
{
    for ( size_t j = 0; j < length; j++ ) {
        cout << idx [ start + j ] << " -> " << val [ start + j ] << endl;
    }
}

Vector ConstSparseVector :: operator+(const Vector &v) const
{
    Vector ret(v);

    for ( size_t j = start; j < length + start; j++ ) {
        size_t index = idx [ j ];
        ret [ index ] += val [ j ];
    }

    return ret;
    //  Vector ret(v) ;
    //
    //  for(size_t j = 0 ; j < length ; j++)
    //  {
    //      ret[(*(idx + j))] += (*(val + j)) ;
    //  }
    //
    //  return ret ;
}


// struct BandSparseVector
// {
// public:
//  Vector & val ;
//  const size_t length ;
//  const size_t start ;
//
//  double zero ;
//
// public:
//  BandSparseVector(Vector & v , const size_t l , const size_t s) ;
//
//  double operator [](const size_t) const ;
//  double & operator [](const size_t) ;
//  double operator *(const Vector&) const ;
//  double operator *(const SparseVector&) const ;
//  Vector operator +(const Vector&) const ;
//
// } ;















CoordinateIndexedSparseMatrix :: CoordinateIndexedSparseMatrix(map< pair< size_t, size_t >, double > &source, unsigned RowCountI, unsigned ColumnCountI) {
    RowCount = RowCountI;
    ColumnCount = ColumnCountI;

    vector< double >temp_array;
    vector< unsigned >temp_column_index;
    vector< unsigned >temp_row_size;
    map< pair< size_t, size_t >, double > :: const_iterator previous =
        source.begin();
    size_t r_s = 0;
    unsigned k = previous->first.first;
    // add initail rows
    while ( k > 0 ) {
        temp_row_size.push_back(r_s);
        k--;
    }
    for ( map< pair< size_t, size_t >, double > :: const_iterator ij =
              source.begin(); ij != source.end(); ++ij ) {
        if ( ij->first.first == previous->first.first ) {
            r_s++;
        } else {
            temp_row_size.push_back(r_s);
            k = ij->first.first - previous->first.first;
            while ( k > 1 ) {
                temp_row_size.push_back(0);
                k--;
            }
            r_s = 1;
        }
        previous = ij;
        temp_array.push_back(ij->second);
        temp_column_index.push_back(ij->first.second);
    }
    temp_row_size.push_back(r_s);

    column_index.resize(temp_column_index.size() );
    copy(temp_column_index.begin(), temp_column_index.end(),
         & column_index [ 0 ]);

    array.resize(temp_array.size() );
    copy(temp_array.begin(), temp_array.end(), & array [ 0 ]);

    row_size.resize(RowCountI);
    copy(temp_row_size.begin(), temp_row_size.end(), & row_size [ 0 ]);

    accumulated_row_size.resize(row_size.size() );
    accumulated_row_size [ 0 ] = 0;
    for ( size_t i = 1; i < accumulated_row_size.size(); i++ ) {
        accumulated_row_size [ i ] += accumulated_row_size [ i - 1 ] + row_size [ i - 1 ];
    }
}

CoordinateIndexedSparseMatrix :: CoordinateIndexedSparseMatrix(map< pair< size_t, size_t >, Matrix > &source, unsigned RowCountI, unsigned ColumnCountI) :
    row_size( ( source.rbegin()->first.first + 1 ) * source.begin()->second.numRows() ) {
    size_t ddl = source.begin()->second.numRows();

    vector< double >temp_array;
    vector< unsigned >temp_column_index;
    map< pair< size_t, size_t >, Matrix > :: const_iterator previous = source.begin();
    vector< vector< double > >to_linerarise(ddl);
    vector< vector< unsigned > >col_to_linerarise(ddl);


    for ( map< pair< size_t, size_t >, Matrix > :: const_iterator ij = source.begin(); ij != source.end(); ) {
        size_t offset = ij->first.first * ddl;
        for ( size_t i = offset; i < ddl + offset; i++ ) {
            row_size [ i ] += ddl;
        }

        if ( ij->first.first != previous->first.first ) {
            for ( size_t i = 0; i < ddl; i++ ) {
                for ( size_t j = 0; j < to_linerarise [ i ].size(); j++ ) {
                    temp_array.push_back(to_linerarise [ i ] [ j ]);
                    temp_column_index.push_back(col_to_linerarise [ i ] [ j ]);
                }
                to_linerarise [ i ].clear();
                col_to_linerarise [ i ].clear();
            }

            for ( size_t i = 0; i < ddl; i++ ) {
                for ( size_t j = 0; j < ddl; j++ ) {
                    to_linerarise [ i ].push_back(ij->second [ i ] [ j ]);
                    col_to_linerarise [ i ].push_back(ij->first.second * ddl + j);
                }
            }

            previous = ij;
        } else {
            for ( size_t i = 0; i < ddl; i++ ) {
                for ( size_t j = 0; j < ddl; j++ ) {
                    to_linerarise [ i ].push_back(ij->second [ i ] [ j ]);
                    col_to_linerarise [ i ].push_back(ij->first.second * ddl + j);
                }
            }

            previous = ij;
        }
        ++ij;
        if ( ij == source.end() ) {
            for ( size_t i = 0; i < ddl; i++ ) {
                for ( size_t j = 0; j < to_linerarise [ i ].size(); j++ ) {
                    temp_array.push_back(to_linerarise [ i ] [ j ]);
                    temp_column_index.push_back(col_to_linerarise [ i ] [ j ]);
                }
                to_linerarise [ i ].clear();
                col_to_linerarise [ i ].clear();
            }
        }
    }

    array.resize(temp_array.size() );
    copy(temp_array.begin(), temp_array.end(), & array [ 0 ]);

    column_index.resize(temp_column_index.size() );
    copy(temp_column_index.begin(), temp_column_index.end(), & column_index [ 0 ]);

    RowCount = RowCountI;
    ColumnCount = ColumnCountI;
}

CoordinateIndexedSparseMatrix :: CoordinateIndexedSparseMatrix(const CoordinateIndexedSparseMatrix &source) {
    this->column_index.resize(source.column_index.size() );
    this->array.resize(source.array.size() );
    this->row_size.resize(source.row_size.size() );
    this->accumulated_row_size.resize(source.accumulated_row_size.size() );

    this->array = source.array;
    this->column_index = source.column_index;
    this->row_size = source.row_size;
    this->accumulated_row_size = source.accumulated_row_size;

    this->RowCount = source.RowCount;
    this->ColumnCount = source.ColumnCount;
}

CoordinateIndexedSparseMatrix :: CoordinateIndexedSparseMatrix(const valarray< unsigned > &rs, const valarray< unsigned > &ci, unsigned RowCountI, unsigned ColumnCountI) : array(0., ci.size() ), column_index(ci), row_size(rs), accumulated_row_size(rs.size() ) {
    //  accumulated_row_size[1] = row_size[0] ;
    for ( size_t i = 1; i < accumulated_row_size.size(); i++ ) {
        accumulated_row_size [ i ] += accumulated_row_size [ i - 1 ] + row_size [ i - 1 ];
    }

    RowCount = RowCountI;
    ColumnCount = ColumnCountI;
}

CoordinateIndexedSparseMatrix :: CoordinateIndexedSparseMatrix() {}

CoordinateIndexedSparseMatrix :: ~CoordinateIndexedSparseMatrix() {}

double CoordinateIndexedSparseMatrix :: froebeniusNorm() const {
    return sqrt(inner_product(& array [ 0 ], & array [ array.size() ], & array [ 0 ], ( double ) ( 0 ) ) );
}

double CoordinateIndexedSparseMatrix :: infinityNorm() const {
    return abs(array).max();
}

void CoordinateIndexedSparseMatrix :: print(unsigned size1, unsigned size2) {
    cout << endl;
    if ( size1 > RowCount ) {
        size1 = RowCount;
    }
    if ( size2 > ColumnCount ) {
        size2 = ColumnCount;
    }

    for ( unsigned i = 0; i < size1; i++ ) {
        if ( i == 0 ) {
            for ( unsigned j = 0; j < size2; j++ ) {
                cout << " " << j << flush;
            }
        }
        cout << endl;
        cout << i << " " << flush;
        for ( unsigned j = 0; j < size2; j++ ) {
            cout << ( * this ) [ i ] [ j ] << " " << flush;
        }
        cout << endl;
    }
    cout << endl;
}

void CoordinateIndexedSparseMatrix :: print() {
    this->print(RowCount, ColumnCount);
}

bool CoordinateIndexedSparseMatrix :: isThereNaN() const {
    for ( auto &k:array ) {
        if ( k != k ) {
            return true;
        }
    }
    return false;
}

SparseVector CoordinateIndexedSparseMatrix :: operator[](const size_t i) {
    if ( i > RowCount ) {
        cerr << "No such row in the matrix" << endl;
    }
    //  size_t start_index = accumulate(&row_size[0], &row_size[i+1], 0) ;
    return SparseVector(array, column_index, row_size [ i ], accumulated_row_size [ i ]);
}

const ConstSparseVector CoordinateIndexedSparseMatrix :: operator[](const size_t i) const {
    if ( i > RowCount ) {
        cerr << "No such row in the matrix" << endl;
    }
    //  size_t start_index = accumulate(&row_size[0], &row_size[i+1], 0) ;
    return ConstSparseVector(array, column_index, row_size [ i ], accumulated_row_size [ i ]);
}


double CoordinateIndexedSparseMatrix :: operator()(const size_t i, const size_t j) const {
    if ( i > RowCount ) {
        cerr << "No such row in the matrix" << endl;
    }
    if ( j > ColumnCount ) {
        cerr << "No such column in the matrix" << endl;
    }

    size_t start_index = accumulated_row_size [ i ];

    for ( size_t k = 0; k < row_size [ i ]; k++ ) {
        if ( column_index [ start_index ] == j ) {
            return array [ start_index ];
        } else {
            start_index++;
        }
    }
    return 0.;
}


Vector CoordinateIndexedSparseMatrix :: operator*(const Vector &v) const {
    if ( v.size() != ColumnCount ) {
        cerr << "Size of sparse matrix did not match the size of the vector" << endl;
    }

    Vector ret(0., RowCount);

    for ( size_t i = 0; i < row_size.size(); i++ ) {
        ret [ i ] += ( * this ) [ i ] * v;
    }

    return ret;
}

void CoordinateIndexedSparseMatrix :: vectorMultiply(const double *vecA, double *vecB) const {
    for ( size_t i = 0; i < row_size.size(); i++ ) {
        vecB [ i ] += ( * this ) [ i ].vectorMultiply(vecA);
    }
}

CoordinateIndexedSparseMatrix CoordinateIndexedSparseMatrix :: operator*(const double d) const {
    CoordinateIndexedSparseMatrix ret(* this);

    ret.array = ret.array * d;

    return ret;
}

Vector CoordinateIndexedSparseMatrix :: inverseDiagonal() const {
    Vector ret(0., min(RowCount, ColumnCount) );

    for ( size_t i = 0; i < min(RowCount, ColumnCount); i++ ) {
        //      double v = (*this)[i][i] ;
        //      if(abs(v) < 1e16)
        ret [ i ] = 1. / ( * this ) [ i ] [ i ];


        //      else
        //          ret[i] = 1 ;
        //      cout << ret[i] << endl ;
    }
    return ret;
}

Vector CoordinateIndexedSparseMatrix :: inverseDiagonalSquared() const {
    Vector ret(0., min(RowCount, ColumnCount) );

    for ( size_t i = 0; i < min(RowCount, ColumnCount); i++ ) {
        ret [ i ] = 1. / ( ( * this ) [ i ] * ( * this ) [ i ] );
        if ( ( ( * this ) [ i ] * ( * this ) [ i ] ) == 0 ) {
            ( * this ) [ i ].print();
        }
    }
    return ret;
}

Vector CoordinateIndexedSparseMatrix :: diagonal() const {
    Vector ret(0., min(RowCount, ColumnCount) );

    for ( size_t i = 0; i < min(RowCount, ColumnCount); i++ ) {
        ret [ i ] = ( * this ) [ i ] [ i ];
        //      if(ret[i] == 0)
        //          cout << i << " is null !!!" << endl ;
    }
    return ret;
}

CoordinateIndexedSparseMatrix &CoordinateIndexedSparseMatrix :: operator=(const CoordinateIndexedSparseMatrix &S) {
    this->column_index.resize(S.column_index.size() );
    this->array.resize(S.array.size() );
    this->row_size.resize(S.row_size.size() );
    this->accumulated_row_size.resize(S.accumulated_row_size.size() );

    this->column_index = S.column_index;
    this->array = S.array;
    this->row_size = S.row_size;
    this->accumulated_row_size = S.accumulated_row_size;

    this->ColumnCount = S.ColumnCount;
    this->RowCount = S.RowCount;

    return * this;
}

CoordinateIndexedSparseMatrix CoordinateIndexedSparseMatrix :: operator*(const CoordinateIndexedSparseMatrix &b) const {
    if ( this->ColumnCount != b.RowCount ) {
        cerr << "Matrix sizes did not match for multiplication" << endl;
    }
    const CoordinateIndexedSparseMatrix *A = this;
    const CoordinateIndexedSparseMatrix *B = & b;
    //if (A->array.size()<B->array.size()){
    //    A = &b;
    //    B = this;
    //}
    map< pair< size_t, size_t >, double >indeces;
    unsigned start_i_A, start_i_B;
    unsigned row_numA = A->row_size.size();
    unsigned row_numB = B->row_size.size();
    unsigned col_numA;
    unsigned col_numB;
    unsigned max_col_numB = B->ColumnCount;
    vector< unsigned >col_index, column_pointer, sum_column;
    vector< double >values, sum;

    for ( unsigned rowA = 0; rowA < row_numA; rowA++ ) {
        values.clear();
        col_index.clear();
        column_pointer.clear();
        column_pointer.resize(max_col_numB);
        sum.clear();
        sum.resize(1);
        sum_column.clear();
        sum_column.resize(1);
        col_numA = A->row_size [ rowA ];
        start_i_A = A->accumulated_row_size [ rowA ];
        for ( unsigned cAi = 0; cAi < col_numA && A->column_index [ start_i_A + cAi ] < row_numB; cAi++ ) {
            col_numB = B->row_size [ A->column_index [ start_i_A + cAi ] ];
            start_i_B = B->accumulated_row_size [ A->column_index [ start_i_A + cAi ] ];
            for ( unsigned i = 0; i < col_numB; i++ ) {
                col_index.push_back(B->column_index [ start_i_B + i ]);
                values.push_back(A->array [ start_i_A + cAi ] * B->array [ start_i_B + i ]);
            }
        }
        for ( unsigned i = 0; i < values.size(); i++ ) {
            if ( column_pointer [ col_index [ i ] ] == 0 ) {
                column_pointer [ col_index [ i ] ] = sum.size();
                sum.push_back(values [ i ]);
                sum_column.push_back(col_index [ i ]);
            } else {
                sum [ column_pointer [ col_index [ i ] ] ] += values [ i ];
            }
        }
        for ( unsigned i = 1; i < sum.size(); i++ ) {
            indeces.insert(pair< pair< size_t, size_t >, double >(pair< size_t, size_t >(rowA, sum_column [ i ]), sum [ i ]) );
        }
    }

    return CoordinateIndexedSparseMatrix(indeces, this->RowCount, b.ColumnCount);
}

CoordinateIndexedSparseMatrix CoordinateIndexedSparseMatrix :: operator+(const CoordinateIndexedSparseMatrix &b) const {
    if ( this->ColumnCount != b.ColumnCount ) {
        cerr << "Matrix sizes did not match for summation" << endl;
    }
    if ( this->RowCount != b.RowCount ) {
        cerr << "Matrix sizes did not match for summation" << endl;
    }
    const CoordinateIndexedSparseMatrix *A = this;
    const CoordinateIndexedSparseMatrix *B = & b;

    map< pair< size_t, size_t >, double >indeces;
    unsigned col, row, ari;

    ari = 0;
    for ( size_t i = 0; i < A->array.size(); i++ ) {
        col = A->column_index [ i ];
        while ( ari < A->row_size.size() - 1 && A->accumulated_row_size [ ari + 1 ] == i ) {
            ari++;
        }
        row = ari;
        indeces.insert(pair< pair< size_t, size_t >, double >(pair< size_t, size_t >(row, col), 0.) );
    }
    ari = 0;
    for ( size_t i = 0; i < B->array.size(); i++ ) {
        col = B->column_index [ i ];
        while ( ari < B->row_size.size() - 1 && B->accumulated_row_size [ ari + 1 ] == i ) {
            ari++;
        }
        row = ari;
        indeces.insert(pair< pair< size_t, size_t >, double >(pair< size_t, size_t >(row, col), 0.) );
    }

    CoordinateIndexedSparseMatrix AplusB(indeces, this->RowCount, this->ColumnCount);

    ari = 0;
    for ( size_t i = 0; i < A->array.size(); i++ ) {
        col = A->column_index [ i ];
        while ( ari < A->row_size.size() - 1 && A->accumulated_row_size [ ari + 1 ] == i ) {
            ari++;
        }
        row = ari;
        AplusB [ row ] [ col ] += A->array [ i ];
    }
    ari = 0;
    for ( size_t i = 0; i < B->array.size(); i++ ) {
        col = B->column_index [ i ];
        while ( ari < B->row_size.size() - 1 && B->accumulated_row_size [ ari + 1 ] == i ) {
            ari++;
        }
        row = ari;
        AplusB [ row ] [ col ] += B->array [ i ];
    }

    return AplusB;
}

CoordinateIndexedSparseMatrix CoordinateIndexedSparseMatrix :: transpose() const {
    map< pair< size_t, size_t >, double >indeces;
    unsigned col, row, ari;

    ari = 0;
    for ( size_t i = 0; i < array.size(); i++ ) {
        col = column_index [ i ];
        while ( ari < row_size.size() - 1 && accumulated_row_size [ ari + 1 ] == i ) {
            ari++;
        }
        row = ari;
        indeces.insert(pair< pair< size_t, size_t >, double >(pair< size_t, size_t >(col, row), array [ i ]) );
    }
    return CoordinateIndexedSparseMatrix(indeces, this->ColumnCount, this->RowCount);
}

CoordinateIndexedSparseMatrix CoordinateIndexedSparseMatrix :: ExtendColumn(const Vector &c) const {
    if ( this->RowCount != c.size() ) {
        cerr << "Matrix and vector sizes did not match" << endl;
    }
    const CoordinateIndexedSparseMatrix *A = this;

    map< pair< size_t, size_t >, double >indeces;
    unsigned col, row, ari;

    ari = 0;
    for ( size_t i = 0; i < A->array.size(); i++ ) {
        col = A->column_index [ i ];
        while ( ari < A->row_size.size() - 1 && A->accumulated_row_size [ ari + 1 ] == i ) {
            ari++;
        }
        row = ari;
        indeces.insert(pair< pair< size_t, size_t >, double >(pair< size_t, size_t >(row, col), A->array [ i ]) );
    }
    for ( unsigned i = 0; i < this->RowCount; i++ ) {
        indeces.insert(pair< pair< size_t, size_t >, double >(pair< size_t, size_t >(this->ColumnCount, i), c [ i ]) );
    }


    return CoordinateIndexedSparseMatrix(indeces, this->RowCount, this->ColumnCount + 1);
}

CoordinateIndexedSparseMatrix CoordinateIndexedSparseMatrix :: ExtendRow(const Vector &c) const {
    if ( this->ColumnCount != c.size() ) {
        cerr << "Matrix and vector sizes did not match" << endl;
    }
    const CoordinateIndexedSparseMatrix *A = this;

    map< pair< size_t, size_t >, double >indeces;
    unsigned col, row, ari;

    ari = 0;
    for ( size_t i = 0; i < A->array.size(); i++ ) {
        col = A->column_index [ i ];
        while ( ari < A->row_size.size() - 1 && A->accumulated_row_size [ ari + 1 ] == i ) {
            ari++;
        }
        row = ari;
        indeces.insert(pair< pair< size_t, size_t >, double >(pair< size_t, size_t >(row, col), A->array [ i ]) );
    }
    for ( size_t i = 0; i < this->ColumnCount; i++ ) {
        indeces.insert(pair< pair< size_t, size_t >, double >(pair< size_t, size_t >(i, this->RowCount), c [ i ]) );
    }


    return CoordinateIndexedSparseMatrix(indeces, this->RowCount + 1, this->ColumnCount);
}

bool ConjGrad(const CoordinateIndexedSparseMatrix &A, Vector &x, const Vector &b, const Vector x0, double precision, double relmaxit) {
    size_t nit = 0;
    size_t Maxit;
    double eps = precision;
    Maxit = b.size() * relmaxit;

    double bnorm = l2_norm(b);

    x.resize(b.size(), 0.);
    if ( x0.size() == b.size() ) {
        x = x0;
    }

    //right hand side is empty
    if ( bnorm < 1E-30 ) {
        x = x * 0.;
        return true;
    }

    //Inverse Diagonal Preconditioner
    Vector preconditioner(A.inverseDiagonal() );

    Vector r = b - A * x;
    double err = l2_norm(r) / bnorm;
    if ( err < eps ) {
        //std::cerr << "b in : " << b.min() << ", " << b.max() << ", err = " << err << std::endl;
        return true;
    }

    Vector z(r);
    for ( size_t i = 0; i < b.size(); i++ ) {
        z [ i ] = r [ i ] * preconditioner [ i ];
    }
    Vector p = z;
    Vector q = A * p;

    double last_rho = std :: inner_product(& r [ 0 ], & r [ r.size() ], & z [ 0 ], ( double ) ( 0 ) );
    double alpha = last_rho / std :: inner_product(& q [ 0 ], & q [ q.size() ], & p [ 0 ], ( double ) ( 0 ) );
    // double rho_0 = last_rho; // unused

    x += p * alpha;
    r -= q * alpha;
    err = l2_norm(r) / bnorm;

    while ( err > eps && nit < Maxit ) {
        for ( size_t i = 0; i < b.size(); i++ ) {
            z [ i ] = r [ i ] * preconditioner [ i ];
        }
        double rho = std :: inner_product(& r [ 0 ], & r [ r.size() ], & z [ 0 ], 0.);
        double beta = rho / last_rho;
        p = z + p * beta;
        q = A * p;
        assert(std :: inner_product(& q [ 0 ], & q [ q.size() ], & p [ 0 ], 0.) != 0);
        alpha = rho / std :: inner_product(& q [ 0 ], & q [ q.size() ], & p [ 0 ], 0.);
        r -= q * alpha;
        x += p * alpha;
        //        if (verbose && nit % 100 == 0) {
        //            r = b - A*x;
        //            std::cerr << "\r iteration : " << nit << " error :" << last_rho << "             " << std::flush;
        //        }
        last_rho = rho;
        nit++;
        err = l2_norm(r) / bnorm;
    }
    r = b - A * x;
    err = l2_norm(r) / bnorm;

    if ( nit == Maxit ) {
        std :: cerr << "\n did not converge after " << nit << " iterations. Error : " << err << ", x norm : " << l2_norm(r) << ", b norm : " << bnorm << std :: endl;
        // for ( unsigned i = 0; i < A.RowCount; i++){
        //   for ( unsigned j = 0; j < A.ColumnCount; j++){
        //     std::cout << "\t" << A(i, j);
        //   }
        //   std::cout << '\n';
        // }
        return false;
    }

    return ( nit < Maxit );
}


bool ConjGrad(const CoordinateIndexedSparseMatrix &A, Vector &x, const Vector &b, const Vector x0) {
    double precision = 1e-16;
    double relmaxit = 0.9;
    if ( b.size() < 50 ) {
        relmaxit = 1.2;               // NOTE using constraint, more iterations are needed (maybe size of original system should be used...?)
    }
    return ConjGrad(A, x, b, x0, precision, relmaxit);
}



bool isMatrixSingular(const CoordinateIndexedSparseMatrix &A) {
    size_t nit = 0;
    size_t Maxit;
    double eps = 1e-26;
    size_t size = A.RowCount;
    Maxit = size * 0.99;

    //check diagonal
    for ( unsigned i = 0; i < Maxit; i++ ) {
        if ( A [ i ] [ i ] < 1e-30 ) {
            return 1;
        }
    }

    Vector b, x;
    b.resize(size, 1.e6);
    x.resize(size, 0.);
    double bnorm = l2_norm(b);

    //Inverse Diagonal Preconditioner
    Vector preconditioner(A.inverseDiagonal() );

    Vector r = b - A * x;
    double err = l2_norm(r) / bnorm;

    Vector z(r);
    for ( size_t i = 0; i < b.size(); i++ ) {
        z [ i ] = r [ i ] * preconditioner [ i ];
    }
    Vector p = z;
    Vector q = A * p;

    double last_rho = std :: inner_product(& r [ 0 ], & r [ r.size() ], & z [ 0 ], ( double ) ( 0 ) );
    double alpha = last_rho / std :: inner_product(& q [ 0 ], & q [ q.size() ], & p [ 0 ], ( double ) ( 0 ) );
    // double rho_0 = last_rho; // unused

    x += p * alpha;
    r -= q * alpha;
    err = l2_norm(r) / bnorm;

    while ( err > eps && nit < Maxit ) {
        for ( size_t i = 0; i < b.size(); i++ ) {
            z [ i ] = r [ i ] * preconditioner [ i ];
        }
        double rho = std :: inner_product(& r [ 0 ], & r [ r.size() ], & z [ 0 ], 0.);
        double beta = rho / last_rho;
        p = z + p * beta;
        q = A * p;
        assert(std :: inner_product(& q [ 0 ], & q [ q.size() ], & p [ 0 ], 0.) != 0);
        alpha = rho / std :: inner_product(& q [ 0 ], & q [ q.size() ], & p [ 0 ], 0.);
        r -= q * alpha;
        x += p * alpha;

        last_rho = rho;
        nit++;
        err = l2_norm(r) / bnorm;
    }
    r = b - A * x;
    err = l2_norm(r) / bnorm;

    if ( nit == Maxit ) {
        return 1;
    }
    return 0;
}


Matrix dyadicProduct(const Vector &a, const Vector &b) {
    Matrix X(a.size(), b.size() );
    for ( unsigned i = 0; i < a.size(); i++ ) {
        for ( unsigned j = 0; j < b.size(); j++ ) {
            X [ i ] [ j ] = a [ i ] * b [ j ];
        }
    }
    return X;
}

double l2_norm(Vector x) {
    return pow(inner_product(& x [ 0 ], & x [ x.size() ], & x [ 0 ], ( double ) ( 0 ) ), 0.5);
}
