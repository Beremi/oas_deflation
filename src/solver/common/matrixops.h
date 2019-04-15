/**
 * @Author: Cyrille Dunant <cyrille.dunant@epfl.ch>, (C) 2005
 * @Date:   2019-04-05T18:20:54+02:00
 * @Last modified by:   jose
 * @Last modified time: 2019-04-05T20:07:00+02:00
 * Copyright: See COPYING file that comes with this distribution
 */


#ifndef __MATRIXOPS_H__
#define __MATRIXOPS_H__

#include <valarray>
#include "sliceiters.h"
#include <iostream>
#include <vector>
#include <numeric>
#include "point.h"

#include <assert.h>

typedef std::valarray<double> Vector;


struct MtM;

class Matrix {

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
    Matrix(size_t rl, size_t k, size_t l, const Matrix & m);


    /** construct a matrix from a given source matrix without specified row and column.
     *
     * @param i
     * @param j
     * @param m
     */
    Matrix(size_t i, size_t j, const Matrix & m);


    Matrix(const Matrix&);

    Matrix(const Point &p);

    virtual ~Matrix() {
        delete v;
    }

    Matrix & operator=(const Matrix &m);

    size_t size() const {
        return r*c;
    }

    size_t numCols() const {
        return c;
    }

    size_t numRows() const {
        return r;
    }

    Slice_iter< double > column(size_t i) {
        return Slice_iter< double >(v, std::slice(i, r, c));
    }

    Cslice_iter< double > column(size_t i) const {
        return Cslice_iter< double >(v, std::slice(i, r, c));
    }

    Slice_iter< double > row(size_t i) {
        return Slice_iter< double >(v, std::slice(i*c, c, 1));
    }

    Cslice_iter< double > row(size_t i) const {
        return Cslice_iter< double >(v, std::slice(i*c, c, 1));
    }

    Matrix transpose() const;

    double& operator()(size_t x, size_t y);
    double operator()(size_t x, size_t y) const;

    Slice_iter< double > operator()(size_t i) {
        return row(i);
    }

    Cslice_iter< double > operator()(size_t i) const {
        return row(i);
    }

    Slice_iter< double > operator[](size_t i) {
        return row(i);
    }

    Cslice_iter< double > operator[](size_t i) const {
        return row(i);
    }

    Matrix & operator *=(double);
    Matrix operator *(double) const;
    Matrix operator /(double) const;
    Matrix & operator /=(double);
    Matrix & operator *=(const Matrix &m);
    Vector & operator *=(const Vector &m);
    // 	const Matrix& operator *(const Matrix &m) const ;
    Matrix & operator +=(const Matrix &m);
    Matrix operator +(const Matrix &m) const;
    Matrix & operator -=(const Matrix &m);
    Matrix operator -(const Matrix &m) const;
    Matrix & operator =(const MtM& m);

    bool operator ==(const Matrix &m);
    bool operator !=(const Matrix &m);

    Vector &array() {
        return *v;
    }

    Vector array() const {
        return *v;
    }

    void print() const;

    double * GetArray() const;
};

struct MtV {
    const Matrix &m;
    const Vector &v;

    MtV(const Matrix &mm, const Vector & vv) : m(mm), v(vv) {
    }

    operator const Vector();
};

struct MtM {
    const Matrix &first;
    const Matrix &second;

    MtM(const Matrix &mm, const Matrix & mmm) : first(mm), second(mmm) {
    }

    operator const Matrix() const;
};

Vector solveSystem(const Matrix & A, const Vector & b, Vector & x);

inline MtV operator*(const Matrix& mm, const Vector& v) {
    return MtV(mm, v);
};

inline MtM operator*(const Matrix& mm, const Matrix& mmm) {
    return MtM(mm, mmm);
};

inline const Matrix matrix_multiply(const Matrix &m0, const Matrix &m1) {
    assert(m0.numCols() == m1.numRows());

    Matrix ret(m0.numRows(), m1.numCols());

    for (size_t i = 0; i < m0.numRows(); i++) {
        for (size_t j = 0; j < m1.numCols(); j++) {
            const Cslice_iter<double>& ri = m0.row(i);
            const Cslice_iter<double>& cj = m1.column(j);
            ret[i][j] = std::inner_product(&ri[0], &ri[m0.numCols()], cj, (double) (0));
        }
    }
    return ret;
}

inline const Vector matrix_vector_multiply(const Matrix &m, const Vector &v) {
    assert(m.numCols() == v.size());

    Vector ret(v.size());

    for (size_t i = 0; i < m.numRows(); i++) {

        const Cslice_iter<double>& ri = m.row(i);
        ret[i] = std::inner_product(ri, ri.end(), &v[0], (double) (0));
    }
    return ret;
}

inline const Vector operator*(const Vector &v, const Matrix &m) {
    assert(m.numRows() == v.size());

    Vector ret(v.size());

    for (size_t i = 0; i < m.numCols(); i++) {

        const Cslice_iter<double>& ri = m.column(i);
        ret[i] = std::inner_product(ri, ri.end(), &v[0], (double) (0));
    }
    return ret;

}

Matrix transfMatrixRot(const Point &x, const Point &y, const Point &z);
Matrix transfMatrixRot(const Point &normal);
Point pointTranformed(const Matrix &rotation, const Point &translation, const Point &point);

#endif  // __MATRIXOPS_H__
