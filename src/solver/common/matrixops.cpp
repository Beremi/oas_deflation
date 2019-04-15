
// Author: Cyrille Dunant <cyrille.dunant@epfl.ch>, (C) 2005-2007
//
// Copyright: See COPYING file that comes with this distribution
//

#include "matrixops.h"

Matrix::Matrix(size_t x, size_t y) {
    r = x;
    c = y;
    v = new Vector(0., x * y);
}

Matrix::Matrix(size_t rl, size_t k, size_t l, const Matrix & m) {
    r = rl;
    c = rl;
    v = new Vector(rl * rl);
    for (size_t i = 0; i < rl; i++) {
        for (size_t j = 0; j < rl; j++) {
            (*v)[i * rl + j] = m[k + i][l + j];
        }
    }
}

Matrix::Matrix(size_t i, size_t j, const Matrix & m) {
    r = m.numRows() - 1;
    c = m.numCols() - 1;
    v = new Vector(r * c);
    for (size_t l = 0; l < m.numRows(); l++) {
        if (l != i) {
            for (size_t k = 0; k < m.numCols(); k++) {
                if (k != j) {
                    size_t r_index = l;
                    if (l > i)
                        r_index--;
                    size_t c_index = k;
                    if (k > j)
                        c_index--;

                    (*this)[r_index][c_index] = m[l][k];
                }
            }
        }
    }
}

Matrix::Matrix(const Point &p) {
    r = 3;
    c = 1;
    v = new Vector(3);
    (*v)[0] = p.getX();
    (*v)[1] = p.getY();
    (*v)[2] = p.getZ();
}

double& Matrix::operator()(size_t x, size_t y) {
    return row(x)[y];
}

double Matrix::operator()(size_t x, size_t y) const {
    return row(x)[y];
}

Matrix Matrix::transpose() const {
    Matrix ret(c, r);
    for (size_t i = 0; i < r; i++) {
        for (size_t j = 0; j < c; j++) {
            ret[j][i] = (*this)[i][j];
        }
    }

    return ret;
}

Matrix &Matrix::operator*=(double d) {
    (*v) *= d;
    return *this;
}

Matrix Matrix::operator*(double d) const {
    Matrix ret(*this);
    ret *= d;

    return ret;
}

Matrix Matrix::operator/(double d) const {
    Matrix ret(*this);
    ret /= d;
    return ret;
}

Matrix &Matrix::operator/=(double d) {
    (*v) /= d;
    return *this;
}

Matrix & Matrix::operator =(const MtM& m) {
    (*this) = (Matrix) m;

    return *this;
}

Matrix &Matrix::operator*=(const Matrix &m) {
    assert(m.numRows() == this->numCols());

    Matrix ret = matrix_multiply((*this), m);

    (*v) = ret.array();
    r = ret.numRows();
    c = ret.numCols();

    return *this;
}

Matrix &Matrix::operator +=(const Matrix &m) {
    (*v) += (m.array());
    return *this;
}

Matrix Matrix::operator +(const Matrix &m) const {
    Matrix ret(*this);
    ret += m;
    ;
    return ret;
}

Matrix &Matrix::operator -=(const Matrix &m) {
    (*v) -= (m.array());
    return *this;
}

Matrix Matrix::operator -(const Matrix &m) const {
    Matrix ret(*this);
    ret -= m;
    ;
    return ret;
}

bool Matrix::operator ==(const Matrix &m) {
    if (v->size() != m.array().size())
        return false;

    for (size_t i = 0; i < v->size(); i++)
        if ((*v)[i] != m.array()[i])
            return false;

    return true;
}

bool Matrix::operator !=(const Matrix &m) {
    if (v->size() != m.array().size())
        return true;
    else
        for (size_t i = 0; i < v->size(); i++)
            if ((*v)[i] != m.array()[i])
                return true;

    return false;
}

Matrix::Matrix(const Matrix& m) : r(m.numRows()), c(m.numCols()) {
    v = new Vector(m.array());
}

Matrix &Matrix::operator =(const Matrix &m) {
    delete v;
    v = new Vector(m.array());
    r = m.numRows();
    c = m.numCols();
    return *this;
}

void Matrix::print() const {
    for (size_t i = 0; i < numRows(); i++) {
        for (size_t j = 0; j < numCols(); j++) {
            std::cout << (*this)[i][j] << " " << std::flush;
        }

        std::cout << std::endl;
    }
}

double * Matrix::GetArray() const {
    return &((*v)[0]);
}


MtM::operator const Matrix() const {
    return matrix_multiply(first, second);
}

MtV::operator const Vector() {
    return matrix_vector_multiply(m, v);
}

Matrix transfMatrixRot(const Point &x, const Point &y, const Point &z){
    /**
     *  to calculate coordinates expressed in basis (x, y, z)
     *  x, y, z must be orthonormal triad !!
     */

    if ((z - cross(x, y)).norm() >= 1e-10){
        std::cout << "x, y, z do not form an orthonormal triad" << "\n";
        exit(1);
    }

    Matrix trans(3, 3);
    trans(0, 0) = x.getX();
    trans(0, 1) = x.getY();
    trans(0, 2) = x.getZ();
    trans(1, 0) = y.getX();
    trans(1, 1) = y.getY();
    trans(1, 2) = y.getZ();
    trans(2, 0) = z.getX();
    trans(2, 1) = z.getY();
    trans(2, 2) = z.getZ();

    return trans;
}

Matrix transfMatrixRot(const Point &norm){
    /**
     *  to transform into coords viewed from direction of vector norm to a plane (n)
     *  multiplilcation of this matrix only gives you a rotation
     *  to transform coords, linear shift is also needed (a' = Q^t a + u) !!
     */
    Point m;
    Point n = norm.normalized();
    if (n.getX() != 0) m = Point((-n.getZ()-n.getY())/n.getX(), 1., 1.).normalized();
    else if (n.getZ() != 0) m = Point(1., 1., (-n.getX()-n.getY())/n.getZ()).normalized();
    else m = Point(1., ((-n.getX()-n.getZ())/n.getY()), 1.).normalized();
    Point l = cross(n, m).normalized();
    /*
    n.print();
    m.print();
    l.print();
    //*/

    Point n2 = cross(m, l);
    if (n2 != n) std::cout << "wrong transformation n" << "\n";
    Point m2 = cross(l, n);
    if (m2 != m) std::cout << "wrong transformation m" << "\n";
    Point l2 = cross(n, m);
    if (l2 != l) std::cout << "wrong transformation m" << "\n";

    Matrix trans(3, 3);
    trans(0, 0) = m.getX();
    trans(0, 1) = m.getY();
    trans(0, 2) = m.getZ();
    trans(1, 0) = l.getX();
    trans(1, 1) = l.getY();
    trans(1, 2) = l.getZ();
    trans(2, 0) = n.getX();
    trans(2, 1) = n.getY();
    trans(2, 2) = n.getZ();

    return trans;
}

Point pointTranformed(const Matrix &rotation, const Point &translation, const Point &original){
    // translation first !!!
    Point point = original + translation;

    std::vector<double> point_coords;
    for (int i = 0; i < 3; i++){
        point_coords.push_back(rotation(i,0)*point.getX() + rotation(i,1)*point.getY() + rotation(i,2)*point.getZ());
    }
    point.setX(point_coords[0]);
    point.setY(point_coords[1]);
    point.setZ(point_coords[2]);

    return point;
}
