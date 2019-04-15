//
// C++ Implementation: sparse_matrix
//
// Description:
//
//
// Author: Cyrille Dunant <cyrille.dunant@epfl.ch>, (C) 2005-2007
//
// Copyright: See COPYING file that comes with this distribution
//
//
#include<map>
#include "sparse_matrix.h"
using namespace std;

CoordinateIndexedSparseMatrix::CoordinateIndexedSparseMatrix(std::map<std::pair<size_t, size_t>, double> &source, int RowCountI, int ColumnCountI) {

    RowCount = RowCountI;
    ColumnCount = ColumnCountI;

    std::vector<double> temp_array;
    std::vector<unsigned int> temp_column_index;
    std::vector<unsigned int> temp_row_size;
    std::map<std::pair<size_t, size_t>, double>::const_iterator previous =
            source.begin();
    size_t r_s = 0;
    int k = previous->first.first;
    // add initail rows
    while (k > 0){
        temp_row_size.push_back(r_s);
        k--;
    }
    for (std::map<std::pair<size_t, size_t>, double>::const_iterator ij =
            source.begin(); ij != source.end(); ++ij) {
        if (ij->first.first == previous->first.first) {
            r_s++;
        } else {
            temp_row_size.push_back(r_s);
            k = ij->first.first - previous->first.first;
            while (k > 1) {
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

    column_index.resize(temp_column_index.size());
    std::copy(temp_column_index.begin(), temp_column_index.end(),
            &column_index[0]);

    array.resize(temp_array.size());
    std::copy(temp_array.begin(), temp_array.end(), &array[0]);

    row_size.resize(RowCountI);
    std::copy(temp_row_size.begin(), temp_row_size.end(), &row_size[0]);

    accumulated_row_size.resize(row_size.size());
    accumulated_row_size[0] = 0;
    for (size_t i = 1; i < accumulated_row_size.size(); i++) {
        accumulated_row_size[i] += accumulated_row_size[i - 1] + row_size[i - 1];
    }

}

CoordinateIndexedSparseMatrix::CoordinateIndexedSparseMatrix(std::map<std::pair<size_t, size_t>, Matrix> &source, int RowCountI, int ColumnCountI) :
row_size((source.rbegin()->first.first + 1) * source.begin()->second.numRows()) {
    size_t ddl = source.begin()->second.numRows();

    std::vector<double > temp_array;
    std::vector<unsigned int> temp_column_index;
    std::map<std::pair<size_t, size_t>, Matrix>::const_iterator previous = source.begin();
    std::vector<std::vector< double > > to_linerarise(ddl);
    std::vector<std::vector< unsigned int > > col_to_linerarise(ddl);


    for (std::map<std::pair<size_t, size_t>, Matrix>::const_iterator ij = source.begin(); ij != source.end();) {
        size_t offset = ij->first.first*ddl;
        for (size_t i = offset; i < ddl + offset; i++) {
            row_size[i] += ddl;
        }

        if (ij->first.first != previous->first.first) {
            for (size_t i = 0; i < ddl; i++) {
                for (size_t j = 0; j < to_linerarise[i].size(); j++) {
                    temp_array.push_back(to_linerarise[i][j]);
                    temp_column_index.push_back(col_to_linerarise[i][j]);
                }
                to_linerarise[i].clear();
                col_to_linerarise[i].clear();
            }

            for (size_t i = 0; i < ddl; i++) {
                for (size_t j = 0; j < ddl; j++) {
                    to_linerarise[i].push_back(ij->second[i][j]);
                    col_to_linerarise[i].push_back(ij->first.second * ddl + j);
                }
            }

            previous = ij;
        } else {
            for (size_t i = 0; i < ddl; i++) {
                for (size_t j = 0; j < ddl; j++) {
                    to_linerarise[i].push_back(ij->second[i][j]);
                    col_to_linerarise[i].push_back(ij->first.second * ddl + j);
                }
            }

            previous = ij;
        }
        ++ij;
        if (ij == source.end()) {

            for (size_t i = 0; i < ddl; i++) {
                for (size_t j = 0; j < to_linerarise[i].size(); j++) {
                    temp_array.push_back(to_linerarise[i][j]);
                    temp_column_index.push_back(col_to_linerarise[i][j]);
                }
                to_linerarise[i].clear();
                col_to_linerarise[i].clear();
            }
        }
    }

    array.resize(temp_array.size());
    std::copy(temp_array.begin(), temp_array.end(), &array[0]);

    column_index.resize(temp_column_index.size());
    std::copy(temp_column_index.begin(), temp_column_index.end(), &column_index[0]);

    RowCount = RowCountI;
    ColumnCount = ColumnCountI;
}

CoordinateIndexedSparseMatrix::CoordinateIndexedSparseMatrix(const CoordinateIndexedSparseMatrix & source) {
    this->column_index.resize(source.column_index.size());
    this->array.resize(source.array.size());
    this->row_size.resize(source.row_size.size());
    this->accumulated_row_size.resize(source.accumulated_row_size.size());

    this->array = source.array;
    this->column_index = source.column_index;
    this->row_size = source.row_size;
    this->accumulated_row_size = source.accumulated_row_size;

    this->RowCount = source.RowCount;
    this->ColumnCount = source.ColumnCount;
}

CoordinateIndexedSparseMatrix::CoordinateIndexedSparseMatrix(const std::valarray<unsigned int> &rs, const std::valarray<unsigned int> &ci, int RowCountI, int ColumnCountI) : array(0., ci.size()), column_index(ci), row_size(rs), accumulated_row_size(rs.size()) {
    // 	accumulated_row_size[1] = row_size[0] ;
    for (size_t i = 1; i < accumulated_row_size.size(); i++) {
        accumulated_row_size[i] += accumulated_row_size[i - 1] + row_size[i - 1];
    }

    RowCount = RowCountI;
    ColumnCount = ColumnCountI;
}

CoordinateIndexedSparseMatrix::CoordinateIndexedSparseMatrix(){
}

CoordinateIndexedSparseMatrix::~CoordinateIndexedSparseMatrix() {
}

double CoordinateIndexedSparseMatrix::froebeniusNorm() const {
    return sqrt(std::inner_product(&array[0], &array[array.size()], &array[0], (double) (0)));
}

double CoordinateIndexedSparseMatrix::infinityNorm() const {
    return std::abs(array).max();
}

void CoordinateIndexedSparseMatrix::print(int size1, int size2) {
    std::cout << std::endl;
    if(size1>RowCount) size1=RowCount;
    if(size2>ColumnCount) size2=ColumnCount;

    for (int i = 0; i < size1; i++) {
        if (i == 0) {
            for (int j = 0; j < size2; j++) {
                std::cout  << " " << j << std::flush;
            }
        }
        std::cout << std::endl;
        std::cout << i << " " << std::flush;
        for (int j = 0; j < size2; j++) {
            std::cout << (*this)[i][j] << " " << std::flush;
        }
        std::cout << std::endl;
    }
    std::cout << std::endl;
}

void CoordinateIndexedSparseMatrix::print() {
	this->print(RowCount,ColumnCount);
}

SparseVector CoordinateIndexedSparseMatrix::operator[](const size_t i) {
    if(i > RowCount) cerr << "No such row in the matrix" << endl;
    // 	size_t start_index = std::accumulate(&row_size[0], &row_size[i+1], 0) ;
    return SparseVector(array, column_index, row_size[i], accumulated_row_size[i]);
}

const ConstSparseVector CoordinateIndexedSparseMatrix::operator[](const size_t i) const {
    if(i > RowCount) cerr << "No such row in the matrix" << endl;
    // 	size_t start_index = std::accumulate(&row_size[0], &row_size[i+1], 0) ;
    return ConstSparseVector(array, column_index, row_size[i], accumulated_row_size[i]);
}


double CoordinateIndexedSparseMatrix::operator()(const size_t i, const size_t j) const {
    if(i > RowCount) cerr << "No such row in the matrix" << endl;
    if(j > ColumnCount) cerr << "No such column in the matrix" << endl;

    size_t start_index = accumulated_row_size[i];

    for (size_t k = 0; k < row_size[i]; k++) {
        if (column_index[start_index] == j)
            return array[start_index];
        else
            start_index++;
    }
    return 0.;
}


Vector CoordinateIndexedSparseMatrix::operator *(const Vector &v) const {

    if (v.size() != ColumnCount) cerr << "Size of sparse matrix did not match the size of the vector" << endl;

    Vector ret(0., RowCount);

    for (size_t i = 0; i < row_size.size(); i++) {
        ret[i] += (*this)[i] * v;
    }

    return ret;
}

CoordinateIndexedSparseMatrix CoordinateIndexedSparseMatrix::operator *(const double d) const {

    CoordinateIndexedSparseMatrix ret(*this);

    ret.array = ret.array*d;

    return ret;
}

Vector CoordinateIndexedSparseMatrix::inverseDiagonal() const {
    Vector ret(0., min(RowCount,ColumnCount));

    for (size_t i = 0; i < min(RowCount,ColumnCount); i++) {
        // 		double v = (*this)[i][i] ;
        // 		if(std::abs(v) < 1e16)
        ret[i] = 1. / (*this)[i][i];


        // 		else
        // 			ret[i] = 1 ;
        // 		std::cout << ret[i] << std::endl ;
    }
    return ret;
}

Vector CoordinateIndexedSparseMatrix::inverseDiagonalSquared() const {
    Vector ret(0., min(RowCount,ColumnCount));

    for (size_t i = 0; i < min(RowCount,ColumnCount); i++) {
        ret[i] = 1. / ((*this)[i]*(*this)[i]);
        if (((*this)[i]*(*this)[i]) == 0)
            (*this)[i].print();
    }
    return ret;
}

Vector CoordinateIndexedSparseMatrix::diagonal() const {
    Vector ret(0., min(RowCount,ColumnCount));

    for (size_t i = 0; i < min(RowCount,ColumnCount); i++) {
        ret[i] = (*this)[i][i];
        // 		if(ret[i] == 0)
        // 			std::cout << i << " is null !!!" << std::endl ;
    }
    return ret;
}

CoordinateIndexedSparseMatrix & CoordinateIndexedSparseMatrix::operator=(const CoordinateIndexedSparseMatrix &S) {
    this->column_index.resize(S.column_index.size());
    this->array.resize(S.array.size());
    this->row_size.resize(S.row_size.size());
    this->accumulated_row_size.resize(S.accumulated_row_size.size());

    this->column_index = S.column_index;
    this->array = S.array;
    this->row_size = S.row_size;
    this->accumulated_row_size = S.accumulated_row_size;

    this->ColumnCount = S.ColumnCount;
    this->RowCount = S.RowCount;

    return *this;
}

CoordinateIndexedSparseMatrix CoordinateIndexedSparseMatrix::operator *(const CoordinateIndexedSparseMatrix & b) const {

    if(this->ColumnCount != b.RowCount) cerr << "Matrix sizes did not match for multiplication" << endl;
    const CoordinateIndexedSparseMatrix *A = this;
    const CoordinateIndexedSparseMatrix *B = &b;
    //if (A->array.size()<B->array.size()){
    //    A = &b;
    //    B = this;
    //}
    std::map<std::pair<size_t, size_t>, double> indeces;
    unsigned int start_i_A, start_i_B;
    unsigned int row_numA = A->row_size.size();
    unsigned int row_numB = B->row_size.size();
    unsigned int col_numA;
    unsigned int col_numB;
    unsigned int max_col_numB = B->ColumnCount;
    vector<unsigned int> col_index, column_pointer, sum_column;
    vector<double> values, sum;

    for (unsigned int rowA = 0; rowA < row_numA; rowA++) {
        values.clear();
        col_index.clear();
        column_pointer.clear();
        column_pointer.resize(max_col_numB);
        sum.clear();
        sum.resize(1);
        sum_column.clear();
        sum_column.resize(1);
        col_numA = A->row_size[rowA];
        start_i_A = A->accumulated_row_size[rowA];
        for (unsigned int cAi = 0; cAi < col_numA && A->column_index[start_i_A + cAi] < row_numB; cAi++) {
            col_numB = B->row_size[A->column_index[start_i_A + cAi]];
            start_i_B = B->accumulated_row_size[A->column_index[start_i_A + cAi]];
            for (unsigned int i = 0; i < col_numB; i++) {
                col_index.push_back(B->column_index[start_i_B + i]);
                values.push_back(A->array[start_i_A + cAi] * B->array[start_i_B + i]);
            }
        }
        for (unsigned int i = 0; i < values.size(); i++) {
            if (column_pointer[col_index[i]] == 0){
                column_pointer[col_index[i]] = sum.size();
                sum.push_back(values[i]);
                sum_column.push_back(col_index[i]);
            }else{
                sum[column_pointer[col_index[i]]]+=values[i];
            }
        }
        for (unsigned int i = 1; i < sum.size(); i++) {
            indeces.insert(std::pair<std::pair<size_t, size_t>, double>(std::pair<size_t, size_t > (rowA, sum_column[i]), sum[i]));
        }
    }

    return CoordinateIndexedSparseMatrix(indeces,this->RowCount, b.ColumnCount);
}

CoordinateIndexedSparseMatrix CoordinateIndexedSparseMatrix::operator +(const CoordinateIndexedSparseMatrix & b) const {

    if(this->ColumnCount != b.ColumnCount) cerr << "Matrix sizes did not match for summation" << endl;
    if(this->RowCount != b.RowCount) cerr << "Matrix sizes did not match for summation" << endl;
    const CoordinateIndexedSparseMatrix *A = this;
    const CoordinateIndexedSparseMatrix *B = &b;

    std::map<std::pair<size_t, size_t>, double> indeces;
    int col, row, ari;

    ari = 0;
    for (int i = 0; i < A->array.size(); i++){
        col = A->column_index[i];
        while (ari<A->row_size.size()-1 && A->accumulated_row_size[ari+1]==i) ari++;
        row = ari;
        indeces.insert(std::pair<std::pair<size_t, size_t>, double>(std::pair<size_t, size_t > (row, col), 0.));
    }
    ari = 0;
    for (int i = 0; i < B->array.size(); i++){
        col = B->column_index[i];
        while (ari<B->row_size.size()-1 && B->accumulated_row_size[ari+1]==i) ari++;
        row = ari;
        indeces.insert(std::pair<std::pair<size_t, size_t>, double>(std::pair<size_t, size_t > (row, col), 0.));
    }

    CoordinateIndexedSparseMatrix AplusB(indeces, this->RowCount, this->ColumnCount);

    ari = 0;
    for (int i = 0; i < A->array.size(); i++){
        col = A->column_index[i];
        while (ari<A->row_size.size()-1 && A->accumulated_row_size[ari+1]==i) ari++;
        row = ari;
        AplusB[row][col] += A->array[i];
    }
    ari = 0;
    for (int i = 0; i < B->array.size(); i++){
        col = B->column_index[i];
        while (ari<B->row_size.size()-1 && B->accumulated_row_size[ari+1]==i) ari++;
        row = ari;
        AplusB[row][col] += B->array[i];
    }

    return AplusB;
}

/*
CoordinateIndexedSparseMatrix CoordinateIndexedSparseMatrix::operator -(const CoordinateIndexedSparseMatrix & b) const {

    if(this->ColumnCount != b.ColumnCount) cerr << "Matrix sizes did not match for summation" << endl;
    if(this->RowCount != b.RowCount) cerr << "Matrix sizes did not match for summation" << endl;
    const CoordinateIndexedSparseMatrix *A = this;
    const CoordinateIndexedSparseMatrix *B = &b;

    std::map<std::pair<size_t, size_t>, double> indeces;
    int col, row, ari;

    ari = 0;
    for (int i = 0; i < A->array.size(); i++){
        col = A->column_index[i];
        while (ari<A->row_size.size()-1 && A->accumulated_row_size[ari+1]==i) ari++;
        row = ari;
        indeces.insert(std::pair<std::pair<size_t, size_t>, double>(std::pair<size_t, size_t > (row, col), 0.));
    }
    ari = 0;
    for (int i = 0; i < B->array.size(); i++){
        col = B->column_index[i];
        while (ari<B->row_size.size()-1 && B->accumulated_row_size[ari+1]==i) ari++;
        row = ari;
        indeces.insert(std::pair<std::pair<size_t, size_t>, double>(std::pair<size_t, size_t > (row, col), 0.));
    }

    CoordinateIndexedSparseMatrix AminusB(indeces, this->RowCount, this->ColumnCount);

    ari = 0;
    for (int i = 0; i < A->array.size(); i++){
        col = A->column_index[i];
        while (ari<A->row_size.size()-1 && A->accumulated_row_size[ari+1]==i) ari++;
        row = ari;
        AminusB[row][col] += A->array[i];
    }
    ari = 0;
    for (int i = 0; i < B->array.size(); i++){
        col = B->column_index[i];
        while (ari<B->row_size.size()-1 && B->accumulated_row_size[ari+1]==i) ari++;
        row = ari;
        AminusB[row][col] -= B->array[i];
    }

    return AminusB;
}
 */

CoordinateIndexedSparseMatrix CoordinateIndexedSparseMatrix::transpose() const {
    std::map<std::pair<size_t, size_t>, double> indeces;
    int col, row, ari;

    ari = 0;
    for (int i = 0; i < array.size(); i++){
        col = column_index[i];
        while (ari<row_size.size()-1 && accumulated_row_size[ari+1]==i) ari++;
        row = ari;
        indeces.insert(std::pair<std::pair<size_t, size_t>, double>(std::pair<size_t, size_t > (col, row), array[i]));
    }
    return CoordinateIndexedSparseMatrix(indeces, this->ColumnCount, this->RowCount);
}

CoordinateIndexedSparseMatrix CoordinateIndexedSparseMatrix::ExtendColumn(const Vector & c) const {

    if(this->RowCount != c.size()) cerr << "Matrix and vector sizes did not match" << endl;
    const CoordinateIndexedSparseMatrix *A = this;

    std::map<std::pair<size_t, size_t>, double> indeces;
    int col, row, ari;

    ari = 0;
    for (int i = 0; i < A->array.size(); i++){
        col = A->column_index[i];
        while (ari<A->row_size.size()-1 && A->accumulated_row_size[ari+1]==i) ari++;
        row = ari;
        indeces.insert(std::pair<std::pair<size_t, size_t>, double>(std::pair<size_t, size_t > (row, col), A->array[i]));
    }
    for(int i=0; i<this->RowCount; i++)
        indeces.insert(std::pair<std::pair<size_t, size_t>, double>(std::pair<size_t, size_t > (this->ColumnCount, i), c[i]));


    return CoordinateIndexedSparseMatrix(indeces, this->RowCount, this->ColumnCount+1);
}

CoordinateIndexedSparseMatrix CoordinateIndexedSparseMatrix::ExtendRow(const Vector & c) const {

    if(this->ColumnCount != c.size()) cerr << "Matrix and vector sizes did not match" << endl;
    const CoordinateIndexedSparseMatrix *A = this;

    std::map<std::pair<size_t, size_t>, double> indeces;
    int col, row, ari;

    ari = 0;
    for (int i = 0; i < A->array.size(); i++){
        col = A->column_index[i];
        while (ari<A->row_size.size()-1 && A->accumulated_row_size[ari+1]==i) ari++;
        row = ari;
        indeces.insert(std::pair<std::pair<size_t, size_t>, double>(std::pair<size_t, size_t > (row, col), A->array[i]));
    }
    for(int i=0; i<this->ColumnCount; i++)
        indeces.insert(std::pair<std::pair<size_t, size_t>, double>(std::pair<size_t, size_t > (i, this->RowCount), c[i]));


    return CoordinateIndexedSparseMatrix(indeces, this->RowCount+1, this->ColumnCount);
}
