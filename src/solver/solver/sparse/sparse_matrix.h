/**
 * @Author: jose
 * @Date:   2019-04-05T19:13:10+02:00
 * @Last modified by:   jose
 * @Last modified time: 2019-04-05T19:44:40+02:00
 */



//
// C++ Interface: sparse_matrix
//
// Description:
//
//
// Author: Cyrille Dunant <cyrille.dunant@epfl.ch>, (C) 2005-2007
//
// Copyright: See COPYING file that comes with this distribution
//
//

#ifndef __SPARSE_MATRIX_H
#define __SPARSE_MATRIX_H

#include<valarray>
#include<map>
#include <complex>
#include "../../common/matrixops.h"
#include "sparse_vector.h"
#include <algorithm>



class CoordinateIndexedSparseMatrix
{
public:
	Vector  array ;
  int RowCount;
  int ColumnCount;
	std::valarray<unsigned int>  column_index ;
	std::valarray<unsigned int>  row_size ;
	std::valarray<unsigned int>  accumulated_row_size ;
public:
	CoordinateIndexedSparseMatrix(std::map<std::pair<size_t, size_t>, double> &source, int RowCountI, int ColumnCountI) ;
	CoordinateIndexedSparseMatrix(std::map<std::pair<size_t, size_t>, Matrix> &source, int RowCountI, int ColumnCountI) ;
	CoordinateIndexedSparseMatrix(const std::valarray<unsigned int> &, const std::valarray<unsigned int> &, int RowCountI, int ColumnCountI) ;
  CoordinateIndexedSparseMatrix(const CoordinateIndexedSparseMatrix & source);
  CoordinateIndexedSparseMatrix();

	~CoordinateIndexedSparseMatrix() ;

	SparseVector operator[](const size_t i) ;
	const ConstSparseVector operator[](const size_t i) const;
	double  operator()(const size_t i, const size_t j) const ;

	Vector operator *(const Vector & v) const ;
	CoordinateIndexedSparseMatrix operator *(const double d) const ;
  CoordinateIndexedSparseMatrix operator *(const CoordinateIndexedSparseMatrix & b) const ;
  CoordinateIndexedSparseMatrix operator +(const CoordinateIndexedSparseMatrix & b) const ;
  //CoordinateIndexedSparseMatrix operator -(const CoordinateIndexedSparseMatrix & b) const ;
	CoordinateIndexedSparseMatrix & operator=(const CoordinateIndexedSparseMatrix &) ;

  CoordinateIndexedSparseMatrix ExtendColumn(const Vector & c) const;
  CoordinateIndexedSparseMatrix ExtendRow(const Vector & c) const;

	Vector inverseDiagonal() const ;
	Vector inverseDiagonalSquared() const ;
	Vector diagonal() const ;
	CoordinateIndexedSparseMatrix transpose() const ;

  double froebeniusNorm() const ;
	double infinityNorm() const ;

  void print(int size1, int size2);
  void print();

} ;

#endif
