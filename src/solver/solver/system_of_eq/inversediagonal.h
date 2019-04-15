/**
 * @Author: jose
 * @Date:   2019-04-05T19:13:10+02:00
 * @Last modified by:   jose
 * @Last modified time: 2019-04-05T19:43:48+02:00
 */



//
// C++ Interface: inversediagonal
//
// Description:
//
//
// Author: Cyrille Dunant <cyrille.dunant@epfl.ch>, (C) 2007
//
// Copyright: See COPYING file that comes with this distribution
//
//

#ifndef __INV_DIAG_H
#define __INV_DIAG_H

#include "preconditionners.h"
#include "../sparse/sparse_matrix.h"

struct InverseDiagonal : public Preconditionner {
    Vector * diagonal;

    virtual ~InverseDiagonal() {
        delete diagonal;
    }
    InverseDiagonal(const CoordinateIndexedSparseMatrix & A);
    virtual void precondition(const Vector &v, Vector &) const;
};

struct InverseLumpedDiagonal : public Preconditionner {
    Vector * diagonal;

    virtual ~InverseLumpedDiagonal() {
        delete diagonal;
    }
    InverseLumpedDiagonal(const CoordinateIndexedSparseMatrix & A);
    virtual void precondition(const Vector &v, Vector &) const;
};

struct InverseDiagonalSquared : public Preconditionner {
    Vector * diagonal;

    virtual ~InverseDiagonalSquared() {
        delete diagonal;
    }
    InverseDiagonalSquared(const CoordinateIndexedSparseMatrix & A);
    virtual void precondition(const Vector &v, Vector &) const;
};

#endif
