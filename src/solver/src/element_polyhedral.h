#ifndef _ELEMENT_POLYHEDRAL_H
#define _ELEMENT_POLYHEDRAL_H

#include <unordered_set>

#include "element.h"


//////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////
class TranspPolygonal : public TransportElement
{
protected:
    double volume;
    Point centroid;
    unsigned nnodes;
    unsigned nfaces;
    vector< vector< unsigned > >faces;
    vector< Point >normals;
    vector< double >surfaces;
    string ip_type;

    void sort2D();
    void findIntegrationPoints();
    Vector WachspressShapeF(Point x) const;
    Matrix WachspressShapeFGrad(Point x) const;
public:
    TranspPolygonal(const unsigned dim);
    ~TranspPolygonal() {};
    virtual void init();
    virtual Vector shapeF(Point x) const { return WachspressShapeF(x); };
    virtual Matrix shapeFGrad(Point x) const { return WachspressShapeFGrad(x); };
    void readFromLine(istringstream &iss, NodeContainer *fullnodes, MaterialContainer *fullmatrs);
    virtual Matrix giveConductivityMatrix(string matrixType) const;
    virtual Matrix giveCapacityMatrix() const;
    virtual Matrix giveSteadyStateMatrix(string matrixType) const { return giveConductivityMatrix(matrixType); };
    virtual Vector giveInternalForces(const Vector &DoFs, bool frozen) const;
};

class TranspVirtPolygonal : public TranspPolygonal
{
private:
    Matrix V1, V2;
public:
    TranspVirtPolygonal(const unsigned dim);
    ~TranspVirtPolygonal() {};
    virtual void init();
    virtual Matrix giveConductivityMatrix(string matrixType) const;
};

//////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////
// TRANSPORT POLYGONAL ELEMENT CONSTRUCTED BY STATIC CONDENSATION OF ISOPARAMETRIC TRIANGLES

class TranspCondensedPolygonal : public TranspPolygonal
{
private:
    unsigned nodeMaxAngle;
    Vector red2full;
    vector< double >angles;
    Vector fullTriShapeF(Point x) const;
    Matrix fullTriShapeFGrad(Point x) const;
    Vector condTriShapeF(Point x) const;
    Matrix condTriShapeFGrad(Point x) const;

    unsigned findFaceNumber(Point x) const;
public:
    TranspCondensedPolygonal(const unsigned dim);
    ~TranspCondensedPolygonal() {};
    virtual Vector shapeF(Point x) const { return condTriShapeF(x); };
    virtual Matrix shapeFGrad(Point x) const { return condTriShapeFGrad(x); };
    virtual void init();
};


class PolyhedralFace : public GeometricalElement
{
protected:

public:
    PolyhedralFace(const unsigned dim);
    ~PolyhedralFace() {};
    void readFromLine(istringstream &iss, NodeContainer *fullnodes, MaterialContainer *fullmatrs);
    void init();
    double giveValue(string code) const {return 0;};
};

//////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////
class TranspPolyhedral : public TranspPolygonal
{
protected:
    vector< double >volumes;
    vector< Point > faceCenters;

    void findIntegrationPoints();
    Vector WachspressShapeF(Point x) const;
    Matrix WachspressShapeFGrad(Point x) const;
public:
    TranspPolyhedral(const unsigned dim);
    ~TranspPolyhedral() {};
    virtual void init();
    virtual Vector shapeF(Point x) const { return WachspressShapeF(x); };
    virtual Matrix shapeFGrad(Point x) const { return WachspressShapeFGrad(x); };
    void readFromLine(istringstream &iss, NodeContainer *fullnodes, MaterialContainer *fullmatrs);
    virtual Matrix giveConductivityMatrix(string matrixType) const;
    virtual Matrix giveCapacityMatrix() const;
    virtual Matrix giveSteadyStateMatrix(string matrixType) const { return giveConductivityMatrix(matrixType); };
    virtual Vector giveInternalForces(const Vector &DoFs, bool frozen) const;
};
/*
class TranspVirtPolyhedral : public TranspPolyhedral
{
private:
    Matrix V1, V2;
public:
    TranspVirtPolyhedral(const unsigned dim);
    ~TranspVirtPolyhedral() {};
    void init();
    Matrix giveConductivityMatrix(string matrixType) const;
};

//////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////
// TRANSPORT POLYHEDRAL ELEMENT CONSTRUCTED BY STATIC CONDENSATION OF ISOPARAMETRIC TRIANGLES

class TranspCondensedPolyhedral : public TranspPolyhedral
{
private:
    unsigned nodeMaxAngle;
    Vector red2full;
    vector< double >angles;
    Vector fullTriShapeF(Point x) const;
    Matrix fullTriShapeFGrad(Point x) const;
    Vector condTriShapeF(Point x) const;
    Matrix condTriShapeFGrad(Point x) const;

    unsigned findFaceNumber(Point x) const;
public:
    TranspCondensedPolyhedral(const unsigned dim);
    ~TranspCondensedPolyhedral() {};
    virtual Vector shapeF(Point x) const { return condTriShapeF(x); };
    virtual Matrix shapeFGrad(Point x) const { return condTriShapeFGrad(x); };
    virtual void init();
};
*/
#endif  /* _ELEMENT_POLYHEDRAL_H */
