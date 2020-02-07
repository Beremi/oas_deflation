#ifndef _ELEMENT_POLYHEDRAL_H
#define _ELEMENT_POLYHEDRAL_H

#include "element.h"


//////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////
class TranspPolyhedral : public TransportElement
{
protected:
    double volume;
    Point centroid;
    unsigned nnodes;
    unsigned nfaces;
    vector< vector< unsigned > >faces;
    vector< Point >normals;
    vector< double >areas;
    string ip_type;

    void sort2D();
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
    virtual Vector giveInternalForces(const Vector &DoFs) const;
};

class TranspVirtPolyhedral : public TranspPolyhedral
{
private:
    Matrix V1, V2;
public:
    TranspVirtPolyhedral(const unsigned dim);
    ~TranspVirtPolyhedral() {};
    virtual void init();
    virtual Matrix giveConductivityMatrix(string matrixType) const;
};

//////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////
// TRANSPORT POLYHEDRAL ELEMENT CONSTRUCTED BY STATIC CONDENSATION OF ISOPARAMETRIC TRIANGLES

class TranspCondensedPolyhedral : public TranspPolyhedral
{
private:
    Vector red2full;
    vector< double >angles;
    Vector fullTriShapeF(Point x) const;
    Matrix fullTriShapeFGrad(Point x) const;
    Vector condTriShapeF(Point x) const;
    Matrix condTriShapeFGrad(Point x) const;
public:
    TranspCondensedPolyhedral(const unsigned dim);
    ~TranspCondensedPolyhedral() {};
    virtual Vector shapeF(Point x) const { return condTriShapeF(x); };
    virtual Matrix shapeFGrad(Point x) const { return condTriShapeFGrad(x); };
    virtual void init();
};

#endif  /* _ELEMENT_POLYHEDRAL_H */
