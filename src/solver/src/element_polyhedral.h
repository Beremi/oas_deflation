#ifndef _ELEMENT_POLYGONAL_H
#define _ELEMENT_POLYGONAL_H

#include "element.h"


//////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////
class TranspPolyhedral : public transportElement
{
private:
    double volume;
    Point centroid;
    unsigned nnodes;
    unsigned nfaces;
    vector<vector < unsigned > > faces;
    vector<Point> normals;
    vector<double> areas;

public:
    TranspPolyhedral(const unsigned dim);
    ~TranspPolyhedral() {};
    void init();
    void readFromLine(istringstream &iss, NodeContainer *fullnodes, MaterialContainer *fullmatrs);
    Matrix giveConductivityMatrix(string matrixType) const;
    Matrix giveCapacityMatrix() const;
    Matrix giveSteadyStateMatrix(string matrixType) const { return giveConductivityMatrix(matrixType); };
    Matrix giveTransientMatrix() const { return giveCapacityMatrix(); };
    virtual Vector giveInternalForces(const Vector &DoFs) const;
    void sort2D();

    Vector WachspressShapeF(Point x) const;
    Matrix WachspressShapeFGrad(Point x) const;
};

#endif  /* _ELEMENT_POLYGONAL_H */
