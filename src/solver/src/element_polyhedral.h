#ifndef _ELEMENT_POLYHEDRAL_H
#define _ELEMENT_POLYHEDRAL_H

#include "element.h"


//////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////
class TranspPolyhedral : public transportElement
{
protected:
    double volume;
    Point centroid;
    unsigned nnodes;
    unsigned nfaces;
    vector<vector < unsigned > > faces;
    vector<Point> normals;
    vector<double> areas;
    string ip_type;

    void sort2D();
    void findIntegrationPoints();
    Vector WachspressShapeF(Point x) const;
    Matrix WachspressShapeFGrad(Point x) const;
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
};

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

#endif  /* _ELEMENT_POLYHEDRAL_H */
