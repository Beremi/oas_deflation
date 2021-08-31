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
    string ip_type;
    vector < vector < unsigned > > faces;
    vector < Point > normals;
    vector < double > surfaces;

    void prepareGeometry();
    virtual void initIntegration();

public:
    TranspPolygonal(const unsigned dim);
    ~TranspPolygonal() {};
    void readFromLine(istringstream &iss, NodeContainer *fullnodes, MaterialContainer *fullmatrs);
    virtual void init();
    virtual Matrix giveBMatrix(const Point *x) const;
    virtual Matrix giveHMatrix(const Point *x) const;
    virtual void setIntegrationPointsAndWeights();
};


//////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////
class TranspVirtPolygonal : public TranspPolygonal
{
protected:
    Matrix V1, V2, W1, W2;
public:
    TranspVirtPolygonal(const unsigned dim);
    ~TranspVirtPolygonal() {};
    virtual void setIntegrationPointsAndWeights();
    virtual Matrix giveStiffnessMatrix(string matrixType) const;
    virtual Matrix giveDampingMatrix() const;
    virtual Vector giveInternalForces(const Vector &DoFs, bool frozen, double timeStep);
};


//////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////
// TRANSPORT POLYGONAL ELEMENT CONSTRUCTED BY STATIC CONDENSATION OF ISOPARAMETRIC TRIANGLES

class TranspCondensedPolygonal : public TranspPolygonal
{
private:
    unsigned findFaceNumber(Point x) const;
public:
    TranspCondensedPolygonal(const unsigned dim);
    ~TranspCondensedPolygonal() {};
    virtual void initIntegration();
};

/*
//////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////
// TRANSPORT POLYGONAL ELEMENT CONSTRUCTED BY STATIC CONDENSATION OF ISOPARAMETRIC TRIANGLES

/*
 * class PolyhedralFace : public GeometricalElement
 * {
 * protected:
 *
 * public:
 *  PolyhedralFace(const unsigned dim);
 *  ~PolyhedralFace() {};
 *  void readFromLine(istringstream &iss, NodeContainer *fullnodes, MaterialContainer *fullmatrs);
 *  void init();
 *  double giveValue(string code) const {return 0;};
 * };
 */

/*
 * //////////////////////////////////////////////////////////
 * //////////////////////////////////////////////////////////
 * class TranspPolyhedral : public TranspPolygonal
 * {
 * protected:
 *  vector< double >volumes;
 *  vector< Point >faceCenters;
 *  vector< vector< int > >faceConnectivity;
 *  vector< double >determinants;
 *
 *  void findIntegrationPoints();
 *  Vector WachspressShapeF(Point x) const;
 *  Matrix WachspressShapeFGrad(Point x) const;
 * public:
 *  TranspPolyhedral(const unsigned dim);
 *  ~TranspPolyhedral() {};
 *  virtual void init();
 *  virtual Vector shapeF(Point x) const { return WachspressShapeF(x); };
 *  virtual Matrix shapeFGrad(Point x) const { return WachspressShapeFGrad(x); };
 *  void readFromLine(istringstream &iss, NodeContainer *fullnodes, MaterialContainer *fullmatrs);
 *  virtual Matrix giveConductivityMatrix(string matrixType) const;
 *  virtual Matrix giveCapacityMatrix() const;
 *  virtual Matrix giveSteadyStateMatrix(string matrixType) const { return giveConductivityMatrix(matrixType); };
 *  virtual Vector giveInternalForces(const Vector &DoFs, bool frozen) const;
 * };
 */

/*
 * class TranspVirtPolyhedral : public TranspPolyhedral
 * {
 * private:
 *  Matrix V1, V2;
 * public:
 *  TranspVirtPolyhedral(const unsigned dim);
 *  ~TranspVirtPolyhedral() {};
 *  void init();
 *  Matrix giveConductivityMatrix(string matrixType) const;
 * };
 *
 * //////////////////////////////////////////////////////////
 * //////////////////////////////////////////////////////////
 * // TRANSPORT POLYHEDRAL ELEMENT CONSTRUCTED BY STATIC CONDENSATION OF ISOPARAMETRIC TRIANGLES
 *
 * class TranspCondensedPolyhedral : public TranspPolyhedral
 * {
 * private:
 *  unsigned nodeMaxAngle;
 *  Vector red2full;
 *  vector< double >angles;
 *  Vector fullTriShapeF(Point x) const;
 *  Matrix fullTriShapeFGrad(Point x) const;
 *  Vector condTriShapeF(Point x) const;
 *  Matrix condTriShapeFGrad(Point x) const;
 *
 *  unsigned findFaceNumber(Point x) const;
 * public:
 *  TranspCondensedPolyhedral(const unsigned dim);
 *  ~TranspCondensedPolyhedral() {};
 *  virtual Vector shapeF(Point x) const { return condTriShapeF(x); };
 *  virtual Matrix shapeFGrad(Point x) const { return condTriShapeFGrad(x); };
 *  virtual void init();
 * };
 */
#endif  /* _ELEMENT_POLYHEDRAL_H */
