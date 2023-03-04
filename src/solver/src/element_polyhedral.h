#ifndef _ELEMENT_POLYHEDRAL_H
#define _ELEMENT_POLYHEDRAL_H

#include <unordered_set>

#include "element.h"


//////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////
class TranspPolygonal : public Element
{
protected:
    double volume;
    Point centroid;
    std :: string ip_type;
    std :: vector< std :: vector< unsigned > >faces;
    std :: vector< Point >normals;
    std :: vector< double >surfaces;

    void prepareGeometry();
    virtual void initIntegration();

public:
    TranspPolygonal(const unsigned dim);
    ~TranspPolygonal() {};
    void readFromLine(std :: istringstream &iss, NodeContainer *fullnodes, MaterialContainer *fullmatrs);
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
    virtual Matrix giveStiffnessMatrix(std :: string matrixType) const;
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
 *  double giveValues(string code, MyVector &result) const {result.resize(0);};
 * };
 *
 *
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
 *  MyVector WachspressShapeF(Point x) const;
 *  MyMatrix WachspressShapeFGrad(Point x) const;
 * public:
 *  TranspPolyhedral(const unsigned dim);
 *  ~TranspPolyhedral() {};
 *  virtual void init();
 *  virtual MyVector shapeF(Point x) const { return WachspressShapeF(x); };
 *  virtual MyMatrix shapeFGrad(Point x) const { return WachspressShapeFGrad(x); };
 *  void readFromLine(istringstream &iss, NodeContainer *fullnodes, MaterialContainer *fullmatrs);
 *  virtual MyMatrix giveConductivityMatrix(string matrixType) const;
 *  virtual MyMatrix giveCapacityMatrix() const;
 *  virtual MyMatrix giveSteadyStateMatrix(string matrixType) const { return giveConductivityMatrix(matrixType); };
 *  virtual MyVector giveInternalForces(const MyVector &DoFs, bool frozen) const;
 * };
 *
 * class TranspVirtPolyhedral : public TranspPolyhedral
 * {
 * private:
 *  MyMatrix V1, V2;
 * public:
 *  TranspVirtPolyhedral(const unsigned dim);
 *  ~TranspVirtPolyhedral() {};
 *  void init();
 *  MyMatrix giveConductivityMatrix(string matrixType) const;
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
 *  MyVector red2full;
 *  vector< double >angles;
 *  MyVector fullTriShapeF(Point x) const;
 *  MyMatrix fullTriShapeFGrad(Point x) const;
 *  MyVector condTriShapeF(Point x) const;
 *  MyMatrix condTriShapeFGrad(Point x) const;
 *
 *  unsigned findFaceNumber(Point x) const;
 * public:
 *  TranspCondensedPolyhedral(const unsigned dim);
 *  ~TranspCondensedPolyhedral() {};
 *  virtual MyVector shapeF(Point x) const { return condTriShapeF(x); };
 *  virtual MyMatrix shapeFGrad(Point x) const { return condTriShapeFGrad(x); };
 *  virtual void init();
 * };
 */
#endif  /* _ELEMENT_POLYHEDRAL_H */
