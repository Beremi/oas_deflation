#ifndef _ELEMENT_DISCRETE_H
#define _ELEMENT_DISCRETE_H

#include "element.h"
#include "simplex.h"

//////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////
// RBSN ELEMENT
class RigidBodyContact : public Element
{
protected:
    std :: vector< Node * >vert;
    double length, area, perimeter;
    Point normal, centroid, userCentroid;
    Point t1, t2;
    Matrix R;
    bool projectArea, userDefinedCentroid;
    std :: string intPoints;
    bool ignoreNegativeAreas;


    Matrix giveRMatrix() const { return R; };
    virtual void checkNodeType() const;
    virtual void setIntegrationPointsAndWeights();

    std :: vector< Simplex * >simplices;


    virtual void computeMassMatrix();
    virtual void computeDampingMatrix();

public:
    RigidBodyContact(const unsigned dim);
    ~RigidBodyContact() {};
    virtual void readFromLine(std :: istringstream &iss, NodeContainer *fullnodes, MaterialContainer *fullmatrs);
    void init();
    virtual Matrix giveBMatrix(const Point *x) const;
    virtual Matrix giveHMatrix(const Point *x) const;
    std :: vector< Node * >giveVertices() const { return vert; };
    Node *giveVertex(unsigned i) const { return vert [ i ]; };
    unsigned giveNumOfVertices() const { return vert.size(); };
    double giveLength() const { return length; }
    double giveArea() const { return area; }
    virtual Vector giveContactStrainNT() const;
    virtual Vector giveContactStrainXYZ() const;
    virtual Vector giveContactStressXYZ();
    virtual Vector transformToLocal(const Vector &DoFs) const;
    virtual Vector transformToGlobal(const Vector &DoFs) const;
    Vector transformVectorToXYZ(Vector &result) const;
    virtual void giveIPValues(std :: string code, unsigned ipnum, Vector &result) const;
    virtual void giveValues(std :: string code, Vector &result) const;
    Vector giveVectorToNode(const unsigned &node_i, const unsigned &ip_id) const;
    Point giveNormal() const { return normal; };
    Point giveT1() const { return t1; };
    Point giveT2() const { return t2; };
    virtual double giveVolumeAssociatedWithNode(unsigned nodenum) const;
    virtual Vector giveStrain(unsigned i, const Vector &DoFs);
    virtual Matrix giveStiffnessMatrix(std :: string matrixType) const;
    virtual Vector giveInternalForces(const Vector &DoFs, bool frozen, double timeStep);
    virtual Vector integrateLoad(BodyLoad *vl, double time) const;
    virtual Vector integrateInternalSources();
    double givePerimeter()const { return perimeter; };
    virtual Vector giveBoundingBox()const;
    virtual Vector giveFacetBoundingBox()const;
    Point giveCentroid() const {return centroid;};

    virtual void extrapolateIPValuesToNodes(std :: string code, std :: vector< Vector > &result, Vector &weights) const;
    virtual bool isPointInside(Point *xn, const Point *x) const;
};

//////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////
// RBSN ELEMENT
class RigidBodyContactWithRotationalStiffness : public RigidBodyContact
{
protected:
    double Inertia;
public:
    RigidBodyContactWithRotationalStiffness(const unsigned dim);
    ~RigidBodyContactWithRotationalStiffness() {};
    virtual Matrix giveBMatrix(const Point *x) const;
    virtual Matrix giveStiffnessMatrix(std :: string matrixType) const;
    virtual void setIntegrationPointsAndWeights();
    double giveMomentOfInertia() const { return Inertia; }
};

//////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////
// RBSN ELEMENT - rotations do not contribute to openings, tractions do not contribute to couples
class RigidBodyContactWithDecoupledRotationsAndTranslations : public RigidBodyContactWithRotationalStiffness
{
protected:
public:
    RigidBodyContactWithDecoupledRotationsAndTranslations(const unsigned dim);
    ~RigidBodyContactWithDecoupledRotationsAndTranslations() {};
    virtual Matrix giveBMatrix(const Point *x) const;
};

//////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////
// COUPLED RBSN ELEMENT
class RigidBodyContact; //forward declaration
class RigidBodyContactCoupled : public RigidBodyContact
{
protected:
    void extractPressureFromSimplices();
public:
    RigidBodyContactCoupled(const unsigned dim);
    ~RigidBodyContactCoupled() {};
    virtual Vector giveStrain(unsigned i, const Vector &DoFs);
};

//////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////
// RBSN BOUNDARY ELEMENT
class RigidBodyBoundary : public RigidBodyContact
{
protected:
    virtual void checkNodeType() const;
    bool active;
    virtual void computeDampingMatrix() { dampC = Matrix :: Zero( ( this->ndim - 1 ) * 3, ( this->ndim - 1 ) * 3); };
public:
    RigidBodyBoundary(const unsigned dim);
    ~RigidBodyBoundary() {};

    virtual Matrix giveBMatrix(const Point *x) const;
    virtual Matrix giveHMatrix(const Point *x) const;

    void init();
    // virtual MyVector giveInternalForces(const MyVector &DoFs, bool frozen, double timeStep);
    virtual Vector giveStrain(unsigned i, const Vector &DoFs);
    virtual void extrapolateIPValuesToNodes(std :: string code, std :: vector< Vector > &result, Vector &weights) const;
    virtual Matrix giveStiffnessMatrix(std :: string matrixType) const;
};


//////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////
// COUPLED RBSN BOUNDARY ELEMENT
// takes boundary pressure (or similar) and transfer into internal force according to biot coeff
class RigidBodyBoundaryCoupled : public RigidBodyContactCoupled
{
protected:
    virtual void checkNodeType() const;
    virtual void computeMassMatrix() { massM = Matrix :: Zero( ( this->ndim - 1 ) * 3, ( this->ndim - 1 ) * 3); };
    virtual void computeDampingMatrix() { dampC = Matrix :: Zero( ( this->ndim - 1 ) * 3, ( this->ndim - 1 ) * 3); };
public:
    RigidBodyBoundaryCoupled(const unsigned dim);
    ~RigidBodyBoundaryCoupled() {};

    virtual Matrix giveBMatrix(const Point *x) const;
    virtual Matrix giveHMatrix(const Point *x) const;

    void init();
    // virtual MyVector giveInternalForces(const MyVector &DoFs, bool frozen, double timeStep);
    virtual Vector giveStrain(unsigned i, const Vector &DoFs);
    virtual void extrapolateIPValuesToNodes(std :: string code, std :: vector< Vector > &result, Vector &weights) const;
};

//////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////
// TRUSS ELEMENT
class Truss : public RigidBodyContact
{
protected:
    virtual void checkNodeType() const;
public:
    Truss(const unsigned dim) : RigidBodyContact(dim) { name = "Truss"; };
    ~Truss() {};
    virtual Vector giveContactStrainNT() const;
    virtual Matrix giveBMatrix(const Point *x) const;
    virtual Matrix giveHMatrix(const Point *x) const;
};


//////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////
// Discrete TRANSPORT ELEMENT
class DiscreteTrsprtElem : public Element
{
protected:
    std :: vector< Node * >vert;
    bool bound;
    Point normal;
    double length, area;
    bool BolanderCapacityMatrix;
    bool ignoreNegativeAreas;
    bool projectArea;

    virtual void checkNodeAndMaterialType() const;
    virtual void setIntegrationPointsAndWeights();

    virtual void computeDampingMatrix();
    virtual void computeMassMatrix();

public:
    DiscreteTrsprtElem(const unsigned dim);
    ~DiscreteTrsprtElem() {};
    void init();
    double giveArea() const { return area; }
    Point giveNormal() const { return normal; }
    double giveLength() const { return length; }
    virtual double giveVolumeAssociatedWithNode(unsigned nodenum) const;
    virtual void readFromLine(std :: istringstream &iss, NodeContainer *fullnodes, MaterialContainer *fullmatrs);
    Vector giveVectorToNode(const unsigned &node_i, const unsigned &ip_id) const;
    virtual Matrix giveBMatrix(const Point *x) const;
    virtual Matrix giveHMatrix(const Point *x) const;
    virtual Matrix giveStiffnessMatrix(std :: string matrixType) const;
    virtual Vector giveInternalForces(const Vector &DoFs, bool frozen, double timeStep);
    virtual Vector integrateLoad(BodyLoad *vl, double time) const;
    virtual Vector integrateInternalSources();
    virtual Vector giveStrain(unsigned i, const Vector &DoFs);
    std :: vector< Node * >giveVertices() const { return vert; };
    virtual void extrapolateIPValuesToNodes(std :: string code, std :: vector< Vector > &result, Vector &weights) const;
    virtual bool isPointInside(Point *xn, const Point *x) const { ( void ) xn; ( void ) x; return false; }; //TODO: discrete elements does not interpolate
    virtual void giveValues(std :: string code, Vector &result) const;
    virtual void giveIPValues(std :: string code, unsigned ipnum, Vector &result) const;
};

//////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////
// 1D TRANSPORT ELEMENT
class DiscreteHeatConductionElem : public DiscreteTrsprtElem
{
protected:

public:
    DiscreteHeatConductionElem(const unsigned dim);
    virtual void checkNodeAndMaterialType() const;
    ~DiscreteHeatConductionElem() {};
};

//////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////
// COUPLED DISCRETE TRANSPORT ELEMENT
class DiscreteTrsprtCoupledElem : public DiscreteTrsprtElem
{
private:
    std :: vector< RigidBodyContact * >friends; //mechanical elements involved in computation
    std :: vector< double >friendsweight;  //weight of mechanical elements

    void findFriendsFromSimplices();

public:
    DiscreteTrsprtCoupledElem(const unsigned dim) : DiscreteTrsprtElem(dim) { name = "DicscreteTrsprtCoupledElem"; solution_order = 1; }; //coupled elements must be solved after all RBSN elements
    ~DiscreteTrsprtCoupledElem() {};
    virtual void giveValues(std :: string code, Vector &result) const;
    void init();
    virtual Vector giveStrain(unsigned i, const Vector &DoFs);
    void addNewFriend(RigidBodyContact *f, double weight);
    unsigned giveNumOfFriends() const { return friends.size(); };
    virtual void collectInformationsFromNeigborhood();
};

//////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////
// RBSN + TEMPERATURE ELEMENT
class RigidBodyContactWithHeatConduction : public RigidBodyContact
{
protected:
    virtual void computeMassMatrix();
    virtual void computeDampingMatrix();
public:
    RigidBodyContactWithHeatConduction(const unsigned dim);
    ~RigidBodyContactWithHeatConduction() {};
    void init();
    virtual Vector giveStrain(unsigned i, const Vector &DoFs);
    virtual Matrix giveBMatrix(const Point *x) const;
    virtual Matrix giveHMatrix(const Point *x) const;
    virtual void giveValues(std :: string code, Vector &result) const;
    virtual Matrix giveStiffnessMatrix(std :: string matrixType) const;
    virtual void checkNodeType() const;
};

#endif  /* _ELEMENT_DISCRETE_H */
