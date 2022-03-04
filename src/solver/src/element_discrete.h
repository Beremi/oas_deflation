#ifndef _ELEMENT_DISCRETE_H
#define _ELEMENT_DISCRETE_H

#include "element.h"
#include "simplex.h"

//////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////
// RBSN ELEMENT
class RigidBodyContact : public MechanicalElement
{
protected:
    vector< Node * >vert;
    double length, area, perimeter;
    Point normal, t1, t2;
    Matrix R;

    Matrix giveRMatrix() const { return R; };
    virtual void checkNodeType() const;
    virtual void setIntegrationPointsAndWeights();

    vector< Simplex * >simplices;

public:
    RigidBodyContact(const unsigned dim);
    virtual Matrix giveAMatrix(Point a, Point x) const;
    ~RigidBodyContact() {};
    void readFromLine(istringstream &iss, NodeContainer *fullnodes, MaterialContainer *fullmatrs);
    void init();
    virtual Matrix giveBMatrix(const Point *x) const;
    virtual Matrix giveHMatrix(const Point *x) const;
    vector< Node * >giveVertices() const { return vert; };
    unsigned giveNumOfVertices() const { return vert.size(); };
    double giveLength() const { return length; }
    double giveArea() const { return area; }
    virtual Vector giveContactStrainNT() const;
    virtual Vector giveContactStrainXYZ() const;
    virtual Vector giveContactStressXYZ();
    virtual Vector transformToLocal(const Vector &DoFs) const;
    virtual Vector transformToGlobal(const Vector &DoFs) const;
    Vector transformVectorToXYZ(Vector &result) const;

    virtual void giveValues(string code, Vector &result) const;
    Vector giveVectorToNode(const unsigned &node_i, const unsigned &ip_id) const;
    Point giveNormal() const { return normal; };
    Point giveT1() const { return t1; };
    Point giveT2() const { return t2; };
    double giveVolumeAssociatedWithNode(unsigned nodenum) const;
    virtual Vector giveStrain(unsigned i, const Vector &DoFs);
    virtual Matrix giveStiffnessMatrix(string matrixType) const;
    virtual Matrix giveDampingMatrix() const;
    virtual Matrix giveMassMatrix() const;
    virtual Vector giveInternalForces(const Vector &DoFs, bool frozen, double timeStep);
    virtual Vector integrateLoad(BodyLoad *vl, double time) const;
    virtual Vector integrateInternalSources();
    double givePerimeter()const { return perimeter; };

    virtual void extrapolateIPValuesToNodes(string code, vector< Vector > &result, Vector &weights) const;
    virtual bool isPointInside(Point *xn, const Point *x) const { ( void ) xn; ( void ) x; return false; }; //TODO: discrete elements does not interpolate
};

// //////////////////////////////////////////////////////////
// //////////////////////////////////////////////////////////
// // RBSN BOUNDARY ELEMENT
// TODO this will be needed in dynamics
// class RigidBodyBoundary : public RigidBodyContact
// {
// protected:
//     virtual void checkNodeType() const;
// public:
//     RigidBodyBoundary(const unsigned dim);
//     ~RigidBodyBoundary() {};
//     virtual Matrix giveStiffnessMatrix(string matrixType) const;
//     virtual Matrix giveDampingMatrix() const;
//     virtual Vector giveInternalForces(const Vector &DoFs, bool frozen, double timeStep);
//     virtual Vector giveStrain(unsigned i, const Vector &DoFs);
// };

//////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////
// COUPLED RBSN ELEMENT
class Transp1D; //forward declaration
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
// COUPLED RBSN BOUNDARY ELEMENT
// takes boundary pressure (or similar) and transfer into internal force according to biot coeff
class RigidBodyBoundaryCoupled : public RigidBodyContactCoupled
{
protected:
    virtual void checkNodeType() const;
public:
    RigidBodyBoundaryCoupled(const unsigned dim);
    ~RigidBodyBoundaryCoupled() {};

    virtual Matrix giveBMatrix(const Point *x) const;
    virtual Matrix giveHMatrix(const Point *x) const;

    void init();

    virtual Matrix giveStiffnessMatrix(string matrixType) const { ( void ) matrixType; return Matrix( ( this->ndim - 1 ) * 3, ( this->ndim - 1 ) * 3 ); };
    virtual Matrix giveDampingMatrix() const { return Matrix( ( this->ndim - 1 ) * 3, ( this->ndim - 1 ) * 3 ); };
    // virtual Vector giveInternalForces(const Vector &DoFs, bool frozen, double timeStep);
    virtual Vector giveStrain(unsigned i, const Vector &DoFs);
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
    virtual Matrix giveAMatrix(Point a, Point x) const;
    virtual Vector giveContactStrainNT(const Vector &DoFs) const;
    virtual Matrix giveBMatrix(const Point *x) const;
    virtual Matrix giveHMatrix(const Point *x) const;
};


//////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////
// 1D TRANSPORT ELEMENT
class Transp1D : public TransportElement
{
protected:
    vector< Node * >vert;
    bool bound;
    Point normal;
    double length, area;
    bool BolanderCapacityMatrix;

    virtual void checkNodeType() const;
    virtual void setIntegrationPointsAndWeights();

public:
    Transp1D(const unsigned dim);
    ~Transp1D() {};
    void init();
    double giveArea() const { return area; }
    Point giveNormal() const { return normal; }
    double giveLength() const { return length; }
    double giveVolumeAssociatedWithNode(unsigned nodenum) const;
    void readFromLine(istringstream &iss, NodeContainer *fullnodes, MaterialContainer *fullmatrs);
    virtual Matrix giveBMatrix(const Point *x) const;
    virtual Matrix giveHMatrix(const Point *x) const;
    virtual Matrix giveDampingMatrix() const;
    virtual Matrix giveStiffnessMatrix(string matrixType) const;
    virtual Vector giveInternalForces(const Vector &DoFs, bool frozen, double timeStep);
    virtual Vector integrateLoad(BodyLoad *vl, double time) const;
    virtual Vector integrateInternalSources();
    virtual Vector giveStrain(unsigned i, const Vector &DoFs);
    vector< Node * >giveVertices() const { return vert; };
    virtual bool isPointInside(Point *xn, const Point *x) const { ( void ) xn; ( void ) x; return false; }; //TODO: discrete elements does not interpolate
};

//////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////
// COUPLED 1D TRANSPORT ELEMENT
class Transp1DCoupled : public Transp1D
{
private:
    vector< RigidBodyContact * >friends; //mechanical elements involved in computation
    vector< double >friendsweight;  //weight of mechanical elements

    void findFriendsFromSimplices();

public:
    Transp1DCoupled(const unsigned dim) : Transp1D(dim) { name = "Transp1DCoupled"; solution_order = 1; }; //coupled elements must be solved after all RBSN elements
    ~Transp1DCoupled() {};
    virtual void giveValues(string code, Vector &result) const;
    void init();
    virtual Vector giveStrain(unsigned i, const Vector &DoFs);
    void addNewFriend(RigidBodyContact *f, double weight);
    unsigned giveNumOfFriends() const { return friends.size(); };
    virtual void collectInformationsFromNeigborhood();
};

#endif  /* _ELEMENT_DISCRETE_H */
