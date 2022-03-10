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
    Point normal;
    Point t1, t2;
    MyMatrix R;

    MyMatrix giveRMatrix() const { return R; };
    virtual void checkNodeType() const;
    virtual void setIntegrationPointsAndWeights();

    vector< Simplex * >simplices;

public:
    RigidBodyContact(const unsigned dim);
    virtual MyMatrix giveAMatrix(unsigned v, Point x) const;
    ~RigidBodyContact() {};
    void readFromLine(istringstream &iss, NodeContainer *fullnodes, MaterialContainer *fullmatrs);
    void init();
    virtual MyMatrix giveBMatrix(const Point *x) const;
    virtual MyMatrix giveHMatrix(const Point *x) const;
    vector< Node * >giveVertices() const { return vert; };
    unsigned giveNumOfVertices() const { return vert.size(); };
    double giveLength() const { return length; }
    double giveArea() const { return area; }
    virtual MyVector giveContactStrainNT() const;
    virtual MyVector giveContactStrainXYZ() const;
    virtual MyVector giveContactStressXYZ();
    virtual MyVector transformToLocal(const MyVector &DoFs) const;
    virtual MyVector transformToGlobal(const MyVector &DoFs) const;
    MyVector transformVectorToXYZ(MyVector &result) const;

    virtual void giveValues(string code, MyVector &result) const;
    MyVector giveVectorToNode(const unsigned &node_i, const unsigned &ip_id) const;
    Point giveNormal() const { return normal; };
    Point giveT1() const { return t1; };
    Point giveT2() const { return t2; };
    double giveVolumeAssociatedWithNode(unsigned nodenum) const;
    virtual MyVector giveStrain(unsigned i, const MyVector &DoFs);
    virtual MyMatrix giveStiffnessMatrix(string matrixType) const;
    virtual MyMatrix giveDampingMatrix() const;
    virtual MyMatrix giveMassMatrix() const;
    virtual MyVector giveInternalForces(const MyVector &DoFs, bool frozen, double timeStep);
    virtual MyVector integrateLoad(BodyLoad *vl, double time) const;
    virtual MyVector integrateInternalSources();
    double givePerimeter()const { return perimeter; };

    virtual void extrapolateIPValuesToNodes(string code, vector< MyVector > &result, MyVector &weights) const;
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
//     virtual MyMatrix giveStiffnessMatrix(string matrixType) const;
//     virtual MyMatrix giveDampingMatrix() const;
//     virtual MyVector giveInternalForces(const MyVector &DoFs, bool frozen, double timeStep);
//     virtual MyVector giveStrain(unsigned i, const MyVector &DoFs);
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
    virtual MyVector giveStrain(unsigned i, const MyVector &DoFs);
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

    virtual MyMatrix giveBMatrix(const Point *x) const;
    virtual MyMatrix giveHMatrix(const Point *x) const;

    void init();

    virtual MyMatrix giveStiffnessMatrix(string matrixType) const { ( void ) matrixType; return MyMatrix::Zero( ( this->ndim - 1 ) * 3, ( this->ndim - 1 ) * 3 ); };
    virtual MyMatrix giveDampingMatrix() const { return MyMatrix :: Zero( ( this->ndim - 1 ) * 3, ( this->ndim - 1 ) * 3 ); };
    // virtual MyVector giveInternalForces(const MyVector &DoFs, bool frozen, double timeStep);
    virtual MyVector giveStrain(unsigned i, const MyVector &DoFs);
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
    virtual MyMatrix giveAMatrix(unsigned v, Point x) const;
    virtual MyVector giveContactStrainNT(const MyVector &DoFs) const;
    virtual MyMatrix giveBMatrix(const Point *x) const;
    virtual MyMatrix giveHMatrix(const Point *x) const;
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
    virtual MyMatrix giveBMatrix(const Point *x) const;
    virtual MyMatrix giveHMatrix(const Point *x) const;
    virtual MyMatrix giveDampingMatrix() const;
    virtual MyMatrix giveStiffnessMatrix(string matrixType) const;
    virtual MyVector giveInternalForces(const MyVector &DoFs, bool frozen, double timeStep);
    virtual MyVector integrateLoad(BodyLoad *vl, double time) const;
    virtual MyVector integrateInternalSources();
    virtual MyVector giveStrain(unsigned i, const MyVector &DoFs);
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
    virtual void giveValues(string code, MyVector &result) const;
    void init();
    virtual MyVector giveStrain(unsigned i, const MyVector &DoFs);
    void addNewFriend(RigidBodyContact *f, double weight);
    unsigned giveNumOfFriends() const { return friends.size(); };
    virtual void collectInformationsFromNeigborhood();
};

#endif  /* _ELEMENT_DISCRETE_H */
