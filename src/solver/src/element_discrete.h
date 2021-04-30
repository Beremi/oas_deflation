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
    double length, area;
    Point normal, t1, t2;
    Matrix R;
    double tempCrackOpening; //needed for coupled analysis;

    Matrix giveRMatrix() const { return R; };
    virtual void checkNodeType() const;
    virtual void setIntegrationPointsAndWeights();

    vector< Simplex * > simplices;
    double volumetricStrain;

public:
    RigidBodyContact(const unsigned dim);
    virtual Matrix giveAMatrix(Point a, Point x) const;
    ~RigidBodyContact() {};
    void readFromLine(istringstream &iss, NodeContainer *fullnodes, MaterialContainer *fullmatrs);
    void init();
    virtual Matrix giveBMatrix(const Point *x) const;
    virtual Matrix giveHMatrix(const Point *x) const;
    vector< Node * >giveVertices() const { return vert; };
    double giveLength() const { return length; }
    double giveArea() const { return area; }
    virtual Vector giveContactStrainNT() const;
    virtual Vector giveContactStrainXYZ() const;
    virtual Vector giveContactStressXYZ();
    virtual Vector transformToLocal(const Vector &DoFs) const;
    virtual Vector transformToGlobal(const Vector &DoFs) const;

    virtual double giveValue(string code) const;
    virtual double giveIPValue(string code, unsigned ipnum) const;
    double giveCrackOpening() { return tempCrackOpening; };
    Vector giveVectorToNode(const unsigned &node_i, const unsigned &ip_id) const;
    Point giveNormal() const { return normal; };
    Point giveT1() const { return t1; };
    Point giveT2() const { return t2; };
    double giveVolume(unsigned nodenum) const;
    double giveVolume() const;
    virtual Vector giveStrain(unsigned i, const Vector &DoFs);
    double giveVolumetricStrain(){return volumetricStrain;};
    virtual Matrix giveDampingMatrix() const;
};

//////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////
// COUPLED RBSN ELEMENT
class Transp1D; //forward declaration
class RigidBodyContactCoupled : public RigidBodyContact
{
protected:
    double averagePressure;
    double volumetricStrainRate;
public:
    RigidBodyContactCoupled(const unsigned dim);
    ~RigidBodyContactCoupled() {};
    virtual Vector giveStrain(unsigned i, const Vector &DoFs);
    double giveVolumetricStrainRate(){return volumetricStrain;};
    double giveAveragePressure(){return averagePressure;};
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
    double averagePressure;

    virtual void checkNodeType() const;
    virtual void setIntegrationPointsAndWeights();

public:
    Transp1D(const unsigned dim);
    ~Transp1D() {};
    void init();
    double giveArea() const { return area; }
    Point giveNormal() const { return normal; }
    double giveLength() const { return length; }
    double giveVolume(unsigned nodenum) const;
    double giveVolume() const;
    void readFromLine(istringstream &iss, NodeContainer *fullnodes, MaterialContainer *fullmatrs);
    virtual Matrix giveBMatrix(const Point *x) const;
    virtual Matrix giveHMatrix(const Point *x) const;
    virtual Matrix giveDampingMatrix() const;
    virtual vector< double >integrateLoad(BodyLoad *vl, double time) const;
    virtual Vector giveStrain(unsigned i, const Vector &DoFs);
    vector < Node* > giveVertices() const { return vert;};
    double giveAveragePressure();
};

//////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////
// COUPLED 1D TRANSPORT ELEMENT
class Transp1DCoupled : public Transp1D
{
private:
    vector< RigidBodyContact * >friends; //mechanical elements involved in computation
    vector< double >friendsweight;  //weight of mechanical elements

    double crackInNeighborhood; ///crack parameter to account for crack opening
    void findFriendsFromSimplices();

public:
    Transp1DCoupled(const unsigned dim) : Transp1D(dim) { name = "Transp1DCoupled"; solution_order = 1; }; //coupled elements must be solved after all RBSN elements
    ~Transp1DCoupled() {};
    virtual double giveValue(string code) const;
    virtual double giveIPValue(string code, unsigned ipnum) const;
    void init();
    virtual Vector giveStrain(unsigned i, const Vector &DoFs);
    void addNewFriend(RigidBodyContact * f, double weight );
    unsigned giveNumOfFriends() const {return friends.size();};
    virtual void collectInformationsFromNeigborhood();
    double giveCrackOpeningInNeigborhood();
    virtual Vector giveInternalForces(const Vector &DoFs, bool frozen, double timeStep);
};

#endif  /* _ELEMENT_DISCRETE_H */
