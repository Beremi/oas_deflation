#ifndef _ELEMENT_H
#define _ELEMENT_H

#include "linear_algebra.h"
#include "node_container.h"
#include "material_container.h"

class ElementContainer; //forward declaration;
class BodyLoad; //forward declaration

//////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////
// BASIC ELEMENT - MASTER CLASS
class Element
{
private:

protected:
    unsigned ndim;
    unsigned idx;
    unsigned solution_order;
    vector< Node * >nodes;
    string name;
    Material *mat;
    vector< Point >ip_locs;
    vector< double >ip_weights;
    vector< Matrix >Bs;
    vector< MaterialStatus * >stats;
    vector< unsigned >DoFids;
    unsigned outDoFs; // for coupled elements, number of input DoFs might be different from number of output DoFs.
public:
    Element() { name = "basic element"; solution_order = 0; }
    virtual ~Element();
    void setID(unsigned i) { idx = i; };
    unsigned giveID() const { return idx; };
    virtual void readFromLine(istringstream &iss, NodeContainer *fullnodes, MaterialContainer *fullmatrs) { ( void ) iss; ( void ) fullnodes; ( void ) fullmatrs; };
    virtual std :: string giveLineToSave(NodeContainer * nodes) const;
    virtual void init();
    void initMaterialStatuses();
    void updateMaterialStatuses();
    virtual Matrix giveSteadyStateMatrix(string matrixType) const;
    virtual Vector giveInternalForces(const Vector &DoFs, bool frozen);
    virtual Matrix giveMassMatrix() const;
    vector< unsigned >giveDoFs() const { return DoFids; };
    vector< unsigned >giveDoFsInDirection(unsigned dir) const;
    unsigned giveNumOutDoFs() const { return outDoFs; };
    virtual double giveValue(string code) const;
    string giveName() const { return name; }
    size_t giveIPNum() const { return ip_locs.size(); };
    Point giveIPLoc(unsigned k) const { return ip_locs [ k ]; };
    virtual double giveIPValue(string code, unsigned ipnum) const;
    MaterialStatus *giveMatStatus(unsigned ipnum) { return stats [ ipnum ]; };
    vector< Node * >giveNodes() const { return nodes; }
    Node *giveNode(unsigned k) const { return nodes [ k ]; }
    Material *giveMaterial() const { return mat; }
    vector< MaterialStatus * >giveMaterialStats() const { return stats; };
    virtual void findElementFriends(ElementContainer *elemcont) { ( void ) elemcont; }
    unsigned giveSolutionOrder() const { return solution_order; }
    virtual void shapeF(const Point *x, Vector &phi) const { ( void ) x; ( void ) phi; };
    virtual double shapeFGrad(const Point *x, Matrix &phiGrad) const { ( void ) x; ( void ) phiGrad; return 0; };
    virtual Matrix giveBMatrix(const Point *x) const { return Matrix(0, 0); };
    Matrix *giveBMatrix(unsigned i) { return & Bs [ i ]; };
    virtual Matrix giveHMatrix(const Point *x) const { return Matrix(0, 0); };
    virtual Vector giveStrain(const Point *x, const Vector &DoFs) { return giveBMatrix(x) * DoFs; };
    virtual Vector giveStrain(unsigned i, const Vector &DoFs) { return Bs [ i ] * DoFs; };
    unsigned giveDimension() const { return ndim; }
    virtual void setIntegrationPointsAndWeights() {};
    virtual vector< double >integrateLoad(BodyLoad *vl, double time) const;
};


//////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////
// GEOMETRICAL ELEMENT - JUST TO REPRESENT GEOMETRICAL ENTITIES
class GeometricalElement : virtual public Element
{
protected:

public:
    GeometricalElement() { mat = nullptr; }
    ~GeometricalElement() {};
    double giveIPValue(string code, unsigned ipnum) const { ( void ) code; ( void ) ipnum; return 0; }
};

//////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////
// TRANSPORT ELEMENT
class TransportElement : public Element
{
protected:

public:
    TransportElement() {}
    ~TransportElement() {};
    virtual Matrix giveConductivityMatrix(string matrixType) const { return giveSteadyStateMatrix(matrixType); };
    virtual Matrix giveCapacityMatrix() const { return giveMassMatrix(); };
};


//////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////
// MECHANICAL ELEMENT
class MechanicalElement : public Element
{
protected:

public:
    MechanicalElement() {}
    ~MechanicalElement() {};
    virtual Matrix giveStiffnessMatrix(string matrixType) const { return giveSteadyStateMatrix(matrixType); };
};

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
public:
    RigidBodyContact(const unsigned dim);
    virtual Matrix giveAMatrix(Point a, Point x) const;
    ~RigidBodyContact() {};
    void readFromLine(istringstream &iss, NodeContainer *fullnodes, MaterialContainer *fullmatrs);
    void init();
    virtual Matrix giveBMatrix(const Point *x) const;
    virtual Matrix giveHMatrix(const Point *x) const;
    virtual void setIntegrationPointsAndWeights();
    vector< Node * >giveVertices() const { return vert; };
    double giveLength() const { return length; }
    double giveArea() const { return area; }
    virtual Vector giveContactStrainNT(const Vector &DoFs) const;
    virtual Vector giveContactStrainXYZ(const Vector &DoFs) const;
    virtual Vector giveContactStressXYZ(const Vector &DoFs);
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
public:
    Transp1D(const unsigned dim);
    ~Transp1D() {};
    void init();
    virtual void setIntegrationPointsAndWeights();
    double giveArea() const { return area; }
    Point giveNormal() const { return normal; }
    double giveLength() const { return length; }
    double giveVolume(unsigned nodenum) const;
    double giveVolume() const;
    void readFromLine(istringstream &iss, NodeContainer *fullnodes, MaterialContainer *fullmatrs);
    virtual Matrix giveBMatrix(const Point *x) const;
    virtual Matrix giveHMatrix(const Point *x) const;
    virtual Matrix giveCapacityMatrix() const;
    virtual vector< double >integrateLoad(BodyLoad *vl, double time) const;
    virtual Vector giveStrain(unsigned i, const Vector &DoFs);
};

//////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////
// COUPLED 1D TRANSPORT ELEMENT
class Transp1DCoupled : public Transp1D
{
private:
    vector< RigidBodyContact * >friends; //mechanical elements involved in computation
    vector< double >friendsweight;  //weight of mechanical elements

    void findFriends2D(ElementContainer *elemcont);
    void findFriends3D(ElementContainer *elemcont);

public:
    Transp1DCoupled(const unsigned dim) : Transp1D(dim) { name = "Transp1DCoupled"; solution_order = 1; }; //coupled elements must be solved after all RBSN elements
    void findElementFriends(ElementContainer *elemcont);
    ~Transp1DCoupled() {};
    virtual double giveValue(string code) const;
    virtual double giveIPValue(string code, unsigned ipnum) const;
    void init();
    virtual Vector giveStrain(unsigned i, const Vector &DoFs);
};

//////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////
// 2D QUADRILATERAL TRANSPORT ELEMENT
class TranspQuad : public TransportElement
{
protected:

public:
    TranspQuad();
    ~TranspQuad() {};
    void setIntegrationPointsAndWeights();
    void readFromLine(istringstream &iss, NodeContainer *fullnodes, MaterialContainer *fullmatrs);
    virtual void shapeF(const Point *x, Vector &phi) const;
    virtual double shapeFGrad(const Point *x, Matrix &phiGrad) const;
    virtual Matrix giveBMatrix(const Point *x) const;
    virtual Matrix giveHMatrix(const Point *x) const;
    virtual Vector giveStrain(unsigned i, const Vector &DoFs);
};

//////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////
// 2D QUADRILATERAL MECHANICAL ELEMENT
class MechanicalQuad : public MechanicalElement
{
protected:

public:
    MechanicalQuad();
    ~MechanicalQuad() {};
    void setIntegrationPointsAndWeights();
    void readFromLine(istringstream &iss, NodeContainer *fullnodes, MaterialContainer *fullmatrs);
    virtual void shapeF(const Point *x, Vector &phi) const;
    virtual double shapeFGrad(const Point *x, Matrix &phiGrad) const;
    virtual Matrix giveBMatrix(const Point *x) const;
    virtual Matrix giveHMatrix(const Point *x) const;
};


//////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////
// 2D QUADRILATERAL COSSERAT MECHANICAL ELEMENT
class CosseratQuad : public MechanicalQuad
{
protected:

public:
    CosseratQuad();
    ~CosseratQuad() {};
    virtual Matrix giveBMatrix(const Point *x) const;
    virtual Matrix giveHMatrix(const Point *x) const;
};
#endif  /* _ELEMENT_STRUCT_H */
