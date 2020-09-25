#ifndef _ELEMENT_H
#define _ELEMENT_H

#include "linear_algebra.h"
#include "node_container.h"
#include "material_container.h"

class ElementContainer; //forward declaration;

//////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////
// BASIC ELEMENT - MASTER CLASS
class Element
{
private:

protected:
    unsigned ndim;
    unsigned solution_order;
    vector< Node * >nodes;
    string name;
    Material *mat;
    vector< Point >ip_locs;
    vector< double >ip_weights;
    vector< MaterialStatus * >stats;
    vector< unsigned >DoFids;
    unsigned outDoFs; // for coupled elements, number of input DoFs might be different from number of output DoFs.
public:
    Element() { name = "basic element"; solution_order = 0; }
    virtual ~Element();
    virtual void readFromLine(istringstream &iss, NodeContainer *fullnodes, MaterialContainer *fullmatrs) { ( void ) iss; ( void ) fullnodes; ( void ) fullmatrs; };
    virtual void init();
    void initMaterialStatuses();
    void updateMaterialStatuses();
    virtual Matrix giveSteadyStateMatrix(string matrixType) const = 0;
    vector< unsigned >giveDoFs() { return DoFids; };
    unsigned giveNumOutDoFs() const { return outDoFs; };
    virtual Vector giveInternalForces(const Vector &DoFs, bool frozen) const = 0;
    virtual double giveValue(string code) const;
    string giveName() const { return name; }
    size_t giveIPNum() const { return ip_locs.size(); };
    virtual double giveIPValue(string code, unsigned ipnum) const;
    vector< Node * >giveNodes() const { return nodes; }
    Node *giveNode(unsigned k) const { return nodes [ k ]; }
    Material *giveMaterial() const { return mat; }
    vector< MaterialStatus * > giveMaterialStats() const { return stats; };
    virtual void findElementFriends(ElementContainer *elemcont) { ( void ) elemcont; }
    unsigned giveSolutionOrder() const { return solution_order; }
    virtual void shapeF(const Point *x, Vector &phi) const { (void) x; (void) phi; };
    virtual double shapeFGrad(const Point *x, Matrix &phiGrad) const { (void) x; (void) phiGrad; return 0;};
    virtual Matrix giveBMatrix(const Point *x) const {return Matrix(0,0);};
    virtual Vector giveStrain(const Point *x, const Vector &DoFs) {return giveBMatrix(x)*DoFs;};
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
    Matrix giveSteadyStateMatrix(string matrixType) const { ( void ) matrixType; return Matrix(0, 0); };
    Vector giveInternalForces(const Vector &DoFs, bool frozen) const { ( void ) DoFs; ( void ) frozen; return Vector(0); };
    double giveIPValue(string code, unsigned ipnum) { ( void ) code; ( void ) ipnum; return 0; }
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
    virtual Matrix giveConductivityMatrix(string matrixType) const = 0;
    virtual Matrix giveCapacityMatrix() const = 0;
    Matrix giveSteadyStateMatrix(string matrixType) const { return giveConductivityMatrix(matrixType); };
};


//////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////
// MECHANICAL ELEMENT
class MechanicalElement : public Element
{
protected:
    Matrix GeomM; //R*B
public:
    MechanicalElement() {}
    ~MechanicalElement() {};
    virtual Matrix giveStiffnessMatrix(string matrixType) const = 0;
    virtual Matrix giveMassMatrix() const = 0;
    Matrix giveSteadyStateMatrix(string matrixType) const { return giveStiffnessMatrix(matrixType); };
    Matrix giveGeomMMatrix() const { return GeomM; };
    virtual Vector giveInternalForces(const Vector &DoFs, bool frozen) const = 0;
};

//////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////
// RBSN ELEMENT
class RigidBodyContact : public MechanicalElement
{
protected:
    vector< Node * >vert;
    double length, area;
    Point normal;
    Matrix R;
    double tempCrackOpening; //needed for coupled analysis;

    virtual void checkNodeType() const;
    virtual Matrix giveBMatrix() const;
public:
    RigidBodyContact(const unsigned dim);
    ~RigidBodyContact() {};
    void readFromLine(istringstream &iss, NodeContainer *fullnodes, MaterialContainer *fullmatrs);
    void init();
    vector< Node * >giveVertices() const { return vert; };
    virtual Matrix giveStiffnessMatrix(string matrixType) const;
    virtual Matrix giveMassMatrix() const;
    Matrix giveRMatrix() const { return R; };
    virtual Matrix giveAMatrix(Point a, Point x) const;
    double giveLength() const { return length; }
    double giveArea() const { return area; }
    virtual Vector giveInternalForces(const Vector &DoFs, bool frozen) const;
    virtual Vector giveContactStrainNT(const Vector &DoFs) const;
    virtual Vector giveContactStrainXYZ(const Vector &DoFs) const;
    virtual double giveValue(string code) const;
    virtual double giveIPValue(string code, unsigned ipnum) const;
    double giveCrackOpening() { return tempCrackOpening; };
};

//////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////
// TRUSS ELEMENT
class Truss : public RigidBodyContact
{
protected:
    virtual void checkNodeType() const;
    virtual Matrix giveBMatrix() const;
public:
    Truss(const unsigned dim) : RigidBodyContact(dim) { name = "Truss"; };
    ~Truss() {};
    virtual Matrix giveAMatrix(Point a, Point x) const;
    virtual Matrix giveStiffnessMatrix(string matrixType) const;
    virtual Vector giveContactStrainNT(const Vector &DoFs) const;
    virtual Matrix giveMassMatrix() const;
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
    bool reducedCapacityMatrix;
public:
    Transp1D(const unsigned dim);
    ~Transp1D() {};
    void init();
    double giveArea() const { return area; }
    Point giveNormal() const { return normal; }
    double giveVolume(unsigned nodenum) const;
    double giveVolume() const;
    void readFromLine(istringstream &iss, NodeContainer *fullnodes, MaterialContainer *fullmatrs);
    Matrix giveConductivityMatrix(string matrixType) const;
    Matrix giveCapacityMatrix() const;
    virtual Vector giveInternalForces(const Vector &DoFs, bool frozen) const;
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
    void init();
    virtual Vector giveInternalForces(const Vector &DoFs, bool frozen) const;
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
    void init();
    void readFromLine(istringstream &iss, NodeContainer *fullnodes, MaterialContainer *fullmatrs);
    virtual void shapeF(const Point *x, Vector &phi) const;
    virtual double shapeFGrad(const Point *x, Matrix &phiGrad) const;
    virtual Matrix giveConductivityMatrix(string matrixType) const;
    virtual Matrix giveCapacityMatrix() const;
    virtual Vector giveInternalForces(const Vector &DoFs, bool frozen) const;
};
#endif  /* _ELEMENT_STRUCT_H */
