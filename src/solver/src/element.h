#ifndef _ELEMENT_H
#define _ELEMENT_H

#include "linear_algebra.h"
#include "node_container.h"
#include "material_container.h"


//////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////
// BASIC ELEMENT - MASTER CLASS
class Element
{
private:

protected:
    unsigned ndim;
    vector< Node * >nodes;
    string name;
    Material *mat;
    vector< Point >ip_locs;
    vector< double >ip_weights;
    vector< MaterialStatus * >stats;
    vector< unsigned >DoFids;
public:
    Element() { name = "basic element"; }
    virtual ~Element();
    virtual void readFromLine(istringstream &iss, NodeContainer *fullnodes, MaterialContainer *fullmatrs) { ( void ) iss; ( void ) fullnodes; ( void ) fullmatrs; };
    virtual void init();
    void initMaterialStatuses();
    void updateMaterialStatuses();
    virtual Matrix giveSteadyStateMatrix(string matrixType) const = 0;
    vector< unsigned >giveDoFs() { return DoFids; };
    virtual Vector giveInternalForces(const Vector &DoFs, bool frozen) const = 0;
    virtual double giveValue(string code) const;
    virtual string giveName() const { return name; }
    virtual size_t giveIPNum() const { return ip_locs.size(); };
    virtual double giveIPValue(string code, unsigned ipnum) const;
    virtual vector< Node * >giveNodes() const { return nodes; }
    virtual Material *giveMaterial() const { return mat; }
};


//////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////
// GEOMETRICAL ELEMENT - JUST TO REPRESENT GEOMETRICAL ENTITIES
class GeometricalElement : public Element
{
private:

public:
    GeometricalElement() {mat = nullptr;}
    ~GeometricalElement() {};
    Matrix giveSteadyStateMatrix(string matrixType) const {(void)matrixType; return Matrix(0, 0);};
    Vector giveInternalForces(const Vector &DoFs, bool frozen) const {(void)DoFs; (void)frozen; return Vector(0);};
    double giveIPValue(string code, unsigned ipnum){(void) code; (void) ipnum; return 0;}
};

//////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////
// TRANSPORT ELEMENT
class TransportElement : public Element
{
private:

public:
    TransportElement() {}
    ~TransportElement() {};
    virtual Matrix giveConductivityMatrix(string matrixType) const = 0;
    virtual Matrix giveCapacityMatrix() const = 0;
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
    virtual double giveValue(string code) const;
    virtual double giveIPValue(string code, unsigned ipnum) const;
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
private:
    vector< Node * >vert;
    bool bound;
    Point normal;
    double length, area;
    bool reducedCapacityMatrix;
public:
    Transp1D(const unsigned dim);
    ~Transp1D() {};
    void init();
    void readFromLine(istringstream &iss, NodeContainer *fullnodes, MaterialContainer *fullmatrs);
    Matrix giveConductivityMatrix(string matrixType) const;
    Matrix giveCapacityMatrix() const;
    Matrix giveSteadyStateMatrix(string matrixType) const { return giveConductivityMatrix(matrixType); };
    virtual Vector giveInternalForces(const Vector &DoFs, bool frozen) const;
};

#endif  /* _ELEMENT_STRUCT_H */
