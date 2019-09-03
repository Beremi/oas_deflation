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
    vector< MaterialStatus * >stats;
    vector< unsigned >DoFids;
public:
    Element() { name = "basic element"; }
    virtual ~Element();
    virtual void readFromLine(istringstream &iss, NodeContainer *fullnodes, MaterialContainer *fullmatrs) {};
    virtual void init();
    void initMaterialStatuses();
    void updateMaterialStatuses();
    virtual Matrix giveSteadyStateMatrix(string matrixType) const = 0;
    virtual Matrix giveTransientMatrix() const = 0;
    vector< unsigned >giveDoFs() { return DoFids; };
    virtual Vector giveInternalForces(const Vector &DoFs) const = 0;
    virtual double giveValue(string code) const;
    virtual string giveName() const { return name; }
    virtual unsigned giveIPNum() const { return ip_locs.size(); };
    virtual double giveIPValue(string code, unsigned ipnum) const;
    virtual vector < Node * > giveNodes() const { return nodes;}
};


//////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////
// TRANSPORT ELEMENT
class transportElement : public Element
{
private:

public:
    transportElement() {}
    ~transportElement() {};
    virtual Matrix giveConductivityMatrix(string matrixType) const = 0;
    virtual Matrix giveCapacityMatrix() const = 0;
};


//////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////
// MECHANICAL ELEMENT
class mechanicalElement : public Element
{
private:

public:
    mechanicalElement() {}
    ~mechanicalElement() {};
    virtual Matrix giveStiffnessMatrix(string matrixType) const = 0;
    virtual Matrix giveInertiaMatrix() const = 0;
};

//////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////
// RBSN ELEMENT
class RigidBodyContact : public mechanicalElement
{
private:
    vector< Node * >vert;
    double length, area;
    Point normal;
    vector< Point > tangs;
    Matrix RB; //R*B
public:
    RigidBodyContact(const unsigned dim);
    ~RigidBodyContact() {};
    void readFromLine(istringstream &iss, NodeContainer *fullnodes, MaterialContainer *fullmatrs);
    void init();
    vector< Node * > giveVertices() const { return vert; };
    Matrix giveStiffnessMatrix(string matrixType) const;
    Matrix giveInertiaMatrix() const;
    Matrix giveSteadyStateMatrix(string matrixType) const { return giveStiffnessMatrix(matrixType); };
    Matrix giveTransientMatrix() const { return giveInertiaMatrix(); };
    //Matrix giveBMatrix(Point x) const {return B;};
    //Matrix giveRMatrix() const {return R;};
    Matrix giveRBMatrix() const { return RB; };
    Matrix giveAMatrix(Point a, Point x) const;
    double giveLength() const { return length; }
    double giveArea() const { return area ; }
    virtual Vector giveInternalForces(const Vector &DoFs) const;
    virtual Vector giveContactStrainNT(const Vector &DoFs) const;
    virtual double giveValue(string code) const;
    virtual double giveIPValue(string code, unsigned ipnum) const;
};

//////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////
// 1D TRANSPORT ELEMENT
class Transp1D : public transportElement
{
private:
    vector< Node * >vert;
    bool bound;
    double length, area;
public:
    Transp1D(const unsigned dim);
    ~Transp1D() {};
    void init();
    void readFromLine(istringstream &iss, NodeContainer *fullnodes, MaterialContainer *fullmatrs);
    Matrix giveConductivityMatrix(string matrixType) const;
    Matrix giveCapacityMatrix() const;
    Matrix giveSteadyStateMatrix(string matrixType) const { return giveConductivityMatrix(matrixType); };
    Matrix giveTransientMatrix() const { return giveCapacityMatrix(); };
    virtual Vector giveInternalForces(const Vector &DoFs) const;
};

#endif  /* _ELEMENT_STRUCT_H */
