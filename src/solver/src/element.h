#ifndef _ELEMENT_H
#define	_ELEMENT_H

#include "linear_algebra.h"
#include "node_container.h"
#include "material_container.h"


//////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////
// BASIC ELEMENT - MASTER CLASS
class Element {
private:

protected:
    unsigned ndim;
    vector<Node *> nodes;
    string name;
    Material *mat;
    vector<Point> ip_locs;
    vector<MaterialStatus *> stats;
    vector<unsigned> DoFids;
public:
  Element(){name = "basic element";}
  virtual ~Element();
  virtual void readFromLine(istringstream &iss, NodeContainer *fullnodes, MaterialContainer *fullmatrs){};
  virtual void init();
  virtual Matrix giveSteadyStateMatrix() const{};  
  virtual Matrix giveTransientMatrix() const{};
  vector<unsigned> giveDoFs(){return DoFids;};
};


//////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////
// TRANSPORT ELEMENT
class transportElement: public Element {
private:

public:
    transportElement (){}
    ~transportElement (){};
    virtual Matrix giveConductivityMatrix()const{};
    virtual Matrix giveCapacityMatrix()const{};
};


//////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////
// MECHANICAL ELEMENT
class mechanicalElement: public Element {
private:

public:
    mechanicalElement (){}
    ~mechanicalElement (){};
    virtual Matrix giveStiffnessMatrix()const{};
    virtual Matrix giveInertiaMatrix()const{};
};

//////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////
// RBSN ELEMENT
class RigidBodyContact: public mechanicalElement {
private:
    vector<Node *> vert;
    double length, area;
    Point normal;
    Matrix B, R;
public:
    RigidBodyContact (const unsigned dim);
    ~RigidBodyContact (){};
    void readFromLine(istringstream &iss, NodeContainer *fullnodes, MaterialContainer *fullmatrs);
    void init();
    Matrix giveStiffnessMatrix() const;
    Matrix giveInertiaMatrix() const;
    Matrix giveSteadyStateMatrix() const {return giveStiffnessMatrix();};  
    Matrix giveTransientMatrix() const {return giveInertiaMatrix();};
    Matrix giveBMatrix(Point x) const {return B;};
    Matrix giveRMatrix() const {return R;};
    Matrix giveAMatrix(Point a, Point x) const;
};

//////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////
// 1D TRANSPORT ELEMENT
class Transp1D: public transportElement {
private:
    vector<Node *> vert;
    bool bound;
    double length, area;
public:
    Transp1D (const unsigned dim);
    ~Transp1D (){};
    void init();
    void readFromLine(istringstream &iss, NodeContainer *fullnodes, MaterialContainer *fullmatrs);
    Matrix giveConductivityMatrix() const;
    Matrix giveCapacityMatrix() const;
    Matrix giveSteadyStateMatrix() const {return giveConductivityMatrix();};  
    Matrix giveTransientMatrix() const {return giveCapacityMatrix();};
};

#endif	/* _ELEMENT_STRUCT_H */
