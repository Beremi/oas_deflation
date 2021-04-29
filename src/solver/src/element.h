#ifndef _ELEMENT_H
#define _ELEMENT_H

#include "linear_algebra.h"
#include "node_container.h"
#include "material_container.h"
#include "shape_functions.h"
#include "integration.h"

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
    virtual void setIntegrationPointsAndWeights();

    ShapeFunc *shafunc;
    IntegrationType *inttype;
    unsigned numOfNodes;

    unsigned vtk_cell_type = 0; //integer detrmining type of cell for VTK plotting,
    //vetrex 1, line 3, triangle 5, polygon 7, quad 9, tetra 10, brick 12, quadratic_triangle 22, quadratic_tetra 24, quadratic_brick 25

public:
    Element() { name = "basic element"; solution_order = 0; }
    virtual ~Element();
    void setID(unsigned i) { idx = i; };
    unsigned giveID() const { return idx; };
    virtual void readFromLine(istringstream &iss, NodeContainer *fullnodes, MaterialContainer *fullmatrs);
    // virtual std :: string giveLineToSave(NodeContainer * nodes) const;
    virtual void init();
    void initMaterialStatuses();
    void updateMaterialStatuses();
    virtual Matrix giveStiffnessMatrix(string matrixType) const;
    virtual Matrix giveDampingMatrix() const;
    virtual Matrix giveMassMatrix() const;
    virtual Vector giveInternalForces(const Vector &DoFs, bool frozen, double timeStep);
    vector< unsigned >giveDoFs() const { return DoFids; };
    vector< unsigned >giveDoFsInDirection(unsigned dir) const;
    unsigned giveNumOutDoFs() const { return outDoFs; };
    virtual double giveValue(string code) const;
    string giveName() const { return name; }
    size_t giveNumIP() const { return inttype->giveNumIP(); };
    Point giveIPLoc(unsigned k) const { return inttype->giveIPLocation(k); };
    virtual double giveIPValue(string code, unsigned ipnum) const;
    MaterialStatus *giveMatStatus(unsigned ipnum) { return stats [ ipnum ]; };
    vector< Node * >giveNodes() const { return nodes; }
    Node *giveNode(unsigned k) const { return nodes [ k ]; }
    Material *giveMaterial() const { return mat; }
    vector< MaterialStatus * >giveMaterialStats() const { return stats; };
    virtual void findElementFriends(ElementContainer *elemcont) { ( void ) elemcont; }
    unsigned giveSolutionOrder() const { return solution_order; }
    virtual Matrix giveBMatrix(const Point *x) const { return Matrix(0, 0); };
    Matrix *giveBMatrix(unsigned i) { return & Bs [ i ]; };
    virtual Matrix giveHMatrix(const Point *x) const { return Matrix(0, 0); };
    virtual Vector giveStrain(const Point *x, const Vector &DoFs) { return giveBMatrix(x) * DoFs; };
    virtual Vector giveStrain(unsigned i, const Vector &DoFs) { return Bs [ i ] * DoFs; };
    unsigned giveDimension() const { return ndim; }
    virtual vector< double >integrateLoad(BodyLoad *vl, double time) const;
    unsigned giveVTKCellType() const { return vtk_cell_type; };
    virtual void changeMaterial(Material *newmat);

    virtual void shapeF(const Point *x, Vector &phi) const { ( void ) x; ( void ) phi; };
    virtual double shapeFGrad(const Point *x, Matrix &phiGrad) const { ( void ) x; ( void ) phiGrad; return 0; };
    
    virtual void collectInformationsFromNeigborhood(){};
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
};

#endif  /* _ELEMENT_STRUCT_H */
