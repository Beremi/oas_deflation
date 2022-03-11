#ifndef _ELEMENT_H
#define _ELEMENT_H

#include "linalg.h"
#include "linalg.h"
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
    double volume;
    std :: vector< Node * >nodes;
    std :: string name;
    Material *mat;
    //vector< Point >ip_locs;
    //vector< double >ip_weights;
    std :: vector< MyMatrix >Bs;     //stored B matrices
    std :: vector< MyMatrix >Hs;     //stored H matrices
    std :: vector< MaterialStatus * >stats;
    std :: vector< unsigned >DoFids;
    unsigned outDoFs; // for coupled elements, number of input DoFs might be different from number of output DoFs.
    virtual void setIntegrationPointsAndWeights();
    virtual void initIntegration();

    ShapeFunc *shafunc;
    IntegrationType *inttype;
    unsigned numOfNodes;

    unsigned vtk_cell_type = 0; //integer detrmining type of cell for VTK plotting,
    //vetrex 1, line 3, triangle 5, polygon 7, quad 9, tetra 10, brick 12, quadratic_triangle 22, quadratic_tetra 24, quadratic_brick 25

public:
    Element() { name = "basic element"; solution_order = 0; volume = 0; }
    virtual ~Element();
    void setID(unsigned i) { idx = i; };
    unsigned giveID() const { return idx; };
    virtual void readFromLine(std :: istringstream &iss, NodeContainer *fullnodes, MaterialContainer *fullmatrs);
    // virtual std :: string giveLineToSave(NodeContainer * nodes) const;
    virtual void init();
    void initMaterialStatuses();
    void updateMaterialStatuses();
    void resetMaterialStatuses();
    virtual MyMatrix giveStiffnessMatrix(std :: string matrixType) const;
    virtual MyMatrix giveDampingMatrix() const;
    virtual MyMatrix giveMassMatrix() const;
    virtual MyVector giveInternalForces(const MyVector &DoFs, bool frozen, double timeStep);
    std :: vector< unsigned >giveDoFs() const { return DoFids; };
    std :: vector< unsigned >giveDoFsInDirection(unsigned dir) const;
    unsigned giveNumOutDoFs() const { return outDoFs; };
    virtual void giveValues(std :: string code, MyVector &result) const;
    std :: string giveName() const { return name; }
    size_t giveNumIP() const { return inttype->giveNumIP(); };
    Point giveIPLoc(unsigned k) const { return inttype->giveIPLocation(k); };
    virtual void giveIPValues(std :: string code, unsigned ipnum, MyVector &result) const;
    std :: vector< Node * >giveNodes() const { return nodes; }
    Node *giveNode(unsigned k) const { return nodes [ k ]; }
    size_t giveNumOfNodes() const { return nodes.size(); }
    Material *giveMaterial() const { return mat; }
    std :: vector< MaterialStatus * >giveMaterialStats() const { return stats; };
    MaterialStatus *giveMatStatus(unsigned ipnum) { return stats [ ipnum ]; };
    virtual void findElementFriends(ElementContainer *elemcont) { ( void ) elemcont; }
    unsigned giveSolutionOrder() const { return solution_order; }
    virtual MyMatrix giveBMatrix(const Point *x) const { ( void ) x; return MyMatrix(0, 0); };
    MyMatrix giveStoredBMatrix(unsigned i) { return Bs [ i ]; };
    virtual MyMatrix giveHMatrix(const Point *x) const { ( void ) x; return MyMatrix(0, 0); };
    MyMatrix giveStoredHMatrix(unsigned i) { return Hs [ i ]; };
    virtual MyVector giveStrain(const Point *x, const MyVector &DoFs) const { return giveBMatrix(x) * DoFs; };
    virtual MyVector giveStrain(unsigned i, const MyVector &DoFs) { return Bs [ i ] * DoFs; };
    unsigned giveDimension() const { return ndim; }
    virtual MyVector integrateLoad(BodyLoad *vl, double time) const;
    unsigned giveVTKCellType() const { return vtk_cell_type; };
    virtual void changeMaterial(Material *newmat);
    virtual MyVector integrateInternalSources();
    //virtual void shapeF(const Point *x, MyVector &phi) const { ( void ) x; ( void ) phi; };
    //virtual double shapeFGrad(const Point *x, MyMatrix &phiGrad) const { ( void ) x; ( void ) phiGrad; return 0; };
    virtual bool giveGlobalCoords(Point *x, const Point *xn) const;
    virtual MyVector giveMasterVariables(const Point *x, const MyVector &DoFs) const { return giveHMatrix(x) * DoFs; };
    MyVector giveElemDoFsFromFullDoFs(const MyVector &FullDoFs) const;
    double giveVolume() const { return volume; };

    virtual void extrapolateIPValuesToNodes(std :: string code, std :: vector< MyVector > &result, MyVector &weights) const;

    virtual void collectInformationsFromNeigborhood() {};
    virtual bool isPointInside(Point *xn, const Point *x) const;
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
    double giveIPValue(std :: string code, unsigned ipnum) const { ( void ) code; ( void ) ipnum; return 0; }
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
