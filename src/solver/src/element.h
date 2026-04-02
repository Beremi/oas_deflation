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


//physicalFields Mechanics, Transport, Thermal, Humidity

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
    std :: vector< Matrix >Bs;     //stored B matrices
    std :: vector< Matrix >Hs;     //stored H matrices
    Matrix massM;           //stored mass matrix
    Matrix dampC;           //stored damping matrix
    std :: vector< MaterialStatus * >stats;
    std :: vector< unsigned >DoFids;
    unsigned outDoFs; // for coupled elements, number of input DoFs might be different from number of output DoFs.
    virtual void setIntegrationPointsAndWeights();
    virtual void initIntegration();
    virtual void computeMassMatrix();
    virtual void computeDampingMatrix();
    bool areIPLocsInNaturalCoords;

    virtual void evaluateStresses(double timeStep);
    virtual void evaluateStressesWithFrozenIntVars(double timeStep);
    

    ShapeFunc *shafunc;
    IntegrationType *inttype;
    unsigned numOfNodes;

    unsigned vtk_cell_type = 0; //integer detrmining type of cell for VTK plotting,
    //vetrex 1, line 3, triangle 5, polygon 7, quad 9, tetra 10, brick 12, quadratic_triangle 22, quadratic_tetra 24, quadratic_brick 25

    std :: vector< bool >physicalFields;
    std :: vector <unsigned> DoFsWhichChangeSignOfInternalForces; //terms in Poisson Equation that are reverted because of negative Flux
    
public:
    Element(unsigned dim);
    virtual ~Element();
    void setID(unsigned i) { idx = i; };
    unsigned giveID() const { return idx; };
    virtual void readFromLine(std :: istringstream &iss, NodeContainer *fullnodes, MaterialContainer *fullmatrs);
    // virtual std :: string giveLineToSave(NodeContainer * nodes) const;
    virtual void init();
    void initMaterialStatuses();
    void updateMaterialStatuses();
    void resetMaterialStatuses();
    virtual Matrix giveStiffnessMatrix(std :: string matrixType) const;
    virtual Matrix giveMassMatrix();
    virtual Matrix giveDampingMatrix();
    virtual Matrix giveLumpedMassMatrix();
    virtual Vector giveInternalForces();
    double giveKineticEnergy(const Vector &velocity) const;
    std :: vector< unsigned >giveDoFs() const { return DoFids; };
    unsigned giveNumDoFs() const { return DoFids.size(); };
    std :: vector< unsigned >giveDoFsInDirection(unsigned dir) const;
    unsigned giveNumOutDoFs() const { return outDoFs; };
    virtual void giveValues(std :: string code, Vector &result) const;
    std :: string giveName() const { return name; }
    size_t giveNumIP() const { return inttype->giveNumIP(); };
    virtual Point giveIPLoc(unsigned k) const;
    double giveIPWeight(unsigned k) const { return inttype->giveIPWeight(k); };
    virtual void giveIPValues(std :: string code, unsigned ipnum, Vector &result) const;
    std :: vector< Node * >giveNodes() const { return nodes; }
    Node *giveNode(unsigned k) const { return nodes [ k ]; }
    size_t giveNumOfNodes() const { return nodes.size(); }
    Material *giveMaterial() const { return mat; }
    std :: vector< MaterialStatus * >giveMaterialStats() const { return stats; };
    MaterialStatus *giveMatStatus(unsigned ipnum);
    virtual void findElementFriends(ElementContainer *elemcont) { ( void ) elemcont; }
    unsigned giveSolutionOrder() const { return solution_order; }
    virtual Matrix giveBMatrix(const Point *x) const { ( void ) x; return Matrix(0, 0); }; //at arbitrary point
    virtual Matrix giveBMatrix(unsigned i) const { return giveBMatrix( inttype->giveIPLocationPointer(i) ); };        //at integration point i
    Matrix giveStoredBMatrix(unsigned i) { return Bs [ i ]; };
    virtual Matrix giveHMatrix(const Point *x) const { ( void ) x; return Matrix(0, 0); };
    virtual Matrix giveHMatrix(unsigned i) const { return giveHMatrix( inttype->giveIPLocationPointer(i) ); };        //at integration point i
    Matrix giveStoredHMatrix(unsigned i) { return Hs [ i ]; };
    virtual Vector giveStrain(const Point *x, const Vector &DoFs) const { return giveBMatrix(x) * DoFs; };
    unsigned giveDimension() const { return ndim; }
    virtual Vector integrateLoad(BodyLoad *vl, double time) const;
    unsigned giveVTKCellType() const { return vtk_cell_type; };
    virtual void changeMaterial(Material *newmat);
    virtual Vector integrateInternalSources();
    //virtual void shapeF(const Point *x, MyVector &phi) const { ( void ) x; ( void ) phi; };
    //virtual double shapeFGrad(const Point *x, MyMatrix &phiGrad) const { ( void ) x; ( void ) phiGrad; return 0; };
    virtual bool giveGlobalCoords(Point *x, const Point *xn) const;
    virtual Vector giveMasterVariables(const Point *x, const Vector &DoFs) const { return giveHMatrix(x) * DoFs; };
    Vector giveElemDoFsFromFullDoFs(const Vector &FullDoFs) const;
    double giveVolume() const { return volume; };
    virtual double giveIPVolume(unsigned i) const { return inttype->giveIPWeight(i); };

    virtual void extrapolateIPValuesToNodes(std :: string code, std :: vector< Vector > &result, Vector &weights) const;
    void setMassMatrix(Matrix Q) { massM = Q; };
    virtual void collectInformationsFromNeigborhood() {};
    virtual bool isPointInside(Point *xn, const Point *x) const;
    virtual Point giveApproxCenter() const;
    Vector giveShapeFunctions(const Point *x) const;
    Matrix giveShapeFunctionsGrad(const Point *x) const;
    Point findNaturalCoords(const Point *x) const;
    virtual double giveVolumeAssociatedWithNode(unsigned nodenum)const { ( void ) nodenum; return volume / nodes.size(); }
    virtual Vector giveBoundingBox()const;
    virtual Vector findIntersectionsWithLine(Point *A, Point *B)const;
    bool doesMechanics()const { return ( physicalFields [ 0 ] > 0 ); };
    bool doesTransport()const { return ( physicalFields [ 1 ] > 0 ); };
    std :: vector< std :: vector< unsigned > >giveTraingulatedFaces()const;
    virtual double giveDissipatedEnergy() const;
    virtual void evaluateStrains(const Vector &DoFs);
    virtual void evaluateStresses(bool frozen, double timeStep);
    void removeEigenStrain();
};


//////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////
// GEOMETRICAL ELEMENT - JUST TO REPRESENT GEOMETRICAL ENTITIES
class GeometricalElement : virtual public Element
{
protected:

public:
    GeometricalElement(unsigned dim) : Element(dim) { mat = nullptr; }
    ~GeometricalElement() {};
    double giveIPValue(std :: string code, unsigned ipnum) const { ( void ) code; ( void ) ipnum; return 0; }
};

//////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////
// MATERIAL TEST ELEMENT - only one material point and virtual loading through prescribed strains
class MaterialTestElement : public Element
{
protected:

public:
    MaterialTestElement(unsigned dim);
    ~MaterialTestElement() {};
    virtual Matrix giveBMatrix(const Point *x) const;
    virtual Matrix giveHMatrix(const Point *x) const;
    virtual void setIntegrationPointsAndWeights();
};

#endif  /* _ELEMENT_H */
