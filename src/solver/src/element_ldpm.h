#ifndef _ELEMENT_LDPMTETRA_H
#define _ELEMENT_LDPMTETRA_H

#include "element.h"
#include "element_discrete.h"
#include "element_ldpm.h"
#include "simplex.h"
#include  <vector>

class ElementContainer; //forward declaration

//////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////
// LDPM ELEMENT
class LDPMTetra : public Element
{
protected:
    std :: vector< Node * >vert; //cent, 0-1, 0-2, 0-3, 1-2, 1-3, 2-3, 0-1-2, 0-1-3, 0-2-3, 1-2-3
    Vector lengths, areas;//, volumes;
    std :: vector< Point >normals, t1s, t2s;
    std :: vector< Matrix >R;

    //Matrix giveRMatrix() const { return R; };
    virtual void checkNodeType() const;
    virtual void setIntegrationPointsAndWeights();

    //std :: vector< Simplex * >simplices;

    std :: vector< unsigned >nodecodes;   //connectivity of facets between nodes
    std :: vector< unsigned >vertcodes;   //connectivity of facets between vertices

    Vector volWeights; //factors to caluclate volumetric strain

    double volumetricStrain;
    virtual void computeMassMatrix();
    virtual void computeDampingMatrix();

public:
    LDPMTetra(unsigned ndim);
    ~LDPMTetra() {};
    virtual void readFromLine(std :: istringstream &iss, NodeContainer *fullnodes, MaterialContainer *fullmatrs);
    virtual void init();
    virtual Matrix giveBMatrix(unsigned k) const;
    virtual Matrix giveHMatrix(const Point *x) const;
    double giveLength(unsigned i) const { return lengths [ i ]; }
    double giveArea(unsigned i) const { return areas [ i ]; }
    Point giveNormal(unsigned i) const { return normals [ i ]; }
    Point giveT1(unsigned i) const { return t1s [ i ]; }
    Point giveT2(unsigned i) const { return t2s [ i ]; }
    std :: vector< Node * >giveVertices() const { return vert; };
    Node *giveVertex(unsigned i) const { return vert [ i ]; };
    Node *giveCentroid() const { return vert [ 0 ]; };
    unsigned giveNumOfVertices() const { return vert.size(); };

    virtual Vector giveStrain(unsigned i, const Vector &DoFs);
    virtual Matrix giveStiffnessMatrix(std :: string matrixType) const;
    virtual Vector giveInternalForces(const Vector &DoFs, bool frozen, double timeStep);
    virtual Vector integrateLoad(BodyLoad *vl, double time) const;
    virtual Vector integrateInternalSources();

    unsigned giveNumOfFacets()const { return 12; };
    std :: vector< unsigned >giveFacetVertCodes(unsigned k) const;
    std :: vector< unsigned >giveFacetNodeCodes(unsigned k) const;
    unsigned giveOppositeSurfaceVertexToNode(unsigned k) const;
    std :: vector< unsigned >giveOppositeFacetsToNode(unsigned k) const;

    double giveVolumetricStrain() const { return volumetricStrain; };
    bool isPointInside(Point *xn, const Point *x) const;

    virtual Vector  giveMasterVariables(const Point *x, const Vector &DoFs) const;
    virtual void giveValues(std :: string code, Vector &result) const;
    //Vector giveVectorToNode(const unsigned &node_i, const unsigned &ip_id) const;
    //double giveVolumeAssociatedWithNode(unsigned nodenum) const;

    //virtual void extrapolateIPValuesToNodes(std :: string code, std :: vector< Vector > &result, Vector &weights) const;
    //virtual bool isPointInside(Point *xn, const Point *x) const { ( void ) xn; ( void ) x; return false; }; //TODO: discrete elements does not interpolate
};

//////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////
// LDPM COUPLED ELEMENT
class LDPMCoupledTetra : public LDPMTetra
{
protected:

public:
    LDPMCoupledTetra();
    ~LDPMCoupledTetra() {};
    virtual Vector giveStrain(unsigned i, const Vector &DoFs);
    virtual void giveValues(std :: string code, Vector &result) const;
};

//////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////
// Discrete LDPM TRANSPORT ELEMENT
class LDPMCoupledTransport : public DiscreteTrsprtCoupledElem
{
protected:
    LDPMTetra *tetA;
    LDPMTetra *tetB;
    unsigned LDPMTetraIDA, LDPMTetraIDB;
    ElementContainer *elems;
    unsigned LDPMsideA, LDPMsideB;

public:
    LDPMCoupledTransport(ElementContainer *allelems);
    ~LDPMCoupledTransport() {};
    virtual void readFromLine(std :: istringstream &iss, NodeContainer *fullnodes, MaterialContainer *fullmatrs);
    virtual void init();
    virtual Vector giveStrain(unsigned i, const Vector &DoFs);
};

//////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////
// Discrete LDPM TRANSPORT BOUNDARY ELEMENT
class LDPMCoupledTransportBoundary : public LDPMCoupledTransport
{
protected:

public:
    LDPMCoupledTransportBoundary(ElementContainer *allelems);
    ~LDPMCoupledTransportBoundary() {};
    virtual void init();
};
#endif  /* _ELEMENT_LDPMTETRA_H */
