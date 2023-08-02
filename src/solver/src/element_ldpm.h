#ifndef _ELEMENT_LDPMTETRA_H
#define _ELEMENT_LDPMTETRA_H

#include "element.h"
#include "simplex.h"
#include  <vector>

//////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////
// RBSN ELEMENT
class LDPMTetra : public Element
{
protected:
    std :: vector< Node * >vert;
    Vector lengths, areas, volumes;
    std :: vector< Point >normals;
    std :: vector< Matrix >R;

    //Matrix giveRMatrix() const { return R; };
    virtual void checkNodeType() const;
    virtual void setIntegrationPointsAndWeights();

    //std :: vector< Simplex * >simplices;

    std :: vector< unsigned >nodecodes;   //coonectivity of facets between nodes
    std :: vector< unsigned >vertcodes;   //coonectivity of facets between vertices

    Vector volWeights; //factors to caluclate volumetric strain

    double volumetricStrain;

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
    std :: vector< Node * >giveVertices() const { return vert; };
    Node *giveVertex(unsigned i) const { return vert [ i ]; };
    unsigned giveNumOfVertices() const { return vert.size(); };

    virtual Vector giveStrain(unsigned i, const Vector &DoFs);
    virtual Matrix giveStiffnessMatrix(std :: string matrixType) const;
    virtual Matrix giveDampingMatrix() const;
    virtual Matrix giveMassMatrix() const;
    virtual Vector giveInternalForces(const Vector &DoFs, bool frozen, double timeStep);
    virtual Vector integrateLoad(BodyLoad *vl, double time) const;
    virtual Vector integrateInternalSources();

    unsigned giveNumOfFacets()const { return 12; };
    std :: vector< unsigned >giveFacetVertCodes(unsigned k) const;
    std :: vector< unsigned >giveFacetNodeCodes(unsigned k) const;

    //virtual void giveValues(std :: string code, Vector &result) const;
    //Vector giveVectorToNode(const unsigned &node_i, const unsigned &ip_id) const;
    //double giveVolumeAssociatedWithNode(unsigned nodenum) const;

    //virtual void extrapolateIPValuesToNodes(std :: string code, std :: vector< Vector > &result, Vector &weights) const;
    //virtual bool isPointInside(Point *xn, const Point *x) const { ( void ) xn; ( void ) x; return false; }; //TODO: discrete elements does not interpolate
};

#endif  /* _ELEMENT_LDPMTETRA_H */
