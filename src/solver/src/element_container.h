#ifndef _ELEMENT_CONTAINER_H
#define _ELEMENT_CONTAINER_H

#include "globals.h"
#include "element.h"
#include "boundary_condition.h"
#include "node_container.h"

//////////////////////////////////////////////////////////
class ElementContainer
{
private:
    std :: vector< Element * >elems;
    NodeContainer *nodes;
    BCContainer *bconds;
    MaterialContainer *materials;
    unsigned max_sol_order = 0; //maximum number of successive rounds of internal force evaluations
    void prepareStructuralMatrix(CoordinateIndexedSparseMatrix &K, unsigned diffType) const;
    void updateStructuralMatrix(CoordinateIndexedSparseMatrix &K, unsigned diffType, std :: string matrixType) const;
    void integrateDampingOrInertiaForces(const Vector &full_v, Vector &full_f, unsigned diffType) const;
    std :: vector< std :: string >file_to_load_from;

public:
    ElementContainer() { nodes = nullptr; bconds = nullptr; materials = nullptr; };
    ~ElementContainer();
    void setContainers(NodeContainer *n, BCContainer *b) { nodes = n; bconds = b; };
    void readFromFile(const std :: string filename, const unsigned ndim, MaterialContainer *matrs);
    // void saveToFile(const std :: string &filepath, std :: vector< unsigned > &elems_to_save) const;
    void saveElemStatsToFile(const std :: string &filepath, const std :: vector< unsigned > &elems_to_save, const double time_now = 0, const unsigned step = 0, const bool saveNodeIds = true, const double idc_time = 0, const double time_step = 1e-4) const;
    void readMatStatsFromFile(double &ini_time, unsigned &ini_step, double &ini_time_step, double &ini_idc_time, const bool get_time_from_file = true);
    void setFileToLoadStatsFrom(const std :: string &str);
    void init();
    void clear();
    size_t giveSize() const { return elems.size(); }
    void findElementFriends();
    void updateMaterialStatuses();
    void resetMaterialStatuses();
    void prepareStiffnessMatrix(CoordinateIndexedSparseMatrix &K) const;
    void updateStiffnessMatrix(CoordinateIndexedSparseMatrix &K, std :: string param) const;
    void prepareDampingMatrix(CoordinateIndexedSparseMatrix &C) const;
    void updateDampingMatrix(CoordinateIndexedSparseMatrix &C) const;
    void prepareMassMatrix(CoordinateIndexedSparseMatrix &M) const;
    void updateMassMatrix(CoordinateIndexedSparseMatrix &M) const;
    void integrateInternalForces(Vector &full_r, Vector &full_f, double timeStep);                        ///< return internal forces with temporary update of internal variables
    void integrateInternalForcesWithFrozenIntVariables(Vector &full_r, Vector &full_f, double timeStep);  ///< return internal forces based on current state of internal variables
    void integrateInternalForces(const Vector &full_r, Vector &full_f, bool frozen, double timeStep);          ///< return internal forces with or without frozen internal variables
    void integrateDampingForces(const Vector &full_v, Vector &full_f) const;
    void integrateInertiaForces(const Vector &full_a, Vector &full_f) const;
    Element *giveElement(unsigned const num) const { return elems [ num ]; }
    Element *giveElementConnectingNodes(std :: vector< unsigned > &node_ids) const;
    unsigned giveElemId(const Element *elem) const;
    bool findElementOwningPoint(Element **elem, Point *xn, const Point *x) const;
    void extrapolateValuesFromIntegrationPointsToNodes(std :: string code, std :: vector< Vector > &results);
    void assignFibersToElems();
    void giveValues(std :: string code, Vector &result) const;
    void sumFromElements(std :: string code, Vector &result) const;

    std :: vector< Element * > :: iterator begin() { return elems.begin(); }
    std :: vector< Element * > :: iterator end() { return elems.end(); }
    std :: vector< Element * > :: const_iterator begin() const { return elems.begin(); }
    std :: vector< Element * > :: const_iterator end() const { return elems.end(); }
};



#endif  /* _ELEMENT_CONTAINER_H */
