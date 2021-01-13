#ifndef _ELEMENT_C_H
#define _ELEMENT_C_H

#include "element.h"
#include "element_polyhedral.h"
#include "boundary_condition.h"
#include "node_container.h"

//////////////////////////////////////////////////////////
class ElementContainer
{
private:
    vector< Element * >elems;
    NodeContainer *nodes;
    BCContainer *bconds;
    unsigned max_sol_order; //maximum number of successive rounds of internal force evaluations

public:
    ElementContainer() {};
    ~ElementContainer();
    void setContainers(NodeContainer *n, BCContainer *b) { nodes = n; bconds = b; };
    void readFromFile(const string filename, const unsigned ndim, MaterialContainer *matrs);
    void init();
    size_t giveSize() const { return elems.size(); }
    void findElementFriends();
    void updateMaterialStatuses();
    void prepareSteadyStateMatrix(CoordinateIndexedSparseMatrix &K, string matrixType) const;
    void prepareSteadyStateMatrix(CoordinateIndexedSparseMatrix &K) const;
    void updateSteadyStateMatrix(CoordinateIndexedSparseMatrix &K, string matrixType) const;
    void prepareCapacityMatrix(CoordinateIndexedSparseMatrix &C) const;
    void updateCapacityMatrix(CoordinateIndexedSparseMatrix &C) const;
    void prepareMassMatrix(CoordinateIndexedSparseMatrix &M) const;
    void updateMassMatrix(CoordinateIndexedSparseMatrix &M) const;
    void giveInternalForces(Vector &full_r, Vector &full_f);                       ///< return internal forces with temporary update of internal variables
    void giveInternalForcesWithFrozenIntVariables(Vector &full_r, Vector &full_f); ///< return internal forces based on current state of internal variables
    void giveInternalForces(Vector &full_r, Vector &full_f, bool frozen);          ///< return internal forces with or without frozen internal variables
    Element *giveElement(unsigned const num) { return elems [ num ]; }

    vector< Element * > :: iterator begin() { return elems.begin(); }
    vector< Element * > :: iterator end() { return elems.end(); }
    vector< Element * > :: const_iterator begin() const { return elems.begin(); }
    vector< Element * > :: const_iterator end() const { return elems.end(); }
};



#endif  /* _ELEMENT_STRUCT_C_H */
