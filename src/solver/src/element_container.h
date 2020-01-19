#ifndef _ELEMENT_C_H
#define _ELEMENT_C_H

#include "element.h"
#include "element_polyhedral.h"
#include "node_container.h"

//////////////////////////////////////////////////////////
class ElementContainer
{
private:
    vector< Element * >elems;
    NodeContainer *nodes;
public:
    ElementContainer() {};
    ~ElementContainer();
    void setNodeContainer(NodeContainer *n) { nodes = n; };
    void readFromFile(const string filename, const unsigned ndim, MaterialContainer *matrs);
    void init();
	size_t giveSize() const { return elems.size(); }
    void updateMaterialStatuses();
    void prepareSteadyStateMatrix(CoordinateIndexedSparseMatrix &K, string matrixType) const;
    void prepareSteadyStateMatrix(CoordinateIndexedSparseMatrix &K) const;
    void updateSteadyStateMatrix(CoordinateIndexedSparseMatrix &K, string matrixType) const;
    void prepareCapacityMatrix(CoordinateIndexedSparseMatrix &C) const;
    void updateCapacityMatrix(CoordinateIndexedSparseMatrix &C) const;
    void prepareMassMatrix(CoordinateIndexedSparseMatrix &M) const;
    void updateMassMatrix(CoordinateIndexedSparseMatrix &M) const;
    void addBodyForces(Vector &R, double time) const;
    void giveInternalForces(Vector &full_r, Vector &full_f);
    Element *giveElement(unsigned const num) { return elems [ num ]; }

    vector< Element * > :: iterator begin(){return elems.begin();}
    vector< Element * > :: iterator end(){return elems.end();}
    vector< Element * > :: const_iterator begin() const {return elems.begin();}
    vector< Element * > :: const_iterator end() const {return elems.end();}
};



#endif  /* _ELEMENT_STRUCT_C_H */
