#ifndef _ELEMENT_C_H
#define _ELEMENT_C_H

#include "element.h"
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
    unsigned giveSize() const { return elems.size(); }
    void updateMaterialStatuses();
    void prepareSteadyStateMatrices(CoordinateIndexedSparseMatrix &K) const;
    void updateSteadyStateMatrices(CoordinateIndexedSparseMatrix &K, string matrixType) const;
    void addBodyForces(Vector &R, double time) const;
    void giveInternalForces(const Vector &full_r, Vector &full_f);
    Element *giveElement(unsigned const num) { return elems [ num ]; }

    vector< Element * > :: iterator begin(){return elems.begin();}
    vector< Element * > :: iterator end(){return elems.end();}
    vector< Element * > :: const_iterator begin() const {return elems.begin();}
    vector< Element * > :: const_iterator end() const {return elems.end();}
};



#endif  /* _ELEMENT_STRUCT_C_H */
