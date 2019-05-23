#ifndef _ELEMENT_C_H
#define	_ELEMENT_C_H

#include "element.h"
#include "node_container.h"

//////////////////////////////////////////////////////////
class ElementContainer {
private:
    vector<Element*> elems; 
    NodeContainer *nodes; 
public:
    ElementContainer(){};    
    ~ElementContainer();
    void setNodeContainer(NodeContainer *n){nodes=n;};
    void readFromFile(const string filename, const unsigned ndim, MaterialContainer *matrs);
    void init();
    void prepareSteadyStateMatrices(CoordinateIndexedSparseMatrix &K11, CoordinateIndexedSparseMatrix &K12) const;
    void updateSteadyStateMatrices(CoordinateIndexedSparseMatrix &K11, CoordinateIndexedSparseMatrix &K12) const;
    void addBodyForces(Vector &R, double time) const;
};



#endif	/* _ELEMENT_STRUCT_C_H */
