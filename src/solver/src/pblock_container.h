#ifndef _PREPROCESSING_BLOCK_CONTAINER_H
#define _PREPROCESSING_BLOCK_CONTAINER_H

#include "linalg.h"
#include <vector>
#include <iostream>
#include <fstream>
#include <algorithm>
#include <typeinfo>

// #include "boundary_condition.h"
// #include "node_container.h"
// #include "element_container.h"
#include "pblock_constraints.h"
#include "pblock_periodic_bc.h"
#include "pblock_rebar.h"
#include "data_exporter.h"
#include "element_discrete.h"
#include "geometry.h"

class Model; //

//////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////
// CONTAINER FOR PREPROCESSOR BLOCKS
class PBlockContainer
{
private:
    std :: vector< PBlock * >blocks;
    Model *model = nullptr;
public:
    PBlockContainer() {};
    virtual ~PBlockContainer();
    void readFromFile(const std :: string filename, unsigned dim);
    void setModel(Model *m);
    void init();
    void clear();
    unsigned giveSize() const { return blocks.size(); };
    PBlock *givePBlock(unsigned i) { return blocks [ i ]; };
protected:
};


#endif /* _PREPROCESSING_BLOCK_CONTAINER_H */
