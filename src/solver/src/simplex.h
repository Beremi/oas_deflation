#ifndef _SIMPLEX_H
#define _SIMPLEX_H

#include "linear_algebra.h"
#include "node.h"
#include <vector>
#include <iostream>

class Simplex
{
protected:
    vector<Particle *> nodes; ///list of nodes
public:
    Simplex(){};
    virtual ~Simplex() {};
    
};


#endif /* _SIMPLEX_H */
