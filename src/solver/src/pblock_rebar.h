#ifndef _REBAR_H
#define _REBAR_H

#include "linalg.h"
#include <vector>
#include <iostream>
#include <fstream>
#include <algorithm>
#include <typeinfo>

#include "pblock_constraints.h"

using namespace std;

//////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////
// REBAR
//////////////////////////////////////////////////////////
class Rebar : public PBlock
{
private:
protected:
    vector<unsigned> nodes;
    double radius;
    unsigned material_id;
public:
    Rebar() {};
    virtual ~Rebar() {};
    virtual void apply(Model *model);
    virtual void readFromLine(std :: istringstream &iss, unsigned d);
    virtual void findIntersectionsWithElements(Model *model);
};

#endif /* _REBAR_H */
