#ifndef _ELEMENT_FIBER_H
#define _ELEMENT_FIBER_H


#include <vector>
#include <iostream>
#include <fstream>
#include <algorithm>

#include "element.h"

class RigidBodyContact; //forward declaration

//////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////
// FIBER
class Fiber : public MechanicalElement
{
protected:
    double diam;
    Point dirVector;
    double length;
    std :: vector< RigidBodyContact * >contacts;
    std :: vector< double >positions;

public:
    Fiber(const unsigned dim);
    virtual ~Fiber() {};
    virtual void readFromLine(std :: istringstream &iss, NodeContainer *fullnodes, MaterialContainer *fullmatrs);
    virtual void init();
    Point giveDirVector()const { return dirVector; }
    double giveLength()const { return length; }
    void createNewCrossing(Point intersec, RigidBodyContact *rbc);
    void setUpCrossings();
    virtual Matrix giveBMatrix(const Point *x) const;
    virtual Matrix giveHMatrix(const Point *x) const;
    double giveDiameter() const {return diam;};
    double giveLeftLength(unsigned id) const { return positions [ id ]; };
    double giveRightLength(unsigned id) const { return length - positions [ id ]; };
};

#endif /* _ELEMENT_FIBER_H */
