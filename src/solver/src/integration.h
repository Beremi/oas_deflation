#ifndef _INTEGRATION_H
#define _INTEGRATION_H

#include "node_container.h"

//////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////
// INTEGRATION POINTS - MASTER CLASS
class IntegrationType
{
private:

protected:
    string name;
    vector< Point >ip_locs;
    vector< double >ip_weights;

public:
    IntegrationType() { name = "basic integration type"; }
    virtual ~IntegrationType() {};
    virtual void init() {};
    unsigned giveNumIP() const { return ip_locs.size(); };
    double giveIPWeight(unsigned i) const { return ip_weights [ i ]; };
    void setIPWeight(unsigned i, double w) { ip_weights [ i ] = w; };
    void setIPLocation(unsigned i, Point p) { ip_locs [ i ].set(p); };
    Point giveIPLocation(unsigned i) const { return ip_locs [ i ]; };
    Point *giveIPLocationPointer(unsigned i) { return & ( ip_locs [ i ] ); };
};

//////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////
// ONE POINT INTEGRATION FOR DISCRETE MODELS
class IntegrDiscrete1 : public IntegrationType
{
public:
    IntegrDiscrete1() { name = "IntegrDiscrete1"; };
    virtual ~IntegrDiscrete1() {};
    virtual void init();
};

//////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////
// FOUR POINT INTEGRATION IN SQUARE
class IntegrQuad4 : public IntegrationType
{
public:
    IntegrQuad4() { name = "IntegrQuad4"; };
    virtual ~IntegrQuad4() {};
    virtual void init();
};


//////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////
// EIGHT POINT INTEGRATION IN CUBE
class IntegrBrick8 : public IntegrationType
{
public:
    IntegrBrick8() { name = "IntegrBrick8"; };
    virtual ~IntegrBrick8() {};
    virtual void init();
};


#endif  /* _INTEGRATION_H */
