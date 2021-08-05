#ifndef _INTEGRATION_H
#define _INTEGRATION_H

#include "node_container.h"
#include "shape_functions.h"

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
    unsigned giveNumIP() const { return ip_locs.size(); };
    double giveIPWeight(unsigned i) const { return ip_weights [ i ]; };
    void setIPWeight(unsigned i, double w) { ip_weights [ i ] = w; };
    void setIPLocation(unsigned i, Point p) { ip_locs [ i ].set(p); };
    Point giveIPLocation(unsigned i) const { return ip_locs [ i ]; };
    Point *giveIPLocationPointer(unsigned i) { return & ( ip_locs [ i ] ); };
    virtual void init();
    virtual void init(const vector< Node* > & nodes);
    virtual void init(const vector< Node* > & nodes, const vector< vector< unsigned > > & faces, Point *centroid);
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
// THREE POINT INTEGRATION IN TRIANGLE
class IntegrTri3 : public IntegrationType
{
public:
    IntegrTri3() { name = "IntegrTri3"; };
    virtual ~IntegrTri3() {};
    virtual void init();
};


//////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////
// THREE POINT INTEGRATION IN POLYGON
class IntegrPolygon : public IntegrationType
{
    string ip_type;
public:
    IntegrPolygon(string type) { name = "IntegrPolygon"; ip_type=type;};
    virtual ~IntegrPolygon() {};
    virtual void init(const vector< Node* > & nodes, const vector< vector< unsigned > > & faces, Point *centroid);
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
