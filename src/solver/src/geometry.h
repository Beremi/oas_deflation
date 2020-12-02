#ifndef _GEOMETRY_H
#define _GEOMETRY_H

#include "linear_algebra.h"
// TODO move all functions regarding geometrical operations etc.

// intersections
// lines perpendicular to plane, line or whatever...


class Region
{
private:

protected:
    string name;
    unsigned dim;

public:
    Region() {};
    virtual ~Region() {};
    virtual bool isInside(const Point &P) const = 0;
};


class RegularRegion : public Region
{
private:

protected:
    Point mainPoint;
    double size;
public:
    RegularRegion() {};
    virtual ~RegularRegion() {};
    virtual bool isInside(const Point &P) const = 0;
    virtual void setMainPoint(const Point &P) { this->mainPoint = P; };
    virtual void setSize(const double &size) { this->size = size; };
};


class Block : public RegularRegion
{
private:
    Point rightTop;

public:
    Block() {};
    virtual ~Block() {};
    Block(const Point &lB, const Point &rT);
    virtual bool isInside(const Point &P) const;
};

class Circle : public RegularRegion
{

public:
    Circle() {};
    virtual ~Circle() {};
    Circle(const Point &c, const double &r);
    virtual bool isInside(const Point &P) const;
};


class Sphere : public Circle
{
public:
    Sphere() {};
    virtual ~Sphere() {};
    Sphere(const Point &c, const double &r);
    virtual bool isInside(const Point &P) const;
};


class Polygon : public Region
{
private:
    std :: vector< Point >vertices;
    void orderClockwise(); // TODO this

public:
    Polygon() {};
    virtual ~Polygon() {};
    Polygon(const std :: vector< Point > &V);
    void addVertex(const Point &P);
    virtual bool isInside(const Point &P) const;
};








bool isInBlock(const Point &P, const Point &leftBottom, const Point &rightTop);
bool isInCircle(const Point &P, const Point &center, const double &radius);

bool onSegment(const Point &p, const Point &q, const Point &r);
int orientation(const Point &p, const Point &q, const Point &r);
bool doIntersect(const Point &p1, const Point &q1, const Point &p2, const Point &q2);
bool isInPolygon(const std :: vector< Point > &polygon, const Point &p);




#endif
