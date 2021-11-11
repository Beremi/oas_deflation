#ifndef _GEOMETRY_H
#define _GEOMETRY_H

#include <memory>
#include "linear_algebra.h"
#include "element.h"
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
    virtual void readFromLine(istringstream &iss) = 0;
    virtual void setMainPoint(const Point &P) = 0;
    virtual void setSize(const double &size) = 0;
    virtual Point giveMainPoint() const = 0;
    virtual double giveSize() const = 0;
    virtual void print() const {};
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
    virtual void setSize(const double &new_size) { this->size = new_size; };
    virtual Point giveMainPoint() const { return this->mainPoint; };
    virtual double giveSize() const { return this->size; };
    virtual void print() const {
        std :: cout << "main point: ";
        this->mainPoint.print();
        std :: cout << "size = " << this->size << '\n';
    };
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
    virtual void readFromLine(istringstream &iss);
};

class Circle : public RegularRegion
{
public:
    Circle() {};
    virtual ~Circle() {};
    Circle(const Point &c, const double &r);
    virtual bool isInside(const Point &P) const;
    virtual void readFromLine(istringstream &iss);
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



bool isInsideRegions(const std :: vector< std :: unique_ptr< Region > > &regions, const Point &p);
void readRegions(const std :: string &filename, std :: vector< std :: unique_ptr< Region > > &regions);
bool isInsideRegions(const std :: vector< std :: unique_ptr< Region > > &regions, const Element *el);
///////////////////////////////////////////////////////////////////////////////////


bool isInBlock(const Point &P, const Point &leftBottom, const Point &rightTop);
bool isInCircle(const Point &P, const Point &center, const double &radius,
                const unsigned dir = 2);

bool onSegment(const Point &p, const Point &q, const Point &r);
int orientation(const Point &p, const Point &q, const Point &r);
bool doIntersect(const Point &p1, const Point &q1, const Point &p2, const Point &q2);
bool isInPolygon(const std :: vector< Point > &polygon, const Point &p);



#endif
