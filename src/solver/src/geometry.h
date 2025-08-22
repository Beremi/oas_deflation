#ifndef _GEOMETRY_H
#define _GEOMETRY_H

#include <memory>
#include "linalg.h"
#include "element.h"
// TODO move all functions regarding geometrical operations etc.

// intersections
// lines perpendicular to plane, line or whatever...


class Region
{
private:

protected:
    std :: string name;
    unsigned dim;

public:
    Region(unsigned d) { dim = d; };
    virtual ~Region() {};
    virtual bool isInside(const Point &P) const = 0;
    virtual void readFromLine(std :: istringstream &iss) = 0;
    virtual void setMainPoint(const Point &P) = 0;
    virtual void setSize(const double &size) = 0;
    virtual Point giveMainPoint() const = 0;
    virtual double giveSize() const = 0;
    virtual void print() const {};
    virtual void init() {};
};


class RegularRegion : public Region
{
private:

protected:
    Point mainPoint;
    double size;
public:
    RegularRegion(unsigned d) : Region(d) {};
    virtual ~RegularRegion() {};
    virtual bool isInside(const Point &P) const = 0;
    virtual void setMainPoint(const Point &P) { this->mainPoint = P; };
    virtual void setSize(const double &new_size) { this->size = new_size; };
    virtual Point giveMainPoint() const { return this->mainPoint; };
    virtual double giveSize() const { return this->size; };
    virtual void print() const {
        std :: cout << "main point: ";
        std :: cout << this->mainPoint.format(VectorSemicolonFmt) << "\n";
        std :: cout << "size = " << this->size << '\n';
    };
};


class Block : public RegularRegion
{
private:
    Point rightTop;

public:
    Block(unsigned d) : RegularRegion(d) {};
    virtual ~Block() {};
    Block(const Point &lB, const Point &rT, unsigned d);
    virtual bool isInside(const Point &P) const;
    virtual void readFromLine(std :: istringstream &iss);
};

class Circle : public RegularRegion
{
private:
    char along = 'z';
public:
    Circle() : RegularRegion(2) {};
    virtual ~Circle() {};
    Circle(const Point &c, const double &r);
    virtual bool isInside(const Point &P) const;
    virtual void readFromLine(std :: istringstream &iss);
};

//////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////
class Sphere : public Circle
{
public:
    Sphere() { dim = 3; };
    virtual ~Sphere() {};
    Sphere(const Point &c, const double &r);
    virtual bool isInside(const Point &P) const;
    virtual void readFromLine(std :: istringstream &iss);
};

//////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////
class Polygon : public Region
{
private:
    std :: vector< Point >vertices;
    void orderClockwise(); // TODO this

public:
    Polygon() : Region(2) {};
    virtual ~Polygon() {};
    Polygon(const std :: vector< Point > &V);
    void addVertex(const Point &P);
    virtual bool isInside(const Point &P) const;
};

//////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////
class Cylinder : public RegularRegion
{
private:
    Point A, B;
    double radius, length;
    Point dir;

public:
    Cylinder() : RegularRegion(3) {};
    virtual ~Cylinder() {};
    Cylinder(const Point &A, const Point &B, double radius);
    virtual bool isInside(const Point &P) const;
    virtual void readFromLine(std :: istringstream &iss);
    virtual void init();
};

//////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////
// CONTAINER FOR REGIONS
class RegionContainer
{
protected:
    std :: vector< Region * >regions;
    unsigned dim;
public:
    RegionContainer() {};
    virtual ~RegionContainer();
    void readFromFile(const std :: string &filename, unsigned dim);
    bool isLocationValid(const Point, const std :: vector< unsigned >in, const std :: vector< unsigned >out) const;
};



bool isInsideRegions(const std :: vector< std :: unique_ptr< Region > > &regions, const Point &p);
void readRegions(const std :: string &filename, std :: vector< std :: unique_ptr< Region > > &regions, unsigned d);
bool isInsideRegions(const std :: vector< std :: unique_ptr< Region > > &regions, const Element *el);
///////////////////////////////////////////////////////////////////////////////////


bool isInBlock(const Point &P, const Point &leftBottom, const Point &rightTop);
bool isInCircle(const Point &P, const Point &center, const double &radius,
                const unsigned &dir);

bool onSegment(const Point &p, const Point &q, const Point &r);
int orientation(const Point &p, const Point &q, const Point &r);
bool doIntersect(const Point &p1, const Point &q1, const Point &p2, const Point &q2);
bool isInPolygon(const std :: vector< Point > &polygon, const Point &p);



#endif /* _GEOMETRY_H */
