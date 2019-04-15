/**
 * @Author: jose
 * @Date:   2019-04-05T19:16:04+02:00
 * @Last modified by:   jose
 * @Last modified time: 2019-04-10T15:46:35+02:00
 */



/*
 * File:   geometry.h
 * Author: jas
 *
 * Created on 11. september 2008, 14:58
 */

#ifndef _POINT_H
#define	_POINT_H

#include <cmath>
#include <iostream>
#include <valarray>
#include "assert.h"

#define POINT_TOLERANCE 1e-6

typedef std::valarray<double> Vector ;

class Point {
private:
    double x;
    double y;
    double z;

public:
    Point();
    Point(double x);
    Point(double x, double y);
    Point(double x, double y, double z);
    Point(const Point & p);
    ~Point(void);

    void setX(double x);
    void setY(double y);
    void setZ(double z);
    double getX() const;
    double getY() const;
    double getZ() const;
    void set(double x, double y, double z);
    void set(const Point & p);
    void set(const Point * p);

    bool operator==(const Point & p) const;
    bool operator!=(const Point & p) const;

    Point operator-(const Point & p) const;
    Point operator+(const Point & p) const;
    Point operator/(const double p) const;
    Point operator*(const double p) const;
    double operator*(const Point &p) const;
    Point normalized() const;
    void normalize();


    void operator+=(const Point &p);
    void operator-=(const Point &p);

    void operator*=(const double d) {
        x *= d;
        y *= d;
        z *= d;
    }

    void operator/=(const double d) {
        x /= d;
        y /= d;
        z /= d;
    }
    double norm() const;
    double sqNorm() const;
    void print() const;
};

Point cross(const Point &p, const Point &q);
double dot(const Point &p, const Point &q);
double squareDist(const Point &v1, const Point & v2);
double squareDist(const Point *v1, const Point *v2);
double determinant(const Point &u, const Point &v, const Point &w);

Vector pointToVector(const Point &p);  // to be able to multiply with matrix
Point vectorToPoint(const Vector &v);  // and back, because of need to count points together

#endif	/* _GEOMETRY_H */
