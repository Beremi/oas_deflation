/**
 * @Author: jose
 * @Date:   2019-04-05T18:20:22+02:00
 * @Last modified by:   jose
 * @Last modified time: 2019-04-05T19:01:51+02:00
 */



#include "point.h"

Point::Point() {
    x = 0;
    y = 0;
    z = 0;
}

Point::Point(double ox) {
    x = ox;
    y = 0;
    z = 0;
}

Point::Point(double ox, double oy) {
    x = ox;
    y = oy;
}

Point::Point(double ox, double oy, double oz) {
    x = ox;
    y = oy;
    z = oz;
}

Point::Point(const Point & p) {
    x = p.x;
    y = p.y;
    z = p.z;
}

Point::~Point(void) {

}

void Point::setX(double ox) {
    x = ox;
}

void Point::setY(double oy) {
    y = oy;
}

void Point::setZ(double oz) {
    z = oz;
}

double Point::getX() const {
    return x;
}

double Point::getY() const {
    return y;
}

double Point::getZ() const {
    return z;
}

void Point::set(double ox, double oy, double oz) {
    x = ox;
    y = oy;
    z = oz;
}

void Point::set(const Point & p) {
    x = p.x;
    y = p.y;
    z = p.z;
}

void Point::set(const Point * p) {
    x = p->x;
    y = p->y;
    z = p->z;
}

bool Point::operator==(const Point &p) const {
    double delta = POINT_TOLERANCE;
    Point a(p);
    a.x += delta;
    a.y += delta;
    a.z += delta;
    Point b(p);
    b.x += delta;
    b.y += delta;
    b.z -= delta;
    Point c(p);
    c.x += delta;
    c.y -= delta;
    c.z += delta;
    Point d(p);
    d.x += delta;
    d.y -= delta;
    d.z -= delta;
    Point e(p);
    e.x -= delta;
    e.y += delta;
    e.z += delta;
    Point f(p);
    f.x -= delta;
    f.y += delta;
    f.z -= delta;
    Point g(p);
    g.x -= delta;
    g.y -= delta;
    g.z += delta;
    Point h(p);
    h.x -= delta;
    h.y -= delta;
    h.z -= delta;

    return squareDist(&p, this) < POINT_TOLERANCE * POINT_TOLERANCE
            || squareDist(&a, this) < POINT_TOLERANCE * POINT_TOLERANCE
            || squareDist(&b, this) < POINT_TOLERANCE * POINT_TOLERANCE
            || squareDist(&c, this) < POINT_TOLERANCE * POINT_TOLERANCE
            || squareDist(&d, this) < POINT_TOLERANCE * POINT_TOLERANCE
            || squareDist(&e, this) < POINT_TOLERANCE * POINT_TOLERANCE
            || squareDist(&f, this) < POINT_TOLERANCE * POINT_TOLERANCE
            || squareDist(&g, this) < POINT_TOLERANCE * POINT_TOLERANCE
            || squareDist(&h, this) < POINT_TOLERANCE*POINT_TOLERANCE;
    return std::abs(x - p.x) < .5 * POINT_TOLERANCE &&
            std::abs(y - p.y) < .5 * POINT_TOLERANCE &&
            std::abs(z - p.z) < .5 * POINT_TOLERANCE
            ;
    return squareDist(&p, this) < POINT_TOLERANCE*POINT_TOLERANCE;
}

bool Point::operator!=(const Point & p) const {
    return squareDist(this, &p) > POINT_TOLERANCE*POINT_TOLERANCE;
}

Point Point::operator-(const Point &p) const {
    Point ret((*this));
    ret.x -= p.x;
    ret.y -= p.y;
    ret.z -= p.z;
    return ret;
}

Point Point::operator+(const Point &p) const {
    Point ret((*this));
    ret.x += p.x;
    ret.y += p.y;
    ret.z += p.z;
    return ret;
}

void Point::operator+=(const Point &p) {
    x += p.x;
    y += p.y;
    z += p.z;
}

void Point::operator-=(const Point &p) {
    x -= p.x;
    y -= p.y;
    z -= p.z;
}

Point Point::operator/(const double p) const {
    Point ret((*this));
    double inv = 1. / p;
    ret.x *= inv;
    ret.y *= inv;
    ret.z *= inv;
    return ret;

}

Point Point::operator*(const double p) const {
    Point ret((*this));
    ret.x *= p;
    ret.y *= p;
    ret.z *= p;
    return ret;

}

double Point::operator*(const Point &p) const {
    return x * p.x + y * p.y + z * p.z;

}

Point Point::normalized() const {
    Point ret((*this));
    ret.normalize();
    return ret;
}

void Point::normalize() {
    double norm = this->norm();
    x = x / norm;
    y = y / norm;
    z = z / norm;
}

double Point::norm() const {
    return sqrt(x * x + y * y + z * z);
}

double Point::sqNorm() const {
    return x * x + y * y + z*z;
}

void Point::print() const {
    std::cout << getX() << std::flush;
    std::cout << "; " << getY() << std::flush;
    std::cout << "; " << getZ() << std::flush << std::endl;
}

Point cross(const Point &p, const Point &q) {
    Point r((p.y * q.z - p.z * q.y), (p.z * q.x - p.x * q.z), (p.x * q.y - p.y * q.x));
    return r;
}

double dot(const Point &p, const Point &q) {
    double r(p.x * q.x + p.y * q.y + p.z * q.z);
    return r;
}

double squareDist(const Point &v1, const Point & v2) {
    double x = v2.getX() - v1.getX();
    double y = v2.getY() - v1.getY();
    double z = v2.getZ() - v1.getZ();

    return x * x + y * y + z*z;
}

double squareDist(const Point *v1, const Point *v2) {
    double x = v2->getX() - v1->getX();
    double y = v2->getY() - v1->getY();
    double z = v2->getZ() - v1->getZ();
    return x * x + y * y + z*z;
}

double determinant(const Point &u, const Point &v, const Point &w) {
    double D = 0.0;
    D = u.x * v.y * w.z - u.x * v.z * w.y + u.y * v.z * w.x - u.y * v.x * w.z + u.z * v.x * w.y - u.z * v.y * w.x;
    return D;
}

Vector pointToVector(const Point &p){
    Vector vect(3);
    vect[0] = p.getX();
    vect[1] = p.getY();
    vect[2] = p.getZ();
    return vect;
}

Point vectorToPoint(const Vector &v){
    Point p;
    p.setX(v[0]);
    p.setY(v[1]);
    p.setZ(v[2]);
    return p;
}
