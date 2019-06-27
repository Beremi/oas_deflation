/*
 * vector.h
 *
 *  Created on: 17-10-2012
 *      Author: Robson
 */

#ifndef ALGEBRA_H_
#define ALGEBRA_H_

#include <model/model.hpp>
#include <exceptions/exceptions.hpp>
#include <3rd/eigen/Eigen/Dense>

namespace dmga{

namespace algebra{

using namespace Eigen;
/** Vector3D is a class that allows for classic vector operation, we use Eigen implementation here */
typedef typename Eigen::Vector3d Vector3D;
/** Matrix3D is a class that allows for classic matrix-matrix-or-vector operation, we use Eigen implementation here */
typedef typename Eigen::Matrix3d Matrix3D;

/**
 * defines a 2D plane in 3D geometry by specifying
 * point in space and normal vector.
 */
class Plane3D{
public:
	Vector3D p;
	Vector3D n;
	double d;

	/** creates new plane in point p with normal n */
	Plane3D(Vector3D p, Vector3D n) : p(p), n(n){
		d = n.dot(p);
	}
	/** creates default plane in (0,0,0) and with normal (1,0,0) */
	Plane3D() : p(0,0,0), n(1,0,0){
		d = n.dot(p);
	}
	/** returns point of the plane */
	Vector3D point(){
		return p;
	}
	/** return normal to the plane */
	Vector3D normal(){
		return n;
	}
	/** output operator */
	friend std::ostream& operator<<(std::ostream& out, const Plane3D& item){
		out << "p = " << item.p << ", n = " << item.n << ", d = " << item.d;
		return out;
	}
};

/**
 * computes parametrization of the power plane given by the
 * power distance of two balls of the radius r - it is necessary for
 * the radical tesselation (voronoi diagram)
 *
 * we assume that coords1 and coords2 holds the information
 * on the position (x, y, z) and r - the radius
 * This is true for the Voro++ implementation
 */
Plane3D ballsToPowerPlane(double* coords1, double* coords2){
	//TODO: optimize it for Voro++, it can be better, and
	//can use r_scale from rad_option.hh etc.
	//r_scale liczy to co chce tutaj policzyc jako d
	//tylko bierze do tego info zawarte w containerze
	DEB(coords1[0] << " " << coords1[1] << " " << coords1[2] << " " << coords2[0] << " " << coords2[1] << " " << coords2[2] << " " << std::flush);
	double xx = coords2[0] - coords1[0];
	double yy = coords2[1] - coords1[1];
	double zz = coords2[2] - coords1[2];
	double r1pow2 = coords1[3] * coords1[3];
	double r2pow2 = coords2[3] * coords2[3];
	double lenpow2 = xx*xx + yy*yy + zz*zz;
	double len = sqrt(lenpow2);
	double d = (lenpow2 + (r1pow2 - r2pow2)) / (2.0 * len);
	Vector3D n(xx/len, yy/len, zz/len);
	Vector3D p = d * n + Vector3D(coords1[0], coords1[1], coords1[2]);
	return Plane3D(p, n);
}

/**
 * computes point at which three 2D Planes in 3D intersects
 * it checks if intersection can be non-unique or non existing and then throws
 * exception.
 */
Vector3D planes3Intersection(Plane3D& pl1, Plane3D& pl2, Plane3D& pl3){
	Vector3D res;
	//P_0 = (d_1 (n_2 x n_3) + d_2 (n_3 x n_1) + d_3 (n_1 x n_2)) / (n_1 * (n_2 x n_3))
	// d_i = p_i * n_i

	Vector3D n_2xn_3 = pl2.n.cross(pl3.n);
	double dot = (pl1.n.dot(n_2xn_3));
	if (dot < 0.00001){ //TODO: choose good constant here!
		throw dmga::exceptions::NotImplementedYet("planes3Intersection(): possible not unique solution, dot < 0.00001");
	}
	return (pl1.d * n_2xn_3 + pl2.d * pl3.n.cross(pl1.n) + pl3.d * pl1.n.cross(pl2.n)) / dot;
}

/**
 * computes normal to the plane defined two vectors an stores it in result
 * it also stores the area of the triangle formed with those 3 vertices in result[3] so result must be at least 4 in size
 */
inline void normal3d_raw(double* v1, double* v2, double* result){
	result[0] = v1[1] * v2[2] - v1[2] * v2[1];
	result[1] = v1[2] * v2[0] - v1[0] * v2[2];
	result[2] = v1[0] * v2[1] - v1[1] * v2[0];
	result[3] = sqrt(result[0] * result[0] + result[1] * result[1] + result[2] *result[2]);
	if (result[3] > 0){
		result[0] /= result[3];
		result[1] /= result[3];
		result[2] /= result[3];
	}
}
/**
 * computes normal to the plane defined with triangle v1, v2, v3, an stores it in result
 * it also stores the area of the triangle formed with those 3 vertices in result[3] so result must be at least 4 in size
 */
inline void normal3d_raw(double* v1, double* v2, double* v3, double* result){
	double a1 = v2[0] - v1[0];
	double a2 = v2[1] - v1[1];
	double a3 = v2[2] - v1[2];
	double b1 = v3[0] - v1[0];
	double b2 = v3[1] - v1[1];
	double b3 = v3[2] - v1[2];
	result[0] = a2 * b3 - a3 * b2;
	result[1] = a3 * b1 - a1 * b3;
	result[2] = a1 * b2 - a2 * b1;
	result[3] = sqrt(result[0] * result[0] + result[1] * result[1] + result[2] *result[2]);
	if (result[3] > 0){
		result[0] /= result[3];
		result[1] /= result[3];
		result[2] /= result[3];
	}
}

/**
 * test if the angle at p in triangle defined by points v1, p, v2 is acute (ostry)
 */
inline bool is_acute_angle(double* v1, double* p, double* v2){
	return ((v1[0] - p[0]) * (v2[0] - p[0]) + (v1[1] - p[1]) * (v2[1] - p[1]) + (v1[2] - p[2]) * (v2[2] - p[2])) > 0 ;
}

/**
 * ccw in 3D - returns
 */
inline bool ccw3d(double* v1, double* p, double* v2){
	throw exceptions::NotImplementedYet("algebra::ccw3d()");

	double* normal; normal3d_raw(p, v1, v2, normal);
	return ((v1[0] - p[0]) * (v2[0] - p[0]) + (v1[1] - p[1]) * (v2[1] - p[1]) + (v1[2] - p[2]) * (v2[2] - p[2])) > 0;
}

/**
 * cw in 3D - returns
 */
inline bool cw3d(double* v1, double* p, double* v2){
	throw exceptions::NotImplementedYet("algebra::cw3d()");

	return ((v1[0] - p[0]) * (v2[0] - p[0]) + (v1[1] - p[1]) * (v2[1] - p[1]) + (v1[2] - p[2]) * (v2[2] - p[2])) > 0;
}

namespace distance{
/**
 * computes euclidean distance squared in 3D, it is faster and more accurate than computing euclid3d()^2
 */
inline double euclid3d_pow_2(double* p1, double* p2){
	return (p1[0] - p2[0]) * (p1[0] - p2[0]) + (p1[1] - p2[1]) * (p1[1] - p2[1]) + (p1[2] - p2[2]) * (p1[2] - p2[2]);
}
/**
 * computes euclidean distance in 3D
 */
inline double euclid3d(double* p1, double* p2){
	return sqrt(euclid3d_pow_2(p1, p2));
}
/**
 * computes distance d from p1 to the point on power plane defined by
 * line connecting balls p1 and p2 (p1, and p2 has 4 coordinates: x,y,z and r - radius)
 * that is the point on the dividing plane is p1 + normal * d
 */
inline double from_power_plane(double* p1, double* p2, double* normal = 0){
	double xx = p2[0] - p1[0];
	double yy = p2[1] - p1[1];
	double zz = p2[2] - p1[2];
	double r1pow2 = p1[3] * p1[3];
	double r2pow2 = p2[3] * p2[3];
	double lenpow2 = xx*xx + yy*yy + zz*zz;
	double len = sqrt(lenpow2);
	if (normal){
		normal[0] = xx / len;
		normal[1] = yy / len;
		normal[2] = zz / len;
	}
	return (lenpow2 + (r1pow2 - r2pow2)) / (2.0 * len);
}
/**
 * computes euclidean distance from point p to plane defined by triangle (v1,v2,v3)
 */
inline double euclid3d(double* p, double* v1, double* v2, double* v3){
	double n[4];
	normal3d_raw(v1, v2, v3, n);
	return Vector3D(n[0], n[1], n[2]).dot(Vector3D(p[0] - v1[0], p[1] - v1[1], p[2] - v1[2]));
}
/**
 * computes squared euclidean distance from point p to plane defined by triangle (v1,v2,v3)
 */
inline double euclid3d_pow_2(double* p, double* v1, double* v2, double* v3){
	double d = euclid3d(p, v1, v2, v3);
	return d * d;
}

/**
 * computes euclidean distance from point p to plane defined by triangle (v1,v2,v3) and returns normal (remember that it also contains the area of the triangle!)
 */
inline double euclid3d(double* p, double* v1, double* v2, double* v3, double* n){
	normal3d_raw(v1, v2, v3, n);
	return Vector3D(n[0], n[1], n[2]).dot(Vector3D(p[0] - v1[0], p[1] - v1[1], p[2] - v1[2]));
}
/**
 * computes squared euclidean distance from point p to plane defined by triangle (v1,v2,v3) and returns normal (remember that it also contains the area of the triangle, so you need n[4] for normal!)
 */
inline double euclid3d_pow_2(double* p, double* v1, double* v2, double* v3, double* n){
	double d = euclid3d(p, v1, v2, v3, n);
	return d * d;
}

} //namespace distance

/**
 * returns the area of the triangle in 3D defined by the points p, v, u
 */
inline double triangle3d_area(double* p, double* v, double* u){
	double n[4];
	normal3d_raw(p, v, u, n);
	return n[3] / 2.0;
}
/**
 * returns dot product of 2 vectors in 3d given in raw coords
 */
double dot3d(double* v, double* u){
	return v[0] * u[0] + v[1] * u[1] + v[2] * u[2];
}
/**
 * returns dot product of a triangle (v, p, u) in the point p in 3d
 */
double dot3d(double* v, double * p, double* u){
	return (v[0] - p[0]) * (u[0] - p[0]) + (v[1] - p[1]) * (u[1] - p[1]) + (v[2] - p[2]) * (u[2] - p[2]);
}
/**
 * returns dot product of 2 vectors in 3d
 */
double dot3d(double x1, double y1, double z1, double x2, double y2, double z2){
	return x1 * x2 + y1 * y2 + z1 * z2;
}
/**
 * convex combination of two points in \mathbb{R}
 */
double convex_combination(double x0, double x1, double t){
	return t * x1 + (1 - t) * x0;
}

/**
 * we assume that input vector n is orthonormal
 *
 * in (x1, y1, z1) and (x2, y2, z2) generates
 * two orthonormal to each other vectors, that are also
 * orthonormal to given vector (nx, ny, nz)
 *
 * this is useful in case we have normal vector to the plane
 * and we want to have basis for this plane
 */
void generate_ortho_basis(double nx, double ny, double nz, double& x1, double& y1, double& z1, double& x2, double& y2, double& z2){
	//we assume we have (nx, ny, nz) vector normal to plane that contains (0,0,0)
	//we want to generate two other points on this plane such that (n, v1, v2) is a orthogonal basis
	x1 = ny; y1 = nz; z1 = nx; //first guess it is usually out of the plane, do not work if (nx, ny, nz) ~= (1,1,1) <- TODO: naprawić
	x2 = nz; y2 = nx; z2 = ny; //first guess it is usually out of the plane, do not work if (nx, ny, nz) ~= (1,1,1) <- TODO: naprawić
	double d;
	//cast it back to the plane
	d = dot3d(nx, ny, nz, x1, y1, z1);
	x1 = x1 - d * nx; y1 = y1 - d * ny; z1 = z1 - d * nz;
	d = dot3d(nx, ny, nz, x2, y2, z2);
	x2 = x2 - d * nx; y2 = y2 - d * ny; z2 = z2 - d * nz;
	//normalize v2
	d = sqrt(dot3d(x2, y2, z2, z2, y2, z2));
	x2 = x2 / d; y2 = y2 / d; z2 = z2 / d;
	//orthogonalize v1 using v2
	d = dot3d(x1, y1, z1, x2, y2, z2);
	x1 = x1 - d * x2; y1 = y1 - d * y2; z1 = z1 - d * z2;
	//normalize v1
	d = sqrt(dot3d(x1, y1, z1, z1, y1, z1));
	x1 = x1 / d; y1 = y1 / d; z1 = z1 / d;
}

/**
 * creates a transformation based on a given plane
 * n is a vector normal to plane, output is a matrix
 * that sends n to (0, 0, 1) and some ortognormal basis
 * on the plane to other two vectors
 */
Matrix3D generate_plane_transform(double nx, double ny, double nz){
	Matrix3D A;
	double v[3][3]; v[0][0] = nx; v[1][0] = ny; v[2][0] = nz;
	generate_ortho_basis(v[0][0], v[1][0], v[2][0], v[0][1], v[1][1], v[2][1], v[0][2], v[1][2], v[2][2]);
	for (int i = 0; i < 3; ++i)
		for (int j = 0; j < 3; ++j)
			A(i,j) = v[i][j];
	return A.inverse();
}
/**
 * helper function
 */
void mul(Matrix3D& A, double* x, double* result){
	Vector3D v; v[0] = x[0]; v[1] = x[1]; v[2] = x[2];
	Vector3D r = A * v; result[0] = r[0]; result[1] = r[1]; result[2] = r[2];
}

bool normalize3d_raw(double* v){
	double d = sqrt(dot3d(v, v));
	if (d > 0){
		v[0] /= d; v[1] /= d; v[2] /= d;
		return true;
	}
	return false;
}

}//namespace algebra

}//namespace dmga


#endif /* VECTOR_H_ */
