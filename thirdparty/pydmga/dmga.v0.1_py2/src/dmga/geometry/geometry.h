/*
 * geometry.h
 *
 *  Created on: 10-11-2012
 *      Author: Robson
 */

#ifndef GEOMETRY_H_
#define GEOMETRY_H_

#include <base.h>
#include <exceptions/exceptions.hpp>
#include <3rd/voro/voro++.hh>

//TYLKO DLA PROBY Z PRESETAMI!
#include <model/model.hpp>
#include <utils/voropp_container_relocator.hpp>

namespace dmga{

namespace geometry{


class GeometryPreset : public DMGAObject{
public:
	std::string show(){ return "dmga::geometry::GeometryPreset"; }

	GeometryPreset(){
	}

	virtual bool on_preset(int id){
		return false;
	}

	virtual ~GeometryPreset(){
	}
};

class SphereGeometryPreset : public voro::wall, public GeometryPreset{
public:
	std::string show(){ return "dmga::geometry::SphereGeometryPreset"; }

	double xc, yc, zc, rc;
	int w_id;

	SphereGeometryPreset(double x, double y, double z, double r) : xc(x), yc(y), zc(z), rc(r), w_id(-99){
	}

	bool point_inside(double x,double y,double z){
		return (x-xc)*(x-xc)+(y-yc)*(y-yc)+(z-zc)*(z-zc)<rc*rc;
	}
	bool cut_cell(voro::voronoicell &c,double x,double y,double z){
		double xd=x-xc,yd=y-yc,zd=z-zc,dq=xd*xd+yd*yd+zd*zd;
		if (dq>1e-5) {
			dq=2*(sqrt(dq)*rc-dq);
			return c.nplane(xd,yd,zd,dq,w_id);
		}
		return true;
	}
	bool cut_cell(voro::voronoicell_neighbor &c,double x,double y,double z){
		double xd=x-xc,yd=y-yc,zd=z-zc,dq=xd*xd+yd*yd+zd*zd;
		if (dq>1e-5) {
			dq=2*(sqrt(dq)*rc-dq);
			return c.nplane(xd,yd,zd,dq,w_id);
		}
		return true;
	}
	bool on_preset(int id){
		return w_id == id;
	}
};

class OneCutCylinderGeometryPreset : public voro::wall, public GeometryPreset{
public:
	std::string show(){ return "dmga::geometry::OneCutCylinderGeometryPreset"; }

	const int w_id;
	const double xc,yc,zc,xa,ya,za,asi,rc;
//	voro::wall_cylinder w;

	OneCutCylinderGeometryPreset(double r, double x, double y, double z, double vx, double vy, double vz) :
			w_id(-99), xc(x), yc(y), zc(z), xa(vx), ya(vy), za(vz),
			asi(1/(vx*vx+vy*vy+vz*vz)), rc(r) {
	}

	bool point_inside(double x,double y,double z){
		double xd=x-xc,yd=y-yc,zd=z-zc;
		double pa=(xd*xa+yd*ya+zd*za)*asi;
		xd-=xa*pa;yd-=ya*pa;zd-=za*pa;
		return xd*xd+yd*yd+zd*zd<rc*rc;
	}
	bool cut_cell(voro::voronoicell &c,double x,double y,double z){
		double xd=x-xc,yd=y-yc,zd=z-zc,pa=(xd*xa+yd*ya+zd*za)*asi;
		xd-=xa*pa;yd-=ya*pa;zd-=za*pa;
		pa=xd*xd+yd*yd+zd*zd;
		if(pa>1e-5) {
			pa=2*(sqrt(pa)*rc-pa);
			return c.nplane(xd,yd,zd,pa,w_id);
		}
		return true;
	}
	bool cut_cell(voro::voronoicell_neighbor &c,double x,double y,double z){
		double xd=x-xc,yd=y-yc,zd=z-zc,pa=(xd*xa+yd*ya+zd*za)*asi;
		xd-=xa*pa;yd-=ya*pa;zd-=za*pa;
		pa=xd*xd+yd*yd+zd*zd;
		if(pa>1e-5) {
			pa=2*(sqrt(pa)*rc-pa);
			return c.nplane(xd,yd,zd,pa,w_id);
		}
		return true;
	}

	bool on_preset(int id){
		return w_id == id;
	}
};

class ApproximateCylinderGeometryPreset : public voro::wall, public GeometryPreset{
public:
	std::string show(){ return "dmga::geometry::ApproximateCylinderGeometryPreset"; }

	const int w_id;
	const double r, H;
	const int n;

	ApproximateCylinderGeometryPreset(double r, double H, int precision = 32) :
			w_id(-99), r(r), H(H), n(precision) {
	}

	bool point_inside(double x,double y,double z){
		return true; //(x*x + y*y <= r * r) && (z >= 0 && z <= H);
	}
	bool cut_cell(voro::voronoicell &c,double x,double y,double z){
		//return c.nplane(-x, 0, 0, 2*x*x, w_id);
		double xp, yp, ca;
		double phi = 2.0 * 3.14159265359 / (double)n;
		for (int i = 0; i < n; ++i){
			xp = cos(phi * (double)i);
			yp = sin(phi * (double)i);
			ca = xp * x + yp * y; //cast on xp, yp vector...
			xp = xp * (r - ca); //rescale the vector
			yp = yp * (r - ca); //rescale the vector
			c.nplane(xp, yp, 0.0, (xp*xp + yp*yp), w_id);
		}
		return true;
	}
	bool cut_cell(voro::voronoicell_neighbor &c,double x,double y,double z){
		//return c.nplane(-x, 0, 0, 2*x*x, w_id);
		double xp, yp, ca;
		double phi = 2.0 * 3.14159265359 / (double)n;
		for (int i = 0; i < n; ++i){
			xp = cos(phi * (double)i);
			yp = sin(phi * (double)i);
			ca = xp * x + yp * y; //cast on xp, yp vector...
			xp = xp * (r - ca); //rescale the vector
			yp = yp * (r - ca); //rescale the vector
			c.nplane(xp, yp, 0.0, 2 * (xp*xp + yp*yp), w_id);
		}
		return true;
	}

	bool on_preset(int id){
		return w_id == id;
	}

	virtual ~ApproximateCylinderGeometryPreset(){
		//to make sure it is virtual
	}
};


class GenericGeometry : public DMGAObject{
public:
	std::string show(){ return "dmga::geometry::GenericGeometry"; }

	/**
	 * return bounding box of this geometry in raw
	 */
	virtual void boundingBox(double &x1,
							double &y1,
							double &z1,
							double &x2,
							double &y2,
							double &z2) = 0;
	/**
	 * return bounding box of this geometry in raw
	 * default impl. delegates to boundingBox(x1,..,z2).
	 */
	virtual void boundingBox(double* raw){
		boundingBox(raw[0], raw[1], raw[2], raw[3], raw[4], raw[5]);
	}
	/**
	 * destructor - because of virtual
	 */
	virtual ~GenericGeometry(){			DEB3("GenericGeometry::__destruct__():");
	}
	/**
	 * we assume that generic geometry can be
	 * enchanced with some geometry presets that
	 * can be used to define borders of geometry
	 * (mainly walls from Voro++ library)
	 */
	virtual bool addPreset(GeometryPreset& preset){
		return false;
	}
	/**
	 * transformaton of coordinates - default do nothing
	 * some geometries may want to transform coordinates
	 * it is used mainly in cast-geometries, where
	 * we want to have atoms aligned near some 2D manifold
	 * we assume that transformaton can act only on first 3
	 * coordinates (position)
	 */
	virtual void transform(int &id, double &x, double &y, double &z) {};
	/**
	 * transformaton of coordinates
	 * some geometries may want to transform coordinates
	 * it is used mainly in cast-geometries, where
	 * we want to have atoms aligned near some 2D manifold
	 * this is for arbitrary things and may be faster in
	 * some settings
	 * we assume that raw contains at least 3 first coordinates
	 * default implementation is that it delegates to transform(id, x, y, z).
	 */
	virtual void transform(int &id, double* raw) { transform(id, raw[0], raw[1], raw[2]); };
	/**
	 * returns true if given neigbour is a boundary
	 * usually it has id < 0, so we use it.
	 * this test is used when CastGeometry is used
	 * to find approximate voronoi diagram of
	 * projection to some 2D manifold.
	 */
	virtual bool on_boundary(int id){
		return id < 0;
	}
};
/**
 * Default base class geometry for Voro++ voronoi implementation
 */
class VoroppGeometry : public GenericGeometry{
public:
	std::string show(){ return "dmga::geometry::VoroppGeometry"; }

	/** this could be done in other ways I think... */
	typedef voro::container_relocator WEIGHTED_STORAGE; //should be poly - relocator is for testing purposes, and if using relocating strategy (which is not better than deleting & recreating)
	/** this could be done in other ways I think... */
	typedef voro::container NORMAL_STORAGE;

	/** here we store presets */
	std::vector<voro::wall*> presets;

	virtual void makeBox(WEIGHTED_STORAGE*& output, int divx, int divy, int divz) = 0;
	virtual void makeBox(NORMAL_STORAGE*& output, int divx, int divy, int divz) = 0;

	~VoroppGeometry(){ 			DEB3("VoroppGeometry::__destruct__():");
		//destructor...
	}
	/**
	 * add geometry presets
	 * now it only handles those that are derived from voro::wall presets
	 */
	virtual bool addPreset(GeometryPreset* preset){						DEB3("VoroppGeometry::addPreset(): ");
		voro::wall* w = dynamic_cast<voro::wall*>(preset);
		if (w) presets.push_back(w);
		DEB3("VoroppGeometry::addPreset(): preset pointer as wall: " << reinterpret_cast<long long>(w));
		return w;
	}
	/**
	 * add voro++ wall object as a preset
	 */
	bool addWall(voro::wall* wall){
		presets.push_back(wall);
		return true;
	}
protected:
	/** TODO: why this is protected? */
	template<typename StorageSpec>
	void applyPresets(StorageSpec& storage){
		std::vector<voro::wall*>::iterator pres_it = presets.begin();
		std::vector<voro::wall*>::iterator pres_end = presets.end();
		for (; pres_it != pres_end; ++pres_it){
			storage.add_wall(*(*pres_it));
		}
	}
};
/**
 * This class represents one of the possible geometries for VoroAtomBank
 * You can use this class
 *
 * VoroGeometry is a class that store information about particle positions
 * and eventually other data associated to it to enhance computational abilities
 */
class OrthogonalGeometry : public VoroppGeometry{
public:
private:
	double x1;
	double y1;
	double z1;
	double x2;
	double y2;
	double z2;
	bool x_periodic;
	bool y_periodic;
	bool z_periodic;
public:
	std::string show(){ return "dmga::geometry::OrthogonalGeometry"; }

	typedef VoroppGeometry super;
	typedef typename super::WEIGHTED_STORAGE WEIGHTED_STORAGE;
	typedef typename super::NORMAL_STORAGE NORMAL_STORAGE;
	/**
	 * orthogonal geometry is defined by box size (3 dimentions) and
	 * definitions of periodicity - 1 for each direction
	 * so OrthogonalGeometry(1,1,1,true,true,true) is a unit box with all
	 * directions periodic
	 *
	 * you must be carefull with not periodic geometry, when you have particles
	 * outside the nonperiodic direction, as this often means that container will
	 * not store them, @see container.
	 */
	OrthogonalGeometry(double x_container_size = 1.0,
						double y_container_size = 1.0,
						double z_container_size = 1.0,
						bool x_periodic = false,
						bool y_periodic = false,
						bool z_periodic = false):
			x1(0),
			y1(0),
			z1(0),
			x2(x_container_size),
			y2(y_container_size),
			z2(z_container_size),
			x_periodic(x_periodic),
			y_periodic(y_periodic),
			z_periodic(z_periodic){										DEB3("OrthogonalGeometry::__construct__(): x=" << x_container_size << ", y=" << y_container_size << ", z=" << z_container_size << ", xper=" << x_periodic << ", yper=" << y_periodic << ", zper=" << z_periodic);
	}
	/**
	 * second constructor - useful for some kinds of geometries
	 *
	 * you must be carefull with not periodic geometry, when you have particles
	 * outside the nonperiodic direction, as this often means that container will
	 * not store them, @see container.
	 */
	OrthogonalGeometry(double x1, double y1, double z1,
						double x2, double y2, double z2,
						bool x_periodic = false,
						bool y_periodic = false,
						bool z_periodic = false):
			x1(x1), y1(y1), z1(z1),
			x2(x2), y2(y2), z2(z2),
			x_periodic(x_periodic),
			y_periodic(y_periodic),
			z_periodic(z_periodic){
		DEB3("OrthogonalGeometry::__construct__(): x1=(" << x1 << "," << x2 << "), y=(" << y1 << "," << y2 << "), z=(" << z1 << "," << z2 << "), xper=" << x_periodic << ", yper=" << y_periodic << ", zper=" << z_periodic);
	}
	/**
	 * reinitializes geometry - useful in some cases
	 *
	 * TODO: maybe we should provide setters for individual values
	 */
	void reinit(double x1, double y1, double z1,
				double x2, double y2, double z2,
				bool x_periodic = false,
				bool y_periodic = false,
				bool z_periodic = false){
		this->x1 = (x1); this->y1 = (y1); this->z1 = (z1);
		this->x2 = (x2); this->y2 = (y2); this->z2 = (z2);
		this->x_periodic = (x_periodic);
		this->y_periodic = (y_periodic);
		this->z_periodic = (z_periodic);
	}
	/** returns size */
	double xSize(){
		return x2 - x1;
	}
	/** returns size */
	double ySize(){
		return y2 - y1;
	}
	/** returns size */
	double zSize(){
		return z2 - z1;
	}
	/** returns if it is periodic */
	double xPeriodic(){
		return x_periodic;
	}
	/** returns if it is periodic */
	double yPeriodic(){
		return y_periodic;
	}
	/** returns if it is periodic */
	double zPeriodic(){
		return z_periodic;
	}
	/**
	 * return bounding box coordinates
	 */
	virtual void boundingBox(double &x1,
							double &y1,
							double &z1,
							double &x2,
							double &y2,
							double &z2){
		x1 = this->x1;
		y1 = this->y1;
		z1 = this->z1;
		x2 = this->x2;
		y2 = this->y2;
		z2 = this->z2;
	}
	/**
	 * creates a new voro++ container for further use.
	 * it is used in our Container
	 * divx, divy, divz are the divisions of container into cubes, which speeds up the computations (see voro++)
	 */
	void makeBox(WEIGHTED_STORAGE*& output, int divx = 1, int divy = 1, int divz = 1){
		output = new WEIGHTED_STORAGE(this->x1, this->x2,
					  this->y1, this->y2,
					  this->z1, this->z2,
					  divx, divy, divz, //for debug only... TODO: estimate?
					  xPeriodic(), yPeriodic(), zPeriodic(),
					  8);
		this->applyPresets(*output);
	}
	/**
	 * creates a new voro++ container for further use.
	 * it is used in our Container
	 * divx, divy, divz are the divisions of container into cubes, which speeds up the computations (see voro++)
	 */
	void makeBox(NORMAL_STORAGE*& output, int divx = 1, int divy = 1, int divz = 1){
		output = new NORMAL_STORAGE(this->x1, this->x2,
					  this->y1, this->y2,
					  this->z1, this->z2,
					  divx, divy, divz, //for debug only... TODO: estimate?
					  xPeriodic(), yPeriodic(), zPeriodic(),
					  8);
		this->applyPresets(*output);
	}

	virtual ~OrthogonalGeometry(){									DEB3("OrthogonalGeometry::__destruct__():");
	}
};

class CastXYZPlaneGeometry : public OrthogonalGeometry {
public:
	//TODO: poki co zaimplementowane tylko w wersji python
};

class CastPlaneGeometry : public OrthogonalGeometry {
public:
	//TODO: poki co zaimplementowane tylko w wersji python
};


/**
 * this class store information needed for non-orthogonal box
 * geometry, that can be used in some Molecular Simulations
 *
 * it creates a pararellogram by definig three vectors.
 * you can also define periodicity.
 */
class ParallelGeometry : public VoroppGeometry{
public:
	//TODO:
	virtual void boundingBox(double &x1,
							double &y1,
							double &z1,
							double &x2,
							double &y2,
							double &z2){
		x1 = 0.0;
		y1 = 0.0;
		z1 = 0.0;
		x2 = 1.0;
		y2 = 1.0;
		z2 = 1.0;
	}

	virtual void makeBox(WEIGHTED_STORAGE*& output){
		output = new WEIGHTED_STORAGE(0, 1.0,
								  0, 1.0,
								  0, 1.0,
								  1,1,1, //for debug only... TODO: estimate?
								  false, false, false,
								  8);
	}

	virtual void makeBox(NORMAL_STORAGE*& output){
		output = new NORMAL_STORAGE(0, 1.0,
								  0, 1.0,
								  0, 1.0,
								  1,1,1, //for debug only... TODO: estimate?
								  false, false, false,
								  8);
	}

	virtual ~ParallelGeometry(){											DEB3("ParallelGeometry::__destruct__():");
	}
};

}//namespace geometry

}//namespace dmga

#endif /* GEOMETRY_H_ */
