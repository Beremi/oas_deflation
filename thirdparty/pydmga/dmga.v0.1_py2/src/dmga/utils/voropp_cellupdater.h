/*
 * voropp_container.h
 *
 *  Created on: 23-10-2012
 *      Author: Robson
 */

#ifndef VOROPP_CELLUPDATER_H_
#define VOROPP_CELLUPDATER_H_

#include "utils.hpp"
#include "algebra/algebra.hpp"

namespace dmga{

/**
 * this is old class for updating voro++ cells to new coordinates - it may produce bad results!
 */
class CellUpdater{
public:
	voro::container& con;
	double** pos_by_id;

	CellUpdater(voro::container& con): con(con){
		int maxId = -1;
		for (int ijk = 0; ijk < con.nxyz; ijk++){
			for (int q = 0; q < con.co[ijk]; q++){
				if (maxId < con.id[ijk][q]){
					maxId = con.id[ijk][q];
				}
			}
		}
		pos_by_id = new double* [maxId + 1];
		for (int ijk = 0; ijk < con.nxyz; ijk++){
			for (int q = 0; q < con.co[ijk]; q++){
				pos_by_id[con.id[ijk][q]] = con.p[ijk] + 3 * q;
			}
		}
	}
	double* pos(int id){
		return pos_by_id[id];
	}
	/**
	 * update the position of single particle
	 * */
	void update(int n, double x, double y, double z){
		throw dmga::exceptions::NotImplementedYet("CellUpdater::update()");
	}
	/**
	 * update position of all particles in the box
	 * */
	void update(FILE *fp=stdin){
		throw dmga::exceptions::NotImplementedYet("CellUpdater::update()");
	}
	/**
	 * update position of all particles in the box
	 * */
	void update(voro::particle_order &vo, FILE *fp=stdin){
		throw dmga::exceptions::NotImplementedYet("CellUpdater::update()");
	}
	/**
	 * update position of all particles in the box
	 * */
	inline void update(const char* filename) {
		throw dmga::exceptions::NotImplementedYet("CellUpdater::update()");
	}
	/**
	 * update position of all particles in the box
	 * */
	inline void update(voro::particle_order &vo, const char* filename) {
		throw dmga::exceptions::NotImplementedYet("CellUpdater::update()");
	}

	/**
	 * update cell given current coordinates
	 * */
	template<class v_cell,class c_loop>
	inline bool update(v_cell &c, c_loop &vl) {
		bool can_simple_update = false;

		/* try to update cell without doing all the job */

		if (!can_simple_update){
			return con.compute_cell(c, vl);
		}
	}
	/**
	 * update cell given current coordinates
	 * */
	template<class v_cell>
	inline bool update(v_cell &c, int ijk, int q) {
		bool can_simple_update = false;

		/* try to update cell without doing all the job */

		if (!can_simple_update){
			return con.compute_cell(c, ijk, q);
		}
	}
};



class CellPolyUpdater{
public:
	voro::container_poly& con;
	double** pos_by_id;

	CellPolyUpdater(voro::container_poly& con) : con(con){
		int maxId = -1;
		for (int ijk = 0; ijk < con.nxyz; ijk++){
			for (int q = 0; q < con.co[ijk]; q++){
				if (maxId < con.id[ijk][q]){
					maxId = con.id[ijk][q];
				}
			}
		}
		pos_by_id = new double* [maxId + 1];
		for (int ijk = 0; ijk < con.nxyz; ijk++){
			for (int q = 0; q < con.co[ijk]; q++){
				pos_by_id[con.id[ijk][q]] = con.p[ijk] + 4 * q;
			}
		}
	}
	/** this isn't safe until recomputePosById isn't called */
	double* pos(int id){
		return pos_by_id[id];
	}

	/**
	 * update the position of single particle
	 * */
	void update(int n, double x, double y, double z){
		throw dmga::exceptions::NotImplementedYet("updatable_container::update()");
	}
	/**
	 * update position of all particles in the box
	 * */
	void update(FILE *fp=stdin){
		throw dmga::exceptions::NotImplementedYet("updatable_container::update()");
	}
	/**
	 * update position of all particles in the box
	 * */
	void update(voro::particle_order &vo, FILE *fp=stdin){
		throw dmga::exceptions::NotImplementedYet("updatable_container::update()");
	}
	/**
	 * update position of all particles in the box
	 * */
	void update(const char* filename) {
		throw dmga::exceptions::NotImplementedYet("updatable_container::update()");
	}
	/**
	 * update position of all particles in the box
	 * */
	void update(voro::particle_order &vo, const char* filename) {
		throw dmga::exceptions::NotImplementedYet("updatable_container::update()");
	}

	/**
	 * update cell given current coordinates
	 * */
	template<class v_cell,class c_loop>
	inline bool update_cell(v_cell &c, c_loop &vl) {
		return update_cell(c, vl.ijk, vl.q);
	}
	/**
	 * update cell given current coordinates
	 * */
	template<class v_cell>
	inline bool update_cell(v_cell &c, int ijk, int q) {
		bool can_simple_update = false;

		/* try to update cell without doing all the job */
		int i, j;
		algebra::Plane3D pl[3];
		for(i=0;i<c.p;i++){
			//compute each point
			if (c.nu[i] > 3){ //for now only info
				DEB3("updateAllCells(): Redundant edges!" << "\n");
			}
			if (c.nu[i] < 3){ //for now only info
				throw throw dmga::exceptions::NotImplementedYet("updateAllCells(): not sufficient information, should never happen");
			}

			bool skip = false;
			for (j=0;j<3;j++){
				if (c.ne[i][j] < 0){
					skip = true;
					continue;
				}

				pl[j] = algebra::ballsToPowerPlane(con.p[ijk] + 4 * q, pos(c.ne[i][j]));
			}

			if (skip) continue;

			//this 2.0 * is because of the way in which Voro++ computes cells (2 times bigger)
			algebra::Vector3D point = 2.0 * planes3Intersection(pl[0], pl[1], pl[2]);
			for(j = 0; j < 3; j++){
				c.pts[i * 3 + j] = point[j];
			}
		}

		if (!can_simple_update){
			return con.compute_cell(c, ijk, q);
		}
	}
};

}//namespace dmga


#endif /* VOROPP_CONTAINER_H_ */
