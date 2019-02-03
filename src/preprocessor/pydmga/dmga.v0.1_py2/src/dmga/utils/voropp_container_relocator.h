/*
 * voropp_container_relocator.h
 *
 *  Created on: Aug 16, 2013
 *      Author: robson
 */

#ifndef VOROPP_CONTAINER_RELOCATOR_H_
#define VOROPP_CONTAINER_RELOCATOR_H_

#include <base.h>
#include <3rd/voro/voro++.hh>
#include <exceptions/exceptions.hpp>
#include <algebra/algebra.hpp>

namespace voro{

/**
 * This class is a kind of wrap on the voro::container
 * that allows updating of positions of perticles
 * inside the container, it maintain the structure of the container
 * (computational blocks). It is supposed to be faster
 * than recreating the whole container once again.
 */
//template<typename VoroppContainerSpec>
typedef container_poly VoroppContainerSpec; //for easier writing of code... (autocomplete)
class container_relocator : public VoroppContainerSpec{
public:
	class NXYZ{
	public:
		int n;
		double raw[4];
		NXYZ(): n(-1) {}
		NXYZ(int n, double x, double y, double z): n(n){
			raw[0] = x; raw[1] = y; raw[2] = z;
		}
		NXYZ(int n, double x, double y, double z, double r): n(n){
			raw[0] = x; raw[1] = y; raw[2] = z; raw[3] = r;
		}
		double& x(){return raw[0];}
		double& y(){return raw[1];}
		double& z(){return raw[2];}
		double& r(){return raw[3];}
	};

	typedef VoroppContainerSpec super;
	/** here for each ijk block we hold list of q which are going to be removed */
	std::vector<std::vector<int> > remove_queue;
	/** here for each ijk block we hold list of raws to be inserted */
	std::vector<std::vector<NXYZ> > insert_queue;

	container_relocator(double ax_,double bx_,double ay_,double by_,double az_,double bz_,
						int nx_,int ny_,int nz_,bool xperiodic_,bool yperiodic_,bool zperiodic_,
						int init_mem) :
			super(ax_, bx_, ay_, by_, az_, bz_, nx_, ny_, nz_, xperiodic_, yperiodic_, zperiodic_, init_mem),
			remove_queue(nx_ * ny_ * nz_),
			insert_queue(nx_ * ny_ * nz_){
	}

	/**
	 * start the batch set - set of positions in multiple particles (best: all)
	 * cleans necessary data structures
	 */
	void batchset_start(){
		auto rem_it = remove_queue.begin();
		auto rem_end = remove_queue.end();
		for (; rem_it != rem_end; ++rem_it) (*rem_it).clear();
		auto ins_it = insert_queue.begin();
		auto ins_end = insert_queue.end();
		for (; ins_it != ins_end; ++ins_it) (*ins_it).clear();
	}

	/**
	 * do all enqueued events - move particles to their final destinations
	 * scans through necessary data structures and swaps things
	 */
	template<typename LocationSpec>
	void batchset_finalize(LocationSpec& some_location){
		for (int ijk = 0; ijk < nxyz; ++ijk){ //process each block
			//sort(remove_queue[ijk].begin(), remove_queue[ijk].end());
			auto ins_it = insert_queue[ijk].begin();
			auto ins_end = insert_queue[ijk].end();
			auto rem_it = remove_queue[ijk].begin();
			auto rem_end = remove_queue[ijk].end();
			int co_resize = ((int)insert_queue[ijk].size() - (int)remove_queue[ijk].size()); //this is the number of added (positive) removed (negative) elements, we need to do (co[ijk] += co_resize) after all operations
			//TODO: sprawdzić przed rozpoczęciem, czy co[ijk] + (insert.size() - remove.size()) nie doprowadzi do
			//      potrzeby zmniejszaenia / zwiekszenia alokacji - jesli tak to najpierw można wykonać alokację,
			//      i do tego przyspieszyć proces wstawiania - podczas kopiowania wywalaj te do remove'a
			//      a potem na końcu dostaw te do dodania.
			while (ins_it != ins_end && rem_it != rem_end){ //while we have both to insert and to remove...
				int q_rem = *rem_it;
				double* new_raw = (*ins_it).raw;
				int new_id = (*ins_it).n;
				DEB("relocating " << new_id << " to " << ijk << " " << q_rem);
				double* pp = dmga::utils::getContainerCoordsRaw(*this, ijk, q_rem); //here we will copy point information
				double* pp_end = pp + ps;
				while (pp != pp_end) *(pp++) = *(new_raw++);
				//now we need to update voro++ internal container architecture
				id[ijk][q_rem] = new_id; //update id of the particle that cames
				some_location[new_id].ijk = ijk;  //update location
				some_location[new_id].q = q_rem; //update location
				//thats all - we do not need to touch co or mem, because we deleted one and added one = everything ok
				++ins_it; ++rem_it; //proceed
			}
			//now we have one or the another but not both
			if (rem_it != rem_end){
				int q_curr = *rem_it; //this is the first empty (marked to remove) place
				int q = q_curr + 1; //this probably is the first thing to shift
				++rem_it; //proceed wit rem_it to the next station
				int q_next = (rem_it == rem_end ? (co[ijk] + 2) : *rem_it); //here we will hold next thing to remove
				while (q != co[ijk]){ //we iterate over all remaining particles
					 if (q == q_next){ //q is at the thing for remove
						 ++rem_it; //skip this and move next to the next thing to remove
						 q_next = (rem_it == rem_end ? (co[ijk] + 2) : *rem_it);
					 }else{	//q is a real thing, copy data from q to q_curr
						 int n = id[ijk][q];
						 id[ijk][q_curr] = n; //copy id
						 some_location[n].ijk = ijk;  //update location
						 some_location[n].q = q_curr; //update location
						 DEB("shifting " << n << " to " << ijk << " " << q_curr);
						 //copy position
						 double* pp_curr = dmga::utils::getContainerCoordsRaw(*this, ijk, q_curr);
						 double* pp_oth = dmga::utils::getContainerCoordsRaw(*this, ijk, q);
						 double* pp_oth_end = pp_oth + ps;
						 while (pp_oth != pp_oth_end) *(pp_curr++) = *(pp_oth++);
						 ++q_curr; //move q_curr to next, because this space is now used
					 }
					 ++q; //shift next
				}
				co[ijk] += co_resize; //set the real size of co[ijk]
			}
			while (ins_it != ins_end){
				//TODO: put this as in put_relocate... (use add_particle_memory if necessary)
				if (co[ijk] == mem[ijk]) add_particle_memory(ijk);
				double* new_raw = (*ins_it).raw;
				id[ijk][co[ijk]] = (*ins_it).n;
				DEB("inserting " << (*ins_it).n << " to " << ijk << " " << co[ijk]);
				some_location[(*ins_it).n].ijk = ijk;   //update location
				some_location[(*ins_it).n].q = co[ijk];	//update location
				double *pp = dmga::utils::getContainerCoordsRaw(*this, ijk, co[ijk]);
				double* pp_end = pp + ps;
				while (pp != pp_end) *(pp++) = *(new_raw++);
				++co[ijk]; //expand
				++ins_it; //proceed
			}
			//we do not add co_resize after inserting, because we increased it along the way... see co[ijk]++
		}
	}

	bool batchset_enqueue(int ijk, int q, double x, double y, double z){
		int old_ijk = ijk;
		int old_q = q;
		int old_id = id[old_ijk][old_q];
		if(!this->put_remap(ijk, x, y, z)){
			//we have moved outside the container - remove old particle from container
			DEB("enqueue remove [outside] " << old_ijk << " " << old_q << " [id = " << id[old_ijk][old_q] << "]");
			remove_queue[old_ijk].push_back(old_q);
			return false;
		}
		if (ijk == old_ijk){
			//this is simple, we have stayed in last box, so we need only to change the coordinates, and thats all!
			DEB("simple relocate " << ijk << " " << q << " [id = " << id[ijk][q] << "]");
			double *pp = dmga::utils::getContainerCoordsRaw(*this, ijk, q);
			pp[0] = x; pp[1] = y; pp[2] = z;
		}else{
			//this is hard, we have moved to other box, mark other as to delete, and add new to queues for addition
			//things planned here will be executed in batchset_finalize() for better performance
			remove_queue[old_ijk].push_back(old_q);
			DEB("enqueue remove " << old_ijk << " " << old_q << " [id = " << id[old_ijk][old_q] << "]");
			DEB("enqueue insert " << old_ijk << " " << old_q << " to " << ijk << " [id = " << id[old_ijk][old_q] << "]");
			if (ps == 4){
				double r = dmga::utils::getContainerCoordsRaw(*this, old_ijk, old_q)[3]; //copy old r
				insert_queue[ijk].push_back(NXYZ(old_id, x, y, z, r));
			}else{
				insert_queue[ijk].push_back(NXYZ(old_id, x, y, z));
			}
		}
		return true;
	}
};

} // namespace voro


#endif /* VOROPP_CONTAINER_RELOCATOR_H_ */
