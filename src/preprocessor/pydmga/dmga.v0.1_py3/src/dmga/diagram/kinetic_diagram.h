/*
 * kinetic_diagram.h
 *
 *  Created on: 15-01-2013
 *      Author: Robson
 */

#ifndef KINETIC_DIAGRAM_H_
#define KINETIC_DIAGRAM_H_

#include "diagram.hpp"

namespace dmga{

namespace diagram{

/**
 * this is basic implementation for Diagram of a moving particles
 *
 * it does move by recomputing all cells that were selected
 *
 * this is handy when we select only few cells and recomputing them will be more efficient
 * than managing some Kinetic Data Structure (KDS) that tries to update cells instead of recomputing.
 */
template<typename ContainerSpec, typename TrajectorySpec, typename CellSpec>
class VoroppBasicKineticDiagram : public virtual VoroppDiagram<ContainerSpec, CellSpec>, public virtual GenericKineticDiagram{
public:
	typedef VoroppDiagram<ContainerSpec, CellSpec> super;
	typedef ContainerSpec ContainerType;
	typedef TrajectorySpec TrajectoryType;
	typedef typename super::GeometryType GeometryType;
	typedef typename super::ParticleType ParticleType;
	typedef typename super::KeyType KeyType;
	typedef typename super::CellType CellType;
	typedef typename super::ResultSetType ResultSetType;
	typedef typename super::SubsetType SubsetType;
	typedef typename super::CacheType CacheType;
	/**
	 * the trajectory of the moving particless
	 * class that allows for updating particle positions
	 */
	TrajectoryType& trajectory;
	using super::container;
	using super::cache;
	/**
	 * creates new KineticDiagram for given container and Trajectory
	 */
	VoroppBasicKineticDiagram(ContainerType& container, TrajectoryType& trajectory) :
			super(container),
			trajectory(trajectory){ 																				DEB3("VoroppBasicKineticDiagram::__construct__()");
	}
	/**
	 * advance to the next frame
	 */
	bool advance(){																									DEB3("VoroppBasicKineticDiagram::advance()");
		bool res = trajectory.advance();
		if (res){
			container.update(trajectory);
			container.updateVoroppStorage(this->storage, this->location);
		}
		return res;
	}
	/**
	 * update all cells stored in cache
	 */
	void updateCells(){																								DEB3("VoroppBasicKineticDiagram::updateCells()");
		auto cell_it = cache.begin();
		auto cell_end = cache.end();
		for (; cell_it < cell_end; ++cell_it){
			//this will do all the work, provided container.update(trajectory)
			//take care of the relocating points as they move.
			//we need not to assign the return value, because ResultSetType
			//holds the same pointers to the cells as cache do.
			if (*cell_it){
				this->computeCell( (*cell_it)->number );
			}
		}
	}
};

/**
 * this is an iplementation of a Kinetic Data Structure as described by Guibas, edelsbrunner and other people
 * we do this computation not on Delaunay Triangulation but directly on Voronoi Diagram
 *
 * this implementation is supossed to work a little beter for computing full Voronoi Diagram (all cells)
 * for now we assume that this is the case - because KDS needs all information on vertices. If not all cells
 * are computed then this implementation will not work...
 *
 * some work is needed for the boundary conditions and periodic containers. Some more work is needed for
 * Presets.
 */
template<typename ContainerSpec, typename TrajectorySpec, typename CellSpec = dmga::model::Cell>
class VoroppClassicKineticDiagram : public virtual VoroppDiagram<ContainerSpec, CellSpec>, public virtual GenericKineticDiagram{
public:
	typedef VoroppDiagram<ContainerSpec, CellSpec> super;
	typedef ContainerSpec ContainerType;
	typedef TrajectorySpec TrajectoryType;
	typedef typename super::GeometryType GeometryType;
	typedef typename super::ParticleType ParticleType;
	typedef typename super::KeyType KeyType;
	typedef typename super::CellType CellType;
	typedef typename super::ResultSetType ResultSetType;
	typedef typename super::SubsetType SubsetType;
	typedef typename super::CacheType CacheType;

	TrajectoryType& trajectory;
	using super::container;
	using super::cache;
	/** newton base for interpolation of trajectory polynomials */
	NewtonBase newton;

	class PreEvent{
	public:
	};

	class ClassicPreEvent : public PreEvent{
	public:
		int s[5];
		ClassicPreEvent(int* data){
			s[0] = data[0];
			s[1] = data[1];
			s[2] = data[2];
			s[3] = data[3];
			s[4] = data[4];
		}
		friend bool operator==(const ClassicPreEvent& a, const ClassicPreEvent&b){
			return a.s[0] == b.s[0] && a.s[1] == b.s[1] && a.s[2] == b.s[2] && a.s[3] == b.s[3] && a.s[4] == b.s[4];
		}
		class PtrComp{
		public:
			bool operator()(const ClassicPreEvent* a, const ClassicPreEvent* b){
				return
						a->s[0] < a->s[0] ||
						(a->s[0] == a->s[0] && a->s[1] < a->s[1]) ||
						(a->s[0] == a->s[0] && a->s[1] == a->s[1] && a->s[2] < a->s[2]) ||
						(a->s[0] == a->s[0] && a->s[1] == a->s[1] && a->s[2] == a->s[2] && a->s[3] < a->s[3]) ||
						(a->s[0] == a->s[0] && a->s[1] == a->s[1] && a->s[2] == a->s[2] && a->s[3] == a->s[3] && a->s[4] < a->s[4]);
			}
		};
	};

	class WallPreEvent : public PreEvent{
	public:
	};

	class PeriodicPreEvent : public PreEvent{
	public:
	};

	class Event{
	public:
		double time;
		class PtrComp{
		public:
			bool operator()(const Event* a, const Event* b){
				return a->time < b->time;
			}
		};
	};

	class ClassicEvent : public Event{
	public:
	};

	class WallEvent : public Event{
	public:
	};

	class PeriodicEvent : public Event{
	public:
	};

	VoroppClassicKineticDiagram(ContainerType& container, TrajectoryType& trajectory) :
			super(container),
			trajectory(trajectory){																						DEB3("VoroppClassicKineticDiagram::__construct__()");
		double interpolation_points[] = {-3.0, -2.0, -1.0, 0.0, 1.0, 2.0};
		newton = NewtonBase(5, interpolation_points);
	}

	bool advance(){																										DEB3("VoroppClassicKineticDiagram::advance()");
		bool res = trajectory.advance();
		if (res){
			//container.update(trajectory);
			//you need to build the KDS event queue
			//and update cells accordingly to it later
			//TODO: build

//			foreach delaunay simplex s do
//				foreach neighbour n of s do
//					create event test T in KDS queue based on 5 vertices of s and n
//					compute positive roots t_k of T <= 1.0 //preposition: there should be only one always
//					put event(n, s, t_k) in the queue

//			zmiana
//			foreach pair of delaunay simplexes n, s that share side
//					create event test T in KDS queue based on 5 vertices of s and n
//					compute positive roots t_k of T <= 1.0 //preposition: there should be only one always
//					put event(n, s, t_k) in the queue

//			the above can be done by
//			foreach cell c in Voronoi do
//				foreach edge e in c do
//					use as s the simplex defined by 3 walls meeting at vertex edge.u
//					use as n the simplex defined by 3 walls meeting at vertex edge.v
//					create event test T in KDS queue based on 5 vertices of s and n
//					compute positive roots t_k of T <= 1.0 //preposition: there should be only one always
//					put event(n, s, t_k) in the queue
//			keep in queue only those events that were created twice
//			mark all particles that formed simplexes deleted in previous step as 'need to rebuild'

			int a, b, s[5];
			int vcw, jcw, ucw, icw;
			std::vector<ClassicPreEvent*> classic_pre_list;
			std::vector<WallPreEvent*> wall_pre_list;
			std::vector<PeriodicPreEvent*> periodic_pre_list;

			//iterate over all cells...
			auto cells_it = cache.begin();
			auto cells_end = cache.end();
			for (; cells_it != cells_end; ++cells_it){
				//TODO: test if all vertex has degree 3 - otherwise mark as 'need reconstruct' - not general position
				CellType& cell = **cells_it;
				DEB2("Ur in da cell " << cell.number);
				auto edges_it = cell.edges.begin();
				auto edges_end = cell.edges.end();
				//iterate over all edges in cell
				for (; edges_it != edges_end; ++edges_it){
					//each edge is a base for one event
					//we need to distinguish 'classic' events, 'wall' events and 'periodic' events

					//you need to perform only one Event creation
					//we only look once per edge, so we ignore backward edges
					if (edges_it.u < edges_it.v){
						if (cell.getEdgesCountAt(edges_it.v) != 3) {DEB2("bad edges count at v"); continue; /*NOT GENERAL*/}
						if (cell.getEdgesCountAt(edges_it.u) != 3) {DEB2("bad edges count at u"); continue; /*NOT GENERAL*/}
						if (cell.getNeighbourId(edges_it.v, 0) < 0 || cell.getNeighbourId(edges_it.u, 0) < 0) { continue; /*WALL EVENT*/}
						if (cell.getNeighbourId(edges_it.v, 1) < 0 || cell.getNeighbourId(edges_it.u, 1) < 0) { continue; /*WALL EVENT*/}
						if (cell.getNeighbourId(edges_it.v, 2) < 0 || cell.getNeighbourId(edges_it.u, 2) < 0) { continue; /*WALL EVENT*/}

						//now take all neighbouring particles that creates vertices at the ends of a given edge
						//this creates two delaunay simplices with common wall
						//we will store this wall is s[1..3] as sorted list
						//we will put the two other particles numbers in s[0] and s[5]
						//min in s[0] and max in s[5]
						//this way we will have a unique description of the event, regardles
						//of which cell we used and which edge to create it as long as it is essentialy the same edge
						s[1] = cell.number;
						s[2] = cell.getNeighbourId(edges_it.u, edges_it.i);
						s[3] = cell.getNeighbourId(edges_it.v, edges_it.j);
						voro::edges::cw(cell.voropp_cell, edges_it.u, edges_it.i, ucw, icw);
						voro::edges::cw(cell.voropp_cell, edges_it.v, edges_it.j, vcw, jcw);
						a = cell.getNeighbourId(ucw, icw);
						b = cell.getNeighbourId(vcw, jcw);
						s[0] = std::min(a, b);
						s[5] = std::max(a, b);
						std::sort(s + 1, s + 4); //TODO: here you can switch to optimal sorting of 3 elements

						DEB2("Event Base: " << s[0] << " " << s[1] << " " << s[2] << " " << s[3] << " " << s[4]);

						//TODO: THINK: how to identify periodic events
						ClassicPreEvent* pre_event = new ClassicPreEvent(s);
						classic_pre_list.push_back(pre_event);
					}
				}
			}

			std::vector<Event*> classic_event_list = prepareEvents(classic_pre_list);

//			THINK: jak zrobi� eventy na scianie
//			THINK: jak zrobi� eventy dla periodic
		}
		return res;
	}

	template<typename PreEventSpec>
	std::vector<Event*> prepareEvents(std::vector<PreEventSpec*>& pre_events){										DEB3("VoroppClassicKineticDiagram::prepareEvents()");
		std::vector<Event*> result;

		typedef typename PreEventSpec::PtrComp PreCompare;
		std::sort(pre_events.begin(), pre_events.end(), PreCompare());

		auto it = pre_events.begin();
		auto end = pre_events.end();
		PreEventSpec* current = 0;
		int count;
		if (it != end) {
			current = *it;
			count = 1;
		}
		for (; it != end; ++it){
			if (*current == **it){
				++count;
			}else{
				if (count == 3){
					//add this to result - it is good, general event,
					//Event* event = current->createEvent(container, trajectory);
					//result.push_back(event);
					count = 1;
					current = *it;
				}else{
					//this is not general event too many cells meet at the defining edge
					//TODO: mark cells as 'to recompute'
				}
			}
		}
		typedef typename Event::PtrComp TimeSortComp;
		std::sort(result.begin(), result.end(), TimeSortComp());
		return result;
	}

	void updateCells(){																											DEB3("VoroppClassicKineticDiagram::updateCells()");
//		while queue not empty do
//			event = queue.front();
//			resolve event, there should be three new simplices in triangulation, s1, s2, s3
//			delete all events from queue that are based on event.n or event.s
//			add to queue events based on s1, s2, s3 and their neighbours - this can be done using only s1, s2, s3 ?
	}
};

/**
 * this is implementation of a Kinetic Data Structure that is simple but we belive it will
 * give a decent boost for the algorithm. It takes into consideration assumption about small
 * moves of the particles, so offten you not need to recompute the cells if neighbouring particles
 * will remain sufficiently close to the particle. If not we recompute the whole cell.
 *
 * This need to be tested if will give some speed boost.
 *
 * This may work better if we do some more work and in recompute can only work with
 * some neighbouring particles - so we need to determine when particles become neighbours.
 */
template<typename ContainerSpec, typename TrajectorySpec, typename CellSpec>
class VoroppKineticDiagram : public virtual VoroppDiagram<ContainerSpec, CellSpec>, public virtual GenericKineticDiagram{
public:
	typedef VoroppDiagram<ContainerSpec, CellSpec> super;
	typedef ContainerSpec ContainerType;
	typedef TrajectorySpec TrajectoryType;
	typedef typename super::GeometryType GeometryType;
	typedef typename super::ParticleType ParticleType;
	typedef typename super::KeyType KeyType;
	typedef typename super::CellType CellType;
	typedef typename super::ResultSetType ResultSetType;
	typedef typename super::SubsetType SubsetType;
	typedef typename super::CacheType CacheType;

	TrajectoryType& trajectory;
	using super::container;

	VoroppKineticDiagram(ContainerType& container, TrajectoryType& trajectory) :
			super(container),
			trajectory(trajectory){																				DEB3("VoroppKineticDiagram::__construct__()");
	}

	bool advance(){																								DEB3("VoroppKineticDiagram::advance()");
		bool res = trajectory.advance();
		if (res){
			container.update(trajectory);
		}
		return res;
	}

	void updateCells(ResultSetType& cells){																		DEB3("VoroppKineticDiagram::updateCells()");
		auto cell_it = cells.begin();
		auto cell_end = cells.end();
		for (; cell_it < cell_end; ++cell_it){
			//you need to build the KDS event queue
			//and update cells accordingly to it later
			//TODO: build
		}
	}
};

} //namespace diagram

} //namespace dmga


#endif /* KINETIC_DIAGRAM_H_ */
