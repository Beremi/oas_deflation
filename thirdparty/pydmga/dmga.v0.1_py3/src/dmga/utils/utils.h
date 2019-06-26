/*
 * utils.h
 *
 *  Created on: 10-11-2012
 *      Author: Robson
 */

#ifndef UTILS_H_
#define UTILS_H_

#include <base.h>
#include <3rd/voro/voro++.hh>
#include <exceptions/exceptions.hpp>
#include <algebra/algebra.hpp>


namespace dmga{

namespace utils{

/**
 * return coords from container
 */
inline double* getContainerCoordsRaw(voro::container_poly& container, int ijk, int q){
	return container.p[ijk] + 4 * q;
}

/**
 * return coords from container
 */
inline double* getContainerCoordsRaw(voro::container& container, int ijk, int q){
	return container.p[ijk] + 3 * q;
}

/**
 * put particle to container
 */
inline void putContainerCoordsRaw(voro::container_poly& container, int id, const double* raw){
	container.put(id, raw[0], raw[1], raw[2], raw[3]);
}

/**
 * put particle to container
 */
inline void putContainerCoordsRaw(voro::container& container, int id, const double* raw){
	container.put(id, raw[0], raw[1], raw[2]);
}

/**
 * put particle to container
 */
inline void putContainerCoordsRaw(voro::container_poly& container, voro::particle_order& order, int id, const double* raw){
	container.put(order, id, raw[0], raw[1], raw[2], raw[3]);
}

/**
 * put particle to container
 */
inline void putContainerCoordsRaw(voro::container& container, voro::particle_order& order, int id, const double* raw){
	container.put(order, id, raw[0], raw[1], raw[2]);
}

/**
 * helper function for safe delete of pointer if not 0 and set to 0 after delete
 */
template<typename T>
void safeDelete(T* &ptr){
	if (ptr){
		delete ptr;
		ptr = 0;
	}
}

/**
 * helper function for safe delete of pointer if not 0 and set to 0 after delete
 */
template<typename T>
void safeDeleteArray(T* &ptr){
	if (ptr){
		delete[] ptr;
		ptr = 0;
	}
}


} //namespace utils

} //namespace dmga


namespace voro{

/** Helper class that can store information on
 * particles in itself nad the parent simultanously
 * holds only one level of hierarchy
 * It;s nice because not depends on internal structer of voro++
 * but uses only their function through inheritance
 */
class multi_particle_order : public particle_order{
	particle_order* orders;
	int orders_count;
	typedef particle_order super;
public:
	/** Adds a record to the order, corresponding to the memory
	 * address of where a particle was placed into the container.
	 * \param[in] ijk the block into which the particle was placed.
	 * \param[in] q the position within the block where the
	 * 		particle was placed. */
	inline void add(int ijk,int q) {
		particle_order* act = orders;
		particle_order* last = orders + orders_count;
		for (; act < last; act++){
			act->add(ijk, q);
		}
	}

	multi_particle_order(int number, particle_order* orders,
						int init_size = init_ordering_size):
		super(init_size), orders(orders), orders_count(number)
	{}
};

class extended_particle_order : public particle_order{
public:
	typedef particle_order super;

	/** returns number of atoms in collection */
	int size(){
		throw dmga::exceptions::NotImplementedYet("extended_particle_order: size() not implemented yet.");
	}
	/**
	 * returns the container coords where Voro++ holds information on
	 * given atom, this should be not used by the users of the library
	 */
	inline void fetchContainerCoords(int i, int& ijk, int& q){
		ijk = o[2*i];
		q = o[2*i+1];
	}
	/**
	 * constructor required by VORO++
	 * @param init_size
	 */
	extended_particle_order(int init_size = init_ordering_size):
		super(init_size)
	{}
};



class fetching_particle_order : public particle_order{
public:
	typedef particle_order super;
	int ijk;
	int q;
	/**
	 * constructor required by VORO++
	 * @param init_size
	 */
	fetching_particle_order(int init_size = init_ordering_size):
		super(2), ijk(-1), q(-1)
	{}

	/** Adds a record to the order, corresponding to the memory
	 * address of where a particle was placed into the container.
	 * \param[in] ijk the block into which the particle was placed.
	 * \param[in] q the position within the block where the
	 * 		particle was placed. */
	inline void add(int ijk, int q) {
		this->ijk = ijk;
		this->q = q;
	}
};


/** Helper class that can store information on atoms
 * (name, id, symbol, molecule) along with the information
 * that VORO++ needs to iterate over this order
 */
template <class AtomSpec>
class dmgalpha_particle_order : public particle_order{
	typedef particle_order super;
public:
	typedef AtomSpec AtomType;
	typedef typename AtomType::AtomPositionType AtomPositionType;
	typedef typename AtomType::AtomCellType AtomCellType;
	typedef typename AtomType::AtomSurfaceType AtomSurfaceType;
	typedef typename AtomType::AtomDescriptionType AtomDescriptionType;

	/** descriptions in the order of putting atoms to this order */
	std::vector<AtomDescriptionType> _atomDescriptions;

	/** returns number of atoms in collection */
	int atomCount(){
		return _atomDescriptions.size();
	}
	/**
	 * get the information about atom
	 * @param i which atom to get
	 * @return description of the atom (name, id, etc)
	 */
	AtomDescriptionType& atomDescription(int i){
		return _atomDescriptions[i];
	}

	/**
	 * get the information about atom
	 * @param i which atom to get
	 * @return description of the atom (name, id, etc)
	 */
	const AtomDescriptionType& atomDescription(int i) const{
		return _atomDescriptions[i];
	}

	/**
	 * returns the container coords where Voro++ holds information on
	 * given atom, this should be not used by the users of the library
	 */
	inline void fetchContainerCoords(int i, int& ijk, int& q){
		ijk = o[2*i];
		q = o[2*i+1];
	}

	//TODO: force add info before add by flags etc
	/**
	 * This routine should be execuded before execution of corresponding
	 * add() function
	 * @param atomDescription
	 */
	inline void addInfo(const AtomDescriptionType& atomDescription){
		_atomDescriptions.push_back(atomDescription);
	}
	/**
	 * constructor required by VORO++
	 * @param init_size
	 */
	dmgalpha_particle_order(int init_size = init_ordering_size):
		super(init_size)
	{}
};

/**
 * renormalize Voro++ cell (original!) that it will hold real coordinates (i.e scaled by 0.5 and moved to (x,y,z))
 */
void normalizeCell(voro::voronoicell_neighbor& pc, double x, double y, double z);
/**
 * denormalize Voro++ cell (previously normalized!) that it holds original coordinates (i.e moved to (0,0,0) and scaled by 2)
 */
void denormalizeCell(voro::voronoicell_neighbor& pc, double x, double y, double z);
/**
 * renormalize Voro++ cell (original!) that it will hold real coordinates (i.e scaled by 0.5 and moved to (x,y,z))
 */
void normalizeCell(voro::voronoicell_neighbor& pc, double* raw);
/**
 * denormalize Voro++ cell (previously normalized!) that it holds original coordinates (i.e moved to (0,0,0) and scaled by 2)
 */
void denormalizeCell(voro::voronoicell_neighbor& pc, double* raw);

namespace edges{
/** returns description of the next edge on the same current side */
inline void next(voro::voronoicell_base& vcell, int v, int j, int& u, int& k);
/** returns description of the previous edge on the current side */
inline void prev(voro::voronoicell_base& vcell, int v, int j, int& z, int& i);
/** returns description of the inverse edge to the current edge of the current side */
inline void inv(voro::voronoicell_base& vcell, int v, int j, int& u, int& l);
/** returns description of the edge clockwise from the current edge on the current side */
inline void cw(voro::voronoicell_base& vcell, int v, int j, int& w, int& m);
/** returns description of the edge counter-clockwise from the current edge on the current side*/
inline void ccw(voro::voronoicell_base& vcell, int v, int j, int& w, int& m);
/** convert v (vertex id) and j (j-th edge at v) to the global number of edge */
int voroToEdge(const int* prefix_nu, int vertex, int jthedge);
}//namespace edges

/** class to describe Edges in Cells */
class EdgeDesc{
public:
	/** the vertex id at which j-th edge is located*/
	int v;
	/** j-th edge of the vertex v */
	int j;
	/** edge number - unique, based on v and j*/
	int edge;
	/** we identify EdgeDesc with the edge number - it's handy in times */
	operator int(){
		return edge;
	}
	/** */
	EdgeDesc(int v = -1, int j = -1, int edge = -1):
			v(v),
			j(j),
			edge(edge){
	}
};

/**
 * this class allows for iteration over voronoicell classes of voro++
 * using only public fields and methods. Moreover this class doesn't change
 * anything in underlying structure, like many methods in voro++ do (eg. they
 * mark visited edges using negative numbers, etc. )
 *
 * it allows for iteration in the way of DCEL data structure
 * t starts from 0-th vertex (first vertex) and its 0-th edge (first edge) and creates
 * all necessary data to allow DCEL structure wandering: using next, prev, inv commands.
 * it also allows for querrying about face that current edge holds
 *
 * Dev-note: Reverse enginerred from Voro++:
 * 	next pair (vertex=u, edge=k) is computed this way:
 * 		u = ed[v][j]
 * 		k = cucle_up(ed[v][j+nu[v]], u)
 * 	previous pair (vertex=z, edge=i) is computed this way:
 * 		z = ed[v][cycle_down(j, v)]
 * 		i = ed[v][nu[v] + cycle_down(j, v)]
 * 	reverse pair (vertex=u, edge=l) is computed this way:
 * 		u = ed[v][j]
 * 		l = ed[v][nu[v]+j]
 */
class CellIterator : public dmga::DMGAObject{
private:
	std::string show(){ return "voro::CellIterator"; }

	/** the cell which we traverse */
	voro::voronoicell_base* vcell;
	/** prefixNu holds offsets in table visited for each edge for each vertex
	 * ie. prefix[i]=k tells that j-th edge of vertex i has its flag stored in visited[k+j]
	 */
	int* prefix_nu;
	/**
	 * @see prefix_nu
	 *
	 * if shared_prefix_nu is set and equal to prefix_nu then the copy constructor will copy
	 * only pointers not whole arrays, moreover destructor will not delete prefix_nu
	 */
	int* shared_prefix_nu;
	/** for each edge we hold here merked/unmarked information */
	std::vector<bool> visited;

	/** voro++ current vertex */
	int v;
	/** voro++ current edge */
	int j;
	/**
	 * @deprecated
	 * for backward compatibility
	 */
	int voroToEdge(int vertex, int jthedge) const{
		return edges::voroToEdge(prefix_nu, vertex, jthedge);
	}
public:
	/** convert v (vertex id) and j (j-th edge at v) to the global number of edge */
	int edge(const EdgeDesc& e) const {
		return voroToEdge(e.v, e.j);
	}
	/** convert v (vertex id) and j (j-th edge at v) to the global number of edge */
	int edge(int vertex, int jthedge) const {
		return voroToEdge(vertex, jthedge);
	}
	/** convert current edge to global number of edge */
	int edge(){
		return voroToEdge(v, j);
	}

	/** with null pointer to pcell, all calls to methods should fail miserably */
	CellIterator();
	/** constructor based on the Voro++ cell class */
	CellIterator(voro::voronoicell_base* pcell, int* precomputed_prefix_nu = 0, int v = 0, int j = 0);
	/** copy constructor - doing only shallow copy, remember! */
	CellIterator(const CellIterator& original);
	/** assign operator */
	CellIterator& operator=(const CellIterator& original);
	/** unmark all edges, prepare for new iteration */
	void reset();
	/** returns description of current edge @see EdgeDesc */
	inline EdgeDesc current() const;
	/** do we seen all edges and vertices? */
	inline bool isFinished();
	/** jump to the next unmarked vertex */
	inline bool jump();
	/** returns description of the next edge on the same current side */
	inline EdgeDesc next() const;
	/** returns description of the previous edge on the current side */
	inline EdgeDesc prev() const;
	/** returns description of the inverse edge to the current edge of the current side */
	inline EdgeDesc inv() const;
	/** returns description of the edge clockwise from the current edge on the current side */
	inline EdgeDesc cw() const;
	/** returns description of the edge counter-clockwise from the current edge on the current side*/
	inline EdgeDesc ccw() const;
	/** move to the next edge on the same side */
	inline CellIterator& forward();
	/** move to the previous edge on the current side */
	inline CellIterator& backward();
	/** switch to the other edge of the current side (also switch current side to the side on the other side of the current edge) */
	inline CellIterator& reverse();
	/** go to the specyfic edge */
	inline CellIterator& go(int v, int j);
	/** go to the specyfic edge */
	inline CellIterator& go(const EdgeDesc& edge);
	/** switch to the clockwise edge on the same side */
	inline CellIterator& rotateCW();
	/** switch to the counter-clockwise edge on the same side */
	inline CellIterator& rotateCCW();
	/** tests if current egde is marked */
	inline bool isMarked() const;
	/** test if given edge is marked */
	inline bool isMarked(int edge) const;
	/** mark given edge */
	inline CellIterator& mark(int edge);
	/** mark current egde */
	inline CellIterator& mark();
	/** unmark given edge */
	inline CellIterator& unmark(int edge);
	/** unmark current edge */
	inline CellIterator& unmark();
	/** returns current vertex at which the current edge is located (i.e. the starting point of current edge) */
	inline dmga::model::Point3D vertex() const;
	/** returns current vertex at which given edge is located (i.e. the starting point of given edge) */
	inline dmga::model::Point3D vertex(const EdgeDesc& edge) const;
	/** returns current vertex at which given edge is located (i.e. the starting point of given edge) */
	inline dmga::model::Point3D vertex(int v) const;
	/** returns current vertex at which the current edge is located (i.e. the starting point of current edge) */
	inline double* vertexRaw() const;
	/** returns current vertex at which given edge is located (i.e. the starting point of given edge) */
	inline double* vertexRaw(const EdgeDesc& edge) const;
	/** returns current vertex at which given edge is located (i.e. the starting point of given edge) */
	inline double* vertexRaw(int v) const;
	/** returns cell object */
	inline voro::voronoicell_base* cell();
	/** returns cell object */
	inline const voro::voronoicell_base* cell() const;
	/**
	 * if this iterator is not shared then it makes it shared.
	 * It returns prefix_nu pointer, you should store it, because you are now responsible for deleting this
	 * if this iterator is already shared then NULL is returned
	 * */
	inline int* makeShared();
	/** destructor */
	~CellIterator();

	bool operator==(const CellIterator& original) const{
		return vcell == original.vcell && v == original.v && j == original.j;
	}
};

class VertexIterator{
public:
	double* current;

	VertexIterator(double* current): current(current){							DEB3("voro::VertexIterator::__construct__(): 1 args");
	}

	double* raw(){																DEB3("voro::VertexIterator::raw()");
		return current;
	}

	dmga::model::Point3D operator*(){											DEB3("voro::VertexIterator::operator*()");
		return dmga::model::Point3D(raw());
	}

	VertexIterator& operator++(){												DEB3("voro::VertexIterator::operator++()");
		current += 3;
		return *this;
	}

	bool operator==(const VertexIterator& original) const{						DEB3("voro::VertexIterator::operator==()");
		return current == original.current;
	}

	bool operator!=(const VertexIterator& original) const{						DEB3("voro::VertexIterator::operator=!=()");
		return current != original.current;
	}
};

class EdgeIterator{
public:
	CellIterator cell_iterator;
	/** next vertex on current side */
	int u;
	/** back edge, pointing to v vertex */
	int i;
	/** current vertex, base of current edge */
	int v;
	/** forward edge, to the next vertex, i.e. u */
	int j;

	static EdgeIterator EDGE_ITERATOR_NULL;

	EdgeIterator() : cell_iterator(){											DEB3("voro::EdgeIterator::__construct__(): 0 args");
		u = i = v = j = 0;
	}

	explicit EdgeIterator(const CellIterator& cell_iterator) :
			cell_iterator(cell_iterator){ 										DEB3("voro::EdgeIterator::__construct__(): 1 arg");
		if (this->cell_iterator.isFinished()){
			this->cell_iterator = CellIterator();
			u = i = v = j = 0;
		}else{
			u = this->cell_iterator.current().v;
			i = this->cell_iterator.current().j;
			v = this->cell_iterator.inv().v;
			j = this->cell_iterator.inv().j;
		}
	}

	explicit EdgeIterator(const CellIterator& cell_iterator, EdgeDesc at_edge, bool marked) :
			cell_iterator(cell_iterator){ 										DEB3("voro::EdgeIterator::__construct__(): 3 arg");
		this->cell_iterator.go(at_edge);
		u = this->cell_iterator.current().v;
		i = this->cell_iterator.current().j;
		v = this->cell_iterator.inv().v;
		j = this->cell_iterator.inv().j;
		if (marked) this->cell_iterator.mark();
	}

	EdgeIterator(const EdgeIterator& original) :
			cell_iterator(original.cell_iterator){								DEB3("voro::EdgeIterator::__construct__(): copy");
		u = original.u;
		i = original.i;
		v = original.v;
		j = original.j;
	}

	EdgeDesc operator*(){														DEB3("voro::EdgeIterator::operator*()");
		return cell_iterator.current();
	}

	EdgeIterator& operator++(){													DEB3("voro::EdgeIterator::operator++()");
		cell_iterator.mark();
		cell_iterator.forward();
		if (!cell_iterator.isFinished() && cell_iterator.isMarked()){ //if we have run of edges on the current side...
			while (!cell_iterator.isFinished() && cell_iterator.isMarked()){ //jump until new unmarked edge is found (or end)
				cell_iterator.jump();
			}
			if(cell_iterator.isFinished()){ //if end return NULL
				cell_iterator = CellIterator();
				return EDGE_ITERATOR_NULL;//NULL
			}
		}
		//we are now on unmarked edge - remember it
		u = cell_iterator.current().v;
		i = cell_iterator.current().j;
		v = cell_iterator.inv().v;
		j = cell_iterator.inv().j;
		return *this;
	}

	bool operator==(const EdgeIterator& original) const{										DEB3("voro::EdgeIterator::operator==()");
		return (cell_iterator.cell() == 0 && original.cell_iterator.cell() == 0) || 		//test if null
				(cell_iterator.cell() == original.cell_iterator.cell() && 					//test if the same cell
						v == original.v && u == original.u );// && 								//test if at the same edge
						//cell_iterator.isMarked() == original.cell_iterator.isMarked());		//test if both edges has the same mark
																							// - important when iterating around a loop
																							//   (for example side)
	}

	bool operator!=(const EdgeIterator& original) const{										DEB3("voro::EdgeIterator::operator!=()");
		return !this->operator==(original);
	}

	inline EdgeIterator& go(int v, int j){														DEB3("voro::EdgeIterator::go(): v, j");
		this->cell_iterator.go(v, j);
		this->u = v;
		this->i = j;
		this->v = cell_iterator.inv().v;
		this->j = cell_iterator.inv().j;
		return *this;
	}
	inline EdgeIterator& go(const EdgeDesc& edge){												DEB3("voro::EdgeIterator::go(): edge");
		this->cell_iterator.go(edge.v, edge.j);
		this->u = edge.v;
		this->i = edge.j;
		this->v = cell_iterator.inv().v;
		this->j = cell_iterator.inv().j;
		return *this;
	}

	~EdgeIterator(){																			DEB3("voro::EdgeIterator::__destruct__()");
	}
};
EdgeIterator EdgeIterator::EDGE_ITERATOR_NULL = EdgeIterator();

class SideDesc{
public:
	CellIterator start, stop;

	SideDesc(){																					DEB3("voro::SideDesc::__construct__(): 0 args");
	}

	SideDesc(const CellIterator& start, const CellIterator& stop) :
			start(start), stop(stop){															DEB3("voro::SideDesc::__construct__(): 1 args");
	}

	SideDesc(const SideDesc& original) :
			start(original.start), stop(original.stop){											DEB3("voro::SideDesc::__construct__(): copy");
	}

	EdgeIterator begin(){																		DEB3("voro::SideDesc::begin()");
		return EdgeIterator(start);
	}
	EdgeIterator end(){																			DEB3("voro::SideDesc::end()");
		return EdgeIterator(stop);
	}

	double area(){																				DEB3("voro::SideDesc::area()");
		EdgeIterator it = this->begin();
		EdgeIterator end = this->end();
		double side_area = 0.0;
		double* p1 = start.vertexRaw(*it); DEB2("edge it = " << (*it).v << " " << (*it).j);
		double* p2 = start.vertexRaw(*(++it)); DEB2("edge it = " << (*it).v << " " << (*it).j);
		double* p3 = 0; ++it;
		for (; it != end; ++it){
			p3 = start.vertexRaw(*it); DEB2("edge it = " << (*it).v << " " << (*it).j);
			side_area += abs(dmga::algebra::dot3d(p1, p2, p3));
			p2 = p3;
		}
		return side_area * 0.5;
	}

	EdgeDesc baseEdge(){																		DEB3("voro::SideDesc::baseEdge()");
		return start.current();
	}
};


class SideIterator{
public:
	CellIterator current;
	CellIterator next;
	int side;

	static SideIterator SIDE_ITERATOR_NULL;

	SideIterator() :
			current(), next(), side(-1){ 														DEB3("voro::SideIterator::__construct__(): 0 args");
	}

	SideIterator(const CellIterator& cell_iterator) :
			current(cell_iterator),
			next(cell_iterator),
			side(0){																			DEB3("voro::SideIterator::__construct__(): 1 arg");
		//fast forward to the last edge of current side
		if (this->next.isFinished()){
			side = -1;
		}else{
			//rewind to the last edge on current side
			this->next.mark();
			this->next.forward();
			while (!this->next.isFinished() && !this->next.isMarked()){
				this->next.mark();
				this->next.forward();
			}
			if (!this->next.isFinished()){
				this->next.jump();
			}
		}
	}

	SideIterator(const SideIterator& original) :
			current(original.current),
			next(original.next),
			side(original.side){																DEB3("voro::SideIterator::__construct__(): copy");
	}

	SideDesc operator*(){																		DEB3("voro::SideIterator::operator*()");
		return SideDesc(current, next);
	}

	SideIterator& operator++(){																	DEB3("voro::SideIterator::operator++()");
		current = next;
		if (this->next.isFinished()){
			side = -1;
		}else{
			side++;
			//rewind to the last edge on current side
			this->next.mark();
			this->next.forward();
			while (!this->next.isFinished() && !this->next.isMarked()){
				this->next.mark();
				this->next.forward();
			}
			if (!this->next.isFinished()){
				this->next.jump();
			}
		}
		return *this;
	}

	bool operator==(const SideIterator& original) const{									DEB3("voro::SideIterator::operator==()");
		return (side == -1 && original.side == -1) ||
					(side != -1 && original.side != -1 &&
							current == original.current);
	}

	bool operator!=(const SideIterator& original) const{									DEB3("voro::SideIterator::operator!=()");
		return !this->operator==(original);
	}

	~SideIterator(){																		DEB3("voro::SideIterator::__destructor__()");
	}
};
SideIterator SideIterator::SIDE_ITERATOR_NULL = SideIterator();

}//namespace voro


#endif /* UTILS_H_ */
