/*
 * cell_collections.h
 *
 *  Created on: 21-03-2013
 *      Author: Robson
 */

#ifndef CELL_COLLECTIONS_H_
#define CELL_COLLECTIONS_H_

#include <3rd/voro/voro++.hh>

#include <vector>
#include <model/primitives.hpp>
#include <model/cell.hpp>
#include <model/cell_components.hpp>
#include <utils/utils.hpp>
#include <io/io.h>

#include "cell_components.hpp"

namespace dmga{

namespace model{

class VertexCollection{
public:
	/** parent cell */
	Cell& parent_cell;
	/** updates internal information. Use after parent cell recomputing */
	void update();
	/** constructor */
	VertexCollection(Cell& parent_cell);
	/** return the first element of points */
	voro::VertexIterator begin();
	/** return the element right after all points */
	voro::VertexIterator end();
	/** returns the number of points */
	int size() const;
	/** returns unique id of the vertex at the vertex v and edge j */
	inline int id(int v, int j) const;
};

class EdgesCollection{
public:
	/** parent cell */
	Cell& parent_cell;
	/** total edges count (forward and backward edges count as different edges) */
	int edges_count;
	/** this allows computing of the unique id of the edge */
	int* prefix_nu;
	/** restores prefix_nu table after recomputing of the cell */
	void update();
	/** constructor */
	EdgesCollection(Cell& parent_cell);
	/** first edge in collection */
	voro::EdgeIterator begin();
	/** returns iterator that says 'NULL' - it indicates the end of collection */
	voro::EdgeIterator end();
	/** returns number of edges (forward and backward edges count as different edges) */
	inline int size() const;
	/** unique id of the j-th edge at vertex v */
	inline int id(int v, int j) const;
	/** destructor */
	~EdgesCollection();
};

class SidesCollection{
public:
	/** parent cell */
	Cell& parent_cell;
	/** holds information edge_id -> side_id */
	std::vector<int> side_by_edge;
	/** updates the side_by_edge table. Use after recomputing of parent cell */
	void update();
	/** constructor */
	SidesCollection(Cell& parent_cell);
	/** returns first element of sides */
	voro::SideIterator begin();
	/** returns element after all sides */
	voro::SideIterator end();
	/** returns number of sides */
	int size() const;
	/** returns unique id of the side clockwise to the given edge */
	int id(int v, int j) const;
};

}//namespace model

}//namespace dmga


#endif /* CELL_COLLECTIONS_H_ */
