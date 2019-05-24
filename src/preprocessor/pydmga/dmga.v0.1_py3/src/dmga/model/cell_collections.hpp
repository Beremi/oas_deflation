/*
 * cell_collections.hpp
 *
 *  Created on: 21-03-2013
 *      Author: Robson
 */

#ifndef CELL_COLLECTIONS_HPP_
#define CELL_COLLECTIONS_HPP_

#include "cell_collections.h"

namespace dmga{

namespace model{

/* VertexCollection inner class */
void VertexCollection::update(){																			DEB3("VertexCollection::update()");
}
/* VertexCollection inner class */
VertexCollection::VertexCollection(Cell& parent_cell) : parent_cell(parent_cell){
}
/* VertexCollection inner class */
voro::VertexIterator VertexCollection::begin(){
	return voro::VertexIterator(parent_cell.voropp_cell.pts);
}
/* VertexCollection inner class */
voro::VertexIterator VertexCollection::end(){
	return voro::VertexIterator(parent_cell.voropp_cell.pts + parent_cell.voropp_cell.p * 3);
}
/* VertexCollection inner class */
int VertexCollection::size() const{
	return parent_cell.getVerticesCount();
}
/* VertexCollection inner class */
int VertexCollection::id(int v, int j) const {
	return v;
}

/* EdgesCollection inner class */
void EdgesCollection::update(){																				DEB3("EdgesCollection::update()");
	utils::safeDelete(prefix_nu);
	int n = parent_cell.voropp_cell.p;
	prefix_nu = new int[n + 1];
	prefix_nu[0] = 0;
	for (int i = 1; i < n + 1; i++){
		prefix_nu[i] = prefix_nu[i-1] + parent_cell.voropp_cell.nu[i-1];
	}
	edges_count = prefix_nu[n];
}

/* EdgesCollection inner class */
EdgesCollection::EdgesCollection(Cell& parent_cell) : parent_cell(parent_cell), edges_count(0), prefix_nu(0){
}
/* EdgesCollection inner class */
voro::EdgeIterator EdgesCollection::begin(){
	return voro::EdgeIterator(parent_cell.getCellIterator());
}
/* EdgesCollection inner class */
voro::EdgeIterator EdgesCollection::end(){
	return voro::EdgeIterator::EDGE_ITERATOR_NULL;
}
/* EdgesCollection inner class */
int EdgesCollection::size() const{
	return edges_count;
}
/* EdgesCollection inner class */
int EdgesCollection::id(int v, int j) const{
	return voro::edges::voroToEdge(prefix_nu, v, j);
}
/* EdgesCollection inner class */
EdgesCollection::~EdgesCollection(){
	utils::safeDelete(prefix_nu);
}

/* SidesCollection inner class */
void SidesCollection::update(){																					DEB3("SidesCollection::update()");
//TODO: zrobic to inaczej (bez edges)
	if (!parent_cell.prefix_nu){
		parent_cell.update(); //make sure we have prefix_nu
		DEB3("prefix_nu updated!");
	}
	voro::CellIterator it = parent_cell.getCellIterator();
	int* prefix_nu = parent_cell.prefix_nu;
	int vertex_count = parent_cell.getVerticesCount();
	DEB3("vertex_count: " << vertex_count);
	int size = prefix_nu[vertex_count];
	if (size > 200){
		for(int k = 0; k < vertex_count; k++){
			DEB("prefix_nu[" << k << "] = " << prefix_nu[k] );
		}
	}
	DEB3("size: " << size);
	side_by_edge.resize(size);
	int side = 0;
	while (!it.isFinished()){
		int v = it.current().v;
		int j = it.current().j;
		side_by_edge[voro::edges::voroToEdge(prefix_nu, v, j)] = side;
		it.mark();
		it.forward();
		while(!it.isMarked()){
			v = it.current().v;
			j = it.current().j;
			side_by_edge[voro::edges::voroToEdge(prefix_nu, v, j)] = side;
			it.mark();
			it.forward();
		}
		//go to the next side!
		side++;
		it.jump();
	}
	DEB3("SidesCollection::update() UPDATED:, size =  " << size << ", side = " << side);
}
/* SidesCollection inner class */
SidesCollection::SidesCollection(Cell& parent_cell) : parent_cell(parent_cell){
}
/* SidesCollection inner class */
voro::SideIterator SidesCollection::begin(){
	return voro::SideIterator(parent_cell.getCellIterator());
}
/* SidesCollection inner class */
voro::SideIterator SidesCollection::end(){
	return voro::SideIterator::SIDE_ITERATOR_NULL;
}
/* SidesCollection inner class */
int SidesCollection::size() const{
	return parent_cell.getSidesCount();
}
/* SidesCollection inner class */
int SidesCollection::id(int v, int j) const{
	return side_by_edge[voro::edges::voroToEdge(parent_cell.prefix_nu, v, j)];
}

}//namespace model

}//namespace dmga


#endif /* CELL_COLLECTIONS_HPP_ */
