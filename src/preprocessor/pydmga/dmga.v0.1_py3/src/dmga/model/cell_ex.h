/*
 * cell_ex.h
 *
 *  Created on: 21-03-2013
 *      Author: Robson
 */

#ifndef CELL_EX_H_
#define CELL_EX_H_

#include "cell.hpp"
#include "cell_components.hpp"
#include "cell_collections.hpp"

namespace dmga{

namespace model{

class CellEx : public Cell{
public:
	typedef Cell super;

	VertexCollection vertices;
	EdgesCollection edges;
	SidesCollection sides;

	CellEx(Cell& base) : super(base), vertices(*this), edges(*this), sides(*this){
		update();
	}
	/**
	 * constructor
	 * creates an default cell for the particle with given number
	 */
	CellEx(int number) : super(number), vertices(*this), edges(*this), sides(*this){
	}
	/**
	 * constructor
	 * creates an default cell for the particle with given number and voro++ cell
	 */
	CellEx(int number, const voro::voronoicell_neighbor& cell) : super(number, cell), vertices(*this), edges(*this), sides(*this){
		update();
	}
	/**
	 * constructor
	 * creates an default cell for the particle with given number
	 */
	CellEx(const CellEx& original) : super(original.number, original.voropp_cell), vertices(*this), edges(*this), sides(*this){
		update();
	}

	void update(){
		super::update();
		vertices.update();
		edges.update();
		sides.update();
	}

};

} //namespace model

} //namespace dmga;

#endif /* CELL_EX_H_ */
