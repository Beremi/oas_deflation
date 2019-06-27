/*
 * kinetic_cell.h
 *
 *  Created on: 16-11-2012
 *      Author: Robson
 */

#ifndef KINETIC_CELL_H_
#define KINETIC_CELL_H_

#include <model/cell.hpp>

namespace dmga{

namespace model{

template<typename CellSpec, typename EventQueueSpec>
class KineticCell{
public:
	typedef CellSpec CellType;
	typedef EventQueueSpec EventQueueType;
	typedef typename EventQueueType::EventType EventType;

	CellType& start_cell;
	EventQueueType& events;

	CellType current_cell;


	KineticCell(CellType& cell, EventQueueType& events) : start_cell(cell), events(events){
		current_cell = start_cell;
	}

	CellType& getCell(double at_time){
		if (at_time < events.time){
			current_cell = start_cell;
			events.reset();
		}
		while (events.next() && events.topological.time <= at_time){
			events.topological.updateCell(current_cell);
		}
		events.update
		return current_cell;
	}
};

} //namespace model

} //namespace dmga

#endif /* KINETIC_CELL_H_ */
