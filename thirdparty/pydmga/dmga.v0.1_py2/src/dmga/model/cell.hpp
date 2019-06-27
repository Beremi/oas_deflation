/*
 * cell.hpp
 *
 *  Created on: 10-11-2012
 *      Author: Robson
 */

#ifndef CELL_HPP_
#define CELL_HPP_

#include <model/cell.h>

namespace dmga{

namespace model{

int Cell::size(){
	if (this->isEmpty()){
		return 0; //we have empty cell
	}else{
		return voropp_cell.number_of_faces();
	}
}

int Cell::getVerticesCount(){
	if (this->isEmpty()){
		return 0; //we have empty cell
	}else{
		return voropp_cell.p;
	}
}

int Cell::getSidesCount(){
	if (this->isEmpty()){
		return 0; //we have empty cell
	}else{
		return voropp_cell.number_of_faces();
	}
}

int Cell::getEdgesCountAt(int vertex_number){
	if (this->isEmpty()){
		return 0;
	}
	return voropp_cell.nu[vertex_number];
}

int Cell::getEdgesCount(){
	if (this->isEmpty()){
		return 0;
	}
	int count = 0;
	for (int i = 0; i < voropp_cell.p; i++){
		count += voropp_cell.nu[i];
	}
	return count;
}

// TODO: te ponizej tylko dlatego ze fajnie widac je w pythonie
// TODO: nalezy usunac i łapac wyjatki oraz przekazywac do pythona!
class A{};
class Aa{};
class Ab{};
class B{};
class C{};
class D{};
class E{};
class F{};

int Cell::getNeighbourId(int vertex_number, int edge_number){
	if (this->isEmpty()){
		throw A(); //TODO:
	}
	return voropp_cell.ne[vertex_number][edge_number];
}

int Cell::getNeighbourId(const voro::CellIterator& cellIterator){
	if (this->isEmpty()){
		throw Aa(); //TODO:
	}
	return voropp_cell.ne[cellIterator.current().v][cellIterator.current().j];
}

int Cell::getNeighbourId(const voro::EdgeDesc& edge){
	if (this->isEmpty()){
		throw Ab(); //TODO:
	}
	return voropp_cell.ne[edge.v][edge.j];
}

voro::CellIterator Cell::getCellIterator(){
	if (this->isEmpty()){
		throw B();
	}
	return voro::CellIterator(&this->voropp_cell, this->prefix_nu);
}

voro::CellIterator* Cell::newCellIterator(){
	if (this->isEmpty()){
		throw C();
	}
	return new voro::CellIterator(&this->voropp_cell, this->prefix_nu);
}

voro::CellIterator Cell::iterator(){
	if (!cell_iterator){
		throw D();
	}
	return *cell_iterator;
}

double* Cell::getVertexRaw(int i){
	if (this->isEmpty()){
		throw E();
	}
	return voropp_cell.pts + 3 * i;
}

Point3D Cell::getVertex(int i){
	if (this->isEmpty()){
		throw F();
	}
	return dmga::model::Point3D(getVertexRaw(i));
}

/**
 * outputs the cell
 */
std::ostream& operator<<(ostream& out, Cell& cell){
	if (cell.isEmpty()){
		out << "0\n";
		out << "0\n";
		return out;
	}
	auto& vcell = cell.voropp_cell;
	voro::CellIterator cell_it = voro::CellIterator(&vcell);
	double* point = vcell.pts;
	double* pend = vcell.pts + vcell.p * 3;
	out << vcell.p << "\n";
	for (; point < pend; point += 3){
		out << (*point) << " " << (point[1]) << " " << (point[2]) << "\n";
	}
	out << cell.size() << "\n";
	while (!cell_it.isFinished()){
		out << cell_it.current().v;
		cell_it.mark();
		cell_it.forward();

		while (!cell_it.isMarked()){
			out << " " << cell_it.current().v;
			cell_it.mark();
			cell_it.forward();
		}
		out << "\n";
		cell_it.jump();
	}
	return out;
}

}// namespace model

}// namespace dmga


#endif /* CELL_HPP_ */
