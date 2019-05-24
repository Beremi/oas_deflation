/*
 * cell.h
 *
 *  Created on: 10-11-2012
 *      Author: Robson
 */

#ifndef CELL_H_
#define CELL_H_

#include <base.h>
#include <3rd/voro/voro++.hh>

#include <vector>
#include <model/primitives.hpp>
#include <utils/utils.hpp>
#include <io/io.h>

namespace dmga{

namespace model{

class BaseCellPreset : public DMGAObject, public voro::wall{
public:
	std::string show(){ return "dmga::model::BaseCellPreset"; }

	voro::voronoicell_base* vcell;

	BaseCellPreset(voro::voronoicell_base* vcell = 0) : vcell(vcell){
	}

	void setTarget(voro::voronoicell_base* vcell){
		this->vcell = vcell;
	}

	virtual bool point_inside(double x,double y,double z) = 0;
	virtual bool cut_cell(voro::voronoicell &c,double x,double y,double z) = 0;
	virtual bool cut_cell(voro::voronoicell_neighbor &c,double x,double y,double z) = 0;

	virtual ~BaseCellPreset(){	DEB3("BaseCellPreset::__destruct__():");
	}
};

class DodecahedronCellPreset : public BaseCellPreset {
public:
	std::string show(){ return "dmga::model::DodecahedronCellPreset"; }

	typedef BaseCellPreset super;

	static double Phi;
	static voro::voronoicell base;
	static voro::voronoicell_neighbor base_neighbor;

	static voro::voronoicell makeBase(){
		DEB2("Phi = " << Phi);
		voro::voronoicell v;
		v.init(-2,2,-2,2,-2,2);
		v.plane(0,Phi,1);v.plane(0,-Phi,1);v.plane(0,Phi,-1);
		v.plane(0,-Phi,-1);v.plane(1,0,Phi);v.plane(-1,0,Phi);
		v.plane(1,0,-Phi);v.plane(-1,0,-Phi);v.plane(Phi,1,0);
		v.plane(-Phi,1,0);v.plane(Phi,-1,0);v.plane(-Phi,-1,0);
		return v;
	}

	static voro::voronoicell_neighbor makeBaseNeighbor(){
		DEB2("Phi = " << Phi);
		voro::voronoicell_neighbor v;
		v.init(-2,2,-2,2,-2,2);
		v.plane(0,Phi,1);v.plane(0,-Phi,1);v.plane(0,Phi,-1);
		v.plane(0,-Phi,-1);v.plane(1,0,Phi);v.plane(-1,0,Phi);
		v.plane(1,0,-Phi);v.plane(-1,0,-Phi);v.plane(Phi,1,0);
		v.plane(-Phi,1,0);v.plane(Phi,-1,0);v.plane(-Phi,-1,0);
		return v;
	}

	bool point_inside(double x,double y,double z) {
		return false;
	}

	bool cut_cell(voro::voronoicell &c,double x,double y,double z) {
		if (this->vcell == 0 || (voro::voronoicell_base*)&c == this->vcell){
			c = base;
			double* p = c.pts;
			double* pend = c.pts + 3 * c.p;
			for (; p < pend; ++p){
				*p *= size;
			}
		}
		return true;
	}

	bool cut_cell(voro::voronoicell_neighbor &c,double x,double y,double z) {
		if (this->vcell == 0 || (voro::voronoicell_base*)&c == this->vcell){
			c = base_neighbor;
			double* p = c.pts;
			double* pend = c.pts + 3 * c.p;
			for (; p < pend; ++p){
				*p *= size;
			}
		}
	    return true;
	}

	double size;

	DodecahedronCellPreset(double size = 1.0, voro::voronoicell_base* vcell = 0) : super(vcell), size(size){
	}

	~DodecahedronCellPreset(){ 			DEB3("DodecahedronCellPreset::__destruct__():");
		//destructor...
	}
};

double DodecahedronCellPreset::Phi = 0.5*(1+sqrt(5.0));
voro::voronoicell DodecahedronCellPreset::base(DodecahedronCellPreset::makeBase());
voro::voronoicell_neighbor DodecahedronCellPreset::base_neighbor(DodecahedronCellPreset::makeBaseNeighbor());


class Cell : public DMGAObject{
public:
	std::string show(){ return "dmga::model::Cell"; }

	int number;
	voro::voronoicell_neighbor voropp_cell;
	voro::CellIterator* cell_iterator;
	int* prefix_nu;
	bool empty;
	/**
	 * constructor
	 * creates an default cell for the particle with given number
	 */
	Cell(int number) : number(number), cell_iterator(0), prefix_nu(0), empty(true){
	}
	/**
	 * constructor
	 * creates an default cell for the particle with given number and voro++ cell
	 */
	Cell(int number, const voro::voronoicell_neighbor& cell) : number(number), voropp_cell(cell), cell_iterator(0), prefix_nu(0), empty(false){
		update();
	}
	/**
	 * constructor
	 * creates an default cell for the particle with given number
	 */
	Cell(const Cell& original) : number(original.number), voropp_cell(original.voropp_cell), cell_iterator(0), prefix_nu(0), empty(original.empty){
		update();
	}

	inline bool isEmpty(){
		return empty;
	}
	inline void setEmpty(bool state){
		empty = state;
	}

	double volume(){
		if (this->isEmpty()){
			return 0.0; //we have empty cell
		}else{
			return this->voropp_cell.volume() * 8.0; //TODO: make it understand if the cell is normalized or not!
		}
	}

	double area(){
		if (this->isEmpty()){
			return 0.0; //we have empty cell
		}else{
			return this->voropp_cell.surface_area() * 4.0; //TODO: make it understand if the cell is normalized or not!
		}
	}
	/**
	 * returns number of sides of this cell
	 */
	inline int size();
	/**
	 * returns the number of vertices of this cell
	 */
	inline int getVerticesCount();
	/**
	 * returns the number of faces in this cell
	 */
	inline int getSidesCount();
	/**
	 * returns the number of faces in this cell
	 */
	inline int getEdgesCount();
	/**
	 * retuns count of the edges at vertex i (in general position should be 3)
	 * vertex_number \in \{0 .. verticesCount()-1\}
	 */
	inline int getEdgesCountAt(int vertex_number);
	/**
	 * returns i-th neighbour of this cell (i \in \{0 .. sidesCount()-1\})
	 * vertex_number \in \{0 .. verticesCount()-1\}
	 * edge_number \in \{0 .. edgesCountAt(vertex_number)-1\}
	 */
	inline int getNeighbourId(int vertex_number, int edge_number);
	/**
	 * returns neighbour to this cell on the side on which the cellIterator actually is
	 */
	inline int getNeighbourId(const voro::CellIterator& cellIterator);
	/**
	 * returns neighbour to this cell on the side on which the given edge actually is
	 */
	inline int getNeighbourId(const voro::EdgeDesc& edge);
	/**
	 * returns most detailed iterator over this cell
	 */
	voro::CellIterator getCellIterator();
	/**
	 * returns most detailed iterator over this cell by creating it with new operator
	 * it is handy when you understands how pointer works - this may be
	 * faster than passing return value by reference. But you are responsible to safely
	 * release it by delete.
	 */
	voro::CellIterator* newCellIterator();
	/**
	 * if necessary creates internal cell_iterator and returns it
	 * this is mostly for internal use inside library, but you can
	 * use this to optimize some code if you understands how stuff works
	 * if you are not sure, you should use getCellIterator() instead
	 */
	voro::CellIterator iterator();
	/**
	 * returns i-th vertex as raw array of doubles (consecutive entries are x, y, z)
	 */
	inline double* getVertexRaw(int i);
	/**
	 * returns i-th vertex as a Point3D class
	 */
	inline dmga::model::Point3D getVertex(int i);
	/**
	 * updates internal data, that helps working with the cell
	 * should be called after recomputing Voro++ cell.
	 */
	void update(){
		utils::safeDeleteArray(prefix_nu);
		utils::safeDelete(cell_iterator);
		this->cell_iterator = new voro::CellIterator(&this->voropp_cell);
		this->prefix_nu = this->cell_iterator->makeShared();
	}
	/**
	 * outputs the cell
	 */
	friend std::ostream& operator<<(ostream& out, Cell& cell);
	/**
	 * destructor
	 */
	virtual ~Cell(){			DEB3("Cell::__destruct__():");
		utils::safeDeleteArray(prefix_nu);
		utils::safeDelete(cell_iterator);
	}
};

}//namespace model;

}//namespace dmga;


#endif /* CELL_H_ */
