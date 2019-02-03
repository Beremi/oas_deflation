/**
 * This file contains command line utility for testing library
 */
#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <map>
#include "dmgalpha.h"
#include "io/io.hpp"

using namespace std;
using namespace dmga;


///////////////////////////////////////////////////////////////////////////////////
// Main testing ///////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////

typedef typename dmga::model::WeightedParticle<int> Particle;
typedef typename dmga::geometry::OrthogonalGeometry Geometry;
typedef typename dmga::container::VectorContainer<Geometry, Particle> Container;
typedef typename dmga::diagram::VoroppDiagram<Container, dmga::model::Cell> Diagram;
typedef typename dmga::diagram::VoroppDiagram<Container, dmga::model::CellEx> DiagramEx;
typedef typename dmga::shape::VoroppAlphaShape<Diagram> Shape;
typedef dmga::io::RepositoryTrajectory Trajectory;
typedef typename dmga::diagram::VoroppBasicKineticDiagram<Container, Trajectory, dmga::model::Cell> KineticDiagram;
typedef typename dmga::diagram::VoroppBasicKineticDiagram<Container, Trajectory, dmga::model::CellEx> KineticDiagramEx;
typedef dmga::model::DodecahedronCellPreset DodecahedronPreset;

//#define NO_OUTPUT 1;

void test_hardcoded_diagram(){															std::cout << "[TESTS][NEW TEST] test_hardcoded_diagram started...\n" << std::flush;
	ofstream file_out1("test_hardcoded_diagram/test_hardcoded_diagram00000000.voro");	std::cout << "[TESTS] Sample outupt opened...\n" << std::flush;
	ofstream file_out0("test_hardcoded_diagram/test_hardcoded_diagram00000000.coord");	std::cout << "[TESTS] Sample outupt opened...\n" << std::flush;
	ofstream file_out2("test_hardcoded_diagram/test_hardcoded_diagram.desc");			std::cout << "[TESTS] Sample outupt opened...\n" << std::flush;

	Geometry geometry(10.0, 10.0, 10.0, true, true, true);								std::cout << "[TESTS] Geometry loaded...\n" << std::flush;
	Container atom_bank(geometry);														std::cout << "[TESTS] Container ready...\n" << std::flush;
	file_out2 << "10.0 10.0 10.0 1 1 1\n";
	atom_bank.add(Particle(1, 5.0, 5.0, 5.0, 1.0));										std::cout << "[TESTS] Particle added...\n" << std::flush;
	file_out2 << "1 A 1 LIP 5.0 5.0 5.0 1.0\n";
	file_out0 << "5.0 5.0 5.0 1.0\n";
	atom_bank.add(Particle(2, 5.0-1.0, 5.0, 5.0, 1.0));									std::cout << "[TESTS] Particle added...\n" << std::flush;
	file_out2 << "2 A 1 LIP 4.0 5.0 5.0 1.0\n";
	file_out0 << "4.0 5.0 0.0 1.0\n";
	atom_bank.add(Particle(3, 5.0, 5.0-1.0, 5.0, 1.0));									std::cout << "[TESTS] Particle added...\n" << std::flush;
	file_out2 << "3 A 1 LIP 5.0 4.0 5.0 1.0\n";
	file_out0 << "5.0 4.0 0.0 1.0\n";
	atom_bank.add(Particle(4, 5.0, 5.0, 5.0-1.0, 1.0));									std::cout << "[TESTS] Particle added...\n" << std::flush;
	file_out2 << "4 A 1 LIP 5.0 5.0 4.0 1.0\n";
	file_out0 << "5.0 5.0 4.0 1.0\n";
	Diagram diagram(atom_bank);															std::cout << "[TESTS] Diagram created...\n" << std::flush;																			std::cout << "[TESTS] Presets added...\n" << std::flush;
	file_out2.close();
	file_out0.close();

	auto cells = diagram.getCells(); 													std::cout << "[TESTS] Cells computed, start updating...\n" << std::flush;
	auto it = cells.begin();															std::cout << "[TESTS] Iterator ready...\n" << std::flush;
	auto iend = cells.end();															std::cout << "[TESTS] Iterator end ready...\n" << std::flush;
	for (; it < iend; ++it){
		file_out1 << (*(*it));															std::cout << "[TESTS] Cell outputed...\n" << std::flush;
	}
	file_out1.close();																	std::cout << "[TESTS] Output closed...\n" << std::flush;

	std::cout << "[TESTS] test_hardcoded_diagram finished...\n" << std::flush;
}

/**
 * Example task - data from �ukasz
 */
void test_presets(){																	std::cout << "[TESTS][NEW TEST] test_presets started...\n" << std::flush;
	ifstream file_in("test_presets/test_presets.desc");									std::cout << "[TESTS] Sample input opened...\n" << std::flush;
	ofstream file_out("test_presets/test_presets00000000.voro");						std::cout << "[TESTS] Sample outupt opened...\n" << std::flush;
	double x, y, z, r;
	int id, molecule_id;

	std::string name, symbol;
	double x_size, y_size, z_size;
	bool x_periodic, y_periodic, z_periodic;

	file_in >> x_size >> y_size >> z_size >> x_periodic >> y_periodic >> z_periodic;	std::cout << "[TESTS] Container description read...\n" << std::flush;
	Geometry geometry(x_size, y_size, z_size, x_periodic, y_periodic, z_periodic);		std::cout << "[TESTS] Geometry loaded...\n" << std::flush;
	Container atom_bank(geometry);														std::cout << "[TESTS] Container ready...\n" << std::flush;
	std::vector<DodecahedronPreset> presets;											std::cout << "[TESTS] Storage for presets ready...\n" << std::flush;
	while (file_in >> id >> name >> molecule_id >> symbol >> x >> y >> z >> r){
		atom_bank.add(Particle(id, x, y, z, r));										std::cout << "[TESTS] Particle added...\n" << std::flush;
		presets.push_back(DodecahedronPreset(2.0));										std::cout << "[TESTS] Preset added...\n" << std::flush;
	}																					std::cout << "[TESTS] Positions loaded...\n" << std::flush;

	Diagram diagram(atom_bank);															std::cout << "[TESTS] Diagram created...\n" << std::flush;
	for (unsigned int i = 0; i < presets.size(); i++){
		diagram.addPreset(i, presets[i]);
	}																					std::cout << "[TESTS] Presets added...\n" << std::flush;

	auto cells = diagram.getCells(); 													std::cout << "[TESTS] Cells computed, start updating...\n" << std::flush;
	auto it = cells.begin();															std::cout << "[TESTS] Iterator ready...\n" << std::flush;
	auto iend = cells.end();															std::cout << "[TESTS] Iterator end ready...\n" << std::flush;
	for (; it < iend; ++it){
		file_out << (*(*it));															std::cout << "[TESTS] Cell outputed...\n" << std::flush;
	}
	file_out.close();																	std::cout << "[TESTS] Output closed...\n" << std::flush;
	file_in.close();																	std::cout << "[TESTS] Input closed...\n" << std::flush;

	std::cout << "[TESTS] test_presets finished...\n" << std::flush;
}

/**
 * Example task - data from �ukasz
 */
void test_diagram_cell(){																std::cout << "[TESTS][NEW TEST] test_diagram_cell started...\n" << std::flush;
	ifstream file_in("test_diagram_cell/test_diagram_cell.desc");						std::cout << "[TESTS] Sample input opened...\n" << std::flush;
	ofstream file_out("test_diagram_cell/test_diagram_cell00000000.voro");				std::cout << "[TESTS] Sample outupt opened...\n" << std::flush;
	double x, y, z, r;
	int id, molecule_id;

	std::string name, symbol;
	double x_size, y_size, z_size;
	bool x_periodic, y_periodic, z_periodic;

	file_in >> x_size >> y_size >> z_size >> x_periodic >> y_periodic >> z_periodic;	std::cout << "[TESTS] Container description read...\n" << std::flush;
	Geometry geometry(x_size, y_size, z_size, x_periodic, y_periodic, z_periodic);		std::cout << "[TESTS] Geometry loaded...\n" << std::flush;
	Container atom_bank(geometry);														std::cout << "[TESTS] Container ready...\n" << std::flush;
	while (file_in >> id >> name >> molecule_id >> symbol >> x >> y >> z >> r){
		atom_bank.add(Particle(id, x, y, z, r));										std::cout << "[TESTS] Particle added...\n" << std::flush;
	}																					std::cout << "[TESTS] Positions loaded...\n" << std::flush;
	Diagram diagram(atom_bank);															std::cout << "[TESTS] Diagram created...\n" << std::flush;																			std::cout << "[TESTS] Presets added...\n" << std::flush;

	auto cells = diagram.getCells(); 													std::cout << "[TESTS] Cells computed, start updating...\n" << std::flush;
	auto it = cells.begin();															std::cout << "[TESTS] Iterator ready...\n" << std::flush;
	auto iend = cells.end();															std::cout << "[TESTS] Iterator end ready...\n" << std::flush;
	for (; it < iend; ++it){
		file_out << (*(*it));															std::cout << "[TESTS] Cell outputed...\n" << std::flush;
	}
	file_out.close();																	std::cout << "[TESTS] Output closed...\n" << std::flush;
	file_in.close();																	std::cout << "[TESTS] Input closed...\n" << std::flush;

	std::cout << "[TESTS] test_diagram_cell finished...\n" << std::flush;
}

/**
 * Example task - data from �ukasz
 */
void test_diagram_cell_ex(){															std::cout << "[TESTS][NEW TEST] test_diagram_cell_ex started...\n" << std::flush;
	ifstream file_in("test_diagram_cell_ex/test_diagram_cell_ex.desc");					std::cout << "[TESTS] Sample input opened...\n" << std::flush;
	ofstream file_out("test_diagram_cell_ex/test_diagram_cell_ex00000000.voro");		std::cout << "[TESTS] Sample outupt opened...\n" << std::flush;
	double x, y, z, r;
	int id, molecule_id;

	std::string name, symbol;
	double x_size, y_size, z_size;
	bool x_periodic, y_periodic, z_periodic;

	file_in >> x_size >> y_size >> z_size >> x_periodic >> y_periodic >> z_periodic;	std::cout << "[TESTS] Container description read...\n" << std::flush;
	Geometry geometry(x_size, y_size, z_size, x_periodic, y_periodic, z_periodic);		std::cout << "[TESTS] Geometry loaded...\n" << std::flush;
	Container atom_bank(geometry);														std::cout << "[TESTS] Container ready...\n" << std::flush;
	while (file_in >> id >> name >> molecule_id >> symbol >> x >> y >> z >> r){
		atom_bank.add(Particle(id, x, y, z, r));										std::cout << "[TESTS] Particle added...\n" << std::flush;
	}																					std::cout << "[TESTS] Positions loaded...\n" << std::flush;
	DiagramEx diagram(atom_bank);															std::cout << "[TESTS] Diagram created...\n" << std::flush;																			std::cout << "[TESTS] Presets added...\n" << std::flush;

	auto cells = diagram.getCells(); 													std::cout << "[TESTS] Cells computed, start updating...\n" << std::flush;
	auto it = cells.begin();															std::cout << "[TESTS] Iterator ready...\n" << std::flush;
	auto iend = cells.end();															std::cout << "[TESTS] Iterator end ready...\n" << std::flush;
	for (; it < iend; ++it){
		file_out << (*(*it));															std::cout << "[TESTS] Cell outputed...\n" << std::flush;
	}
	file_out.close();																	std::cout << "[TESTS] Output closed...\n" << std::flush;
	file_in.close();																	std::cout << "[TESTS] Input closed...\n" << std::flush;

	std::cout << "[TESTS] test_diagram_cell_ex finished...\n" << std::flush;
}

/**
 * Example task - standard data from generator
 */
void test_kinetic_diagram_cell(){														std::cout << "[TESTS][NEW TEST] test_kinetic_diagram_cell started...\n" << std::flush;
	ifstream file_in("test_kinetic_diagram_cell/test_kinetic_diagram_cell.desc");		std::cout << "[TESTS] Sample input opened...\n" << std::flush;
	double x, y, z, r;
	int id, molecule_id;
	std::string name, symbol;

	double x_size, y_size, z_size;
	bool x_periodic, y_periodic, z_periodic;

	file_in >> x_size >> y_size >> z_size >> x_periodic >> y_periodic >> z_periodic;	std::cout << "[TESTS] Container description read...\n" << std::flush;
	Geometry geometry(x_size, y_size, z_size, x_periodic, y_periodic, z_periodic);		std::cout << "[TESTS] Geometry loaded...\n" << std::flush;
	Container atom_bank(geometry);														std::cout << "[TESTS] Container ready...\n" << std::flush;
	while (file_in >> id >> name >> molecule_id >> symbol >> x >> y >> z >> r){
		atom_bank.add(Particle(id, x, y, z, r));										std::cout << "[TESTS] Particle added...\n" << std::flush;
	}																					std::cout << "[TESTS] Positions loaded...\n" << std::flush;

	Trajectory trajectory("test_kinetic_diagram_cell");									std::cout << "[TESTS] Trajectory loaded...\n" << std::flush;
	KineticDiagram kinetic_diagram(atom_bank, trajectory);								std::cout << "[TESTS] Diagram created...\n" << std::flush;

	auto kinetic_cells = kinetic_diagram.getCells(); 									std::cout << "[TESTS] Cells computed, start updating...\n" << std::flush;

	while (kinetic_diagram.advance()){													std::cout << "[TESTS] Next frame...\n" << std::flush;
		kinetic_diagram.updateCells();													std::cout << "[TESTS] Cells updated...\n" << std::flush;
		auto it = kinetic_cells.begin();												std::cout << "[TESTS] Iterator created...\n" << std::flush;
		auto iend = kinetic_cells.end();												std::cout << "[TESTS] Iterator end created...\n" << std::flush;
		for (; it < iend; ++it){
			trajectory.ostream() << (*(*it));											std::cout << "[TESTS] Cell outputed...\n" << std::flush;
		}
	}
	std::cout << "[TESTS] test_kinetic_diagram_cell finished...\n" << std::flush;
}
/**
 * Example task - standard data from generator
 */
void test_kinetic_diagram_cell_ex(){														std::cout << "[TESTS][NEW TEST] test_kinetic_diagram_cell started...\n" << std::flush;
	ifstream file_in("test_kinetic_diagram_cell_ex/test_kinetic_diagram_cell_ex.desc");		std::cout << "[TESTS] Sample input opened...\n" << std::flush;
	double x, y, z, r;
	int id, molecule_id;
	std::string name, symbol;

	double x_size, y_size, z_size;
	bool x_periodic, y_periodic, z_periodic;

	file_in >> x_size >> y_size >> z_size >> x_periodic >> y_periodic >> z_periodic;	std::cout << "[TESTS] Container description read...\n" << std::flush;
	Geometry geometry(x_size, y_size, z_size, x_periodic, y_periodic, z_periodic);		std::cout << "[TESTS] Geometry loaded...\n" << std::flush;
	Container atom_bank(geometry);														std::cout << "[TESTS] Container ready...\n" << std::flush;
	while (file_in >> id >> name >> molecule_id >> symbol >> x >> y >> z >> r){
		atom_bank.add(Particle(id, x, y, z, r));										std::cout << "[TESTS] Particle added...\n" << std::flush;
	}																					std::cout << "[TESTS] Positions loaded...\n" << std::flush;

	Trajectory trajectory("test_kinetic_diagram_cell_ex");								std::cout << "[TESTS] Trajectory loaded...\n" << std::flush;
	KineticDiagramEx kinetic_diagram(atom_bank, trajectory);							std::cout << "[TESTS] Diagram created...\n" << std::flush;

	auto kinetic_cells = kinetic_diagram.getCells(); 									std::cout << "[TESTS] Cells computed, start updating...\n" << std::flush;

	while (kinetic_diagram.advance()){													std::cout << "[TESTS] Next frame...\n" << std::flush;
		kinetic_diagram.updateCells();													std::cout << "[TESTS] Cells updated...\n" << std::flush;
		auto it = kinetic_cells.begin();												std::cout << "[TESTS] Iterator created...\n" << std::flush;
		auto iend = kinetic_cells.end();												std::cout << "[TESTS] Iterator end created...\n" << std::flush;
		for (; it < iend; ++it){
			trajectory.ostream() << (*(*it));											std::cout << "[TESTS] Cell outputed...\n" << std::flush;
		}
	}
	std::cout << "[TESTS] test_kinetic_diagram_cell finished...\n" << std::flush;
}
/**
 * Example task - standard data from generator
 *
 * TODO: osobne dane do testów
 */
void test_kinetic_shape(){																std::cout << "[TESTS][NEW TEST] test_kinetic_shape started...\n" << std::flush;
	ifstream file_in("test_kinetic_diagram_cell/test_kinetic_diagram_cell.desc");		std::cout << "[TESTS] Sample input opened...\n" << std::flush;
	double x, y, z, r;
	int id, molecule_id;
	std::string name, symbol;

	double x_size, y_size, z_size;
	bool x_periodic, y_periodic, z_periodic;

	file_in >> x_size >> y_size >> z_size >> x_periodic >> y_periodic >> z_periodic;	std::cout << "[TESTS] Container description read...\n" << std::flush;
	Geometry geometry(x_size, y_size, z_size, x_periodic, y_periodic, z_periodic);		std::cout << "[TESTS] Geometry loaded...\n" << std::flush;
	Container atom_bank(geometry);														std::cout << "[TESTS] Container ready...\n" << std::flush;
	while (file_in >> id >> name >> molecule_id >> symbol >> x >> y >> z >> r){
		atom_bank.add(Particle(id, x, y, z, r));										std::cout << "[TESTS] Particle added...\n" << std::flush;
	}																					std::cout << "[TESTS] Positions loaded...\n" << std::flush;

	Trajectory trajectory("test_kinetic_diagram_cell");									std::cout << "[TESTS] Trajectory loaded...\n" << std::flush;
	KineticDiagram kinetic_diagram(atom_bank, trajectory);								std::cout << "[TESTS] Diagram created...\n" << std::flush;

	auto kinetic_cells = kinetic_diagram.getCells(); 									std::cout << "[TESTS] Cells computed, start updating...\n" << std::flush;

	while (kinetic_diagram.advance()){													std::cout << "[TESTS] Next frame...\n" << std::flush;
		kinetic_diagram.updateCells();													std::cout << "[TESTS] Cells updated...\n" << std::flush;
		Shape shape(kinetic_diagram);													std::cout << "[TESTS] Shape updated...\n" << std::flush;
		auto it = kinetic_cells.begin();												std::cout << "[TESTS] Iterator created...\n" << std::flush;
		auto iend = kinetic_cells.end();												std::cout << "[TESTS] Iterator end created...\n" << std::flush;
		for (; it < iend; ++it){
			trajectory.ostream() << *shape.getCellShape(*(*it));						std::cout << "[TESTS] Cell outputed...\n" << std::flush;
		}
	}
	std::cout << "[TESTS] test_kinetic_shape finished...\n" << std::flush;
}




int vectcmp(const std::vector<int>& a, const std::vector<int>& b){
	std::vector<int>::const_iterator a_it = a.begin();
	std::vector<int>::const_iterator a_end = a.end();
	std::vector<int>::const_iterator b_it = b.begin();
	std::vector<int>::const_iterator b_end = b.end();
	for ( ; a_it != a_end && a_it != a_end; a_it++, b_it++){
		if (*a_it < *b_it){
			return 1;  // a is smaller
		}else if (*a_it > *b_it){
			return -1; // b is smaller
		}
	}
	//one and/or the other ended, so we compare sizes
	if (a.size() < b.size()){ // a is smaller or equal
		return 1;
	}else if (a.size() > b.size()){
		return -1; //b is smaller
	}

	return 0; //equal
}
/**
 * Example task - standard data from generator
 */
void test_change_count_kinetic(){														std::cout << "test_change_count_kinetic started...\n" << std::flush;
	ifstream file_in("test_change_count_kinetic_3/test_change_count_kinetic_3.desc");		std::cout << "Sample input opened...\n" << std::flush;
	double x, y, z, r;
	int id, molecule_id;
	std::string name, symbol;

	double x_size, y_size, z_size;
	bool x_periodic, y_periodic, z_periodic;

	file_in >> x_size >> y_size >> z_size >> x_periodic >> y_periodic >> z_periodic;	std::cout << "Container description read...\n" << std::flush;
	Geometry geometry(x_size, y_size, z_size, x_periodic, y_periodic, z_periodic);		std::cout << "Geometry loaded...\n" << std::flush;
	Container atom_bank(geometry);														std::cout << "Container ready...\n" << std::flush;
	while (file_in >> id >> name >> molecule_id >> symbol >> x >> y >> z >> r){
		atom_bank.add(Particle(id, x, y, z, r));										std::cout << "Particle added...\n" << std::flush;
	}																					std::cout << "Positions loaded...\n" << std::flush;

	Trajectory trajectory("test_change_count_kinetic_3");								std::cout << "Trajectory loaded...\n" << std::flush;
	KineticDiagramEx kinetic_diagram(atom_bank, trajectory);							std::cout << "Diagram created...\n" << std::flush;

	auto kinetic_cells = kinetic_diagram.getCells(); 									std::cout << "Cells computed, start updating...\n" << std::flush;
	//compute statistic for comparision
	auto it = kinetic_cells.begin();													std::cout << "Iterator created...\n" << std::flush;
	auto iend = kinetic_cells.end();													std::cout << "Iterator end created...\n" << std::flush;
	int cell_num = 0;
	std::vector<std::vector<int> > current_list_of_neighbours(atom_bank.size());
	for (; it < iend; ++it, ++cell_num){
		//we have cell, we need to compute its characteristic...
		//first characteristic: neighbours - list all neighbours
		dmga::model::CellEx& cell = **it;
		auto sides_it = cell.sides.begin();
		auto sides_end = cell.sides.end();
		for (; sides_it != sides_end; ++sides_it){
			current_list_of_neighbours[cell_num].push_back(cell.getNeighbourId((*sides_it).baseEdge()));
		}
		sort(current_list_of_neighbours[cell_num].begin(), current_list_of_neighbours[cell_num].end());
		//second characteristic: list all sides, and sort neighbours on edges (TODO: think?)
	}
	int diffcount = 0;
	int samecount = 0;
	while (kinetic_diagram.advance()){													std::cout << "Next frame...\n" << std::flush;
		continue;
		kinetic_diagram.updateCells();													std::cout << "Cells updated...\n" << std::flush;
		auto it = kinetic_cells.begin();												std::cout << "Iterator created...\n" << std::flush;
		auto iend = kinetic_cells.end();												std::cout << "Iterator end created...\n" << std::flush;
		std::vector<std::vector<int> > last_list_of_neighbours = current_list_of_neighbours;
		cell_num = 0;
		for (; it < iend; ++it, ++cell_num){
			//std::cout << "comp cell " << cell_num << "... ";
			trajectory.ostream() << (*(*it));
			current_list_of_neighbours[cell_num].clear();
			dmga::model::CellEx& cell = **it;
			auto sides_it = cell.sides.begin();
			auto sides_end = cell.sides.end();
			for (; sides_it != sides_end; ++sides_it){
				current_list_of_neighbours[cell_num].push_back(cell.getNeighbourId((*sides_it).baseEdge()));
			}
			sort(current_list_of_neighbours[cell_num].begin(), current_list_of_neighbours[cell_num].end());
			int result = vectcmp(current_list_of_neighbours[cell_num], last_list_of_neighbours[cell_num]);
			if (result == 0){
				//std::cout << cell_num << " SAME AS PREVIOUS\n";
				samecount++;
			}else{
				//std::cout << cell_num << "NOT SAME AS PREVIOUS\n";
				diffcount++;
			}
		}
		std::cout << "Samecount = " << samecount << "\n";
		std::cout << "Diffcount = " << diffcount << "\n";
	}
	std::cout << "test_change_count_kinetic finished...\n" << std::flush;
}

int main(int argc, char* argv[]) {
	int all = 0;
	int ok = 0;

	//NOW IM TESTING SPEED
	//test_change_count_kinetic();
	//return 0;

	try{ all++; test_hardcoded_diagram(); ok++; }catch(...){ std::cout << "[TESTS] some error...";}
	try{ all++; test_diagram_cell(); ok++; }catch(...){ std::cout << "[TESTS] some error...";}
	try{ all++; test_diagram_cell_ex(); ok++; }catch(...){ std::cout << "[TESTS] some error...";}
	try{ all++; test_kinetic_diagram_cell(); ok++; }catch(...){ std::cout << "[TESTS] some error...";}
	try{ all++; test_kinetic_diagram_cell_ex(); ok++; }catch(...){ std::cout << "[TESTS] some error...";}
	try{ all++; test_kinetic_shape(); ok++; }catch(...){ std::cout << "[TESTS] some error...";}
	try{ all++; test_presets(); ok++; }catch(...){ std::cout << "[TESTS] some error...";}


	std::cout << "[TESTS] Summary: " << ok << "/" << all << "\n";
	if (all != ok) std::cout << "[TESTS] There were some errors... \n";
	return 0;
}




