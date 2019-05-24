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
typedef typename dmga::container::VoroppContainer<Geometry, Particle> Container;
typedef typename dmga::diagram::VoroppDiagram<Container, dmga::model::Cell> Diagram;
typedef typename dmga::diagram::VoroppDiagram<Container, dmga::model::CellEx> DiagramEx;
typedef typename dmga::shape::VoroppAlphaShape<Diagram> Shape;
typedef dmga::io::RepositoryTrajectory Trajectory;
typedef typename dmga::diagram::VoroppBasicKineticDiagram<Container, Trajectory, dmga::model::Cell> KineticDiagram;
typedef dmga::model::DodecahedronCellPreset DodecahedronPreset;

/**
 * compute on standard data from generator
 */
void compute(std::string repname){																			std::cout << "Computation started...\n" << std::flush;
	ifstream file_in(repname + "/" + repname + ".desc"); //open sample file
	double x, y, z, r;
	int id, molecule_id;
	std::string name, symbol;

	double x_size, y_size, z_size;
	bool x_periodic, y_periodic, z_periodic;

	file_in >> x_size >> y_size >> z_size >> x_periodic >> y_periodic >> z_periodic;
	Geometry geometry(x_size, y_size, z_size, x_periodic, y_periodic, z_periodic);		std::cout << "Geometry loaded...\n" << std::flush;
	Container atom_bank(geometry);
	while (file_in >> id >> name >> molecule_id >> symbol >> x >> y >> z >> r){
		atom_bank.add(Particle(id, x, y, z, r));
	}																					std::cout << "Positions loaded...\n" << std::flush;

	Trajectory trajectory(repname);														std::cout << "Trajectory loaded...\n" << std::flush;
	KineticDiagram kinetic_diagram(atom_bank, trajectory);								std::cout << "Diagram created...\n" << std::flush;

	auto kinetic_cells = kinetic_diagram.getCells(); 									std::cout << "Cells computed, start updating...\n" << std::flush;

	while (kinetic_diagram.advance()){													std::cout << "Next frame...\n" << std::flush;
		kinetic_diagram.updateCells();
		Shape shape(kinetic_diagram);
		auto it = kinetic_cells.begin();
		auto iend = kinetic_cells.end();
		for (; it < iend; ++it){
			//trajectory.ostream() << (*(*it));
			trajectory.ostream() << shape.getCellShape(*(*it));
		}
	}																					std::cout << "Computation finished...\n" << std::flush;
}

int main(int argc, char* argv[]) {
	std::string filename = "data";
	if (argc > 1){
		filename = std::string(argv[1]);
	}
	compute(filename);
	return 0;
}




