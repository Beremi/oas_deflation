/**
 * This file contains command line utility for testing library
 */

#define DEBUG 3

#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <map>
#include "dmgalpha.h"
#include "shape/sasa_shape.hpp"
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
typedef typename dmga::shape::VoroppSASAShape<Diagram> Shape;

//#define NO_OUTPUT 1;

int main(int argc, char* argv[]) {
	Geometry geometry(10.0, 10.0, 10.0, false, false, false);
	Container atom_bank(geometry);
//	atom_bank.add(Particle(0, 5.0, 5.0, 5.0, 1.0));
//	atom_bank.add(Particle(1, 5.0-1.0, 5.0, 5.0, 0.75));
//	atom_bank.add(Particle(2, 5.0+1.0, 5.0, 5.0, 0.75));
//	atom_bank.add(Particle(3, 5.0+0.5, 5.0+0.86, 5.0, 0.75));
//	atom_bank.add(Particle(4, 5.0+0.5, 5.0-0.86, 5.0, 0.75));
//	atom_bank.add(Particle(5, 5.0-0.5, 5.0+0.86, 5.0, 0.75));
//	atom_bank.add(Particle(6, 5.0-0.5, 5.0-0.86, 5.0, 0.75));
//	atom_bank.add(Particle(7, 5.0, 5.0, 5.0+1.5, 0.25));
//	atom_bank.add(Particle(8, 5.0, 5.0, 5.0-1.5, 0.9));
	atom_bank.add(Particle(0, 5.0, 5.0, 5.0, 1.0));
	atom_bank.add(Particle(1, 5.0-1.0, 5.0, 5.0, 1.0));
	atom_bank.add(Particle(2, 5.0+1.0, 5.0, 5.0, 1.0));
	atom_bank.add(Particle(1, 5.0, 5.0-1.0, 5.0, 1.0));
	atom_bank.add(Particle(2, 5.0, 5.0+1.0, 5.0, 1.0));
	Diagram diagram(atom_bank);
	DEB1("[SASA_TEST] Starting");
	Shape sasa_shape(diagram);																								DEB1("[SASA_TEST] sasa_shape initialized");
	dmga::shape::SASAContours* mesh = sasa_shape.computeOne(0);																	DEB1("[SASA_TEST] mesh computed");
	ofstream file_out("sasa_test_1.txt");																					DEB1("[SASA_TEST] file opened");
	//draw points
	file_out << "sasa = {\n";
	file_out << "\t'centers': " << "[]" << ",\n";			DEB1("[SASA_TEST] centers...");
	file_out << "\t'vertices': " << dmga::io::python::asList(mesh->vertices.begin(), mesh->vertices.end()) << ",\n";			DEB1("[SASA_TEST] vertices...");
	file_out << "\t'borders': [";
	auto border_it = mesh->borderIterator();
	bool coma = false;
	if (border_it.ready()) do {
		if (coma) file_out << ","; coma = true;
		file_out << "["
				 << border_it.currentPoint() << ","
				 << border_it.centerPointOnPlane() << ","
				 << border_it.nextPoint()
				 << "]";
	} while (border_it.iterate() != -1);																					DEB1("[SASA_TEST] borders...");
	file_out << "],\n";
	//draw all on polygons
	file_out << "\t'polygons': [";
	auto polygon_it = mesh->polygonIterator();
	coma = false;
	if (polygon_it.ready()) do {
		if (coma) file_out << ","; coma = true;
		file_out << "["
				 << polygon_it.currentPoint() << ","
				 << polygon_it.centerPointOnSphere() << ","
				 << polygon_it.nextPoint()
				 << "]";
	} while (polygon_it.iterate() != -1);
	file_out << "],\n";																										DEB1("[SASA_TEST] polygons DONE...");
	file_out << "\t'caps': [";																								DEB1("[SASA_TEST] output caps...");
	coma = false;
	for (auto cit = mesh->centers.begin(); cit != mesh->centers.end(); ++cit){
		if (cit->cap_radius_sq > 0.0){
			if (coma) file_out << ","; coma = true;
			double *pp = cit->point_on_plane.getRaw();
			file_out << "(" << pp[0] << ", " << pp[1] << "," << pp[2] << "," << sqrt(cit->cap_radius_sq) << ")";
		}
	}
	file_out << "],\n";																										DEB1("[SASA_TEST] caps DONE...");
	file_out << "\t'voronoi': ";																							DEB1("[SASA_TEST] output voronoi...");
	file_out << dmga::io::python::cellAsSides(diagram[0]);
	file_out << "\n";																										DEB1("[SASA_TEST] voronoi DONE...");
	file_out << "}\n";
	file_out.close();																										DEB1("[SASA_TEST] DONE");
	return 0;
}




