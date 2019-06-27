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
#include "io/pdb.hpp"

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
typedef dmga::model::DodecahedronCellPreset DodecahedronPreset;

#define NO_OUTPUT 1;

class DPPCParser : public dmga::io::SimplePDBDataParser{
	bool heavy_atoms_only;
	std::string weighting;
	int current_molecule;
	int last_part;
public:
	typedef dmga::io::SimplePDBDataParser super;

	vector<int> molecules;
	vector<std::string> molecules_desc;
	std::vector<dmga::model::Point3D> orig_coords;
	std::vector<bool> heavy_atom;
	double radius_lookup[256];

	DPPCParser(bool heavy_atoms_only = true, std::string weighting = "none") :
			heavy_atoms_only(heavy_atoms_only),
			weighting(weighting),
			current_molecule(0),
			last_part(3){
		for (int i = 0; i < 256; i++) radius_lookup[i] = 1.0;
		for (int i = 0; i < 256; i++) radius_lookup[i] = 1.0;
		//Van Der Waals radii:
		if (weighting == "vanderwaals"){
			std::cout << "Using Van Der Waals radii weighting\n" << std::flush;
			radius_lookup['H'] = 1.09;
			radius_lookup['C'] = 1.70;
			radius_lookup['N'] = 1.55;
			radius_lookup['O'] = 1.52;
			radius_lookup['F'] = 1.47;
			radius_lookup['P'] = 1.80;
			radius_lookup['S'] = 1.80;
		}

		//Van Der Waals radii^2:
		if (weighting == "vanderwaals2"){
			std::cout << "Using Van Der Waals radii squares weighting\n" << std::flush;
			radius_lookup['H'] = 1.09 * 1.09;
			radius_lookup['C'] = 1.70 * 1.70;
			radius_lookup['N'] = 1.55 * 1.55;
			radius_lookup['O'] = 1.52 * 1.52;
			radius_lookup['F'] = 1.47 * 1.47;
			radius_lookup['P'] = 1.80 * 1.80;
			radius_lookup['S'] = 1.80 * 1.80;
		}

		//Atomic Mass:
		if (weighting == "mass"){
			std::cout << "Atomic mass weighting\n" << std::flush;
			radius_lookup['H'] = 1.0;
			radius_lookup['C'] = 12.0;
			radius_lookup['N'] = 14.0;
			radius_lookup['O'] = 16.0;
			radius_lookup['F'] = 19.0;
			radius_lookup['P'] = 31.0;
			radius_lookup['S'] = 32.0;
		}

		//Number of electron orbits:
		if (weighting == "orbit"){
			std::cout << "Number of electron orbits weighting\n" << std::flush;
			radius_lookup['H'] = 1.0;
			radius_lookup['C'] = 2.0;
			radius_lookup['N'] = 2.0;
			radius_lookup['O'] = 2.0;
			radius_lookup['F'] = 2.0;
			radius_lookup['P'] = 3.0;
			radius_lookup['S'] = 3.0;
		}
	}

	virtual bool accept(int serial,
			const std::string& name,
			char alt_loc,
			const std::string& residue,
			char chain,
			int res_seq,
			char i_code,
			double x, double y, double z, double r,
			double temp_factor,
			const std::string& seq_id,
			const std::string& element,
			const std::string& charge){
		super::accept(serial, name, alt_loc, residue, chain, res_seq, i_code, x, y, z, r, temp_factor, seq_id, element, charge);
		if ( name[1] != 'H' || !heavy_atoms_only){
			heavy_atom.push_back(name[1] != 'H');
			orig_coords.push_back(dmga::model::Point3D(x, y, z));
			this->m_r = radius_lookup[(unsigned int)name[1]];

			if (residue == "PCH" || residue == "PAL"){
				if (last_part == 3 && res_seq != 3){
					++current_molecule;
				}
				last_part = res_seq;
				molecules.push_back(current_molecule);
				molecules_desc.push_back(io::trim(name) + " " + io::trim(residue) + " " + io::paddedNumber(res_seq, 1));
			}
			return true;
		}
		return false;
	}
};

void computation_graph_test(DiagramEx& diagram){
	std::cout << "computation_graph_test started\n" << std::flush;
	dmga::graph::ComputationGraph<DiagramEx, dmga::graph::FIFOQueue, dmga::graph::AllExpander> cg(diagram);
	cg.append(1);
	auto& g = diagram.container.geometry;
	double test_volume = g.xSize() * g.ySize() * g.zSize();
	double computed_volume = 0;
	while (auto* cell = cg.next()){
		//std::cout << "Computing cell: " << cell->number << "\n" << std::flush;
		computed_volume += cell->volume();
		dmga::utils::safeDelete(cell);
	}
	std::cout << "VOLUME TEST: " << test_volume << " vs " << computed_volume << "\n" << std::flush;
}

class IntervalOnlyExpander : public dmga::graph::AcceptExpander{
public:
	int start;
	int end;

	IntervalOnlyExpander(int start, int end) : start(start), end(end){
	}

	bool accept(int cell_number){
		return start <= cell_number && cell_number <= end;
	}
};

//BAD BAD BAD GLOBAL VARIABLE BUT I NEED TO DO THIS FAST!
std::vector<bool> to_cast;

double water_groups_test(DiagramEx& diagram, int waterStartsAt, bool with_output = false, std::string output_file_base = ""){
	std::cout << "water_groups_test started\n" << std::flush;
	IntervalOnlyExpander expander(waterStartsAt, diagram.container.size());
	dmga::graph::ComputationGraph<DiagramEx, dmga::graph::FIFOQueue, IntervalOnlyExpander> cg(diagram, expander);
	int groups = 0;
	int out_file_num = 0;
	double water_volume = 0.0;
	std::ofstream file_out;
	for (int i = waterStartsAt; i < diagram.container.size(); i++){
		cg.append(i);
		if (cg.hasNext()){
			groups++;
			if (with_output){
				file_out.open((output_file_base + "_surface_" + dmga::io::paddedNumber(groups, 5) + " " + dmga::io::paddedNumber(out_file_num, 4) + ".pdb").c_str());
			}
			double this_group_volume = 0.0;
			double this_group_surface_area = 0.0;
			int base_vertex_id = 1;
			try{
				while (auto* cell = cg.next()){
					if (cell->isEmpty()){
						std::cout << ">>>>> Empty Water Cell\n" << std::flush;
						dmga::utils::safeDelete(cell);
						continue;
					}
					this_group_volume += cell->volume();
					std::vector<bool> printed_vertex(cell->getVerticesCount(), false);
					if (base_vertex_id + cell->getVerticesCount() >= 100000){
						file_out.close();
						out_file_num++;
						file_out.open((output_file_base + "_surface_" + dmga::io::paddedNumber(groups, 5) + "_" + dmga::io::paddedNumber(out_file_num, 4) + ".pdb").c_str());
						base_vertex_id = 1;
					}
					auto sit = cell->sides.begin();
					auto send = cell->sides.end();
					for (; sit != send; ++sit){
						auto side = (*sit);
						int neighbour = cell->getNeighbourId(side.baseEdge().v, side.baseEdge().j);
						if (neighbour < waterStartsAt){

							//$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
							while (to_cast.size() <= (size_t)neighbour) to_cast.push_back(false);
							to_cast[neighbour] = true;
							//$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$

						//if (neighbour < 160 * 50){
							std::cout << "Side " << sit.side << " is on surface with neighbour = " << neighbour << "\n" << std::flush;
							this_group_surface_area += side.area();

							auto eit = side.begin();
							auto eend = side.end();
							int first_v = (*eit).v;
							if (!printed_vertex[first_v]){
								double* pos = cell->getVertexRaw(first_v);
								if (with_output){
									file_out << dmga::io::pdb::hetatm(base_vertex_id + first_v, pos[0], pos[1], pos[2]);
								}
								printed_vertex[first_v] = true;
							}
							int last_v = first_v;
							++eit;
							for (; eit != eend; ++eit){
								int v = (*eit).v;
								if (!printed_vertex[v]){
									double* pos = cell->getVertexRaw(v);
									if (with_output){
										file_out << dmga::io::pdb::hetatm(base_vertex_id + v, pos[0], pos[1], pos[2]);
									}
									printed_vertex[v] = true;
								}
								if (with_output){
									file_out << dmga::io::pdb::conect(base_vertex_id + last_v, base_vertex_id + v);
								}
								last_v = v;
							}
							file_out << dmga::io::pdb::conect(base_vertex_id + last_v, base_vertex_id + first_v);
						}else{
							std::cout << "Side " << sit.side << " is not on surface\n" << std::flush;
						}
					}
					base_vertex_id += cell->getVerticesCount();
					dmga::utils::safeDelete(cell);
				}
			}catch(dmga::exceptions::BaseException& ex){
				std::cout << "Uhu... \n" << ex << std::flush;
				throw ex;
			}catch(int a){
				std::cout << "INT (" << a << ") exception in WATER TEST\n" << std::flush;
				throw a;
			}catch(...){
				std::cout << "Some exception in WATER TEST\n" << std::flush;
				throw;
			}
			water_volume += this_group_volume;
			std::cout << "new group found with volume = " << this_group_volume << ", surface area = " << this_group_surface_area << "\n" << std::flush;
			if (with_output){
				file_out.close();
			}
		}
	}
	std::cout << "THERE ARE " << groups << " groups of water in this file...\n" << std::flush;
	return water_volume;
}

/**
 * Example task from dr. Murzyn
 */
void compute_dppc(std::string filename, DPPCParser& parser){							std::cout << "Task 01 started...\n" << std::flush;
	ifstream file_in(filename + ".pdb"); //open sample file
	ofstream file_out(filename + ".dist"); //open sample file
	#ifndef NO_OUTPUT
	ofstream file_out_pml(filename + ".pml"); //open sample file
	ofstream file_out_pdb(filename + "_0000.pdb"); //open sample file
	ofstream file_out_surface_pml(filename + "_surf.pml"); //open sample file
	ofstream file_out_surface_pdb(filename + "_surf.pdb"); //open sample file
	#endif

	Geometry* geometry = 0;
	Container* atom_bank = 0;
	int test_size = dmga::io::pdb_load(file_in, parser, geometry, atom_bank);	std::cout << "Positions loaded...\n" << std::flush;
	double test_volume = geometry->xSize() * geometry->ySize() * geometry->zSize();

	DiagramEx diagram(*atom_bank, Diagram::CACHE_OFF);										std::cout << "Diagram created...\n" << std::flush;

	std::cout << "TEST SIZE: " << atom_bank->size() << " " << test_size << "\n";
	try{
		double all_volume = 0.0;
		double water_volume = water_groups_test(diagram, parser.molecules.size(), true, filename);
		std::cout << "Water Volume = " << water_volume << "\n";

		int last_molecule = -1;
		#ifndef NO_OUTPUT
		int pml_file = 0;
		#endif
		int current_atom_in_molecule = 0;
		file_out << "molecule" << " " << "atom" << " " << "volume" << " " << "is_heavy" << " " << "name" << " " << "symbol" << " " << "part" "\n";
		unsigned int i = 0;
		for (; i < parser.molecules.size(); ++i){
			if (parser.molecules[i] != last_molecule){
				current_atom_in_molecule = 0;
				#ifndef NO_OUTPUT
				file_out_pml.close();
				file_out_pdb.close();
				file_out_pml.open( (filename + "_" + dmga::io::paddedNumber(pml_file, 4) + ".pml").c_str() ); //open sample file
				file_out_pdb.open( (filename + "_" + dmga::io::paddedNumber(pml_file, 4) + ".pdb").c_str() ); //open sample file
				pml_file++;
				#endif
			}
			last_molecule = parser.molecules[i];
			current_atom_in_molecule++;
			auto* cell = diagram.getCell(i);
			double volume = cell->volume();
			if (cell->isEmpty()) std::cout << ">>>>>>> Empty cell!\n";
			all_volume += volume;
			file_out << parser.molecules[i] << " " << current_atom_in_molecule << " " << volume << " " << (parser.heavy_atom[i] ? "1" : "0") << " " << parser.molecules_desc[i] << "\n";
			#ifndef NO_OUTPUT
			if (!cell->isEmpty()){
				voro::denormalizeCell(cell->voropp_cell, atom_bank->getRawCoords(i));
				voro::normalizeCell(cell->voropp_cell, parser.orig_coords[i].x, parser.orig_coords[i].y, parser.orig_coords[i].z);
				file_out_pml << dmga::io::PymolCell<dmga::model::CellEx>(*cell, "br" + dmga::io::paddedNumber(parser.molecules[i] % 10, 1));
				file_out_pdb << dmga::io::PDBCell<dmga::model::CellEx>(*cell, "br" + dmga::io::paddedNumber(parser.molecules[i] % 10, 1), "dppc" + dmga::io::paddedNumber(parser.molecules[i], 4));
			}
			#endif
			delete cell;
		}
		std::cout << "Lipid volume = " << all_volume << "\n";

		//all_volume += water_groups_test(diagram, i);
		std::cout << "VOLUME TEST := " << water_volume + all_volume << " vs " << test_volume << "\n" << std::flush;

	}catch(...){
		std::cout << "Some exception...\n" << std::flush;

		file_out.close();
		file_in.close();
		#ifndef NO_OUTPUT
		file_out_pml.close();
		file_out_pdb.close();
		#endif
	}
	file_out.close();
	file_in.close();
	#ifndef NO_OUTPUT
	file_out_pml.close();
	file_out_pdb.close();
	file_out_surface_pml.close();
	file_out_surface_pdb.close();
	#endif
	delete atom_bank;																		std::cout << "Task 01 finished...\n" << std::flush;
}

//$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
//$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$


class DPPC2DCastParser : public dmga::io::SimplePDBDataParser{
public:
	typedef dmga::io::SimplePDBDataParser super;

	DPPC2DCastParser(){
	}

	virtual bool accept(int serial,
			const std::string& name,
			char alt_loc,
			const std::string& residue,
			char chain,
			int res_seq,
			char i_code,
			double x, double y, double z, double r,
			double temp_factor,
			const std::string& seq_id,
			const std::string& element,
			const std::string& charge){
		super::accept(serial, name, alt_loc, residue, chain, res_seq, i_code, x, y, z, r, temp_factor, seq_id, element, charge);
		if ( (name[1] != 'H') && (to_cast[serial-1]) && (serial-1 < 160 * 130)){
			this->m_z = 0.0; //cast!
			return true;
		}
		return false;
	}
};

void cast_test(std::string filename){
	ifstream file_in(filename + ".pdb"); //open sample file
//	ofstream file_out(filename + ".cast"); //open sample file
	ofstream file_out("_cast_" + filename + ".pml"); //open sample file
	// int pml_file = 0;

	DPPC2DCastParser parser;
	Geometry* insert_geometry = 0;
	Container* insert_atom_bank = 0;
	int test_size = dmga::io::pdb_load(file_in, parser, insert_geometry, insert_atom_bank);	std::cout << test_size << " positions loaded...\n" << std::flush;

	Geometry* geometry = new Geometry(insert_geometry->xSize(), insert_geometry->ySize(), insert_geometry->zSize(), true, true, false);
	Container* atom_bank = new Container(*geometry);
	for (int i = 0; i < insert_atom_bank->size(); i++){
		atom_bank->add(insert_atom_bank->get(i));
	}

	DiagramEx diagram(*atom_bank, Diagram::CACHE_OFF);
	int pymol_id = 1;
	double testArea = geometry->xSize() * geometry->ySize();
	double area = 0.0;
	for (int i = 0; i < atom_bank->size(); i++){
		auto* cell = diagram.getCell(i);
		auto s_it = cell->sides.begin();
		auto s_end = cell->sides.end();
		for (; s_it != s_end; ++s_it){
			auto s_desc = *s_it;
			int neigh = cell->getNeighbourId(s_desc.baseEdge().v, s_desc.baseEdge().j);
			if (neigh < 0){
				auto it = s_desc.begin();
				auto end = s_desc.end();
				area += s_desc.area();
				std::cout << i << " += " << s_desc.area() << "\n";

				int start_id = pymol_id;
				double* p = cell->getVertexRaw(it.u);
				file_out << "pseudoatom " << "dupa" << ",pos=[" << (*p) << "," << (p[1]) << "," << (p[2]) << "]" << "\n";
				pymol_id++;
				++it;
				for (; it != end; ++it){
					double* p = cell->getVertexRaw(it.u);
					file_out << "pseudoatom " << "dupa" << ",pos=[" << (*p) << "," << (p[1]) << "," << (p[2]) << "]" << "\n";
					file_out << "bond " << "dupa" << "/PSDO/P/PSD`1/PS" << (pymol_id-1) <<"," << "dupa" << "/PSDO/P/PSD`1/PS" << pymol_id << "\n";
					pymol_id++;
				}
				file_out << "bond " << "dupa" << "/PSDO/P/PSD`1/PS" << (pymol_id-1) <<"," << "dupa" << "/PSDO/P/PSD`1/PS" << start_id << "\n";
			}
		}
	}
	std::cout << "testArea = " << testArea << " vs " << area << " vs " << area / 2.0 << "\n";
}

//$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
//$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$

int main(int argc, char* argv[]) {
	std::string filename = "data";
	if (argc > 1){
		filename = std::string(argv[1]);
	}
	std::string have_atoms_options = "all";
	std::string weighting_option = "none";
	if (argc > 2){
		have_atoms_options = std::string(argv[2]);
	}
	if (argc > 3){
		weighting_option = std::string(argv[3]);
		weighting_option = io::trim(weighting_option);
	}
	std::cout << "working on " << filename << ".pdb... \n" << std::flush;
	std::cout << "weighting flag: " << weighting_option << "\n";
	std::cout << "heavy atoms only flag: " << have_atoms_options << "\n";
	DPPCParser parser((have_atoms_options == "hao"), weighting_option);
	compute_dppc(filename, parser); /* weighting option \in none|orbit|mass|vanderwaals|vanderwaals2 */
	//cast_test(filename);
	return 0;
}




