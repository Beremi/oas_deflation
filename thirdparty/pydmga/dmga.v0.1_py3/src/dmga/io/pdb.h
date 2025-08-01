/*
 * pdb.h
 *
 *  Created on: 25-03-2013
 *      Author: Robson
 */

#ifndef PDB_H_
#define PDB_H_

#include <fstream>
#include <string>

namespace dmga{

namespace io{

class PDBDataParser{
public:
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
						const std::string& charge) = 0;

	virtual ~PDBDataParser(){
	}
};

class SimplePDBDataParser : public PDBDataParser{
protected:
	double m_x, m_y, m_z, m_r;
	int m_id;
public:
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
		m_id = serial;
		m_x = x;
		m_y = y;
		m_z = z;
		m_r = r;
		return true;
	}

	int id(){
		return m_id;
	}

	double x(){
		return m_x;
	}

	double y(){
		return m_y;
	}

	double z(){
		return m_z;
	}

	double r(){
		return m_r;
	}
};


template<typename InserterSpec>
int pdb_load_positions(std::istream& file_in, InserterSpec& ins){
	std::string cmd = "";
	double x, y, z, r, temp_factor;
	int serial, res_seq;
	char chain, alt_loc, i_code;
	std::string name, residue, seq_id, element, charge;
	double x_size, y_size, z_size;
	bool x_periodic, y_periodic, z_periodic;

	int added_count = 0;
	while(!file_in.eof()){
		char buff[256];
		file_in.getline(buff, 7); file_in.clear();
		cmd = buff;
		if (cmd == "CRYST1"){
			x_periodic = y_periodic = z_periodic = false;
			file_in.getline(buff, 10); file_in.clear(); { std::istringstream ss(buff); ss >> x_size; }//atom id
			file_in.getline(buff, 10); file_in.clear(); { std::istringstream ss(buff); ss >> y_size; }//atom id
			file_in.getline(buff, 10); file_in.clear(); { std::istringstream ss(buff); ss >> z_size; }//atom id
			file_in.getline(buff, 8); file_in.clear(); //[dump] (alpha angle)
			file_in.getline(buff, 8); file_in.clear(); //[dump] (beta angle)
			file_in.getline(buff, 8); file_in.clear(); //[dump] (gamma angle)
			file_in.getline(buff, 256); //[dump] everything!

			//ignore this data!

		} else if (cmd == "ATOM  "){
//			if (!atom_bank){
//				file_in.getline(buff, 256); //[dump] anything else!
//				continue;
//			}

			//I need to add +1 to each size because of the way the readline works - see documentation
			file_in.getline(buff, 6); file_in.clear(); { std::istringstream ss(buff); ss >> serial; }	//atom id (int)
			file_in.getline(buff, 2); file_in.clear(); 													//[dump] space
			file_in.getline(buff, 5); file_in.clear(); { name = buff; } 								//atom name (string)
			file_in.getline(buff, 2); file_in.clear(); { alt_loc = buff[0]; } 							//alternate location indicator (char)
			file_in.getline(buff, 4); file_in.clear(); { residue = buff; } 								//residue symbol (string)
			file_in.getline(buff, 2); file_in.clear(); 													//[dump] (space)
			file_in.getline(buff, 2); file_in.clear(); { chain = buff[0]; } 							//chain id (char)
			file_in.getline(buff, 5); file_in.clear(); { std::istringstream ss(buff); ss >> res_seq; }	//residue sequence no
			file_in.getline(buff, 5); file_in.clear(); { i_code = buff[0]; } 							//iCode + [dump] 3 spaces
			file_in.getline(buff, 9); file_in.clear(); { std::istringstream ss(buff); ss >> x; }		//x
			file_in.getline(buff, 9); file_in.clear(); { std::istringstream ss(buff); ss >> y; }		//y
			file_in.getline(buff, 9); file_in.clear(); { std::istringstream ss(buff); ss >> z; }		//z
			file_in.getline(buff, 7); file_in.clear(); { std::istringstream ss(buff); ss >> r; }		//r
			file_in.getline(buff, 7); file_in.clear(); {temp_factor = 0.0; /* TODO: */} 				//[dump] temperature factor
			file_in.getline(buff, 5); file_in.clear(); {seq_id = ""; /* TODO: */} 						//[dump]
			file_in.getline(buff, 3); file_in.clear(); {element = ""; /* TODO: */} 						//[dump]
			file_in.getline(buff, 3); file_in.clear(); {charge = ""; /* TODO: */} 						//[dump]
			file_in.getline(buff, 256); //[dump] anything else!

			bool test = ins.insert(serial,
									name,
									alt_loc,
									residue,
									chain,
									res_seq,
									i_code,
									x, y, z, r,
									temp_factor,
									seq_id,
									element,
									charge);
			if (test){
				added_count++;
			}
		} else {
			file_in.getline(buff, 255); //[dump] anything left!
		}
	}

	return added_count;
}

template<typename GeometrySpec, typename ContainerSpec, typename DataParserSpec>
int pdb_load(std::istream& file_in, DataParserSpec& data_parser, GeometrySpec* &geometry, ContainerSpec* &atom_bank){
	dmga::utils::safeDelete(atom_bank);
	dmga::utils::safeDelete(geometry);
	typedef typename ContainerSpec::ElementType ParticleType;

	std::string cmd = "";
	double x, y, z, r, temp_factor;
	int serial, res_seq;
	char chain, alt_loc, i_code;
	std::string name, residue, seq_id, element, charge;
	double x_size, y_size, z_size;
	bool x_periodic, y_periodic, z_periodic;

	int added_count = 0;
	while(!file_in.eof()){
		char buff[256];
		file_in.getline(buff, 7); file_in.clear();
		cmd = buff;
		if (cmd == "CRYST1"){
			x_periodic = y_periodic = z_periodic = false;
			file_in.getline(buff, 10); file_in.clear(); { std::istringstream ss(buff); ss >> x_size; }//atom id
			file_in.getline(buff, 10); file_in.clear(); { std::istringstream ss(buff); ss >> y_size; }//atom id
			file_in.getline(buff, 10); file_in.clear(); { std::istringstream ss(buff); ss >> z_size; }//atom id
			file_in.getline(buff, 8); file_in.clear(); //[dump] (alpha angle)
			file_in.getline(buff, 8); file_in.clear(); //[dump] (beta angle)
			file_in.getline(buff, 8); file_in.clear(); //[dump] (gamma angle)
			file_in.getline(buff, 256); //[dump] everything!

			//hardcoded- nie za bardzo rozumiem ta notacj� (Hermann-Mauguin space group symbol)
			x_periodic = true;
			y_periodic = true;
			z_periodic = true;

			geometry = new GeometrySpec(x_size, y_size, z_size, x_periodic, y_periodic, z_periodic);
			atom_bank = new ContainerSpec(*geometry);
		} else if (cmd == "ATOM  "){
			if (!atom_bank){
				file_in.getline(buff, 256); //[dump] anything else!
				continue;
			}

			//I need to add +1 to each size because of the way the readline works - see documentation
			file_in.getline(buff, 6); file_in.clear(); { std::istringstream ss(buff); ss >> serial; }	//atom id (int)
			file_in.getline(buff, 2); file_in.clear(); 													//[dump] space
			file_in.getline(buff, 5); file_in.clear(); { name = buff; } 								//atom name (string)
			file_in.getline(buff, 2); file_in.clear(); { alt_loc = buff[0]; } 							//alternate location indicator (char)
			file_in.getline(buff, 4); file_in.clear(); { residue = buff; } 								//residue symbol (string)
			file_in.getline(buff, 2); file_in.clear(); 													//[dump] (space)
			file_in.getline(buff, 2); file_in.clear(); { chain = buff[0]; } 							//chain id (char)
			file_in.getline(buff, 5); file_in.clear(); { std::istringstream ss(buff); ss >> res_seq; }	//residue sequence no
			file_in.getline(buff, 5); file_in.clear(); { i_code = buff[0]; } 							//iCode + [dump] 3 spaces
			file_in.getline(buff, 9); file_in.clear(); { std::istringstream ss(buff); ss >> x; }		//x
			file_in.getline(buff, 9); file_in.clear(); { std::istringstream ss(buff); ss >> y; }		//y
			file_in.getline(buff, 9); file_in.clear(); { std::istringstream ss(buff); ss >> z; }		//z
			file_in.getline(buff, 7); file_in.clear(); { std::istringstream ss(buff); ss >> r; }		//r
			file_in.getline(buff, 7); file_in.clear(); {temp_factor = 0.0; /* TODO: */} 				//[dump] temperature factor
			file_in.getline(buff, 5); file_in.clear(); {seq_id = ""; /* TODO: */} 						//[dump]
			file_in.getline(buff, 3); file_in.clear(); {element = ""; /* TODO: */} 						//[dump]
			file_in.getline(buff, 3); file_in.clear(); {charge = ""; /* TODO: */} 						//[dump]
			file_in.getline(buff, 256); //[dump] anything else!

			bool test = data_parser.accept(serial,
											name,
											alt_loc,
											residue,
											chain,
											res_seq,
											i_code,
											x, y, z, r,
											temp_factor,
											seq_id,
											element,
											charge);
			if (test){
				atom_bank->add(ParticleType(data_parser.id(), data_parser.x(), data_parser.y(), data_parser.z(), data_parser.r()));
				added_count++;
			}
		} else {
			file_in.getline(buff, 255); //[dump] anything left!
		}
	}

	return added_count;
}

} //namespace io

} //namespace dmga


#endif /* PDB_H_ */
