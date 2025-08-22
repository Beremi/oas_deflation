/*
 * input.h
 *
 *  Created on: 10-08-2012
 *      Author: Robson
 */

#ifndef INPUT_H_
#define INPUT_H_

#include <iostream>
#include <fstream>
#include <string>
#include <trajectory/trajectory.hpp>
#include <exceptions/exceptions.hpp>
#include <model/model.hpp>
#include <map>

namespace dmga{

namespace io{

std::string trim(std::string s, char symbol = ' '){
	unsigned int i = 0, j = 0;
	while (s[i] == symbol) i++;
	while (i + j < s.length() && s[i + j] != symbol) j++;
	return s.substr(i, j);
}

std::string paddedNumber(unsigned int number, unsigned int width, char fill = '0'){
	std::string strNumber = "";
	int n = number;
	while (n > 0){
		strNumber = (char)((n % 10)+'0') + strNumber;
		n /= 10;
	}
	while (strNumber.length() < width){
		strNumber = fill + strNumber;
	}
	return strNumber;
}

std::string paddedDouble(double number, unsigned int left, unsigned int right, char fill = '0'){
	std::string strNumber = "";
	std::string strRight = "";
	bool sign = number < 0; if (number < 0 ) number = -number;
	int n = (int)number;
	double mm = (number - (double)n); mm = (mm < 0 ? -mm : mm);
	for(unsigned int i = 0; i < right; i++){
		mm *= 10;
	}
	int m = (int)mm;
	if (n == 0){
		strNumber == "0";
	}
	while (n > 0){
		strNumber = (char)((n % 10)+'0') + strNumber;
		n /= 10;
	}
	if (sign){
		strNumber = "-" + strNumber;
	}
	while (strNumber.length() < left){
		strNumber = fill + strNumber;
	}
	while (m > 0){
		strRight = (char)((m % 10)+'0') + strRight;
		m /= 10;
	}
	while (strRight.length() < right){
		strRight = strRight + "0";
	}
	strNumber = strNumber + '.' + strRight;
	return strNumber;
}

/**
 * this kinetic input is more for debug purposes as it maintains
 * structure that is easy to generate, browse, read and understand
 * (see python tools voro_render_gui and voro_gen_gui)
 *
 * in future we will probably have a KineticInput that is suited for
 * existing biological file formats such as .PDB with multiple frames
 */
class RepositoryTrajectory : public dmga::trajectory::Trajectory{
private:
	int first_frame_number;
	int current_frame_number;
	std::string repository_name;
	std::string separator;
	std::ifstream input_file;
	std::ofstream voro_file;
	int frame_size;
	double* raw_data;

	std::string createBaseFilename(){
		return repository_name + separator + repository_name;
	}

	std::string createCoordFilename(){
		return createBaseFilename() + paddedNumber(current_frame_number, 8) + ".coord";
	}

	std::string createDescFilename(){
		return createBaseFilename() + ".desc";
	}

	std::string createVoroFilename(){
		return createBaseFilename() + paddedNumber(current_frame_number, 8) + ".voro";
	}

	bool reloadCurrent(){
		std::ifstream new_file( createCoordFilename().c_str() );
		if (!new_file){
			return false;
		}
		input_file.close();
		voro_file.close();
		input_file.open(  createCoordFilename().c_str() );
		voro_file.open( createVoroFilename().c_str() );



		return true;
	}
public:
	RepositoryTrajectory(std::string repository_name):
			first_frame_number(0),
			current_frame_number(-1), //before simulation
			repository_name(repository_name),
			separator("/"),
			frame_size(0),
			raw_data(0){
		input_file.open( (createBaseFilename() + ".desc").c_str());
		if (!input_file){
			throw exceptions::InvalidInput(std::string("RepositoryTrajectory::constructor(): cannot open description file in ") + repository_name);
		}
		std::string data;
		if (std::getline(input_file, data)){
			while (std::getline(input_file, data)){
				++frame_size;
			}
		}else{
			throw exceptions::InvalidInput(std::string("RepositoryTrajectory::constructor(): invalid description file format in ") + repository_name);
		}
		raw_data = new double[3 * frame_size];
		voro_file.open( (createVoroFilename()).c_str() );
	}

	~RepositoryTrajectory(){
		input_file.close();
		voro_file.close();
	}

	std::istream& istream(){
		return input_file;
	}

	std::ostream& ostream(){
		return voro_file;
	}

	bool advance(){
		this->current_frame_number++;
		if (!reloadCurrent()){
			this->current_frame_number--;
			return false;
		}
		double x, y, z, r;

		double* data_ptr = raw_data;
		while (input_file >> x >> y >> z >> r){
			*(data_ptr++) = x;
			*(data_ptr++) = y;
			*(data_ptr++) = z;
		}

		return true;
	}

	double* getRaw(){
		return raw_data;
	}

	double* getRaw(int i){
		return raw_data + 3 * i;
	}

	dmga::model::Point3D getPos(int i){
		return dmga::model::Point3D(raw_data + 3 * i);
	}
};

template<typename CellSpec>
class PymolCell{
public:
	typedef CellSpec CellType;
	CellType& cell;
	std::string color;
	std::string group;
	static std::map<std::string, int> base_points_by_group;

	PymolCell(CellType& item, std::string color = "br0", std::string group = ""): cell(item), color(color), group(group){
	}

	friend std::ostream& operator<<(std::ostream& out, const PymolCell& item){
		auto& vcell = item.cell.voropp_cell;
		int id = item.cell.number;
		voro::CellIterator cell_it = voro::CellIterator(&vcell);
		double* point = vcell.pts;
		double* pend = vcell.pts + vcell.p * 3;

		std::string pymol_id;
		int base_point_id;
		if(item.group == ""){
			pymol_id = "cell_" + paddedNumber(id, 6);
			base_point_id = 0;
		}else{
			pymol_id = item.group;
			if (PymolCell<CellType>::base_points_by_group.find(item.group) == PymolCell<CellType>::base_points_by_group.end()){
				PymolCell<CellType>::base_points_by_group.insert(std::make_pair(item.group, 0));
			}
			base_point_id = PymolCell<CellType>::base_points_by_group[item.group];
			PymolCell<CellType>::base_points_by_group[item.group] = base_point_id + vcell.p;
		}

		auto sid = paddedNumber(id, 6);
		for (; point < pend; point += 3){
			out << "pseudoatom " << pymol_id << ",pos=[" << (*point) << "," << (point[1]) << "," << (point[2]) << "],color=" << item.color << "\n";
		}
		int last_v;
		int current_v;
		while (!cell_it.isFinished()){
			last_v = base_point_id + cell_it.current().v + 1;
			cell_it.mark();
			cell_it.forward();

			while (!cell_it.isMarked()){
				current_v = base_point_id + cell_it.current().v + 1;
				out << "bond " << pymol_id << "/PSDO/P/PSD`1/PS" << last_v <<"," << pymol_id << "/PSDO/P/PSD`1/PS" << current_v << "\n";
				cell_it.mark();
				cell_it.forward();
				last_v = current_v;
			}
			//add last edge!
			current_v = base_point_id + cell_it.current().v + 1;
			out << "bond " << pymol_id << "/PSDO/P/PSD`1/PS" << last_v <<"," << pymol_id << "/PSDO/P/PSD`1/PS" << current_v << "\n";
			cell_it.jump();
		}
		return out;
	}
};

template<typename CellSpec>
std::map<std::string, int> PymolCell<CellSpec>::base_points_by_group;

namespace pdb{

std::string hetatm(int id, double x, double y, double z){
	std::ostringstream out;
	out << "HETATM" << dmga::io::paddedNumber(id, 5, ' ') << " PS1  PSD P   1    " << dmga::io::paddedDouble(x , 3, 4, ' ') << dmga::io::paddedDouble(y , 3, 4, ' ') << dmga::io::paddedDouble(z, 3, 4, ' ') << "  0.00  0.00      PSDOPS  \n";
	return out.str();
}

std::string conect(int id1, int id2){
	std::ostringstream out;
	out << "CONECT" << dmga::io::paddedNumber(id1, 5, ' ') << dmga::io::paddedNumber(id2, 5, ' ') << "\n";
	return out.str();
}

}//namespace pdb

namespace pml{

}//namespace pml


template<typename CellSpec>
class PDBCell{
public:
	typedef CellSpec CellType;
	CellType& cell;
	std::string color;
	std::string group;
	static std::map<std::string, int> base_points_by_group;

	PDBCell(CellType& item, std::string color = "br0", std::string group = ""): cell(item), color(color), group(group){
	}

	friend std::ostream& operator<<(std::ostream& out, const PDBCell& item){
		auto& vcell = item.cell.voropp_cell;
		int id = item.cell.number;
		voro::CellIterator cell_it = voro::CellIterator(&vcell);
		double* point = vcell.pts;
		double* pend = vcell.pts + vcell.p * 3;

		std::string pymol_id;
		int base_point_id;
		if(item.group == ""){
			pymol_id = "cell_" + paddedNumber(id, 6);
			base_point_id = 0;
		}else{
			pymol_id = item.group;
			if (PDBCell<CellType>::base_points_by_group.find(item.group) == PDBCell<CellType>::base_points_by_group.end()){
				PDBCell<CellType>::base_points_by_group.insert(std::make_pair(item.group, 0));
			}
			base_point_id = PDBCell<CellType>::base_points_by_group[item.group];
			PDBCell<CellType>::base_points_by_group[item.group] = base_point_id + vcell.p;
		}

		auto sid = paddedNumber(id, 6);
		int i = 0;
		for (; point < pend; point += 3, ++i){
			out << "HETATM" << dmga::io::paddedNumber(base_point_id + i + 1, 5, ' ') << " PS1  PSD P   1    " << dmga::io::paddedDouble(point[0] , 4, 3, ' ') << dmga::io::paddedDouble(point[1] , 3, 4, ' ') << dmga::io::paddedDouble(point[2] , 3, 4, ' ') << "  0.00  0.00      PSDOPS  \n";
		}
		int last_v;
		int current_v;
		while (!cell_it.isFinished()){
			last_v = base_point_id + cell_it.current().v + 1;
			cell_it.mark();
			cell_it.forward();

			while (!cell_it.isMarked()){
				current_v = base_point_id + cell_it.current().v + 1;
				out << "CONECT" << dmga::io::paddedNumber(current_v, 5, ' ') << dmga::io::paddedNumber(last_v, 5, ' ') << "\n";
				cell_it.mark();
				cell_it.forward();
				last_v = current_v;
			}
			//add last edge!
			current_v = base_point_id + cell_it.current().v + 1;
			out << "CONECT" << dmga::io::paddedNumber(current_v, 5, ' ') << dmga::io::paddedNumber(last_v, 5, ' ') << "\n";
			cell_it.jump();
		}
		return out;
	}
};

template<typename CellSpec>
std::map<std::string, int> PDBCell<CellSpec>::base_points_by_group;

namespace python{

	template<typename IteratorSpec>
	std::string asList(IteratorSpec start, IteratorSpec stop){
		std::ostringstream out;
		bool coma = false;
		out << "[";
		for (; start != stop; ++start){
			if (coma) out << ", "; coma = true;
			out << *start;
		}
		out << "]";
		return out.str();
	}

	template<typename CellSpec>
	std::string cellAsSides(CellSpec& item){
		std::ostringstream out;
		auto& vcell = item.voropp_cell;
		voro::CellIterator cell_it = voro::CellIterator(&vcell);

		out << "[";
		bool coma_list = false;
		while (!cell_it.isFinished()){
			if (coma_list) out << ", ";
			//we are at the new side
			bool coma = false;
			out << "[";
			while (!cell_it.isMarked()){
				double* p = cell_it.vertex().getRaw();
				if (coma) out << ", ";
				out << "(" << p[0] << "," << p[1] << "," << p[2] << ")";
				cell_it.mark();
				cell_it.forward();
				coma = true;
			}
			out << "]";
			cell_it.jump();
			coma_list = true;
		}
		out << "]";
		return out.str();
	}

}

} //namespace io

} //namespace dmga


#endif /* INPUT_H_ */
