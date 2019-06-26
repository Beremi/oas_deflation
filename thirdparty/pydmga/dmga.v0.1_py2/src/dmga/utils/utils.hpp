/*
 * utils.hpp
 *
 *  Created on: 10-11-2012
 *      Author: Robson
 */

#ifndef UTILS_HPP_
#define UTILS_HPP_

#include <utils/utils.h>

namespace voro{

void normalizeCell(voro::voronoicell_neighbor& pc, double x, double y, double z){
	double* point = pc.pts;
	double* pend = pc.pts + 3 * pc.p;
	for (; point < pend; point += 3){
		point[0] = 0.5 * point[0] + x;
		point[1] = 0.5 * point[1] + y;
		point[2] = 0.5 * point[2] + z;
	}
}

void denormalizeCell(voro::voronoicell_neighbor& pc, double x, double y, double z){
	double* point = pc.pts;
	double* pend = pc.pts + 3 * pc.p;
	for (; point < pend; point += 3){
		point[0] = (point[0] - x) * 2.0;
		point[1] = (point[1] - y) * 2.0;
		point[2] = (point[2] - z) * 2.0;
	}
}

void normalizeCell(voro::voronoicell_neighbor& pc, double* raw){
	normalizeCell(pc, raw[0], raw[1], raw[2]);
}

void denormalizeCell(voro::voronoicell_neighbor& pc, double* raw){
	denormalizeCell(pc, raw[0], raw[1], raw[2]);
}

namespace edges{
inline void next(voro::voronoicell_base& vcell, int v, int j, int& u, int& k){
	u = vcell.ed[v][j];
	k = vcell.cycle_up(vcell.ed[v][j + vcell.nu[v]], u); //6.12.2013 ... bylo na koncu v); zamiast u);
}

inline void prev(voro::voronoicell_base& vcell, int v, int j, int& z, int& i){
	int rev = vcell.cycle_down(j, v);
	z = vcell.ed[v][rev];
	i = vcell.ed[v][vcell.nu[v] + rev];
}

inline void inv(voro::voronoicell_base& vcell, int v, int j, int& u, int& l){
	u = vcell.ed[v][j];
	l = vcell.ed[v][vcell.nu[v] + j];
}

inline void cw(voro::voronoicell_base& vcell, int v, int j, int& w, int& m){
	w = v, m = vcell.cycle_down(j, v);
}

inline void ccw(voro::voronoicell_base& vcell, int v, int j, int& w, int& m){
	w = v, m = vcell.cycle_up(j, v);
}

int voroToEdge(const int* prefix_nu, int vertex, int jthedge) {
	return prefix_nu[vertex] + jthedge;
}

} //namespace edges

CellIterator::CellIterator(): vcell(0), prefix_nu(0), shared_prefix_nu(0), v(0), j(0){
	DEB3("CellIterator::__construct__(): 0 args");
}

CellIterator::CellIterator(voro::voronoicell_base* pcell, int* precomputed_prefix_nu, int v, int j):
		vcell(pcell),
		prefix_nu(0), //make prefix array for each vertex, it will hold offset in visited array for edges for that vertex
		visited(pcell->number_of_edges() << 1), // *2 because we need for each edge two half-edges (forward and backward)
		v(v),
		j(j){
	if (precomputed_prefix_nu){																		DEB3("CellIterator::__construct__(): base*, int*");
		prefix_nu = precomputed_prefix_nu;
		shared_prefix_nu = precomputed_prefix_nu;
	}else{
		prefix_nu = new int[pcell->p + 1];
		prefix_nu[0] = 0;
		for (int i = 1; i < vcell->p + 1; i++){
			prefix_nu[i] = prefix_nu[i-1] + vcell->nu[i-1];
		}
		shared_prefix_nu = 0;
	}
}

CellIterator::CellIterator(const CellIterator& original) :
		vcell(original.vcell),
		prefix_nu(0), //make prefix array for each vertex, it will hold offset in visited array for edges for that vertex
		//visited(vcell ? (vcell->number_of_edges() << 1) : 1), // *2 because we need for each edge two half-edges (forward and backward)
		visited(original.visited),
		v(original.v),
		j(original.j){																				DEB3("CellIterator::__construct__(): copy...");
	if (original.vcell == 0){
		shared_prefix_nu = 0;
	}else{
		shared_prefix_nu = original.shared_prefix_nu;
		if (original.shared_prefix_nu != original.prefix_nu){
			prefix_nu = new int[vcell->p + 1];
			prefix_nu[0] = 0;
			for (int i = 1; i < vcell->p + 1; i++){
				prefix_nu[i] = prefix_nu[i-1] + vcell->nu[i-1];
			}
		}else{
			prefix_nu = shared_prefix_nu;
		}
	}
}

CellIterator& CellIterator::operator=(const CellIterator& original){								DEB3("CellIterator::operator=(): assignment!");
	vcell = original.vcell;
	visited = original.visited;
	v = original.v;
	j = original.j;
	if (shared_prefix_nu != prefix_nu){
		dmga::utils::safeDeleteArray(prefix_nu);
	}
	if (original.vcell == 0){
		shared_prefix_nu = 0;
		prefix_nu = 0;
	}else{
		shared_prefix_nu = original.shared_prefix_nu;
		if (original.shared_prefix_nu != original.prefix_nu){
			prefix_nu = new int[vcell->p + 1];
			prefix_nu[0] = 0;
			for (int i = 1; i < vcell->p + 1; i++){
				prefix_nu[i] = prefix_nu[i-1] + vcell->nu[i-1];
			}
		}else{
			prefix_nu = shared_prefix_nu;
		}
	}
	return *this;
}

void CellIterator::reset(){
	for (size_t i = 0; i < visited.size(); i++){
		visited[i] = false;
	}
	v = 0;
	j = 0;
}

inline voro::EdgeDesc CellIterator::current() const{
	return EdgeDesc(v, j, prefix_nu[v] + j);
}

inline bool CellIterator::isFinished(){
	return (current() >= (int)visited.size());
}

inline bool CellIterator::jump(){
	if (isFinished()){
		return false;
	}
	++j;
	if (j >= vcell->nu[v]){ //no edges, jump to next vertex
		j = 0;//first edge of new vertex
		++v;//new vertex
	}
	while (!isFinished() && isMarked()){
		++j;
		if (j >= vcell->nu[v]){ //no edges, jump to next vertex
			j = 0;//first edge of new vertex
			++v;//new vertex
		}
	}
	return !isFinished();
}

inline voro::EdgeDesc CellIterator::next() const{
	int u = vcell->ed[v][j];
	int k = vcell->cycle_up(vcell->ed[v][j + vcell->nu[v]], u); //6.12.2013 ... bylo na koncu v); zamiast u);
	return EdgeDesc(u, k, voroToEdge(u, k));
}

inline voro::EdgeDesc CellIterator::prev() const{
	int rev = vcell->cycle_down(j, v);
	int z = vcell->ed[v][rev];
	int i = vcell->ed[v][vcell->nu[v] + rev];
	return EdgeDesc(z, i, voroToEdge(z, i));
}

inline voro::EdgeDesc CellIterator::inv() const{
	int u = vcell->ed[v][j];
	int l = vcell->ed[v][vcell->nu[v] + j];
	return EdgeDesc(u, l, voroToEdge(u, l));
}

inline voro::EdgeDesc CellIterator::cw() const{
	int l = vcell->cycle_down(j, v);
	return EdgeDesc(v, l, voroToEdge(v, l));
}

inline voro::EdgeDesc CellIterator::ccw() const{
	int l = vcell->cycle_up(j, v);
	return EdgeDesc(v, l, voroToEdge(v, l));
}

inline CellIterator& CellIterator::forward(){
	EdgeDesc desc = next();
	v = desc.v;
	j = desc.j;
	return *this;
}

inline CellIterator& CellIterator::backward(){
	EdgeDesc desc = prev();
	v = desc.v;
	j = desc.j;
	return *this;
}

inline CellIterator& CellIterator::reverse(){
	EdgeDesc desc = inv();
	v = desc.v;
	j = desc.j;
	return *this;
}

inline CellIterator& CellIterator::go(int v, int j){
	this->v = v;
	this->j = j;
	return *this;
}

inline CellIterator& CellIterator::go(const EdgeDesc& edge){
	this->v = edge.v;
	this->j = edge.j;
	return *this;
}

inline CellIterator& CellIterator::rotateCW(){
	EdgeDesc desc = cw();
	v = desc.v;
	j = desc.j;
	return *this;
}

inline CellIterator& CellIterator::rotateCCW(){
	EdgeDesc desc = ccw();
	v = desc.v;
	j = desc.j;
	return *this;
}

inline bool CellIterator::isMarked() const{
	return visited[voroToEdge(v, j)];
}

inline bool CellIterator::isMarked(int edge) const{
	return visited[edge];
}

inline CellIterator& CellIterator::mark(int edge){
	visited[edge] = true;
	return *this;
}

inline CellIterator& CellIterator::mark(){
	visited[voroToEdge(v, j)] = true;
	return *this;
}

inline CellIterator& CellIterator::unmark(int edge){
	visited[edge] = true;
	return *this;
}

inline CellIterator& CellIterator::unmark(){
	visited[voroToEdge(v, j)] = false;
	return *this;
}

inline dmga::model::Point3D CellIterator::vertex() const{
	double* point = vcell->pts + v * 3;
	return dmga::model::Point3D(*point, point[1], point[2]);
}

inline dmga::model::Point3D CellIterator::vertex(const EdgeDesc& edge) const{
	double* point = vcell->pts + edge.v * 3;
	return dmga::model::Point3D(*point, point[1], point[2]);
}

inline dmga::model::Point3D CellIterator::vertex(int v) const{
	double* point = vcell->pts + v * 3;
	return dmga::model::Point3D(*point, point[1], point[2]);
}

inline double* CellIterator::vertexRaw() const{
	return vcell->pts + v * 3;
}

inline double* CellIterator::vertexRaw(const EdgeDesc& edge) const{
	return vcell->pts + edge.v * 3;
}

inline double* CellIterator::vertexRaw(int v) const{
	return vcell->pts + v * 3;
}

inline voro::voronoicell_base* CellIterator::cell(){
	return this->vcell;
}

/** returns cell object */
inline const voro::voronoicell_base* CellIterator::cell() const {
	return this->vcell;
}

inline int* CellIterator::makeShared(){
	if (prefix_nu == shared_prefix_nu){
		return 0;
	}
	shared_prefix_nu = prefix_nu;
	return prefix_nu;
}

CellIterator::~CellIterator(){				DEB3("CellIterator::__destruct__():");
	if(prefix_nu != shared_prefix_nu)
		dmga::utils::safeDeleteArray(prefix_nu);
	vcell = 0;
	shared_prefix_nu = 0;
}

}//namespace voro

#endif /* UTILS_HPP_ */
