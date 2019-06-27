/*
 * alpha_shape.h
 *
 *  Created on: Oct 9, 2013
 *      Author: robson
 */

#ifndef ALPHA_SHAPE_H_
#define ALPHA_SHAPE_H_

#include <shape/shape.h>
#include <model/cell.hpp>
#include <model/cell_ex.hpp>
#include <algorithm>
#include <ostream>

//TODO: poprawić kod, żeby był bardziej zgodny z nową strukturą shapeów

namespace dmga{

namespace shape{

/**
 * This is a basic building block for alpha-complexes (and alpha shapes)
 * it is strongly connected with the implementation of Voro++ cells
 */
class AlphaSimplex{
public:
	/** simplex dimension, may be 1, 2, 3 */
	signed char dim;
	/** vertex of the cell */
	int v;
	/** edge at the vertex - not used when 3-d simplex */
	int j;
	/** alpha threshold at which it appears in the complex */
	double alpha_threshold;
	/** TODO: ??? is this true: id in apropriate not-sorted table inside CellComplex */
	int id;
	/** ONLY IN DEBUG */
	bool is_killed;

	AlphaSimplex() :
			dim(-1),
			v(-1),
			j(-1),
			alpha_threshold(-1.0),
			id(-1),
			is_killed(false){
	}

	AlphaSimplex(char dim, int v, int j, double alpha_threshold, int back_id, bool is_killed = false) :
		dim(dim),
		v(v),
		j(j),
		alpha_threshold(alpha_threshold),
		id(back_id),
		is_killed(is_killed){
	}

	/**
	 * for sorting
	 */
	friend bool operator<(const AlphaSimplex& a, const AlphaSimplex& b){
		return (a.alpha_threshold < b.alpha_threshold) ||
				(a.alpha_threshold == b.alpha_threshold && a.dim < b.dim);
	}
	/**
	 * for searching
	 */
	friend bool operator<(double alpha, const AlphaSimplex& b){
		return (alpha < b.alpha_threshold);
	}
	/**
	 * for searching
	 */
	friend bool operator<(const AlphaSimplex& b, double alpha){
		return (b.alpha_threshold < alpha);
	}
	/**
	 * output simplex in format: {dim,v,j,alpha,id}
	 */
	friend std::ostream& operator<<(std::ostream& out, const AlphaSimplex& simplex){
		out << "{" << (int)simplex.dim << "," << simplex.v << "," << simplex.j << "," << simplex.alpha_threshold << ","<< simplex.id << "}";
		return out;
	}
};

/**
 *
 */
class AlphaComplex : DMGAObject{
public:
	std::string show(){ return "dmga::shape::AlphaComplex"; }

	typedef std::vector<AlphaSimplex> StorageType;
	typedef typename StorageType::iterator IteratorType;
	StorageType simplices;
	std::vector<int> back_ref[4];
	double max_alpha_threshold;

	AlphaComplex():
		max_alpha_threshold(0.0){										DEB3("AlphaComplex::__construct__(): ");
	}
	/**
	 * this should be private
	 */
	void makeBackRefs(int verticesMax, int edgesMax, int sidesMax){		DEB3("AlphaComplex::makeBackRefs(): vertices = " << verticesMax << ", edges = "<< edgesMax << ", sides = "<< sidesMax);
		back_ref[1].resize(sidesMax); //dim 1
		back_ref[2].resize(edgesMax); //dim 2
		back_ref[3].resize(verticesMax); //dim 3
		auto it = simplices.begin();
		auto end = simplices.end();
		int n = 0;
		for (; it != end; ++it, ++n){
			back_ref[(*it).dim][(*it).id] = n;
		}
	}

	AlphaSimplex& getSimplex(int dim, int id){							DEB3("AlphaComplex::getSimplex(): dimension = " << dim << ", id = " << id);
		return simplices[back_ref[dim][id]];
	}

	AlphaSimplex& operator()(int dim, int id){							DEB3("AlphaComplex::operator()(): dim = " << dim << ", id = " << id);
		return getSimplex(dim, id);
	}

	void add(const AlphaSimplex& new_simplex){							DEB3("AlphaComplex::add(): " << new_simplex);
		if (new_simplex.alpha_threshold > max_alpha_threshold){
			max_alpha_threshold = new_simplex.alpha_threshold;			DEB3("AlphaComplex::add(): Max alpha shape changed... ");
		}
		simplices.push_back(new_simplex);
	}

	AlphaSimplex& item(int number){										DEB3("AlphaComplex::item(): ");
		return simplices[number];
	}

	AlphaSimplex& operator[](int number){								DEB3("AlphaComplex::operator[](): ");
		return item(number);
	}

	IteratorType begin(){												DEB3("AlphaComplex::begin(): ");
		return simplices.begin();
	}

	IteratorType end(){													DEB3("AlphaComplex::end(): ");
		return simplices.end();
	}

	IteratorType end(double alpha){										DEB3("AlphaComplex::end(): alpha = " << alpha);
		return std::upper_bound(simplices.begin(), simplices.end(), alpha);
	}

	friend std::ostream& operator<<(std::ostream& out, const AlphaComplex& complex){
		auto it = complex.simplices.begin();
		auto end = complex.simplices.end();
		if (it != end){
			out << (*it);
			++it;
		}
		for (; it != end; ++it){
			out << " " << (*it);
		}
		out << "\n";
		return out;
	}

	~AlphaComplex(){ DEB3("AlphaComplex::__destruct__():");
	}
};


/**
 * w tej klasie umieścic wszystko co bedzie umiała robić jedna komórka
 * AlphaShape - sprawdzać wolne krawędzie, wierzchołki, ściany, wypisywac sie na stdout itp
 */
template<typename CellSpec>
class AlphaShapeCell : public DMGAObject{
public:
	std::string show(){ return "dmga::shape::AlphaShapeCell"; }

	typedef CellSpec CellType;

	CellType* cell;
	AlphaComplex* complex;

	AlphaShapeCell(CellType* cell = 0, AlphaComplex* complex = 0) :
			cell(cell), complex(complex){
	}

	AlphaShapeCell(CellType& cell, AlphaComplex& complex) :
			cell(&cell), complex(&complex){
	}

	AlphaShapeCell(CellType& cell) :
			cell(*cell), complex(0){
	}

	/**
	 * serializuj do JSONa
	 */
	std::string serialize(){
		std::ostringstream out;
		return out.str();
	}

	friend std::ostream& operator<<(std::ostream& out, const AlphaShapeCell& shape){
		//TODO: przepisac render, aby potrafił korzystac z JSONa
		auto& cell = *shape.cell;
		auto& vcell = cell.voropp_cell;
		auto& complex = *shape.complex;
		voro::CellIterator cell_it = voro::CellIterator(&vcell);
		double* point = vcell.pts;
		double* pend = vcell.pts + vcell.p * 3;
		out << vcell.p << "\n";
		int v = 0, j = 0;
		int dim3, dim2, dim1;
		double alpha_threshold = 0.0;
		bool is_killed;
		for (; point < pend; point += 3, ++v){
			dim3 = v;
			alpha_threshold = complex.getSimplex(3, dim3).alpha_threshold;
			out << (*point) << " " << (point[1]) << " " << (point[2]) << " {'alpha':" << alpha_threshold << "}\n";
		}
		out << cell.size() << "\n";
		int current_side = 0;
		while (!cell_it.isFinished()){
			v = cell_it.current().v;
			j = cell_it.current().j;
			dim1 = current_side;
			alpha_threshold = complex.getSimplex(1, dim1).alpha_threshold;
			out << "{'alpha':" << alpha_threshold << ",'killed':" << is_killed << "}"; //side threshold

			while (!cell_it.isMarked()){
				v = cell_it.current().v;
				j = cell_it.current().j;
				dim2 = voro::edges::voroToEdge(cell.prefix_nu, v, j);
				alpha_threshold = complex.getSimplex(2, dim2).alpha_threshold;
				is_killed = complex.getSimplex(2, dim2).is_killed;

				out << " " << cell_it.current().v  << " {'alpha':" << alpha_threshold << ",'killed':" << is_killed << "}";
				cell_it.mark();
				cell_it.forward();
			}
			out << "\n";
			current_side++;
			cell_it.jump();
		}
		return out;
	}

	~AlphaShapeCell(){ DEB3("AlphaShapeCell::__destruct__():");
	}
};

/**
 * non-kinetic shape of the set of balls in 3D (for power diagram)
 * it uses Voronoi (Power) diagram directly, not as standard implementations the dual - delaunay (regular) triangulation
 * this way it is simpler to deal with non-general positions of balls, and in our problems we pay much more attention to
 * voronoi diagrams than in normal applications where only proximity is of importance (thus delaunay is better)
 */
template<typename DiagramSpec>
class VoroppAlphaShape : public GenericShape{
public:
	std::string show(){ return "dmga::shape::VoroppAlphaShape"; }

	typedef DiagramSpec DiagramType;
	typedef typename DiagramType::ContainerType ContainerType;
	typedef typename DiagramType::GeometryType GeometryType;
	typedef typename DiagramType::ParticleType ParticleType;
	typedef typename DiagramType::KeyType KeyType;
	typedef typename DiagramType::CellType CellType;
	typedef typename DiagramType::ResultSetType ResultSetType;
	typedef typename DiagramType::SubsetType SubsetType;
	typedef typename DiagramType::CacheType CacheType;
	typedef AlphaShapeCell<CellType> CellShapeType;

	DiagramType& diagram;
	std::vector<AlphaComplex*> cell_complexes;
	double max_alpha_threshold;

	VoroppAlphaShape(DiagramType& diagram) : diagram(diagram), cell_complexes(diagram.container.size()), max_alpha_threshold(-1) {		DEB3("VoroppAlphaShape::__construct__(): ");
		recompute();
	}

	AlphaComplex* computeOne(int i){					DEB3("VoroppAlphaShape::computeOne(): compute cell " << i);
		CellType& cell = *(diagram.getCell(i));

		//iterate over all sides in cell - each side defines edge (1d-simplex) in dual triangulation
		//iterate over all edges in cell - each edge defines side (2d-simplex) in dual triangulation
		//iterate over all vertices in cell - each vertex defines tetrahedron (3d-simplex) in dual triangulation
		//all the above things you can do when iterating over the skeleton (all edges) of the cell, using bare CellIterator
		//this way it is faster than using 'nice' collection .sides, .vertices, .edges from CellType

		auto cell_it = cell.getCellIterator();
		std::vector<bool> created_vertices(cell.getVerticesCount(), false); //here we mark which vertices alpha_thresholds we already computed
		std::vector<double> vertex_alpha_thresholds(cell.getVerticesCount(), 0.0); // here we will cache alpha_thresholds for vertices to avoid recomputing
		int current_side = 0; //current side, it's used to generate id for each side.
		cell_complexes[cell.number] = new AlphaComplex();
		while (!cell_it.isFinished()){
			DEB3("VoroppAlphaShape::computeOne(): \tstart new side");

			// we will use those single letter (bad coding... sorry) variables to describe:
			int v = cell_it.current().v;	// v - current vertex on the side
			int j = cell_it.current().j;	// j - current edge on the side
			int u = cell_it.inv().v; 		// u - next vertex on the side
			int k = cell_it.next().j;		// k - the next edge on the side (from u)
			int w = cell_it.prev().v;		// w - the previous vertex on the side

			// now we are at new side of the voronoi diagram. We need to gather some information about
			// this side in order to compute it correctly. We now compute alpha_threshold, but we remember it only
			// we can use it only if no other higher dimensional simplex 'kills' the side. We need some values
			// to test if something is killing our side before we add simplex.
			//
			// So you can build now the data needed for the side simplex (1d) - two neighbouring particles p_1 and p_2
			// the distance is computed simply by finding weighted point at which line connecting p_1 and p_2
			// crosses this side. It is simply computable as it is the scaling factor by which we move
			// the plane to cut the cell. The computation of the scaling factor is (p_1 = cell.number, p_2 the other vertex):
			// dist(p_1, p) = (dist(p_1, p_2)^2 + (r_1^2 - r_2^2)) / (2.0 * dist(p_1, p_2));
			// d = dist(p_1, p)^2 - r_1^2
			//
			// we of course need to consider 'killing' of the given alpha_threshold by other, higher dimensional simplices.
			// By 'killing' we mean, that the alpha-ball defined by the two vertices which we use to compute alpha-threshold
			// is not the real alpha-ball which expose the two particles...
			// for a side it means that it can be killed in two ways:
			// 1) by single other particle, that will lie inside the alpha-(power)-ball defined by two particles
			//    we say that the ball is killed by the voronoi edge (defined by three particles)
			// 2) by two particles forming with the our two particles the 3d-simplex
			// REMARK: other cases are impossible as it would violate the definition of the current side.
			// the two cases may be blend to one repair procedure, by seeing that if we are in case 2 then it
			// must be that alpha-ball at the some edge was also killed by those 4 particles. Now if we first
			// resolve killing of the alpha-balls at edges then repair the alpha ball at the side adjecent to that edge
			// then we will need only to care about the first case
			//
			// we can detect if the repair is needed if we consider the point p of intersection of
			// the plane containing the side and the line connecting two particles. This point can
			// Be easily computed (see comments below)
			// If this point lies inside the side then no repair is needed. Otherwise we need to find
			// the closest edge to that point and use its alpha-threshold.
			// The test may be done using the signed area of the triangle defined by the edges of the side and the point p.
			// if all signs are the same (assume positive) then the point lies inside the side. Otherwise there should be
			// some edges that form triangle with 'negative' area with point p. There should be exactly one such an edge e
			// that the projection of the p onto the line containing e belongs to e
			// (in other words the triangle formed with e and p is 'nice', not rozwartokatny, we can even check rozwartokatnosc ;) ).
			// we choose the alpha-threshold of this edge e.
			// the above discussion is in exact match with the definition of alpha-complex by growing balls, as this edge e is the edge
			// which belong to the side and is first hit by the growing balls in the two defining points, when we reach the alpha-threshold.

			DEB3("VoroppAlphaShape::computeOne(): \tcomputing side simplex (edge in Delaunay)");
			double point_on_side[3]; // this point will be used to test if side alpha-threshold needs to be repaired (is killed by something)
			double side_alpha_threshold; //we need to store it for adding at the end
			// we need normal to the side, to compare it in the various tests
			// for example in testing if crossing point of delaunay edge is hitting plane inside side or outside
			// or to test if the current side is created by periodic or non-periodic particle
			double side_normal[4];
			//we compute alpha_threshold using the distance from the side - this way we will automatically deal with the periodic sides
			double* particle_coords = diagram.container.getRawCoords(cell.number);
			// we will use normal to compute the point on the plane defined by side where the alpha-ball first hit the plane,
			// it will be used to test if alpha-threshold need to pe repaired
			// we get a side_normal as a side effect of the computation of the euclidean distance to the plane
			double eucl_dist = algebra::distance::euclid3d(particle_coords,
															cell.getVertexRaw(v),
															cell.getVertexRaw(u),
															cell.getVertexRaw(w),
															side_normal); //here the sequence of v, u, w is important, the same as in computing side_normal
			side_alpha_threshold = eucl_dist * eucl_dist - particle_coords[3] * particle_coords[3];
			//remember point on the side
			point_on_side[0] = particle_coords[0] - side_normal[0] * eucl_dist;
			point_on_side[1] = particle_coords[1] - side_normal[1] * eucl_dist;
			point_on_side[2] = particle_coords[2] - side_normal[2] * eucl_dist;
			//and prepare some variables used to test 'killing'
			bool side_is_killed = false;	//we assume at first that no edge is 'killing' the side
			double dist_from_killing_edge = 0; // we need only to care about edges in Voronoi that can kill side, because if vertex kills side then it also kills edge
			double killing_edge_threshold = 0; // here we will remember the killer alpha-threshold to replace previously computed alpha_threshold for the side

			// we compute the alpha threshold for the first vertex on this side (if it is neccessary):
			if (!created_vertices[v]){
				DEB3("VoroppAlphaShape::computeOne(): \tcomputing vertex simplex (tetrahedron in Delaunay) [place: first vertex on side]");
				vertex_alpha_thresholds[v] = computeVertexAlphaThreshold(cell.getVertexRaw(v), diagram.container.getRawCoords(cell.number));
				created_vertices[v] = true;
				cell_complexes[cell.number]->add(AlphaSimplex(3, v, j, vertex_alpha_thresholds[v], v));
			}

			//start iteration on the curent side
			while (!cell_it.isMarked()){
				v = cell_it.current().v; // the v-th vertex
				j = cell_it.current().j; // j-th edge in vertex v
				u = cell_it.inv().v; // the end of the j-th edge in v is vertex u
				k = cell_it.next().j;

				// we have now the u-vertex, we can check if the tetrahedron simplex for u was created before
				// (when we were on the different side of the cell) if no, we create it. We always create the latter vertex of the
				// current edge, to assure we had all necessary information when we test for edge killing (se below)
				if (!created_vertices[u]){
					DEB3("VoroppAlphaShape::computeOne(): \tcomputing vertex simplex (tetrahedron in Delaunay) [place: in loop]");
					vertex_alpha_thresholds[u] = computeVertexAlphaThreshold(cell.getVertexRaw(u), diagram.container.getRawCoords(cell.number));
					created_vertices[u] = true;
					cell_complexes[cell.number]->add(AlphaSimplex(3, u, k, vertex_alpha_thresholds[u], u));
				}

				// we need to compute the edge twice, as it is on the two sides - to ensure both sides will compute their 'killing' posiible edges
				DEB3("VoroppAlphaShape::computeOne(): \tcomputing edge simplex (side triangle in Delaunay)");
				// those are the particles that form the 2d-simplex (triangle)

				// to create simplex we should gather below information, but, by duality we identify this with the edge
				// so we do not explicitely store this information in our AlphaComplex...
				// s[0] = cell.number;
				// s[1] = cell.getNeighbourId(u, i);
				// s[2] = cell.getNeighbourId(v, j);
				// std::sort(s, s + 3);
				// instead we remember the edge description (v, j), and compute the alpha_threshold
				// alpha_threshold can be computed as a power distance from the any of the
				// particles to the point of interseciton of the line containing edge
				// and a plane containing three particles
				// but this is the same as computing:
				// d = dist(p_i, l)^2 - r_i^2 for any particle p_i from the triangle and line l containing edge
				// by the definition of the edge and the power distance
				// computing dist(p_i, l) is easy in our setting because
				// 1) we have line given by the edge = interval in 3d given by two points (v, u)
				// 2) dist(p_i, l) = h - the height of the triangle on the base (v, u)
				// 3) area_of_the_triangle = h * len(v, u) / 2 =>
				//    h = 2 * area_of_the_triangle / ||u - v||
				//    h = |(p_i - u) x (p_i - v)| / ||u - v|| :)

				// THINK: maybe from the points p_1, p_2, p_3 choose smaller index - always get the same distance, whatever you compute for cell p_1 or p_2 or p_3
				// THINK: if change the particle, then choose the biggest number - it will automatically solve the problem of the walls (which have negative ids)

				// We need also to test if one of the vertices 'kills' the alpha-ball. It is when the two vertices u and v lies on the same side of the plane defined by
				// three particles (or more generally - for walls - the plane perpendicular to the edge). However, we only need to do simple test on the one particle and two
				// vertices. Namely, we need to test if the two angles between the edge and lines connecting vertices to one of the particles are smaller than pi/2. It can be
				// measured as the sign of the dot product of the appropriate vectors, as the dot product is cos(angle) * [something positive] and the cosine is positive
				// when the angle is lower than pi/2 or greater than 3*pi/2 which is perfect for our purpose.

				int choosen_particle = cell.number;
				double* v_coords = cell.getVertexRaw(v);
				double* u_coords = cell.getVertexRaw(u);
				double* p_coords = diagram.container.getRawCoords(choosen_particle); // we can safely assume that, if we choose real particle, then this coordinates will be good for the angles test for 'killing'
				double alpha_threshold;
				// first test if the two ends of the edge (i.e. v and u) lies on the opposite sides of the 2d-simplex defined by
				// the three particles. If not, then we need to get the alpha-threshold from the bad particle (the one on the bad
				// side of the plane, that is closer to the plane. Below there must be one or the other on none, but not both, so we have only three cases.
				bool edge_is_killed = false;
				if (!algebra::is_acute_angle(p_coords, u_coords, v_coords)){
					// need to use alpha-threshold for u, we have them computed and stored in vertex_alpha_thresholds
					alpha_threshold = vertex_alpha_thresholds[u];
					edge_is_killed = true;
				}else if (!algebra::is_acute_angle(p_coords, v_coords, u_coords)){
					// need to use alpha-threshold for v, we have them computed and stored in vertex_alpha_thresholds
					alpha_threshold = vertex_alpha_thresholds[v];
					edge_is_killed = true;
				}else{
					// we need to simply compute the alpha_threshold
					double n[4];
					algebra::normal3d_raw(p_coords, v_coords, u_coords, n); //in n[3] we have the doubled area of the triangle
					alpha_threshold = ((n[3] * n[3]) / algebra::distance::euclid3d_pow_2(v_coords, u_coords)) - p_coords[3] * p_coords[3];
				}
				// now we have real alpha threshold computed, we can add simplices
				cell_complexes[cell.number]->add(AlphaSimplex(2, v, j, alpha_threshold, voro::edges::voroToEdge(cell.prefix_nu, v, j), edge_is_killed));
				//TODO: przemyslec czy nie zrobic cache'u jak w przypadku wierzcholkow... //cell_complexes[cell.number].add(AlphaSimplex(2, u, i, alpha_threshold, voro::edges::voroToEdge(cell.prefix_nu, u, i)));

				// now, we have real alpha-threshold for the edge. We must test if the edge is not killing the side simplex.
				// we have remembered those variables:
				// double point_on_side[3] - this point (on the plane defined by side) will be used to test if side alpha-threshold needs to be repaired
				// side_is_killed - if this is true then the variables below are defined, and we need to compare current edge to it
				// killing_v, killing_j - the description of the killing edge
				// dist_from_killing_edge - distance from the current killing edge, if our edge is closer then we must replace the current edge
				double test_normal[4];
				algebra::normal3d_raw(v_coords, u_coords, point_on_side, test_normal); //in n[3] we have the doubled area of the triangle
				// TODO: THINK, this is standard dot product, should I rewrite it using some inline function?
				if (algebra::dot3d(side_normal, test_normal) < 0){
					// this is 'killing' edge
					// in test_normal[3] we have the doubled area of the triangle we need
					// to compute its height - the distance from the point on the side!
					// and the test-distance is the height of the triangle created with v, u and point_on_side, that is perpendicular to (v, u)
					// THIS TEST IS CERTAINLY GOOD - because we computed test_normal putting particles in the same order as in side_normal
					double test_dist = (test_normal[3] * test_normal[3]) / algebra::distance::euclid3d_pow_2(v_coords, u_coords);
					if (!side_is_killed || test_dist < dist_from_killing_edge){
						dist_from_killing_edge = test_dist;
						killing_edge_threshold = alpha_threshold;
						side_is_killed = true;
					}
				}

				cell_it.mark();
				cell_it.forward();
			}// for all edges in this side

			// before we jump to the next side, we check if we need repair 'killed' alpha-threshold for this side
			if (side_is_killed){
				side_alpha_threshold = killing_edge_threshold;
			}
			// finally ad the AlphaSimplex for the side
			// we can use as the base vertex and the base edge the last remembered v and j
			cell_complexes[cell.number]->add(AlphaSimplex(1, v, j, side_alpha_threshold, current_side, side_is_killed));

			//now we can jump to the next side
			current_side++;
			cell_it.jump();
		}
		// now we constructed all simplices and it is time to sort them for this particular cell
		std::sort(cell_complexes[cell.number]->begin(), cell_complexes[cell.number]->end());
		int vertices_count = cell.getVerticesCount();
		cell_complexes[cell.number]->makeBackRefs(vertices_count, cell.prefix_nu[vertices_count], current_side);
		double complex_max_alpha_threshold = cell_complexes[cell.number]->max_alpha_threshold;
		if (complex_max_alpha_threshold > max_alpha_threshold) max_alpha_threshold = complex_max_alpha_threshold;
		DEB3("VoroppAlphaShape::computeOne(): side " << i << " finished.");
		return cell_complexes[cell.number];
	}
	/**
	 * recomputes the whole alpha-shape for a diagram
	 * it iterates over all cells currently computed and stored in diagram.cache
	 */
	void recompute(){										DEB3("VoroppAlphaShape::recompute(): ");
//		for each cell c in diagram
//			for each 'pseudo-simplex' s in c
//				//pseudo-simplex is each hypothetic simplex defined by vertexes, edges and sides of c
//				compute alpha_threshold t for s
//				store description of simplex s with t
//			sort simplexes with respect to t and assign this list to cell c
//			now for alpha > 0 the prefix of this sorted array with t <= alpha is alpha-complex
//		for each c we have a sorted list of simplexes - thus we can for each cell retrieve all
		auto cache_it = diagram.cache.begin();
		auto cells_end = diagram.cache.end();
		for (; cache_it != cells_end; ++cache_it){
			if (*cache_it){
				computeOne((*cache_it)->number);
			}
		} // foreach cell c
		// now we constructed all simplices for all cells, we have them in the AlphaComplex
		DEB3("VoroppAlphaShape::recompute(): finished");
	}

	AlphaComplex* getAlphaComplex(int i){					DEB3("VoroppAlphaShape::getAlphaComplex(): i = " << i);
		if (!cell_complexes[i]){
			//TODO: test if cell computed
			this->computeOne(i);
		}
		return cell_complexes[i];
	}

	AlphaComplex* getAlphaComplex(const CellType& cell){	DEB3("VoroppAlphaShape::getAlphaComplex(): cell.number = " << cell.number);
		if (!cell_complexes[cell.number]){
			//TODO: test if cell belongs to this diagram????
			this->computeOne(cell.number);
		}
		return getAlphaComplex(cell.number);
	}

	CellShapeType* getCellShape(int i){						DEB3("VoroppAlphaShape::getCellShape(): i = " << i);
		CellType* cell = diagram.cache[i];	//TODO: diagram.computeCell(i)? or cache?
		AlphaComplex* complex = this->getAlphaComplex(i);
		return new CellShapeType(cell, complex);
	}

	CellShapeType* getCellShape(CellType& cell){			DEB3("VoroppAlphaShape::getCellShape(): cell.number = " << cell.number);
		return new CellShapeType(cell, *this->getAlphaComplex(cell.number));
	}

	~VoroppAlphaShape(){ DEB3("VoroppAlphaShape::__destruct__():");
	}

private:
	/**
	 * each vertex defines tetrahedron (3d-simplex) in dual triangulation
	 */
	double computeVertexAlphaThreshold(double* vertex_coords, double* particle_coords){
		// compute the power distance from v do any of p_i's (one of them is a parameter to this function, we always can simply use cell.number)
		// d = dist(v, p_i)^2 - r_i^2
		return algebra::distance::euclid3d_pow_2(vertex_coords, particle_coords) - particle_coords[3] * particle_coords[3];
	}
};

template<typename DiagramSpec>
class VoroppKineticShape : public GenericKineticShape{
public:
	std::string show(){ return "dmga::shape::VoroppKineticShape"; }

	typedef DiagramSpec DiagramType;

	DiagramType& diagram;

	VoroppKineticShape(DiagramType& diagram) : diagram(diagram){

	}

	~VoroppKineticShape(){ DEB3("VoroppKineticShape::__destruct__():");
	}
};


} //namespace shape

} //namespace dmga


#endif /* ALPHA_SHAPE_H_ */
