/*
 * sasa_shape.h
 *
 *  Created on: Oct 9, 2013
 *      Author: robson
 */

#ifndef SASA_SHAPE_H_
#define SASA_SHAPE_H_

#include <vector>
#include <model/primitives.hpp>
#define _USE_MATH_DEFINES 1 //for PI
#include <math.h> //for PI and cos/acos functions

namespace dmga{

namespace shape{

using namespace model;

/**
 * this class holds the description of the surface of the sphere
 * laying outside a voronoi cell. It can compute its area then
 * we can find the SASArea of inside by substracting this from 4 * \pi * R^2
 */
class SASAContours : DMGAObject{
public:
	std::string show(){ return "dmga::shape::SASAContours"; }

	/**
	 * here we store information of the intersection of the line connecting
	 * two particles with appropriate elements (sphere, side, etc.)
	 * it is used to store cap area and volume if there is no intersection
	 * of sphere with side border.
	 */
	class CenterDesc{
	public:
		int neighbour_id;
		Point3D point_on_sphere;
		Point3D point_on_plane;
		double cap_intersection_area; // > 0 only if we have clear intersection with side - not side edges
		double cap_intersection_volume; // > 0 only if we have clear intersection with side - not side edges
		double cap_radius_sq;
		CenterDesc():
				neighbour_id(-1),
				point_on_sphere(0,0,0),
				point_on_plane(0,0,0),
				cap_intersection_area(0.0),
				cap_intersection_volume(0.0),
				cap_radius_sq(0.0){
		}
	};
	/**
	 * this class represents data of the next element on border of
	 * the cap area or polygon area. We also store appropriate center
	 * id, as polygons and border arcs always are based on some center.
	 * moreover, separating centers allows for very simple representation
	 * of borders and polygons as each vertex on border has exactly one next on
	 * border and exactly one on polygon neighbour (thats not true for centers, as
	 * they can be part of any number of polygons/borders)
	 */
	class VertexLink{
	public:
		int vertex;
		int center;
		VertexLink(): vertex(-1), center(-1) {}
	};
	/**
	 * this class represents a single vertex on the boundary of the cap region
	 * as we mentioned in @see VertexLink, we can only remember one neighbour
	 * on border of the region and one neighbour on the polygon.
	 */
	class BorderVertex : public Point3D{
	public:
		typedef Point3D super;
		VertexLink next_on_border;
		VertexLink next_on_polygon;
		bool visited; //TODO: this shoud be dumped, and a iterator should be created for the purpose of iterating thru this structure
		BorderVertex() : super(0,0,0), visited(false){
		}
		BorderVertex(double x, double y, double z): super(x, y, z), visited(false){
		}
	};
	class EdgePoints{
	public:
		enum EdgeType{
			EDGE_TYPE_UNKNOWN, 	// not computed yet
			EDGE_TYPE_BURRIED, 	// whole edge under the dome
			EDGE_TYPE_FREE, 	// whole edge outside the dome
			EDGE_TYPE_DOUBLE, 	// two intersection points
			EDGE_TYPE_FAR, 		// only one  intersection, base vertex burried under dome
			EDGE_TYPE_NEAR, 	// only one  intersection, base vertex outside the dome
		};

		int near_point;
		int far_point;
		EdgeType edge_type;
		// if point1 or point2 < 0 then appropriate point not exist
		EdgePoints(): near_point(-1), far_point(-1), edge_type(EDGE_TYPE_UNKNOWN) {}
	};
	class SASAContourIterator{
	public:
		SASAContours* parent;
		int vertex_count;
		std::vector<bool> visited;
		int current;
		int next;
		int center;
		SASAContourIterator():
				parent(0),
				vertex_count(0),
				current(-1),
				next(-1),
				center(-1){
			//creates empty contour iterator
		}
		SASAContourIterator(SASAContours* parent):
				parent(parent),
				vertex_count(parent->vertices.size()),
				visited(vertex_count, false),
				current(0),
				next(-1),
				center(-1){
		}
		virtual int nextVertex(int current) = 0;
		virtual int nextCenter(int current) = 0;
		bool isMarked()					{ return visited[current]; }
		bool isMarked(int vertex_id)	{ return visited[vertex_id]; }
		void mark()						{ visited[current]=true; }
		void unmark()					{ visited[current]=false; }
		void mark(int vertex_id)		{ visited[vertex_id]=true; }
		void unmark(int vertex_id)		{ visited[vertex_id]=false; }

		const Point3D& currentPoint(){
			return parent->vertices[current];
		}
		const Point3D& nextPoint(){
			return parent->vertices[next];
		}
		const SASAContours::CenterDesc& centerDesc(){
			return parent->centers[center];
		}
		const Point3D& centerPointOnPlane(){
			return parent->centers[center].point_on_plane;
		}
		const Point3D& centerPointOnSphere(){
			return parent->centers[center].point_on_sphere;
		}
		/*
		 * test if this conotur is empty, one should test this prior to use
		 * this iterator
		 * if return true then not empty and ready for iteration (current vertex at 0)
		 * if false then empty or current vertex not at 0 and should not iterate
		 * you may call reset() if not empty to start new iteration
		 */
		bool ready(){
			if (parent == 0 || current != 0 || vertex_count == 0){
				return false;
			}
			next = nextVertex(current);
			center = nextCenter(current);
			return true;
		}
		/**
		 * this resets iteration and returns true as in ready
		 * (so it is safe to call reset instead of ready, but it is slower (resets visisted flags and so on))
		 * so you may like to call ready first time, then use reset in subsequent calls
		 */
		bool reset(){
			for (unsigned int i = 0; i < visited.size(); i++) visited[i] = false;
			current = 0;
			return ready();
		}
		/**
		 * go to next vertex on current contour
		 *
		 * (beware: it cycle thru this contour indefinitely if you do not use jump!)
		 *
		 * returns current vertex
		 */
		int forward(){
			current = next;
			next = nextVertex(current);
			center = nextCenter(current);
			return current;
		}
		/**
		 * jump to the next unvisited vertex
		 *
		 * returns current vertex, or -1 if end of container
		 */
		int jump(){
			while (current < vertex_count && isMarked()) ++current;
			if (current >= vertex_count) return -1;
			next = nextVertex(current);
			center = nextCenter(current);
			return current;
		}
		/**
		 * iterate over all contours of this SASAContours if new contour started
		 * then it returns 1, if no new centour then 0, if end of container is reached then -1
		 *
		 * it should be used in do ... while loop like:
		 * do { smth } while (iterator.iterate() != -1)
		 */
		int iterate(){
			mark();
			if (isMarked(next)){
				return jump() == -1 ? -1 : 1;
			}else{
				forward();
				return 0;
			}
		}
	};
	//ABANDONED CONCEPT AFTER PLACING CAPS ALONG WITH THE NORMAL BORDER ARCS
//	/**
//	 * this class simply iterates over all arcs in this Contours container
//	 * at each step it returns four points: first, on_plane, on_sphere, second (all Point3D classes)
//	 * if first == second then we have cap (case 0), otherwise we have two arcs going from
//	 * first to on_sphere and from on_sphere to second
//	 * please note that those arcs are in
//	 * 1) some border arc (that is arc from first to second with radius dist(first, on_plane)
//	 *                     and lying on the plane given by first, on_plane, second)
//	 * 2) the polygon arcs (themselves - two great arcs from first to on_sphere, and from on-sphere to second).
//	 *  if you wan to draw outline of of the intersection of the voronoi cell and sphere then you
//	 *  need only to draw caps (case 0) and case 1)
//	 *
//	 *  it is DMGAObject as we will use it to iterate over contours in pydmga
//	 */
//	class SASAArcIterator : public DMGAObject{
//	public:
//		std::string show(){ return "dmga::shape::SASAArcIterator"; }
//		~SASAArcIterator(){ DEB3("SASAArcIterator::__destruct__():");
//		}
//		SASAContours& contours;
//		SASAArcIterator(SASAContours& parent) : contours(parent){
//			prepare();
//		}
//		struct ArcData{
//			Point3D first;
//			Point3D on_plane;
//			Point3D on_sphere;
//			Point3D second;
//			ArcData(): first(0,0,0), on_plane(0,0,0), on_sphere(0,0,0), second(0,0,0) {}
//			ArcData(const Point3D& fir, const Point3D& pla, const Point3D& sph, const Point3D& sec) : first(fir), on_plane(pla), on_sphere(sph), second(sec){}
//			ArcData(double* fir, double*  pla, double*  sph, double*  sec) : first(fir), on_plane(pla), on_sphere(sph), second(sec){}
//		};
//		virtual void prepare() {};
//		virtual bool hasNext() = 0;
//		virtual ArcData next() = 0;
//	};
//	class SortedSASAArcIterator : public SASAArcIterator{
//		//if is sorted then we will sort Centers by neighbour_id
//		//and we will sort BorderVertex by next_on_border.center.neighbour_id
//		//thus we will have all elements for given neighbour at a time
//		//morever we will sort borderVertices clockwise, so one will be
//		//able to connect them to form skeleton of a cell only using contours
//		//this may be however slow, so we will give option to skip sorting, and
//		//use unsorted things.
//		size_t i, j;
//		bool hasNext(){
//			if (i < contours.vertices.size()) return true;
//			while (j < contours.centers.size() && contours.centers[j].cap_radius_sq == 0) ++j;
//			return (j < contours.centers.size() && contours.centers[j].cap_radius_sq > 0);
//		}
//		ArcData next(){
//			if (i < contours.vertices.size()){
//				//we have Arc
//				BorderVertex first = contours.vertices[i];
//				CenterDesc center = contours.centers[first.next_on_border.center];
//				BorderVertex second = contours.vertices[first.next_on_border.vertex];
//				return ArcData(first, center.point_on_plane, center.point_on_sphere, second);
//			}
//			if (j < contours.centers.size() && contours.centers[j].cap_radius_sq > 0){
//				//we have cap
//				// CenterDesc center = contours.centers[j];
//				//we need to generate first == second on the plane in distance sqrt(center.rsq)
//			}
//			throw dmga::exceptions::StopIteration();
//		}
//		void prepare(){
//			i = 0;
//			j = 0;
//		}
//	};
//	class UnsortedSASAArcIterator : public SASAArcIterator{
//		bool hasNext(){
//			return false;
//		}
//		ArcData next(){
//			return ArcData();
//		}
//		void prepare(){
//
//		}
//	};

	class SASABorderIterator : public SASAContourIterator{
	public:
		typedef SASAContourIterator super;
		SASABorderIterator(): super() {};
		SASABorderIterator(SASAContours* parent): super(parent) {}
		int nextVertex(int current) { return parent->vertices[current].next_on_border.vertex; }
		int nextCenter(int current) { return parent->vertices[current].next_on_border.center; }
	};

	class SASAPolygonIterator : public SASAContourIterator{
	public:
		typedef SASAContourIterator super;
		SASAPolygonIterator(): super() {};
		SASAPolygonIterator(SASAContours* parent): super(parent) {}
		int nextVertex(int current) { return parent->vertices[current].next_on_polygon.vertex; }
		int nextCenter(int current) { return parent->vertices[current].next_on_polygon.center; }
	};
	std::vector<BorderVertex> vertices;
	std::vector<CenterDesc> centers;
	std::vector<EdgePoints> points_by_edge; // vertices[points_by_edge[edge_id].near_point] is a point near the vertex v of edge edge_id
	double sphere[4];
	double excluded_area;
	double excluded_volume;

	SASAContours():
			excluded_area(0.0),
			excluded_volume(0.0){
	}
	/**
	 * adds new data to a given side, it defines two points: on sphere, on plane (on side). If on_plane is further than on_side
	 * then on this side there exists the outside region that needs to be computed
	 */
	int addCenter(size_t side_id, int neighbour_id, double x_sphere, double y_sphere, double z_sphere, double x_plane, double y_plane, double z_plane){		DEB3("SASAContours::addCenter(): side_id = " << side_id << ", neighbour_id = " << neighbour_id);
		while (centers.size() <= side_id) centers.push_back(CenterDesc()); //TODO: you can initialize this outright in constructor!
		centers[side_id].neighbour_id = neighbour_id;
		centers[side_id].point_on_plane = Point3D(x_plane, y_plane, z_plane);
		centers[side_id].point_on_sphere = Point3D(x_sphere, y_sphere, z_sphere);
		return side_id;
	}
	/**
	 * returns data connected with a given edge (edge_id is an Voro++ edge id as defined in CellIterator)
	 */
	EdgePoints& getEdgePoints(size_t edge_id){																													DEB3("SASAContours::getEdgePoints(): edge_id = " << edge_id);
		while (points_by_edge.size() <= edge_id) points_by_edge.push_back(EdgePoints()); //TODO: you can initialize this outright in constructor!
		return points_by_edge[edge_id];
	}
	/**
	 * adds new near vertex (on forward edge). Near vertex is a vertex at which sphere goes outside of the voronoi cell
	 */
	int addNearVertex(size_t forward_edge_id, size_t backward_edge_id, double x, double y, double z){																DEB3("SASAContours::addNearVertex(): forward_edge = " << forward_edge_id << ", backward_edge = " << backward_edge_id);
		while (points_by_edge.size() <= forward_edge_id || points_by_edge.size() <= backward_edge_id) points_by_edge.push_back(EdgePoints()); //TODO: you can initialize this outright in constructor!
		int current_vertex_num = vertices.size();
		points_by_edge[forward_edge_id].near_point = current_vertex_num;
		points_by_edge[backward_edge_id].far_point = current_vertex_num;
		if (points_by_edge[forward_edge_id].edge_type == EdgePoints::EDGE_TYPE_FAR){
			points_by_edge[forward_edge_id].edge_type = EdgePoints::EDGE_TYPE_DOUBLE;
		}else{
			points_by_edge[forward_edge_id].edge_type = EdgePoints::EDGE_TYPE_NEAR;
		}
		if (points_by_edge[backward_edge_id].edge_type == EdgePoints::EDGE_TYPE_NEAR){
			points_by_edge[backward_edge_id].edge_type = EdgePoints::EDGE_TYPE_DOUBLE;
		}else{
			points_by_edge[backward_edge_id].edge_type = EdgePoints::EDGE_TYPE_FAR;
		}
		vertices.push_back(BorderVertex(x,y,z));
		return current_vertex_num;
	}
	/**
	 * adds new near vertex to this side. This vertex is a base for full cap on this side
	 * those vertices arise when no intersection with side edges is detected
	 * OnSide vertices are allweys connected to themselves
	 */
	int addOnSideVertex(int side_id, double x, double y, double z){																							DEB3("SASAContours::addOnSideVertex(): side_id = " << side_id);
		int current_vertex_num = vertices.size();
		vertices.push_back(BorderVertex(x,y,z));
		conditionalConnectBorder(current_vertex_num, current_vertex_num, side_id);
		conditionalConnectPolygon(current_vertex_num, current_vertex_num, side_id);
		return current_vertex_num;
	}
	/**
	 * adds new far vertex (on forward edge). Far vertex is a vertex at which sphere comes back to inside of the voronoi cell
	 */
	int addFarVertex(int forward_edge_id, int backward_edge_id, double x, double y, double z){																DEB3("SASAContours::addFarVertex(): forward_edge = " << forward_edge_id << ", backward_edge = " << backward_edge_id);
		return this->addNearVertex(backward_edge_id, forward_edge_id, x, y, z); //oh yeee I'm so cleva!
	}
	/**
	 * marks given edge as Burried edge (inside the sphere)
	 *
	 * Currently we do not use this in any way, but we remember this as a side effect of computation
	 */
	void addBurriedEdge(size_t forward_edge_id, size_t backward_edge_id){																							DEB3("SASAContours::addBurriedEdge(): forward_edge = " << forward_edge_id << ", backward_edge = " << backward_edge_id);
		while (points_by_edge.size() <= forward_edge_id || points_by_edge.size() <= backward_edge_id) points_by_edge.push_back(EdgePoints()); //TODO: you can initialize this outright in constructor!
		points_by_edge[forward_edge_id].edge_type = EdgePoints::EDGE_TYPE_BURRIED;
		points_by_edge[backward_edge_id].edge_type = EdgePoints::EDGE_TYPE_BURRIED;
	}
	/**
	 * marks given edge as Free edge (outside the sphere)
	 *
	 * Currently we do not use this in any way, but we remember this as a side effect of computation
	 */
	void addFreeEdge(size_t forward_edge_id, size_t backward_edge_id){																							DEB3("SASAContours::addFreeVertex(): forward_edge = " << forward_edge_id << ", backward_edge = " << backward_edge_id);
		while (points_by_edge.size() <= forward_edge_id || points_by_edge.size() <= backward_edge_id) points_by_edge.push_back(EdgePoints()); //TODO: you can initialize this outright in constructor!
		points_by_edge[forward_edge_id].edge_type = EdgePoints::EDGE_TYPE_FREE;
		points_by_edge[backward_edge_id].edge_type = EdgePoints::EDGE_TYPE_FREE;
	}
	/**
	 * adds vertex2 and side_id center as on_polygon VertexLink, only if it wasn't set earlier
	 */
	bool conditionalConnectPolygon(int vertex_1, int vertex_2, int side_id){																				DEB3("SASAContours::conditionalConnectPolygon(): v_1 = " << vertex_1 << ", v_2 = " << vertex_2 << ", c = " << side_id);
		if (vertices[vertex_1].next_on_polygon.vertex < 0){
			vertices[vertex_1].next_on_polygon.center = side_id;
			vertices[vertex_1].next_on_polygon.vertex = vertex_2;
			return true;
		}else{				DEB3("SASAContours::conditionalConnectPolygon(): Already connected...");
			return false;
		}
	}
	/**
	 * adds vertex2 and side_id center as on_border VertexLink, only if it wasn't set earlier
	 */
	bool conditionalConnectBorder(int vertex_1, int vertex_2, int side_id){																					DEB3("SASAContours::conditionalConnectBorder(): v_1 = " << vertex_1 << ", v_2 = " << vertex_2 << ", c = " << side_id);
		if (vertices[vertex_1].next_on_border.vertex < 0){
			vertices[vertex_1].next_on_border.center = side_id;
			vertices[vertex_1].next_on_border.vertex = vertex_2;
			return true;
		}else{				DEB3("SASAContours::conditionalConnectBorder(): Already connected...");
			return false;
		}
	}
	/**
	 * iterator thru all border contours
	 */
	SASABorderIterator borderIterator(){
		return SASABorderIterator(this);
	}
	/**
	 * iterator thru all inside polygon contours
	 */
	SASAPolygonIterator polygonIterator(){
		return SASAPolygonIterator(this);
	}

	/*
	Iterator begin(){

	}

	Iterator end(){

	}
	*/

	~SASAContours(){ DEB3("SASAContours::__destruct__():");
	}
};

/**
 * w tej klasie umieścic wszystko co bedzie umiała robić jedna komórka
 * SASA - liczyc powierzchnie, objętości, wypisywac sie na cout, itp
 */
template<typename CellSpec>
class SASAShapeCell : public DMGAObject{
public:
	std::string show(){ return "dmga::shape::SASAShapeCell"; }

	typedef CellSpec CellType;

	CellType* cell;
	SASAContours* contours;

	SASAShapeCell(CellType* cell = 0, SASAContours* contours = 0) :
			cell(cell), contours(contours){
	}

	SASAShapeCell(CellType& cell, SASAContours& contours) :
			cell(&cell), contours(&contours){
	}
	/**
	 * serializuj do JSONa
	 */
	std::string serialize(){
		std::ostringstream out;
		return out.str();
	}
	double area(){
		double r = this->contours->sphere[3];
		return 4.0 * r*r * M_PI - this->contours->excluded_area;
	}
	double volume(){
		throw dmga::exceptions::NotImplementedYet("SASAShapeCell::getVolume()");
		return 0.0;
	}
	CellType* getCell(){
		return cell;
	}
	SASAContours* getContours(){
		return contours;
	}
	friend std::ostream& operator<<(std::ostream& out, const SASAShapeCell& shape){
		out << shape.serialize();
		return out;
	}

	~SASAShapeCell(){ DEB3("SASAShapeCell::__destruct__():");
	}
};

//TODO: dodac sprawdzanie czy DiagramSpec działa na WeightedParticle...
template<typename DiagramSpec>
class VoroppSASAShape : public GenericShape{
public:
	typedef DiagramSpec DiagramType;
	typedef typename DiagramType::ContainerType ContainerType;
	typedef typename DiagramType::GeometryType GeometryType;
	typedef typename DiagramType::ParticleType ParticleType;
	typedef typename DiagramType::KeyType KeyType;
	typedef typename DiagramType::CellType CellType;
	typedef typename DiagramType::ResultSetType ResultSetType;
	typedef typename DiagramType::SubsetType SubsetType;
	typedef typename DiagramType::CacheType CacheType;
	typedef SASAShapeCell<CellType> CellShapeType;

	DiagramType* diagram;
	ContainerType* container;
	std::vector<SASAContours*> contours;
	bool internal_diagram;
	double solvent_radius; //this radius will be added to the

	virtual ~VoroppSASAShape(){								DEB3("VoroppSASAShape::__destruct__():");
		if (internal_diagram){
			DEB3("VoroppSASAShape::__destruct__(): freeing internal diagram object");
			dmga::utils::safeDelete(this->diagram);
		}
	}

	VoroppSASAShape(DiagramType& diagram):
				diagram(&diagram),
				container(&(diagram.container)),
				contours(diagram.container.size()),
				internal_diagram(false),
				solvent_radius(0.0){
		recompute();
	}

	//TODO: zmianiemy container, to co jest w srodku, czy tak dobrze?
	VoroppSASAShape(ContainerType& container, double solvent_radius = 0.0):
				diagram(new DiagramType(container)),
				container(&container),
				contours(container.size()),
				internal_diagram(true),
				solvent_radius(solvent_radius){
		for (int i = 0; i < container.size(); ++i){
			//TODO: przemyśleć, czy tak jest oki...
			container[i]->getRaw()[4] += solvent_radius;
		}
		recompute();
	}

	//TODO: zmianiemy container, to co jest w srodku, czy tak dobrze?
	VoroppSASAShape(ContainerType& container, SubsetType& subset, double solvent_radius = 0.0):
					solvent_radius(solvent_radius),
					internal_diagram(true){
		for (int i = 0; i < container.size(); ++i){
			//TODO: przemyśleć, czy tak jest oki...
			container[i]->getRaw()[4] += solvent_radius;
		}
		diagram = new DiagramType(container, subset);
		this->container = &container;
		recompute();
	}

	SASAContours* computeOne(int i){
		CellType* cell = diagram->getCell(i); 				DEB3("(SASAShape) compute one cell " << cell->number);
		voro::CellIterator cell_it = cell->getCellIterator();
		int current_side = 0; //current side, we do standard iteration, so we get the right ordering

		//get the data of the ball, it can be done once, because this will be the same for all computations
		double* c = container->getRawCoords(i); //center of the ball
		double R = c[3]; //radius of the ball
		double RSQ = R*R; //R^2

		if (cell->isEmpty()){
			//nic nie trzeba robić - całą powierzchnię kuli wciągnęły inne kule?... TODO: przemyśleć
			//zapisz, że kula pusta, wszystko puste...
			return 0;
		}

		//stwórz miejsce do trzymania danych o tym co mamy SASAContours
		SASAContours* cell_mesh = new SASAContours();
		while (!cell_it.isFinished()){						DEB3("(new side)" << current_side);
			//start iteration on the current side
			int v, j, u, k, w; //current vertex, its forward edge on this side, next vertex, edge dual to j (from u to v)...
			// we need triangle w -> v -> u -> w to compute side normal and point on the side and on the sphere
			v = cell_it.current().v; // the v-th vertex
			u = cell_it.inv().v; // the end of the j-th edge in v is vertex u
			w = cell_it.prev().v;

			//get the neigbour id
			int neighbour_id = cell->getNeighbourId(cell_it); 	DEB3("VoroppSASAShape::computeOne(): This side neighbour is " << neighbour_id);
			//double *nc = this->container->getRawCoords(neighbour_id); 	DEB3("VoroppSASAShape::computeOne(): Neighbour pos is " << nc[0] << "," << nc[1] << "," << nc[2] << ", r = " << nc[3]);

			double on_sphere[3]; //here we have point on sphere on line connecting this sphere center and neighbour sphere center
			double on_plane[3]; //here we have point on power plane defined by this sphere and its neighbour
			// this is the first time we visit a side, so we can add center (one per each side)
			// double* s = container->getRawCoords(neighbour_id);	DEB3("VoroppSASAShape::computeOne(): Neighbour coords are " << s[0] << ", " << s[1] << ", " << s[2]);
			// double n[3];
			// double pow_d = algebra::distance::from_power_plane(c, s, n);		DEB3("VoroppSASAShape::computeOne(): Distance to power plane is " << pow_d);
			// DEB3("VoroppSASAShape::computeOne(): Normal is " << n[0] << ", " << n[1] << ", " << n[2]);

			// to overcome problem with sides which do not have neighbours we use triangle w -> v -> u -> w (the same as in AlphaShape)
			double side_normal[4];
			// we will use normal to compute the point on the plane defined by side where the alpha-ball first hit the plane,
			// it will be used to test if alpha-threshold need to pe repaired
			// we get a side_normal as a side effect of the computation of the euclidean distance to the plane
			double to_side_eucl_dist = algebra::distance::euclid3d(c,
															cell_it.vertexRaw(v),
															cell_it.vertexRaw(u),
															cell_it.vertexRaw(w),
															side_normal); //here the sequence of v, u, w is important, the same as in computing side_normal
			//Add center of this side to the mesh
			// the point coordinates are: c + n * R
			on_sphere[0] = c[0] - R * side_normal[0];
			on_sphere[1] = c[1] - R * side_normal[1];
			on_sphere[2] = c[2] - R * side_normal[2];				DEB3("VoroppSASAShape::computeOne(): On sphere point is " << on_sphere[0] << ", " << on_sphere[1] << ", " << on_sphere[2]);
			on_plane[0] = c[0] - to_side_eucl_dist * side_normal[0];
			on_plane[1] = c[1] - to_side_eucl_dist * side_normal[1];
			on_plane[2] = c[2] - to_side_eucl_dist * side_normal[2];		DEB3("VoroppSASAShape::computeOne(): On plane point is " << on_plane[0] << ", " << on_plane[1] << ", " << on_plane[2]);
			//add the center to the collection of centers
			cell_mesh->addCenter(current_side, neighbour_id, on_sphere[0], on_sphere[1], on_sphere[2], on_plane[0], on_plane[1], on_plane[2]);		DEB3("VoroppSASAShape::computeOne(): center added...");

			if (to_side_eucl_dist >= R){
				DEB3("VoroppSASAShape::computeOne(): Side not touched by ball. Going to next side. ");
				//we have situation when sphere doesn't touch the dividing plane between particles - we can skip whole computation
				//TODO: add free edges... ?
				while (!cell_it.isMarked()){
					cell_it.mark();
					cell_it.forward();
				}
				current_side++; //increase current side id
				cell_it.jump(); //go to next side
				continue; //go to the beginning of the loop
			}

			DEB3("VoroppSASAShape::computeOne(): Side touched by ball. Computing intersections");
			int first_near = -1;
			int first_far = -1;
			int current_near = -1;
			int current_far = -1;
			while (!cell_it.isMarked()){
				v = cell_it.current().v; // the v-th vertex
				j = cell_it.current().j; // j-th edge in vertex v
				u = cell_it.inv().v; // the end of the j-th edge in v is vertex u
				k = cell_it.inv().j; // k - the next edge on the side (from u)

				// we have data for forward and backward edge - we can store it and use it later, so we skip half of the edges...
				int edge_id_forward = cell_it.edge(v, j);		DEB3("VoroppSASAShape::computeOne(): forward edge id is " << edge_id_forward);
				int edge_id_backward = cell_it.edge(u, k);		DEB3("VoroppSASAShape::computeOne(): backward edge id is " << edge_id_backward);
				//if not computed then compute...
				int near_vertex_id = -1;
				int far_vertex_id = -1;
				if (cell_mesh->getEdgePoints(edge_id_forward).edge_type == SASAContours::EdgePoints::EDGE_TYPE_UNKNOWN){		DEB3("VoroppSASAShape::computeOne(): Edge not computed...");
					// compute itersetion of the current edge with the ball of a given radii (ball_data, ball_r)
					double* p1 = cell_it.vertexRaw(v); //first point on the edge
					double* p2 = cell_it.vertexRaw(u); //second point on the edge
														DEB3("VoroppSASAShape::computeOne(): p1 is " << p1[0] << "," << p1[1] << "," << p1[2]);
														DEB3("VoroppSASAShape::computeOne(): p2 is " << p2[0] << "," << p2[1] << "," << p2[2]);
														DEB3("VoroppSASAShape::computeOne(): c " << c[0] << "," << c[1] << "," << c[2]);
														DEB3("VoroppSASAShape::computeOne(): R^2 " << RSQ);
					//let p1 = (x1, y1, z1), p2 = (x2,y2,z2), c = (xc,yc,zc), R = raduis of the ball
					//then each point on the line between p1 and p2 is given by: p(t) = t*p2 + (1-t) * p1 for t \in [0,1]
					//we need to solve A*t^2 + B*t + C for t \in [0,1]
					//where A = (x2-x1)^2 + (y2-y1)^2 + (z2-z1)^2
					//      B = (x1-xc)*(x2-x1) + (y1-yc)*(y2-y1) + (z1-zc)*(z2-z1)
					//      C = (x1-xc)^2 + (y1-yc)^2 + (z1-zc)^2 - R^2
					// this comes from the equation (p(t)-c).(p(t)-c) = R^2 for point at R from center of ball
					// REMARK: this may be boosted, or we simply may assume that compiler do their best do optimize it!
					double A = (p2[0]-p1[0])*(p2[0]-p1[0]) + (p2[1]-p1[1])*(p2[1]-p1[1]) + (p2[2]-p1[2])*(p2[2]-p1[2]);
					double B = 2.0 * ((p1[0]-c[0])*(p2[0]-p1[0]) + (p1[1]-c[1])*(p2[1]-p1[1]) + (p1[2]-c[2])*(p2[2]-p1[2]));
					double C = (p1[0]-c[0])*(p1[0]-c[0]) + (p1[1]-c[1])*(p1[1]-c[1]) + (p1[2]-c[2])*(p1[2]-c[2]) - RSQ;			DEB3("VoroppSASAShape::computeOne(): Points on edge equation = " << A << " x^2 + " << B << " x + " << C);
					double delta = B*B - 4.0 * A * C;																			DEB3("VoroppSASAShape::computeOne(): Delta = " << delta);
					if (delta > 0.0){ //we do not want to waste time for delta == 0, when ball is only 'touching' the edge - in that case we can treat this as whole circle on the side
						//moreover we also merge points if they are too close to each other, because they may then cause problems when computing tangent vectors etc.
						//if A too small then errors - if p1 and p2 too close to each other...
						//TODO: if A or delta, or dist(t1, t2) small then HEURISTIC: test if p1, p2 inside ball then threat p1, p2 as inside, otherwise threat this as ball 'touching' edge - should work?
						double t1, t2; //we would like to have t1 <= t2
						t1 = (-B - sqrt(delta)) / (2.0 * A); //t1 is closer to v
						t2 = (-B + sqrt(delta)) / (2.0 * A); //t2 is closer to u
																					DEB3("VoroppSASAShape::computeOne(): t1 = " << t1 << ", t2 = " << t2);
						//compute points pc1 and pc2 - on the line
						if (0.0 < t1 && t1 < 1.0){ // t1 is smaller than t1 -> so we call the point p(t = t1) the near point
							//addYYYVertex takes care of remembering edge type in this case
							near_vertex_id = cell_mesh->addNearVertex(edge_id_forward, edge_id_backward, t1 * p2[0] + (1-t1) * p1[0], t1 * p2[1] + (1-t1) * p1[1], t1 * p2[2] + (1-t1) * p1[2]);
						}
						if (0.0 < t2 && t2 < 1.0){ // t2 is greater than t1 -> so we call the point p(t = t2) the far point
							//addYYYVertex takes care of remembering edge type in this case
							far_vertex_id = cell_mesh->addFarVertex(edge_id_forward, edge_id_backward, t2 * p2[0] + (1-t2) * p1[0], t2 * p2[1] + (1-t2) * p1[1], t2 * p2[2] + (1-t2) * p1[2]);
						}
						//if after two above if-s we do not have edge type, then edge is burried
						if (cell_mesh->getEdgePoints(edge_id_forward).edge_type == SASAContours::EdgePoints::EDGE_TYPE_UNKNOWN){
							cell_mesh->addBurriedEdge(edge_id_forward, edge_id_backward);		DEB3("VoroppSASAShape::computeOne(): Burried edge...");
						}
					}else{ // sphere doesn't cross the edge
						cell_mesh->addFreeEdge(edge_id_forward, edge_id_backward);		DEB3("VoroppSASAShape::computeOne(): Free edge...");
					}
				}else{		DEB3("VoroppSASAShape::computeOne(): Already computed...");
					//we have computed vertices - we need to retrive them and use them
					auto edge_points = cell_mesh->getEdgePoints(edge_id_forward);
					if (edge_points.edge_type == SASAContours::EdgePoints::EDGE_TYPE_DOUBLE ||
							edge_points.edge_type == SASAContours::EdgePoints::EDGE_TYPE_NEAR){
						near_vertex_id = edge_points.near_point;
					}
					if (edge_points.edge_type == SASAContours::EdgePoints::EDGE_TYPE_DOUBLE ||
							edge_points.edge_type == SASAContours::EdgePoints::EDGE_TYPE_FAR){
						far_vertex_id = edge_points.far_point;
					}
				}

				//update near vertex links if it exists
				if (near_vertex_id >= 0){ //we have new near_vertex
					if (first_near < 0) first_near = near_vertex_id; //remember first near vertex...
					current_near = near_vertex_id;
					if (current_far >= 0){
						//connect both elements - we know that this is far->near relation = boundary connection
						cell_mesh->conditionalConnectBorder(current_far, current_near, current_side);
					}
				}

				if (far_vertex_id >= 0){ //we have new far_vertex
					//update far vertex links if it exists
					if (first_far < 0) first_far = far_vertex_id; //remember the first far if it leaves queue
					current_far = far_vertex_id;
					if (current_near >= 0){
						//connect both elements - we know that this is near->far relation = polygon connection
						cell_mesh->conditionalConnectPolygon(current_near, current_far, current_side);
					}
				}
				cell_it.mark();
				cell_it.forward();
			}
			//TODO: przemyslec zamykanie loops - poki co ok (conditionalConnect), ale moze da sie lepiej
			if (first_near >= 0 && current_far >= 0){
				//close loop in far->near relation = boundary connection
				cell_mesh->conditionalConnectBorder(current_far, first_near, current_side);
			}
			if (first_far >= 0 && current_near >= 0){
				//close loop in near->far relation = polygon connection
				cell_mesh->conditionalConnectPolygon(current_near, first_far, current_side);
			}
			if (first_near < 0 && first_far < 0){ //it suffices to check only one, but for clarity I check both - maybe in release version remove
				//we can have only crossing on the side
				double *pp = cell_mesh->centers[current_side].point_on_plane.getRaw();
				double *ps = cell_mesh->centers[current_side].point_on_sphere.getRaw();
				double dist_pp = algebra::distance::euclid3d_pow_2(pp, c);
				double dist_ps = algebra::distance::euclid3d_pow_2(ps, c);
				if (dist_pp < dist_ps){
					//move around the side and find pair of edges around the closest vertex
					//one of them is the closest edge to this point (pp - on plane)
					//because only two points may be the closest points if the point is on the outside
					//if the point is inside we do not care
					int prev_v = cell_it.prev().v;
					int v = cell_it.current().v;
					int first_v = v;
					int next_v = cell_it.next().v;
					double* v_coords = cell_it.vertexRaw();
					double min_dist_sq = algebra::distance::euclid3d_pow_2(pp, v_coords);
					int min_prev_v = prev_v;
					int min_v = v;
					int min_next_v = next_v;
					cell_it.forward();
					while (cell_it.current().v != first_v){
						prev_v = v;
						v = next_v;
						next_v = cell_it.next().v;
						v_coords = cell_it.vertexRaw();
						double d = algebra::distance::euclid3d_pow_2(pp, v_coords);
						if (d < min_dist_sq){
							min_dist_sq = d;
							min_prev_v = prev_v;
							min_v = v;
							min_next_v = next_v;
						}
						cell_it.forward();
					}
					// ok we have the edges: prev_v -> v -> next_v
					// do standard test if point is inside
					v_coords = cell_it.vertexRaw(min_v);
					double* prev_v_coords = cell_it.vertexRaw(min_prev_v);
					double* next_v_coords = cell_it.vertexRaw(min_next_v);
					//those are vectors prev->v and v->next, we will use them to compute normals to the plane given
					// by them and side_normal - those will poin outwards on the plane of the cutting plane - we will
					// use it to test if point is inside or outside
					double n_prev_v_v[3] = {v_coords[0] - prev_v_coords[0], v_coords[1] - prev_v_coords[1], v_coords[2] - prev_v_coords[2]};
					double n_v_next_v[3] = {next_v_coords[0] - v_coords[0], next_v_coords[1] - v_coords[1], next_v_coords[2] - v_coords[2]};
					bool test_outside = false;
					double test_normal[4];
					//here is the vector from v to pp that we test if is outside the plane
					double test_vector[3] = {pp[0] - v_coords[0], pp[1] - v_coords[1], pp[2] - v_coords[2]};

					algebra::normal3d_raw(side_normal, n_prev_v_v, test_normal);
					if (algebra::dot3d(test_normal, test_vector) < 0) test_outside = true;
					algebra::normal3d_raw(side_normal, n_v_next_v, test_normal);
					if (algebra::dot3d(test_normal, test_vector) < 0) test_outside = true;

					//bool test_outside = false;
					if (!test_outside){
						//we have intersection because point on side is closer than a point on sphere...
						// V = (\pi * h^2 / 3.0) * (3.0 * r - h)
						// A = 2 * \pi * r * h
						//where
						//   h = dist(ps, pp) \in [0, 2r] (h \in [r, 2r]) is if the hemisphere is greater than half the ball
						//   r = R (radius of the ball)
						double h = algebra::distance::euclid3d(ps, pp);
	//					cell_mesh->centers[current_side].cap_intersection_area = 2.0 * M_PI * R * h;
	//					cell_mesh->centers[current_side].cap_intersection_volume = M_PI * h * h * (3.0 * R - h) / 3.0;
	//					cell_mesh->centers[current_side].cap_radius_sq = 2 * h * R - h * h;
						double cap_radius_sq = 2 * h * R - h * h;
						//we need to find a point on side that lies on intersection with the sphere. Any point will be good.
						//so we take arbitrary vertex on the side and place some vertex at distance sqrt(cap_radius_sq) on the line from the center to this point
						double *some_point = cell_it.vertexRaw();
						double dist_sq = algebra::distance::euclid3d_pow_2(pp, some_point);
						double scale = sqrt(cap_radius_sq / dist_sq); //to use only one sqrt
						//behold (the evil is strong!): here we know that center exists (centers are always added) so we can trust
						//addOnSideVertex to set appropriate centers in connectBorder and connectPolygon (see addOnSideVertex)
						cell_mesh->addOnSideVertex(	current_side,
													pp[0] + (some_point[0] - pp[0]) * scale,
													pp[1] + (some_point[1] - pp[1]) * scale,
													pp[2] + (some_point[2] - pp[2]) * scale);
						DEB3("VoroppSASAShape::computeOne(): cap_intersection_area at cell " << i << " side " << current_side << " is " << (2.0 * M_PI * R * h));
					}else{
						DEB3("VoroppSASAShape::computeOne(): cap is outside the side " << i << " side " << current_side);
					}
				}
			}

			current_side++; //increase current side id
			cell_it.jump();
		}

		//OK we have the mesh - SASAContours - we need to compute some more things to get the area and volume!
		DEB3("VoroppSASAShape::computeOne(): Computing arcs finished. Computing area & volume.");
		double sphere_angle = 4.0 * M_PI; //this is sphere_area / R^2
		double sasa_angle = 0;

		//przejdź po polygonach i policz kąty
		auto pit = cell_mesh->polygonIterator();
		double test_normal[4]; //here we store normal for test if the angle is greater than 180*
		double angle_normal[4]; //here we store normal for test if the angle is greater than 180*
		double prev_normal[4]; //here we store normal to the plane given by (origin, prev_center, current)
		double curr_normal[4]; //here we store normal to the plane defined by (origin, centrum, current)
		if (pit.ready()){
			while (true){ //for all polygon contours
				//remark on orientation of computing normals - normal3d_raw should be feed coordinates
				// the same order each time, ie following rule of thumb (left/right hand) - do the picture if you need
				if (pit.current == pit.next){ //we have cap. Caps doesn't have polygons inside, so we skip this iteration and jump to the next vertex
					pit.mark(); //to make things clearer, but we can skip this.
					DEB3("VoroppSASAShape::computeOne(): skipping polygon at cap at cell " << i << " side " << pit.center);
					if (pit.jump() == -1) break; //we jump and we break when there is no other polygons
				}
				int prev_center = pit.center;
				int first = pit.current;
				pit.mark(); //mark it
				pit.forward(); //and go further
				// remark: origin is in c variable
				double *prev_center_raw = cell_mesh->centers[prev_center].point_on_sphere.getRaw(); // it should be ok to use anyone point
				double *current_raw = cell_mesh->vertices[pit.current].getRaw();
				double *center_raw, *next_raw; //we will use this later
				dmga::algebra::normal3d_raw(c, current_raw, prev_center_raw, prev_normal);
				//now we start iteration
				int tes = 2;
				int n = 0;
				double angle_sum = 0.0;
				while (tes > 0){ //we know that this will eventually happen
					DEB3("VoroppSASAShape::computeOne(): polygon arcs: " << pit.current << ", " << pit.center << ", " << pit.next);
					DEB3("VoroppSASAShape::computeOne(): polygon arcs on side with neighbour: " << cell_mesh->centers[pit.center].neighbour_id);
					//we have those elements:
					// prev_center
					// pit.current (vertex)
					// pit.center
					// pit.next (vertex)
					//we need to compute two angles
					// alpha - prev_center, pit.current, pit.center
					// beta - pit.current, pit.center, pit.next
					// to compute these angles we need also origin of the sphere
					// because we need to compute angles between two planes given
					// for alpha we have planes:
					// 		(origin, prev_center, pit.current)
					//	and
					// 		(origin, pit.current, pit.center)
					// remark: origin is in c variable
					current_raw = cell_mesh->vertices[pit.current].getRaw();
					center_raw = cell_mesh->centers[pit.center].point_on_sphere.getRaw(); //we assuming that R > 0 so the point os sphere != origin (c) - it is ok? It should be that if R == 0 then we should not be here to compute this
					next_raw = cell_mesh->vertices[pit.next].getRaw();
					dmga::algebra::normal3d_raw(c, current_raw, center_raw, curr_normal); //current normal
//					std::cout << "[";
//					std::cout << "(" << prev_normal[0] + current_raw[0] << "," << prev_normal[1]+current_raw[1] << "," << prev_normal[2]+current_raw[2] << "),";
//					std::cout << "(" << current_raw[0] << "," << current_raw[1] << "," << current_raw[2] << "),";
//					std::cout << "(" << curr_normal[0] + current_raw[0] << "," << curr_normal[1]+current_raw[1] << "," << curr_normal[2]+current_raw[2] << ")";
//					std::cout << "]\n";
					double alpha = acos(dmga::algebra::dot3d(curr_normal, prev_normal));
					//we also need to check if the angle is greater than pi ...
					test_normal[0] = current_raw[0] - c[0]; //we need to justify that
					test_normal[1] = current_raw[1] - c[1]; //in beta we will have exact matching along normals
					test_normal[2] = current_raw[2] - c[2]; //here it may not be the case, it should be close, but still...
					dmga::algebra::normal3d_raw(curr_normal, prev_normal, angle_normal);
					if (dmga::algebra::dot3d(test_normal, angle_normal) < 0){ //kąt rozwarty - nalezy dodac PI do wyniku
						alpha += M_PI;
					}

//					std::cout << "alpha: " << alpha << " (" << (alpha * 180.0 / M_PI) <<  ")" << "\n";

					//ok we have computed angle alpha, prev_normal is not needed anymore, we may use it, and
					//automaticaly we get prev_normal for next iteration (for current = next) - see it clear!
					dmga::algebra::normal3d_raw(c, next_raw, center_raw, prev_normal); //current normal
//					std::cout << "dot: " << dmga::algebra::dot3d(curr_normal, prev_normal) << "\n";
					double beta = acos(dmga::algebra::dot3d(curr_normal, prev_normal));
//					std::cout << "[";
//					std::cout << "(" << prev_normal[0] + center_raw[0] << "," << prev_normal[1] + center_raw[1] << "," << prev_normal[2] + center_raw[2] << "),";
//					std::cout << "(" << center_raw[0] << "," << center_raw[1] << "," << center_raw[2] << "),";
//					std::cout << "(" << curr_normal[0] + center_raw[0] << "," << curr_normal[1] + center_raw[1] << "," << curr_normal[2] + center_raw[2] << ")";
//					std::cout << "]\n";
					//we also need to check if the angle is greater than pi
					test_normal[0] = center_raw[0] - c[0];
					test_normal[1] = center_raw[1] - c[1];
					test_normal[2] = center_raw[2] - c[2];
					dmga::algebra::normal3d_raw(prev_normal, curr_normal, angle_normal); //order makes difference - it is used that dot3 < 0 for angles > pi
					if (dmga::algebra::dot3d(test_normal, angle_normal) < 0){ //kąt rozwarty - nalezy dodac PI do wyniku
						beta += M_PI;
					}

					DEB3("VoroppSASAShape::computeOne(): beta: " << beta << " (" << (beta * 180.0 / M_PI) <<  ")" << "\n");

					n += 2;
					angle_sum += alpha + beta;

					//ok we computed angles, we have prev_normal, we can go further
					prev_center = pit.center; //not necessary anymore to remember this?
					pit.mark();
					pit.forward();
					if (tes < 2) tes--; if (pit.current == first) tes--; //this is to ensure that we will add last angles (after returning to starting point on this contour
				}
				angle_sum -= M_PI * (n-2);
				angle_sum *= RSQ;
				DEB3("VoroppSASAShape::computeOne(): This side angles: " << angle_sum << "\n");
				while (angle_sum > sphere_angle) angle_sum -= sphere_angle;
				DEB3("VoroppSASAShape::computeOne(): This side angles (adjusted): " << angle_sum << "\n");
				sasa_angle += angle_sum;
				//now we finished with this part, we can go to the next part
				if (pit.jump() == -1) break; //we break when there is no other polygons
				DEB3("VoroppSASAShape::computeOne(): next contour...\n");
			}
		}
		while (sasa_angle > sphere_angle) sasa_angle -= sphere_angle; //adjust it to not exceed allowed max angle

		double sasa_excluded_area = sasa_angle * RSQ; //this is area excluded for now from the sphere (sum of areas of polygons)

		//przejdź po borderze i policz katy
		//this is simpler, as we only have a collection of arcs, and the angle in each arc is given by all points in bit iterator
		auto bit = cell_mesh->borderIterator();
		double *current_raw;
		double *sphere_center_raw;
		double *center_raw;
		double *next_raw;
		double v1[3], v2[3]; //for tests
		double gamma; //here we will store the angle of the cap part in radians gamma \in 0 .. 2*\pi
		if (bit.ready()) while (true){ //for all polygon contours
			while (!bit.isMarked()){ //we know that this will eventually happen
				bit.mark(); //mark it

				if (bit.current == bit.next){
					gamma = 2 * M_PI;
					//we need to heve those points to compute h - height of the cap.
					sphere_center_raw = cell_mesh->centers[bit.center].point_on_sphere.getRaw();
					center_raw = cell_mesh->centers[bit.center].point_on_plane.getRaw();
					DEB3("VoroppSASAShape::computeOne(): we have cap here at cell " << i << " side " << bit.center);
					DEB3("VoroppSASAShape::computeOne(): gamma: " << gamma << " (" << (gamma * 180.0 / M_PI) <<  ")" << "\n");
				}else{
					current_raw = cell_mesh->vertices[bit.current].getRaw();
					sphere_center_raw = cell_mesh->centers[bit.center].point_on_sphere.getRaw();
					center_raw = cell_mesh->centers[bit.center].point_on_plane.getRaw();
					next_raw = cell_mesh->vertices[bit.next].getRaw();
					v1[0] = current_raw[0] - center_raw[0];
					v1[1] = current_raw[1] - center_raw[1];
					v1[2] = current_raw[2] - center_raw[2];
					v2[0] = next_raw[0] - center_raw[0];
					v2[1] = next_raw[1] - center_raw[1];
					v2[2] = next_raw[2] - center_raw[2];
					dmga::algebra::normalize3d_raw(v1);
					dmga::algebra::normalize3d_raw(v2);
					gamma = acos(dmga::algebra::dot3d(v1, v2));
					DEB3("VoroppSASAShape::computeOne(): current_raw = " << Point3D(current_raw) << ", center_raw = " << Point3D(center_raw) << ", next_raw = " << Point3D(next_raw) << ";");
					DEB3("VoroppSASAShape::computeOne(): v1 = " << Point3D(&v1[0]) << ", v2 = " << Point3D(&v2[0]));
					DEB3("VoroppSASAShape::computeOne(): gamma: " << gamma << " (" << (gamma * 180.0 / M_PI) <<  ")" << "\n");
				}

				//compute this cap area and scale it with gamma / (2 * PI)
				double h = algebra::distance::euclid3d(sphere_center_raw, center_raw);
				double arc_area = gamma * R * h; // cap = 2.0 * M_PI * R * h, arc = 2.0 * M_PI * R * h * gamma / (2 * M_PI) = gamma * R * h;
				if (gamma == 2.0 * M_PI) DEB3("VoroppSASAShape::computeOne(): we have cap here at cell " << i << " side " << bit.center << " with area " << gamma * R * h);
				//double arc_volume = gamma * ((3.0 * R - h) * h*h) / 6.0; // cap = (M_PI / 3.0) * (3r-h)*h^2, just multiply by gamma / (2 * M_PI)
				//arc_volume not used yet - we must re-think how to compute volume...
				DEB3("VoroppSASAShape::computeOne(): arc_area: " << arc_area << "\n");
				sasa_excluded_area += arc_area;
				//sasa_volume -= arc_volume;

				bit.forward(); //and go further
			}
			if (bit.jump() == -1) break;
			DEB3("VoroppSASAShape::computeOne(): next contour...\n");
		}
		//THIS IS OLD PART WHERE WE DIDN'T STORE FULL CAPS (i.e. when no intersection on side edge is detected). Now we hold those caps as Arcs with some artificial vertex on side.
//		int centers_size = cell_mesh->centers.size();
//		for (int j = 0; j < centers_size; ++j){
//			if (cell_mesh->centers[j].cap_intersection_area > 0){
//				// we have cap_intersection_area only if sphere clearly intersected this side as cap
//				sasa_excluded_area += cell_mesh->centers[j].cap_intersection_area;
//			}
//		}

		if (sasa_excluded_area > sphere_angle * RSQ) sasa_excluded_area = sphere_angle * RSQ; //should not happen but just in case...
		cell_mesh->excluded_area = sasa_excluded_area;
		cell_mesh->excluded_volume = 0.0; //TODO: obliczyc i zapisac
		cell_mesh->sphere[0] = c[0];
		cell_mesh->sphere[1] = c[1];
		cell_mesh->sphere[2] = c[2];
		cell_mesh->sphere[3] = c[3]; //you can use memcpy or something, but for now it was faster to code ;)
		return cell_mesh;
	}

	void recompute(){										DEB3("VoroppSASAShape::recompute(): ");
		auto cache_it = diagram->cache.begin();
		auto cells_end = diagram->cache.end();
		for (; cache_it != cells_end; ++cache_it){
			if (*cache_it){
				computeOne((*cache_it)->number);
			}
		}
		DEB3("VoroppSASAShape::recompute(): finished");
	}

	double area(){
		std::vector<SASAContours*>::iterator it = contours.begin();
		std::vector<SASAContours*>::iterator end = contours.end();
		double sasa_area = 0;
		for (; it != end; ++it){
			if (*it){
				double r = (*it)->sphere[3];
				sasa_area += 4.0 * r*r * M_PI - (*it)->excluded_area;
			}
		}
		return sasa_area;
	}

	SASAContours* getContours(int i){					DEB3("VoroppSASAShape::getContours(): i = " << i);
		if (!contours[i]){
			//TODO: test if cell computed
			contours[i] = this->computeOne(i);
		}
		return contours[i];
	}

	AlphaComplex* getContours(const CellType& cell){	DEB3("VoroppSASAShape::getContours(): cell.number = " << cell.number);
		return getContours(cell.number);
	}

	CellShapeType* getCellShape(int i){						DEB3("VoroppSASAShape::getCellShape(): i = " << i);
		//CellType* cell = diagram->cache[i];
		CellType* cell = diagram->getCell(i);
		SASAContours* contours = this->getContours(i);
		return new CellShapeType(cell, contours);
	}

	CellShapeType* getCellShape(CellType& cell){			DEB3("VoroppSASAShape::getCellShape(): cell.number = " << cell.number);
		return new CellShapeType(cell, this->getContours(cell.number));
	}


};

} //namespace shape

} // namespace dmga


#endif /* SASA_SHAPE_H_ */
