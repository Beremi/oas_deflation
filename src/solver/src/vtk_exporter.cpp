#include "vtk_exporter.h"

//////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////
// VTK EXPORTERS

void VTKExporter :: giveFileName(unsigned step, char *buffer) const {
  sprintf(buffer, "%s_%05d.vtu", filename.c_str(), step);
}

//////////////////////////////////////////////////////////
void VTKExporter :: readFromLine(istringstream &iss, unsigned dimension){
  string dummy;  // to skip the specifier
  iss >> filename >> dummy >> time_each;
  unsigned num;
  iss >> num;
  codes.resize(num);
  for ( unsigned i = 0; i < num; i++ ) {
      iss >> codes [ i ];
  }
  time_last = 0.;
}
//////////////////////////////////////////////////////////
bool VTKExporter :: doExportNow(const double &time) {
  if (time < time_last + time_each) {
    return false;
  } else {
    time_last = time;
    return true;
  }
}



//////////////////////////////////////////////////////////
// ELEMENTS TO VTU FILE
//////////////////////////////////////////////////////////
void VTKElementExporter :: exportData(unsigned step, const Vector &DoFs, const Vector &reactions) const{
  // Export of elements into vtu xml file format (vtu = vtk for unstructured grid)
  // NOTE this is messy construction of xml file, will be remade using some of xml libraries for cpp
  char buffer [ 100 ];
  // Point P;
  // Element *ee;

  vector< int > points_id;
  vector< vector< int > > all_points_id;
  vector<int> cell_types;
  vector<int> offsets;
  vector<Point> displ;
  vector<double> damage;  // test version, this and more will be specified on the exporter input
  int offset = 0;
  for (auto const &el : *elems){
    for (auto const &n : el->giveNodes()){
      auto res = std::find(begin(*nodes), end(*nodes), n);
      points_id.push_back(std::distance(begin(*nodes), res));
      // points_id.push_back(n - *nodes->begin());
    }
    all_points_id.push_back(points_id);
    cell_types.push_back(points_id.size()*2 - 1);  // NOTE this works for line (type 3), triangle (type 5), be careful with quad (type 9), but closed polygon is type 7, needs to be enhanced for bricks etc...
    offset += points_id.size();
    offsets.push_back(offset);
    damage.push_back(el->giveValue("damage"));
    points_id.clear();
  }
  // for (unsigned e; e < elems->giveSize(); e++){
  // }

  giveFileName(step, buffer);
  ofstream outputfile((GlobPaths::RESULTDIR / buffer).string());

  if ( outputfile.is_open() ) {
      outputfile << std :: scientific;
      outputfile << "<?xml version=\"1.0\"?>\n";
      outputfile << "<VTKFile type=\"UnstructuredGrid\" version=\"0.1\" byte_order=\"LittleEndian\">" << '\n';
      outputfile << "<UnstructuredGrid>" << '\n';
      outputfile << "<Piece NumberOfPoints=\"" << nodes->giveSize() << "\" NumberOfCells=\"" << elems->giveSize() << "\">" << '\n';


      outputfile << "<Points>" << '\n';
      outputfile << "<DataArray type=\"Float32\" NumberOfComponents=\"3\" format=\"ascii\">" << '\n';
      for (auto const &n : *nodes){
        displ.push_back(Point(n->giveDoFBasedValue("ux", DoFs),
                              n->giveDoFBasedValue("uy", DoFs),
                              n->giveDoFBasedValue("uz", DoFs)));  // this is not correct for 2D, because it gives the rotation
        outputfile << n->givePoint().getX() << "\t" << n->givePoint().getY() << "\t"<< n->givePoint().getZ() << '\n';
      }
      // for (unsigned n; n < nodes->giveSize(); n++) {
      //   P = nodes->giveNode(n)->givePoint();
      // }
      outputfile << "</DataArray>" << "\n";
      outputfile << "</Points>" << "\n";
      // /*
      outputfile << "<Cells>" << '\n';
      outputfile << "<DataArray type=\"Int32\" Name=\"connectivity\" format=\"ascii\">" << '\n';
      for (auto const & value : all_points_id){
        for (auto const & id : value){
          outputfile << "\t" << id;
        }
        outputfile << '\n';
      }
      outputfile << "</DataArray>" << '\n';
      outputfile << "<DataArray type=\"Int32\" Name=\"offsets\" format=\"ascii\">" << '\n';
      for (auto const & value : offsets){
        outputfile << value << '\n';
      }
      outputfile << "</DataArray>" << '\n';
      outputfile << "<DataArray type=\"UInt8\" Name=\"types\" format=\"ascii\">" << '\n';
      for (auto const & value : cell_types){
        outputfile << value << '\n';
      }
      outputfile << "</DataArray>" << '\n';
      outputfile << "</Cells>" << '\n';
      // */
      outputfile << "<PointData Scalars=\"scalars\">" << "\n";
      outputfile << "<DataArray type=\"Float32\" Name=\"displacement\" NumberOfComponents=\"3\" format=\"ascii\">" << '\n';
      for (auto const & p : displ){
        outputfile << p.getX() << '\t' << p.getY() << '\t' << p.getZ() << '\n';
      }
      outputfile << "</DataArray>" << '\n';
      outputfile << "</PointData>" << '\n';
      // /*
      outputfile << "<CellData Scalars=\"scalars\">" << "\n";
      outputfile << "<DataArray type=\"Float32\" Name=\"damage\" format=\"ascii\">" << '\n';
      for (auto const & value : damage){
        outputfile << value << '\n';
      }
      outputfile << "</DataArray>" << '\n';
      outputfile << "</CellData>" << '\n';
      // */
      outputfile << "</Piece>" << '\n';
      outputfile << "</UnstructuredGrid>" << '\n';
      outputfile << "</VTKFile>" << '\n';
      outputfile.close();
  }
}

//////////////////////////////////////////////////////////
// function tahat calculates displacement of any point of rigid body from its rotations and
Point calculateVertexDisplacement2D(const RigidBodyContact &rbc, const Node &v, const Node &a, const Vector &DoFs){
  Matrix A = rbc.giveAMatrix(a.givePoint(), v.givePoint());
  Matrix U(3, 1);
  U [ 0 ] [ 0 ] = a.giveDoFBasedValue("ux", DoFs);
  U [ 1 ] [ 0 ] = a.giveDoFBasedValue("uy", DoFs);
  U [ 2 ] [ 0 ] = a.giveDoFBasedValue("uz", DoFs);  // in 2D, this DOF stays for rotation

  Matrix P = A * U;

  Point p = Point(P[0][0], P[1][0], 0.);  // zero displacement in Z direction

  return p;
}
//////////////////////////////////////////////////////////

// RIGID POLYGONS TO VTU FILE
//////////////////////////////////////////////////////////
void VTKRBExporter :: exportData(unsigned step, const Vector &DoFs, const Vector &reactions) const{
  // Export of elements into vtu xml file format (vtu = vtk for unstructured grid)
  // NOTE this is messy construction of xml file, will be remade using some of xml libraries for cpp
  char buffer [ 100 ];
  // Point P;
  // Element *ee;

  vector< Point > all_vertices_twice;  // because the vertices are needed twice (on each side of the contact)
  vector< Point > vertices_displ;  // displaceemnt of the contact
  vector< int > points_id;
  vector< int > node_id;
  vector< vector< int > > all_points_id;
  vector<int> cell_types;
  vector<int> offsets;
  vector<Point> displ;
  vector<double> damage;  // test version, this and more will be specified on the exporter input
  int offset = 0;
  for (auto const &el : *elems){
    RigidBodyContact *rbc = static_cast< RigidBodyContact * >( el );
    for (auto const &n : rbc->giveNodes()){
      auto nod_id_ptr = std::find(begin(*nodes), end(*nodes), n);
      for (auto const &v : rbc->giveVertices()){
        all_vertices_twice.push_back(v->givePoint());
        points_id.push_back(nodes->giveSize()+all_vertices_twice.size()-1);
        vertices_displ.push_back(calculateVertexDisplacement2D(*rbc, *v, *n, DoFs));
      }
      points_id.push_back(std::distance(begin(*nodes), nod_id_ptr));
      node_id.push_back(std::distance(begin(*nodes), nod_id_ptr));
      all_points_id.push_back(points_id);
      cell_types.push_back(points_id.size()*2 - 1);  // NOTE this works for line (type 3), triangle (type 5), be careful with quad (type 9), but closed polygon is type 7, needs to be enhanced for bricks etc...
      offset += points_id.size();
      offsets.push_back(offset);
      damage.push_back(el->giveValue("damage"));
      points_id.clear();
    }
  }
  // for (unsigned e; e < elems->giveSize(); e++){
  // }

  giveFileName(step, buffer);
  ofstream outputfile((GlobPaths::RESULTDIR / buffer).string());

  if ( outputfile.is_open() ) {
      outputfile << std :: scientific;
      outputfile << "<?xml version=\"1.0\"?>\n";
      outputfile << "<VTKFile type=\"UnstructuredGrid\" version=\"0.1\" byte_order=\"LittleEndian\">" << '\n';
      outputfile << "<UnstructuredGrid>" << '\n';
      outputfile << "<Piece NumberOfPoints=\"" << nodes->giveSize()+all_vertices_twice.size() << "\" NumberOfCells=\"" << elems->giveSize()*2 << "\">" << '\n';


      outputfile << "<Points>" << '\n';
      outputfile << "<DataArray type=\"Float32\" NumberOfComponents=\"3\" format=\"ascii\">" << '\n';
      for (auto const &n : *nodes){
        displ.push_back(Point(n->giveDoFBasedValue("ux", DoFs),
                              n->giveDoFBasedValue("uy", DoFs),
                              0));  // for 2D this DOF stays for rotation
        outputfile << n->givePoint().getX() << "\t" << n->givePoint().getY() << "\t"<< n->givePoint().getZ() << '\n';
      }
      for (auto const &p : all_vertices_twice){
        outputfile << p.getX() << "\t" << p.getY() << "\t"<< p.getZ() << '\n';
      }
      // for (unsigned n; n < nodes->giveSize(); n++) {
      //   P = nodes->giveNode(n)->givePoint();
      // }
      outputfile << "</DataArray>" << "\n";
      outputfile << "</Points>" << "\n";
      // /*
      outputfile << "<Cells>" << '\n';
      outputfile << "<DataArray type=\"Int32\" Name=\"connectivity\" format=\"ascii\">" << '\n';
      for (auto const & value : all_points_id){
        for (auto const & id : value){
          outputfile << "\t" << id;
        }
        outputfile << '\n';
      }
      outputfile << "</DataArray>" << '\n';
      outputfile << "<DataArray type=\"Int32\" Name=\"offsets\" format=\"ascii\">" << '\n';
      for (auto const & value : offsets){
        outputfile << value << '\n';
      }
      outputfile << "</DataArray>" << '\n';
      outputfile << "<DataArray type=\"UInt8\" Name=\"types\" format=\"ascii\">" << '\n';
      for (auto const & value : cell_types){
        outputfile << value << '\n';
      }
      outputfile << "</DataArray>" << '\n';
      outputfile << "</Cells>" << '\n';
      // */
      outputfile << "<PointData Scalars=\"scalars\">" << "\n";
      outputfile << "<DataArray type=\"Float32\" Name=\"displacement\" NumberOfComponents=\"3\" format=\"ascii\">" << '\n';
      for (auto const & p : displ){
        outputfile << p.getX() << '\t' << p.getY() << '\t' << p.getZ() << '\n';
      }
      for (auto const &p : vertices_displ){
        outputfile << p.getX() << '\t' << p.getY() << '\t' << p.getZ() << '\n';
      }
      outputfile << "</DataArray>" << '\n';
      outputfile << "</PointData>" << '\n';
      // /*
      outputfile << "<CellData Scalars=\"scalars\">" << "\n";
      outputfile << "<DataArray type=\"Float32\" Name=\"damage\" format=\"ascii\">" << '\n';
      for (auto const & value : damage){
        outputfile << value << '\n';
      }
      outputfile << "</DataArray>" << '\n';

      outputfile << "<DataArray type=\"Float32\" Name=\"node_id\" format=\"ascii\">" << '\n';
      for (auto const & value : node_id){
        outputfile << value << '\n';
      }
      outputfile << "</DataArray>" << '\n';

      outputfile << "</CellData>" << '\n';
      // */
      outputfile << "</Piece>" << '\n';
      outputfile << "</UnstructuredGrid>" << '\n';
      outputfile << "</VTKFile>" << '\n';
      outputfile.close();
  }
}
