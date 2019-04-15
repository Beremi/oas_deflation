/**
 * @Author: jose
 * @Date:   2019-04-05T19:01:08+02:00
 * @Last modified by:   jose
 * @Last modified time: 2019-04-10T19:34:26+02:00
 */


#include "loadGeometry.h"
#include <fstream>
using std::ifstream;

using namespace std;


void loadInput(Superintendent &SUP, std::vector<std::vector<Element *> > &ELEMENTS, std::vector<Node> &Nodes, std::vector<Vertex> &Vertices){
  string path_to_files;

  // first, load nodes
  path_to_files = SUP.get_fol();
  path_to_files.append("/nodes_out.txt");
  readNodes(path_to_files, Nodes);
  //  and vertices
  path_to_files = SUP.get_fol();
  path_to_files.append("/vertices_out.txt");
  readVertices(path_to_files, Vertices);

  // now materials
  // TODO this is more complicated, same as elements, need vector of pointers to different materials , which should not be a problem, just not sure if it works correctly (vector of pointers declared as to parents but in reality to a derrived class instances - virtual constructors in parent class needed)




}

// vector<node> readNodes(const string &label){
//   vector<node> Nodes;
//
//   string inputfile = label;
//   inputfile.append("/nodes_out.txt");
//
//   ifstream indata; // indata is like cin
//   string readstr, rf_name;
//   char t[inputfile.length()];
//   strcpy(t, inputfile.c_str());
//   indata.open(t); // opens the file
//   if (!indata) { // file couldn't be opened
//       cerr << "Error: input file " << inputfile << " could not be opened" << endl;
//       exit(1);
//   }
//   double x, y, z, r;
//   node A;
//   while (!indata.eof()) {
//     indata >> readstr;
//     if (readstr.compare("#") == 0){
//       continue;
//     }
//     x = std::stod(readstr);
//     indata >> y >> z >> r;
//     A.center = Point(x, y, z);
//     A.radius = r;
//
//     Nodes.push_back(A);
//   }
//   return Nodes;
// }
//
// vector<Point> readVertices(const string &label){
//   vector<Point> Vertices;
//
//   string inputfile = label;
//   inputfile.append("/veritces_out.txt");
//
//   ifstream indata; // indata is like cin
//   string readstr, rf_name;
//   char t[inputfile.length()];
//   strcpy(t, inputfile.c_str());
//   indata.open(t); // opens the file
//   if (!indata) { // file couldn't be opened
//       cerr << "Error: input file " << inputfile << " could not be opened" << endl;
//       exit(1);
//   }
//   double x, y, z, r;
//   Point P;
//   while (!indata.eof()) {
//     indata >> readstr;
//     if (readstr.compare("#") == 0){
//       continue;
//     }
//     x = std::stod(readstr);
//     indata >> y >> z;
//     P = Point(x, y, z);
//
//     Vertices.push_back(P);
//   }
//   return Vertices;
// }
//
// void readElems(const string &label, const vector<Point> &Vertices, const vector<node> &Nodes, vector<Element> &Elems){
//
//     string inputfile = label;
//     inputfile.append("/elems_out.txt");
//
//     ifstream indata; // indata is like cin
//     string readstr, rf_name;
//     char t[inputfile.length()];
//     strcpy(t, inputfile.c_str());
//     indata.open(t); // opens the file
//     if (!indata) { // file couldn't be opened
//         cerr << "Error: input file " << inputfile << " could not be opened" << endl;
//         exit(1);
//     }
//     double x, y, z, r;
//     int a, b, c;
//     int id = 0;
//     string mat;
//     Element E;
//     while (!indata.eof()) {
//       indata >> readstr;
//       if (readstr.compare("#") == 0){
//         continue;
//       }
//       if (readstr.compare("M_ELEM") == 0){
//         // toto schovat do samostatnýho readeru pro daný typ elementu
//         // TODO what for T_ELEM? (transport elements) or any other?
//         indata >> a >> b;
//         E = Element(id, a, b, Nodes[a].center, Nodes[b].center, Nodes[a].radius, Nodes[b].radius);
//         indata >> c;
//         for (int ci = 0 ; ci < c; c++){
//           indata >> a;
//           E.vertices.push_back(Vertices[a]);
//         }
//         indata >> mat;
//         E.material = mat;
//         id ++;
//       }
//       Elems.push_back(E);
//     }
// }
//
// void ArrangeAggregates(const Experiment &EXP, vector<Element> &E, vector<Aggregate> &AGG){
//   // TODO
// }
//
// void readPrepData(const Experiment &EXP, vector<Element> &E, vector<Aggregate> &AGG){
//   vector<node> Nodes = readNodes(EXP.label);
//   vector<Point> Vertices = readVertices(EXP.label);
//   readElems(EXP.label, Vertices, Nodes, E);
//
// }
