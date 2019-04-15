/**
 * @Author: jose
 * @Date:   2019-04-10T16:21:24+02:00
 * @Last modified by:   jose
 * @Last modified time: 2019-04-10T17:57:04+02:00
 */

#include "vertex.h"
#include <sstream>
#include <fstream>

using namespace std;

Vertex::Vertex(const string &str){
  istringstream iss(str);
  double x, y, z;
  // NOTE how to distinguish between 2D and 3D?
  iss >> x >> y >> z;
  V.set(x, y, z);
}

void readVertices(const std::string &inputfile, std::vector<Vertex> &Vertices){
  ifstream indata(inputfile); // indata is like cin

  string line;
  int line_num = 0;
  Vertex V;
  for (string line; std::getline(indata, line); line_num++){
    if (line_num == 0) continue;

    V = Vertex(line);
    Vertices.push_back(V);
  }
}
