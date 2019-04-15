/**
 * @Author: jose
 * @Date:   2019-04-10T15:48:34+02:00
 * @Last modified by:   jose
 * @Last modified time: 2019-04-10T17:56:43+02:00
 */

#include "node.h"
#include <sstream>
#include <fstream>

using namespace std;

Node::Node(const string &str){
  istringstream iss(str);
  double x, y, z, rad;
  // NOTE how to distinguish between 2D and 3D?
  iss >> x >> y >> z >> rad;
  N.set(x, y, z);
  r = rad;
}

Node::Node(istringstream &iss){
  double x, y, z, rad;
  // NOTE how to distinguish between 2D and 3D?
  iss >> x >> y >> z >> rad;
  N.set(x, y, z);
  r = rad;
}

int Node::get_id() const{
  return id;
}

void Node::set_id(const int &i){
  id = i;
}

Point Node::get_node() const{
  return N;
}
void Node::set_node(const Point &P){
  N = P;
}

double Node::get_radius() const{
  return r;
}

void Node::set_radius(const double &a){
  r = a;
}


void readNodes(const std::string &inputfile, std::vector<Node> &Nodes){
  ifstream indata(inputfile); // indata is like cin

  string line;
  int line_num = 0;
  Node N;
  for (string line; std::getline(indata, line); line_num++){
    // if (line.front() == '#') continue; // to skip # comments
    // NOTE but it is more convenient to just skip the first line in this case, since the only comment should be there, otherwise it would check the beginning of every line
    // istringstream iss(line);
    // N = Node(iss);  // other type of constructor
    N = Node(line);
    Nodes.push_back(N);
  }
}
