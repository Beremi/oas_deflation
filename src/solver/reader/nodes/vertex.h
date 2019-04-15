/**
 * @Author: jose
 * @Date:   2019-04-10T16:21:36+02:00
 * @Last modified by:   jose
 * @Last modified time: 2019-04-10T16:58:55+02:00
 */

#include "../../common/point.h"
#include <vector>

#ifndef _VERTEX_H
#define _VERTEX_H


class Vertex {
private:
  int id;
  Point V;

public:
  Vertex ();
  Vertex(const std::string &str);
  virtual ~Vertex ();

  int get_id() const{
    return id;
  }
  void set_id(const int &i){
    id = i;
  }

  Point get_vertex() const {
    return V;
  }
  void set_vertex(const Point &P){
    V = P;
  }
};

void readVertices(const std::string &inputfile, std::vector<Vertex> &Vertices);

#endif  /* _VERTEX_H */
