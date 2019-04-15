/**
 * @Author: jose
 * @Date:   2019-04-05T15:25:52+02:00
 * @Last modified by:   jose
 * @Last modified time: 2019-04-11T10:37:53+02:00
 */

#include "element_struct.h"
#include "../../common/matrixops.h"
#include "../../common/point.h"
#include <sstream>

#ifndef _ELEMENT_DISCRETE_H
#define	_ELEMENT_DISCRETE_H

// TODO nebude elastic, ale obecně, ostatní bude řešit až materiál
class Element_discrete_elastic: public Element_struct {
private:
  std::vector<int> node_ids; // ids of particle centers it connects
  // particles are coincidental with nodes only in case of pure discrete model
  // (no FEM etc.)
  std::vector<int> vertices; // ids of particle centers it connects

  // Point mid;  // midpoint between two connecting nodes
  Point c;  // centroid of the facet
  Point n, t1, t2;  // normal and two tangential direction forming orthonormal basis at the facet
  double length;  // legth of element
  double area;

  // double ela_mod;
  // double alpha;
  // TODO pointer na material
  // TODO polke pointrů na statusy

  Vector eps;

  Vector Bn1;  // slice of matrix B (A1*n)
  Vector Bn2;  // slice of matrix B (A2*n)
  Vector Bt11;  // slice of matrix B (A1*t1)
  Vector Bt12;  // slice of matrix B (A1*t2)
  Vector Bt21;  // slice of matrix B (A2*t1)
  Vector Bt22;  // slice of matrix B (A2*t2)

public:
  Element_discrete_elastic () ;
  Element_discrete_elastic (const std::string &str);
  Element_discrete_elastic (std::istringstream &iss);
  ~Element_discrete_elastic ();

  std::vector<int> get_node_ids() const ;
  void set_node_ids(const std::vector<int> &ids);

  std::vector<int> get_vertices() const ;
  void set_vertices(const std::vector<int> &ids);

  // Point get_point(const char &a) const ;
  // void set_point(const char &a, const Point &P);

  Point get_c() const ;
  void set_c(const Point &P);

  Point get_n() const ;
  void set_n(const Point &P);

  Point get_t1() const ;
  void set_t1(const Point &P);

  Point get_t2() const ;
  void set_t2(const Point &P);

  double get_area() const ;
  void set_area(const double &a);

  double get_length() const ;
  void set_length(const double &l);

  double get_ela_mod() const ;
  void set_ela_mod(const double &s);

  double get_alpha() const ;
  void set_alpha(const double &a);

  Vector get_epsilon() const ;
  void set_epsilon(const Vector &e);

  Vector get_vectorB(const int &b) const ;
  void set_vectorB(const int &b, const Vector &B);

};

#endif	/* _ELEMENT_DISCRETE_H */
