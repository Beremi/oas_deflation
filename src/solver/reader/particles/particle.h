/**
 * @Author: jose
 * @Date:   2019-04-05T17:54:51+02:00
 * @Last modified by:   jose
 * @Last modified time: 2019-04-10T12:08:54+02:00
 */

#include <string>
#include <vector>
#include "../../common/matrixops.h"
#include "../../common/point.h"

#ifndef _PARTICLE_H
#define	_PARTICLE_H


class Particle {
private:
  int id;  // does it need to be here?
  double radius;
  std::string label;
  int node_id; // Voronoi point - alll the nodes and vertices will be in one
  double volume;
  double linear_momentum; //added for inertia tensor and mass matrix
  Point centroid; // center of gravity
  Point static_momentum;
  Matrix rotational_momentum; //added for inertia tensor and mass matrix

public:
  Particle ();
  virtual ~Particle ();


  int get_id() const;
  void set_id(const int &i);

  int get_node_id() const;
  void set_node_id(const int &i);

  double get_radius() const;
  void set_radius(const double &r);

  std::string get_label() const;
  void set_label(const std::string &s);

  double get_volume() const;
  void set_volume(const double &v);

  // template <class T> T get_info(const std::string &str) const ;
  // template <class T> void set_info(const std::string &str, const T &parameter);

  double get_lin_mom() const;
  void set_lin_mom(const double &l);

  Point get_centroid() const;
  void set_centroid(const Point &P);

  Point get_static_mom() const;
  void set_static_mom(const Point &P);

  Matrix get_rot_mom() const;
  void set_rot_mom(const Matrix &M);

};



#endif /* _PARTICLE_H */
