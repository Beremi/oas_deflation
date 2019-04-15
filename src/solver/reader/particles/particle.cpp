/**
 * @Author: jose
 * @Date:   2019-04-05T17:54:28+02:00
 * @Last modified by:   jose
 * @Last modified time: 2019-04-10T12:12:51+02:00
 */

#include "particle.h"

// template <class T> T Particle::get_info(const std::string &str) const {
//   if (str.compare("id")==0) return id;
//   else if (str.compare("radius")==0) return radius;
//   else if (str.compare("label")==0) return label;
//   else if (str.compare("node_id")==0) return node_id;
//   else if (str.compare("volume")==0) return volume;
//   else {
//     std::cerr << "no such variable in class Particle: " << str << '\n';
//   }
// }


// template <class T> void Particle::set_info(const std::string &str, const T &parameter){
//   if (str.compare("id")==0) id = parameter;
//   else if (str.compare("radius")==0) radius = parameter;
//   else if (str.compare("label")==0) label = parameter;
//   else if (str.compare("node_id")==0) node_id = parameter;
//   else if (str.compare("volume")==0) volume = parameter;
//   else {
//     std::cerr << "no such variable in class Particle: " << str << '\n';
//   }
// }


int Particle::get_id() const{
  return id;
}

void Particle::set_id(const int &i){
  id = i;
}

int Particle::get_node_id() const{
  return node_id;
}

void Particle::set_node_id(const int &i){
  node_id = i;
}

double Particle::get_radius() const{
  return radius;
}

void Particle::set_radius(const double &r){
  radius = r;
}

std::string Particle::get_label() const{
  return label;
}

void Particle::set_label(const std::string &s){
  label = s;
}

double Particle::get_volume() const{
  return volume;
}

void Particle::set_volume(const double &v){
  volume = v;
}

double Particle::get_lin_mom() const {
  return linear_momentum;
}

void Particle::set_lin_mom(const double &l){
  linear_momentum = l;
}

Point Particle::get_centroid() const{
  return centroid;
}
void Particle::set_centroid(const Point &P){
  centroid = P;
}

Point Particle::get_static_mom() const{
  return static_momentum;
}

void Particle::set_static_mom(const Point &P){
  static_momentum = P;
}

Matrix Particle::get_rot_mom() const{
  return rotational_momentum;
}

void Particle::set_rot_mom(const Matrix &M){
  rotational_momentum = M;
}
