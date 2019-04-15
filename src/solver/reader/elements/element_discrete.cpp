/**
 * @Author: jose
 * @Date:   2019-04-05T15:21:45+02:00
 * @Last modified by:   jose
 * @Last modified time: 2019-04-10T18:06:30+02:00
 */

#include "element_discrete.h"
using namespace std;


Element_discrete_elastic::Element_discrete_elastic (const string &str){

}

Element_discrete_elastic::Element_discrete_elastic (istringstream &iss){
  int id, num;
  for (int i = 0; i < 2; i++){
    iss >> id;
    node_ids.push_back(id);
  }
  iss >> num;
  for (int i = 0; i < num; i++){
    iss >> id;
    vertices.push_back(id);
  }
}

std::vector<int> Element_discrete_elastic::get_node_ids() const {
  return node_ids;
}

void Element_discrete_elastic::set_node_ids(const std::vector<int> &ids){
  node_ids = ids;
}

std::vector<int> Element_discrete_elastic::get_vertices() const {
  return vertices;
}

void Element_discrete_elastic::set_vertices(const std::vector<int> &ids){
  vertices = ids;
}

// Point Element_discrete_elastic::get_point(const char &a) const {
//   if (a == 'c') return c;
//   else if (a == 'n') return n;
//   else if (a == 'm') return t1;
//   else if (a == 'l') return t2;
//   else {
//     std::cerr << "no such point " << a << '\n';
//     exit(1);
//   }
// }
//
// void Element_discrete_elastic::set_point(const char &a, const Point &P){
//   if (a == 'c') c = P;
//   else if (a == 'n') n = P;
//   else if (a == 'm') t1 = P;
//   else if (a == 'l') t2 = P;
//   else {
//     std::cerr << "no such point " << a << '\n';
//     exit(1);
//   }
// }

Point Element_discrete_elastic::get_c() const {
  return c;
}

void Element_discrete_elastic::set_c(const Point &P){
  c = P;
}

Point Element_discrete_elastic::get_n() const {
  return n;
}

void Element_discrete_elastic::set_n(const Point &P){
  n = P;
}

Point Element_discrete_elastic::get_t1() const {
  return t1;
}

void Element_discrete_elastic::set_t1(const Point &P){
  t1 = P;
}

Point Element_discrete_elastic::get_t2() const {
  return t2;
}

void Element_discrete_elastic::set_t2(const Point &P){
  t2 = P;
}

double Element_discrete_elastic::get_area() const {
  return area;
}
void Element_discrete_elastic::set_area(const double &a){
  area = a;
}

double Element_discrete_elastic::get_length() const {
  return length;
}

void Element_discrete_elastic::set_length(const double &l){
  length = l;
}

double Element_discrete_elastic::get_ela_mod() const {
  return ela_mod;
}

void Element_discrete_elastic::set_ela_mod(const double &s){
  ela_mod = s;
}

double Element_discrete_elastic::get_alpha() const {
  return alpha;
}

void Element_discrete_elastic::set_alpha(const double &a){
  alpha = a;
}

Vector Element_discrete_elastic::get_epsilon() const {
  // TODO does this work? eps je array double
  return eps;
}

void Element_discrete_elastic::set_epsilon(const Vector &e){
  eps = e;
}

Vector Element_discrete_elastic::get_vectorB(const int &b) const {
  if (b == 1) return Bn1;
  else if (b == 2) return Bn2;
  else if (b == 11) return Bt11;
  else if (b == 12) return Bt12;
  else if (b == 21) return Bt21;
  else if (b == 22) return Bt22;
  else {
    std::cerr << "no such vector B " << b << '\n';
    exit(1);
  }

}

void Element_discrete_elastic::set_vectorB(const int &b, const Vector &B){
  if (b == 1) Bn1 = B;
  else if (b == 2) Bn2 = B;
  else if (b == 11) Bt11 = B;
  else if (b == 12) Bt12 = B;
  else if (b == 21) Bt21 = B;
  else if (b == 22) Bt22 = B;
  else {
    std::cerr << "no such vector B " << b << '\n';
    exit(1);
  }
}
