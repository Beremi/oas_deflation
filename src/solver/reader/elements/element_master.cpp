/**
 * @Author: jose
 * @Date:   2019-04-05T15:24:29+02:00
 * @Last modified by:   jose
 * @Last modified time: 2019-04-10T18:10:11+02:00
 */

 #include "element_master.h"

std::string Element::get_type() const{
  return type;
}

std::vector<int> Element::get_nodes() const{
  return nodes;
}


void readElements(const std::string &inputfile);
