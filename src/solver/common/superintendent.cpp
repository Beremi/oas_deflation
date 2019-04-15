/**
 * @Author: jose
 * @Date:   2019-04-05T20:22:31+02:00
 * @Last modified by:   jose
 * @Last modified time: 2019-04-10T19:08:55+02:00
 */

#include "superintendent.h"

std::string Superintendent::get_fol() const{
  return fol_name;
}

void Superintendent::set_fol(const std::string &s){
  fol_name = s;
}
