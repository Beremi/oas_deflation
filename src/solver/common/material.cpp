/**
 * @Author: jose
 * @Date:   2019-04-10T19:20:12+02:00
 * @Last modified by:   jose
 * @Last modified time: 2019-04-10T19:27:45+02:00
 */

#include "material.h"

LinearElastic Materials::get_linel_data() const{
  return linel;
}
void Materials::set_linel_data(const LinearElastic &le){
  linel = le;
}
