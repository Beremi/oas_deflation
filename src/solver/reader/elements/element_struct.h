/**
 * @Author: jose
 * @Date:   2019-04-05T15:25:52+02:00
 * @Last modified by:   jose
 * @Last modified time: 2019-04-11T10:22:34+02:00
 */

#include "element_master.h"
#include "../../common/matrixops.h"

#ifndef _ELEMENT_STRUCT_H
#define	_ELEMENT_STRUCT_H

class Element_struct: public Element {
private:
  //  tady možná jen funkce get_matici
  Matrix* K;  // for different elements can be different matrix size
  Matrix* M;

public:
  Element_struct ();
  ~Element_struct ();
  Matrix* get_matrix(const char &a) const;
};

#endif	/* _ELEMENT_STRUCT_H */
