/**
 * @Author: jose
 * @Date:   2019-04-10T19:19:51+02:00
 * @Last modified by:   jose
 * @Last modified time: 2019-04-10T19:30:39+02:00
 */

#include <string>

#ifndef _MATERIAL_H
#define	_MATERIAL_H

typedef struct {
  double ela_mod;
  double alpha;
  double conductivity;
  double capacity;
  double density;
} LinearElastic;


class Materials {
private:
  LinearElastic linel;

public:
  Materials ();
  Materials (const std::string &filename);
  virtual ~Materials ();

  LinearElastic get_linel_data() const;
  void set_linel_data(const LinearElastic &le);
};

#endif
