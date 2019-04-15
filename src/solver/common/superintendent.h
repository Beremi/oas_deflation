/**
 * @Author: jose
 * @Date:   2019-04-05T20:22:48+02:00
 * @Last modified by:   jose
 * @Last modified time: 2019-04-10T19:07:51+02:00
 */

#include <stdio.h>
#include <string>

#ifndef _SUPERINTENDENT_H
#define	_SUPERINTENDENT_H


class Superintendent {
private:
  std::string fol_name;

public:
  Superintendent ();
  virtual ~Superintendent ();

  std::string get_fol() const;
  void set_fol(const std::string &s);
};


#endif /* _SUPERINTENDENT_H */
