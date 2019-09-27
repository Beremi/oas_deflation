#include "node.h"
#include "boundary_condition.h"

//////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////
// BASIC NODE - ONLY MASTER CLASS
unsigned Node :: giveNumberOfFreeDoFs() const {
    if ( !bc ) {
        return nDoFs;
    } else {
        return nDoFs - bc->giveNumberOfBlockedDoFs();
    }
}

//////////////////////////////////////////////////////////
void Node :: readFromLine(istringstream &iss, int dim) {
    double x, y, z;
    if ( dim == 2 ) {
        iss >> x >> y;
        point = Point(x, y);
    } else if ( dim == 3 )      {
        iss >> x >> y >> z;
        point = Point(x, y, z);
    }
}

//////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////
// MASTER DOF - GOVERN DEPENDENT DOFs
void MasterDoF :: readFromLine(istringstream &iss, int dim){
  double x, y, z;
  if ( dim == 2 ) {
      iss >> x >> y;
      point = Point(x, y);
  } else if ( dim == 3 )      {
      iss >> x >> y >> z;
      point = Point(x, y, z);
  }
  unsigned i;
  iss >> i;
  nDoFs = i;
}

//////////////////////////////////////////////////////////
// MASTER DOF - GOVERN DEPENDENT DOFs
MasterDoF :: MasterDoF(Point c, unsigned n){
  point = c;
  nDoFs = n;
  name = "Master DoF";
}

//////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////
// MASTER Node - GOVERN multiple dependent DOFs
void MasterNode :: readFromLine(istringstream &iss, int dim){
  double x, y, z;
  if ( dim == 2 ) {
      iss >> x >> y;
      point = Point(x, y);
  } else if ( dim == 3 )      {
      iss >> x >> y >> z;
      point = Point(x, y, z);
  }
}

//////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////
// MECHANICAL NODE - translational DoFs
double MechNode :: giveDoFBasedValue(string code, const Vector &DoFs) const {
    if ( code.compare("ux") == 0 ) {
        return DoFs [ firstDoF ];
    } else if ( code.compare("uy") == 0 )   {
        return DoFs [ firstDoF + 1 ];
    } else if ( code.compare("uz") == 0 )                                                                         {
        return DoFs [ firstDoF + 2 ];
    } else                                                                                                                                          {
        return Node :: giveDoFBasedValue(code, DoFs);
    }
};

//////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////
// TRANSPORT NODE - pressure DoF
double TrsNode :: giveDoFBasedValue(string code, const Vector &DoFs) const {
    if ( code.compare("pressure") == 0 ) {
        return DoFs [ firstDoF ];
    } else {
        return Node :: giveDoFBasedValue(code, DoFs);
    }
};

//////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////
// MECHANICAL NODE - translational and rotational DoFs
void Particle :: readFromLine(istringstream &iss, int dim) {
    double x, y, z;
    if ( dim == 2 ) {
        iss >> x >> y >> r;
        point = Point(x, y);
    } else if ( dim == 3 )      {
        iss >> x >> y >> z >> r;
        point = Point(x, y, z);
    }
}

//////////////////////////////////////////////////////////
double Particle :: giveDoFBasedValue(string code, const Vector &DoFs) const {
    if ( code.compare("rotx") == 0 ) {
        return DoFs [ firstDoF + 3 ];
    } else if ( code.compare("roty") == 0 )   {
        return DoFs [ firstDoF + 4 ];
    } else if ( code.compare("rotz") == 0 )                                                                           {
        return DoFs [ firstDoF + 5 ];
    } else                                                                                                                                              {
        return MechNode :: giveDoFBasedValue(code, DoFs);
    }
};
