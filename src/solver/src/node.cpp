#include "node.h"
#include "boundary_condition.h"

//////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////
// BASIC NODE - ONLY MASTER CLASS

void Node :: readFromLine(istringstream &iss) {
    double x, y, z;
    if ( dim == 2 ) {
        iss >> x >> y;
        point = Point(x, y, 0.);
    } else if ( dim == 3 ) {
        iss >> x >> y >> z;
        point = Point(x, y, z);
    }
}


std :: string Node :: giveLineToSave() const {
  std::ostringstream out;
  out.precision(10);
  out << this->name << '\t';

  out << std::fixed << this->point.getX();
  out << '\t';
  out << std::fixed << this->point.getY();
  if ( this->dim == 3 ){
    out << '\t';
    out << std::fixed << this->point.getZ();
  }
  return out.str();
}

//////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////
// MASTER DOF - GOVERN DEPENDENT DOFs
void MechDoF :: readFromLine(istringstream &iss) {
    ( void ) iss;
}

//////////////////////////////////////////////////////////
MechDoF :: MechDoF(unsigned dimension) : MechNode(dimension) {
    point = Point(0, 0, 0);
    nDoFs = dimension;
    dim = dimension;
    name = "MechDoF";
}

//////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////
// MECHANICAL NODE - translational DoFs
double MechNode :: giveDoFBasedValue(string code, const Vector &DoFs) const {
    if ( code.compare("ux") == 0 ) {
        return DoFs [ firstDoF ];
    } else if ( dim > 1 && code.compare("uy") == 0 ) {
        return DoFs [ firstDoF + 1 ];
    } else if ( dim > 2 && code.compare("uz") == 0 ) {
        return DoFs [ firstDoF + 2 ];
    } else {
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
// MASTER DOF - GOVERN DEPENDENT DOFs
void TrsDoF :: readFromLine(istringstream &iss) {
    ( void ) iss;
}

//////////////////////////////////////////////////////////
// MASTER DOF - GOVERN DEPENDENT DOFs
TrsDoF :: TrsDoF(unsigned dimension) : TrsNode(dimension) {
    point = Point(0, 0, 0);
    nDoFs = dimension;
    dim = dimension;
    name = "TrsDoF";
}


//////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////
// MECHANICAL NODE - translational and rotational DoFs
void Particle :: readFromLine(istringstream &iss) {
    double x, y, z;
    if ( dim == 2 ) {
        iss >> x >> y >> r;
        point = Point(x, y, 0);
    } else if ( dim == 3 ) {
        iss >> x >> y >> z >> r;
        point = Point(x, y, z);
    }
}

//////////////////////////////////////////////////////////
double Particle :: giveDoFBasedValue(string code, const Vector &DoFs) const {
    if ( dim > 1 && code.compare("rotx") == 0 ) {
        return DoFs [ firstDoF + 3 ];
    } else if ( dim > 2 && code.compare("roty") == 0 ) {
        return DoFs [ firstDoF + 4 ];
    } else if ( dim > 2 && code.compare("rotz") == 0 ) {
        return DoFs [ firstDoF + 5 ];
    } else {
        return MechNode :: giveDoFBasedValue(code, DoFs);
    }
};


//////////////////////////////////////////////////////////
std :: string Particle :: giveLineToSave() const {
  std::ostringstream out;
  out.precision(10);
  out << Node :: giveLineToSave() << "\t";
  out << std :: fixed << this->r;
  return out.str();
};
