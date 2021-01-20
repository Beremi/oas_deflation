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
  std :: string  str = this->name + "\t" + to_string(this->point.getX()) +
      "\t" + to_string(this->point.getY());
  if ( this->dim == 3 ){
    str += "\t" + to_string(this->point.getZ());
  }
  return str;
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
  return Node :: giveLineToSave() + "\t" + to_string(this->r);
};
