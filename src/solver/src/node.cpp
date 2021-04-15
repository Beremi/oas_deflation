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
    std :: ostringstream out;
    out.precision(10);
    out << this->name << '\t';

    out << std :: fixed << this->point.getX();
    out << '\t';
    out << std :: fixed << this->point.getY();
    if ( this->dim == 3 ) {
        out << '\t';
        out << std :: fixed << this->point.getZ();
    }
    return out.str();
}

//////////////////////////////////////////////////////////
unsigned Node :: giveOrderOfForceCode(string code) const {
    char* p;
    long converted = strtol(code.c_str(), &p, 10);
    if (not *p) return converted;
    else{
        cerr << name << " error: there is no force corresponding to code "<< code << endl;
        exit(1);
    }
    return 0;
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
unsigned MechNode :: giveOrderOfForceCode(string code) const {    
    if ( code.compare("fx") == 0 ) return 0;
    else if ( code.compare("fy") == 0 ) return 1;
    else if ( dim >2 && code.compare("fz") == 0 ) return 2;
    else return Node::giveOrderOfForceCode(code);
}


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
unsigned TrsNode :: giveOrderOfForceCode(string code) const {
    if ( code.compare("flux") == 0 || code.compare("fx") == 0 ) return 0;
    else return Node::giveOrderOfForceCode(code);
}

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
// PARTICLE - translational and rotational DoFs
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
unsigned Particle :: giveOrderOfForceCode(string code) const {
    if ( dim==3 && code.compare("mx") == 0 ) return 3;
    else if ( dim==3 && code.compare("my") == 0 ) return 4;
    else if (code.compare("mz") == 0 ) {
        if (dim==2) return 2;
        else return 5;
    }
    else return MechNode::giveOrderOfForceCode(code);
}

//////////////////////////////////////////////////////////
std :: string Particle :: giveLineToSave() const {
    std :: ostringstream out;
    out.precision(10);
    out << Node :: giveLineToSave() << "\t";
    out << std :: fixed << this->r;
    return out.str();
};

//////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////
// COUPLED PARTICLE - translational and rotational and transport DoFs
//////////////////////////////////////////////////////////
double CoupledParticle :: giveDoFBasedValue(string code, const Vector &DoFs) const {
    if ( code.compare("pressure") == 0 ) {
        return DoFs [ firstDoF + 6 ];
    } else {
        return Particle :: giveDoFBasedValue(code, DoFs);
    }
};

//////////////////////////////////////////////////////////
unsigned CoupledParticle :: giveOrderOfForceCode(string code) const {
    if (code.compare("flux") == 0 ) {
        if (dim==2) return 3;
        else return 6;
    }
    else return Particle::giveOrderOfForceCode(code);
}
