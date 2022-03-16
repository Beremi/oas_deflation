#include "node.h"
#include "boundary_condition.h"
#include "element_discrete.h"
#include "simplex.h"

using namespace std;

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

    out << std :: fixed << this->point.x();
    out << '\t';
    out << std :: fixed << this->point.y();
    if ( this->dim == 3 ) {
        out << '\t';
        out << std :: fixed << this->point.z();
    }
    return out.str();
}

//////////////////////////////////////////////////////////
unsigned Node :: giveOrderOfForceCode(string code) const {
    char *p;
    long converted = strtol(code.c_str(), & p, 10);
    if ( !* p ) {
        return converted;
    } else {
        cerr << name << " error: there is no force corresponding to code " << code << endl;
        exit(1);
    }
    return 0;
}

//////////////////////////////////////////////////////////
Simplex *Node :: addElementToSimplex(RigidBodyContact *rbc) {
    if ( !simplex ) {
        simplex = new Simplex(this);
    }
    simplex->addElement(rbc);
    return simplex;
}

//////////////////////////////////////////////////////////
void Node :: initSimplex() {
    if ( simplex ) {
        simplex->init(dim);
    }
}

//////////////////////////////////////////////////////////
void Node :: updateSimplexVolumetricStrain(const Vector &fullDoFs) {
    if ( simplex ) {
        simplex->computeVolumetricStrain(fullDoFs);
    }
}

//////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////
// MASTER DOF - GOVERN DEPENDENT DOFs
void MechDoF :: readFromLine(istringstream &iss) {
    double x, y, z;
    if ( dim == 2 ) {
        iss >> x >> y;
        point = Point(x, y, 0.);
    } else if ( dim == 3 ) {
        iss >> x >> y >> z;
        point = Point(x, y, z);
    }
    iss >> nDoFs;
    if ( iss.fail() ) {
        cerr << "MechDoF reading failed!" << endl;
        exit(1);
    }
}

//////////////////////////////////////////////////////////
MechDoF :: MechDoF(unsigned dimension, unsigned numDoFs) : MechNode(dimension) {
    point = Point(0, 0, 0);
    nDoFs = numDoFs;
    dim = dimension;
    name = "MechDoF";
}

//////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////
// MECHANICAL NODE - translational DoFs
void MechNode :: giveDoFBasedValues(string code, const Vector &DoFs, Vector &result) const {
    if ( code.compare("displacements") == 0 ) {
        result.resize(dim);
        for ( unsigned i = 0; i < dim; i++ ) {
            result [ i ] = DoFs [ firstDoF + i ];
        }
    } else if ( code.compare("ux") == 0 ) {
        result.resize(1);
        result [ 0 ] = DoFs [ firstDoF ];
    } else if ( dim > 1 && code.compare("uy") == 0 ) {
        result.resize(1);
        result [ 0 ] = DoFs [ firstDoF + 1 ];
    } else if ( dim > 2 && code.compare("uz") == 0 ) {
        result.resize(1);
        result [ 0 ] = DoFs [ firstDoF + 2 ];
    } else {
        Node :: giveDoFBasedValues(code, DoFs, result);
    }
};

//////////////////////////////////////////////////////////
unsigned MechNode :: giveOrderOfForceCode(string code) const {
    if ( code.compare("fx") == 0 ) {
        return 0;
    } else if ( code.compare("fy") == 0 ) {
        return 1;
    } else if ( dim > 2 && code.compare("fz") == 0 ) {
        return 2;
    } else {
        return Node :: giveOrderOfForceCode(code);
    }
}


//////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////
// TRANSPORT NODE - pressure DoF
void TrsNode :: giveDoFBasedValues(string code, const Vector &DoFs, Vector &result) const {
    if ( code.compare("pressure") == 0 ) {
        result.resize(1);
        result [ 0 ] = DoFs [ firstDoF ];
    } else {
        Node :: giveDoFBasedValues(code, DoFs, result);
    }
};

//////////////////////////////////////////////////////////
unsigned TrsNode :: giveOrderOfForceCode(string code) const {
    if ( code.compare("flux") == 0 || code.compare("fx") == 0 ) {
        return 0;
    } else {
        return Node :: giveOrderOfForceCode(code);
    }
}

//////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////
// TRANSPORT + TEMPERATURE NODE
void TrsTemprtrCoupledNode :: giveDoFBasedValues(string code, const Vector &DoFs, Vector &result) const {
    if ( code.compare("humidity") == 0 ) {
        result.resize(1);
        result [ 0 ] = DoFs [ firstDoF ];
    } else if ( code.compare("temperature") == 0 ) {
        result.resize(1);
        result [ 0 ] = DoFs [ firstDoF + 1 ];
    } else {
        Node :: giveDoFBasedValues(code, DoFs, result);
    }
};

//////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////
// MASTER DOF - GOVERN DEPENDENT DOFs
void TrsDoF :: readFromLine(istringstream &iss) {
    double x, y, z;
    if ( dim == 2 ) {
        iss >> x >> y;
        point = Point(x, y, 0.);
    } else if ( dim == 3 ) {
        iss >> x >> y >> z;
        point = Point(x, y, z);
    }
    iss >> nDoFs;
    if ( iss.fail() ) {
        cerr << "TrsDoF reading failed!" << endl;
        exit(1);
    }
}

//////////////////////////////////////////////////////////
// MASTER DOF - GOVERN DEPENDENT DOFs
TrsDoF :: TrsDoF(unsigned dimension, unsigned numDoFs) : TrsNode(dimension) {
    point = Point(0, 0, 0);
    nDoFs = numDoFs;
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
void Particle :: giveDoFBasedValues(string code, const Vector &DoFs, Vector &result) const {
    if ( dim > 1 && code.compare("rotx") == 0 ) {
        result.resize(1);
        result [ 0 ] = DoFs [ firstDoF + 3 ];
    } else if ( dim > 2 && code.compare("roty") == 0 ) {
        result.resize(1);
        result [ 0 ] = DoFs [ firstDoF + 4 ];
    } else if ( dim > 2 && code.compare("rotz") == 0 ) {
        result.resize(1);
        result [ 0 ] = DoFs [ firstDoF + 5 ];
    } else {
        MechNode :: giveDoFBasedValues(code, DoFs, result);
    }
};

//////////////////////////////////////////////////////////
unsigned Particle :: giveOrderOfForceCode(string code) const {
    if ( dim == 3 && code.compare("mx") == 0 ) {
        return 3;
    } else if ( dim == 3 && code.compare("my") == 0 ) {
        return 4;
    } else if ( code.compare("mz") == 0 ) {
        if ( dim == 2 ) {
            return 2;
        } else {
            return 5;
        }
    } else {
        return MechNode :: giveOrderOfForceCode(code);
    }
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
void CoupledParticle :: giveDoFBasedValues(string code, const Vector &DoFs, Vector &result) const {
    if ( code.compare("pressure") == 0 ) {
        result.resize(1);
        result [ 0 ] = DoFs [ firstDoF + 6 ];
    } else {
        Particle :: giveDoFBasedValues(code, DoFs, result);
    }
};

//////////////////////////////////////////////////////////
unsigned CoupledParticle :: giveOrderOfForceCode(string code) const {
    if ( code.compare("flux") == 0 ) {
        if ( dim == 2 ) {
            return 3;
        } else {
            return 6;
        }
    } else {
        return Particle :: giveOrderOfForceCode(code);
    }
}
