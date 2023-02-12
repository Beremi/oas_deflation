#include "node.h"
#include "boundary_condition.h"
#include "element_discrete.h"
#include "simplex.h"

using namespace std;

//////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////
// BASIC NODE - ONLY MASTER CLASS

Node :: Node(unsigned dimension) {
    name = "generic node";
    physicalFields.resize(4, false); //mechanical, transport, thermal, humidity
    physicalFieldsDoFNum.resize(4, 0);
    point = Point(0, 0, 0);
    dim = dimension;
    simplex = nullptr;
}

//////////////////////////////////////////////////////////
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

//////////////////////////////////////////////////////////
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
bool Node :: giveDoFBasedValues(string code, const Vector &DoFs, Vector &result) const {
    if ( code.compare("id") == 0 ) {
        result.resize(1);
        result [ 0 ] = id;
        return true;
    }

    //mechanics
    unsigned initialDoF = firstDoF;
    if ( physicalFields [ 0 ] ) {
        if ( code.compare("displacements") == 0  || code.compare("displacement") == 0  ) {
            result.resize(3); //always 3 components to show it in VTK
            unsigned i;
            for ( i = 0; i < dim; i++ ) {
                result [ i ] = DoFs [ initialDoF + i ];
            }
            for ( ; i < 3; i++ ) {
                result [ i ] = 0;
            }
            return true;
        } else if ( code.compare("ux") == 0 ) {
            result.resize(1);
            result [ 0 ] = DoFs [ initialDoF ];
            return true;
        } else if ( dim > 1 && code.compare("uy") == 0 ) {
            result.resize(1);
            result [ 0 ] = DoFs [ initialDoF + 1 ];
            return true;
        } else if ( dim > 2 && code.compare("uz") == 0 ) {
            result.resize(1);
            result [ 0 ] = DoFs [ initialDoF + 2 ];
            return true;
        }
    }

    //transport
    initialDoF += physicalFieldsDoFNum [ 0 ];
    if ( physicalFields [ 1 ] ) {
        if ( code.compare("pressure") == 0 ) {
            result.resize(1);
            result [ 0 ] = DoFs [ initialDoF ];
            return true;
        }
    }

    //thermal
    initialDoF += physicalFieldsDoFNum [ 1 ];
    if ( physicalFields [ 2 ] ) {
        if ( code.compare("temperature") == 0 ) {
            result.resize(1);
            result [ 0 ] = DoFs [ initialDoF ];
            return true;
        }
    }

    //humidity
    initialDoF += physicalFieldsDoFNum [ 2 ];
    if ( physicalFields [ 3 ] ) {
        if ( code.compare("humidity") == 0 ) {
            result.resize(1);
            result [ 0 ] = DoFs [ initialDoF ];
            return true;
        }
    }

    char *pEnd;
    unsigned converted = strtol(code.c_str(), & pEnd, 10);
    if ( !* pEnd ) {
        if ( converted < nDoFs ) {
            result.resize(1);
            result [ 0 ] = DoFs [ firstDoF + converted ];
            return true;
        } else  {
            cerr << name << "Error: Requested DoFid exceeded number of DoFs: " << code << endl;
            exit(1);
        }
    }

    result.resize(0);
    return false;
};

//////////////////////////////////////////////////////////
unsigned Node :: giveOrderOfEnergyConjugateCode(string code) const {
    //mechanics
    unsigned initialDoF = firstDoF;
    if ( physicalFields [ 0 ] ) {
        if ( code.compare("fx") == 0 ) {
            return initialDoF;
        } else if ( dim > 1 && code.compare("fy") == 0 ) {
            return initialDoF + 1;
        } else if ( dim > 2 && code.compare("fz") == 0 ) {
            return initialDoF + 2;
        }
    }

    //transport
    initialDoF += physicalFieldsDoFNum [ 0 ];
    if ( physicalFields [ 1 ] ) {
        if ( code.compare("fluid_flux") == 0 ) {
            return initialDoF;
        }
    }

    //thermal
    initialDoF += physicalFieldsDoFNum [ 1 ];
    if ( physicalFields [ 1 ] ) {
        if ( code.compare("heat_flux") == 0 ) {
            return initialDoF;
        }
    }

    //humidity
    initialDoF += physicalFieldsDoFNum [ 2 ];
    if ( physicalFields [ 1 ] ) {
        if ( code.compare("humidity_flux") == 0 ) {
            return initialDoF;
        }
    }

    char *pEnd;
    unsigned converted = strtol(code.c_str(), & pEnd, 10);
    if ( !* pEnd ) {
        if ( converted < nDoFs ) {
            return firstDoF + converted;
        } else  {
            cerr << name << "Error: Requested DoFid exceeded number of DoFs: " << code << endl;
            exit(1);
        }
    }

    cerr << "Required code not found: " << code << endl;
    exit(1);
};

//////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////
// MECHANICAL NODE - translational DoFs
void Node :: init() {
    nDoFs = std :: reduce( physicalFieldsDoFNum.begin(), physicalFieldsDoFNum.end() );
};

//////////////////////////////////////////////////////////
unsigned Node :: giveRelativeDoFPhysicalFieldNum(unsigned k) const {
    unsigned h = 0;
    for ( size_t p = 0; p < physicalFields.size(); p++ ) {
        h += physicalFieldsDoFNum [ p ];
        if ( k < h ) {
            return p;
        }
    }
    cerr << name << ": DoF number exceeded number of DoFs" << endl;
    exit(1);
    return 0;
};

//////////////////////////////////////////////////////////
unsigned Node :: giveAbsoluteDoFPhysicalFieldNum(unsigned k) const {
    return giveRelativeDoFPhysicalFieldNum(k - firstDoF);
};

//////////////////////////////////////////////////////////
vector< unsigned >Node :: givePhysicalFieldNumForAllDoFs() const {
    vector< unsigned >pf(nDoFs);
    unsigned h = 0;
    for ( size_t p = 0; p < physicalFields.size(); p++ ) {
        for ( unsigned k = 0; k < physicalFieldsDoFNum [ p ]; k++ ) {
            pf [ h + k ] = p;
        }
        h += physicalFieldsDoFNum [ p ];
    }
    return pf;
};

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
void FreeDoF :: readFromLine(istringstream &iss) {
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
        cerr << name << " Error: reading failed!" << endl;
        exit(1);
    }
}

//////////////////////////////////////////////////////////
FreeDoF :: FreeDoF(unsigned dimension, unsigned numDoFs) : Node(dimension) {
    nDoFs = numDoFs;
    name = "FreeDoF";
}

//////////////////////////////////////////////////////////
void FreeDoF :: init() {
    Node :: init();
}

//////////////////////////////////////////////////////////
MechDoF :: MechDoF(unsigned dimension, unsigned numDoFs) : FreeDoF(dimension, numDoFs) {
    name = "MechDoF";
    physicalFields [ 0 ] = true;
}

//////////////////////////////////////////////////////////
void MechDoF :: init() {
    physicalFieldsDoFNum [ 0 ] = nDoFs;
    FreeDoF :: init();
}

//////////////////////////////////////////////////////////
TrsDoF :: TrsDoF(unsigned dimension, unsigned numDoFs) : FreeDoF(dimension, numDoFs) {
    name = "TrsDoF";
    physicalFields [ 1 ] = true;
}

//////////////////////////////////////////////////////////
void TrsDoF :: init() {
    physicalFieldsDoFNum [ 1 ] = nDoFs;
    FreeDoF :: init();
}
//////////////////////////////////////////////////////////
TempDoF :: TempDoF(unsigned dimension, unsigned numDoFs) : FreeDoF(dimension, numDoFs) {
    name = "TempDoF";
    physicalFields [ 2 ] = true;
}

//////////////////////////////////////////////////////////
void TempDoF :: init() {
    physicalFieldsDoFNum [ 2 ] = nDoFs;
    FreeDoF :: init();
}
//////////////////////////////////////////////////////////
HumidityDoF :: HumidityDoF(unsigned dimension, unsigned numDoFs) : FreeDoF(dimension, numDoFs) {
    physicalFields [ 3 ] = true;
    name = "HumdityDoF";
}

//////////////////////////////////////////////////////////
void HumidityDoF :: init() {
    physicalFieldsDoFNum [ 3 ] = nDoFs;
    FreeDoF :: init();
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
Vector Particle :: calculateRigidBodyMotionVector(const Point *x, const Vector &DoFs) const {
    unsigned DofsPerNode = ( dim - 1 ) * 3;
    Vector u = Vector :: Zero(DofsPerNode);
    for ( unsigned i = 0; i < DofsPerNode; i++ ) {
        u [ i ]  = DoFs [ firstDoF + i ];
    }
    return giveRigidBodyMotionMatrix(x) * u;
}

//////////////////////////////////////////////////////////
Point Particle :: calculateRigidBodyMotionPoint(const Point *x, const Vector &DoFs) const {
    Vector u = calculateRigidBodyMotionVector(x, DoFs);
    return Point(u [ 0 ], u [ 1 ], dim > 2 ? u [ 2 ] : 0);
}

//////////////////////////////////////////////////////////
Matrix Particle :: giveRigidBodyMotionMatrix(const Point *x) const {
    Matrix A = Matrix :: Zero( dim, 3 * ( dim - 1 ) );
    if ( dim == 3 ) {
        A(0, 0) = A(1, 1) = A(2, 2) = 1;
        A(1, 3) = point.z() - x->z();
        A(0, 4) = -A(1, 3);
        A(2, 3) = x->y() - point.y();
        A(0, 5) = -A(2, 3);
        A(2, 4) = point.x() - x->x();
        A(1, 5) = -A(2, 4);
    } else if ( dim == 2 ) {
        A(0, 0) = A(1, 1) = 1;
        A(0, 2) = point.y() - x->y();
        A(1, 2) = x->x() - point.x();
    } else {
        cerr << "Error - Particle: dimension " << dim << "not implemented" << endl;
        exit(EXIT_FAILURE);
    }
    return A;
}

//////////////////////////////////////////////////////////
bool Particle :: giveDoFBasedValues(string code, const Vector &DoFs, Vector &result) const {
    if ( code.compare("rotation") == 0  || code.compare("rotations") == 0  ) {
        result.resize(2 * dim - 3);
        for ( unsigned i = 0; i < 2 * dim - 3; i++ ) {
            result [ i ] = DoFs [ firstDoF + dim + i ];
        }
        return true;
    } else if ( dim == 2 && code.compare("rotz") == 0 ) {
        result.resize(1);
        result [ 0 ] = DoFs [ firstDoF + 2 ];
        return true;
    } else if ( dim == 3 && code.compare("rotx") == 0 ) {
        result.resize(1);
        result [ 0 ] = DoFs [ firstDoF + 3 ];
        return true;
    } else if ( dim == 3 && code.compare("roty") == 0 ) {
        result.resize(1);
        result [ 0 ] = DoFs [ firstDoF + 4 ];
        return true;
    } else if ( dim == 3 && code.compare("rotz") == 0 ) {
        result.resize(1);
        result [ 0 ] = DoFs [ firstDoF + 5 ];
        return true;
    }
    return Node :: giveDoFBasedValues(code, DoFs, result);
}

//////////////////////////////////////////////////////////
unsigned Particle :: giveOrderOfEnergyConjugateCode(string code) const {
    //mechanics
    if ( dim == 3 && code.compare("mx") == 0 ) {
        return firstDoF + 3;
    } else if ( dim == 3 && code.compare("my") == 0 ) {
        return firstDoF + 4;
    } else if ( code.compare("mz") == 0 ) {
        if ( dim == 2 ) {
            return firstDoF + 2;
        } else {
            return firstDoF + 5;
        }
    }

    return Node :: giveOrderOfEnergyConjugateCode(code);
}

//////////////////////////////////////////////////////////
std :: string Particle :: giveLineToSave() const {
    std :: ostringstream out;
    out.precision(10);
    out << Node :: giveLineToSave() << "\t";
    out << std :: fixed << this->r;
    return out.str();
};
