#include "element.h"
#include "element_container.h"
#include "boundary_condition.h"

using namespace std;

//////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////
// BASIC ELEMENT - MASTER CLASS
Element :: ~Element() {
    for ( std :: vector< MaterialStatus * > :: iterator e = stats.begin(); e != stats.end(); ++e ) {
        delete * e;
    }
    delete shafunc;
    delete inttype;
}

// std :: string Element :: giveLineToSave(NodeContainer * nodes_all) const {
//   std :: string str = this->giveName() + "\t";
//   if ( this->nodes.size() > 2 ){
//     // for polygonal elems etc
//     str += to_string(this->nodes.size()) + "\t";
//   }
//   for ( auto const &node : this->nodes ){
//     str += to_string(nodes_all->giveNodeId(node)) + "\t";
//   }
//   return str;
// }

//////////////////////////////////////////////////////////
void Element :: initIntegration() {
    shafunc->init(nodes);
    inttype->init();
};

//////////////////////////////////////////////////////////
void Element :: setIntegrationPointsAndWeights() {
    stats.resize( inttype->giveNumIP() );
    for ( unsigned k = 0; k < inttype->giveNumIP(); k++ ) {
        stats [ k ] = mat->giveNewMaterialStatus(this, k);
        inttype->setIPWeight( k, inttype->giveIPWeight(k) * shafunc->giveJacobian(inttype->giveIPLocationPointer(k) ) );
    }
};


//////////////////////////////////////////////////////////
void Element :: readFromLine(std :: istringstream &iss, NodeContainer *fullnodes, MaterialContainer *fullmatrs) {
    unsigned num;
    nodes.resize(numOfNodes);
    for ( unsigned k = 0; k < numOfNodes; k++ ) {
        iss >> num;
        nodes [ k ] = fullnodes->giveNode(num);
    }
    iss >> num;
    mat = fullmatrs->giveMaterial(num);
}

//////////////////////////////////////////////////////////
void Element :: init() {
    //delete possible previous statuses
    for ( std :: vector< MaterialStatus * > :: iterator e = stats.begin(); e != stats.end(); ++e ) {
        delete * e;
    }

    unsigned totalDoFs = 0;
    for ( std :: vector< Node * > :: const_iterator n = nodes.begin(); n != nodes.end(); ++n ) {
        totalDoFs += ( * n )->giveNumberOfDoFs();
    }

    initIntegration();
    setIntegrationPointsAndWeights();

    DoFids.resize(totalDoFs);
    unsigned i = 0;
    unsigned k;
    for ( std :: vector< Node * > :: const_iterator n = nodes.begin(); n != nodes.end(); ++n ) {
        k = ( * n )->giveStartingDoF();
        for ( unsigned s = 0; s < ( * n )->giveNumberOfDoFs(); s++, i++ ) {
            DoFids [ i ] = k + s;
        }
    }
    outDoFs = totalDoFs; //basic elems will alway have input = output

    Bs.resize( inttype->giveNumIP() );
    Hs.resize( inttype->giveNumIP() );
    for ( k = 0; k < inttype->giveNumIP(); k++ ) {
        Bs [ k ] = giveBMatrix(k);
        Hs [ k ] = giveHMatrix(k);
    }

    //set stress and strain vectors at integration points
    for ( k = 0; k < inttype->giveNumIP(); k++ ) {
        stats [ k ]->initializeStressAndStrainVector(Bs [ k ].rows() );
    }
}

//////////////////////////////////////////////////////////
std :: vector< unsigned >Element :: giveDoFsInDirection(unsigned dir) const {
    std :: vector< unsigned >DoFinDir( nodes.size() );
    for ( unsigned i = 0; i < nodes.size(); i++ ) {
        DoFinDir [ i ] = nodes [ i ]->giveStartingDoF() + dir;
    }
    return DoFinDir;
}

//////////////////////////////////////////////////////////
void Element :: initMaterialStatuses() {
    unsigned num = 0;
    for ( std :: vector< MaterialStatus * > :: iterator m = stats.begin(); m != stats.end(); ++m, num++ ) {
        ( * m )->init();
    }
}

//////////////////////////////////////////////////////////
void Element :: updateMaterialStatuses() {
    for ( std :: vector< MaterialStatus * > :: iterator m = stats.begin(); m != stats.end(); ++m ) {
        ( * m )->update();
    }
}

//////////////////////////////////////////////////////////
void Element :: resetMaterialStatuses() {
    for ( std :: vector< MaterialStatus * > :: iterator m = stats.begin(); m != stats.end(); ++m ) {
        ( * m )->resetTemporaryVariables();
    }
}


//////////////////////////////////////////////////////////
void Element :: giveIPValues(std :: string code, unsigned ipnum, Vector &result) const {
    if ( ipnum >= inttype->giveNumIP() ) {
        std :: cerr << name <<  " Error: intergration point number " << ipnum << " exceeds number of integration points" << std :: endl;
        exit(1);
    }
    if ( code.compare("x") == 0 ) {
        result.resize(1);
        result [ 0 ] = inttype->giveIPLocationPointer(ipnum)->x();
    } else if ( code.compare("y") == 0 ) {
        result.resize(1);
        result [ 0 ] = inttype->giveIPLocationPointer(ipnum)->y();
    } else if ( code.compare("z") == 0 ) {
        result.resize(1);
        result [ 0 ] = inttype->giveIPLocationPointer(ipnum)->z();
    } else {        
        stats [ ipnum ]->giveValues(code, result);
    }
};

//////////////////////////////////////////////////////////
MaterialStatus * Element :: giveMatStatus(unsigned ipnum){
    if ( ipnum >= inttype->giveNumIP() ) {
        std :: cerr << name <<  " Error: intergration point number " << ipnum << " exceeds number of integration points" << std :: endl;
        exit(1);
    }        
    return stats[ipnum];
}

//////////////////////////////////////////////////////////
void Element :: giveValues(std :: string code, Vector &result) const {
    if ( code.compare("id") == 0 ) {
        result.resize(1);
        result [ 0 ] = idx;
    } else if ( code.compare("strain_energy") == 0 ) {
        Element :: giveValues("strain_energy_density", result);
        result *= giveVolume() * ndim; //ndim because of discrete elements
        //cout << volume << endl; exit(1);
    } else if ( code.compare("total_energy") == 0 ) {
        Element :: giveValues("total_energy_density", result);
        result *= giveVolume() * ndim; //ndim because of discrete elements
    } else if ( code.compare("dissipated_energy") == 0 ) {
        Element :: giveValues("dissipated_energy_density", result);
        result *= giveVolume() * ndim; //ndim because of discrete elements
    } else if ( code.compare("dissipated_energy_inc") == 0 ) {
        Element :: giveValues("dissipated_energy_density_inc", result);
        result *= giveVolume() * ndim; //ndim because of discrete elements
    } else if ( code.compare("kinetic_energy") == 0 ) {
        result.resize(1);
        result [ 0 ] = 0;
    } else if ( code.compare("material_ID") == 0 ) {
        result.resize(1);
        result [ 0 ] = mat->giveId();
    } else {  //TODO: should be weighted average
        //average values from IP
        if ( inttype->giveNumIP() > 0 ) {
            stats [ 0 ]->giveValues(code, result);
            Vector res2(result.size() );
            for ( unsigned i = 1; i < inttype->giveNumIP(); i++ ) {
                stats [ i ]->giveValues(code, res2);
                result += res2;
            }
            result /= inttype->giveNumIP();
        } else {
            result.resize(0);
        }
    }
}

//////////////////////////////////////////////////////////
Matrix Element :: giveStiffnessMatrix(std :: string matrixType) const {
    unsigned nDoFs = DoFids.size();
    Matrix K = Matrix :: Zero(nDoFs, nDoFs);
    Matrix D;
    for ( unsigned i = 0; i < inttype->giveNumIP(); i++ ) {
        D = stats [ i ]->giveStiffnessTensor(matrixType);
        K += Bs [ i ].transpose() * D * ( Bs [ i ] * inttype->giveIPWeight(i) );
    }
    return K;
}

//////////////////////////////////////////////////////////
Vector Element :: giveStrain(unsigned i, const Vector &DoFs) {
    return Bs [ i ] * DoFs;
}


//////////////////////////////////////////////////////////
Vector Element :: giveInternalForces(const Vector &DoFs, bool frozen, double timeStep) {
    Vector intF = Vector :: Zero( DoFids.size() );
    Vector stress;
    for ( unsigned i = 0; i < inttype->giveNumIP(); i++ ) {
        if ( frozen ) {
            stress = stats [ i ]->giveStressWithFrozenIntVars(giveStrain(i, DoFs), timeStep);  //frozen internal variables
        } else {
            stress = stats [ i ]->giveStress(giveStrain(i, DoFs), timeStep); //full evaluation of stress including change of state variables
        }
        intF  += Bs [ i ].transpose() * (  stress * inttype->giveIPWeight(i) );
    }
    // Remove cout
    // cout << name << endl;
    // cout << Bs [ 0 ].transpose() << endl;
    // cout << intF << endl;

    //add internal sources
    if ( mat->isProducingInternalSources() ) {
        Vector intS = integrateInternalSources();
        for ( unsigned i = 0; i < intF.size(); i++ ) {
            intF [ i ] += intS [ i ];
        }
    }
    return intF;
}

//////////////////////////////////////////////////////////
Vector Element :: integrateInternalSources() {
    Vector intS = Vector :: Zero( DoFids.size() );
    Vector intmats;
    for ( unsigned i = 0; i < inttype->giveNumIP(); i++ ) {
        intmats = stats [ i ]->giveInternalSource();
        intS += Hs [ i ].transpose() * ( intmats * inttype->giveIPWeight(i) );
    }

    return intS;
}


//////////////////////////////////////////////////////////
Matrix Element :: giveDampingMatrix() const {
    unsigned nDoFs = DoFids.size();
    Matrix M = Matrix :: Zero(nDoFs, nDoFs);
    Matrix c;
    for ( unsigned i = 0; i < inttype->giveNumIP(); i++ ) {
        c = stats [ i ]->giveDampingTensor();
        M += Hs [ i ].transpose() * ( c * inttype->giveIPWeight(i) ) * Hs [ i ];
    }
    return M;
}

//////////////////////////////////////////////////////////
Matrix Element :: giveMassMatrix() const {
    unsigned nDoFs = DoFids.size();
    Matrix M = Matrix :: Zero(nDoFs, nDoFs);
    Matrix c;
    for ( unsigned i = 0; i < inttype->giveNumIP(); i++ ) {
        c = stats [ i ]->giveMassTensor();
        M += Hs [ i ].transpose() * ( c * inttype->giveIPWeight(i) ) * Hs [ i ];
    }
    return M;
}

//////////////////////////////////////////////////////////
Vector Element :: integrateLoad(BodyLoad *vl, double time) const {
    unsigned nDoFs = DoFids.size();
    Vector load = Vector :: Zero(nDoFs);
    double fvalue;
    unsigned dir = vl->giveDirection();
    for ( unsigned i = 0; i < inttype->giveNumIP(); i++ ) {
        fvalue = vl->giveValue(inttype->giveIPLocationPointer(i), time);
        for ( unsigned j = 0; j < nDoFs; j++ ) {
            load [ j ] += Hs [ i ](dir, j) * fvalue * inttype->giveIPWeight(i);
        }
    }
    return load;
}

//////////////////////////////////////////////////////////
void Element :: changeMaterial(Material *newmat) {
    this->mat = newmat;
    for ( unsigned p = 0; p < stats.size(); p++ ) {
        if ( stats [ p ] != nullptr ) {
            delete stats [ p ];
        }
        stats [ p ] = this->mat->giveNewMaterialStatus(this, p);
    }
    this->initMaterialStatuses();
}

//////////////////////////////////////////////////////////
bool Element :: giveGlobalCoords(Point *x, const Point *xn) const {
    Vector phi = Vector :: Zero( nodes.size() );
    shafunc->giveShapeF(xn, phi);
    * x = Point(0, 0, 0);
    for ( unsigned n = 0; n < nodes.size(); n++ ) {
        * x += nodes [ n ]->givePoint() * phi [ n ];
    }
    return true;
}


//////////////////////////////////////////////////////////
bool Element :: isPointInside(Point *xn, const Point *x) const {
    //initial screening
    Point maxc(-1e10, -1e10, -1e10);
    Point minc(1e10, 1e10, 1e10);
    Point *p;
    for ( auto &n: nodes ) {
        p = n->givePointPointer();
        for ( unsigned c = 0; c < ndim; c++ ) {
            maxc(c) = std :: max( maxc(c), ( * p )(c) );
            minc(c) = std :: min( minc(c), ( * p )(c) );
        }
    }
    for ( unsigned c = 0; c < ndim; c++ ) {
        if ( ( * x )(c) > maxc(c) || ( * x )(c) < minc(c) ) {
            return false;
        }
    }

    //find natural coordinates
    //initial estimation
    Point center(0, 0, 0);
    * xn = Point(0, 0, 0);
    Point aux(0, 0, 0), diff(0, 0, 0), diffC(0, 0, 0);
    Point size = maxc - minc;
    giveGlobalCoords(& center, xn);

    //initial estimation
    aux = ( ( * x ) - center ) * 2.;
    for ( unsigned c = 0; c < ndim; c++ ) {
        ( * xn )(c) = aux(c) / size(c);
    }
    giveGlobalCoords(& aux, xn);

    int i = 0;
    int max_i = 500;
    double maxerror = 1;
    diff = aux - ( * x );
    while ( maxerror > 1e-4 && i < max_i ) {
        for ( unsigned c = 0; c < ndim; c++ ) {
            if ( abs( diff(c) / size(c) ) < 1e-8 ) {
                continue;
            }
            if ( diffC(c) > 1e-16 ) {
                ( * xn )(c) = ( * xn )(c) - ( * xn )(c) * diff(c) / diffC(c);
            } else {
                ( * xn )(c) = ( * xn )(c) - 2. * diff(c) / size(c);
            }
        }

        giveGlobalCoords(& aux, xn);
        diff = aux - ( * x );
        diffC = aux - center;
        maxerror = 0.;
        for ( unsigned c = 0; c < ndim; c++ ) {
            maxerror = std :: max( maxerror, abs( diff(c) / size(c) ) );
        }
        i++;
    }
    if ( i == max_i ) {
        return false;
    }

    //check natural coordinates are inside limits
    //TODO: works only for brick and quadrilateral
    for ( unsigned c = 0; c < ndim; c++ ) {
        if ( abs( ( * xn )(c) ) > 1. ) {
            return false;
        }
    }
    return true;
}

//////////////////////////////////////////////////////////
Vector Element :: giveElemDoFsFromFullDoFs(const Vector &FullDoFs) const {
    Vector elemDoFs = Vector :: Zero( DoFids.size() );
    for ( unsigned i = 0; i < DoFids.size(); i++ ) {
        elemDoFs [ i ] = FullDoFs [ DoFids [ i ] ];
    }
    return elemDoFs;
}

//////////////////////////////////////////////////////////
void Element :: extrapolateIPValuesToNodes(std :: string code, vector< Vector > &result, Vector &weights) const {
    Vector phi = Vector :: Zero( nodes.size() );
    Vector res = Vector :: Zero( nodes.size() );
    double jacobian;
    Vector ipres;
    Matrix M = Matrix :: Zero( nodes.size(), nodes.size() );

    if ( inttype->giveNumIP() == 0 ) {
        std :: cerr << "Error in function extrapolateIPValuesToNodes: zero number of integration points" << std :: endl;
        exit(1);
    }

    giveIPValues(code, 0, ipres);
    unsigned reslen = ipres.size();
    result.resize(reslen);

    std :: vector< Vector >rhs(reslen);
    for ( unsigned h = 0; h < reslen; h++ ) {
        rhs [ h ] = Vector :: Zero(nodes.size() );
        result [ h ] = Vector :: Zero(nodes.size() );
    }
    weights.resize( nodes.size() );
    weights.setOnes(); //for(auto &h: weights) h = 1;


    for ( unsigned i = 0; i < inttype->giveNumIP(); i++ ) {
        shafunc->giveShapeF(inttype->giveIPLocationPointer(i), phi);
        jacobian = shafunc->giveJacobian(inttype->giveIPLocationPointer(i) );
        giveIPValues(code, i, ipres);
        for ( unsigned k = 0; k < nodes.size(); k++ ) {
            for ( unsigned h = 0; h < reslen; h++ ) {
                rhs [ h ] [ k ] += phi [ k ] * jacobian * ipres [ h ];
            }
            for ( unsigned l = 0; l < nodes.size(); l++ ) {
                M(k, l) += phi [ k ] * phi [ l ] * jacobian;
            }
        }
    }

    Matrix M_inv = M.inverse();

    for ( unsigned h = 0; h < reslen; h++ ) {
        result [ h ] = M_inv * rhs [ h ];
    }
}

//////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////
// MATERIAL TEST ELEMENT - only one material point and virtual loading through prescribed strains
//////////////////////////////////////////////////////////
MaterialTestElement :: MaterialTestElement(unsigned dim) : Element(dim) {
    ndim = dim;
    numOfNodes = 1;
    name = "MaterialTestElement";
    vtk_cell_type = 1;
    shafunc = new Linear1DLineShapeF();
    inttype = new IntegrDiscrete1();
}


//////////////////////////////////////////////////////////
void MaterialTestElement :: setIntegrationPointsAndWeights() {
    stats.resize(1);
    inttype->setIPLocation(0, Point(0., 0., 0.) );
    inttype->setIPWeight(0, 1);
    stats [ 0 ] = mat->giveNewMaterialStatus(this, 0);
}

//////////////////////////////////////////////////////////
Matrix MaterialTestElement :: giveBMatrix(const Point *x) const {
    ( void ) x;
    Matrix B = Matrix :: Identity(DoFids.size(), DoFids.size() );
    return B;
}

//////////////////////////////////////////////////////////
Matrix MaterialTestElement :: giveHMatrix(const Point *x) const {
    ( void ) x;
    unsigned numOfIntSources = 7;  //TODO: THIS IS WRONG, NEEDS TO BE TREATED AUTOMATICALLY
    Matrix H = Matrix :: Zero(numOfIntSources, DoFids.size() );
    return H;
}
