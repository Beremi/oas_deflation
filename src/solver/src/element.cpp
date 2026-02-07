#include "element.h"
#include "element_container.h"
#include "boundary_condition.h"
#include <Eigen/Dense>

using namespace std;

//////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////
// BASIC ELEMENT - MASTER CLASS

//////////////////////////////////////////////////////////
Element :: Element(unsigned dim) {
    name = "basic element";
    solution_order = 0;
    volume = 0;
    ndim = dim;
    physicalFields.resize(4, false); //mechanical, transport, thermal, humidity
    areIPLocsInNaturalCoords = true;
}

//////////////////////////////////////////////////////////
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

    volume  = 0;
    for ( k = 0; k < inttype->giveNumIP(); k++ ) {
        volume += giveIPVolume(k);
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
void Element :: removeEigenStrain() {
    for ( auto &s: stats ) {
        s->removeEigenStrain();
    }
}


//////////////////////////////////////////////////////////
void Element :: giveIPValues(std :: string code, unsigned ipnum, Vector &result) const {
    if ( ipnum >= inttype->giveNumIP() ) {
        std :: cerr << name <<  " Error: intergration point number " << ipnum << " exceeds number of integration points" << std :: endl;
        exit(1);
    }
    if ( code.compare("location") == 0 ) {
        result.resize(ndim);
        fill(result.begin(), result.end(), 0);
        if ( areIPLocsInNaturalCoords ) {
            Vector phi(this->giveNumOfNodes() );
            shafunc->giveShapeF(inttype->giveIPLocationPointer(ipnum), phi);
            for ( unsigned n = 0; n < this->giveNumOfNodes(); n++ ) {
                Point *nn = nodes [ n ]->givePointPointer();
                for ( unsigned k = 0; k < ndim; k++ ) {
                    result [ k ] += ( * nn ) [ k ] * phi [ n ];
                }
            }
        } else {
            for ( unsigned k = 0; k < ndim; k++ ) {
                result [ k ] = ( * inttype->giveIPLocationPointer(ipnum) ) [ k ];
            }
        }
    } else if ( code.compare("weight") == 0 ) {
        result.resize(1);
        result [ 0 ] = inttype->giveIPWeight(ipnum);
    } else if ( code.compare("id") == 0 || code.compare("element_id") == 0 ) {
        result.resize(1);
        result [ 0 ] = idx;
    } else if ( code.compare("x") == 0 ) {
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
MaterialStatus *Element :: giveMatStatus(unsigned ipnum) {
    if ( ipnum >= inttype->giveNumIP() ) {
        std :: cerr << name <<  " Error: intergration point number " << ipnum << " exceeds number of integration points" << std :: endl;
        exit(1);
    }
    return stats [ ipnum ];
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
    } else if ( code.compare("material_ID") == 0 or code.compare("materialID") == 0 ) {
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
    /*
     * // This check works only if flag "-ffast-math" is removed from CMake
     * if (std::isnan(K.template maxCoeff<Eigen::PropagateNaN>())){
     *  cout << "ELEMENT " << idx << " ("<< name <<") has NaN in stiffness matrix" << endl;
     * }
     */
    return K;
}

//////////////////////////////////////////////////////////
void Element :: evaluateStrains(const Vector &DoFs) {
    for ( unsigned i = 0; i < inttype->giveNumIP(); i++ ) {
        stats [ i ]->setTotalTempStrain(Bs [ i ] * DoFs);
    }
}

//////////////////////////////////////////////////////////
void Element :: evaluateStresses(bool frozen, double timeStep) {
    if ( frozen ) {
        evaluateStressesWithFrozenIntVars(timeStep);
    } else {
        evaluateStresses(timeStep);
    }
}

//////////////////////////////////////////////////////////
void Element :: evaluateStresses(double timeStep) {
    for ( auto &s: stats) {
        s->computeStress(timeStep);
    }
}

//////////////////////////////////////////////////////////
void Element :: evaluateStressesWithFrozenIntVars(double timeStep) {
    for ( auto &s: stats) {
        s->computeStressWithFrozenIntVars(timeStep);
    }
}

//////////////////////////////////////////////////////////
Vector Element :: giveInternalForces() {
    Vector intF = Vector :: Zero( DoFids.size() );
    Vector stress;
    for ( unsigned i = 0; i < inttype->giveNumIP(); i++ ) {
        stress = stats [ i ]->giveTempStress();  //frozen internal variables
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
double Element :: giveKineticEnergy(const Vector &velocity) const {
    return 0.5 * ( velocity.dot(massM * velocity) );
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
void Element :: computeDampingMatrix() {
    unsigned nDoFs = DoFids.size();
    dampC = Matrix :: Zero(nDoFs, nDoFs);
    Matrix c;
    for ( unsigned i = 0; i < inttype->giveNumIP(); i++ ) {
        c = stats [ i ]->giveDampingTensor();
        dampC += Hs [ i ].transpose() * ( c * inttype->giveIPWeight(i) ) * Hs [ i ];
    }
}

//////////////////////////////////////////////////////////
Matrix Element :: giveDampingMatrix() {
    if ( mat->requiresDampingsMatrixUpdate()  || dampC.rows() == 0 ) {
        computeDampingMatrix();
    }
    return dampC;
}

//////////////////////////////////////////////////////////
void Element :: computeMassMatrix() {
    unsigned nDoFs = DoFids.size();
    massM = Matrix :: Zero(nDoFs, nDoFs);
    Matrix c;
    for ( unsigned i = 0; i < inttype->giveNumIP(); i++ ) {
        c = stats [ i ]->giveMassTensor();
        massM += Hs [ i ].transpose() * ( c * inttype->giveIPWeight(i) ) * Hs [ i ];
    }
}

//////////////////////////////////////////////////////////
Matrix Element :: giveMassMatrix() {
    if ( mat->requiresMassMatrixUpdate() || massM.rows() == 0 ) {
        computeMassMatrix();
    }
    return massM;
}

//////////////////////////////////////////////////////////
Matrix Element :: giveLumpedMassMatrix() {
    if ( mat->requiresMassMatrixUpdate() || massM.rows() == 0 ) {
        computeMassMatrix();
    }

    //cannot mix rotations and translations
    unsigned nDoFs = DoFids.size();
    vector< bool >indicateRot(nDoFs);
    unsigned k = 0;
    for ( auto &nn:nodes ) {
        for ( unsigned p = 0; p < nn->giveNumberOfDoFs(); p++ ) {
            if ( physicalFields [ k ] == 0 && p > ndim ) {
                indicateRot [ k ] = true;
            } else {
                indicateRot [ k ] = false;
            }
            k += 1;
        }
    }

    Matrix lumpedMassM = Matrix :: Zero(nDoFs, nDoFs);
    for ( unsigned i = 0; i < nDoFs; i++ ) {
        lumpedMassM(i, i) += massM(i, i);
        for ( unsigned j = i + 1; j < nDoFs; j++ ) {
            if ( physicalFields [ i ] == physicalFields [ j ] ) {
                //if (physicalFields [ i ]!=0 || indicateRot[i]==indicateRot[j]) lumpedMassM[i] += massM.coeff(i,j);
                if ( physicalFields [ i ] != 0 || ( !indicateRot [ i ] && !indicateRot [ j ] ) ) {
                    lumpedMassM(i, i) += massM(i, j);                                                                //rotations only at the diagonal
                }
            }
        }
    }
    return lumpedMassM;
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
Point Element :: findNaturalCoords(const Point *x) const {
    Point natcoords;
    if ( shafunc->isInNaturalCoords() ) {
        natcoords = Point(0, 0, 0);
        Point testglocoords, glocoords, testnatcoords;
        Vector errcoordsvec(ndim);
        Matrix grad = Matrix :: Zero(ndim, ndim);
        Vector natcoordsvec;
        giveGlobalCoords(& glocoords, & natcoords);
        for (unsigned i = 0; i < ndim; i++) {
            errcoordsvec [ i ] = glocoords [ i ] - ( * x ) [ i ];
        }
        double err = errcoordsvec.norm();
        double step = 1e-3;
        unsigned maxit = 10;
        unsigned it = 0;
        double xnorm = max(x->norm(), 1e-5);
        while ( err / xnorm > 1e-4 && it < maxit ) {
            testnatcoords = natcoords;
            for ( unsigned dim = 0; dim < ndim; dim++ ) {
                testnatcoords [ dim ] += step;
                giveGlobalCoords(& testglocoords, & testnatcoords);
                for (unsigned i = 0; i < ndim; i++) {
                    grad(i, dim) = ( testglocoords [ i ] - glocoords [ i ] ) / step;
                }
                testnatcoords [ dim ] -= step;
            }
            natcoordsvec = grad.householderQr().solve(errcoordsvec); //maybe grad.norm square
            for (unsigned i = 0; i < ndim; i++) {
                natcoords [ i ] -= natcoordsvec [ i ];
            }
            giveGlobalCoords(& glocoords, & natcoords);
            for (unsigned i = 0; i < ndim; i++) {
                errcoordsvec [ i ] = glocoords [ i ] - ( * x ) [ i ];
            }
            err = errcoordsvec.norm();
            it++;
        }
        if ( it == maxit ) {
            cerr << "Error in " << name << ": Natural coordinates were not found, error " << err << endl;
            cerr << "coords: " << ( * x ) [ 0 ] << " " << ( * x ) [ 1 ] << " " << ( * x ) [ 2 ] << endl;
            cerr << "natcoords: " << natcoords [ 0 ] << " " << natcoords [ 1 ] << " " << natcoords [ 2 ] << endl;
            cerr << "found coords: " << glocoords [ 0 ] << " " << glocoords [ 1 ] << " " << glocoords [ 2 ] << endl;
            cerr << "elem nodes: " << endl;
            for ( auto &n:nodes ) {
                cerr << n->givePoint() [ 0 ] << " " << n->givePoint() [ 1 ] << " " << n->givePoint() [ 2 ] << endl;
            }
            exit(1);
        }
    }
    return natcoords;
}

//////////////////////////////////////////////////////////
vector< vector< unsigned > >Element :: giveTraingulatedFaces()const {
    //vetrex 1, line 3, triangle 5, polygon 7, quad 9, tetra 10, brick 12, quadratic_triangle 22, quadratic_tetra 24, quadratic_brick 25
    vector< vector< unsigned > >tf;
    if ( vtk_cell_type == 9 ) {//triangle
        tf = { { 0, 1 }, { 1, 2 }, { 2, 0 } };
    } else if ( vtk_cell_type == 9 ) {     //quad
        tf = { { 0, 1 }, { 1, 2 }, { 2, 3 }, { 3, 0 } };
    } else if ( vtk_cell_type == 10 ) {     //tetra
        tf = { { 0, 1, 2 }, { 0, 2, 3 }, { 0, 1, 3 }, { 1, 2, 3 } };
    } else if ( vtk_cell_type == 12 ) {     //brick
        tf = { { 1, 0, 2 }, { 3, 2, 0 }, { 4, 5, 6 }, { 7, 6, 5 }, { 5, 0, 1 }, { 4, 1, 0 }, { 2, 1, 6 }, { 5, 6, 1 }, { 3, 2, 7 }, { 6, 7, 2 }, { 0, 3, 4 }, { 7, 4, 3 } };
    }
    return tf;
}

//////////////////////////////////////////////////////////
Vector Element :: findIntersectionsWithLine(Point *A, Point *B)const {
    double t;
    vector< double >intersections;
    vector< vector< unsigned > >tf = giveTraingulatedFaces();
    for (unsigned k = 0; k < tf.size(); k++) {
        t = find_intesection_of_segment_and_triangle(A, B, nodes [ tf [ k ] [ 0 ] ]->givePointPointer(), nodes [ tf [ k ] [ 1 ] ]->givePointPointer(), nodes [ tf [ k ] [ 2 ] ]->givePointPointer() );
        if ( t >= 0 ) {
            intersections.push_back(t);
        }
    }
    if ( intersections.size() == 0 ) {
        return Vector :: Zero(0);
    } else if ( intersections.size() == 1 ) {
        Vector extint(1);
        extint [ 0 ] = intersections [ 0 ];
        return extint;
    } else {
        Vector extint = Vector :: Zero(2);
        sort(intersections.begin(), intersections.end() );
        extint [ 0 ] = intersections [ 0 ];
        extint [ 1 ] = intersections.back();
        return extint;
    }
}

//////////////////////////////////////////////////////////
Vector Element :: giveShapeFunctions(const Point *x) const {
    Vector phi = Vector :: Zero( nodes.size() );
    shafunc->giveShapeF(x, phi);
    return phi;
}

//////////////////////////////////////////////////////////
Matrix Element :: giveShapeFunctionsGrad(const Point *x) const {
    Matrix phi = Matrix :: Zero(nodes.size(), ndim);
    shafunc->giveShapeFGrad(x, phi);
    return phi;
}

//////////////////////////////////////////////////////////
double Element :: giveDissipatedEnergy() const {
    double E = 0;
    for ( unsigned i = 0; i < inttype->giveNumIP(); i++ ) {
        E = stats [ i ]->giveDissipatedEnergyDensity() * inttype->giveIPWeight(i);
    }
    return E;
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
    IntegrDiscrete1 *it = dynamic_cast< IntegrDiscrete1 * >( inttype );
    it->setNumIP(1);
    it->setIPLocation(0, Point(0, 0, 0) );
    it->setIPWeight(0, 1);
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

//////////////////////////////////////////////////////////
Point Element :: giveApproxCenter() const {
    Point c = Point(0, 0, 0);
    for ( auto &k: nodes ) {
        c += k->givePoint();
    }
    if ( nodes.size() > 0 ) {
        c /= nodes.size();
    }
    return c;
}

//////////////////////////////////////////////////////////
Point Element :: giveIPLoc(unsigned k) const {
    Point nc = inttype->giveIPLocation(k);
    if ( shafunc->isInNaturalCoords() ) {
        Vector phi = Vector :: Zero( nodes.size() );
        shafunc->giveShapeF(& nc, phi);
        Point tc = Point(0, 0, 0);
        for ( unsigned i = 0; i < nodes.size(); i++ ) {
            tc += phi [ i ] * nodes [ i ]->givePoint();
        }
        return tc;
    } else {
        return nc;
    }
}

//////////////////////////////////////////////////////////
Vector Element :: giveBoundingBox()const {
    Vector bbox = Vector :: Zero(2 * ndim);
    Point t;
    for (unsigned k = 0; k < nodes.size(); k++) {
        t = nodes [ k ]->givePoint();
        for (unsigned i = 0; i < ndim; i++) {
            if ( k == 0 ) {
                bbox [ 2 * i ] = t [ i ];
                bbox [ 2 * i + 1 ] = t [ i ];
            } else {
                bbox [ 2 * i ] = min(bbox [ 2 * i ], t [ i ]);
                bbox [ 2 * i + 1 ] = max(bbox [ 2 * i + 1 ], t [ i ]);
            }
        }
    }
    return bbox;
}
