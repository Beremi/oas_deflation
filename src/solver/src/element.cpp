#include "element.h"
#include "element_container.h"
#include "boundary_condition.h"

//////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////
// BASIC ELEMENT - MASTER CLASS
Element :: ~Element() {
    for ( vector< MaterialStatus * > :: iterator e = stats.begin(); e != stats.end(); ++e ) {
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
    stats.resize(inttype->giveNumIP() );
    for ( unsigned k = 0; k < inttype->giveNumIP(); k++ ) {
        stats [ k ] = mat->giveNewMaterialStatus(this, k);
        inttype->setIPWeight(k, inttype->giveIPWeight(k) * shafunc->giveJacobian( inttype->giveIPLocationPointer(k) ) );
    }
};


//////////////////////////////////////////////////////////
void Element :: readFromLine(istringstream &iss, NodeContainer *fullnodes, MaterialContainer *fullmatrs) {
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
    unsigned totalDoFs = 0;
    for ( vector< Node * > :: const_iterator n = nodes.begin(); n != nodes.end(); ++n ) {
        totalDoFs += ( * n )->giveNumberOfDoFs();
    }

    initIntegration();
    setIntegrationPointsAndWeights();

    DoFids.resize(totalDoFs);
    unsigned i = 0;
    unsigned k;
    for ( vector< Node * > :: const_iterator n = nodes.begin(); n != nodes.end(); ++n ) {
        k = ( * n )->giveStartingDoF();
        for ( unsigned s = 0; s < ( * n )->giveNumberOfDoFs(); s++, i++ ) {
            DoFids [ i ] = k + s;
        }
    }
    outDoFs = totalDoFs; //basic elems will alway have input = output

    Bs.resize(inttype->giveNumIP() );
    Hs.resize(inttype->giveNumIP() );
    for ( k = 0; k < inttype->giveNumIP(); k++ ) {
        Bs [ k ] = giveBMatrix(inttype->giveIPLocationPointer(k) );
        Hs [ k ] = giveHMatrix(inttype->giveIPLocationPointer(k) );
    }
}

//////////////////////////////////////////////////////////
vector< unsigned >Element :: giveDoFsInDirection(unsigned dir) const {
    vector< unsigned >DoFinDir(nodes.size() );
    for ( unsigned i = 0; i < nodes.size(); i++ ) {
        DoFinDir [ i ] = nodes [ i ]->giveStartingDoF() + dir;
    }
    return DoFinDir;
}

//////////////////////////////////////////////////////////
void Element :: initMaterialStatuses() {
    unsigned num = 0;
    for ( vector< MaterialStatus * > :: iterator m = stats.begin(); m != stats.end(); ++m, num++ ) {
        ( * m )->init();
    }
}

//////////////////////////////////////////////////////////
void Element :: updateMaterialStatuses() {
    for ( vector< MaterialStatus * > :: iterator m = stats.begin(); m != stats.end(); ++m ) {
        ( * m )->update();
    }
}

//////////////////////////////////////////////////////////
double Element :: giveValue(string code) const {
    ( void ) code;
    return 0;
};

//////////////////////////////////////////////////////////
double Element :: giveIPValue(string code, unsigned ipnum) const {
    if ( code.compare("x") == 0 ) {
        return inttype->giveIPLocationPointer(ipnum)->getX();
    } else if ( code.compare("y") == 0 ) {
        return inttype->giveIPLocationPointer(ipnum)->getY();
    } else if ( code.compare("z") == 0 ) {
        return inttype->giveIPLocationPointer(ipnum)->getZ();
    } else if ( code.compare("materialID") == 0 || code.compare("materialId") == 0 ) {
        return stats [ ipnum ]->giveMaterial()->giveId();
    } else {
        return stats [ ipnum ]->giveValue(code);
    }
};

//////////////////////////////////////////////////////////
Matrix Element :: giveStiffnessMatrix(string matrixType) const {
    unsigned nDoFs = DoFids.size();
    Matrix K(nDoFs, nDoFs);
    Matrix D(0, 0);
    for ( unsigned i = 0; i < inttype->giveNumIP(); i++ ) {
        D = stats [ i ]->giveStiffnessTensor(matrixType, ndim);
        K += Bs [ i ].transpose() * D * ( Bs [ i ] * inttype->giveIPWeight(i) );
    }
    return K;
}

//////////////////////////////////////////////////////////
Vector Element :: giveInternalForces(const Vector &DoFs, bool frozen, double timeStep) {
    Vector intF(DoFids.size() );
    Vector stress;
    for ( unsigned i = 0; i < inttype->giveNumIP(); i++ ) {
        if ( frozen ) {
            stress = stats [ i ]->giveStressWithFrozenIntVars(giveStrain(i, DoFs), timeStep);  //frozen internal variables
        } else {
            stress = stats [ i ]->giveStress(giveStrain(i, DoFs), timeStep); //full evaluation of stress including change of state variables
        }
        intF  += Bs [ i ].transpose() * (  stress * inttype->giveIPWeight(i) );
    }

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
    Vector intS(DoFids.size() );
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
    Matrix M(nDoFs, nDoFs);
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
    Matrix M(nDoFs, nDoFs);
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
    Vector load(nDoFs);
    double fvalue;
    unsigned dir = vl->giveDirection();
    for ( unsigned i = 0; i < inttype->giveNumIP(); i++ ) {
        fvalue = vl->giveValue(inttype->giveIPLocationPointer(i), time);
        for ( unsigned j = 0; j < nDoFs; j++ ) {
            load [ j ] += Hs [ i ] [ dir ] [ j ] *fvalue *inttype->giveIPWeight(i);
        }
    }
    return load;
}

//////////////////////////////////////////////////////////
void Element :: changeMaterial(Material *newmat) {
    this->mat = newmat;
    for ( unsigned p = 0; p < stats.size(); p++ ) {
        if (stats [ p ] != nullptr) delete stats [ p ];
        stats [ p ] = this->mat->giveNewMaterialStatus(this, p);
    }
    this->initMaterialStatuses();
}

//////////////////////////////////////////////////////////
bool Element :: giveGlobalCoords(Point *x, const Point *xn) const{
    Vector phi(nodes.size());
    shafunc->giveShapeF(xn, phi);
    *x = Point(0,0,0);
    for (unsigned n=0; n<nodes.size(); n++){
        *x += nodes[n]->givePoint()*phi[n];  
    }    
    return true;
}


//////////////////////////////////////////////////////////
bool Element :: isPointInside(Point* xn, const Point *x) const {
    //initial screening
    Point maxc(-1e10,-1e10,-1e10);
    Point minc( 1e10, 1e10, 1e10);
    Point *p;
    for(auto &n: nodes){
        p = n->givePointPointer();
        for(unsigned c=0; c<ndim; c++){
            maxc.setCoord(c,max(maxc.giveCoord(c),p->giveCoord(c)));
            minc.setCoord(c,min(minc.giveCoord(c),p->giveCoord(c)));
        }
    }   
    for(unsigned c=0; c<ndim; c++){
        if (x->giveCoord(c)>maxc.giveCoord(c) || x->giveCoord(c)<minc.giveCoord(c)) return false;
    }
    
    //find natural coordinates
    //initial estimation
    Point center;
    *xn = Point(0,0,0);    
    Point aux, diff, diffC;
    Point size = maxc-minc;
    giveGlobalCoords(&center, xn);
        
    //initial estimation
    aux = ((*x) - center)*2.;    
    for(unsigned c=0; c<ndim; c++) xn->setCoord(c, aux.giveCoord(c)/size.giveCoord(c));
    giveGlobalCoords(&aux, xn);      

    int i = 0;
    int max_i = 500;
    double maxerror = 1;
    diff = aux-(*x);
    while(maxerror>1e-4 && i<max_i){
        for(unsigned c=0; c<ndim; c++){
            if (abs(diff.giveCoord(c)/size.giveCoord(c))<1e-8) continue;
            if (diffC.giveCoord(c)>1e-16)  xn->setCoord(c, xn->giveCoord(c) - xn->giveCoord(c)*diff.giveCoord(c)/diffC.giveCoord(c)); 
            else xn->setCoord(c, xn->giveCoord(c) - 2.*diff.giveCoord(c)/size.giveCoord(c));
        }

        giveGlobalCoords(&aux, xn);
        diff = aux-(*x);
        diffC = aux - center;
        maxerror = 0.;
        for(unsigned c=0; c<ndim; c++) maxerror = max(maxerror, abs(diff.giveCoord(c)/size.giveCoord(c)) ); 
        i++;
    }
    if(i==max_i) return false;
    
    //check natural coordinates are inside limits
    //TODO: works only for brick and quadrilateral
    for(unsigned c=0; c<ndim; c++) if (abs(xn->giveCoord(c))>1.) return false;
    return true;
}

//////////////////////////////////////////////////////////
Vector Element :: giveElemDoFsFromFullDoFs( const Vector &FullDoFs) const{
    Vector elemDoFs(DoFids.size());
    for (unsigned i=0; i<DoFids.size(); i++){
        elemDoFs[i] = FullDoFs[DoFids[i]];
    }
    return elemDoFs;
}
