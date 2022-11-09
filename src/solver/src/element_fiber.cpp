#include "element_fiber.h"
#include "element_discrete.h"
#include "element_container.h"
#include "model.h"

using namespace std;
//////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////
// FIBER

Fiber :: Fiber(const unsigned dim) {
    ndim = dim;
    name = "Fiber";
    diam = 0;
    solution_order = 0;
    volume = 0;
    numOfNodes = 2;
    nodes.resize(2);
    vtk_cell_type = 3;
    shafunc = new Linear1DLineShapeF();
    inttype = new IntegrFiber();
}

//////////////////////////////////////////////////////////
void Fiber :: init() {
    dirVector = nodes [ 1 ]->givePoint() - nodes [ 0 ]->givePoint();
    length = dirVector.norm();
    dirVector /= length;
}

//////////////////////////////////////////////////////////
void Fiber :: readFromLine(istringstream &iss, NodeContainer *fullnodes, MaterialContainer *fullmatrs) {
    unsigned num;
    iss >> num;
    nodes [ 0 ] = fullnodes->giveNode(num);
    iss >> num;
    nodes [ 1 ] = fullnodes->giveNode(num);
    iss >> diam >> num;
    mat = fullmatrs->giveMaterial(num);
}

//////////////////////////////////////////////////////////
void Fiber :: createNewCrossing(Point intersec, RigidBodyContact *rbc) {
    stats.push_back(mat->giveNewMaterialStatus(this, stats.size() ) );
    IntegrFiber *intf = static_cast< IntegrFiber * >( inttype );
    intf->addNewIP(intersec);
    contacts.push_back(rbc);
    positions.push_back( ( intersec - nodes [ 0 ]->givePoint() ).norm() );
}

//////////////////////////////////////////////////////////
void Fiber :: setUpCrossings() {
    std :: vector< Node * >rbcnodes;
    std :: vector< Node * > :: iterator pos;
    Node *kn;
    vector< unsigned >np(2 * stats.size() );
    unsigned totalDoFs = 0;
    unsigned i = 0;
    for ( auto &r:contacts ) {
        for ( unsigned p = 0; p < 2; p++ ) {
            kn = r->giveNode(p);
            pos = std :: find(rbcnodes.begin(), rbcnodes.end(), kn);
            if ( pos != rbcnodes.end() ) {
                np [ 2 * i + p ] = pos - rbcnodes.begin();
            } else {
                np [ 2 * i + p ] = rbcnodes.size();
                rbcnodes.push_back(kn);
                totalDoFs += kn->giveNumberOfDoFs();
            }
        }
        i++;
    }

    DoFids.resize(totalDoFs);
    i = 0;
    unsigned k;
    for ( auto &rn:rbcnodes ) {
        k = rn->giveStartingDoF();
        for ( unsigned s = 0; s < rn->giveNumberOfDoFs(); s++, i++ ) {
            DoFids [ i ] = k + s;
        }
    }
    outDoFs = totalDoFs; //basic elems will always have input = output

    unsigned nodedof = 3 * ( ndim - 1 );
    Bs.resize( inttype->giveNumIP() );
    Hs.resize( inttype->giveNumIP() );
    for ( k = 0; k < inttype->giveNumIP(); k++ ) {
        Hs [ k ] = Element :: giveHMatrix(k);

        Bs [ k ] = Matrix :: Zero(ndim, DoFids.size() );
        Matrix rbcB = contacts [ k ]->giveBMatrix(inttype->giveIPLocationPointer(i) ) * contacts [ k ]->giveLength() ;
        for ( unsigned cc = 0; cc < nodedof; cc++ ) {
            for ( unsigned rr = 0; rr < ndim; rr++ ) {
                Bs [ k ](rr, nodedof *np [ 2 * k ] + cc)   = rbcB(rr, cc);
                Bs [ k ](rr, nodedof *np [ 2 * k + 1 ] + cc) = rbcB(rr, cc + nodedof);
            }
        }
    }

    //set stress and strain vectors at integration points
    for ( k = 0; k < inttype->giveNumIP(); k++ ) {
        stats [ k ]->initializeStressAndStrainVector(Bs [ k ].rows() );
        stats [k] ->init();
    }
}


//////////////////////////////////////////////////////////
Matrix Fiber :: giveBMatrix(const Point *x) const {
    ( void ) x;
    return Matrix :: Zero(DoFids.size(), ndim);
};

//////////////////////////////////////////////////////////
Matrix Fiber :: giveHMatrix(const Point *x) const {
    ( void ) x;
    return Matrix(DoFids.size(), DoFids.size() );
};
