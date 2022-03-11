#include "integration.h"

using namespace std;

//////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////
// INTEGRATION POINTS - MASTER CLASS
//////////////////////////////////////////////////////////
void IntegrationType :: init() {
    cerr << "IntegrationType Error: you are using incorrect initialization function" << endl;
    exit(1);
}

//////////////////////////////////////////////////////////
void IntegrationType :: init(const vector< Node * > &nodes) {
    ( void ) nodes;
    IntegrationType :: init();
}

//////////////////////////////////////////////////////////
void IntegrationType :: init(const vector< Node * > &nodes, const vector< vector< unsigned > > &faces, Point *centroid) {
    ( void ) nodes;
    ( void ) faces;
    ( void ) centroid;
    IntegrationType :: init();
}

//////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////
// SINGLE POINT INTEGRATION FOR DISCRETE MODELS
//////////////////////////////////////////////////////////
void IntegrDiscrete1 :: init() {
    ip_locs.resize(1);
    ip_locs [ 0 ] = Point(0.5, 0, 0); //this will be rewritten by actual centroid of discrete facet
    ip_weights.resize(1);
    ip_weights [ 0 ] = 1.;
}

//////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////
// FOUR POINT INTEGRATION IN SQUARE
//////////////////////////////////////////////////////////
void IntegrQuad4 :: init() {
    ip_locs.resize(4);
    double s = 1. / sqrt(3);
    ip_locs [ 0 ] = Point(-s, -s, 0);
    ip_locs [ 1 ] = Point(s, -s, 0);
    ip_locs [ 2 ] = Point(s, s, 0);
    ip_locs [ 3 ] = Point(-s, s, 0);
    ip_weights.resize(4);
    for ( unsigned i = 0; i < 4; i++ ) {
        ip_weights [ i ] = 1.;
    }
}

//////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////
// EIGHT POINT INTEGRATION IN CUBE
//////////////////////////////////////////////////////////
void IntegrBrick8 :: init() {
    ip_locs.resize(8);
    double s = 1. / sqrt(3);
    ip_locs [ 0 ] = Point(-s, -s, -s);
    ip_locs [ 1 ] = Point(s, -s, -s);
    ip_locs [ 2 ] = Point(s, s, -s);
    ip_locs [ 3 ] = Point(-s, s, -s);
    ip_locs [ 4 ] = Point(-s, -s, s);
    ip_locs [ 5 ] = Point(s, -s, s);
    ip_locs [ 6 ] = Point(s, s, s);
    ip_locs [ 7 ] = Point(-s, s, s);
    ip_weights.resize(8);
    for ( unsigned i = 0; i < 8; i++ ) {
        ip_weights [ i ] = 1.;
    }
}


//////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////
// THREE POINT INTEGRATION IN TRIANGLE
//
//////////////////////////////////////////////////////////
void IntegrTri3 :: init() {
    ip_locs.resize(3);
    double s = 1. / 6.;
    ip_locs [ 0 ] = Point(s, s, 0.);
    ip_locs [ 1 ] = Point(4 * s, s, 0.);
    ip_locs [ 2 ] = Point(s, 4 * s, 0);
    ip_weights.resize(3);
    for ( unsigned i = 0; i < 3; i++ ) {
        ip_weights [ i ] = 1. / 6.;
    }
}


//////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////
// INTEGRATION IN POLYGON
void IntegrPolygon :: init(const vector< Node * > &nodes, const vector< vector< unsigned > > &faces, Point *centroid) {
    //based on isoparametric elements
    unsigned nnodes = nodes.size();

    if ( ip_type.compare("tri") == 0 ) {
        IntegrTri3 localINT; //integration type
        localINT.init();
        Linear2DTriShapeF localSF; //shape functions
        unsigned nIP = localINT.giveNumIP();
        MyVector phi = MyVector :: Zero(nIP);
        ip_locs.resize(nIP * nnodes);
        ip_weights.resize(nIP * nnodes);
        Point a, b, c;
        a = Point(0, 0, 0);
        b = Point(1, 0, 0);
        c = Point(0, 1, 0);
        vector< Point * >pp(nIP);
        pp [ 0 ] = & a;
        pp [ 1 ] = & b;
        pp [ 2 ] = & c;
        localSF.init(pp);
        double area;
        for ( unsigned i = 0; i < nnodes; i++ ) {
            area = triArea3D(nodes [ faces [ i ] [ 0 ] ]->givePointPointer(), nodes [ faces [ i ] [ 1 ] ]->givePointPointer(), centroid);
            for ( unsigned r = 0; r < nIP; r++ ) {
                localSF.giveShapeF(localINT.giveIPLocationPointer(r), phi);
                ip_locs [ 3 * i + r ] = nodes [ faces [ i ] [ 0 ] ]->givePoint() * phi [ 0 ] + nodes [ faces [ i ] [ 1 ] ]->givePoint() * phi [ 1 ] + ( * centroid ) * phi [ 2 ];
                ip_weights [ 3 * i + r ] = localINT.giveIPWeight(r) * area * 2.;
            }
        }
    } else if ( ip_type.compare("quad") == 0 ) {
        IntegrQuad4 localINT; //integration type
        localINT.init();
        Linear2DQuadShapeF localSF; //shape functions
        unsigned nIP = localINT.giveNumIP();

        MyVector phi = MyVector :: Zero(nIP);
        ip_locs.resize(nIP * nnodes);
        ip_weights.resize(nIP * nnodes);
        Point a = ( nodes [ faces [ nnodes - 1 ] [ 0 ] ]->givePoint() + nodes [ faces [ nnodes - 1 ] [ 1 ] ]->givePoint() ) / 2;
        Point b, c, d;
        d = ( * centroid );
        vector< Point * >pp(nIP);

        for ( unsigned i = 0; i < nnodes; i++ ) {
            c = ( nodes [ faces [ i ] [ 0 ] ]->givePoint() + nodes [ faces [ i ] [ 1 ] ]->givePoint() ) / 2;
            b = nodes [ i ]->givePoint();
            pp [ 0 ] = & a;
            pp [ 1 ] = & b;
            pp [ 2 ] = & c;
            pp [ 3 ] = & d;
            localSF.init(pp);

            for ( unsigned r = 0; r < nIP; r++ ) {
                localSF.giveShapeF(localINT.giveIPLocationPointer(r), phi);
                ip_locs [ 4 * i + r ] = Point(0, 0, 0);
                for ( unsigned s = 0; s < nIP; s++ ) {
                    ip_locs [ 4 * i + r ] += ( * pp [ s ] ) * phi [ s ];
                }
                ip_weights [ 4 * i + r ] = localINT.giveIPWeight(r) * localSF.giveJacobian( localINT.giveIPLocationPointer(r) );
            }
            a = c;
        }
    } else {
        cerr << "Error in " << name << ": ip_type '" << ip_type << "' not implemented" << endl;
        exit(1);
    }
}
