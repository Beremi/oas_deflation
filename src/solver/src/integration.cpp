#include "integration.h"

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
