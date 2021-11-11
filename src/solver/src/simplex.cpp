#include "simplex.h"
#include "node.h"
#include "element_discrete.h"

//////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////
// SIMPLEX
//////////////////////////////////////////////////////////
void Simplex :: addElement(RigidBodyContact *rbc) {
    Particle *a;
    for ( unsigned i = 0; i < 2; i++ ) {
        a = static_cast< Particle * >( rbc->giveNode(i) );
        if ( find(nodes.begin(), nodes.end(), a) == nodes.end() ) {
            nodes.push_back(a);
        }
    }

    elems.push_back(rbc);
}

//////////////////////////////////////////////////////////
void Simplex :: init(unsigned ndim) {
    TrsNode *trs = dynamic_cast< TrsNode * >( center );
    if ( trs ) {
        transport = true;
        pressureDoF = trs->giveStartingDoF();
    }


    if ( ndim != 3 ) {
        return;             //TODO:extend for 2D models, need to take into account plane stress and strain
    }
    if ( nodes.size() == ndim + 1 ) {
        valid = true;

        //calculate DoF array
        DoFs.resize( ndim * nodes.size() );
        DoFweights.resize( ndim * nodes.size() );
        unsigned initDoF;
        unsigned pos = 0;
        Point volumeChangeWeights;
        unsigned j, k, l;
        double sign;
        double averageSide = 0;
        for ( unsigned i = 0; i < ndim + 1; i++ ) {
            j = ( i + 1 ) % ( ndim + 1 );
            k = ( i + 2 ) % ( ndim + 1 );
            l = ( i + 3 ) % ( ndim + 1 );
            averageSide += ( nodes [ j ]->givePoint() - nodes [ l ]->givePoint() ).norm();
            volumeChangeWeights = cross( nodes [ j ]->givePoint() - nodes [ l ]->givePoint(), nodes [ k ]->givePoint() - nodes [ l ]->givePoint() ) / 6.;
            volume = dot(nodes [ i ]->givePoint() - nodes [ l ]->givePoint(), volumeChangeWeights);
            sign = volume / abs(volume);
            initDoF = nodes [ i ]->giveStartingDoF();
            for ( unsigned v = 0; v < ndim; v++ ) {
                DoFs [ pos + v ] = initDoF + v;
                DoFweights [ pos + v ] = sign * volumeChangeWeights.giveCoord(v) / ndim; //divided by ndim, othewise total trace of strain vector would be returned
            }
            pos += ndim;
        }
        volume = abs(volume);

        averageSide /= 4;
        if ( volume < 0.05 * sqrt(2) / 12. * pow(averageSide, 3) || volume < 1e-25 ) {
            valid = false;                                                              //cancel simplexes that are in wrong shape, by comparing to 5% of volume of regular Tet
        }
    }
}

//////////////////////////////////////////////////////////
void Simplex :: findNeighbors() {
    //for non-valid simplices, they steal volumetric deformation from neighborhood
    if ( valid ) {
        return;
    }

    Simplex *s;
    vector< Node * >vertices;
    for ( auto &rbc: elems ) {
        vertices = rbc->giveVertices();
        for ( auto &v: vertices ) {
            s = v->giveSimplex();
            if ( s->isValid() ) {
                neighbors.insert(s);
            }
        }
    }
}

//////////////////////////////////////////////////////////
void Simplex :: computeVolumetricStrain(const Vector &fullDoFs) {
    volstrain = 0;
    if ( valid ) { //valid simplex taking care of its own
        for ( unsigned i = 0; i < DoFs.size(); i++ ) {
            volstrain += fullDoFs [ DoFs [ i ] ] * DoFweights [ i ];
        }
        volstrain /= volume; //mechanical volumetric stress, one third of strain tensor strace
    } else if ( neighbors.size() > 0 ) {  //steal volumetric strain from neigborhood
        for ( auto &simn: neighbors ) {
            volstrain += simn->giveVolumetricStrain();
        }
        volstrain /= neighbors.size();
    }

    if ( transport ) {
        pressure = fullDoFs [ pressureDoF ];
    }
}

//////////////////////////////////////////////////////////
double Simplex :: giveVolumetricStrain() const {
    return volstrain;
}

//////////////////////////////////////////////////////////
double Simplex :: givePressure() const {
    return pressure;
}
