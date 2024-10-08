#include "element_ldpm.h"
#include "element_container.h"
#include "boundary_condition.h"
#include "material_vectorial.h"
#include "model.h"

using namespace std;

//////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////
// LDPM COUPLED TETRA
LDPMTetra :: LDPMTetra(unsigned dim) : Element{dim} {
    if ( ndim != 3 ) {
        cerr << "LDPMTetra implemented only in 3D" << endl;
        exit(1);
    }
    numOfNodes = 4;
    nodes.resize(4);
    name = "LDPMTetra";
    vtk_cell_type = 10;
    shafunc = new NullShapeF(3);
    inttype = new IntegrLDPM12();

    vert.resize(12);
    lengths.resize(12);
    volumes.resize(12);
    areas.resize(12);
    normals.resize(12);
    t1s.resize(12);
    t2s.resize(12);
    R.resize(12);

    nodecodes = { 0, 1, 0, 1, 0, 2, 0, 2, 0, 3, 0, 3,  1, 2, 1, 2, 1,  3, 1, 3,  2, 3, 2, 3 };
    vertcodes = { 8, 1, 1, 7, 2, 9, 7, 2, 9, 3, 3, 8, 10, 4, 4, 7, 5, 10, 8, 5, 10, 6, 6, 9 }; //last point is always centroid at position 0
    //only for help
    //opposed_verts = { 10, 9, 8, 7}

    //OLD NUMBERING
    //nodecodes = { 1, 2, 1, 3, 2, 3, 0, 2, 0, 3, 2, 3, 0, 1, 0, 3, 1, 3, 0, 1, 0, 2, 1, 2 };
    //vertcodes = { 1, 2, 3, 1, 1, 4, 6, 5, 5, 7, 4, 5, 8, 9, 7, 8, 8, 3, 9, 10, 10, 6, 2, 10 }; //last point is always centroid at position 0
}

//////////////////////////////////////////////////////////
unsigned LDPMTetra :: giveOppositeSurfaceVertexToNode(unsigned k) const{
    return 10-k;
}

//////////////////////////////////////////////////////////
vector<unsigned> LDPMTetra :: giveOppositeFacetsToNode(unsigned k) const{
    vector<unsigned> f(3);
    if (k==0){
        f[0] = 6; f[1] = 8; f[2] = 10;
    }else if (k==1){
        f[0] = 2; f[1] = 4; f[2] = 11;
    }else if (k==2){
        f[0] = 0; f[1] = 5; f[2] = 9;
    }else if (k==3){
        f[0] = 1; f[1] = 3; f[2] = 7;
    }else{
        cerr << "Error " << name << ": node does not exist" << endl;
        exit(1);
    }
    return f;
}

/////////////////////////////////////////////////////////
void LDPMTetra :: readFromLine(istringstream &iss, NodeContainer *fullnodes, MaterialContainer *fullmatrs) {
    unsigned num;
    for ( unsigned i = 0; i < 4; i++ ) {
        iss >> num;
        nodes [ i ]  = fullnodes->giveNode(num);
    }
    for ( unsigned i = 0; i < 11; i++ ) {
        iss >> num;
        vert [ i ]  = fullnodes->giveNode(num);
    }
    iss >> num;
    mat = fullmatrs->giveMaterial(num);
}




//////////////////////////////////////////////////////////

void LDPMTetra :: checkNodeType() const {
    //check that nodes are particles
    for ( unsigned i = 0; i < 4; i++ ) {
        Particle *p = dynamic_cast< Particle * >( nodes [ i ] );
        if ( !p ) {
            cerr << "Error in " << name << ": nodes must be inherited from Particle, " << nodes [ i ]->giveName() << " provided" << endl;
            exit(1);
        }
    }
}

//////////////////////////////////////////////////////////
void LDPMTetra :: setIntegrationPointsAndWeights() {
    stats.resize(12);

    for ( unsigned i = 0; i < 12; i++ ) {
        //true face normal
        Point n = ( vert [ vertcodes [ 2 * i ] ]->givePoint() - vert [ vertcodes [ 2 * i + 1 ] ]->givePoint() ).cross( vert [ 0 ]->givePoint() - vert [ vertcodes [ 2 * i ] ]->givePoint() );
        n /= n.norm();
        //contact vector
        normals [ i ] = nodes [ nodecodes [ 2 * i + 1 ] ]->givePoint() - nodes [ nodecodes [ 2 * i ] ]->givePoint();
        lengths [ i ] = normals [ i ].norm();
        normals [ i ] /= lengths [ i ];
        inttype->setIPLocation(i, ( vert [ vertcodes [ 2 * i  ] ]->givePoint() + vert [ vertcodes [ 2 * i + 1 ] ]->givePoint() + vert [ 0 ]->givePoint() ) / 3.);
        areas [ i ] = triArea3D(vert [ vertcodes [ 2 * i  ] ]->givePointPointer(), vert [ vertcodes [ 2 * i + 1 ] ]->givePointPointer(), vert [ 0 ]->givePointPointer() );
        if ( n.norm() == n.norm() ) { //NaN test
            areas [ i ] *= abs(n.dot(normals [ i ]) );    //projection of area
        }

        //t1 = inttype->giveIPLocation(i)-vert [ 0 ]->givePoint();   this is wrong for irregular TET
        // coordinate swap for tangential vector according to https://orbit.dtu.dk/files/126824972/onb_frisvad_jgt2012_v2.pdf
        Point arbit(sqrt(2.), -sqrt(3.), M_PI);
        if ( ( normals [ i ] - arbit ).norm() < 1e-3 ) {
            t1s[ i ] = arbit.cross(normals [ i ]);
        } else {
            // the following results in zeros in stiffness matrix in case of normal in direction of any of global base axes
            if ( abs( normals [ i ].x() ) > 1e-3 ) {
                t1s[ i ] = Point(-normals [ i ].y() / normals [ i ].x(), 1, 0);
            } else if ( abs( normals [ i ].y() ) > 1e-3 ) {
                t1s[ i ] = Point(0, -normals [ i ].z() / normals [ i ].y(), 1);
            } else {
                t1s[ i ] = Point( 1, 0, -normals [ i ].x() / normals [ i ].z() );
            }
        }
        t1s[ i ].normalize();
        t2s[ i ] = normals [ i ].cross(t1s [i] );
        R [ i ] = Matrix :: Zero(3, 3);
        R [ i ](0, 0) = normals [ i ].x();
        R [ i ](0, 1) = normals [ i ].y();
        R [ i ](0, 2) = normals [ i ].z();
        R [ i ](1, 0) = t1s [ i ].x();
        R [ i ](1, 1) = t1s [ i ].y();
        R [ i ](1, 2) = t1s [ i ].z();
        R [ i ](2, 0) = t2s [ i ].x();
        R [ i ](2, 1) = t2s [ i ].y();
        R [ i ](2, 2) = t2s [ i ].z();

        volumes [ i ] = areas [ i ] * lengths [ i ] / ndim;
        inttype->setIPWeight(i, lengths [ i ] * areas [ i ] / ndim);
        stats [ i ] = mat->giveNewMaterialStatus(this, i);
    }
}

//////////////////////////////////////////////////////////
void LDPMTetra :: init() {
    Element :: init(); //calling base class method;
    checkNodeType();

    //check that material is VectMechMat
    VectMechMaterial *p = dynamic_cast< VectMechMaterial * >( mat );
    if ( !p ) {
        cerr << "Error in " << name << ": material must be inherited from VectMechMaterial, " << mat->giveName() << " provided" << endl;
        exit(1);
    }


    //weights for volumetric calculations
    volWeights.resize(12);
    Point volumeChangeWeights;
    unsigned j, k, l;
    double sign;
    double averageSide = 0;
    for ( unsigned i = 0; i < 4; i++ ) {
        j = ( i + 1 ) % 4;
        k = ( i + 2 ) % 4;
        l = ( i + 3 ) % 4;
        averageSide += ( nodes [ j ]->givePoint() - nodes [ l ]->givePoint() ).norm();
        volumeChangeWeights = ( nodes [ j ]->givePoint() - nodes [ l ]->givePoint() ).cross( nodes [ k ]->givePoint() - nodes [ l ]->givePoint() ) / 6.;
        volume = ( nodes [ i ]->givePoint() - nodes [ l ]->givePoint() ).dot(volumeChangeWeights);
        sign = volume / abs(volume);
        for ( unsigned v = 0; v < 3; v++ ) {
            volWeights [ 3 * i + v ] = sign * volumeChangeWeights(v)  / 3; //divided by ndim, othewise total trace of strain vector would be returned
        }
    }

    volume = abs(volume);
    averageSide /= 4;
    if ( volume < 1e-25 ) {
        cerr << name << "Error: wrong geometry" << endl;
        exit(1);
    }

    double faceVolume = 0;
    for ( unsigned i = 0; i < 12; i++ ) {
        faceVolume += volumes [ i ];
    }
    if ( abs(faceVolume - volume) / max(faceVolume, volume) > 1e-4 ) {
        cerr << "LDPM Tetra Warning: " << idx << " - total volume is " << volume << ", volume from faces is " << faceVolume << ", error " << abs(faceVolume - volume) / max(faceVolume, volume)  << endl;
        //for(unsigned i=0; i<12; i++) cout << areas[i]*lengths[i]/ndim << endl;
        //exit(1);
    }
}

//////////////////////////////////////////////////////////

Matrix LDPMTetra :: giveBMatrix(unsigned k) const {
    unsigned nA = nodecodes [ 2 * k ];
    unsigned nB = nodecodes [ 2 * k + 1 ];
    Matrix B = Matrix :: Zero(3, 24);
    Particle *a = static_cast< Particle * >( nodes [ nA ] );
    Matrix Aa = a->giveRigidBodyMotionMatrix(inttype->giveIPLocationPointer(k) );
    a = static_cast< Particle * >( nodes [ nB ] );
    Matrix Ab = a->giveRigidBodyMotionMatrix(inttype->giveIPLocationPointer(k) );

    for ( unsigned i = 0; i < 3; i++ ) {
        for ( unsigned j = 0; j < 6; j++ ) {
            B(i, nA * 6 + j) = -Aa(i, j);
            B(i, nB * 6 + j) = Ab(i, j);
        }
    }
    return ( R [ k ] * B ) / lengths [ k ];
}

//////////////////////////////////////////////////////////
Matrix LDPMTetra :: giveHMatrix(const Point *x) const {
    ( void ) x;
    return Matrix :: Zero(3, 24);  // NOTE JK: this should be based on ndim
}

//////////////////////////////////////////////////////////
Vector LDPMTetra :: giveStrain(unsigned i, const Vector &DoFs) {
    //compute volumetric strain

    if ( i == 0 ) { //first IP
        volumetricStrain = 0;
        unsigned r = 0;
        for ( unsigned k = 0; k < 4; k++ ) {
            for ( unsigned p = 0; p < 3; p++ ) {
                volumetricStrain += DoFs [ r + p ] * volWeights [ 3 * k + p ];
            }
            r += nodes [ k ]->giveNumberOfDoFs();
        }
        volumetricStrain /= volume; //mechanical volumetric stress, one third of strain tensor strace
    }
    stats [ i ]->setParameterValue("volumetric_strain", volumetricStrain);

    return Element :: giveStrain(i, DoFs);
};

//////////////////////////////////////////////////////////
Matrix LDPMTetra :: giveStiffnessMatrix(string matrixType) const {
    return Element :: giveStiffnessMatrix(matrixType) * ndim; //ndim needs to be included here for discrete elements
}

//////////////////////////////////////////////////////////
void LDPMTetra :: computeDampingMatrix() {
    dampC = giveStiffnessMatrix("elastic") * 1e-15;           //rough fix of zeros, here can be anything
}

//////////////////////////////////////////////////////////
Vector LDPMTetra :: giveInternalForces(const Vector &DoFs, bool frozen, double timeStep) {
    return Element :: giveInternalForces(DoFs, frozen, timeStep) * ndim; //ndim needs to be included here for discrete elements
}

//////////////////////////////////////////////////////////
Vector LDPMTetra :: integrateLoad(BodyLoad *vl, double time) const {
    return Element :: integrateLoad(vl, time) / ndim;
}

//////////////////////////////////////////////////////////
Vector LDPMTetra :: integrateInternalSources() {
    return Element :: integrateInternalSources() / ndim;
}

//////////////////////////////////////////////////////////
void LDPMTetra :: computeMassMatrix() {
    massM = Matrix :: Zero(24, 24);
    VectMechMaterialStatus *mechstat;
    double density;
    double tetvol;
    Matrix tetI;
    unsigned nodeID;
    Point tetcentr, diff;
    vector< Point * >tetnodes(4);
    vector< Point >reltetnodes(4);
    tetnodes [ 3 ] = vert [ 0 ]->givePointPointer();

    for ( unsigned i = 0; i < 12; i++ ) {
        mechstat = static_cast< VectMechMaterialStatus * >( stats [ 0 ] );
        density = mechstat->giveDensity();
        for ( unsigned j = 0; j < 2; j++ ) {
            nodeID = nodecodes [ 2 * i + j ];
            tetnodes [ 0 ] = nodes [ nodeID ]->givePointPointer();
            tetnodes [ 1 ] = vert [ vertcodes [ 2 * i + j ] ]->givePointPointer();
            tetnodes [ 2 ] = vert [ vertcodes [ 2 * i + 1 - j ] ]->givePointPointer();
            tetcentr = ( ( * tetnodes [ 0 ] ) + ( * tetnodes [ 1 ] ) + ( * tetnodes [ 2 ] ) + ( * tetnodes [ 3 ] ) ) / 4;
            for ( unsigned k = 0; k < 4; k++ ) {
                reltetnodes [ k ] = ( * tetnodes [ k ] ) - tetcentr;
            }
            tetvol = tetraVolumeSigned(& reltetnodes [ 0 ], & reltetnodes [ 1 ], & reltetnodes [ 2 ], & reltetnodes [ 3 ]);
            tetI = tetraInertia3D(& reltetnodes [ 0 ], & reltetnodes [ 1 ], & reltetnodes [ 2 ], & reltetnodes [ 3 ]);
            if ( tetvol < 0 ) {
                tetvol *= -1;           //should not happen for correctly generated mesh
            }
            if ( tetI(0, 0) < 0 ) {
                tetI *= -1;              //should not happen for correctly generated mesh
            }
            for ( unsigned k = 0; k < 3; k++ ) {
                massM(nodeID * 6 + k, nodeID * 6 + k) += density * tetvol;
            }
            diff = tetcentr - ( * tetnodes [ 0 ] );
            massM(nodeID * 6, nodeID * 6 + 4) += density * tetvol * diff [ 2 ];
            massM(nodeID * 6, nodeID * 6 + 5) -= density * tetvol * diff [ 1 ];
            massM(nodeID * 6 + 1, nodeID * 6 + 3) -= density * tetvol * diff [ 2 ];
            massM(nodeID * 6 + 1, nodeID * 6 + 5) += density * tetvol * diff [ 0 ];
            massM(nodeID * 6 + 2, nodeID * 6 + 3) += density * tetvol * diff [ 1 ];
            massM(nodeID * 6 + 2, nodeID * 6 + 4) -= density * tetvol * diff [ 0 ];
            massM(nodeID * 6 + 3, nodeID * 6 + 3) += density * ( tetI(0, 0) + tetvol * ( pow( ( diff [ 1 ] ), 2 ) + pow( ( diff [ 2 ] ), 2 ) ) );
            massM(nodeID * 6 + 4, nodeID * 6 + 4) += density * ( tetI(1, 1) + tetvol * ( pow( ( diff [ 0 ] ), 2 ) + pow( ( diff [ 2 ] ), 2 ) ) );
            massM(nodeID * 6 + 5, nodeID * 6 + 5) += density * ( tetI(2, 2) + tetvol * ( pow( ( diff [ 0 ] ), 2 ) + pow( ( diff [ 1 ] ), 2 ) ) );
            massM(nodeID * 6 + 3, nodeID * 6 + 4) += density * ( tetI(0, 1) - tetvol * ( ( diff [ 0 ] ) * ( diff [ 1 ] ) ) );
            massM(nodeID * 6 + 3, nodeID * 6 + 5) += density * ( tetI(0, 2) - tetvol * ( ( diff [ 0 ] ) * ( diff [ 2 ] ) ) );
            massM(nodeID * 6 + 4, nodeID * 6 + 5) += density * ( tetI(1, 2) - tetvol * ( ( diff [ 1 ] ) * ( diff [ 2 ] ) ) );
        }
    }
    for ( nodeID = 0; nodeID < 4; nodeID++ ) {
        massM(nodeID * 6 + 4, nodeID * 6) = massM(nodeID * 6, nodeID * 6 + 4);
        massM(nodeID * 6 + 5, nodeID * 6) = massM(nodeID * 6, nodeID * 6 + 5);
        massM(nodeID * 6 + 3, nodeID * 6 + 1) = massM(nodeID * 6 + 1, nodeID * 6 + 3);
        massM(nodeID * 6 + 5, nodeID * 6 + 1) = massM(nodeID * 6 + 1, nodeID * 6 + 5);
        massM(nodeID * 6 + 3, nodeID * 6 + 2) = massM(nodeID * 6 + 2, nodeID * 6 + 3);
        massM(nodeID * 6 + 4, nodeID * 6 + 2) = massM(nodeID * 6 + 2, nodeID * 6 + 4);
        massM(nodeID * 6 + 4, nodeID * 6 + 3) = massM(nodeID * 6 + 3, nodeID * 6 + 4);
        massM(nodeID * 6 + 5, nodeID * 6 + 3) = massM(nodeID * 6 + 3, nodeID * 6 + 5);
        massM(nodeID * 6 + 5, nodeID * 6 + 4) = massM(nodeID * 6 + 4, nodeID * 6 + 5);
    }
}


//////////////////////////////////////////////////////////
vector< unsigned >LDPMTetra :: giveFacetVertCodes(unsigned k) const {
    vector< unsigned >v;
    v.resize(3);
    v [ 0 ] = vertcodes [ 2 * k ];
    v [ 1 ] = vertcodes [ 2 * k + 1 ];
    v [ 2 ] = 0;
    return v;
}

//////////////////////////////////////////////////////////
vector< unsigned >LDPMTetra :: giveFacetNodeCodes(unsigned k) const {
    vector< unsigned >v;
    v.resize(2);
    v [ 0 ] = nodecodes [ 2 * k ];
    v [ 1 ] = nodecodes [ 2 * k + 1 ];
    return v;
}

//////////////////////////////////////////////////////////
void LDPMTetra :: giveValues(string code, Vector &result) const {
    if ( code.compare("volumetric_strain") == 0 ) {
        result.resize(0);
        result [ 0 ] = volumetricStrain;
    } else {
        Element :: giveValues(code, result);
    }
};

/*
 * //////////////////////////////////////////////////////////
 * void LDPMTetra :: extrapolateIPValuesToNodes(string code, vector< Vector > &result, Vector &weights) const {
 *  Vector ipres;
 *  giveIPValues(code, 0, ipres);
 *  Vector A = giveVectorToNode(0, 0);
 *  Vector B = giveVectorToNode(1, 0);
 *  size_t d;
 *
 *  weights.resize(2);
 *  weights [ 0 ] = giveVolumeAssociatedWithNode(0);
 *  weights [ 1 ] = giveVolumeAssociatedWithNode(1);
 *
 *  if ( ipres.size() == 0 ) {   //empty answer
 *      result.resize(0);
 *  } else if ( ipres.size() == 1 ) {   //scalar times vector //needs to be checked, probably not theoretically correct
 *      result.resize( ndim );
 *      for ( d = 0; d < ndim; d++ ) {
 *          result [ d ].resize(2);
 *          result [ d ] [ 0 ] =  area * ipres [ 0 ] * abs(A [ d ]);
 *          result [ d ] [ 1 ] =  area * ipres [ 0 ] * abs(B [ d ]);
 *      }
 *  } else if ( ipres.size() == A.size() ) { //vector times vector of same length, symmetrization
 *      //transform result to xyz
 *      Vector ipresglobal = transformVectorToXYZ(ipres);
 *
 *      //dyadic product
 *      unsigned k = A.size();
 *      result.resize( ( k * ( k - 1 ) ) / 2 + k );
 *      for ( d = 0; d < ( k * ( k - 1 ) ) / 2 + k; d++ ) {
 *          result [ d ].resize(2);
 *      }
 *      //diagonal
 *      for ( d = 0; d < k; d++ ) {
 *          result [ d ] [ 0 ] =  area * ipresglobal [ d ] * A [ d ];
 *          result [ d ] [ 1 ] =  -area * ipresglobal [ d ] * B [ d ];
 *      }
 *      //off diagonal
 *      if ( k == 2 ) {
 *          result [ 2 ] [ 0 ] =  area * ( ipresglobal [ 1 ] * A [ 0 ] + ipresglobal [ 0 ] * A [ 1 ] ) / 2.;
 *          result [ 2 ] [ 1 ] = -area * ( ipresglobal [ 1 ] * B [ 0 ] + ipresglobal [ 0 ] * B [ 1 ] ) / 2.;
 *      } else if ( k == 3 ) {
 *          result [ 3 ] [ 0 ] =  area * ( ipresglobal [ 1 ] * A [ 2 ] + ipresglobal [ 2 ] * A [ 1 ] ) / 2.;
 *          result [ 3 ] [ 1 ] = -area * ( ipresglobal [ 1 ] * B [ 2 ] + ipresglobal [ 2 ] * B [ 1 ] ) / 2.;
 *          result [ 4 ] [ 0 ] =  area * ( ipresglobal [ 2 ] * A [ 0 ] + ipresglobal [ 0 ] * A [ 2 ] ) / 2.;
 *          result [ 4 ] [ 1 ] = -area * ( ipresglobal [ 2 ] * B [ 0 ] + ipresglobal [ 0 ] * B [ 2 ] ) / 2.;
 *          result [ 5 ] [ 0 ] =  area * ( ipresglobal [ 1 ] * A [ 0 ] + ipresglobal [ 0 ] * A [ 1 ] ) / 2.;
 *          result [ 5 ] [ 1 ] = -area * ( ipresglobal [ 1 ] * B [ 0 ] + ipresglobal [ 0 ] * B [ 1 ] ) / 2.;
 *      } else {
 *          cerr << "Error in " << name << ": transformation of matrix of size " << k << " to vector not implemented" << endl;
 *          exit(1);
 *      }
 *  } else {
 *      cerr << "Error in " << name << ": dyadic product of vectors of different length in function extrapolateIPValuesToNodes" << endl;
 *      exit(1);
 *  }
 * }
 *
 */



//////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////
// LDPM COUPLED TETRA
LDPMCoupledTetra :: LDPMCoupledTetra() : LDPMTetra{3} {
    name = "LDPMCoupledTetra";
}

//////////////////////////////////////////////////////////
Vector LDPMCoupledTetra :: giveStrain(unsigned i, const Vector &DoFs) {
    Vector res;
    vert[0]->giveDoFBasedValues("pressure",masterModel->giveSolver(),res);
    double pressure = 0;
    if (res.size()==1) pressure=res[0];
    stats [ i ]->setParameterValue("pressure",pressure);

    return LDPMTetra :: giveStrain(i, DoFs);
};

//////////////////////////////////////////////////////////
void LDPMCoupledTetra :: giveValues(string code, Vector &result) const {
    LDPMTetra :: giveValues(code, result);
};

//////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////
// LDPM TRANSPORT
//////////////////////////////////////////////////////////
LDPMCoupledTransport::LDPMCoupledTransport(ElementContainer* allelems):DiscreteTrsprtCoupledElem(3){
    elems = allelems;
    name = "LDPMCoupledTransport";
};

//////////////////////////////////////////////////////////
Vector LDPMCoupledTransport :: giveStrain(unsigned i, const Vector &DoFs) {
    //crack opening
    /*
    double crackInNeighborhood = 0;
    double crackVolume = 0.;
    double elem_crack_opening;

    for ( auto &f: friends ) {
        elem_crack_opening = 0.;
        for ( unsigned k = 0; k < f->giveNumIP(); k++ ) {
            f->giveIPValues("crack_opening", k, res);
            if ( res.size() == 1 ) {
                elem_crack_opening += abs(res [ 0 ]);
            }
        }
        crackInNeighborhood += pow(elem_crack_opening / f->giveNumIP(), 3) * friendsweight [ m ];
        if ( ndim == 3 ) {
            crackVolume += ( elem_crack_opening / f->giveNumIP() ) * f->giveArea() * length / f->givePerimeter();
        } else if ( ndim == 2 ) {
            crackVolume += ( elem_crack_opening / f->giveNumIP() ) * f->giveArea();
        }
        m++;
    }
    stats [ 0 ]->setParameterValue("crack_opening", crackInNeighborhood);
    stats [ 0 ]->setParameterValue("crack_volume", crackVolume);
    */

    //Biot effect
    double volStrain = tetA->giveVolumetricStrain();
    if (tetB) { volStrain += tetB->giveVolumetricStrain();
        volStrain /=2.;
    }
    stats [ 0 ]->setParameterValue("volumetric_strain", volStrain); //trace divided by dimension to obtain mechanical volumetric strain

    return DiscreteTrsprtElem :: giveStrain(i, DoFs);
};

//////////////////////////////////////////////////////////
void LDPMCoupledTransport::readFromLine(std :: istringstream &iss, NodeContainer *fullnodes, MaterialContainer *fullmatrs){
    (void) fullnodes;
    iss >> LDPMTetraIDA;
    iss >> LDPMTetraIDB;
    if (LDPMTetraIDA==LDPMTetraIDB){
        cerr << "Error " << name << ": two identical tetrahedras supplied, " << LDPMTetraIDA << endl;
    }
    unsigned num;
    iss >> num;
    mat = fullmatrs->giveMaterial(num);

    
    string code;
    while ( iss >> code ) {
        if ( code.compare("BolanderCapacityMatrix") == 0 ) {
            BolanderCapacityMatrix = true;
        }
    }
}

//////////////////////////////////////////////////////////
void LDPMCoupledTransport::init(){
    tetA = dynamic_cast<LDPMTetra*>(elems->giveElement(LDPMTetraIDA));
    tetB = dynamic_cast<LDPMTetra*>(elems->giveElement(LDPMTetraIDB));
    if ((!tetA) or (!tetB)){
        cout << "Error " << name << ": the element numbers are not LDPM tetras"<< endl;
        exit(1);
    } 
    nodes.resize(2);
    nodes[0]=tetA->giveCentroid();
    nodes[1]=tetB->giveCentroid();

    LDPMsideA = 5;
    LDPMsideB = 5;

    vector<Node*> nodesA = tetA->giveNodes();
    vector<Node*> nodesB = tetB->giveNodes();
    bool found;
    for (unsigned k=0; k<4; k++){
        found=false;
        for (unsigned l=0; l<4; l++){
            if (nodesA[k]==nodesB[l]){
                found=true;
                break;
            }
        }
        if (!found){
            if (LDPMsideA!=5){
                cerr << "Error " << name << ": LDPM tetras are not neighbors, "<< LDPMTetraIDA << " " << LDPMTetraIDB << endl;
                for (auto a:nodesA) cout << a->giveID() << " ";
                cout << endl;
                for (auto a:nodesB) cout << a->giveID() << " ";
                cout << endl;
                exit(1);
            } else LDPMsideA=k;
        }
    }
    for (unsigned k=0; k<4; k++){
        found=false;
        for (unsigned l=0; l<4; l++){
            if (nodesB[k]==nodesA[l]){
                found=true;
                break;
            }
        }
        if (!found){
            if (LDPMsideB!=5){
                cerr << "Error " << name << ": LDPM tetras are not neighbors, "<< LDPMTetraIDA << " " << LDPMTetraIDA << endl;
                for (auto a:nodesA) cout << a->giveID() << " ";
                cout << endl;
                for (auto a:nodesB) cout << a->giveID() << " ";
                cout << endl;
                exit(1);
            } else LDPMsideB=k;
        }
    }
    vert.resize(3);
    unsigned p = 0;
    for (unsigned k=0; k<4; k++){   
        if (k==LDPMsideA) continue;
        vert[p] = nodesA[k];
        p++;
    }

    DiscreteTrsprtCoupledElem::init();
}


//////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////
// LDPM TRANSPORT BOUNDARY
//////////////////////////////////////////////////////////
LDPMCoupledTransportBoundary :: LDPMCoupledTransportBoundary(ElementContainer* allelems):LDPMCoupledTransport(allelems){
    name = "LDPMCoupledTransportBoundary";
};

//////////////////////////////////////////////////////////
void LDPMCoupledTransportBoundary::init(){
    tetA = dynamic_cast<LDPMTetra*>(elems->giveElement(LDPMTetraIDA));
    LDPMsideA = LDPMTetraIDB;
    
    nodes.resize(2);
    nodes[0]=tetA->giveCentroid();
    nodes[1]=tetA->giveVertex(tetA->giveOppositeSurfaceVertexToNode(LDPMsideA));

    vector<Node*> nodesA = tetA->giveNodes();    
    vert.resize(3);
    unsigned p = 0;
    for (unsigned k=0; k<4; k++){   
        if (k==LDPMsideA) continue;
        vert[p] = nodesA[k];
        p++;
    }

    DiscreteTrsprtElem::init();
}




