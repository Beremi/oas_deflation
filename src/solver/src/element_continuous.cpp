#include "element_continuous.h"
#include "element_container.h"
#include "boundary_condition.h"
#include "material_tensorial.h"
#include "material_rve.h"

using namespace std;

//////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////
// 2D TRIANGULAR TRANSPORT ELEMENT
TrsprtTriangle :: TrsprtTriangle() : Element(2) {
    numOfNodes = 3;
    name = "TrsprtTriangle";
    vtk_cell_type = 5;
    shafunc = new Linear2DTriShapeF();
    inttype = new IntegrTri3();
    physicalFields [ 1 ] = true; //transport
}

//////////////////////////////////////////////////////////
Matrix TrsprtTriangle :: giveBMatrix(const Point *x) const {
    Matrix phiG = Matrix :: Zero(ndim, numOfNodes);
    shafunc->giveShapeFGrad(x, phiG);
    return phiG;
}

//////////////////////////////////////////////////////////
Matrix TrsprtTriangle :: giveHMatrix(const Point *x) const {
    Vector phi = Vector :: Zero( DoFids.size() );
    shafunc->giveShapeF(x, phi);
    Matrix H = Matrix :: Zero( 1, DoFids.size() );
    for ( unsigned k = 0; k < DoFids.size(); k++ ) {
        H(0, k) = phi(k);
    }
    return H;
}


//////////////////////////////////////////////////////////
Vector TrsprtTriangle :: giveStrain(unsigned i, const Vector &DoFs) {
    Vector strain = Element :: giveStrain(i, DoFs);

    //pressure at integration point for material model
    Vector pf = Hs [ i ] * DoFs;
    stats [ i ]->setParameterValue("pressure", pf [ 0 ]);

    return strain;
}

//////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////
// 2D QUADRILATERAL TRANSPORT ELEMENT
TrsprtQuad :: TrsprtQuad() {
    numOfNodes = 4;
    name = "TrsprtQuad";
    vtk_cell_type = 9;
    shafunc = new Linear2DQuadShapeF();
    inttype = new IntegrQuad4();
}

//////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////
// 3D TETRA TRANSPORT ELEMENT
TrsprtTetra :: TrsprtTetra() {
    ndim = 3;
    name = "TrsprtTetra";
    numOfNodes = 4;
    vtk_cell_type = 10;
    shafunc = new Linear3DTetraShapeF();
    inttype = new IntegrTetra4();
    physicalFields [ 1 ] = true; //transport

}


//////////////////////////////////////////////////////////
void TrsprtTetra :: init() {
    //check orientation of vertices
    if ( tetraVolumeSigned( nodes [ 0 ]->givePointPointer(), nodes [ 1 ]->givePointPointer(), nodes [ 2 ]->givePointPointer(), nodes [ 3 ]->givePointPointer() ) < 0. ) {
        Node *a = nodes [ 3 ];
        nodes [ 3 ] = nodes [ 2 ];
        nodes [ 2 ] = a;
    }
    TrsprtTriangle :: init();
}

//////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////
// 3D BRICK TRANSPORT ELEMENT
TrsprtBrick :: TrsprtBrick() {
    ndim = 3;
    name = "TrsprtBrick";
    numOfNodes = 8;
    vtk_cell_type = 12;
    shafunc = new Linear3DBrickShapeF();
    inttype = new IntegrBrick8();
    physicalFields [ 1 ] = true; //transport
}

//////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////
// 3D BRICK  TRANSPORT + TEMPERATURE COUPLED ELEMENT
TrsprtTemprtrCoupledBrick :: TrsprtTemprtrCoupledBrick() {
    ndim = 3;
    numOfNodes = 8;
    name = "TrsprtTemprtrCoupledBrick";
    vtk_cell_type = 12;
    shafunc = new Linear3DBrickShapeF();
    inttype = new IntegrBrick8();
}

//////////////////////////////////////////////////////////
Matrix TrsprtTemprtrCoupledBrick :: giveBMatrix(const Point *x) const {
    Matrix phiG = Matrix :: Zero(ndim, numOfNodes);
    shafunc->giveShapeFGrad(x, phiG);
    Matrix B = Matrix :: Zero(2 * ndim, numOfNodes * 2);
    for ( unsigned i = 0; i < numOfNodes; i++ ) {
        for ( unsigned v = 0; v < ndim; v++ ) {
            B(v, 2 * i) = B(ndim + v, 2 * i + 1) = phiG(v, i);
        }
    }
    return B;
}

//////////////////////////////////////////////////////////
Matrix TrsprtTemprtrCoupledBrick :: giveHMatrix(const Point *x) const {
    Vector phi = Vector :: Zero( DoFids.size() );
    shafunc->giveShapeF(x, phi);
    Matrix H = Matrix :: Zero(2,  numOfNodes * 2);
    for ( unsigned k = 0; k < numOfNodes; k++ ) {
        H(0, 2 * k) = H(1, 2 * k + 1) = phi(k);
    }
    return H;
}

//////////////////////////////////////////////////////////
Vector TrsprtTemprtrCoupledBrick :: giveStrain(unsigned i, const Vector &DoFs) {
    Vector strain = Element :: giveStrain(i, DoFs);

    Vector masters = Hs [ i ] * DoFs;
    stats [ i ]->setParameterValue("humidity", masters [ 0 ]);
    stats [ i ]->setParameterValue("temperature", masters [ 1 ]);
    return strain;
}

//////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////
// 2D TIANGULAR MECHANICAL ELEMENT
MechanicalTriangle :: MechanicalTriangle() : Element(2) {
    numOfNodes = 3;
    name = "MechanicalTriangle";
    vtk_cell_type = 5;
    shafunc = new Linear2DTriShapeF();
    inttype = new IntegrTri3(); //ONE IP suffices, but it does not work in extrapolation routine
    physicalFields [ 0 ] = true; //mechanics
    b_bar_integration_split = false;
}

//////////////////////////////////////////////////////////
void MechanicalTriangle :: readFromLine(std :: istringstream &iss, NodeContainer *fullnodes, MaterialContainer *fullmatrs) {
    Element :: readFromLine(iss, fullnodes, fullmatrs);

    //iss.clear(); // clear string stream
    //iss.seekg(0, iss.beg); //reset position in string stream

    string param;
    while (  iss >> param ) {
        if ( param.compare("b_bar_int") == 0 ) {
            b_bar_integration_split = true;
        }
    }
}

//////////////////////////////////////////////////////////
void MechanicalTriangle :: setIntegrationPointsAndWeights() {
    Element :: setIntegrationPointsAndWeights();
    if ( b_bar_integration_split ) {
        if ( ndim == 2 ) {
            TensMechMaterial *tmm = dynamic_cast< TensMechMaterial * >( mat );
            if ( tmm && tmm->isPlaneStress() ) {
                cerr << "WARNING: The B-bar integration for plane stress problem is not supported and it is therefore skipped" << endl;
            }
        }
        computeAverageBVolumeMatrix();
    }
}

//////////////////////////////////////////////////////////
void MechanicalTriangle :: computeAverageBVolumeMatrix() {
    Matrix phiG = Matrix :: Zero( ndim, nodes.size() );
    averageVolumeB = Matrix :: Zero( ndim, nodes.size() );
    double sumw = 0;
    for ( unsigned ipnum = 0; ipnum < inttype->giveNumIP(); ipnum++ ) {
        shafunc->giveShapeFGrad(inttype->giveIPLocationPointer(ipnum), phiG);
        averageVolumeB += phiG * inttype->giveIPWeight(ipnum);
        sumw += inttype->giveIPWeight(ipnum);
    }
    averageVolumeB /= ndim * sumw;
}

//////////////////////////////////////////////////////////
void MechanicalTriangle :: applyAverageVolumeB(Matrix &B, const Matrix phiG) const {
    int s = 0;
    for ( unsigned i = 0; i < numOfNodes; i++ ) {
        for ( unsigned d = 0; d < ndim; d++ ) {
            for ( unsigned dd = 0; dd < ndim; dd++ ) {
                B(dd, s + d) -= phiG(d, i) / ndim - averageVolumeB(d, i);
            }
        }
        s += nodes [ i ]->giveNumberOfDoFs();
    }
}


//////////////////////////////////////////////////////////
Matrix MechanicalTriangle :: giveBMatrix(const Point *x) const {
    Matrix phiG = Matrix :: Zero( ndim, nodes.size() );
    shafunc->giveShapeFGrad(x, phiG);
    Matrix B = Matrix :: Zero( 3, DoFids.size() );

    for ( unsigned i = 0; i < numOfNodes; i++ ) {
        B(0, 2 * i)     =   B(2, 2 * i + 1) =   phiG(0, i);
        B(1, 2 * i + 1) =   B(2, 2 * i)     =   phiG(1, i);
    }
    if ( b_bar_integration_split ) {
        applyAverageVolumeB(B, phiG);
    }
    return B;
}

//////////////////////////////////////////////////////////
Matrix MechanicalTriangle :: giveHMatrix(const Point *x) const {
    Vector phi = Vector :: Zero( nodes.size() );
    shafunc->giveShapeF(x, phi);
    Matrix H = Matrix :: Zero( ndim, DoFids.size() );
    for ( unsigned i = 0; i < ndim; i++ ) {
        for ( unsigned j = 0; j < numOfNodes; j++ ) {
            H(i, ndim * j + i) = phi(j);
        }
    }
    return H;
}

//////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////
// 2D QUADRILATERAL MECHANICAL ELEMENT
MechanicalQuad :: MechanicalQuad() {
    numOfNodes = 4;
    name = "MechanicalQuad";
    vtk_cell_type = 9;
    shafunc = new Linear2DQuadShapeF();
    inttype = new IntegrQuad4();
    b_bar_integration_split = false;
}

//////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////
// 3D TETRAHEDRAL MECHANICAL ELEMENT
MechanicalTetra :: MechanicalTetra() {
    ndim = 3;
    name = "MechanicalTetra";
    numOfNodes = 4;
    vtk_cell_type = 10;
    shafunc = new Linear3DTetraShapeF();
    inttype = new IntegrTetra4();
    physicalFields [ 0 ] = true; //mechanics
    b_bar_integration_split = false;
}

//////////////////////////////////////////////////////////
void MechanicalTetra :: init() {
    //check orientation of vertices
    if ( tetraVolumeSigned( nodes [ 0 ]->givePointPointer(), nodes [ 1 ]->givePointPointer(), nodes [ 2 ]->givePointPointer(), nodes [ 3 ]->givePointPointer() ) < 0. ) {
        Node *a = nodes [ 3 ];
        nodes [ 3 ] = nodes [ 2 ];
        nodes [ 2 ] = a;
    }
    Element :: init();
}

//////////////////////////////////////////////////////////
Matrix MechanicalTetra :: giveBMatrix(const Point *x) const {
    Matrix phiG = Matrix :: Zero( ndim, nodes.size() );
    shafunc->giveShapeFGrad(x, phiG);
    Matrix B = Matrix :: Zero( 6, DoFids.size() );
    for ( unsigned i = 0; i < numOfNodes; i++ ) {
        B(0, 3 * i)       =   B(4, 3 * i + 2)  =   B(5, 3 * i + 1) =   phiG(0, i);
        B(1, 3 * i + 1)   =   B(3, 3 * i + 2)  =   B(5, 3 * i)     =   phiG(1, i);
        B(2, 3 * i + 2)   =   B(3, 3 * i + 1)  =   B(4, 3 * i)     =   phiG(2, i);
    }
    if ( b_bar_integration_split ) {
        applyAverageVolumeB(B, phiG);
    }
    return B;
}

//////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////
// 3D BRICK MECHANICAL ELEMENT
MechanicalBrick :: MechanicalBrick() {
    ndim = 3;
    name = "MechanicalBrick";
    numOfNodes = 8;
    vtk_cell_type = 12;
    shafunc = new Linear3DBrickShapeF();
    inttype = new IntegrBrick8();
    b_bar_integration_split = false;
}

//////////////////////////////////////////////////////////
void MechanicalBrick :: init() {
    Element :: init();
}

//////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////
// 2D QUADRILATERAL COSSERAT MECHANICAL ELEMENT
CosseratQuad :: CosseratQuad() {
    ndim = 2;
    name = "CosseratQuad";
    numOfNodes = 4;
    vtk_cell_type = 9;
    shafunc = new Linear2DQuadShapeF();
    inttype = new IntegrQuad4();
    physicalFields [ 0 ] = true; //mechanics
    b_bar_integration_split = false;
}


//////////////////////////////////////////////////////////
Matrix CosseratQuad :: giveBMatrix(const Point *x) const {
    Matrix phiG = Matrix :: Zero(ndim, numOfNodes);
    shafunc->giveShapeFGrad(x, phiG);
    Vector phi = Vector :: Zero( nodes.size() );
    shafunc->giveShapeF(x, phi);
    Matrix B = Matrix :: Zero( 6, DoFids.size() );
    for ( unsigned i = 0; i < numOfNodes; i++ ) {
        //strains
        // 00 03
        // 02 01
        B(0, 3 * i) = B(3, 3 * i + 1) =  phiG(0, i);
        B(2, 3 * i) = B(1, 3 * i + 1) =  phiG(1, i);
        //effect of rotation = - LeviCivita * rotation
        //  0 -1
        //  1  0
        B(3, 3 * i + 2) = -phi(i);
        B(2, 3 * i + 2) =  phi(i);
        //curvatures
        // 04 05
        B(4, 3 * i + 2) =  phiG(0, i);
        B(5, 3 * i + 2) =  phiG(1, i);
    }
    if ( b_bar_integration_split ) {
        applyAverageVolumeB(B, phiG);
    }
    return B;
}

//////////////////////////////////////////////////////////
Matrix CosseratQuad :: giveHMatrix(const Point *x) const {
    Vector phi = Vector :: Zero( nodes.size() );
    shafunc->giveShapeF(x, phi);
    Matrix H = Matrix :: Zero( 3, DoFids.size() );         //2 transl, 1 rot
    for ( unsigned j = 0; j < numOfNodes; j++ ) {
        H(0, 3 * j) = H(1, 3 * j + 1) = H(1, 3 * j + 2) = phi(j);
    }
    return H;
}


//////////////////////////////////////////////////////////
void CosseratQuad :: computeMassMatrix() {
    MechanicalQuad :: computeMassMatrix();
    TensCosseratMechMaterial *cosmat = static_cast< TensCosseratMechMaterial * >( mat );
    for ( unsigned i = 0; i < 4; i++ ) {
        for ( unsigned j = 0; j < 4; j++ ) {
            massM(3 * i + 2, 3 * j + 2) *= 0.5 * pow(cosmat->giveParticleSize(), 2);
            massM(3 * i + 1, 3 * j + 2) = 0;
            massM(3 * i + 2, 3 * j + 1) = 0;
            massM(3 * i + 1, 3 * j) = 0;
            massM(3 * i + 2, 3 * j) = 0;
        }
    }
}

//////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////
// 3D BRICK COSSERAT MECHANICAL ELEMENT
CosseratBrick :: CosseratBrick() {
    ndim = 3;
    name = "CosseratBrick";
    numOfNodes = 8;
    vtk_cell_type = 12;
    shafunc = new Linear3DBrickShapeF();
    inttype = new IntegrBrick8();
    physicalFields [ 0 ] = true; //mechanics
}


//////////////////////////////////////////////////////////
Matrix CosseratBrick :: giveBMatrix(const Point *x) const {
    Matrix phiG = Matrix :: Zero(ndim, numOfNodes);
    shafunc->giveShapeFGrad(x, phiG);
    Vector phi = Vector :: Zero( nodes.size() );
    shafunc->giveShapeF(x, phi);
    Matrix B = Matrix :: Zero( 18, DoFids.size() );         //9 strains, 9 curvatures
    for ( unsigned i = 0; i < numOfNodes; i++ ) {
        //strains
        // 00 08 06
        // 07 01 04
        // 05 03 02
        B(0, 6 * i) = B(8, 6 * i + 1) = B(6, 6 * i + 2) = phiG(0, i);
        B(7, 6 * i) = B(1, 6 * i + 1) = B(4, 6 * i + 2) = phiG(1, i);
        B(5, 6 * i) = B(3, 6 * i + 1) = B(2, 6 * i + 2) = phiG(2, i);
        //effect of rotation = - LeviCivita * rotations
        //  0 -Z  Y
        //  Z  0 -X
        // -Y  X  0
        B(8, 6 * i + 5) = B(4, 6 * i + 3) = B(5, 6 * i + 4) = -phi(i);
        B(6, 6 * i + 4) = B(7, 6 * i + 5) = B(3, 6 * i + 3) =  phi(i);
        //curvatures
        // 09 17 15
        // 16 10 13
        // 14 12 11
        B(9, 6 * i + 3) = B(17, 6 * i + 4) = B(15, 6 * i + 5) = phiG(0, i);
        B(16, 6 * i + 3) = B(10, 6 * i + 4) = B(13, 6 * i + 5) = phiG(1, i);
        B(14, 6 * i + 3) = B(12, 6 * i + 4) = B(11, 6 * i + 5) = phiG(2, i);
    }
    return B;
}

//////////////////////////////////////////////////////////
Matrix CosseratBrick :: giveHMatrix(const Point *x) const {
    Vector phi = Vector :: Zero( nodes.size() );
    shafunc->giveShapeF(x, phi);
    Matrix H = Matrix :: Zero( 6, DoFids.size() );         //3 transl, 3 rot
    for ( unsigned v = 0; v < ndim; v++ ) {
        for ( unsigned j = 0; j < numOfNodes; j++ ) {
            H(v, 6 * j + v) = H(3 + v, 6 * j + 3 + v) = phi(j);
        }
    }
    return H;
}

//////////////////////////////////////////////////////////
void CosseratBrick :: computeMassMatrix() {
    MechanicalBrick :: computeMassMatrix();
    TensCosseratMechMaterial *cosmat = static_cast< TensCosseratMechMaterial * >( mat );
    unsigned k = 2 * ndim - 3; //number of rotations
    for ( unsigned i = 0; i < 4; i++ ) {
        for ( unsigned j = 0; j < 4; j++ ) {
            for (unsigned d = 0; d < k; d++) {
                for (unsigned e = 0; e < k; e++) {
                    if ( e == d ) {
                        massM(3 * i + e + ndim, 3 * j + d + ndim) *= 0.5 * pow(cosmat->giveParticleSize(), 2);
                    } else {
                        massM(3 * i + e + ndim, 3 * j + d + ndim) *= 0;
                    }
                }
                for (unsigned e = 0; e < ndim; e++) {
                    massM(3 * i + d + ndim, 3 * j + e) = 0;
                    massM(3 * i + e, 3 * j + d + ndim) = 0;
                }
            }
        }
    }
}

//////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////
// 2D QUADRILATERAL COSSERAT COUPLED MECHANICAL-TRANSPORT ELEMENT
CoupledCosseratTransportQuad :: CoupledCosseratTransportQuad() {
    ndim = 2;
    name = "CoupledCosseratQuad";
    numOfNodes = 4;
    vtk_cell_type = 9;
    shafunc = new Linear2DQuadShapeF();
    inttype = new IntegrQuad4();
    physicalFields [ 0 ] = true; //mechanics
    physicalFields [ 1 ] = true; //transport
}


//////////////////////////////////////////////////////////
Matrix CoupledCosseratTransportQuad :: giveBMatrix(const Point *x) const {
    Matrix phiG = Matrix :: Zero(ndim, numOfNodes);
    shafunc->giveShapeFGrad(x, phiG);
    Vector phi = Vector :: Zero( nodes.size() );
    shafunc->giveShapeF(x, phi);
    Matrix B = Matrix :: Zero( 8, DoFids.size() );         //4 strains, 2 curvatures, 2 pressure gradients
    for ( unsigned i = 0; i < numOfNodes; i++ ) {
        //strains
        // 00 03
        // 02 01
        B(0, 4 * i) = B(3, 4 * i + 1) =  phiG(0, i);
        B(2, 4 * i) = B(1, 4 * i + 1) =  phiG(1, i);
        //effect of rotation = - LeviCivita * rotation
        //  0 -1
        //  1  0
        B(3, 4 * i + 2) = -phi(i);
        B(2, 4 * i + 2) =  phi(i);
        //curvatures
        // 04 05
        B(4, 4 * i + 2) =  phiG(0, i);
        B(5, 4 * i + 2) =  phiG(1, i);
        //pressure gradient
        // 06 07
        B(6, 4 * i + 3) =  phiG(0, i);
        B(7, 4 * i + 3) =  phiG(1, i);
    }
    return B;
}

//////////////////////////////////////////////////////////
Matrix CoupledCosseratTransportQuad :: giveHMatrix(const Point *x) const {
    Vector phi = Vector :: Zero( nodes.size() );
    shafunc->giveShapeF(x, phi);
    Matrix H = Matrix :: Zero( 4, DoFids.size() );         //2 transl, 1 rot, 1 pressure
    for ( unsigned j = 0; j < numOfNodes; j++ ) {
        H(0, 4 * j) = H(1, 4 * j + 1) = H(2, 4 * j + 2) = H(3, 4 * j + 3) = phi(j);
    }
    return H;
}


//////////////////////////////////////////////////////////
void CoupledCosseratTransportQuad :: init() {
    Element :: init();

    ParticleWithTransport *cp;
    for ( auto &n: nodes ) {
        cp = dynamic_cast< ParticleWithTransport * >( n );
        if ( cp == nullptr ) {
            cerr << name << " requires nodes to be inhereted from ParticleWithTransport" << endl;
            exit(1);
        }
    }
}

//////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////
// 3D BRICK COSSERAT COUPLED MECHANICAL-TRANSPORT ELEMENT
CoupledCosseratTransportBrick :: CoupledCosseratTransportBrick() {
    ndim = 3;
    name = "CoupledCosseratBrick";
    numOfNodes = 8;
    vtk_cell_type = 12;
    shafunc = new Linear3DBrickShapeF();
    inttype = new IntegrBrick8();
    physicalFields [ 1 ] = true; //transport
}

//////////////////////////////////////////////////////////
Vector CoupledCosseratTransportBrick :: giveStrain(unsigned i, const Vector &DoFs) {
    Vector strain = CosseratBrick :: giveStrain(i, DoFs);

    //pressure at integration point for material model
    Vector pf = Hs [ i ] * DoFs;
    stats [ i ]->setParameterValue("pressure", pf [ 6 ]);
    //volumetric strain at integration point for material model
    stats [ i ]->setParameterValue("volumetricStrain", ( strain [ 0 ] + strain [ 1 ] + strain [ 2 ] ) / 3.);

    return strain;
}

/*
 * //////////////////////////////////////////////////////////
 * MyVector CoupledCosseratTransportBrick :: giveInternalForces(const MyVector &DoFs, bool frozen, double timeStep) {
 *  MyVector intF = CosseratBrick :: giveInternalForces(DoFs, frozen, timeStep);
 *
 *
 *  //Biot effect
 *  DiscreteCoupledRVEMaterialStatus *couplm;
 *  double intTrsprtSource;
 *  MyMatrix H;
 *  for ( unsigned k = 0; k < inttype->giveNumIP(); k++ ) {
 *      couplm = static_cast< DiscreteCoupledRVEMaterialStatus * >( stats [ k ] );
 *      intTrsprtSource = couplm->giveInternalSource();
 *      for ( unsigned j = 0; j < 8; j++ ) {
 *          intF [ 7*j+6 ] -= Hs [ k ] [ 6 ] [ 7*j+6 ] *intTrsprtSource *inttype->giveIPWeight(k);
 *      }
 *  }
 *
 *  return intF;
 * }
 */

//////////////////////////////////////////////////////////
Matrix CoupledCosseratTransportBrick :: giveBMatrix(const Point *x) const {
    Matrix phiG = Matrix :: Zero(ndim, numOfNodes);
    shafunc->giveShapeFGrad(x, phiG);
    Vector phi = Vector :: Zero( nodes.size() );
    shafunc->giveShapeF(x, phi);
    Matrix B = Matrix :: Zero( 21, DoFids.size() );         //9 strains, 9 curvatures, 3 pressure gradients
    for ( unsigned i = 0; i < numOfNodes; i++ ) {
        //strains
        // 00 08 06
        // 07 01 04
        // 05 03 02
        B(0, 7 * i) = B(8, 7 * i + 1) = B(6, 7 * i + 2) = phiG(0, i);
        B(7, 7 * i) = B(1, 7 * i + 1) = B(4, 7 * i + 2) = phiG(1, i);
        B(5, 7 * i) = B(3, 7 * i + 1) = B(2, 7 * i + 2) = phiG(2, i);
        //effect of rotation = - LeviCivita * rotations
        //  0 -Z  Y
        //  Z  0 -X
        // -Y  X  0
        B(8, 7 * i + 5) = B(4, 7 * i + 3) = B(5, 7 * i + 4) = -phi(i);
        B(6, 7 * i + 4) = B(7, 7 * i + 5) = B(3, 7 * i + 3) =  phi(i);
        //curvatures
        // 09 17 15
        // 16 10 13
        // 14 12 11
        B(9, 7 * i + 3) = B(17, 7 * i + 4) = B(15, 7 * i + 5) = phiG(0, i);
        B(16, 7 * i + 3) = B(10, 7 * i + 4) = B(13, 7 * i + 5) = phiG(1, i);
        B(14, 7 * i + 3) = B(12, 7 * i + 4) = B(11, 7 * i + 5) = phiG(2, i);
        //pressure gradient
        // 18 19 20
        B(18, 7 * i + 6) = phiG(0, i);
        B(19, 7 * i + 6) = phiG(1, i);
        B(20, 7 * i + 6) = phiG(2, i);
    }
    return B;
}

//////////////////////////////////////////////////////////
Matrix CoupledCosseratTransportBrick :: giveHMatrix(const Point *x) const {
    Vector phi = Vector :: Zero( nodes.size() );
    shafunc->giveShapeF(x, phi);
    Matrix H = Matrix :: Zero( 7, DoFids.size() );         //3 transl, 3 rot, 1 pressure
    for ( unsigned j = 0; j < numOfNodes; j++ ) {
        for ( unsigned v = 0; v < ndim; v++ ) {
            H(v, 7 * j + v) = H(3 + v, 7 * j + 3 + v) = phi(j);
        }
        H(6, 7 * j + 6) = phi(j); //pressure gradient
    }
    return H;
}

//////////////////////////////////////////////////////////
void CoupledCosseratTransportBrick :: computeDampingMatrix() {
    Element :: computeDampingMatrix();
}

//////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////
// 3D BRICK COSSERAT MECHANICAL ELEMENT
CoupledCosseratBrickWithDependentUpperZLayer :: CoupledCosseratBrickWithDependentUpperZLayer() {
    name = "CoupledCosseratBrickWithDependentUpperZLayer";
    bindlayers = false; // TRUE not working, need fix
}


//////////////////////////////////////////////////////////
void CoupledCosseratBrickWithDependentUpperZLayer :: setIntegrationPointsAndWeights() {
    if ( bindlayers ) {
        stats.resize(inttype->giveNumIP() / 2.);
    } else {
        stats.resize( inttype->giveNumIP() );
    }
    for ( unsigned k = 0; k < inttype->giveNumIP(); k++ ) {
        inttype->setIPWeight( k, inttype->giveIPWeight(k) * shafunc->giveJacobian(inttype->giveIPLocationPointer(k) ) );

        if ( k < 4 || !bindlayers ) {
            stats [ k ] = mat->giveNewMaterialStatus(this, k);
            DiscreteCoupledRVEMaterialStatus *dcrve = dynamic_cast< DiscreteCoupledRVEMaterialStatus * >( stats [ k ] );
            if ( dcrve != nullptr ) {
                Matrix r = Matrix :: Zero(3, 3);
                r(2, 2) = 1.;
                double xcent = 0;
                double ycent = 0;
                for ( unsigned i = 0; i < 4; i++ ) {
                    xcent += nodes [ i ]->givePoint().x();
                    ycent += nodes [ i ]->givePoint().y();
                }
                double alpha = atan2(ycent / 4., xcent / 4.);
                r(0, 0) = r(1, 1) = cos(alpha);
                r(0, 1) = sin(alpha);
                r(1, 0) = -sin(alpha);

                dcrve->setReferenceSystemDirections(r);
            }
        }
    }
};

//////////////////////////////////////////////////////////
Matrix CoupledCosseratBrickWithDependentUpperZLayer :: giveStiffnessMatrix(string matrixType) const {
    unsigned nDoFs = DoFids.size();
    Matrix K = Matrix :: Zero(nDoFs, nDoFs);
    Matrix D(0, 0);

    unsigned q;
    for ( unsigned i = 0; i < inttype->giveNumIP(); i++ ) {
        if ( i < 4 || !bindlayers ) {
            q = i;
        } else {
            q = i - 4;
        }
        D = stats [ q ]->giveStiffnessTensor(matrixType);
        K += Bs [ i ].transpose() * D * ( Bs [ i ] * inttype->giveIPWeight(i) );
    }
    return K;
}

//////////////////////////////////////////////////////////
Vector CoupledCosseratBrickWithDependentUpperZLayer :: giveInternalForces(const Vector &DoFs, bool frozen, double timeStep) {
    Vector intF = Vector :: Zero( DoFids.size() );
    Vector stress;
    for ( unsigned i = 0; i < inttype->giveNumIP(); i++ ) {
        if ( i < 4 || !bindlayers ) {
            if ( frozen ) {
                stress = stats [ i ]->giveStressWithFrozenIntVars(giveStrain(i, DoFs), timeStep);  //frozen internal variables
            } else {
                stress = stats [ i ]->giveStress(giveStrain(i, DoFs), timeStep); //full evaluation of stress including change of state variables
            }
        } else {
            stress = stats [ i - 4 ]->giveStressWithFrozenIntVars(giveStrain(i, DoFs), timeStep);  //frozen internal variables
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
Vector CoupledCosseratBrickWithDependentUpperZLayer :: integrateInternalSources() {
    Vector intS = Vector :: Zero( DoFids.size() );
    Vector intmats;
    unsigned q;
    for ( unsigned i = 0; i < inttype->giveNumIP(); i++ ) {
        if ( i < 4 || !bindlayers ) {
            q = i;
        } else {
            q = i - 4;
        }
        intmats = stats [ q ]->giveInternalSource();
        intS += Hs [ i ].transpose() * ( intmats * inttype->giveIPWeight(i) );
    }

    return intS;
}


//////////////////////////////////////////////////////////
void CoupledCosseratBrickWithDependentUpperZLayer :: computeDampingMatrix() {
    unsigned nDoFs = DoFids.size();
    dampC = Matrix :: Zero(nDoFs, nDoFs);
    Matrix c;
    unsigned q;
    for ( unsigned i = 0; i < inttype->giveNumIP(); i++ ) {
        if ( i < 4 || !bindlayers ) {
            q = i;
        } else {
            q = i - 4;
        }
        c = stats [ q ]->giveDampingTensor();
        dampC += Hs [ i ].transpose() * ( c * inttype->giveIPWeight(i) ) * Hs [ i ];
    }
}

//////////////////////////////////////////////////////////
void CoupledCosseratBrickWithDependentUpperZLayer :: computeMassMatrix() {
    unsigned nDoFs = DoFids.size();
    massM = Matrix :: Zero(nDoFs, nDoFs);
    Matrix c;
    unsigned q;
    for ( unsigned i = 0; i < inttype->giveNumIP(); i++ ) {
        if ( i < 4 || !bindlayers ) {
            q = i;
        } else {
            q = i - 4;
        }
        c = stats [ q ]->giveMassTensor();
        massM += Hs [ i ].transpose() * ( c * inttype->giveIPWeight(i) ) * Hs [ i ];
    }
}

//////////////////////////////////////////////////////////
Vector CoupledCosseratBrickWithDependentUpperZLayer :: giveStrain(unsigned i, const Vector &DoFs) {
    Vector strain = CosseratBrick :: giveStrain(i, DoFs);

    //pressure at integration point for material model
    Vector pf = Hs [ i ] * DoFs;
    unsigned q;
    if ( i < 4 || !bindlayers ) {
        q = i;
    } else {
        q = i - 4;
    }
    stats [ q ]->setParameterValue("pressure", pf [ 6 ]);
    //volumetric strain at integration point for material model
    stats [ q ]->setParameterValue("volumetricStrain", ( strain [ 0 ] + strain [ 1 ] + strain [ 2 ] ) / 3.);

    return strain;
}
