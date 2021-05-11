#include "element_continuous.h"
#include "element_container.h"
#include "boundary_condition.h"
#include "material_RVE.h"

//////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////
// 2D QUADRILATERAL TRANSPORT ELEMENT
TrsprtQuad :: TrsprtQuad() {
    ndim = 2;
    numOfNodes = 4;
    name = "TrsprtQuad";
    vtk_cell_type = 9;
    shafunc = new Linear2DQuadShapeF();
    inttype = new IntegrQuad4();
}

//////////////////////////////////////////////////////////
Matrix TrsprtQuad :: giveBMatrix(const Point *x) const {
    Matrix phiG(ndim, numOfNodes);
    shafunc->giveShapeFGrad(x, nodes, phiG);
    return phiG;
}

//////////////////////////////////////////////////////////
Matrix TrsprtQuad :: giveHMatrix(const Point *x) const {
    Vector phi(DoFids.size() );
    shafunc->giveShapeF(x, phi);
    Matrix H(1, DoFids.size() );
    for ( unsigned k = 0; k < DoFids.size(); k++ ) {
        H [ 0 ] [ k ] = phi [ k ];
    }
    return H;
}


//////////////////////////////////////////////////////////
Vector TrsprtQuad :: giveStrain(unsigned i, const Vector &DoFs) {
    Vector strain0 = Element :: giveStrain(i, DoFs);
    Vector strain(strain0.size() + 1);
    for ( unsigned i = 0; i < strain0.size(); i++ ) {
        strain [ i ] = strain0 [ i ];
    }

    //evaluate pressure at gauss point to account for nonlinearity
    strain [ strain0.size() ] = matrix_vector_multiply(giveHMatrix(inttype->giveIPLocationPointer(i) ), DoFs) [ 0 ];

    return strain;
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
    Matrix phiG(ndim, numOfNodes);
    shafunc->giveShapeFGrad(x, nodes, phiG);
    Matrix B(2*ndim, numOfNodes*2);
    for(unsigned i=0; i<numOfNodes; i++){
        for(unsigned v=0; v<ndim; v++){
            B[v][2*i] = B[ndim+v][2*i+1] = phiG[v][i];
        }
    }
    return B;
}

//////////////////////////////////////////////////////////
Matrix TrsprtTemprtrCoupledBrick :: giveHMatrix(const Point *x) const {
    Vector phi(DoFids.size() );
    shafunc->giveShapeF(x, phi);
    Matrix H(2,  numOfNodes*2 );
    for ( unsigned k = 0; k < numOfNodes; k++ ) {
        H [ 0 ] [ 2*k ] = H [ 1 ] [ 2*k+1 ] = phi [ k ];
    }
    return H;
}

//////////////////////////////////////////////////////////
Vector TrsprtTemprtrCoupledBrick :: giveStrain(unsigned i, const Vector &DoFs){
    Vector strain = Element :: giveStrain(i, DoFs);
    
    Vector masters;
    for ( unsigned i = 0; i < inttype->giveNumIP(); i++ ) {
        masters = Hs[i]*DoFs;
        stats [ i ]->setParameterValue("humidity", masters[0]);
        stats [ i ]->setParameterValue("temperature", masters[1]);
    }    
    return strain;
}

//////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////
// 2D QUADRILATERAL MECHANICAL ELEMENT
MechanicalQuad :: MechanicalQuad() {
    ndim = 2;
    numOfNodes = 4;
    name = "MechanicalQuad";
    vtk_cell_type = 9;
    shafunc = new Linear2DQuadShapeF();
    inttype = new IntegrQuad4();
}

//////////////////////////////////////////////////////////
Matrix MechanicalQuad :: giveBMatrix(const Point *x) const {
    Matrix phiG(ndim, nodes.size() );
    shafunc->giveShapeFGrad(x, nodes, phiG);
    Matrix B(3, DoFids.size() );

    for ( unsigned i = 0; i < numOfNodes; i++ ) {
        B [ 0 ] [ 2 * i ]     =   B [ 2 ] [ 2 * i + 1 ] =   phiG [ 0 ] [ i ];
        B [ 1 ] [ 2 * i ]     =   B [ 0 ] [ 2 * i + 1 ] =   0.;
        B [ 1 ] [ 2 * i + 1 ]   =   B [ 2 ] [ 2 * i ]   =   phiG [ 1 ] [ i ];
    }
    return B;
}

//////////////////////////////////////////////////////////
Matrix MechanicalQuad :: giveHMatrix(const Point *x) const {
    Vector phi(nodes.size() );
    shafunc->giveShapeF(x, phi);
    Matrix H(ndim, DoFids.size() );
    for ( unsigned i = 0; i < ndim; i++ ) {
        for ( unsigned j = 0; j < numOfNodes; j++ ) {
            H [ i ] [ ndim * j + i ] = phi [ j ];
        }
    }
    return H;
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
}

//////////////////////////////////////////////////////////
Matrix MechanicalBrick :: giveBMatrix(const Point *x) const {
    Matrix phiG(ndim, nodes.size() );
    shafunc->giveShapeFGrad(x, nodes, phiG);
    Matrix B(6, DoFids.size() );
    for ( unsigned i = 0; i < numOfNodes; i++ ) {
        B [ 0 ] [ 3 * i ]       =   B [ 4 ] [ 3 * i + 2 ]  =   B [ 5 ] [ 3 * i + 1 ] =   phiG [ 0 ] [ i ];
        B [ 1 ] [ 3 * i + 1 ]   =   B [ 3 ] [ 3 * i + 2 ]  =   B [ 5 ] [ 3 * i ]     =   phiG [ 1 ] [ i ];
        B [ 2 ] [ 3 * i + 2 ]   =   B [ 3 ] [ 3 * i + 1 ]  =   B [ 4 ] [ 3 * i ]     =   phiG [ 2 ] [ i ];
    }
    return B;
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
}


//////////////////////////////////////////////////////////
Matrix CosseratQuad :: giveBMatrix(const Point *x) const {
    Matrix phiG(ndim, numOfNodes);
    shafunc->giveShapeFGrad(x, nodes, phiG);
    Vector phi(nodes.size() );
    shafunc->giveShapeF(x, phi);
    Matrix B(6, DoFids.size() );
    for ( unsigned i = 0; i < numOfNodes; i++ ) {
        //strains
        // 00 03
        // 02 01
        B [ 0 ] [ 3 * i ] = B [ 3 ] [ 3 * i + 1 ] =  phiG [ 0 ] [ i ];
        B [ 2 ] [ 3 * i ] = B [ 1 ] [ 3 * i + 1 ] =  phiG [ 1 ] [ i ];
        //effect of rotation = - LeviCivita * rotation
        //  0 -1
        //  1  0
        B [ 3 ] [ 3 * i + 2 ] = -phi [ i ];
        B [ 2 ] [ 3 * i + 2 ] =  phi [ i ];
        //curvatures
        // 04 05
        B [ 4 ] [ 3 * i + 2 ] =  phiG [ 0 ] [ i ];
        B [ 5 ] [ 3 * i + 2 ] =  phiG [ 1 ] [ i ];
    }
    return B;
}

//////////////////////////////////////////////////////////
Matrix CosseratQuad :: giveHMatrix(const Point *x) const {
    Vector phi(nodes.size() );
    shafunc->giveShapeF(x, phi);
    Matrix H(3, DoFids.size() ); //2 transl, 1 rot
    for ( unsigned j = 0; j < numOfNodes; j++ ) {
        H [ 0 ] [ 3 * j ] = H [ 1 ] [ 3 * j + 1 ] = H [ 1 ] [ 3 * j + 2 ] = phi [ j ];
    }
    return H;
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
}


//////////////////////////////////////////////////////////
Matrix CosseratBrick :: giveBMatrix(const Point *x) const {
    Matrix phiG(ndim, numOfNodes);
    shafunc->giveShapeFGrad(x, nodes, phiG);
    Vector phi(nodes.size() );
    shafunc->giveShapeF(x, phi);
    Matrix B(18, DoFids.size() ); //9 strains, 9 curvatures
    for ( unsigned i = 0; i < numOfNodes; i++ ) {
        //strains
        // 00 08 06
        // 07 01 04
        // 05 03 02
        B [ 0 ] [ 6 * i ] = B [ 8 ] [ 6 * i + 1 ] = B [ 6 ] [ 6 * i + 2 ] = phiG [ 0 ] [ i ];
        B [ 7 ] [ 6 * i ] = B [ 1 ] [ 6 * i + 1 ] = B [ 4 ] [ 6 * i + 2 ] = phiG [ 1 ] [ i ];
        B [ 5 ] [ 6 * i ] = B [ 3 ] [ 6 * i + 1 ] = B [ 2 ] [ 6 * i + 2 ] = phiG [ 2 ] [ i ];
        //effect of rotation = - LeviCivita * rotations
        //  0 -Z  Y
        //  Z  0 -X
        // -Y  X  0
        B [ 8 ] [ 6 * i + 5 ] = B [ 4 ] [ 6 * i + 3 ] = B [ 5 ] [ 6 * i + 4 ] = -phi [ i ];
        B [ 6 ] [ 6 * i + 4 ] = B [ 7 ] [ 6 * i + 5 ] = B [ 3 ] [ 6 * i + 3 ] =  phi [ i ];
        //curvatures
        // 09 17 15
        // 16 10 13
        // 14 12 11
        B [  9 ] [ 6 * i + 3 ] = B [ 17 ] [ 6 * i + 4 ] = B [ 15 ] [ 6 * i + 5 ] = phiG [ 0 ] [ i ];
        B [ 16 ] [ 6 * i + 3 ] = B [ 10 ] [ 6 * i + 4 ] = B [ 13 ] [ 6 * i + 5 ] = phiG [ 1 ] [ i ];
        B [ 14 ] [ 6 * i + 3 ] = B [ 12 ] [ 6 * i + 4 ] = B [ 11 ] [ 6 * i + 5 ] = phiG [ 2 ] [ i ];
    }
    return B;
}

//////////////////////////////////////////////////////////
Matrix CosseratBrick :: giveHMatrix(const Point *x) const {
    Vector phi(nodes.size() );
    shafunc->giveShapeF(x, phi);
    Matrix H(6, DoFids.size() ); //3 transl, 3 rot
    for ( unsigned v = 0; v < ndim; v++ ) {
        for ( unsigned j = 0; j < numOfNodes; j++ ) {
            H [ v ] [ 6 * j + v ] = H [ 3 + v ] [ 6 * j + 3 + v ] = phi [ j ];
        }
    }
    return H;
}

//////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////
// 2D QUADRILATERAL COSSERAT COUPLED MECHANICAL-TRANSPORT ELEMENT
CoupledCosseratQuad :: CoupledCosseratQuad() {
    ndim = 2;
    name = "CoupledCosseratQuad";
    numOfNodes = 4;
    vtk_cell_type = 9;
    shafunc = new Linear2DQuadShapeF();
    inttype = new IntegrQuad4();
}


//////////////////////////////////////////////////////////
Matrix CoupledCosseratQuad :: giveBMatrix(const Point *x) const {
    Matrix phiG(ndim, numOfNodes);
    shafunc->giveShapeFGrad(x, nodes, phiG);
    Vector phi(nodes.size() );
    shafunc->giveShapeF(x, phi);
    Matrix B(8, DoFids.size() ); //4 strains, 2 curvatures, 2 pressure gradients
    for ( unsigned i = 0; i < numOfNodes; i++ ) {
        //strains
        // 00 03
        // 02 01
        B [ 0 ] [ 4 * i ] = B [ 3 ] [ 4 * i + 1 ] =  phiG [ 0 ] [ i ];
        B [ 2 ] [ 4 * i ] = B [ 1 ] [ 4 * i + 1 ] =  phiG [ 1 ] [ i ];
        //effect of rotation = - LeviCivita * rotation
        //  0 -1
        //  1  0
        B [ 3 ] [ 4 * i + 2 ] = -phi [ i ];
        B [ 2 ] [ 4 * i + 2 ] =  phi [ i ];
        //curvatures
        // 04 05
        B [ 4 ] [ 4 * i + 2 ] =  phiG [ 0 ] [ i ];
        B [ 5 ] [ 4 * i + 2 ] =  phiG [ 1 ] [ i ];
        //pressure gradient
        // 06 07
        B [ 6 ] [ 4 * i + 3 ] =  phiG [ 0 ] [ i ];
        B [ 7 ] [ 4 * i + 3 ] =  phiG [ 1 ] [ i ];
    }
    return B;
}

//////////////////////////////////////////////////////////
Matrix CoupledCosseratQuad :: giveHMatrix(const Point *x) const {
    Vector phi(nodes.size() );
    shafunc->giveShapeF(x, phi);
    Matrix H(4, DoFids.size() ); //2 transl, 1 rot, 1 pressure
    for ( unsigned j = 0; j < numOfNodes; j++ ) {
        H [ 0 ] [ 4 * j ] = H [ 1 ] [ 4 * j + 1 ] = H [ 2 ] [ 4 * j + 2 ] = H [ 3 ] [ 4 * j + 3 ] = phi [ j ];
    }
    return H;
}


//////////////////////////////////////////////////////////
void CoupledCosseratQuad :: init() {
    MechanicalElement :: init();

    CoupledParticle *cp;
    for ( auto &n: nodes ) {
        cp = dynamic_cast< CoupledParticle * >( n );
        if ( cp == nullptr ) {
            cerr << name << " requires nodes to be inhereted from CoupledParticle" << endl;
            exit(1);
        }
    }
}

//////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////
// 3D BRICK COSSERAT COUPLED MECHANICAL-TRANSPORT ELEMENT
CoupledCosseratBrick :: CoupledCosseratBrick() {
    ndim = 3;
    name = "CoupledCosseratBrick";
    numOfNodes = 8;
    vtk_cell_type = 12;
    shafunc = new Linear3DBrickShapeF();
    inttype = new IntegrBrick8();
}

//////////////////////////////////////////////////////////
Vector CoupledCosseratBrick :: giveStrain(unsigned i, const Vector &DoFs){

    Vector strain = CosseratBrick::giveStrain(i,DoFs);

    //save pressure at integration point for material model
    IPpressures.resize(numOfNodes);
    Vector IPvalues = giveHMatrix(inttype->giveIPLocationPointer(i))*DoFs; 
    IPpressures[i] = IPvalues[6];

    //save volumetric strain at integration point for material model
    IPvolstrains.resize(numOfNodes);
    IPvolstrains[i] = (strain[0] + strain[1] + strain[2])/3.;

    return strain;
}

//////////////////////////////////////////////////////////
Vector CoupledCosseratBrick :: giveInternalForces(const Vector &DoFs, bool frozen, double timeStep) {
    Vector intF = CosseratBrick :: giveInternalForces(DoFs, frozen, timeStep);

    
    //Biot effect
    DiscreteCoupledRVEMaterialStatus *couplm;
    double volumetricBiotPart;
    Matrix H;
    for ( unsigned k = 0; k < inttype->giveNumIP(); k++ ) {
        couplm = static_cast< DiscreteCoupledRVEMaterialStatus * >( stats [ k ] );
        volumetricBiotPart = couplm->computeBiotEffect();
        cout << volumetricBiotPart << endl;
        H = giveHMatrix(inttype->giveIPLocationPointer(k));
        for ( unsigned j = 0; j < 7*8; j++ ) intF[j] -= H[6][j]*volumetricBiotPart*inttype->giveIPWeight(k);
    }

    return intF;
}

//////////////////////////////////////////////////////////
Matrix CoupledCosseratBrick :: giveBMatrix(const Point *x) const {
    Matrix phiG(ndim, numOfNodes);
    shafunc->giveShapeFGrad(x, nodes, phiG);
    Vector phi(nodes.size() );
    shafunc->giveShapeF(x, phi);
    Matrix B(21, DoFids.size() ); //9 strains, 9 curvatures, 3 pressure gradients
    for ( unsigned i = 0; i < numOfNodes; i++ ) {
        //strains
        // 00 08 06
        // 07 01 04
        // 05 03 02
        B [ 0 ] [ 7 * i ] = B [ 8 ] [ 7 * i + 1 ] = B [ 6 ] [ 7 * i + 2 ] = phiG [ 0 ] [ i ];
        B [ 7 ] [ 7 * i ] = B [ 1 ] [ 7 * i + 1 ] = B [ 4 ] [ 7 * i + 2 ] = phiG [ 1 ] [ i ];
        B [ 5 ] [ 7 * i ] = B [ 3 ] [ 7 * i + 1 ] = B [ 2 ] [ 7 * i + 2 ] = phiG [ 2 ] [ i ];
        //effect of rotation = - LeviCivita * rotations
        //  0 -Z  Y
        //  Z  0 -X
        // -Y  X  0
        B [ 8 ] [ 7 * i + 5 ] = B [ 4 ] [ 7 * i + 3 ] = B [ 5 ] [ 7 * i + 4 ] = -phi [ i ];
        B [ 6 ] [ 7 * i + 4 ] = B [ 7 ] [ 7 * i + 5 ] = B [ 3 ] [ 7 * i + 3 ] =  phi [ i ];
        //curvatures
        // 09 17 15
        // 16 10 13
        // 14 12 11
        B [  9 ] [ 7 * i + 3 ] = B [ 17 ] [ 7 * i + 4 ] = B [ 15 ] [ 7 * i + 5 ] = phiG [ 0 ] [ i ];
        B [ 16 ] [ 7 * i + 3 ] = B [ 10 ] [ 7 * i + 4 ] = B [ 13 ] [ 7 * i + 5 ] = phiG [ 1 ] [ i ];
        B [ 14 ] [ 7 * i + 3 ] = B [ 12 ] [ 7 * i + 4 ] = B [ 11 ] [ 7 * i + 5 ] = phiG [ 2 ] [ i ];
        //pressure gradient
        // 18 19 20
        B [ 18 ] [ 7 * i + 6 ] = phiG [ 0 ] [ i ];
        B [ 19 ] [ 7 * i + 6 ] = phiG [ 1 ] [ i ];
        B [ 20 ] [ 7 * i + 6 ] = phiG [ 2 ] [ i ];
    }
    return B;
}

//////////////////////////////////////////////////////////
Matrix CoupledCosseratBrick :: giveHMatrix(const Point *x) const {
    Vector phi(nodes.size() );
    shafunc->giveShapeF(x, phi);
    Matrix H(7, DoFids.size() ); //3 transl, 3 rot, 1 pressure
    for ( unsigned j = 0; j < numOfNodes; j++ ) {
        for ( unsigned v = 0; v < ndim; v++ ) {
            H [ v ] [ 7 * j + v ] = H [ 3 + v ] [ 7 * j + 3 + v ] = phi [ j ];
        }
        H [ 6 ] [ 7 * j + 6 ] = phi [ j ]; //pressure gradient
    }
    return H;
}

//////////////////////////////////////////////////////////
Matrix CoupledCosseratBrick :: giveDampingMatrix() const {    
    //return giveStiffnessMatrix("elastic") * 1e-15;           //rough fix of zeros, here can be anything
    return Element :: giveDampingMatrix();
}
