#include "element_ldpm.h"
#include "element_container.h"
#include "boundary_condition.h"
#include "material_vectorial.h"

using namespace std;

//////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////
// RBSN ELEMENT
LDPMTetra :: LDPMTetra(unsigned dim) : Element{dim} {
    cout << "LDPMTetra :: LDPMTetra(unsigned dim)" << endl; cout.flush();

    if ( ndim != 3 ) {
        cerr << "LDPMTetra implemented only in 3D" << endl;
        exit(1);
    }
    numOfNodes = 4;
    nodes.resize(4);
    name = "LDPMTetra";
    vtk_cell_type = 10;
    shafunc = new Linear1DLineShapeF();
    inttype = new IntegrLDPM12();

    vert.resize(12);
    lengths.resize(12);
    volumes.resize(12);
    areas.resize(12);
    normals.resize(12);
    R.resize(12);

    nodecodes = { 1, 2, 1, 3, 2, 3, 0, 2, 0, 3, 2, 3, 0, 1, 0, 3, 1, 3, 0, 1, 0, 2, 1, 2};
    vertcodes = { 1, 2, 3, 1, 1, 4, 6, 5, 5, 7, 4, 5, 8, 9, 7, 8, 8, 3, 9, 10, 10, 6, 2, 10}; //last point is always centroid at position 0
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
        Point n = ( vert [ vertcodes [ 2 * i ] ]->givePoint() - vert [ vertcodes [ 2 * i + 1] ]->givePoint() ).cross( vert [ 0 ]->givePoint() - vert [ vertcodes [ 2 * i ] ]->givePoint() );
        n /= n.norm();
        //contact vector
        normals [ i ] = nodes [ nodecodes [ 2 * i + 1 ] ]->givePoint() - nodes [ nodecodes [ 2 * i ] ]->givePoint();
        lengths [ i ] = normals [ i ].norm();
        normals [ i ] /= lengths [ i ];
        inttype->setIPLocation(i, ( vert [ vertcodes [ 2 * i  ] ]->givePoint() + vert [ vertcodes [ 2 * i + 1 ] ]->givePoint() + vert [ 0 ]->givePoint() ) / 3.);
        areas [ i ] = triArea3D(vert [ vertcodes [ 2 * i  ] ]->givePointPointer(), vert [ vertcodes [ 2 * i + 1 ] ]->givePointPointer(), vert [ 0 ]->givePointPointer() );
        areas [ i ] *= n.dot(normals [ i ]); //projection of area

        Point t1, t2;
        t1 = vert [ 0 ]->givePoint() - inttype->giveIPLocation(i);                
        t1.normalize();
        t2 = normals [ i ].cross(t1);
        R [ i ] = Matrix :: Zero(3, 3);
        R [ i ](0, 0) = normals [ i ].x();
        R [ i ](0, 1) = normals [ i ].y();
        R [ i ](0, 2) = normals [ i ].z();
        R [ i ](1, 0) = t1.x();
        R [ i ](1, 1) = t1.y();
        R [ i ](1, 2) = t1.z();
        R [ i ](2, 0) = t2.x();
        R [ i ](2, 1) = t2.y();
        R [ i ](2, 2) = t2.z();

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
    if ( i == 0 ) {
        volumetricStrain = 0;
    }
    stats [ i ]->setParameterValue("volumetric_strain", volumetricStrain);

    return Element :: giveStrain(i, DoFs);
};

//////////////////////////////////////////////////////////
Matrix LDPMTetra :: giveStiffnessMatrix(string matrixType) const {
    return Element :: giveStiffnessMatrix(matrixType) * ndim; //ndim needs to be included here for discrete elements
}

//////////////////////////////////////////////////////////

Matrix LDPMTetra :: giveDampingMatrix() const {
    return giveStiffnessMatrix("elastic") * 1e-15;           //rough fix of zeros, here can be anything
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
Matrix LDPMTetra :: giveMassMatrix() const {
    Matrix M = Matrix :: Zero(24, 24);
    return M;
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

/*
 * //////////////////////////////////////////////////////////
 * void RigidBodyContact :: giveValues(string code, Vector &result) const {
 *  if ( code.compare("damage") == 0 ) {
 *      stats [ 0 ]->giveValues(code, result);
 *  } else if ( code.compare("normal") == 0 ) {
 *      result.resize(ndim);
 *      for ( unsigned i = 0; i < ndim; i++ ) {
 *          result [ i ] = normal(i);
 *      }
 *  } else if ( code.compare("t1") == 0 ) {
 *      result.resize(ndim);
 *      for ( unsigned i = 0; i < ndim; i++ ) {
 *          result [ i ] = R(1, i);
 *      }
 *
 *  } else if ( code.compare("t2") == 0 ) {
 *      result.resize(ndim);
 *      for ( unsigned i = 0; i < ndim; i++ ) {
 *          result [ i ] = R(1, i);
 *      }
 *  } else if ( code.compare("volume") == 0 ) {
 *      result.resize(1);
 *      result [ 0 ] = area * length / ndim;
 *  } else if ( code.compare("area") == 0 ) {
 *      result.resize(1);
 *      result [ 0 ] = area;
 *  } else if ( code.compare("length") == 0 ) {
 *      result.resize(1);
 *      result [ 0 ] = length;
 *  } else {
 *      MechanicalElement :: giveValues(code, result);
 *  }
 * }
 *
 *
 *
 *
 *
 * //////////////////////////////////////////////////////////
 * Matrix RigidBodyContact::giveMassMatrix() const {
 *  Matrix M = Matrix::Zero(6 * (ndim - 1), 6 * (ndim - 1));
 *  VectMechMaterialStatus *mechstat = static_cast<VectMechMaterialStatus *>(stats[0]);
 *  double density = mechstat->giveDensity();
 *  double m0 = giveVolumeAssociatedWithNode(0) * density; ///mass
 *  double m1 = giveVolumeAssociatedWithNode(1) * density; ///mass
 *  if (ndim == 2) {
 *          // Define points
 *          Point *A = nodes[0]->givePointPointer();
 *          Point *B = nodes[1]->givePointPointer();
 *          Point *C = vert[0]->givePointPointer();
 *          Point *D = vert[1]->givePointPointer();
 *          Point cg0 = ( (*A) + (*C) + (*D) ) / 3.;
 *          Point cg1 = ( (*B) + (*C) + (*D) ) / 3.;
 *          Point null(0,0,0);
 *          Point C_ = (*C) - (*A);
 *          Point D_ = (*D) - (*A);
 *          Point C__ = (*C) - (*B);
 *          Point D__ = (*D) - (*B);
 *          // MassMatrix
 *          M(0, 0) = M(1, 1) = m0;
 *          M(3, 3) = M(4, 4) = m1;
 *          M(2, 2) = density * ( triInertia2D(&null, &C_, &D_) ); // Inertia relative to the node A[0,0]
 *          M(5, 5) = density * ( triInertia2D(&null, &C__, &D__) ); // Inertia relative to the node B[0,0]
 *          M(0, 2) = M(2, 0) = - m0 * ( cg0.y() - A->y() );
 *          M(1, 2) = M(2, 1) = + m0 * ( cg0.x() - A->x() );
 *          M(3, 5) = M(5, 3) = - m1 * ( cg1.y() - B->y() );
 *          M(4, 5) = M(5, 4) = + m1 * ( cg1.x() - B->x() );
 *  }
 *  else if (ndim == 3) {
 *          // Define points
 *          Point *A = nodes[0]->givePointPointer();
 *          Point *B = nodes[1]->givePointPointer();
 *          Point *FaceCentroid = inttype->giveIPLocationPointer(0);
 *          // Mass
 *          M(0, 0) = M(1, 1) = M(2, 2) = m0;
 *          M(6, 6) = M(7, 7) = M(8, 8) = m1;
 *          size_t n = vert.size();
 *          for ( unsigned i = 0 ; i < n ; i++ ) {
 *              Point *C = vert[i]->givePointPointer();
 *              Point *D;
 *              if (i == 0) {
 *              D = vert[n-1]->givePointPointer();
 *              } else {
 *              D = vert[i-1]->givePointPointer();
 *              }
 *              Point cg0 = ( (*A) + (*C) + (*D) + (*FaceCentroid) ) / 4.;
 *              Point cg1 = ( (*B) + (*C) + (*D) + (*FaceCentroid) ) / 4.;
 *          // Volume of tetrahedron
 *              Point A_ = (*A) - (*FaceCentroid);
 *              Point B_ = (*B) - (*FaceCentroid);
 *              Point C_ = (*C) - (*FaceCentroid);
 *              Point D_ = (*D) - (*FaceCentroid);
 *              double tetraVolume0 = A_.dot( C_.cross( D_ ) ) / 6.;
 *              double tetraVolume1 = B_.dot( C_.cross( D_ ) ) / 6.;
 * //cout << tetraVolume0 << endl;
 *          // Inertia matrix relative to the centroid [0,0,0]
 *              Point A__ = (*A) - cg0;
 *              Point C__ = (*C) - cg0;
 *              Point D__ = (*D) - cg0;
 *              Point FaceCentroid__ = (*FaceCentroid) - cg0;
 *              Point B___ = (*B) - cg1;
 *              Point C___ = (*C) - cg1;
 *              Point D___ = (*D) - cg1;
 *              Point FaceCentroid___ = (*FaceCentroid) - cg1;
 *              Matrix Icg0 = tetraInertia3D(&A__, &C__, &D__, &FaceCentroid__);
 *              Matrix Icg1 = tetraInertia3D(&B___, &C___, &D___, &FaceCentroid___);
 *
 *          // MassMatrix
 *              M(3, 3) += density * ( Icg0(0, 0) + tetraVolume0 * ( pow((cg0.y() - A->y()),2) + pow((cg0.z() - A->z()),2) ) );
 *              M(4, 4) += density * ( Icg0(1, 1) + tetraVolume0 * ( pow((cg0.x() - A->x()),2) + pow((cg0.z() - A->z()),2) ) );
 *              M(5, 5) += density * ( Icg0(2, 2) + tetraVolume0 * ( pow((cg0.x() - A->x()),2) + pow((cg0.y() - A->y()),2) ) );
 *              M(3, 4) = M(4, 3) += density * ( Icg0(0, 1) - tetraVolume0 * ( (cg0.x() - A->x()) * (cg0.y() - A->y()) ) );
 *              M(3, 5) = M(5, 3) += density * ( Icg0(0, 2) - tetraVolume0 * ( (cg0.x() - A->x()) * (cg0.z() - A->z()) ) );
 *              M(4, 5) = M(5, 4) += density * ( Icg0(1, 2) - tetraVolume0 * ( (cg0.y() - A->y()) * (cg0.z() - A->z()) ) );
 *
 *              M(9, 9) += density * ( Icg1(0, 0) + tetraVolume1 * ( pow((cg1.y() - B->y()),2) + pow((cg1.z() - B->z()),2) ) );
 *              M(10, 10) += density * ( Icg1(1, 1) + tetraVolume1 * ( pow((cg1.x() - B->x()),2) + pow((cg1.z() - B->z()),2) ) );
 *              M(11, 11) += density * ( Icg1(2, 2) + tetraVolume1 * ( pow((cg1.x() - B->x()),2) + pow((cg1.y() - B->y()),2) ) );
 *              M(9, 10) = M(10, 9) += density * ( Icg1(0, 1) - tetraVolume1 * ( (cg1.x() - B->x()) * (cg1.y() - B->y()) ) );
 *              M(9, 11) = M(11, 9) += density * ( Icg1(0, 2) - tetraVolume1 * ( (cg1.x() - B->x()) * (cg1.z() - B->z()) ) );
 *              M(10, 11) = M(11, 10) += density * ( Icg1(1, 2) - tetraVolume1 * ( (cg1.y() - B->y()) * (cg1.z() - B->z()) ) );
 *
 *              M(0, 4) = M(4, 0) += + ( tetraVolume0 * density ) * ( cg0.z() - A->z() );
 *              M(0, 5) = M(5, 0) += - ( tetraVolume0 * density ) * ( cg0.y() - A->y() );
 *              M(1, 3) = M(3, 1) += - ( tetraVolume0 * density ) * ( cg0.z() - A->z() );
 *              M(1, 5) = M(5, 1) += + ( tetraVolume0 * density ) * ( cg0.x() - A->x() );
 *              M(2, 3) = M(3, 2) += + ( tetraVolume0 * density ) * ( cg0.y() - A->y() );
 *              M(2, 4) = M(4, 2) += - ( tetraVolume0 * density ) * ( cg0.x() - A->x() );
 *
 *              M(6, 10) = M(10, 6) += + ( tetraVolume1 * density ) * ( cg0.z() - B->z() );
 *              M(6, 11) = M(11, 6) += - ( tetraVolume1 * density ) * ( cg0.y() - B->y() );
 *              M(7, 9) = M(9, 7) += - ( tetraVolume1 * density ) * ( cg0.z() - B->z() );
 *              M(7, 11) = M(11, 7) += + ( tetraVolume1 * density ) * ( cg0.x() - B->x() );
 *              M(8, 9) = M(9, 8) += + ( tetraVolume1 * density ) * ( cg0.y() - B->y() );
 *              M(8, 10) = M(10, 8) += - ( tetraVolume1 * density ) * ( cg0.x() - B->x() );
 *          }
 *  }
 * //    cout << M << endl;
 * //    exit(1);
 *  return M;
 * }
 *
 *
 *
 *
 *
 *
 * //////////////////////////////////////////////////////////
 * Matrix RigidBodyContact::giveMassMatrix() const {
 *  Matrix M = Matrix::Zero(6 * (ndim - 1), 6 * (ndim - 1));
 *  VectMechMaterialStatus *mechstat = static_cast<VectMechMaterialStatus *>(stats[0]);
 *  double density = mechstat->giveDensity();
 *  double m0 = giveVolumeAssociatedWithNode(0) * density; ///mass
 *  double m1 = giveVolumeAssociatedWithNode(1) * density; ///mass
 *  if (ndim == 2) {
 *          // Define points
 *          Point *A = nodes[0]->givePointPointer();
 *          Point *B = nodes[1]->givePointPointer();
 *          Point *C = vert[0]->givePointPointer();
 *          Point *D = vert[1]->givePointPointer();
 *          Point cg0 = ( (*A) + (*C) + (*D) ) / 3.;
 *          Point cg1 = ( (*B) + (*C) + (*D) ) / 3.;
 *          Point null(0,0,0);
 *          Point C_ = (*C) - (*A);
 *          Point D_ = (*D) - (*A);
 *          Point C__ = (*C) - (*B);
 *          Point D__ = (*D) - (*B);
 *          // MassMatrix
 *          M(0, 0) = M(1, 1) = m0;
 *          M(3, 3) = M(4, 4) = m1;
 *          M(2, 2) = density * ( triInertia2D(&null, &C_, &D_) ); // Inertia relative to the node A[0,0]
 *          M(5, 5) = density * ( triInertia2D(&null, &C__, &D__) ); // Inertia relative to the node B[0,0]
 *          M(0, 2) = M(2, 0) = - m0 * ( cg0.y() - A->y() );
 *          M(1, 2) = M(2, 1) = + m0 * ( cg0.x() - A->x() );
 *          M(3, 5) = M(5, 3) = - m1 * ( cg1.y() - B->y() );
 *          M(4, 5) = M(5, 4) = + m1 * ( cg1.x() - B->x() );
 *  }
 *  else if (ndim == 3) {
 *          // Define points
 *          Point *A = nodes[0]->givePointPointer();
 *          Point *B = nodes[1]->givePointPointer();
 *          Point *FaceCentroid = inttype->giveIPLocationPointer(0);
 *          // Mass
 *          M(0, 0) = M(1, 1) = M(2, 2) = m0;
 *          M(6, 6) = M(7, 7) = M(8, 8) = m1;
 *          size_t n = vert.size();
 *          for ( unsigned i = 0 ; i < n ; i++ ) {
 *              Point *C = vert[i]->givePointPointer();
 *              Point *D;
 *              if (i == 0) {
 *              D = vert[n-1]->givePointPointer();
 *              } else {
 *              D = vert[i-1]->givePointPointer();
 *              }
 *              Point cg0 = ( (*A) + (*C) + (*D) + (*FaceCentroid) ) / 4.;
 *              Point cg1 = ( (*B) + (*C) + (*D) + (*FaceCentroid) ) / 4.;
 *          // Volume of tetrahedron
 *              Point A_ = (*A) - (*FaceCentroid);
 *              Point B_ = (*B) - (*FaceCentroid);
 *              Point C_ = (*C) - (*FaceCentroid);
 *              Point D_ = (*D) - (*FaceCentroid);
 *              double tetraVolume0 = A_.dot( C_.cross( D_ ) ) / 6.;
 *              double tetraVolume1 = B_.dot( C_.cross( D_ ) ) / 6.;
 * //cout << tetraVolume0 << endl;
 *          // Inertia matrix relative to the centroid [0,0,0]
 *              Point A__ = (*A) - cg0;
 *              Point C__ = (*C) - cg0;
 *              Point D__ = (*D) - cg0;
 *              Point FaceCentroid__ = (*FaceCentroid) - cg0;
 *              Point B___ = (*B) - cg1;
 *              Point C___ = (*C) - cg1;
 *              Point D___ = (*D) - cg1;
 *              Point FaceCentroid___ = (*FaceCentroid) - cg1;
 *              Matrix Icg0 = tetraInertia3D(&A__, &C__, &D__, &FaceCentroid__);
 *              Matrix Icg1 = tetraInertia3D(&B___, &C___, &D___, &FaceCentroid___);
 *
 *          // MassMatrix
 *              M(3, 3) += density * ( Icg0(0, 0) + tetraVolume0 * ( pow((cg0.y() - A->y()),2) + pow((cg0.z() - A->z()),2) ) );
 *              M(4, 4) += density * ( Icg0(1, 1) + tetraVolume0 * ( pow((cg0.x() - A->x()),2) + pow((cg0.z() - A->z()),2) ) );
 *              M(5, 5) += density * ( Icg0(2, 2) + tetraVolume0 * ( pow((cg0.x() - A->x()),2) + pow((cg0.y() - A->y()),2) ) );
 *              M(3, 4) = M(4, 3) += density * ( Icg0(0, 1) - tetraVolume0 * ( (cg0.x() - A->x()) * (cg0.y() - A->y()) ) );
 *              M(3, 5) = M(5, 3) += density * ( Icg0(0, 2) - tetraVolume0 * ( (cg0.x() - A->x()) * (cg0.z() - A->z()) ) );
 *              M(4, 5) = M(5, 4) += density * ( Icg0(1, 2) - tetraVolume0 * ( (cg0.y() - A->y()) * (cg0.z() - A->z()) ) );
 *
 *              M(9, 9) += density * ( Icg1(0, 0) + tetraVolume1 * ( pow((cg1.y() - B->y()),2) + pow((cg1.z() - B->z()),2) ) );
 *              M(10, 10) += density * ( Icg1(1, 1) + tetraVolume1 * ( pow((cg1.x() - B->x()),2) + pow((cg1.z() - B->z()),2) ) );
 *              M(11, 11) += density * ( Icg1(2, 2) + tetraVolume1 * ( pow((cg1.x() - B->x()),2) + pow((cg1.y() - B->y()),2) ) );
 *              M(9, 10) = M(10, 9) += density * ( Icg1(0, 1) - tetraVolume1 * ( (cg1.x() - B->x()) * (cg1.y() - B->y()) ) );
 *              M(9, 11) = M(11, 9) += density * ( Icg1(0, 2) - tetraVolume1 * ( (cg1.x() - B->x()) * (cg1.z() - B->z()) ) );
 *              M(10, 11) = M(11, 10) += density * ( Icg1(1, 2) - tetraVolume1 * ( (cg1.y() - B->y()) * (cg1.z() - B->z()) ) );
 *
 *              M(0, 4) = M(4, 0) += + ( tetraVolume0 * density ) * ( cg0.z() - A->z() );
 *              M(0, 5) = M(5, 0) += - ( tetraVolume0 * density ) * ( cg0.y() - A->y() );
 *              M(1, 3) = M(3, 1) += - ( tetraVolume0 * density ) * ( cg0.z() - A->z() );
 *              M(1, 5) = M(5, 1) += + ( tetraVolume0 * density ) * ( cg0.x() - A->x() );
 *              M(2, 3) = M(3, 2) += + ( tetraVolume0 * density ) * ( cg0.y() - A->y() );
 *              M(2, 4) = M(4, 2) += - ( tetraVolume0 * density ) * ( cg0.x() - A->x() );
 *
 *              M(6, 10) = M(10, 6) += + ( tetraVolume1 * density ) * ( cg0.z() - B->z() );
 *              M(6, 11) = M(11, 6) += - ( tetraVolume1 * density ) * ( cg0.y() - B->y() );
 *              M(7, 9) = M(9, 7) += - ( tetraVolume1 * density ) * ( cg0.z() - B->z() );
 *              M(7, 11) = M(11, 7) += + ( tetraVolume1 * density ) * ( cg0.x() - B->x() );
 *              M(8, 9) = M(9, 8) += + ( tetraVolume1 * density ) * ( cg0.y() - B->y() );
 *              M(8, 10) = M(10, 8) += - ( tetraVolume1 * density ) * ( cg0.x() - B->x() );
 *          }
 *  }
 * //    cout << M << endl;
 * //    exit(1);
 *  return M;
 * }
 *
 * //////////////////////////////////////////////////////////
 * Vector RigidBodyContact :: giveContactStrainNT() const {
 *  return stats [ 0 ]->giveTempStrain();
 * };
 *
 * //////////////////////////////////////////////////////////
 * Vector RigidBodyContact :: giveContactStrainXYZ() const {
 *  return this->R.transpose() * this->giveContactStrainNT();
 * };
 *
 * //////////////////////////////////////////////////////////
 * Vector RigidBodyContact :: giveContactStressXYZ() {
 *  return this->R.transpose() * stats [ 0 ]->giveTempStress();
 * };
 *
 * //////////////////////////////////////////////////////////
 * Vector RigidBodyContact :: transformVectorToXYZ(Vector &result) const {
 *  return this->R.transpose() * result;
 * };
 *
 * //////////////////////////////////////////////////////////
 * Vector RigidBodyContact :: transformToGlobal(const Vector &DoFs) const {
 *  return this->R * DoFs;
 * }
 *
 * //////////////////////////////////////////////////////////
 * Vector RigidBodyContact :: giveVectorToNode(const unsigned &node_i, const unsigned &ip_id) const {
 *  ( void ) ip_id;
 *  Point distance = inttype->giveIPLocation(0) - nodes [ node_i ]->givePoint();
 *  Vector dst = Vector :: Zero(ndim);
 *  for ( unsigned i = 0; i < ndim; i++ ) {
 *      if ( i == 0 ) {
 *          dst [ i ] = distance.x();
 *      } else if ( i == 1 ) {
 *          dst [ i ] = distance.y();
 *      } else if ( i == 2 ) {
 *          dst [ i ] = distance.z();
 *      }
 *  }
 *  return dst;
 * }
 *
 * //////////////////////////////////////////////////////////
 * double RigidBodyContact :: giveVolumeAssociatedWithNode(unsigned nodenum) const {
 *  if ( nodenum == 0 ) {
 *      return ( vert [ 0 ]->givePoint() - nodes [ 0 ]->givePoint() ).dot(normal) * area / ndim;
 *  } else if ( nodenum == 1 ) {
 *      return -( vert [ 0 ]->givePoint() - nodes [ 1 ]->givePoint() ).dot(normal) * area / ndim;
 *  } else {
 *      cerr << "Error in " << name << ": attempting to reach node number different form 0 or 1." << endl;
 *      exit(1);
 *  }
 * };
 *
 *
 *
 *
 *
 *
 *
 * //////////////////////////////////////////////////////////
 * //////////////////////////////////////////////////////////
 * void RigidBodyContact :: extrapolateIPValuesToNodes(string code, vector< Vector > &result, Vector &weights) const {
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
