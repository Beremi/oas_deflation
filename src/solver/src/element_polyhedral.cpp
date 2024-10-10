#include "element_polyhedral.h"
#include "material_tensorial.h"

using namespace std;

//////////////////////////////////////////////////////////
// 2 D
//////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////
// TRANSPORT POLYHEDRAL ELEMENT
TranspPolygonal :: TranspPolygonal(const unsigned dim) : Element(dim) {
    name = "TranspPolygonal";
    ip_type = "quad";
    numOfNodes = 0;
    vtk_cell_type = 7;
    shafunc = new Wachspress2DShapeF();
    inttype = new IntegrPolygon(ip_type);
    physicalFields [ 1 ] = true; //transport
}

//////////////////////////////////////////////////////////
void TranspPolygonal :: readFromLine(istringstream &iss, NodeContainer *fullnodes, MaterialContainer *fullmatrs) {
    unsigned num, num2;

    iss >> numOfNodes;
    nodes.resize(numOfNodes);
    for ( unsigned i = 0; i < numOfNodes; i++ ) {
        iss >> num2;
        nodes [ i ] = fullnodes->giveNode(num2);
    }
    iss >> num;
    mat = fullmatrs->giveMaterial(num);
}

//////////////////////////////////////////////////////////
void TranspPolygonal :: prepareGeometry() {
    //estimate centroid
    Point cpoint = Point(0, 0, 0);
    for ( const auto n: nodes ) {
        cpoint += n->givePoint();
    }
    cpoint /= numOfNodes;

    //compute angles
    vector< pair< double, unsigned > >angles;
    angles.resize(numOfNodes);
    for ( unsigned i = 0; i < numOfNodes; i++ ) {
        angles [ i ].second = i;
        angles [ i ].first = atan2(nodes [ i ]->givePoint().y() - cpoint.y(), nodes [ i ]->givePoint().x() - cpoint.x() );
    }

    //sort to have counterclockwise direction
    sort(angles.begin(), angles.end() );
    vector< Node * >newnodes;
    newnodes.resize(numOfNodes);
    for ( unsigned i = 0; i < numOfNodes; i++ ) {
        newnodes [ i ] = nodes [ angles [ i ].second ];
    }
    nodes = newnodes;


    volume = 0;
    double triVolume;
    faces.resize(numOfNodes);
    normals.resize(numOfNodes);
    surfaces.resize(numOfNodes);
    centroid = Point(0., 0., 0.);
    Point diff;
    //face X must start at node X and end at node X+1
    for ( unsigned i = 0; i < numOfNodes; i++ ) {
        faces [ i ].resize(2);
        faces [ i ] [ 0 ] = i;
        if ( i == numOfNodes - 1 ) {
            faces [ i ] [ 1 ] = 0;
        } else {
            faces [ i ] [ 1 ] = i + 1;
        }
        diff = nodes [ faces [ i ] [ 1 ] ]->givePoint() - nodes [ faces [ i ] [ 0 ] ]->givePoint();
        surfaces [ i ] = diff.norm();
        normals [ i ] = Point(diff.y() / surfaces [ i ], -diff.x() / surfaces [ i ], 0);
        triVolume = triArea2D(nodes [ faces [ i ] [ 0 ] ]->givePointPointer(), nodes [ faces [ i ] [ 1 ] ]->givePointPointer(), & cpoint);
        centroid += ( nodes [ faces [ i ] [ 0 ] ]->givePoint() + nodes [ faces [ i ] [ 1 ] ]->givePoint() + cpoint ) * triVolume;
        volume += triVolume;
    }
    centroid /= volume * 3.;
}

//////////////////////////////////////////////////////////
void TranspPolygonal :: initIntegration() {
    prepareGeometry();

    Wachspress2DShapeF *wsf = static_cast< Wachspress2DShapeF * >( shafunc );
    wsf->setFacesAndNormals(faces, normals);
    inttype->init(nodes, faces, & centroid);
    shafunc->init(nodes);
}


//////////////////////////////////////////////////////////
void TranspPolygonal :: setIntegrationPointsAndWeights() {
    stats.resize(inttype->giveNumIP() );
    for ( unsigned k = 0; k < inttype->giveNumIP(); k++ ) {
        stats [ k ] = mat->giveNewMaterialStatus(this, k);
    }
}

//////////////////////////////////////////////////////////
void TranspPolygonal :: init() {
    //check that nodes are TrsNodes
    for ( const auto n: nodes ) {
        TrsNode *p = dynamic_cast< TrsNode * >( n );
        if ( !p ) {
            cerr << "Error in " << name << ": nodes must be inherited from TrsNode, " << n->giveName() << " provided" << endl;
            exit(1);
        }
    }
    //check that material is VectMechMat
    TensTrsprtMaterial *p = dynamic_cast< TensTrsprtMaterial * >( mat );
    if ( !p ) {
        cerr << "Error in " << name << ": material must be inherited from TensTrsprtMaterial, " << mat->giveName() << " provided" << endl;
        exit(1);
    }

    Element :: init(); //calling base class method;
}

//////////////////////////////////////////////////////////
Matrix TranspPolygonal :: giveBMatrix(const Point *x) const {
    Matrix B = Matrix :: Zero(ndim, nodes.size() );
    shafunc->giveShapeFGrad(x, B);
    return B;
}

//////////////////////////////////////////////////////////
Matrix TranspPolygonal :: giveHMatrix(const Point *x) const {
    Vector X = Vector :: Zero(numOfNodes);
    shafunc->giveShapeF(x, X);
    Matrix H = Matrix :: Zero(1, numOfNodes);
    for ( unsigned k = 0; k < numOfNodes; k++ ) {
        H(0, k) = X [ k ];
    }
    return H;
}


//////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////
// TRANSPORT VIRTUAL POLYGONAL ELEMENT

TranspVirtPolygonal :: TranspVirtPolygonal(const unsigned dim) : TranspPolygonal(dim) {
    name = "TranspVirtPolygonal";
}

//////////////////////////////////////////////////////////
void TranspVirtPolygonal :: setIntegrationPointsAndWeights() {
    TranspPolygonal :: setIntegrationPointsAndWeights(); //calling base class method;

    normals.resize(numOfNodes);
    surfaces.resize(numOfNodes);
    Point diff;
    unsigned i, j, v;
    //face X must start at node X and end at node X+1
    for ( i = 0; i < numOfNodes; i++ ) {
        if ( i == numOfNodes - 1 ) {
            j = 0;
        } else {
            j = i + 1;
        }
        diff = nodes [ j ]->givePoint() - nodes [ i ]->givePoint();
        surfaces [ i ] = diff.norm();
        normals [ i ] = Point(diff.y() / surfaces [ i ], -diff.x() / surfaces [ i ], 0);
    }

    double radius = pow(volume / M_PI, 0.5);
    Matrix D = Matrix :: Zero(numOfNodes, ndim + 1);
    Point x;
    for ( i = 0; i < numOfNodes; i++ ) {
        D(i, 0) = 1.;
        x = nodes [ i ]->givePoint();
        for ( v = 0; v < ndim; v++ ) {
            D(i, v + 1) = ( x(v) - centroid(v) ) / radius;
        }
    }

    Matrix B = Matrix :: Zero(ndim + 1, numOfNodes);
    j = numOfNodes - 1;
    for (  i = 0; i < numOfNodes; i++ ) {
        B(0, i) = 1. / numOfNodes;
        for ( v = 0; v < ndim; v++ ) {
            B(v + 1, i) = ( normals [ i ](v) * surfaces [ i ] + normals [ j ](v) * surfaces [ j ] );
        }
        j = i;
    }
    B /= 2. * radius;

    Matrix G = B * D;
    Matrix Gtilde = G;
    for ( i = 0; i < ndim + 1; i++ ) {
        Gtilde(0, i) = 0.;
    }

    double Gdet = G(0, 0) * G(1, 1) * G(2, 2) + G(0, 2) * G(1, 0) * G(2, 1) + G(2, 0) * G(0, 1) * G(1, 2) - G(0, 0) * G(2, 1) * G(1, 2) - G(0, 1) * G(1, 0) * G(2, 2) - G(0, 2) * G(1, 1) * G(2, 0);
    Matrix Ginv = Matrix :: Zero(ndim + 1, ndim + 1);
    Ginv(0, 0) = G(1, 1) * G(2, 2) - G(1, 2) * G(2, 1);
    Ginv(1, 0) = -G(1, 0) * G(2, 2) - G(1, 2) * G(2, 0);
    Ginv(2, 0) = G(1, 0) * G(2, 1) - G(1, 1) * G(2, 0);
    Ginv(0, 1) = -G(0, 1) * G(2, 2) - G(0, 2) * G(2, 1);
    Ginv(1, 1) = G(0, 0) * G(2, 2) - G(0, 2) * G(2, 0);
    Ginv(2, 1) = -G(0, 0) * G(2, 1) - G(0, 1) * G(2, 0);
    Ginv(0, 2) = G(0, 1) * G(1, 2) - G(0, 2) * G(1, 1);
    Ginv(1, 2) = -G(0, 0) * G(1, 2) - G(0, 2) * G(1, 0);
    Ginv(2, 2) = G(0, 0) * G(1, 1) - G(0, 1) * G(1, 0);
    Ginv /= Gdet;

    V1 = Ginv * B;
    V2 = D * V1 * ( -1. );
    for ( i = 0; i < numOfNodes; i++ ) {
        V2(i, i) += 1.;
    }
    V1 = ( V1.transpose() * Gtilde ) * V1;

    Matrix H = Matrix :: Zero(ndim + 1, ndim + 1);
    Vector m = Vector :: Zero(ndim + 1);
    m [ 0 ] = 1;
    for ( i = 0; i < inttype->giveNumIP(); i++ ) {
        for ( v = 0; v < ndim; v++ ) {
            m [ v + 1 ] = ( inttype->giveIPLocation(i)(v) - centroid(v) ) / radius;
        }
        H += dyadicProduct( m, m * inttype->giveIPWeight(i) );
    }

    double Hdet = H(0, 0) * H(1, 1) * H(2, 2) + H(0, 2) * H(1, 0) * H(2, 1) + H(2, 0) * H(0, 1) * H(1, 2) - H(0, 0) * H(2, 1) * H(1, 2) - H(0, 1) * H(1, 0) * H(2, 2) - H(0, 2) * H(1, 1) * H(2, 0);
    Matrix Hinv = Matrix :: Zero(ndim + 1, ndim + 1);
    Hinv(0, 0) = H(1, 1) * H(2, 2) - H(1, 2) * H(2, 1);
    Hinv(1, 0) = -H(1, 0) * H(2, 2) - H(1, 2) * H(2, 0);
    Hinv(2, 0) = H(1, 0) * H(2, 1) - H(1, 1) * H(2, 0);
    Hinv(0, 1) = -H(0, 1) * H(2, 2) - H(0, 2) * H(2, 1);
    Hinv(1, 1) = H(0, 0) * H(2, 2) - H(0, 2) * H(2, 0);
    Hinv(2, 1) = -H(0, 0) * H(2, 1) - H(0, 1) * H(2, 0);
    Hinv(0, 2) = H(0, 1) * H(1, 2) - H(0, 2) * H(1, 1);
    Hinv(1, 2) = -H(0, 0) * H(1, 2) - H(0, 2) * H(1, 0);
    Hinv(2, 2) = H(0, 0) * H(1, 1) - H(0, 1) * H(1, 0);
    Hinv /= Hdet;

    W2 = Ginv * B;
    Matrix C = H * W2;
    W2 = ( D * W2 ) * ( -1. );
    for ( i = 0; i < numOfNodes; i++ ) {
        W2(i, i) += 1.;
    }
    W1 = ( C.transpose() * Hinv ) * C;
}

//////////////////////////////////////////////////////////
Matrix TranspVirtPolygonal :: giveStiffnessMatrix(string matrixType) const {
    Matrix C = TranspPolygonal :: giveStiffnessMatrix(matrixType);
    double cond = 0;
    for ( unsigned i = 0; i < inttype->giveNumIP(); i++ ) {
        cond += inttype->giveIPWeight(i) * stats [ i ]->giveStiffnessTensor("elastic")(0, 0);
    }
    cond /= volume;

    return V1 * cond + ( ( V2.transpose() * C ) * V2 );
}

//////////////////////////////////////////////////////////
void TranspVirtPolygonal :: computeDampingMatrix() {
    double cap = 0;
    for ( unsigned i = 0; i < inttype->giveNumIP(); i++ ) {
        cap += inttype->giveIPWeight(i) * stats [ i ]->giveDampingTensor()(0, 0);
    }
    cap /= volume;

    /*
     * cout << "------------" << endl;
     * M.print();
     * cout << "------------" << endl;
     * ( W1 * cap + matrix_multiply(matrix_multiply(W2.transpose(), M), W2) ).print();
     * cout << "------------" << endl;
     * ( W1 * cap + matrix_multiply(matrix_multiply(W2.transpose(), M), W2)*volume ).print();
     * cout << "------------" << endl;
     * ( (W1*cap + matrix_multiply(W2.transpose(), W2)*volume*volume )* cap ).print();
     * exit(1);
     * return M;
     */

    dampC =  ( W1 + ( W2.transpose() * W2 ) * volume ) * cap;
}

//////////////////////////////////////////////////////////
Vector TranspVirtPolygonal :: giveInternalForces(const Vector &DoFs, bool frozen, double timeStep) {
    ( void ) frozen;
    ( void ) timeStep;
    //return Element::giveInternalForces(DoFs, frozen); //incorrect integration
    return giveStiffnessMatrix("elastic") * DoFs;  //using VEM integration, only elastic material!
}


//////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////
// TRANSPORT POLYGONAL ELEMENT CONSTRUCTED BY STATIC CONDENSATION OF ISOPARAMETRIC TRIANGLES

TranspCondensedPolygonal :: TranspCondensedPolygonal(const unsigned dim) : TranspPolygonal(dim) {
    name = "TranspCondensedPolygonal";
    ip_type = "tri";
    delete shafunc;
    delete inttype;
    shafunc = new Linear2DPolygonShapeF();
    inttype = new IntegrPolygon(ip_type);
}

//////////////////////////////////////////////////////////
void TranspCondensedPolygonal :: initIntegration() {
    prepareGeometry();

    Linear2DPolygonShapeF *lds = static_cast< Linear2DPolygonShapeF * >( shafunc );
    inttype->init(nodes, faces, & centroid);
    lds->setFacesCentroidAndIntegration(faces, centroid, inttype);
    shafunc->init(nodes);
}


/*
 * //////////////////////////////////////////////////////////
 * // 3 D
 * //////////////////////////////////////////////////////////
 * //////////////////////////////////////////////////////////
 * // POLYHEDRAL FACE
 * PolyhedralFace :: PolyhedralFace(const unsigned dim) {
 *  ndim = 3;
 *  if (dim!=3){
 *      cerr << "Polyhedral Face error: this element can be used only in 3D models" << endl;
 *      exit(1);
 *  }
 *  name = "PolyhedralFace";
 * }
 *
 * //////////////////////////////////////////////////////////
 * void PolyhedralFace :: readFromLine(istringstream &iss, NodeContainer *fullnodes, MaterialContainer *fullmatrs){
 *  (void) fullmatrs;
 *
 *  unsigned num, num2;
 *
 *  iss >> num;
 *  nodes.resize(num);
 *  for ( unsigned i = 0; i < num; i++ ) {
 *      iss >> num2;
 *      nodes [ i ] = fullnodes->giveNode(num2);
 *  }
 * }
 *
 * //////////////////////////////////////////////////////////
 * void PolyhedralFace :: init(){
 * }
 */

/*
 * //////////////////////////////////////////////////////////
 * //////////////////////////////////////////////////////////
 * // TRANSPORT POLYHEDRAL ELEMENT
 * TranspPolyhedral :: TranspPolyhedral(const unsigned dim)  : TranspPolygonal(2) {
 *  ndim = dim;
 *  name = "TranspPolyhedral";
 *  ip_type = "quad";
 * }
 *
 * //////////////////////////////////////////////////////////
 * void TranspPolyhedral :: readFromLine(istringstream &iss, NodeContainer *fullnodes, MaterialContainer *fullmatrs) {
 *  unsigned num, num2;
 *  unordered_set< unsigned >nn;
 *
 *  iss >> nfaces;
 *  faces.resize(nfaces);
 *  for ( unsigned i = 0; i < nfaces; i++ ) {
 *      iss >> num;
 *      faces [ i ].resize(num);
 *      for ( unsigned j = 0; j < num; j++ ) {
 *          iss >> num2;
 *          nn.insert(num2);
 *          faces [ i ] [ j ] = num2;
 *      }
 *  }
 *  iss >> num;
 *  mat = fullmatrs->giveMaterial(num);
 *
 *  numOfNodes = nn.size();
 *  vector< unsigned >nv;
 *  nv.resize(numOfNodes);
 *  nodes.resize(numOfNodes);
 *  faceConnectivity.resize(numOfNodes);
 *
 *  unsigned i = 0;
 *  for ( auto &n: nn ) {
 *      nv [ i ] = n;
 *      nodes [ i ] = fullnodes->giveNode(n);
 *      i++;
 *  }
 *  for ( unsigned f = 0; f < nfaces; f++ ) {
 *      for ( unsigned n = 0; n < faces [ f ].size(); n++ ) {
 *          i = 0;
 *          while ( nv [ i ] != faces [ f ] [ n ] && i < nv.size() ) {
 *              i++;
 *          }
 *          if ( i == nv.size() ) {
 *              cerr << "Inconzistence in Polyhedral input file" << endl;
 *              exit(1);
 *          }
 *          faces [ f ] [ n ] = i;
 *          faceConnectivity [ i ].push_back(f);
 *      }
 *  }
 * }
 *
 *
 * //////////////////////////////////////////////////////////
 * void TranspPolyhedral :: WachspressShapeF(const Point x, MyMatrix phi) const {
 *  //compute hf
 *  MyVector h(nfaces);
 *  for ( unsigned i = 0; i < nfaces; i++ ) {
 *      h [ i ] = dot(nodes [ faces [ i ] [ 0 ] ]->givePoint() - x, normals [ i ]);
 *  }
 *
 *  //compute wf
 *  double sumW = 0;
 *  phi *= 0;
 *  for ( unsigned i = 0; i < numOfNodes; i++ ) {
 *      phi[ i ] = 1;
 *      for ( unsigned j = 0; j < 3; j++ ) {
 *          w [0] [ i ] *= h [ faceConnectivity [ i ] [ j ] ];
 *      }
 *      phi[ i ] = determinants [ i ] / w [ i ];
 *      sumW += phi[ i ];
 *  }
 *  return phi / sumW;
 * }
 *
 * //////////////////////////////////////////////////////////
 * double TranspPolyhedral :: WachspressShapeFGrad(const Point x, MyMatrix phiGrad) const {
 *  MyVector phi(numOfNodes);
 *  WachspressShapeF(x, phi);
 *  MyMatrix R(numOfNodes, 3);
 *
 *  MyVector h(nfaces);
 *  for ( unsigned i = 0; i < nfaces; i++ ) {
 *      h [ i ] = dot(nodes [ faces [ i ] [ 0 ] ]->givePoint() - x, normals [ i ]);
 *  }
 *
 *  unsigned n;
 *  MyVector phiR(3);
 *  phiR [ 0 ] = 0;
 *  phiR [ 1 ] = 0;
 *  phiR [ 2 ] = 0;
 *  for ( unsigned i = 0; i < numOfNodes; i++ ) {
 *      R [ i ] [ 0 ] = 0;
 *      R [ i ] [ 1 ] = 0;
 *      R [ i ] [ 2 ] = 0;
 *      for ( unsigned j = 0; j < 3; j++ ) {
 *          n = faceConnectivity [ i ] [ 0 ];
 *          R [ i ] [ 0 ] += normals [ n ].getX() / h [ n ];
 *          R [ i ] [ 1 ] += normals [ n ].getY() / h [ n ];
 *          R [ i ] [ 2 ] += normals [ n ].getZ() / h [ n ];
 *      }
 *      phiR [ 0 ] += R [ i ] [ 0 ] * phi[ i ];
 *      phiR [ 1 ] += R [ i ] [ 1 ] * phi[ i ];
 *      phiR [ 2 ] += R [ i ] [ 2 ] * phi[ i ];
 *  }
 *
 *  for ( unsigned i = 0; i < numOfNodes; i++ ) {
 *      for ( unsigned j = 0; j < 3; j++ ) {
 *          phiGrad [ j ] [ i ] = phi[ i ] * ( R [ i ] [ j ] - phiR [ j ] );
 *      }
 *  }
 *
 *  return 1;
 * }
 *
 * //////////////////////////////////////////////////////////
 * void TranspPolyhedral :: findIntegrationPoints() {
 *  if ( ip_type.compare("quad") == 0 ) {
 *      //based on quadrilateral isoparametric elements
 *      ip_locs.resize(4 * numOfNodes);
 *      ip_weights.resize(4 * numOfNodes);
 *      stats.resize(4 * numOfNodes);
 *      Point a = ( nodes [ faces [ numOfNodes - 1 ] [ 0 ] ]->givePoint() + nodes [ faces [ numOfNodes - 1 ] [ 1 ] ]->givePoint() ) / 2;
 *      Point b, c, d, derxyr, derxys;
 *      double detJ;
 *      d = centroid;
 *      double r = 1. / sqrt(3.);
 *
 *      for ( unsigned i = 0; i < numOfNodes; i++ ) {
 *          c = ( nodes [ faces [ i ] [ 0 ] ]->givePoint() + nodes [ faces [ i ] [ 1 ] ]->givePoint() ) / 2;
 *          b = nodes [ i ]->givePoint();
 *          for ( int k = -1; k < 2; k = k + 2 ) {
 *              for ( int l = -1; l < 2; l = l + 2 ) {
 *                  ip_locs [ 4 * i + ( k + 1 ) + ( l + 1 ) / 2 ] = ( a * ( ( 1 + k * r ) * ( 1 + l * r ) ) + b * ( ( 1 - k * r ) * ( 1 + l * r ) ) + c * ( ( 1 - k * r ) * ( 1 - l * r ) ) + d * ( ( 1 + k * r ) * ( 1 - l * r ) ) ) / 4.;
 *                  derxyr = ( a * ( 1 + l * r ) - b * ( 1 + l * r ) - c * ( 1 - l * r ) + d * ( 1 - l * r ) ) / 4.;
 *                  derxys = ( a * ( 1 + k * r ) + b * ( 1 - k * r ) - c * ( 1 - k * r ) - d * ( 1 + k * r ) ) / 4.;
 *                  detJ = derxyr.x * derxys.y - derxyr.y * derxys.x;
 *                  ip_weights [ 4 * i + ( k + 1 ) + ( l + 1 ) / 2 ] = detJ;
 *                  stats [ 4 * i + ( k + 1 ) + ( l + 1 ) / 2 ] = mat->giveNewMaterialStatus(this);
 *              }
 *          }
 *          a = c;
 *      }
 *  } else if ( ip_type.compare("tri") == 0 ) {
 *      //based on triangular isoparametric elements
 *      ip_locs.resize(3 * numOfNodes);
 *      ip_weights.resize(3 * numOfNodes);
 *      stats.resize(3 * numOfNodes);
 *      Point a, b;
 *
 *      for ( unsigned i = 0; i < numOfNodes; i++ ) {
 *          a = nodes [ faces [ i ] [ 0 ] ]->givePoint();
 *          b = nodes [ faces [ i ] [ 1 ] ]->givePoint();
 *          double tarea = triArea2D(a, b, centroid);
 *          ip_locs [ 3 * i  ] = b * ( 1. - 1. / 6. - 1. / 6. ) + a / 6. + centroid / 6.;
 *          ip_locs [ 3 * i + 1 ] = b * ( 1. - 2. / 3. - 1. / 6. ) + a * 2. / 3. + centroid / 6.;
 *          ip_locs [ 3 * i + 2 ] = b * ( 1. - 1. / 6. - 2. / 3. ) + a / 6. + centroid * 2. / 3.;
 *          for ( unsigned t = 0; t < 3; t++ ) {
 *              ip_weights [ 3 * i + t ] = tarea / 3.;
 *              stats [ 3 * i + t ] = mat->giveNewMaterialStatus(this);
 *          }
 *      }
 *     } else {
 *      cerr << "Error in " << name << ": ip_type '" << ip_type << "' not implemented" << endl;
 *     }
 * }
 *
 * //////////////////////////////////////////////////////////
 * void TranspPolyhedral :: init() {
 *  Element :: init(); //calling base class method;
 *
 *  //check that nodes are TrsNodes
 *  for ( const auto n: nodes ) {
 *      TrsNode *p = dynamic_cast< TrsNode * >( n );
 *      if ( !p ) {
 *          cerr << "Error in " << name << ": nodes must be inherited from TrsNode, " << n->giveName() << " provided" << endl;
 *          exit(1);
 *      }
 *  }
 *  //check that material is TrsprtMaterial
 *  TrsprtMaterial *p = dynamic_cast< TrsprtMaterial * >( mat );
 *  if ( !p ) {
 *      cerr << "Error in " << name << ": material must be inherited from TrsprtMaterial, " << mat->giveName() << " provided" << endl;
 *      exit(1);
 *  }
 *
 *  //centroid estimation
 *  Point estcentroid = Point(0, 0, 0);
 *  for ( auto &n: nodes ) {
 *      estcentroid += n->givePoint();
 *  }
 *  estcentroid /= numOfNodes;
 *
 *  //centers, areas and volumes of faces
 *  surfaces.resize(nfaces);
 *  faceCenters.resize(nfaces);
 *  volumes.resize(nfaces);
 *
 *  Point *last;
 *  Point *current;
 *  Point estcenter;
 *  unsigned i = 0;
 *  double fa, fv;
 *  volume = 0.;
 *  centroid = Point(0, 0, 0);
 *  for ( auto &f: faces ) {
 *      //estimate face center
 *      estcenter = Point(0, 0, 0);
 *      for ( auto &n: f ) {
 *          estcenter += nodes [ n ]->givePoint();
 *      }
 *      estcenter /= f.size();
 *
 *
 *      last = nodes [ f [ f.size() - 1 ] ]->givePointPointer();
 *      surfaces [ i ] = 0;
 *      volumes [ i ] = 0;
 *      faceCenters [ i ] = Point(0, 0, 0);
 *      for ( auto &n: f ) {
 *          current = nodes [ n ]->givePointPointer(); \
 *          fa = triArea3D(last, current, & estcenter);
 *          fv = -tetVolume3D(last, current, & estcenter, & estcentroid); //normal oriented outside
 *          surfaces [ i ] += fa;
 *          volumes [ i ] += fv;
 *          faceCenters [ i ] = ( * last + * current ) * fa;
 *          centroid = ( * last + * current + estcenter ) * fv;
 *          last = current;
 *      }
 *      faceCenters [ i ] /= 3. * surfaces [ i ];
 *      faceCenters [ i ] += estcenter / 3.;
 *      volume += volumes [ i ];
 *      i++;
 *  }
 *  centroid /= 4. * volume;
 *  centroid += estcentroid / 4;
 *
 *  i = 0;
 *  //correct volumes of each face
 *  for ( auto &f: faces ) {
 *      last = nodes [ f [ f.size() - 1 ] ]->givePointPointer();
 *      volumes [ i ] = 0;
 *      for ( auto &n: f ) {
 *          current = nodes [ n ]->givePointPointer(); \
 *          volumes [ i ] = tetVolume3D(last, current, & faceCenters [ i ], & centroid);
 *          last = current;
 *      }
 *      i++;
 *  }
 *
 *  //normals
 *  i = 0;
 *  for ( auto &f: faces ) {
 *      Point a = nodes [ f [ 1 ] ]->givePoint() - faceCenters [ i ];
 *      Point b = nodes [ f [ 0 ] ]->givePoint() - faceCenters [ i ];
 *      normals [ i ].setX(a.getY() * b.getZ() - a.getZ() * b.getY() );
 *      normals [ i ].setX(a.getZ() * b.getX() - a.getX() * b.getZ() );
 *      normals [ i ].setX(a.getX() * b.getY() - a.getY() * b.getX() );
 *      normals [ i ] /= normals [ i ].norm();
 *      i++;
 *  }
 *
 *  //check face connectivity and compute determinants from normals
 *  i = 0;
 *  determinants.resize(numOfNodes);
 *  Point a, b, c;
 *  for ( auto &fc: faceConnectivity ) {
 *      if ( fc.size() != 3 ) {
 *          cerr << "Error in TranspPolyhedral: each node must have 3 connected faces, " << fc.size() << " found" << endl;
 *          exit(EXIT_FAILURE);
 *      }
 *      a = nodes [ fc [ 0 ] ]->givePoint();
 *      b = nodes [ fc [ 1 ] ]->givePoint();
 *      c = nodes [ fc [ 2 ] ]->givePoint();
 *      determinants [ i ] = a.getX() * ( b.getY() * c.getZ() - b.getZ() * c.getY() ) - b.getX() * ( a.getY() * c.getZ() - a.getZ() * c.getY() ) + c.getX() * ( a.getY() * b.getZ() - a.getZ() * b.getY() );
 *      if ( determinants [ i ] < 0 ) {
 *          determinants [ i ] *= -1;
 *          unsigned p = fc [ 1 ];
 *          fc [ 1 ] = fc [ 2 ];
 *          fc [ 2 ] = p;
 *      }
 *      i++;
 *  }
 *
 *  findIntegrationPoints();
 * }
 */

/*
 * //////////////////////////////////////////////////////////
 * //////////////////////////////////////////////////////////
 * // TRANSPORT VIRTUAL POLYHEDRAL ELEMENT
 *
 * TranspVirtPolyhedral :: TranspVirtPolyhedral(const unsigned dim) : TranspPolyhedral(dim) {
 *  name = "TranspVirtPolyhedral";
 * }
 *
 * //////////////////////////////////////////////////////////
 * void TranspVirtPolyhedral :: init() {
 *  TranspPolyhedral :: init(); //calling base class method;
 *
 *  MyMatrix R(numOfNodes, ndim);
 *  unsigned j = numOfNodes - 1;
 *  for ( unsigned i = 0; i < numOfNodes; i++ ) {
 *      R [ i ] [ 0 ] = ( normals [ i ].x * surfaces [ i ] + normals [ j ].x * surfaces [ j ] ) / 2.;
 *      R [ i ] [ 1 ] = ( normals [ i ].y * surfaces [ i ] + normals [ j ].y * surfaces [ j ] ) / 2.;
 *      j = i;
 *  }
 *
 *  MyMatrix N(numOfNodes, ndim);
 *  Point x;
 *  for ( unsigned i = 0; i < numOfNodes; i++ ) {
 *      x = nodes [ i ]->givePoint();
 *      N [ i ] [ 0 ] = x.x;
 *      N [ i ] [ 1 ] = x.y;
 *  }
 *
 *  MyMatrix I(numOfNodes, numOfNodes);
 *  for ( unsigned i = 0; i < numOfNodes; i++ ) {
 *      I [ i ] [ i ] = 1.;
 *  }
 *
 *  MyMatrix P0(numOfNodes, numOfNodes);
 *  for ( unsigned i = 0; i < numOfNodes; i++ ) {
 *      for ( unsigned j = 0; j < numOfNodes; j++ ) {
 *          P0 [ i ] [ j ] = 1. / numOfNodes;
 *      }
 *  }
 *
 *  MyMatrix H = matrix_multiply(N, R.transpose() ) / volume;
 *  MyMatrix P =  H +  matrix_multiply(P0, I - H);
 *
 *  V1 = matrix_multiply(R, R.transpose() ) / volume;
 *  V2 = I - P;
 * }
 *
 * //////////////////////////////////////////////////////////
 * MyMatrix TranspVirtPolyhedral :: giveConductivityMatrix(string matrixType) const {
 *  MyMatrix C = TranspPolyhedral :: giveConductivityMatrix(matrixType);
 *  TrsprtMaterial *tmat = static_cast< TrsprtMaterial * >( mat );
 *  return V1 * tmat->giveConductivity() + matrix_multiply(matrix_multiply(V2.transpose(), C), V2);
 * }
 *
 *
 * //////////////////////////////////////////////////////////
 * //////////////////////////////////////////////////////////
 * // TRANSPORT POLYHEDRAL ELEMENT CONSTRUCTED BY STATIC CONDENSATION OF ISOPARAMETRIC TRIANGLES
 *
 * TranspCondensedPolyhedral :: TranspCondensedPolyhedral(const unsigned dim) : TranspPolyhedral(dim) {
 *  name = "TranspCondensedPolyhedral";
 *  ip_type = "tri";
 * }
 *
 * //////////////////////////////////////////////////////////
 * MyVector TranspCondensedPolyhedral :: fullTriShapeF(Point x) const {
 *  unsigned face = findFaceNumber(x);
 *  MyVector phi(0., numOfNodes + 1); //include the centroid
 *  double tarea = triArea2D(centroid, nodes [ faces [ face ] [ 0 ] ]->givePoint(), nodes [ faces [ face ] [ 1 ] ]->givePoint() );
 *  phi [ faces [ face ] [ 1 ] ] = triArea2D(x, centroid, nodes [ faces [ face ] [ 0 ] ]->givePoint() ) / tarea;
 *  phi [ faces [ face ] [ 0 ] ] = triArea2D(x, nodes [ faces [ face ] [ 1 ] ]->givePoint(), centroid) / tarea;
 *  phi [ numOfNodes ] = triArea2D(x, nodes [ faces [ face ] [ 0 ] ]->givePoint(), nodes [ faces [ face ] [ 1 ] ]->givePoint() ) / tarea;
 *  return phi;
 * }
 *
 * //////////////////////////////////////////////////////////
 * MyMatrix TranspCondensedPolyhedral :: fullTriShapeFGrad(Point x) const {
 *  unsigned face = findFaceNumber(x);
 *  MyMatrix phiGrad(ndim, numOfNodes + 1); //include the centroid
 *  double tarea = triArea2D(centroid, nodes [ faces [ face ] [ 0 ] ]->givePoint(), nodes [ faces [ face ] [ 1 ] ]->givePoint() );
 *  phiGrad [ 0 ] [ faces [ face ] [ 1 ] ] = 0.5 * ( centroid.getY() - nodes [ faces [ face ] [ 0 ] ]->givePoint().getY() ) / tarea;
 *  phiGrad [ 1 ] [ faces [ face ] [ 1 ] ] = 0.5 * ( nodes [ faces [ face ] [ 0 ] ]->givePoint().getX() - centroid.getX() ) / tarea;
 *  phiGrad [ 0 ] [ faces [ face ] [ 0 ] ] = 0.5 * ( nodes [ faces [ face ] [ 1 ] ]->givePoint().getY() - centroid.getY() ) / tarea;
 *  phiGrad [ 1 ] [ faces [ face ] [ 0 ] ] = 0.5 * ( centroid.getX() - nodes [ faces [ face ] [ 1 ] ]->givePoint().getX() ) / tarea;
 *  phiGrad [ 0 ] [ numOfNodes ]         = 0.5 * ( nodes [ faces [ face ] [ 0 ] ]->givePoint().getY() - nodes [ faces [ face ] [ 1 ] ]->givePoint().getY() ) / tarea;
 *  phiGrad [ 1 ] [ numOfNodes ]         = 0.5 * ( nodes [ faces [ face ] [ 1 ] ]->givePoint().getX() - nodes [ faces [ face ] [ 0 ] ]->givePoint().getX() ) / tarea;
 *  return phiGrad;
 * }
 *
 * //////////////////////////////////////////////////////////
 * unsigned TranspCondensedPolyhedral :: findFaceNumber(Point x) const {
 *  double alpha = atan2(x.getY() - centroid.getY(), x.getX() - centroid.getX() );
 *  unsigned face;
 *  for ( face = 0; face < numOfNodes; face++ ) {
 *      if ( alpha < angles [ faces [ face ] [ 1 ] ] &&  alpha > angles [ faces [ face ] [ 0 ] ] ) {
 *          return face;
 *      }
 *  }
 *  return nodeMaxAngle;
 * }
 *
 * //////////////////////////////////////////////////////////
 * MyVector TranspCondensedPolyhedral :: condTriShapeF(Point x) const {
 *  MyVector full = fullTriShapeF(x);
 *  MyVector reduced(numOfNodes);
 *  for ( unsigned i = 0; i < numOfNodes; i++ ) {
 *      reduced [ i ] = full [ i ] + full [ numOfNodes ] * red2full [ i ];
 *  }
 *  return reduced;
 * }
 *
 * //////////////////////////////////////////////////////////
 * MyMatrix TranspCondensedPolyhedral :: condTriShapeFGrad(Point x) const {
 *  MyMatrix full = fullTriShapeFGrad(x);
 *  MyMatrix reduced(ndim, numOfNodes);
 *  for ( unsigned d = 0; d < ndim; d++ ) {
 *      for ( unsigned i = 0; i < numOfNodes; i++ ) {
 *          reduced [ d ] [ i ] = full [ d ] [ i ] + full [ d ] [ numOfNodes ] * red2full [ i ];
 *      }
 *  }
 *  return reduced;
 * }
 *
 * //////////////////////////////////////////////////////////
 * void TranspCondensedPolyhedral :: init() {
 *  TranspPolyhedral :: init(); //calling base class method;
 *  angles.resize(numOfNodes);
 *
 *  nodeMaxAngle = 0;
 *  double maxAngle = -2;
 *  for ( unsigned i = 0; i < numOfNodes; i++ ) {
 *      angles [ i ] = atan2(nodes [ i ]->givePoint().getY() - centroid.getY(), nodes [ i ]->givePoint().getX() - centroid.getX() );
 *      if (maxAngle<angles [ i ]){
 *          maxAngle = angles [ i ];
 *          nodeMaxAngle = i;
 *      }
 *  }
 *
 *  //build transformation matrix allowing to calculate inner degree of freedom
 *  MyMatrix FullK(numOfNodes + 1, numOfNodes + 1);
 *  MyMatrix phiGrad;
 *  for ( size_t i = 0; i < ip_weights.size(); i++ ) {
 *      phiGrad = fullTriShapeFGrad(ip_locs [ i ]);
 *      FullK += matrix_multiply(phiGrad.transpose(), phiGrad) * ip_weights [ i ];
 *  }
 *
 *  red2full.resize(numOfNodes);
 *  for ( unsigned i = 0; i < numOfNodes; i++ ) {
 *      red2full [ i ] = -FullK [ i ] [ numOfNodes ] / FullK [ numOfNodes ] [ numOfNodes ];
 *  }
 * }
 */
