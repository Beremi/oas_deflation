#include "element_polyhedral.h"

//////////////////////////////////////////////////////////
double triArea(const Point a, const Point b, const Point c) { //points in counter clockwise direction
    return 0.5 * ( a.getX() * ( b.getY() - c.getY() ) + b.getX() * ( c.getY() - a.getY() ) + c.getX() * ( a.getY() - b.getY() ) );
}

//////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////
// TRANSPORT POLYHEDRAL ELEMENT
TranspPolyhedral :: TranspPolyhedral(const unsigned dim) {
    ndim = dim;
    name = "TranspPolyhedral";
    ip_type = "quad";
}

//////////////////////////////////////////////////////////
void TranspPolyhedral :: readFromLine(istringstream &iss, NodeContainer *fullnodes, MaterialContainer *fullmatrs) {
    unsigned num, num2;

    iss >> nnodes;
    nodes.resize(nnodes);
    for ( unsigned i = 0; i < nnodes; i++ ) {
        iss >> num2;
        nodes [ i ] = fullnodes->giveNode(num2);
    }
    iss >> num;
    mat = fullmatrs->giveMaterial(num);
}

//////////////////////////////////////////////////////////
void TranspPolyhedral :: sort2D() {
    //estimate centroid
    for ( const auto n: nodes ) {
        centroid += n->givePoint();
    }
    centroid /= nnodes;

    //compute angles
    vector< pair< double, unsigned > >angles;
    angles.resize(nnodes);
    for ( unsigned i = 0; i < nnodes; i++ ) {
        angles [ i ].second = i;
        angles [ i ].first = atan2(nodes [ i ]->givePoint().getY() - centroid.getY(), nodes [ i ]->givePoint().getX() - centroid.getX() );
    }
    ;

    //sort to have counterclockwise direction
    sort(angles.begin(), angles.end() );
    vector< Node * >newnodes;
    newnodes.resize(nnodes);
    //do NOT change node ordering, the fisrt face should run through +-pi axis
    for ( unsigned i = 0; i < nnodes - 1; i++ ) {
        newnodes [ i + 1 ] = nodes [ angles [ i ].second ];
    }
    newnodes [ 0 ] = nodes [ angles [ nnodes - 1 ].second ];
    nodes = newnodes;
}

//////////////////////////////////////////////////////////
Vector TranspPolyhedral :: WachspressShapeF(Point x) const {
    Vector h(nnodes);
    Vector w(nnodes);
    Point oldNormal = normals [ nnodes - 1 ];
    for ( unsigned i = 0; i < nnodes; i++ ) {
        h [ i ] = abs(dot(nodes [ faces [ i ] [ 0 ] ]->givePoint() - x, normals [ i ]) );
        w [ i ] = abs(oldNormal.x * normals [ i ].y - oldNormal.y * normals [ i ].x);
        oldNormal = normals [ i ];
    }
    double oldh = h [ nnodes - 1 ];
    double sumW = 0;
    for ( unsigned i = 0; i < nnodes; i++ ) {
        w [ i ] /= oldh * h [ i ];
        sumW += w [ i ];
        oldh = h [ i ];
    }
    return w / sumW;
}

//////////////////////////////////////////////////////////
Matrix TranspPolyhedral :: WachspressShapeFGrad(Point x) const {
    Vector phi = WachspressShapeF(x);
    Matrix R(nnodes, 2);
    Matrix phiGrad(2, nnodes);

    Vector h(nnodes);
    for ( unsigned i = 0; i < nnodes; i++ ) {
        h [ i ] = abs(dot(nodes [ faces [ i ] [ 0 ] ]->givePoint() - x, normals [ i ]) );
    }
    unsigned oldi = nnodes - 1;
    Vector phiR(2);
    phiR [ 0 ] = 0;
    phiR [ 1 ] = 0;
    for ( unsigned i = 0; i < nnodes; i++ ) {
        R [ i ] [ 0 ] = normals [ oldi ].getX() / h [ oldi ] + normals [ i ].getX() / h [ i ];
        R [ i ] [ 1 ] = normals [ oldi ].getY() / h [ oldi ] + normals [ i ].getY() / h [ i ];
        phiR [ 0 ] += R [ i ] [ 0 ] * phi [ i ];
        phiR [ 1 ] += R [ i ] [ 1 ] * phi [ i ];
        oldi = i;
    }
    for ( unsigned i = 0; i < nnodes; i++ ) {
        for ( unsigned j = 0; j < 2; j++ ) {
            phiGrad [ j ] [ i ] = phi [ i ] * ( R [ i ] [ j ] - phiR [ j ] );
        }
    }
    return phiGrad;
}

//////////////////////////////////////////////////////////
void TranspPolyhedral :: findIntegrationPoints() {
    if ( ip_type.compare("quad") == 0 ) {
        //based on quadrilateral isoparametric elements
        ip_locs.resize(4 * nnodes);
        ip_weights.resize(4 * nnodes);
        stats.resize(4 * nnodes);
        Point a = ( nodes [ faces [ nnodes - 1 ] [ 0 ] ]->givePoint() + nodes [ faces [ nnodes - 1 ] [ 1 ] ]->givePoint() ) / 2;
        Point b, c, d, derxyr, derxys;
        double detJ;
        d = centroid;
        double r = 1. / sqrt(3.);

        for ( unsigned i = 0; i < nnodes; i++ ) {
            c = ( nodes [ faces [ i ] [ 0 ] ]->givePoint() + nodes [ faces [ i ] [ 1 ] ]->givePoint() ) / 2;
            b = nodes [ i ]->givePoint();
            for ( int k = -1; k < 2; k = k + 2 ) {
                for ( int l = -1; l < 2; l = l + 2 ) {
                    ip_locs [ 4 * i + ( k + 1 ) + ( l + 1 ) / 2 ] = ( a * ( ( 1 + k * r ) * ( 1 + l * r ) ) + b * ( ( 1 - k * r ) * ( 1 + l * r ) ) + c * ( ( 1 - k * r ) * ( 1 - l * r ) ) + d * ( ( 1 + k * r ) * ( 1 - l * r ) ) ) / 4.;
                    derxyr = ( a * ( 1 + l * r ) - b * ( 1 + l * r ) - c * ( 1 - l * r ) + d * ( 1 - l * r ) ) / 4.;
                    derxys = ( a * ( 1 + k * r ) + b * ( 1 - k * r ) - c * ( 1 - k * r ) - d * ( 1 + k * r ) ) / 4.;
                    detJ = derxyr.x * derxys.y - derxyr.y * derxys.x;
                    ip_weights [ 4 * i + ( k + 1 ) + ( l + 1 ) / 2 ] = detJ;
                    stats [ 4 * i + ( k + 1 ) + ( l + 1 ) / 2 ] = mat->giveNewMaterialStatus(this);
                }
            }
            a = c;
        }
    } else if ( ip_type.compare("tri") == 0 )       {
        //based on triangular isoparametric elements
        ip_locs.resize(3 * nnodes);
        ip_weights.resize(3 * nnodes);
        stats.resize(3 * nnodes);
        Point a, b;

        for ( unsigned i = 0; i < nnodes; i++ ) {
            a = nodes [ faces [ i ] [ 0 ] ]->givePoint();
            b = nodes [ faces [ i ] [ 1 ] ]->givePoint();
            double tarea = triArea(a, b, centroid);
            ip_locs [ 3 * i  ] = b * ( 1. - 1. / 6. - 1. / 6. ) + a / 6. + centroid / 6.;
            ip_locs [ 3 * i + 1 ] = b * ( 1. - 2. / 3. - 1. / 6. ) + a * 2. / 3. + centroid / 6.;
            ip_locs [ 3 * i + 2 ] = b * ( 1. - 1. / 6. - 2. / 3. ) + a / 6. + centroid * 2. / 3.;
            for ( unsigned t = 0; t < 3; t++ ) {
                ip_weights [ 3 * i + t ] = tarea / 3.;
                stats [ 3 * i + t ] = mat->giveNewMaterialStatus(this);
            }
        }
    } else  {
        cerr << "Error in " << name << ": ip_type '" << ip_type << "' not implemented" << endl;
    }
}

//////////////////////////////////////////////////////////
void TranspPolyhedral :: init() {
    //reorder nodes before calling base calls initialization
    sort2D();
    Element :: init(); //calling base class method;

    //check that nodes are TrsNodes
    for ( const auto n: nodes ) {
        TrsNode *p = dynamic_cast< TrsNode * >( n );
        if ( !p ) {
            cerr << "Error in " << name << ": nodes must be inherited from TrsNode, " << n->giveName() << " provided" << endl;
            exit(1);
        }
    }
    //check that material is DisMechMat
    TrsprtMaterial *p = dynamic_cast< TrsprtMaterial * >( mat );
    if ( !p ) {
        cerr << "Error in " << name << ": material must be inherited from TrsprtMaterial, " << mat->giveName() << " provided" << endl;
        exit(1);
    }

    volume = 0;
    double triVolume;
    if ( ndim == 2 ) {
        if ( nodes.size() < 3 ) {
            cerr << "Error: more than 2 nodes must be involved, " << nodes.size() << " provided" << endl;
            exit(1);
        }

        //area,centroid, faces, ...
        faces.resize(nnodes);
        normals.resize(nnodes);
        areas.resize(nnodes);
        Point cpoint = centroid;
        centroid = Point(0., 0., 0.);
        Point diff;
        for ( unsigned i = 0; i < nnodes; i++ ) {
            faces [ i ].resize(2);
            faces [ i ] [ 0 ] = i;
            faces [ i ] [ 1 ] = ( i == nnodes - 1 ) ? 0 : i + 1;
            diff = nodes [ faces [ i ] [ 1 ] ]->givePoint() - nodes [ faces [ i ] [ 0 ] ]->givePoint();
            areas [ i ] = diff.norm();
            normals [ i ] = Point(diff.y / areas [ i ], -diff.x / areas [ i ], 0);
            triVolume = triArea(nodes [ faces [ i ] [ 0 ] ]->givePoint(), nodes [ faces [ i ] [ 1 ] ]->givePoint(), cpoint);
            centroid += ( nodes [ faces [ i ] [ 0 ] ]->givePoint() + nodes [ faces [ i ] [ 1 ] ]->givePoint() + cpoint ) * triVolume;
            volume += triVolume;
        }
        centroid /= volume * 3.;
    } else if ( ndim == 3 )  {
        cerr << name << ": 3rd dimension not implemented yet" << endl;
        exit(1);
    }

    findIntegrationPoints();
}

//////////////////////////////////////////////////////////
Matrix TranspPolyhedral :: giveConductivityMatrix(string matrixType) const {
    ( void ) matrixType;
    Matrix C(nnodes, nnodes);
    Matrix phiGrad;
    TrsprtMaterialStatus *tstats;
    for ( size_t i = 0; i < ip_weights.size(); i++ ) {
        phiGrad = shapeFGrad(ip_locs [ i ]);
        tstats = static_cast< TrsprtMaterialStatus * >( stats [ 0 ] );
        C += matrix_multiply(phiGrad.transpose(), phiGrad) * ( ip_weights [ i ] * tstats->giveConductivity() );
    }
    return C;
}

//////////////////////////////////////////////////////////
Matrix TranspPolyhedral :: giveCapacityMatrix() const {
    Matrix S(2, 2);
    /*
     * TrsprtMaterialStatus *tstats = static_cast< TrsprtMaterialStatus * >( stats [ 0 ] );
     * double s = area * tstats->giveCapacity() * length / 6.;
     * S [ 0 ] [ 0 ] = S [ 1 ] [ 1 ] = 2 * s;
     * S [ 1 ] [ 0 ] = S [ 0 ] [ 1 ] = s;
     */
    return S;
}

//////////////////////////////////////////////////////////
Vector TranspPolyhedral :: giveInternalForces(const Vector &DoFs) const {
    return giveConductivityMatrix("secant") * DoFs;
};

//////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////
// TRANSPORT VIRTUAL POLYHEDRAL ELEMENT

TranspVirtPolyhedral :: TranspVirtPolyhedral(const unsigned dim) : TranspPolyhedral(dim) {
    name = "TranspVirtPolyhedral";
}

//////////////////////////////////////////////////////////
void TranspVirtPolyhedral :: init() {
    TranspPolyhedral :: init(); //calling base class method;

    Matrix R(nnodes, ndim);
    unsigned j = nnodes - 1;
    for ( unsigned i = 0; i < nnodes; i++ ) {
        R [ i ] [ 0 ] = ( normals [ i ].x * areas [ i ] + normals [ j ].x * areas [ j ] ) / 2.;
        R [ i ] [ 1 ] = ( normals [ i ].y * areas [ i ] + normals [ j ].y * areas [ j ] ) / 2.;
        j = i;
    }

    Matrix N(nnodes, ndim);
    Point x;
    for ( unsigned i = 0; i < nnodes; i++ ) {
        x = nodes [ i ]->givePoint();
        N [ i ] [ 0 ] = x.x;
        N [ i ] [ 1 ] = x.y;
    }

    Matrix I(nnodes, nnodes);
    for ( unsigned i = 0; i < nnodes; i++ ) {
        I [ i ] [ i ] = 1.;
    }

    Matrix P0(nnodes, nnodes);
    for ( unsigned i = 0; i < nnodes; i++ ) {
        for ( unsigned j = 0; j < nnodes; j++ ) {
            P0 [ i ] [ j ] = 1. / nnodes;
        }
    }

    Matrix H = matrix_multiply(N, R.transpose() ) / volume;
    Matrix P =  H +  matrix_multiply(P0, I - H);

    V1 = matrix_multiply(R, R.transpose() ) / volume;
    V2 = I - P;
}

//////////////////////////////////////////////////////////
Matrix TranspVirtPolyhedral :: giveConductivityMatrix(string matrixType) const {
    Matrix C = TranspPolyhedral :: giveConductivityMatrix(matrixType);
    TrsprtMaterialStatus *tstats = static_cast< TrsprtMaterialStatus * >( stats [ 0 ] );
    return V1 * tstats->giveConductivity() + matrix_multiply(matrix_multiply(V2.transpose(), C), V2);
}


//////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////
// TRANSPORT POLYHEDRAL ELEMENT CONSTRUCTED BY STATIC CONDENSATION OF ISOPARAMETRIC TRIANGLES

TranspCondensedPolyhedral :: TranspCondensedPolyhedral(const unsigned dim) : TranspPolyhedral(dim) {
    name = "TranspCondensedPolyhedral";
    ip_type = "tri";
}

//////////////////////////////////////////////////////////
Vector TranspCondensedPolyhedral :: fullTriShapeF(Point x) const {
    double alpha = atan2(x.getY() - centroid.getY(), x.getX() - centroid.getX() );
    unsigned face;
    for ( face = 0; face < nnodes; face++ ) {
        if ( alpha < angles [ faces [ face ] [ 1 ] ] ) {
            break;
        }
    }
    if ( face == nnodes ) {
        face = 0;
    }
    Vector phi(0., nnodes + 1); //include the centroid
    double tarea = triArea(centroid, nodes [ faces [ face ] [ 0 ] ]->givePoint(), nodes [ faces [ face ] [ 1 ] ]->givePoint() );
    phi [ faces [ face ] [ 1 ] ] = triArea(x, centroid, nodes [ faces [ face ] [ 0 ] ]->givePoint() ) / tarea;
    phi [ faces [ face ] [ 0 ] ] = triArea(x, nodes [ faces [ face ] [ 1 ] ]->givePoint(), centroid) / tarea;
    phi [ nnodes ] = triArea(x, nodes [ faces [ face ] [ 0 ] ]->givePoint(), nodes [ faces [ face ] [ 1 ] ]->givePoint() ) / tarea;
    return phi;
}

//////////////////////////////////////////////////////////
Matrix TranspCondensedPolyhedral :: fullTriShapeFGrad(Point x) const {
    double alpha = atan2(x.getY() - centroid.getY(), x.getX() - centroid.getX() );
    unsigned face;
    for ( face = 0; face < nnodes; face++ ) {
        if ( alpha < angles [ faces [ face ] [ 1 ] ] ) {
            break;
        }
    }
    if ( face == nnodes ) {
        face = 0;
    }
    Matrix phiGrad(ndim, nnodes + 1); //include the centroid
    double tarea = triArea(centroid, nodes [ faces [ face ] [ 0 ] ]->givePoint(), nodes [ faces [ face ] [ 1 ] ]->givePoint() );
    phiGrad [ 0 ] [ faces [ face ] [ 1 ] ] = 0.5 * ( centroid.getY() - nodes [ faces [ face ] [ 0 ] ]->givePoint().getY() ) / tarea;
    phiGrad [ 1 ] [ faces [ face ] [ 1 ] ] = 0.5 * ( nodes [ faces [ face ] [ 0 ] ]->givePoint().getX() - centroid.getX() ) / tarea;
    phiGrad [ 0 ] [ faces [ face ] [ 0 ] ] = 0.5 * ( nodes [ faces [ face ] [ 1 ] ]->givePoint().getY() - centroid.getY() ) / tarea;
    phiGrad [ 1 ] [ faces [ face ] [ 0 ] ] = 0.5 * ( centroid.getX() - nodes [ faces [ face ] [ 1 ] ]->givePoint().getX() ) / tarea;
    phiGrad [ 0 ] [ nnodes ]         = 0.5 * ( nodes [ faces [ face ] [ 0 ] ]->givePoint().getY() - nodes [ faces [ face ] [ 1 ] ]->givePoint().getY() ) / tarea;
    phiGrad [ 1 ] [ nnodes ]         = 0.5 * ( nodes [ faces [ face ] [ 1 ] ]->givePoint().getX() - nodes [ faces [ face ] [ 0 ] ]->givePoint().getX() ) / tarea;
    return phiGrad;
}

//////////////////////////////////////////////////////////
Vector TranspCondensedPolyhedral :: condTriShapeF(Point x) const {
    Vector full = fullTriShapeF(x);
    Vector reduced(nnodes);
    for ( unsigned i = 0; i < nnodes; i++ ) {
        reduced [ i ] = full [ i ] + full [ nnodes ] * red2full [ i ];
    }
    return reduced;
}

//////////////////////////////////////////////////////////
Matrix TranspCondensedPolyhedral :: condTriShapeFGrad(Point x) const {
    Matrix full = fullTriShapeFGrad(x);
    Matrix reduced(ndim, nnodes);
    for ( unsigned d = 0; d < ndim; d++ ) {
        for ( unsigned i = 0; i < nnodes; i++ ) {
            reduced [ d ] [ i ] = full [ d ] [ i ] + full [ d ] [ nnodes ] * red2full [ i ];
        }
    }
    return reduced;
}

//////////////////////////////////////////////////////////
void TranspCondensedPolyhedral :: init() {
    TranspPolyhedral :: init(); //calling base class method;
    angles.resize(nnodes);

    for ( unsigned i = 0; i < nnodes; i++ ) {
        angles [ i ] = atan2(nodes [ i ]->givePoint().getY() - centroid.getY(), nodes [ i ]->givePoint().getX() - centroid.getX() );
    }

    //build transformation matrix allowing to calculate inner degree of freedom
    Matrix FullK(nnodes + 1, nnodes + 1);
    Matrix phiGrad;
    for ( size_t i = 0; i < ip_weights.size(); i++ ) {
        phiGrad = fullTriShapeFGrad(ip_locs [ i ]);
        FullK += matrix_multiply(phiGrad.transpose(), phiGrad) * ip_weights [ i ];
    }

    red2full.resize(nnodes);
    for ( unsigned i = 0; i < nnodes; i++ ) {
        red2full [ i ] = -FullK [ i ] [ nnodes ] / FullK [ nnodes ] [ nnodes ];
    }
}
