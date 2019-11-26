#include "element_polyhedral.h"

//////////////////////////////////////////////////////////
double triArea(const Point a, const Point b, const Point c){ //points in counter clockwise direction
        return 0.5* ( a.getX()*(b.getY()-c.getY()) + b.getX()*(c.getY()-a.getY()) + c.getX()*(a.getY()-b.getY()) );
}

//////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////
// TRANSPORT POLYHEDRAL ELEMENT
TranspPolyhedral :: TranspPolyhedral(const unsigned dim) {
    ndim = dim;
    ip_locs.resize(1);
    stats.resize(1);
    name = "TranspPolyhedral";
}

//////////////////////////////////////////////////////////
void TranspPolyhedral :: readFromLine(istringstream &iss, NodeContainer *fullnodes, MaterialContainer *fullmatrs) {
    unsigned num, num2;

    iss >> nnodes;
    nodes.resize(nnodes);
    for(unsigned i=0; i<nnodes; i++) {
        iss >> num2;
        nodes [ i ] = fullnodes->giveNode(num2);
    }
    iss >> num;
    mat = fullmatrs->giveMaterial(num);
}

//////////////////////////////////////////////////////////
void TranspPolyhedral :: sort2D() {

    //estimate centroid
    for ( const auto n: nodes) centroid += n->givePoint();
    centroid /= nnodes;

    //compute angles
    vector<pair <double,unsigned> > angles;
    angles.resize(nnodes); 
    for ( unsigned i=0; i<nnodes; i++){
        angles[i].second = i;
        angles[i].first = atan2(nodes[i]->givePoint().getY() - centroid.getY(),nodes[i]->givePoint().getX() - centroid.getX());        
    };

    //sort to have counterclockwise direction
    sort(angles.begin(), angles.end());  
    vector< Node * > newnodes;
    newnodes.resize(nnodes);
    for ( unsigned i=0; i<nnodes; i++) newnodes[i] = nodes[angles[i].second];
    nodes = newnodes;    
}

//////////////////////////////////////////////////////////
Vector TranspPolyhedral :: WachspressShapeF(Point x) const{
    Vector h = Vector(nnodes);
    Vector w = Vector(nnodes);    
    Point oldNormal = normals[nnodes-1]; 
    for(unsigned i=0; i<nnodes; i++){
        h[i] = abs(dot(nodes[faces[i][0]]->givePoint()-x,normals[i]));
        w[i] = abs(oldNormal.x*normals[i].y-oldNormal.y*normals[i].x);                
        oldNormal = normals[i];
    }
    double oldh = h[nnodes-1];
    double sumW = 0;
    for(unsigned i=0; i<nnodes; i++){
        w[i] /= oldh*h[i];
        sumW += w[i];
        oldh = h[i];
    }    
    return w/sumW;   
}

//////////////////////////////////////////////////////////
Matrix TranspPolyhedral :: WachspressShapeFGrad(Point x) const{
    Vector phi = WachspressShapeF(x);
    Matrix R;

    /*
        R = np.divide(np.vstack((self.normals[-1,:],self.normals[:-1,:])),np.tile(np.hstack((h[-1],h[:-1])),(2,1)).T) + np.divide(self.normals,np.tile(h,(2,1)).T)
        phiR = np.multiply(np.tile(phi,(2,1)).T,R)
        return np.multiply(np.tile(phi,(2,1)).T,R-np.tile(np.sum(phiR,axis=0),(self.nn,1)))


    vector < double > h, w;
    Point oldNormal = normals[nnodes-1]; 
    for(unsigned i=0; i<nnodes; i++){
        h[i] = abs(dot(nodes[faces[i][0]]->givePoint()-x,normals[i]));
        w[i] = abs(oldNormal.x*normals[i].y-oldNormal.y*normals[i].x);                
        oldNormal = normals[i];
    }
    double oldh = h[nnodes-1];
    double sumW = 0;
    for(unsigned i=0; i<nnodes; i++){
        w[i] /= oldh*h[i];
        sumW += w[i];
        oldh = h[i];
    }
    for(unsigned i=0; i<nnodes; i++) w[i] /= sumW;
    return w;   
    */
    return R;
}

//////////////////////////////////////////////////////////
void TranspPolyhedral :: init() {
    Element :: init(); //calling base class method;

    //check that nodes are TrsNodes
    for ( const auto n: nodes) {
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
    centroid = Point(0.,0.,0.);
    if ( ndim == 2 ) {
        if ( nodes.size() <3 ) {
            cerr << "Error: more than 2 nodes must be involved, " << nodes.size() << " provided" << endl;
            exit(1);
        }

        //reorder nodes
        sort2D();

        //area,centroid, faces, ...
        faces.resize(nnodes);
        normals.resize(nnodes);
        areas.resize(nnodes);
        Point cpoint = centroid;
        Point diff;
        for (unsigned i = 0; i< nnodes; i++){
            faces[i].resize(2);
            faces[i][0] = (i == 0) ? -1 : i-1;
            faces[i][1] = i; 
            diff = nodes[i]->givePoint() - nodes[faces[i][0]]->givePoint(); 
            areas[i] = diff.norm();
            normals[i] = Point(diff.y/areas[i], -diff.x/areas[i],0);  
            triVolume = triArea(nodes[faces[i][0]]->givePoint(),nodes[i]->givePoint(), cpoint);
            centroid += (nodes[faces[i][0]]->givePoint() + nodes[i]->givePoint() + cpoint) * triVolume;
            volume += triVolume;
        }
        centroid /= volume*3.;

    }else if ( ndim == 3 ) {  
        cerr << name << ": 3rd dimension not implemented yet" << endl; 
        exit(1);
    }
  

    /*
    stats [ 0 ] = mat->giveNewMaterialStatus(this);
    normal = nodes [ 1 ]->givePoint() - nodes [ 0 ]->givePoint();
    length = normal.norm();
    normal = normal / length;

    if ( abs(normal * t) > 1e-4) {
        cout << v0.x << " " <<  v0.y <<  " X " << v1.x << " " <<  v1.y << endl;
        cout << nodes [ 0 ]->givePoint().x << " " <<  nodes [ 0 ]->givePoint().y <<  " X " << nodes [ 1 ]->givePoint().x << " " <<  nodes [ 1 ]->givePoint().y << endl;
        cerr << "Error: normal and contact vector are not parallel, error " << normal * t << " normal v." << normal.x << " " << normal.y << " contact v. " << t.x << " " << t.y << endl;
        exit(1);
    }
    */
}

//////////////////////////////////////////////////////////
Matrix TranspPolyhedral :: giveConductivityMatrix(string matrixType) const {    
    Matrix C(2, 2);
    /*
    TrsprtMaterialStatus *tstats = static_cast< TrsprtMaterialStatus * >( stats [ 0 ] );
    double c = area * tstats->giveConductivity() / length;
    C [ 0 ] [ 0 ] = C [ 1 ] [ 1 ] = c;
    C [ 1 ] [ 0 ] = C [ 0 ] [ 1 ] = -c;
    */
    return C;
}

//////////////////////////////////////////////////////////
Matrix TranspPolyhedral :: giveCapacityMatrix() const {    
    Matrix S(2, 2);
    /*
    TrsprtMaterialStatus *tstats = static_cast< TrsprtMaterialStatus * >( stats [ 0 ] );
    double s = area * tstats->giveCapacity() * length / 6.;
    S [ 0 ] [ 0 ] = S [ 1 ] [ 1 ] = 2 * s;
    S [ 1 ] [ 0 ] = S [ 0 ] [ 1 ] = s;
    */
    return S;
}

//////////////////////////////////////////////////////////
Vector TranspPolyhedral :: giveInternalForces(const Vector &DoFs) const {
    return giveConductivityMatrix("elastic") * DoFs;
};
