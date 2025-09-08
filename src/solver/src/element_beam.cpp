#include "element_beam.h"
#include "boundary_condition.h"
#include "material_beam.h"

using namespace std;

//////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////
// TIMOSHENKO BEAM ELEMENT
/////////////////////////////////////////////////////////
TimoshenkoBeam3D::TimoshenkoBeam3D(unsigned dim):Element(dim){
    ndim = dim;
    name = "TimoshenkoBeam3D";
    numOfNodes = 2;
    nodes.resize(numOfNodes);
    vtk_cell_type = 3;
    shafunc = new NullShapeF(1);
    inttype = new IntegrLine();
    zdir = Point(0,0,1);
}

/////////////////////////////////////////////////////////
TimoshenkoBeam3D::TimoshenkoBeam3D(Node* a, Node* b, Material * m, CrossSection *cs, Point zrefpoint):Element(3){
    ndim = 3;
    name = "TimoshenkoBeam3D";
    numOfNodes = 2;
    nodes.resize(numOfNodes);
    nodes[0] = a;
    nodes[1] = b;
    vtk_cell_type = 3;
    shafunc = new NullShapeF(1);
    inttype = new IntegrLine();
    zdir = zrefpoint;
    CS = cs;
    mat = m;
}

/////////////////////////////////////////////////////////
void TimoshenkoBeam3D::setIntegrationPointsAndWeights(){
        IntegrLine *it = static_cast< IntegrLine * >( inttype );
        it->setNumIP(3);
        stats.resize(inttype->giveNumIP());
        for ( unsigned k = 0; k < inttype->giveNumIP(); k++ ) {
            stats [ k ] = mat->giveNewMaterialStatus(this, k);
            BeamMaterialStatus* tm = dynamic_cast<BeamMaterialStatus*>(stats [ k ]);
            if (tm) tm -> setCrossSection(CS);
            it->setIPWeight(k, it->giveIPWeight(k)*length/2);
            it->setIPLocation(k, (it->giveIPLocation(k)+Point(1,0,0))/2);
        }
}

/////////////////////////////////////////////////////////
void TimoshenkoBeam3D::computeDampingMatrix(){
    dampC = Matrix :: Zero(DoFids.size(), DoFids.size());
}

/////////////////////////////////////////////////////////
void TimoshenkoBeam3D::readFromLine(std :: istringstream &iss, NodeContainer *fullnodes, MaterialContainer *fullmatrs, CrossSectionContainer *csc){
    
    bool csb, zb;
    csb = zb = false;

    unsigned num;
    iss >> num;
    nodes [ 0 ] = fullnodes->giveNode(num);
    iss >> num;
    nodes [ 1 ] = fullnodes->giveNode(num);
    iss >> num;
    mat = fullmatrs->giveMaterial(num);
    
    string code;
    
    string param;
    while ( iss >> code ) {
        if ( code.compare("crosssection") == 0 ) {
            iss >> num;
            CS = csc->giveCrossSection(num);
            csb = true;
        } else if  ( code.compare("zrefpoint") == 0 ) {
            double x,y,z;
            iss >> x >> y >> z;
            zdir = Point(x,y,z);            
            zb = true;
        }
    }
    
    if (!csb) {
        cerr << "Error in TimoshenkoBeam3D: cross section was not defined, use parameter 'crosssection'" << endl;
        exit(1);
    }
    if (!zb) {
        cerr << "Error in TimoshenkoBeam3D: local reference system was not defined, use parameter 'zrefpoint'" << endl;
        exit(1);
    }
}

/////////////////////////////////////////////////////////
void TimoshenkoBeam3D::init(){
    normal = nodes [ 1 ]->givePoint() - nodes [ 0 ]->givePoint();
    length = normal.norm();
    normal /= length;

    Point t2 = zdir - normal*(zdir.dot(normal));
    t2.normalize();
    Point t1 = t2.cross(normal);
    R = Matrix :: Zero( 12, 12 );
    for (unsigned i=0; i<4; i++){
        R(3*i+0,3*i+0) = normal.x();
        R(3*i+0,3*i+1) = normal.y();
        R(3*i+0,3*i+2) = normal.z();
        R(3*i+1,3*i+0) = t1.x();
        R(3*i+1,3*i+1) = t1.y();
        R(3*i+1,3*i+2) = t1.z();
        R(3*i+2,3*i+0) = t2.x();
        R(3*i+2,3*i+1) = t2.y();
        R(3*i+2,3*i+2) = t2.z();
    }

    Element :: init(); //calling base class method;
    BeamMaterial* tm = dynamic_cast<BeamMaterial*>(mat);

    if(!tm){
        cerr << "Error in TimoshenkoBeam3D: material must be instance of BeamMaterial" << endl;
        exit(1);
    }
    checkNodeType();
}

/////////////////////////////////////////////////////////
void TimoshenkoBeam3D::checkNodeType() const{
    //check that nodes are particles
    for ( unsigned i = 0; i < 2; i++ ) {
        Particle *p = dynamic_cast< Particle * >( nodes [ i ] );
        if ( !p ) {
            cerr << "Error in " << name << ": nodes must be inherited from Particle, " << nodes [ i ]->giveName() << " provided" << endl;
            exit(1);
        }
    }

    //check that material is TensMechMat
    TensMechMaterial *p = dynamic_cast< TensMechMaterial * >( mat );
    if ( !p ) {
        cerr << "Error in " << name << ": material must be inherited from TensMechMaterial, " << mat->giveName() << " provided" << endl;
        exit(1);
    }
}
    

/////////////////////////////////////////////////////////
Matrix TimoshenkoBeam3D::giveBMatrix(const Point *x) const {
    //EN234: Three-dimentional Timoshenko beam element undergoing axial, torsional and bending deformations, Wenqiang Fang, 2015
    double ksi  = x->x();
    BeamMaterial* tm = static_cast<BeamMaterial*>(mat);
    double FiZ = 12*tm->giveElasticModulus()*CS->giveIz()/(CS->giveKappaY()*tm->giveShearModulus()*CS->giveArea()*pow(length,2));
    double FiY = 12*tm->giveElasticModulus()*CS->giveIy()/(CS->giveKappaZ()*tm->giveShearModulus()*CS->giveArea()*pow(length,2));
    double ksi2 = pow(ksi,2);
    Matrix B = Matrix :: Zero( 6, 12 );
    //x strain.
    B(0, 0) = -1./length;
    B(0, 6) = 1./length;
    //y strain = diff v - theta z
    B(1, 1) = 1./(1.+FiZ)/length*(-6.*ksi+6.*ksi2-FiZ);
    B(1, 5) = 1./(1.+FiZ)*(1. - 4.*ksi + 3.*ksi2 + FiZ/2.*(1.-2.*ksi));
    B(1, 7) = 1./(1.+FiZ)/length*(6.*ksi-6.*ksi2+FiZ);
    B(1, 11) = 1./(1.+FiZ)*(- 2.*ksi + 3.*ksi2 + FiZ/2.*(-1.+2.*ksi));
    // rotations z
    B(1, 1) -= 6./(1.+FiZ)/length*(-ksi+ksi2);
    B(1, 5) -= 1./(1.+FiZ)*(1.-4.*ksi+3.*ksi2+FiZ*(1.-ksi));
    B(1, 7) -= -6./(1.+FiZ)/length*(-ksi+ksi2);
    B(1, 11) -= 1./(1.+FiZ)*(-2.*ksi+3.*ksi2+FiZ*ksi);
    //z strain = diff w + theta y
    B(2, 2) = 1./(1.+FiY)/length*(-6.*ksi+6.*ksi2-FiY);
    B(2, 4) = -1./(1.+FiY)*(1 - 4*ksi + 3*ksi2 + FiY/2.*(1-2.*ksi));
    B(2, 8) = 1./(1.+FiY)/length*(6.*ksi-6.*ksi2+FiY);
    B(2, 10) = -1./(1.+FiY)*(- 2*ksi + 3*ksi2 + FiY/2.*(-1+2.*ksi));
    // rotations y
    B(2, 2) -= 6./(1.+FiY)/length*(-ksi+ksi2);
    B(2, 4) -= -1./(1.+FiY)*(1.-4*ksi+3*ksi2+FiY*(1.-ksi));
    B(2, 8) -= -6./(1.+FiY)/length*(-ksi+ksi2);
    B(2, 10) -= -1./(1.+FiY)*(-2.*ksi+3.*ksi2+FiY*ksi);
    // x rot diff
    B(3, 3) = -1./length;
    B(3, 9) = 1./length;
    // y curvature
    B(4, 2) = -6./(1.+FiY)/pow(length,2)*(-1.+2.*ksi);
    B(4, 4) = 1./(1.+FiY)/length*(-4.+6.*ksi-FiY);
    B(4, 8) = 6./(1.+FiY)/pow(length,2)*(-1.+2.*ksi);
    B(4, 10) = 1./(1.+FiY)/length*(-2.+6.*ksi+FiY);
    // z curvature
    B(5, 1) = 6./(1.+FiZ)/pow(length,2)*(-1.+2.*ksi);
    B(5, 5) = 1./(1.+FiZ)/length*(-4.+6.*ksi-FiZ);
    B(5, 7) = -6./(1.+FiZ)/pow(length,2)*(-1.+2.*ksi);
    B(5, 11) = 1./(1.+FiZ)/length*(-2.+6.*ksi+FiZ);
    return B*R;
}

/////////////////////////////////////////////////////////
Matrix TimoshenkoBeam3D::giveHMatrix(const Point *x) const{
    //EN234: Three-dimentional Timoshenko beam element undergoing axial, torsional and bending deformations, Wenqiang Fang, 2015
    double ksi  = x->x();
    TensMechMaterial* tm = static_cast<TensMechMaterial*>(mat);
    double FiZ = 12*tm->giveElasticModulus()*CS->giveIz()/(CS->giveKappaY()*tm->giveShearModulus()*CS->giveArea()*pow(length,2));
    double FiY = 12*tm->giveElasticModulus()*CS->giveIy()/(CS->giveKappaZ()*tm->giveShearModulus()*CS->giveArea()*pow(length,2));
    double ksi2 = pow(ksi,2);
    double ksi3 = pow(ksi,3);

    Matrix H = Matrix :: Zero( 6, 12 );
    //x displ.
    H(0, 0) = 1. - ksi;
    H(0, 6) = ksi;
    //y displ.
    H(1, 1) = 1./(1.+FiZ)*(1-3*ksi2+2*ksi3+FiZ*(1-ksi));
    H(1, 5) = length/(1.+FiZ)*(ksi - 2*ksi2 + ksi3 + FiZ/2.*(ksi-ksi2));
    H(1, 7) = 1./(1.+FiZ)*(3*ksi2-2*ksi3+FiZ*ksi);
    H(1, 11) = length/(1.+FiZ)*(- ksi2 + ksi3 + FiZ/2.*(-ksi+ksi2));
    //z displ.
    H(2, 2) = 1./(1.+FiY)*(1-3*ksi2+2*ksi3+FiY*(1-ksi));
    H(2, 4) = -length/(1.+FiY)*(ksi - 2*ksi2 + ksi3 + FiY/2.*(ksi-ksi2));
    H(2, 8) = 1./(1.+FiY)*(3*ksi2-2*ksi3+FiY*ksi);
    H(2, 10) = -length/(1.+FiY)*(- ksi2 + ksi3 + FiY/2.*(-ksi+ksi2));
    // rotations x
    H(3, 3) = 1. - ksi;
    H(3, 9) = ksi;
    // rotations y
    H(4, 2) = 6./(1.+FiY)/length*(-ksi+ksi2);
    H(4, 4) = -1./(1.+FiY)*(1.-4*ksi+3*ksi2+FiY*(1-ksi));
    H(4, 8) = -6./(1.+FiY)/length*(-ksi+ksi2);
    H(4, 10) = -1./(1.+FiY)*(-2*ksi+3*ksi2+FiY*ksi);
    // rotations z
    H(5, 1) = 6./(1.+FiZ)/length*(-ksi+ksi2);
    H(5, 5) = 1./(1.+FiZ)*(1.-4*ksi+3*ksi2+FiZ*(1-ksi));
    H(5, 7) = -6./(1.+FiZ)/length*(-ksi+ksi2);
    H(5, 11) = 1./(1.+FiZ)*(-2*ksi+3*ksi2+FiZ*ksi);

    return H*R;
}

/////////////////////////////////////////////////////////
void TimoshenkoBeam3D::giveValues(std :: string code, Vector &result) const{
    if ( code.compare("length") == 0 ) {
        result.resize(1);
        result [ 0 ] = length;
    } else if ( code.compare("internal_forces") == 0 ) {
        Vector IF = Vector::Zero(6);
        for ( unsigned k = 0; k < inttype->giveNumIP(); k++ ) {
            IF += stats [ k ]->giveTempStress()*inttype->giveIPWeight(k);            
        }        
        result.resize(IF.size());
        for(unsigned i=0; i<IF.size(); i++){
            result[i] = IF[i]/length;
        }
    } else if ( code.compare("diameter") == 0 ) {
        CircularCrossSection* circ = dynamic_cast<CircularCrossSection*>(CS);
        if (circ){
          result.resize(1);
          result[0] = circ->giveDiameter();
        }
    }
    
}

/////////////////////////////////////////////////////////
Vector TimoshenkoBeam3D::integrateLoad(BodyLoad *vl, double time) const{
    (void) vl; (void) time;
    return Vector :: Zero( DoFids.size() );
}

/////////////////////////////////////////////////////////
Vector TimoshenkoBeam3D::integrateInternalSources(){
    return Vector :: Zero( DoFids.size() );
} 
