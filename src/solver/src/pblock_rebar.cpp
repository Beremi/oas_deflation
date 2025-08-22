#include "pblock_rebar.h"
#include "model.h"
#include "element_beam.h"
#include "material_beam.h"
#include "cross_section.h"

using namespace std;


bool sortPairsA(const pair<double,Element*>& p1, const pair<double,Element*>& p2)
{
   return p1.first > p2.first;
}

bool sortPairsB(const pair<double,Node*>& p1, const pair<double,Node*>& p2)
{
   return p1.first > p2.first;
}

//////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////
// REBAR
//////////////////////////////////////////////////////////
void Rebar :: readFromLine(istringstream &iss, unsigned d) {
    dim = d;
    string param;
    unsigned num;
    while (  iss >> param ) {
        if ( param.compare("nodes") == 0 ) {
            iss >> num;
            nodes.resize(num);
            for(unsigned i=0; i<nodes.size(); i++){            
                iss >> nodes[i];
            }
        } else if ( param.compare("diameter") == 0 ) {
            iss >> radius;
            radius = radius/2;
        } else if ( param.compare("material") == 0 ) {
            iss >> material_id;
        } else {
            std :: cerr << "parameter '" << param << "' not implemented yet" << '\n';
            exit(EXIT_FAILURE);
        }
    }
}

//////////////////////////////////////////////////////////
void Rebar :: apply(Model *model) {

    BeamMaterial * bm = dynamic_cast<BeamMaterial *> (model->giveMaterials()->giveMaterial(material_id));
    if (!bm) {
        cerr << "Error in rebar preprocessing block: material must be derived from Beam Material" << endl;
    }

    findIntersectionsWithElements(model);
}


//////////////////////////////////////////////////////////
void Rebar :: findIntersectionsWithElements(Model *model) {
    RigidBodyContact *rbc;
    Point normal, dirvec, intersec;
    Point a, b, minp, maxp;
    Point * s, *r;
    Point auxA, auxB, centroid;
    double t, length, d=0;
    bool bintersect;
    double vol1, vol2, vol3;
    vector< Node * >verts;
    vector< pair<double, Element* > > intersections; 
    

    BeamMaterial * bm = dynamic_cast<BeamMaterial *> (model->giveMaterials()->giveMaterial(material_id));
    CircularCrossSection * cs = new CircularCrossSection(radius,bm->givePoissonsRatio());  
    cs -> init();
    model->giveCrossSections()->addCrossSection(cs); 

    //find elements for reinforcement nodes
    Vector bbox;
    Point loc;
    Point coord = Point(0,0,0);
    bool found, inside;
    for(unsigned i=0; i<nodes.size(); i++){
        loc = model->giveNodes()->giveNode(nodes[i])->givePoint();
        
        //inster particle instead of AuxNode if needed
        if(dynamic_cast<AuxNode*>(model->giveNodes()->giveNode(nodes[i]))){
            Particle* newp = new Particle(dim, 0, loc);
            model->giveNodes()->addNode(newp); 
            nodes[i] = newp->giveID();
        }        

        found = false;
        for ( vector< Element * > :: iterator ee = model->giveElements()->begin(); ee!= model->giveElements()->end() && !found; ++ee) {    
            if (not (*ee)->doesMechanics()) continue; 
            bbox = (*ee)->giveBoundingBox();
            inside=true;
            //test bounding box
            for(unsigned v=0; v<dim; v++){
                if(loc[v]<bbox[2*v] || loc[v]>bbox[2*v+1]) inside=false;                
            }  
            //true test             
            if(inside){
                inside=(*ee)->isPointInside(&coord, &loc);
            }
            if(inside){
                found = true;
                rbc = dynamic_cast< RigidBodyContact * >( *ee );
                if (rbc){
                    model->giveConstraints()->addRigidArmConstraint(dim,model->giveNodes()->giveNode(nodes[i]),rbc->giveNode(unsigned(coord[0]+0.2)),false);
                }else{                    
                    model->giveConstraints()->addHangingNodeConstraint(model->giveNodes()->giveNode(nodes[i]), *ee);
                }
            }
        }
        if(!found){
            cout << "Warning: rebar node " << nodes[i] << " at coordinates " << loc[0] << " " << loc[1] << " " << loc[2] << " is not connected to any mechanical element" << endl;
        }
    }

    //find centers of intersection with elements
    for(unsigned k=1; k<nodes.size(); k++){
        intersections.clear();
        a = model->giveNodes()->giveNode(nodes[k-1])->givePoint();
        b = model->giveNodes()->giveNode(nodes[k])->givePoint();
        dirvec = (b-a);
        length = dirvec.norm();
        dirvec /= length;
        minp = Point(min(a[0],b[0]),min(a[1],b[1]),min(a[2],b[2]));        
        maxp = Point(max(a[0],b[0]),max(a[1],b[1]),max(a[2],b[2]));        
        for ( vector< Element * > :: iterator ee = model->giveElements()->begin(); ee!= model->giveElements()->end(); ++ee) {
            if (!( *ee )->doesMechanics()) continue;
            rbc = dynamic_cast< RigidBodyContact * >( *ee );
            if (rbc){
                bbox = rbc->giveFacetBoundingBox();
                inside=true;
                //test bounding box
                for(unsigned v=0; v<dim; v++){
                    if(maxp[v]<bbox[2*v]) inside=false;
                    if(minp[v]>bbox[2*v+1]) inside=false;            
                }   
                if(!inside) continue;
                //compute intersection with facet plane
                normal = rbc->giveNormal();
                centroid = rbc->giveCentroid();
                d = -( centroid [ 0 ] * normal [ 0 ] + centroid [ 1 ] * normal [ 1 ] + centroid [ 2 ] * normal [ 2 ] );
                t = -( normal [ 0 ] * a [ 0 ] + normal [ 1 ] * a[ 1 ] + normal [ 2 ] * a [ 2 ] + d ) / ( normal [ 0 ] * dirvec [ 0 ] + normal [ 1 ] * dirvec [ 1 ] + normal [ 2 ] * dirvec [ 2 ] );
                if ( t <= 0 || t >= length ) {
                    continue;
                }
                intersec = a + t * dirvec;
                //check that it is inside the facet
                //according to https://stackoverflow.com/questions/42740765/intersection-between-line-and-triangle-in-3d
                bintersect = false;
                auxA = intersec + dirvec * ( 5 * sqrt(rbc->giveArea() ) );
                auxB = intersec - dirvec * ( 5 * sqrt(rbc->giveArea() ) );
                verts = rbc->giveVertices();
                for ( unsigned i = 0; i < verts.size() && !bintersect; i++ ) {
                    s = verts [ i ]->givePointPointer();
                    if ( i == 0 ) {
                        r = verts [ verts.size() - 1 ]->givePointPointer();
                    } else {
                        r = verts [ i - 1 ]->givePointPointer();
                    }
                    vol1 = tetraVolumeSigned(& auxA, & auxB, & centroid, r);
                    vol2 = tetraVolumeSigned(& auxA, & auxB, r, s);
                    vol3 = tetraVolumeSigned(& auxA, & auxB, s, & centroid);
                    if ( ( vol1 < 0 && vol2 < 0 && vol3 < 0 ) || ( vol1 > 0 && vol2 > 0 && vol3 > 0 ) ) {
                        bintersect = true;
                    }
                }
                if ( !bintersect ) {
                    continue;
                }
                intersections.push_back(pair(t,rbc));
            //elems that are not RigidBodyContact
            } else {                
                Vector tt = ( *ee )->findIntersectionsWithLine(&a, &b);    
                for(unsigned i=0; i<tt.size(); i++){
                    intersections.push_back(pair(tt[i],(*ee)));
                }
            }
        }
    
        //find centers of intersection with particles and elements
        vector<Node*>primaryParticles;
        vector<Element*>primaryElements;
        vector<pair<double,double>>minmaxt;
        sort(intersections.begin(), intersections.end(), sortPairsA);
        unsigned l;
        
        for(unsigned i=0; i<intersections.size(); i++){
            if(dynamic_cast<RigidBodyContact*>(intersections[i].second)){
                for(unsigned j=0; j<2; j++){
                    auto m = find(primaryParticles.begin(), primaryParticles.end(), intersections[i].second->giveNode(j));
                    if (m!=primaryParticles.end()){
                        l = m-primaryParticles.begin();
                        minmaxt[l].first=intersections[i].first;        
                    } else {
                        primaryParticles.push_back(intersections[i].second->giveNode(j));
                        primaryElements.push_back(nullptr);
                        minmaxt.push_back(pair(-1,intersections[i].first));
                    }
                }
            }else{
                auto m = find(primaryElements.begin(), primaryElements.end(), intersections[i].second);
                if (m!=primaryElements.end()){
                    l = m-primaryElements.begin();
                    minmaxt[l].first=intersections[i].first;        
                } else {
                    primaryParticles.push_back(nullptr);
                    primaryElements.push_back(intersections[i].second);
                    minmaxt.push_back(pair(-1,intersections[i].first));
                }                
            }
        }

        vector<pair<double, Node*>> newts;
        double newt;
        //add new nodes and constraints
        for(unsigned i=0; i<primaryParticles.size(); i++){
            if(minmaxt[i].first>0){
                newt = (minmaxt[i].first+minmaxt[i].second)/2.;
                if (newts.size()==0 || abs(newts.back().first-newt)>length/1000.){
                    Point newloc = a + dirvec*newt;
                    Particle* newp = new Particle(dim, 0, newloc);
                    model->giveNodes()->addNode(newp);
                    if(primaryParticles[i]){
                        model->giveConstraints()->addRigidArmConstraint(dim,newp,primaryParticles[i],false);
                    }else{
                        model->giveConstraints()->addHangingNodeConstraint(newp,primaryElements[i]);
                    }
                    newts.push_back(pair(newt,newp));
                }  
            }                       
        }
        //create rebar beams
        sort(newts.begin(), newts.end(), sortPairsB);
        Node* second = model->giveNodes()->giveNode(nodes[k]);
        Node* first;
        for(vector<pair<double,Node*>>::const_iterator tt = newts.begin(); tt!=newts.end(); ++tt){
            first = tt->second;
            TimoshenkoBeam3D * tb = new TimoshenkoBeam3D(first , second, bm, cs, Point(100,100,100));
            model->giveElements()->addElement(tb);            
            second = first;
        }        
        first = model->giveNodes()->giveNode(nodes[k-1]);  
        TimoshenkoBeam3D * tb = new TimoshenkoBeam3D(first , second, bm, cs, Point(100,100,100));
        model->giveElements()->addElement(tb);  

    }  
}
