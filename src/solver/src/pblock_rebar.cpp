#include "pblock_rebar.h"
#include "model.h"
#include "element_beam.h"
#include "material_beam.h"
#include "cross_section.h"
#include "element_ldpm.h"

using namespace std;


bool sortPairsA(const pair< double, Node * > &p1, const pair< double, Node * > &p2)
{
    return p1.first > p2.first;
}

bool sortPairsB(const pair< double, int > &p1, const pair< double, int > &p2)
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
            for (unsigned i = 0; i < nodes.size(); i++) {
                iss >> nodes [ i ];
            }
        } else if ( param.compare("diameter") == 0 ) {
            iss >> radius;
            radius = radius / 2;
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
    BeamMaterial *bm = dynamic_cast< BeamMaterial * >( model->giveMaterials()->giveMaterial(material_id) );
    if ( !bm ) {
        cerr << "Error in rebar preprocessing block: material must be derived from Beam Material" << endl;
    }

    findIntersectionsWithElements(model);
}


//////////////////////////////////////////////////////////
void Rebar :: findIntersectionsWithElements(Model *model) {
    RigidBodyContact *rbc;
    LDPMTetra *tet;
    Point normal, dirvec, intersec;
    Point a, b, minp, maxp;
    Point *s, *r;
    Point auxA, auxB, centroid;
    double t, length, d = 0;
    bool bintersect;
    double vol1, vol2, vol3;
    vector< Node * >verts;
    vector< pair< double, Node * > >intersections;


    BeamMaterial *bm = dynamic_cast< BeamMaterial * >( model->giveMaterials()->giveMaterial(material_id) );
    CircularCrossSection *cs = new CircularCrossSection( radius, bm->givePoissonsRatio() );
    cs->init();
    model->giveCrossSections()->addCrossSection(cs);

    //find elements for reinforcement nodes
    Vector bbox;
    Point loc;
    Point coord = Point(0, 0, 0);
    bool found, inside;
    for (vector< unsigned > :: iterator nn = nodes.begin(); nn != nodes.end(); ++nn) {
        //one node might be repeated
        vector< unsigned > :: iterator same = find(nodes.begin(), nn, * nn);
        if ( same != nn ) {
            continue;
        }

        loc = model->giveNodes()->giveNode(* nn)->givePoint();

        //insert particle instead of AuxNode if needed
        if ( dynamic_cast< AuxNode * >( model->giveNodes()->giveNode(* nn) ) ) {
            Particle *newp = new Particle(dim, 0, loc);
            model->giveNodes()->addNode(newp);
            replace( nodes.begin(), nodes.end(), * nn, newp->giveID() );
        }

        found = false;
        for ( vector< Element * > :: iterator ee = model->giveElements()->begin(); ee != model->giveElements()->end() && !found; ++ee) {
            if ( not ( * ee )->doesMechanics() ) {
                continue;
            }
            bbox = ( * ee )->giveBoundingBox();
            inside = true;
            //test bounding box
            for (unsigned v = 0; v < dim; v++) {
                if ( loc [ v ] < bbox [ 2 * v ] || loc [ v ] > bbox [ 2 * v + 1 ] ) {
                    inside = false;
                }
            }
            //true test
            if ( inside ) {
                inside = ( * ee )->isPointInside(& coord, & loc);
            }
            if ( inside ) {
                found = true;
                rbc = dynamic_cast< RigidBodyContact * >( * ee );
                tet = dynamic_cast< LDPMTetra * >( * ee );
                if ( rbc || tet ) {
                    model->giveConstraints()->addRigidArmConstraint(dim, model->giveNodes()->giveNode(* nn), ( * ee )->giveNode( unsigned( coord [ 0 ] + 0.2 ) ), false);
                } else  {
                    model->giveConstraints()->addHangingNodeConstraint(model->giveNodes()->giveNode(* nn), * ee);
                }
            }
        }
        if ( !found ) {
            cout << "Warning: rebar node " << ( * nn ) << " at coordinates " << loc [ 0 ] << " " << loc [ 1 ] << " " << loc [ 2 ] << " is not connected to any mechanical element" << endl;
        }
    }

    //find centers of intersection with elements
    for (unsigned k = 1; k < nodes.size(); k++) {
        vector< pair< double, int > >newts;
        vector< Element * >primElems;
        intersections.clear();
        a = model->giveNodes()->giveNode(nodes [ k - 1 ])->givePoint();
        b = model->giveNodes()->giveNode(nodes [ k ])->givePoint();
        dirvec = ( b - a );
        length = dirvec.norm();
        dirvec /= length;
        minp = Point( min(a [ 0 ], b [ 0 ]), min(a [ 1 ], b [ 1 ]), min(a [ 2 ], b [ 2 ]) );
        maxp = Point( max(a [ 0 ], b [ 0 ]), max(a [ 1 ], b [ 1 ]), max(a [ 2 ], b [ 2 ]) );
        for ( vector< Element * > :: iterator ee = model->giveElements()->begin(); ee != model->giveElements()->end(); ++ee) {
            if ( !( * ee )->doesMechanics() ) {
                continue;
            }
            rbc = dynamic_cast< RigidBodyContact * >( * ee );
            tet = dynamic_cast< LDPMTetra * >( * ee );
            if ( rbc ) {
                bbox = rbc->giveFacetBoundingBox();
                inside = true;
                //test bounding box
                for (unsigned v = 0; v < dim; v++) {
                    if ( maxp [ v ] < bbox [ 2 * v ] ) {
                        inside = false;
                    }
                    if ( minp [ v ] > bbox [ 2 * v + 1 ] ) {
                        inside = false;
                    }
                }
                if ( !inside ) {
                    continue;
                }
                //compute intersection with facet plane
                normal = rbc->giveNormal();
                centroid = rbc->giveCentroid();
                d = -( centroid [ 0 ] * normal [ 0 ] + centroid [ 1 ] * normal [ 1 ] + centroid [ 2 ] * normal [ 2 ] );
                t = -( normal [ 0 ] * a [ 0 ] + normal [ 1 ] * a [ 1 ] + normal [ 2 ] * a [ 2 ] + d ) / ( normal [ 0 ] * dirvec [ 0 ] + normal [ 1 ] * dirvec [ 1 ] + normal [ 2 ] * dirvec [ 2 ] );
                if ( t <= 0 || t >= length ) {
                    continue;
                }
                intersec = a + t * dirvec;
                //check that it is inside the facet
                //according to https://stackoverflow.com/questions/42740765/intersection-between-line-and-triangle-in-3d
                bintersect = false;
                auxA = intersec + dirvec * ( 5 * sqrt( rbc->giveArea() ) );
                auxB = intersec - dirvec * ( 5 * sqrt( rbc->giveArea() ) );
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
                intersections.push_back( pair( t, rbc->giveNode(0) ) );
                intersections.push_back( pair( t, rbc->giveNode(1) ) );
            } else if ( tet ) {
                Vector tt = ( * ee )->findIntersectionsWithLine(& a, & b);
                if ( tt.size() > 1 ) {
                    sort( tt.begin(), tt.end() );
                    Point B = a + dirvec * tt [ 0 ];
                    Point E = a + dirvec * tt [ tt.size() - 1 ];
                    double distB = ( tet->giveNode(0)->givePoint() - B ).squaredNorm();
                    unsigned indB = 0;
                    double distE = ( tet->giveNode(0)->givePoint() - E ).squaredNorm();
                    unsigned indE = 0;
                    for (unsigned i = 1; i < 4; i++) {
                        if ( ( tet->giveNode(i)->givePoint() - B ).squaredNorm() < distB ) {
                            distB = ( tet->giveNode(i)->givePoint() - B ).squaredNorm();
                            indB = i;
                        }
                        if ( ( tet->giveNode(i)->givePoint() - E ).squaredNorm() < distE ) {
                            distE = ( tet->giveNode(i)->givePoint() - E ).squaredNorm();
                            indE = i;
                        }
                    }
                    intersections.push_back( pair( tt [ 0 ], tet->giveNode(indB) ) );
                    intersections.push_back( pair( tt [ tt.size() - 1 ], tet->giveNode(indE) ) );
                    double q = ( tt [ 0 ] + tt [ tt.size() - 1 ] ) / 2;
                    intersections.push_back( pair( q, tet->giveNode(indB) ) );
                    intersections.push_back( pair( q, tet->giveNode(indE) ) );
                }
                //elems that are not RigidBodyContact
            } else {
                Vector tt = ( * ee )->findIntersectionsWithLine(& a, & b);
                if ( tt.size() > 1 ) {
                    sort( tt.begin(), tt.end() );
                    primElems.push_back(* ee);
                    newts.push_back({ ( tt [ 0 ] + tt [ tt.size() - 1 ] ) / 2., primElems.size() });
                }
            }
        }

        //find centers of intersection with particles
        vector< Node * >primParticles;
        vector< pair< double, double > >minmaxt;
        sort(intersections.begin(), intersections.end(), sortPairsA);
        unsigned l;
        for (unsigned i = 0; i < intersections.size(); i++) {
            if ( intersections [ i ].second == model->giveNodes()->giveNode(nodes [ k - 1 ]) || intersections [ i ].second == model->giveNodes()->giveNode(nodes [ k - 1 ]) ) {
                continue;
            }
            auto m = find(primParticles.begin(), primParticles.end(), intersections [ i ].second);
            if ( m != primParticles.end() ) {
                l = m - primParticles.begin();
                minmaxt [ l ].first = intersections [ i ].first;
            } else {
                primParticles.push_back(intersections [ i ].second);
                minmaxt.push_back({ -1., intersections [ i ].first });
            }
        }
        double newt;
        //add new nodes and constraints
        for (unsigned i = 0; i < primParticles.size(); i++) {
            if ( minmaxt [ i ].first >= 0 ) {
                newt = ( minmaxt [ i ].first + minmaxt [ i ].second ) / 2.;
                if ( newts.size() == 0 || abs(newts.back().first - newt) > length / 10000. ) {
                    newts.push_back({ newt, -( i + 1 ) });
                }
            }
        }
        //add new nodes, constraints and
        Particle *newp;
        sort(newts.begin(), newts.end(), sortPairsB);
        Node *second = model->giveNodes()->giveNode(nodes [ k ]);
        for (vector< pair< double, int > > :: const_iterator tt = newts.begin(); tt != newts.end(); ++tt) {
            newp = new Particle(dim, 0, a + dirvec * tt->first);
            model->giveNodes()->addNode(newp);
            TimoshenkoBeam3D *tb = new TimoshenkoBeam3D( newp, second, bm, cs, Point(100, 100, 100) );
            model->giveElements()->addElement(tb);
            second = newp;
            if ( tt->second > 0 ) {
                model->giveConstraints()->addHangingNodeConstraint(newp, primElems [ tt->second - 1 ]);
                auto nnn = primElems [ tt->second - 1 ]->giveNodes();
            } else  {
                model->giveConstraints()->addRigidArmConstraint(dim, newp, primParticles [ -tt->second - 1 ], false);
            }
        }
        TimoshenkoBeam3D *tb = new TimoshenkoBeam3D( model->giveNodes()->giveNode(nodes [ k - 1 ]), second, bm, cs, Point(100, 100, 100) );
        model->giveElements()->addElement(tb);
    }
}
