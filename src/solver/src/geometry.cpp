// TODO move all functions regarding geometrical operations etc.
// TODO change thsi to class region -> inherited classes circle, polygon etc with same fns isinside etc
#include "geometry.h"
#include <fstream>
#include <memory>

using namespace std;

// Define Infinite (Using INT_MAX caused overflow problems)
#define INF 10000.0
#define TOL 1e-9


Block :: Block(const Point &lB, const Point &rT, unsigned d) : RegularRegion(d) {
    dim = d;
    this->mainPoint = lB;
    this->rightTop = rT;
}

void Block :: readFromLine(istringstream &iss) {
    double x, y, z;
    string param;
    bool ba = false;
    bool bb = false;
    while (  iss >> param ) {
        if ( param.compare("min_point") == 0 ) {
            iss >> x >> y;
            if ( dim == 3 ) {
                iss >> z;
            } else {
                z = 0;
            }
            mainPoint = Point(x, y, z);
            ba = true;
        }
        if ( param.compare("max_point") == 0 ) {
            iss >> x >> y;
            if ( dim == 3 ) {
                iss >> z;
            } else {
                z = 0;
            }
            rightTop = Point(x, y, z);
            bb = true;
        }
    }
    if ( !ba ) {
        cerr << name << " Error: parameter min_point not specified" << endl;
        exit(1);
    }
    if ( !bb ) {
        cerr << name << " Error: parameter max_point not specified" << endl;
        exit(1);
    }
}

bool Block :: isInside(const Point &P) const {
    if ( ( this->mainPoint.x() - P.x() ) < TOL && ( P.x() - this->rightTop.x() ) < TOL ) {
        if ( ( this->mainPoint.y() - P.y() ) < TOL && ( P.y() - this->rightTop.y() ) < TOL ) {
            if ( dim == 3 ) {
                if ( ( this->mainPoint.z() - P.z() ) < TOL && ( P.z() - this->rightTop.z() ) < TOL ) {
                    return true;
                }
            } else {
                return true;
            }
        }
    }
    return false;
}

///////////////////////////////////////////////////////////////////////////////
Circle :: Circle(const Point &c, const double &r) : RegularRegion(2) {
    this->mainPoint = c;
    this->size = r;
}

void Circle :: readFromLine(istringstream &iss) {
    double x, y;
    string param;
    bool bcenter = false;
    bool brad = false;
    while (  iss >> param ) {
        if ( param.compare("center") == 0 ) {
            iss >> x >> y;
            mainPoint = Point(x, y, 0);
            bcenter = true;
        }
        if ( param.compare("radius") == 0 ) {
            iss >> size;
            brad = true;
        }
    }
    if ( !bcenter ) {
        cerr << name << " Error: parameter cente not specified" << endl;
        exit(1);
    }
    if ( !brad ) {
        cerr << name << " Error: parameter radius not specified" << endl;
        exit(1);
    }
    init();
}

bool Circle :: isInside(const Point &P) const {
    // JK in 2D circle == sphere
    // in 3D circle = cylinder with axis along z axis
    // TODO make this more universal to be able to switch cylinder axis to x or y
    Point P_to_compare = P;
    if ( this->along == 'x' ) {
        P_to_compare.x() = 0;
    } else if ( this->along == 'x' ) {
        P_to_compare.y() = 0;
    } else if ( this->along == 'z' ) {
        P_to_compare.z() = 0;
    }
    if ( ( P - this->mainPoint ).norm() < this->size ) {
        return true;
    }
    return false;
}


///////////////////////////////////////////////////////////////////////////////
Sphere :: Sphere(const Point &c, const double &r) {
    dim = 3;

    this->mainPoint = c;
    this->size = r;
}

void Sphere :: readFromLine(istringstream &iss) {
    double x, y, z;
    string param;
    bool bcenter = false;
    bool brad = false;
    while (  iss >> param ) {
        if ( param.compare("center") == 0 ) {
            iss >> x >> y >> z;
            mainPoint = Point(x, y, z);
            bcenter = true;
        }
        if ( param.compare("radius") == 0 ) {
            iss >> size;
            brad = true;
        }
    }
    if ( !bcenter ) {
        cerr << name << " Error: parameter cente not specified" << endl;
        exit(1);
    }
    if ( !brad ) {
        cerr << name << " Error: parameter radius not specified" << endl;
        exit(1);
    }

    init();
}

bool Sphere :: isInside(const Point &P) const {
    if ( ( P - this->mainPoint ).norm() < this->size ) {
        return true;
    }
    return false;
}

///////////////////////////////////////////////////////////////////////////////
Polygon :: Polygon(const std :: vector< Point > &V) : Region(3) {
    this->vertices = V;
}

void Polygon :: addVertex(const Point &P) {
    this->vertices.push_back(P);
    if ( this->vertices.size() < 3 ) {
        // There must be at least 3 vertices in polygon[]
        std :: cerr << "Error: there must be at least 3 vertices in Polygon" << '\n';
        exit(1);
    }
}

bool Polygon :: isInside(const Point &P) const {
    // Create a point for line segment from p to infinite (ray)
    Point extreme(INF, P.y(), 0.0);

    // Count intersections of the above line with sides of polygon
    int count = 0;
    for ( size_t i = 0; i < this->vertices.size(); i++ ) {
        int next = ( i == this->vertices.size() - 1 ) ? 0 : i + 1;

        // Check if the line segment from 'p' to 'extreme' intersects
        // with the line segment from 'polygon[i]' to 'polygon[next]'
        if ( doIntersect(this->vertices [ i ], this->vertices [ next ], P, extreme) ) {
            // If the point 'p' is colinear with line segment 'i-next',
            // then check if it lies on segment. If it lies, return true,
            // otherwise false
            if ( orientation(this->vertices [ i ], P, this->vertices [ next ]) == 0 ) {
                return onSegment(this->vertices [ i ], P, this->vertices [ next ]);
            }

            count++;
        }
    }

    // Return true if count is odd, false otherwise
    return count & 1; // Same as (count%2 == 1)
}

///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////
Cylinder :: Cylinder(const Point &a, const Point &b, double r) : RegularRegion(3) {
    A = a;
    B = b;
    radius = r;
    init();
}

///////////////////////////////////////////////////////////////////////////////
void Cylinder :: readFromLine(std :: istringstream &iss) {
    double x, y, z;
    string param;
    bool bbcenter = false;
    bool btcenter = false;
    bool brad = false;
    while (  iss >> param ) {
        if ( param.compare("center_bottom") == 0 ) {
            iss >> x >> y >> z;
            A = Point(x, y, z);
            bbcenter = true;
        } else if ( param.compare("center_top") == 0 ) {
            iss >> x >> y >> z;
            B = Point(x, y, z);
            btcenter = true;
        }
        if ( param.compare("radius") == 0 ) {
            iss >> radius;
            brad = true;
        }
    }
    if ( !bbcenter ) {
        cerr << name << " Error: parameter center_bottom not specified" << endl;
        exit(1);
    }
    if ( !btcenter ) {
        cerr << name << " Error: parameter center_top not specified" << endl;
        exit(1);
    }
    if ( !brad ) {
        cerr << name << " Error: parameter radius not specified" << endl;
        exit(1);
    }

    init();
}

///////////////////////////////////////////////////////////////////////////////
void Cylinder :: init() {
    dir = B - A;
    length = dir.norm();
    dir /= length;
}

///////////////////////////////////////////////////////////////////////////////
bool Cylinder :: isInside(const Point &P) const {
    Point vec = P - A;
    double pos = dir [ 0 ] * vec [ 0 ] + dir [ 1 ] * vec [ 1 ] + dir [ 2 ] * vec [ 2 ];
    if ( pos< 0 or pos >length ) {
        return false;
    }
    vec = vec - pos * dir;
    if ( vec.norm() > radius ) {
        return false;
    }
    return true;
}


///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////
bool isInBlock(const Point &P, const Point &leftBottom, const Point &rightTop) {
    if ( ( leftBottom.x() - P.x() ) < TOL && ( P.x() - rightTop.x() ) < TOL ) {
        if ( ( leftBottom.y() - P.y() ) < TOL && ( P.y() - rightTop.y() ) < TOL ) {
            if ( ( leftBottom.z() - P.z() ) < TOL && ( P.z() - rightTop.z() ) < TOL ) {
                return true;
            }
        }
    }
    return false;
}

bool isInCircle(const Point &P, const Point &center, const double &radius,
                const unsigned &dir) {
    Point P1, C1;
    P1 = Point( ( dir == 0 ) ? 0 : P.x(), ( dir == 1 ) ? 0 : P.y(), ( dir == 2 ) ? 0 : P.z() );
    C1 = Point( ( dir == 0 ) ? 0 : center.x(), ( dir == 1 ) ? 0 : center.y(), ( dir == 2 ) ? 0 : center.z() );
    if ( ( P1 - C1 ).norm() < radius ) {
        return true;
    }
    return false;
}


// Given three colinear points p, q, r, the function checks if
// point q lies on line segment 'pr'
bool onSegment(const Point &p, const Point &q, const Point &r)
{
    if ( q.x() <= max( p.x(), r.x() ) && q.x() >= min( p.x(), r.x() ) &&
         q.y() <= max( p.y(), r.y() ) && q.y() >= min( p.y(), r.y() ) ) {
        return true;
    }
    return false;
}

// To find orientation of ordered triplet (p, q, r).
// The function returns following values
// 0 --> p, q and r are colinear
// 1 --> Clockwise
// 2 --> Counterclockwise
int orientation(const Point &p, const Point &q, const Point &r)
{
    double val = ( q.y() - p.y() ) * ( r.x() - q.x() ) -
                 ( q.x() - p.x() ) * ( r.y() - q.y() );

    if ( val == 0 ) {
        return 0;            // colinear
    }
    return ( val > 0 ) ? 1 : 2; // clock or counterclock wise
}

// The function that returns true if line segment 'p1q1'
// and 'p2q2' intersect.
bool doIntersect(const Point &p1, const Point &q1, const Point &p2, const Point &q2)
{
    // Find the four orientations needed for general and
    // special cases
    int o1 = orientation(p1, q1, p2);
    int o2 = orientation(p1, q1, q2);
    int o3 = orientation(p2, q2, p1);
    int o4 = orientation(p2, q2, q1);

    // General case
    if ( o1 != o2 && o3 != o4 ) {
        return true;
    }

    // Special Cases
    // p1, q1 and p2 are colinear and p2 lies on segment p1q1
    if ( o1 == 0 && onSegment(p1, p2, q1) ) {
        return true;
    }

    // p1, q1 and p2 are colinear and q2 lies on segment p1q1
    if ( o2 == 0 && onSegment(p1, q2, q1) ) {
        return true;
    }

    // p2, q2 and p1 are colinear and p1 lies on segment p2q2
    if ( o3 == 0 && onSegment(p2, p1, q2) ) {
        return true;
    }

    // p2, q2 and q1 are colinear and q1 lies on segment p2q2
    if ( o4 == 0 && onSegment(p2, q1, q2) ) {
        return true;
    }

    return false; // Doesn't fall in any of the above cases
}

// Returns true if the point p lies inside the polygon[] with n vertices (vertices must be ordered)
bool isInPolygon(const std :: vector< Point > &polygon, const Point &p)
{
    // There must be at least 3 vertices in polygon[]
    if ( polygon.size() < 3 ) {
        return false;
    }

    // Create a point for line segment from p to infinite (ray)
    Point extreme = { INF, p.y(), 0.0 };

    // Count intersections of the above line with sides of polygon
    int count = 0;
    for ( unsigned i = 0; i < polygon.size(); i++ ) {
        int next = ( i == polygon.size() - 1 ) ? 0 : i + 1;

        // Check if the line segment from 'p' to 'extreme' intersects
        // with the line segment from 'polygon[i]' to 'polygon[next]'
        if ( doIntersect(polygon [ i ], polygon [ next ], p, extreme) ) {
            // If the point 'p' is colinear with line segment 'i-next',
            // then check if it lies on segment. If it lies, return true,
            // otherwise false
            if ( orientation(polygon [ i ], p, polygon [ next ]) == 0 ) {
                return onSegment(polygon [ i ], p, polygon [ next ]);
            }

            count++;
        }
    }

    // Return true if count is odd, false otherwise
    return count & 1;  // Same as (count%2 == 1)
}


///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////
// CONTAINER FOR REGIONS
///////////////////////////////////////////////////////////////////////////////
RegionContainer :: ~RegionContainer() {
    for ( auto &r: regions ) {
        if ( r != nullptr ) {
            delete r;
        }
    }
}

///////////////////////////////////////////////////////////////////////////////
void RegionContainer :: readFromFile(const std :: string &filename, unsigned d) {
    dim = d;
    size_t origsize = regions.size();
    string line, regionType;
    ifstream inputfile( filename.c_str() );
    if ( inputfile.is_open() ) {
        while ( getline(inputfile >> std :: ws, line) ) {
            if ( line.empty() || ( line.at(0) == '#' ) ) {
                continue;
            }
            istringstream iss(line);
            iss >> std :: ws >> regionType;
            if ( !( regionType.rfind("#", 0) == 0 ) ) {
                if ( regionType.compare("Box") == 0 || regionType.compare("Rectangle") == 0 ) {
                    Block *newregion = new Block(dim);
                    newregion->readFromLine(iss);
                    regions.push_back(newregion);
                } else if ( regionType.compare("Circle") == 0 || regionType.compare("Sphere") == 0 ) {
                    Sphere *newregion = new Sphere();
                    newregion->readFromLine(iss);
                    regions.push_back(newregion);
                } else if ( regionType.compare("Cylinder") == 0 ) {
                    Cylinder *newregion = new Cylinder();
                    newregion->readFromLine(iss);
                    regions.push_back(newregion);
                } else {
                    cerr << "Error: region type '" <<  regionType <<  "' does not exists" << endl;
                    exit(EXIT_FAILURE);
                }
            }
        }
        inputfile.close();
        cout << "Input file '" <<  filename << "' succesfully loaded; " << regions.size() - origsize << " regions found" << endl;
    } else {
        cerr << "Error: unable to open input file '" <<  filename <<  "'" << endl;
        exit(EXIT_FAILURE);
    }
}

///////////////////////////////////////////////////////////////////////////////
bool RegionContainer :: isLocationValid(const Point p, const vector< unsigned >in, const vector< unsigned >out) const {
    for ( auto &k: in ) {
        if ( k >= regions.size() ) {
            cerr << "Region number " << k << " does not exists" << endl;
            exit(1);
        }
        if ( !regions [ k ]->isInside(p) ) {
            return false;
        }
    }
    for ( auto &k: out ) {
        if ( k >= regions.size() ) {
            cerr << "Region number " << k << " does not exists" << endl;
            exit(1);
        }
        if ( regions [ k ]->isInside(p) ) {
            return false;
        }
    }
    return true;
}

//TO BE USE IN ADAPTIVITY
///////////////////////////////////////////////////////////////////////////////
void readRegions(const std :: string &filename, std :: vector< std :: unique_ptr< Region > > &regions, unsigned dim) {
    size_t origsize = regions.size();
    string line, regionType;
    ifstream inputfile( filename.c_str() );
    if ( inputfile.is_open() ) {
        while ( getline(inputfile >> std :: ws, line) ) {
            if ( line.empty() || ( line.at(0) == '#' ) ) {
                continue;
            }
            istringstream iss(line);
            iss >> std :: ws >> regionType;
            if ( !( regionType.rfind("#", 0) == 0 ) ) {
                if ( regionType.compare("box") == 0 || regionType.compare("rectangle") == 0 ) {
                    auto newregion = std :: make_unique< Block >(dim);
                    newregion->readFromLine(iss);
                    regions.push_back( std :: move(newregion) );
                } else if ( regionType.compare("circle") == 0 || regionType.compare("sphere") == 0 ) {
                    auto newregion = std :: make_unique< Sphere >();
                    newregion->readFromLine(iss);
                    regions.push_back( std :: move(newregion) );
                } else {
                    cerr << "Error: region type '" <<  regionType <<  "' does not exists" << endl;
                    exit(EXIT_FAILURE);
                }
            }
        }
        inputfile.close();
        cout << "Input file '" <<  filename << "' succesfully loaded; " << regions.size() - origsize << " regions found" << endl;
    } else {
        cerr << "Error: unable to open input file '" <<  filename <<  "'" << endl;
        exit(EXIT_FAILURE);
    }
}



///////////////////////////////////////////////////////////////////////////////
bool isInsideRegions(const std :: vector< std :: unique_ptr< Region > > &regions, const Point &p) {
    for ( auto const &reg : regions ) {
        if ( reg->isInside(p) ) {
            return true;
        }
    }
    return false;
}


bool isInsideRegions(const std :: vector< std :: unique_ptr< Region > > &regions, const Element *el) {
    unsigned inside = 0;
    for ( auto const &reg : regions ) {
        inside = 0;
        for ( auto const &n : el->giveNodes() ) {
            if ( reg->isInside( n->givePoint() ) ) {
                inside++;  // must be in the same region, not in two neighboring
            }
        }
        if ( inside == el->giveNodes().size() ) {
            // return true only if all nodes of the particular elem are in one of the regions
            return true;
        }
    }
    return false;
}
