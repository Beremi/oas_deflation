// TODO move all functions regarding geometrical operations etc.
// TODO change thsi to class region -> inherited classes circle, polygon etc with same fns isinside etc
#include "geometry.h"
#include <fstream>
// Define Infinite (Using INT_MAX caused overflow problems)
#define INF 10000
#define TOL 1e-9


Block :: Block(const Point &lB, const Point &rT) {
    this->mainPoint = lB;
    this->rightTop = rT;
}

void Block :: readFromLine(istringstream &iss) {
    double x, y, z;
    iss >> x >> y >> z;
    this->mainPoint = Point(x, y, z);
    iss >> x >> y >> z;
    this->rightTop = Point(x, y, z);
}

bool Block :: isInside(const Point &P) const {
    if ( ( this->mainPoint.getX() - P.getX() ) < TOL && ( P.getX() - this->rightTop.getX() ) < TOL ) {
        if ( ( this->mainPoint.getY() - P.getY() ) < TOL && ( P.getY() - this->rightTop.getY() ) < TOL ) {
            if ( ( this->mainPoint.getZ() - P.getZ() ) < TOL && ( P.getZ() - this->rightTop.getZ() ) < TOL ) {
                return true;
            }
        }
    }
    return false;
}

///////////////////////////////////////////////////////////////////////////////
Circle :: Circle(const Point &c, const double &r) {
    this->mainPoint = c;
    this->size = r;
}

void Circle :: readFromLine(istringstream &iss) {
    double x, y, z, r;
    iss >> x >> y >> z >> r;
    this->mainPoint = Point(x, y, z);
    this->size = r;
}

bool Circle :: isInside(const Point &P) const {
    if ( ( P - this->mainPoint ).norm() < this->size ) {
        return true;
    }
    return false;
}


///////////////////////////////////////////////////////////////////////////////
Sphere :: Sphere(const Point &c, const double &r) {
    this->mainPoint = c;
    this->size = r;
}

bool Sphere :: isInside(const Point &P) const {
    if ( ( P - this->mainPoint ).norm() < this->size ) {
        return true;
    }
    return false;
}

///////////////////////////////////////////////////////////////////////////////
Polygon :: Polygon(const std :: vector< Point > &V) {
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
    Point extreme = { INF, P.getY() };

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
bool isInBlock(const Point &P, const Point &leftBottom, const Point &rightTop) {
    if ( ( leftBottom.getX() - P.getX() ) < TOL && ( P.getX() - rightTop.getX() ) < TOL ) {
        if ( ( leftBottom.getY() - P.getY() ) < TOL && ( P.getY() - rightTop.getY() ) < TOL ) {
            if ( ( leftBottom.getZ() - P.getZ() ) < TOL && ( P.getZ() - rightTop.getZ() ) < TOL ) {
                return true;
            }
        }
    }
    return false;
}

bool isInCircle(const Point &P, const Point &center, const double &radius) {
    // TODO make fn isInCilinder (for arbitrary cilinder orientation)
    if ( ( P - center ).norm() < radius ) {
        return true;
    }
    return false;
}


// Given three colinear points p, q, r, the function checks if
// point q lies on line segment 'pr'
bool onSegment(const Point &p, const Point &q, const Point &r)
{
    if ( q.getX() <= max( p.getX(), r.getX() ) && q.getX() >= min( p.getX(), r.getX() ) &&
         q.getY() <= max( p.getY(), r.getY() ) && q.getY() >= min( p.getY(), r.getY() ) ) {
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
    double val = ( q.getY() - p.getY() ) * ( r.getX() - q.getX() ) -
                 ( q.getX() - p.getX() ) * ( r.getY() - q.getY() );

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
    Point extreme = { INF, p.getY() };

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


void readRegions(const std :: string &filename, std :: vector< Region * > &regions) {
    size_t origsize = regions.size();
    string line, regionType;
    ifstream inputfile( filename.c_str() );
    if ( inputfile.is_open() ) {
        while ( getline(inputfile >> std :: ws, line) ) {
            if ( line.empty() ) {
                continue;
            }
            if ( line.at(0) == '#' ) {
                continue;
            }
            istringstream iss(line);
            iss >> std :: ws >> regionType;
            if ( !regionType.rfind("#", 0) == 0 ) {
                if ( regionType.compare("block") == 0 || regionType.compare("rectangle") == 0 ) {
                    Block *newregion = new Block();
                    newregion->readFromLine(iss);
                    regions.push_back(newregion);
                } else if ( regionType.compare("circle") == 0 || regionType.compare("sphere") == 0 ) {
                    Sphere *newregion = new Sphere();
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


bool isInsideRegions(const std :: vector< Region * > &regions, const Point &p) {
    for ( auto const &reg : regions ) {
        if ( reg->isInside(p) ) {
            return true;
        }
    }
    return false;
}


bool isInsideRegions(const std :: vector< Region * > &regions, const Element *el) {
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
