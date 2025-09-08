#ifndef _ELEMENT_BEAM_H
#define _ELEMENT_BEAM_H

#include "element.h"
#include "cross_section.h"
using namespace std;

//////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////
// TIMOSHENKO BEAM ELEMENT
class TimoshenkoBeam3D : public Element
{
protected:
    CrossSection *CS;
    double length;
    unsigned csnum;
    Point normal;
    Point zdir; //reference point of local z axis - local z axis is perpendicular to local x axis and runs through this point (it is the standard local coordinate system where cross section is given by yz - y points right, z points up).
    Matrix R;

    virtual void checkNodeType() const;
    virtual void setIntegrationPointsAndWeights();
    virtual void computeDampingMatrix();
    Matrix giveDMatrix(unsigned ipnum, std :: string matrixType) const;

public:
    TimoshenkoBeam3D(unsigned dim);
    TimoshenkoBeam3D(Node *a, Node *b, Material *m, CrossSection *cs, Point zrefpoint);
    virtual ~TimoshenkoBeam3D() { };
    virtual void readFromLine(std :: istringstream &iss, NodeContainer *fullnodes, MaterialContainer *fullmatrs, CrossSectionContainer *csc);
    void init();
    virtual Matrix giveBMatrix(const Point *x) const;
    virtual Matrix giveHMatrix(const Point *x) const;
    double giveLength() const { return length; }
    //virtual Vector transformToLocal(const Vector &DoFs) const;
    //virtual Vector transformToGlobal(const Vector &DoFs) const;
    virtual void giveValues(std :: string code, Vector &result) const;
    virtual Vector integrateLoad(BodyLoad *vl, double time) const;
    virtual Vector integrateInternalSources();
    Point giveNormal() const { return normal; }
    virtual bool isPointInside(Point *xn, const Point *x) const { ( void ) xn; ( void ) x; return false; };
};

#endif  /* _ELEMENT_BEAM_H */
