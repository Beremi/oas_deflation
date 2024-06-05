#ifndef _ELEMENT_SUPERELEM_H
#define _ELEMENT_SUPERELEM_H

#include "element.h"

class ElementContainer; //forward declaration;

//////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////
// AI BASED ELEMENT
class MLMechElement : public Element
{
protected:
    unsigned poly_degree;
    Matrix stiffmat;

    Matrix readMatrixFromFile(std::string filepath) const;

public:
    MLMechElement(unsigned dim);
    virtual ~MLMechElement() {};
    virtual Vector giveStrain(unsigned i, const Vector &DoFs);
    virtual void readFromLine(std :: istringstream &iss, NodeContainer *fullnodes, MaterialContainer *fullmatrs);
    virtual void init();
    virtual Matrix giveStiffnessMatrix(std :: string matrixType) const;
    virtual Matrix giveMassMatrix();
    virtual Matrix giveDampingMatrix();
    virtual Vector giveLumpedMassMatrix();
    virtual Vector giveInternalForces(const Vector &DoFs, bool frozen, double timeStep);
    virtual void giveValues(std :: string code, Vector &result) const;
    virtual Vector integrateInternalSources();
    double giveVolume() const { return volume; };    
};

#endif  /* _ELEMENT_SUPERELEM_H */
