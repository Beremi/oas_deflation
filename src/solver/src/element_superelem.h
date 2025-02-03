#ifndef _ELEMENT_SUPERELEM_H
#define _ELEMENT_SUPERELEM_H

#include "element.h"
#include <torch/script.h>

class ElementContainer; //forward declaration;

//////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////
// AI BASED ELEMENT
class MLMechElement : public Element
{
protected:
    unsigned poly_degree;
    Matrix stiffmat_elastic;
    Matrix stiffmat;
    double max_elastic_strain_energy = 1e4;
    fs :: path sm_path;
    fs :: path nm_path;
    fs :: path nmL_path = "";
    fs :: path ml_path;
    Matrix readStiffMatrixFromFile() const;
    Matrix readDataNormalizationMatrix(int size, fs :: path matrix_path) const;
    torch :: jit :: script :: Module module;
    Matrix norm;
    Matrix normL;
    std :: string Ftype = "Fp"; // F, Fp
    std :: string Ktype = "L_tangent"; // L_tangent, Lp_tangent, K, Kp
    std :: string normalizationType = "default"; // default, bisymLog



public:
    MLMechElement(unsigned dim);
    virtual ~MLMechElement() {};
    virtual Vector giveStrain(unsigned i, const Vector &DoFs);
    virtual void readFromLine(std :: istringstream &iss, NodeContainer *fullnodes, MaterialContainer *fullmatrs);
    virtual void init();
    virtual Matrix giveStiffnessMatrix(std :: string matrixType) const;
    virtual Matrix giveMassMatrix();
    virtual Matrix giveDampingMatrix();
    virtual Matrix giveLumpedMassMatrix();
    virtual Vector giveInternalForces(const Vector &DoFs, bool frozen, double timeStep);
    virtual void giveValues(std :: string code, Vector &result) const;
    virtual Vector integrateInternalSources();
    double giveVolume() const { return volume; };
};

#endif  /* _ELEMENT_SUPERELEM_H */
