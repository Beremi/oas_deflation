#ifndef _MATERIAL_CONTAINER_H
#define _MATERIAL_CONTAINER_H

#include "material.h"
#include <vector>
#include <iostream>
#include <fstream>


class ElementContainer; //forward declaration
class Model; //forward declaration

//////////////////////////////////////////////////////////
class MaterialContainer
{
private:
    std :: vector< Material * >matrs;
    Model *masterModel;
public:
    MaterialContainer() {};
    ~MaterialContainer();

    void readFromFile(const std :: string filename, unsigned dim);
    void init();
    Material *giveMaterial(unsigned const mat);
    void runPreparationForStressEvaluation(ElementContainer *elems);
    Model * giveModel() const;
    void setModel(Model *m);
    bool requestTetrahedralBackgroundMesh() const;
protected:
};

#endif /* _MATERIAL_CONTAINER_H */
