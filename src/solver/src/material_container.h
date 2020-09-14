#ifndef _MATERIAL_C_H
#define _MATERIAL_C_H

#include "material.h"
#include "material_RVE.h"
#include "mars_material.h"
#include "fatigue_material.h"
#include <vector>
#include <iostream>
#include <fstream>


//////////////////////////////////////////////////////////
class MaterialContainer
{
private:
    vector< Material * >matrs;
public:
    MaterialContainer() {};
    ~MaterialContainer();

    void readFromFile(const string filename);
    void init();
    Material *giveMaterial(unsigned const mat) { return matrs [ mat ]; }
protected:
};

#endif /* _MATERIAL_C_H */
