#ifndef _MATERIAL_C_H
#define _MATERIAL_C_H

#include "material.h"
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
    Material *giveMaterial(unsigned const mat);
protected:
};

#endif /* _MATERIAL_C_H */
