#ifndef IDC_H
#define IDC_H

#include <vector>
#include <iostream>
#include <fstream>
#include <iomanip>      // std::setw
#include "linear_algebra.h"
#include "element_container.h"
#include "node_container.h"

//////////////////////////////////////////////////////////
class IndirectDC
{
protected:
    string name;
    ElementContainer *elems;
    NodeContainer *nodes;
    Function *func;
    unsigned nummaxunit;
    vector <unsigned> c_nodes;
    vector <unsigned> c_dirs;
    vector <unsigned> c_DoFs;
    vector <double> c_weights;

public:
    IndirectDC() { name = "indirect displacement controller"; };
    virtual ~IndirectDC() {};
    virtual void init();
    virtual void readFromStream(ifstream &inputfile);
};

#endif
