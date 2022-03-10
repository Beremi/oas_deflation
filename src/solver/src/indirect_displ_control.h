#ifndef IDC_H
#define IDC_H

#include <vector>
#include <iostream>
#include <fstream>
#include <iomanip>      // std::setw
#include "linalg.h"
#include "element_container.h"
#include "node_container.h"
#include "function.h"

//////////////////////////////////////////////////////////
class IndirectDC
{
protected:
    string name;
    Function *func;
    double target_value;
    int funcnum;
    unsigned nummaxunit;
    vector< bool >coords_active;
    vector< bool >nodes_active;
    vector< vector< unsigned > >c_nodes;
    vector< vector< unsigned > >c_dirs;
    vector< vector< unsigned > >c_DoFs;
    vector< vector< double > >xcoords;
    vector< vector< double > >ycoords;
    vector< vector< double > >zcoords;
    vector< vector< double > >c_weights;

    double givePrescribedDisplacement(double time);

public:
    IndirectDC();
    ~IndirectDC() {};
    void init(NodeContainer *nodes, FunctionContainer *funcs, bool initial = true);
    void readFromStream(unsigned num, ifstream &inputfile);
    double giveMultiplierCorrection(MyVector &prev_displ, MyVector &displ_f, double time);
    double giveControlValue(MyVector &displ);
};

#endif
