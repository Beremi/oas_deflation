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
    std :: string name;
    Function *func;
    double target_value;
    int funcnum;
    unsigned nummaxunit;
    std :: vector< bool >coords_active;
    std :: vector< bool >nodes_active;
    std :: vector< std :: vector< unsigned > >c_nodes;
    std :: vector< std :: vector< unsigned > >c_dirs;
    std :: vector< std :: vector< unsigned > >c_DoFs;
    std :: vector< std :: vector< double > >xcoords;
    std :: vector< std :: vector< double > >ycoords;
    std :: vector< std :: vector< double > >zcoords;
    std :: vector< std :: vector< double > >c_weights;

    double givePrescribedDisplacement(double time);

public:
    IndirectDC();
    ~IndirectDC() {};
    void init(NodeContainer *nodes, FunctionContainer *funcs, bool initial = true);
    void readFromStream(unsigned num, std :: ifstream &inputfile);
    double giveMultiplierCorrection( Vector &prev_displ, Vector &displ_f, double time);
    double giveControlValue( Vector &displ);
};

#endif
