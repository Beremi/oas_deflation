#ifndef _INDIRECT_CONTROL_H
#define _INDIRECT_CONTROL_H

#include <vector>
#include <iostream>
#include <fstream>
#include <iomanip>      // std::setw
#include "linalg.h"
#include "element_container.h"
#include "node_container.h"
#include "function.h"

//forward declaration 
class Model;

//////////////////////////////////////////////////////////
class IndirectControl
{
protected:
    std :: string name;
    Function *func;
    double target_value;
    int funcnum;
    unsigned nummaxunit;
    std :: vector< bool >coords_active;
    std :: vector< bool >nodes_active;
    std :: vector< std :: vector< unsigned > >nodes;
    std :: vector< std :: vector< unsigned > >dirs;
    std :: vector< std :: vector< unsigned > >DoFs;
    std :: vector< std :: vector< double > >xcoords;
    std :: vector< std :: vector< double > >ycoords;
    std :: vector< std :: vector< double > >zcoords;
    std :: vector< std :: vector< double > >r_weights;
    std :: vector< std :: vector< double > >f_weights;

    double givePrescribedValue(double time);

public:
    IndirectControl();
    ~IndirectControl() {};
    void init(NodeContainer *mnodes, FunctionContainer *funcs, bool initial = true);
    void readFromStream(unsigned num, std :: ifstream &inputfile);
    virtual double giveMultiplierCorrection(Vector &prev_displ, Vector &prev_force, Vector &diff_displ, Vector &diff_force, double time);
};

//////////////////////////////////////////////////////////
class IndirectControlSumOfSquares: public IndirectControl
{
protected:

public:
    IndirectControlSumOfSquares();
    ~IndirectControlSumOfSquares() {};
    virtual double giveMultiplierCorrection(Vector &prev_displ, Vector &prev_force, Vector &diff_displ, Vector &diff_force, double time);
};
#endif /* _INDIRECT_CONTROL_H */
