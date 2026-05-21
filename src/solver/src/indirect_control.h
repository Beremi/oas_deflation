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
    bool requiref;
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

    double givePrescribedValue(double time) const;

public:
    IndirectControl();
    virtual ~IndirectControl() {};
    void init(NodeContainer *mnodes, FunctionContainer *funcs, bool initial = true);
    void readFromStream(unsigned num, std :: ifstream &inputfile);
    virtual double giveMultiplierCorrection(Vector &prev_displ, Vector &prev_force, Vector &diff_displ, Vector &diff_force, double time);
    double giveTargetValue(double time) const;
    double giveControlValue(const Vector &displ, const Vector &force, unsigned unit = 0) const;
    double giveControlDifference(const Vector &diff_displ, const Vector &diff_force, unsigned unit = 0) const;
    bool requireForces()const;
};

//////////////////////////////////////////////////////////
class IndirectControlSumOfSquares : public IndirectControl
{
protected:

public:
    IndirectControlSumOfSquares();
    virtual ~IndirectControlSumOfSquares() {};
    virtual double giveMultiplierCorrection(Vector &prev_displ, Vector &prev_force, Vector &diff_displ, Vector &diff_force, double time);
};
#endif /* _INDIRECT_CONTROL_H */
