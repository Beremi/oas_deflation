#ifndef _BC_H
#define _BC_H


#include <vector>
#include <iostream>
#include <fstream>
#include <algorithm>
#include <typeinfo>

#include "function.h"
#include "linear_algebra.h"

class Node; //forward declaration
class NodeContainer; //forward declaration

//////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////
// DIRICHLET AND NEUMANN BOUNDARY CONDITION
class BoundaryCondition
{
private:
    Node *node;
    vector< int >dirichBC; //kinematic - pressure BC
    vector< int >neumannBC; //static - flux BC
    vector<double>multipliers; //multipliers of functions
    unsigned blockedDoFNum, loadedDoFNum;
public:
    BoundaryCondition() {};
    BoundaryCondition(Node *n, vector< int >dBC, vector< int >nBC, vector< double >m) { node = n; dirichBC = dBC; neumannBC = nBC; multipliers = m;};
    BoundaryCondition(Node *n, vector< int >dBC, vector< int >nBC) { node = n; dirichBC = dBC; neumannBC = nBC; multipliers.resize(dBC.size(),1.);};
    ~BoundaryCondition() {};
    void init();
    unsigned giveNumberOfBlockedDoFs() const { return blockedDoFNum; };
    unsigned giveNumberOfLoadedDoFs() const { return loadedDoFNum; };
    vector< unsigned >giveBlockedDoFs() const;
    vector< unsigned >giveLoadedDoFs() const;
    vector< unsigned >giveBlockedFunctions() const;
    vector< unsigned >giveLoadedFunctions() const;
    vector< double >giveBlockedMultipliers() const;
    vector< double >giveLoadedMultipliers() const;
    Node *giveNode() { return node; };

protected:
};

//////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////
// CONTAINER FOR BOUNDARY CONDITIONS
class BCContainer
{
private:
    FunctionContainer *functions;
    vector< BoundaryCondition * >BC;
    vector< unsigned >dirichF;
    vector< unsigned >neumannF;
    vector< unsigned >dirichDoFs;
    vector< unsigned >neumannDoFs;
    vector< double >dirichMultipliers;
    vector< double >neumannMultipliers;


public:
    BCContainer(FunctionContainer *f) { functions = f; };
    virtual ~BCContainer();
    void init();
    void readFromFile(const string filename, NodeContainer *nodes);
    vector< unsigned >giveArrayOfBlockedDoFs() const { return dirichDoFs; };
    vector< unsigned >giveArrayOfLoadedDoFs() const { return neumannDoFs; };
    vector< double >giveBlockedDoFValues(double time) const;
    vector< double >giveLoadedDoFValues(double time) const;
    BoundaryCondition *giveBC(unsigned i) { return BC [ i ]; };
    void calculateDoFfields();  
    unsigned giveSize(){return BC.size();}
    void addBoundaryCondition(BoundaryCondition *bc){BC.push_back(bc);}
protected:
};


#endif /* _BC_H */
