#ifndef _BC_H
#define _BC_H


#include <vector>
#include <iostream>
#include <fstream>
#include <algorithm>
#include <typeinfo>

#include "function.h"
#include "linalg.h"

class Node; //forward declaration
class NodeContainer; //forward declaration
class Element; //forward declaration
class ElementContainer; //forward declaration

//////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////
// DIRICHLET AND NEUMANN BOUNDARY CONDITION
class BoundaryCondition
{
protected:
    Node *node;
    vector< int >dirichBC; //kinematic - pressure BC
    vector< int >neumannBC; //static - flux BC
    vector< Function * >dirichF; //kinematic - pressure BC
    vector< Function * >neumannF; //static - flux BC
    vector< double >multipliers; //multipliers of functions
    unsigned blockedDoFNum, loadedDoFNum;
public:
    BoundaryCondition() {};
    BoundaryCondition(Node *n, vector< int > &dBC, vector< int > &nBC, vector< double > &m) { node = n; dirichBC = dBC; neumannBC = nBC; multipliers = m; };
    BoundaryCondition(Node *n, vector< int > &dBC, vector< int > &nBC) { node = n; dirichBC = dBC; neumannBC = nBC; multipliers.resize(dBC.size(), 1.); };
    virtual ~BoundaryCondition() {};
    void replaceDirichBC(vector< int > &newdBC) { dirichBC = newdBC; };
    void replaceNeumannBC(vector< int > &newnBC) { neumannBC = newnBC; };
    void readFromLine(istringstream &iss, NodeContainer *nodes);
    void init(FunctionContainer *funcs);
    unsigned giveNumberOfBlockedDoFs() const { return blockedDoFNum; };
    unsigned giveNumberOfLoadedDoFs() const { return loadedDoFNum; };
    vector< unsigned >giveBlockedDoFs() const;
    vector< unsigned >giveLoadedDoFs() const;
    virtual vector< double >giveBlockedDoFValues(double t) const;
    virtual vector< double >giveLoadedDoFValues(double t) const;
    void setMultipliers(vector< double > &m) { multipliers = m; };
    Node *giveNode() { return node; };
};


//////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////
// VOLUME LOAD
class BodyLoad
{
protected:
    vector< Element * >els;
    unsigned timeFunctionNum;
    Function *timeFunction;
    unsigned spatialFunctionNum;
    Function *spatialFunction;
    unsigned dir;

public:
    BodyLoad() {};
    virtual ~BodyLoad() {};
    void readFromLine(istringstream &iss, ElementContainer *elems);
    void init(FunctionContainer *funcs);
    double giveValue(const Point *xyz, double time);
    vector< double >giveBodyForceDoFValues(double t);
    vector< unsigned >giveArrayOfBodyForceDoFs() const;
    unsigned giveDirection() const { return dir; };
};


//////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////
// CONTAINER FOR BOUNDARY CONDITIONS
class BCContainer
{
private:
    FunctionContainer *functions;
    vector< BoundaryCondition * >BC;
    vector< unsigned >dirichDoFs;
    vector< unsigned >neumannDoFs;

    vector< BodyLoad * >loads;

public:
    BCContainer() { functions = nullptr; };
    virtual ~BCContainer();
    void setContainers(FunctionContainer *f) { functions = f; };
    void init();
    void clear();
    void readFromFile(const string filename, NodeContainer *nodes, ElementContainer *elements);
    vector< unsigned >giveArrayOfBlockedDoFs() const { return dirichDoFs; };
    vector< unsigned >giveArrayOfLoadedDoFs() const { return neumannDoFs; };
    vector< unsigned >giveArrayOfBodyForceDoFs() const;
    unsigned giveNumBlockedDoFs() const { return dirichDoFs.size(); };//todo: conversion from 'size_t' to 'unsigned int', possible loss of data
    vector< double >giveBlockedDoFValues(double time) const;
    vector< double >giveLoadedDoFValues(double time) const;
    vector< double >giveBodyForceDoFValues(double time);
    BoundaryCondition *giveBC(unsigned i) { return BC [ i ]; };
    void calculateDoFfields();
    size_t giveSize() { return BC.size(); }
    void addBoundaryCondition(BoundaryCondition *bc) { BC.push_back(bc); }
    void removeBoundaryCondition(unsigned i);
protected:
};


#endif /* _BC_H */
