#ifndef _BOUNDARY_CONDITION_H
#define _BOUNDARY_CONDITION_H


#include <vector>
#include <iostream>
#include <fstream>
#include <algorithm>
#include <typeinfo>
#include <set>

#include "function.h"
#include "linalg.h"

class Node; //forward declaration
class NodeContainer; //forward declaration
class Element; //forward declaration
class ElementContainer; //forward declaration
class Solver; //forward declaration

//////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////
// DIRICHLET AND NEUMANN BOUNDARY CONDITION
class BoundaryCondition
{
protected:
    Node *node;
    std :: vector< int >dirichBC; //kinematic - pressure BC
    std :: vector< int >neumannBC; //static - flux BC
    std :: vector< Function * >dirichF; //kinematic - pressure BC
    std :: vector< double >initialState;  //to be added to the function values
    std :: vector< Function * >neumannF; //static - flux BC
    std :: vector< double >multipliers; //multipliers of functions
    unsigned blockedDoFNum, loadedDoFNum;
    double beginTime, endTime;
    bool active, addInitialDoFVals;
public:
    BoundaryCondition() { beginTime = -INFINITY, endTime = INFINITY; addInitialDoFVals = false; };
    BoundaryCondition(Node *n, std :: vector< int > &dBC, std :: vector< int > &nBC) : BoundaryCondition() { node = n; dirichBC = dBC; neumannBC = nBC; multipliers.resize(dBC.size(), 1.); };
    BoundaryCondition(Node *n, std :: vector< int > &dBC, std :: vector< int > &nBC, std :: vector< double > &m) : BoundaryCondition(n, dBC, nBC) { multipliers = m; };

    virtual ~BoundaryCondition() {};
    void replaceDirichBC(std :: vector< int > &newdBC) { dirichBC = newdBC; };
    void replaceNeumannBC(std :: vector< int > &newnBC) { neumannBC = newnBC; };
    void readFromLine(std :: istringstream &iss, NodeContainer *nodes);
    void init(FunctionContainer *funcs, double time);
    unsigned giveNumberOfBlockedDoFs() const { return blockedDoFNum; };
    unsigned giveNumberOfLoadedDoFs() const { return loadedDoFNum; };
    std :: vector< unsigned >giveBlockedDoFs() const;
    std :: vector< unsigned >giveLoadedDoFs() const;
    virtual std :: vector< double >giveBlockedDoFValues(double t) const;
    virtual std :: vector< double >giveLoadedDoFValues(double t) const;
    void setMultipliers(std :: vector< double > &m) { multipliers = m; };
    Node *giveNode() { return node; };
    bool isActive() const { return active; };
    double giveBeginTime() const { return beginTime; };
    double giveEndTime() const { return endTime; };
    void setInitialDoFFields(Solver *s);
    void setInitialDoFFields(Vector &DoFs);
};


//////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////
// VOLUME LOAD
class BodyLoad
{
protected:
    std :: vector< Element * >els;
    unsigned timeFunctionNum;
    Function *timeFunction;
    unsigned spatialFunctionNum;
    Function *spatialFunction;
    unsigned dir;
    double beginTime, endTime;
    bool active;

public:
    BodyLoad() { beginTime = -INFINITY, endTime = INFINITY; };
    virtual ~BodyLoad() {};
    virtual void readFromLine(std :: istringstream &iss, ElementContainer *elems);
    virtual void init(FunctionContainer *funcs, double time);
    virtual double giveValue(const Point *xyz, double time);
    virtual std :: vector< double >giveBodyForceDoFValues(double t);
    std :: vector< unsigned >giveArrayOfBodyForceDoFs() const;
    unsigned giveDirection() const { return dir; };
    bool isActive() const { return active; };
    double giveBeginTime() const { return beginTime; };
    double giveEndTime() const { return endTime; };
};

//////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////
// SELF WEIGHT
class SelfWeight : public BodyLoad
{
protected:
    double g;
public:
    SelfWeight() { };
    virtual ~SelfWeight() {};
    virtual void readFromLine(std :: istringstream &iss, ElementContainer *elems);
    virtual void init(FunctionContainer *funcs, double time);
    virtual double giveValue(const Point *xyz, double time);    
    virtual std :: vector< double >giveBodyForceDoFValues(double t);
};


//////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////
// CONTAINER FOR BOUNDARY CONDITIONS
class BCContainer
{
private:
    FunctionContainer *functions;
    std :: vector< BoundaryCondition * >BC;
    
    std :: vector< unsigned >dirichDoFs;
    std :: vector< unsigned >neumannDoFs;

    std :: vector< BodyLoad * >loads;
    std :: set< double >stageTimes;

public:
    BCContainer() { functions = nullptr; };
    virtual ~BCContainer();
    void setContainers(FunctionContainer *f) { functions = f; };
    void setInitialDoFFields(Solver *s);
    void init(double time);
    void clear();
    void readFromFile(const std :: string filename, NodeContainer *nodes, ElementContainer *elements);
    std :: vector< unsigned >giveArrayOfBlockedDoFs() const { return dirichDoFs; };
    std :: vector< unsigned >giveArrayOfLoadedDoFs() const { return neumannDoFs; };
    std :: vector< unsigned >giveArrayOfBodyForceDoFs() const;
    unsigned giveNumBlockedDoFs() const { return dirichDoFs.size(); };//todo: conversion from 'size_t' to 'unsigned int', possible loss of data
    std :: vector< double >giveBlockedDoFValues(double time) const;
    std :: vector< double >giveLoadedDoFValues(double time) const;
    std :: vector< double >giveBodyForceDoFValues(double time);
    BoundaryCondition *giveBC(unsigned i) { return BC [ i ]; };
    void calculateDoFfields();
    size_t giveSize() { return BC.size(); }
    void addBoundaryCondition(BoundaryCondition *bc) { BC.push_back(bc); }
    void removeBoundaryCondition(unsigned i);
    double giveTimeOfNextChange(double time);
protected:
};


#endif /* _BOUNDARY_CONDITION_H */
