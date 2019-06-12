#ifndef _BC_H
#define _BC_H

#include "linear_algebra.h"
#include <vector>
#include <iostream>
#include <fstream>
#include <algorithm>

class Node; //forward declaration
class NodeContainer; //forward declaration

//////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////
// BASIC FUNCTION - MASTER CLASS
class Function
{
private:
public:
    Function() {};
    virtual ~Function() {};
    virtual double giveY(double t) const  = 0;
    virtual void readFromLine(istringstream &iss) = 0;
protected:
};

//////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////
// PIECE-WISE LINEAR FUNCTION
class PieceWiseLinearFunction : public Function
{
private:
    vector< double >x;
    vector< double >y;
public:
    PieceWiseLinearFunction() {};
    virtual ~PieceWiseLinearFunction() {};
    void readFromLine(istringstream &iss);
    double giveY(double t) const;
protected:
};

//////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////
// CONTAINER FOR FUNCTIONS
class FunctionContainer
{
private:
    vector< Function * >functions;
public:
    FunctionContainer() {};
    virtual ~FunctionContainer();
    void readFromFile(const string filename);
    double giveY(unsigned f, double t) const;
protected:
};


//////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////
// DIRRICHLET AND NEUMANN BOUNDARY CONDITION
class BoundaryCondition
{
private:
    Node *node;
    vector< int >dirrichBC; //kinematic - pressure BC
    vector< int >neumannBC; //static - flux BC
    unsigned blockedDoFNum, loadedDoFNum;
public:
    BoundaryCondition() {};
    BoundaryCondition(Node *n, vector< int >dBC, vector< int >nBC) { node = n; dirrichBC = dBC; neumannBC = nBC; };
    ~BoundaryCondition() {};
    void init();
    unsigned giveNumberOfBlockedDoFs() const { return blockedDoFNum; };
    unsigned giveNumberOfLoadedDoFs() const { return loadedDoFNum; };
    vector< unsigned >giveBlockedDoFs() const;
    vector< unsigned >giveLoadedDoFs() const;
    vector< unsigned >giveBlockedFunctions() const;
    vector< unsigned >giveLoadedFunctions() const;
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
    vector< unsigned >dirrichF;
    vector< unsigned >neumannF;
    vector< unsigned >dirrichDoFs;
    vector< unsigned >neumannDoFs;

public:
    BCContainer(FunctionContainer *f) { functions = f; };
    virtual ~BCContainer();
    void init();
    void readFromFile(const string filename, NodeContainer *nodes);
    vector< unsigned >giveArrayOfBlockedDoFs() const { return dirrichDoFs; };
    vector< unsigned >giveArrayOfLoadedDoFs() const { return neumannDoFs; };
    vector< double >giveBlockedDoFValues(double time) const;
    vector< double >giveLoadedDoFValues(double time) const;
    BoundaryCondition *giveBC(unsigned i) { return BC [ i ]; };
    void calculateDoFfields();
protected:
};


#endif /* _BC_H */
