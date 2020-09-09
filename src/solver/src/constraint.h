#ifndef _CONSTRAINT_H
#define _CONSTRAINT_H

#include "linear_algebra.h"


class Node; //forward declaration
class NodeContainer; //forward declaration
class ElementContainer;  // forward declaration
class ConstraintContainer;  // forward declaration


class JointDoF
{
protected:
    Node *slaveNode;
    unsigned direction;
    std :: vector< Node * >masters;
    std :: vector< unsigned >directions;
    std :: vector< double >multipliers;
public:
    JointDoF() {};
    ~JointDoF() {};
    JointDoF(Node *s, unsigned &dir, std :: vector< Node * > &m, std :: vector< unsigned > &dirs, std :: vector< double > &mult);
    void readFromLine(istringstream &iss, NodeContainer *nodes);
    void print();
    virtual void init();
    unsigned giveSlaveDoF() const;
    Node *giveSlaveNode() const { return slaveNode; };
    unsigned giveSlaveDir() const { return direction; };
    std :: vector< Node * >giveMasterNodes() { return masters; };
    unsigned giveNumOfMasters() const { return masters.size(); };
    Node *giveMasterNode(unsigned k) { return masters[k]; };
    unsigned giveMasterDoF(unsigned k) const;
    std :: vector< unsigned >giveMasterDirs() { return directions; };
    unsigned giveMasterDir(unsigned k) const {return directions[k];}
    std :: vector< double >giveMasterMultipliers() { return multipliers; };
    double giveMasterMultiplier(unsigned k) const {return multipliers[k];}
};

//////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////
class VolumetricAverage : public JointDoF
{
protected:
    ElementContainer *elems;
    ConstraintContainer *constraints;
    vector< Node * > nodes;
    vector< unsigned > dirs;
    Node *masternode;
    unsigned masterdir;
public:
    VolumetricAverage(vector< Node * > &n, std :: vector< unsigned > &d, Node *mn, unsigned md, ElementContainer *ec, ConstraintContainer *cc);
    ~VolumetricAverage() {};
    void init();
};

//////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////
class ConstraintContainer
{
private:
    std :: vector< JointDoF * >constraints;
    CoordinateIndexedSparseMatrix X;

public:
    ConstraintContainer() {};
    ~ConstraintContainer() {};
    void readFromFile(const string filename, const unsigned ndim, NodeContainer *nodes);
    // void calculateSlaveDoFfield(NodeContainer *nodes);
    void init(NodeContainer *nodes); // here matrix X will be created
    void transformToConstraintSpace(CoordinateIndexedSparseMatrix &K);
    void calculateDependentDoFs(Vector &fullDoFs);
    void calculateMasterForces(Vector &fullForces);
    JointDoF *giveConstraint(const unsigned &i) { return constraints [ i ]; };
    void addConstraint(JointDoF *jd) { constraints.push_back(jd); };
    void connectSlaveMaster(Node *slave, Node *master, unsigned const &ndim, const string &which);
    size_t giveSize() { return constraints.size(); };
    bool isActive() const { return !constraints.empty(); }

    std :: vector< JointDoF * > :: iterator begin() { return constraints.begin(); }
    std :: vector< JointDoF * > :: iterator end() { return constraints.end(); }
    std :: vector< JointDoF * > :: const_iterator begin() const { return constraints.begin(); }
    std :: vector< JointDoF * > :: const_iterator end() const { return constraints.end(); }

protected:
};

bool isInBlock(const Point &P, const Point &leftBottom, const Point &rightTop);

#endif /* _CONSTRAINT_H */
