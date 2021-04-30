#ifndef _CONSTRAINT_H
#define _CONSTRAINT_H

#include "linear_algebra.h"
#include "function.h"


class Node; //forward declaration
class NodeContainer; //forward declaration
class ElementContainer;  // forward declaration
class ConstraintContainer;  // forward declaration
class BCContainer;  // forward declaration


class JointDoF
{
protected:
    Node *slaveNode;
    unsigned direction;
    std :: vector< Node * >masters;
    std :: vector< unsigned >directions;
    std :: vector< double >multipliers;
    std :: vector< Function * >time_fns;
    std :: vector< double >additional_term;
public:
    JointDoF() {};
    ~JointDoF() {};
    JointDoF(Node *s, const unsigned &dir, const std :: vector< Node * > &m, const std :: vector< unsigned > &dirs, const std :: vector< double > &mult, const std :: vector< Function * > &fns={}, const std :: vector< double > &time_mult={} );
    void readFromLine(istringstream &iss, NodeContainer *nodes);
    void print();
    virtual void init();
    unsigned giveSlaveDoF() const;
    Node *giveSlaveNode() const { return slaveNode; };
    unsigned giveSlaveDir() const { return direction; };
    std :: vector< Node * >giveMasterNodes() { return masters; };
    unsigned giveNumOfMasters() const { return masters.size(); }; //TODO: warning C4267: 'return': conversion from 'size_t' to 'unsigned int', possible loss of data
    Node *giveMasterNode(unsigned k) { return masters [ k ]; };
    unsigned giveMasterDoF(unsigned k) const;
    std :: vector< unsigned >giveMasterDirs() { return directions; };
    unsigned giveMasterDir(unsigned k) const { return directions [ k ]; };
    std :: vector< double >giveMasterMultipliers() { return multipliers; };
    double giveMasterMultiplier(unsigned k, const double time_now=0.0) const { return multipliers [ k ]; };
    double giveFnDepPart(unsigned k, const double time_now=0.0) const {
        return (time_fns[k] == nullptr) ? 0.0 : additional_term[k] * time_fns[k]->giveY(time_now);
    };
    std :: vector< Function * >giveTimeFns() { return time_fns; };
    Function * giveTimeFn(unsigned k) const { return time_fns [ k ]; };
    bool isTimeDependent() { return !time_fns.empty(); };
};

//////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////
class VolumetricAverage : public JointDoF
{
protected:
    ElementContainer *elems;
    ConstraintContainer *constraints;
    vector< Node * >nodes;
    vector< unsigned >dirs;
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
    NodeContainer *nodes;
    BCContainer *bconds;
    std :: vector< JointDoF * >constraints;
    CoordinateIndexedSparseMatrix X;  // for connection due to geometry
    bool time_dependent = false;

public:
    ConstraintContainer() {};
    ~ConstraintContainer() {};
    void readFromFile(const string filename, const unsigned ndim, NodeContainer *nodes);
    // void calculateSlaveDoFfield(NodeContainer *nodes);
    void init(NodeContainer *nodes, BCContainer *bconds); // here matrix X will be created
    void transformToConstraintSpace(CoordinateIndexedSparseMatrix &K, const double time_now=0);
    void calculateDependentDoFs(Vector &fullDoFs, const double time_now=0.0, const bool all=false);
    void calculateMasterForces(Vector &fullForces);
    JointDoF *giveConstraint(const unsigned &i) { return constraints [ i ]; };
    void addConstraint(JointDoF *jd) { constraints.push_back(jd); };
    size_t giveSize() { return constraints.size(); };
    bool isActive() const { return !constraints.empty(); }

    std :: vector< JointDoF * > :: iterator begin() { return constraints.begin(); }
    std :: vector< JointDoF * > :: iterator end() { return constraints.end(); }
    std :: vector< JointDoF * > :: const_iterator begin() const { return constraints.begin(); }
    std :: vector< JointDoF * > :: const_iterator end() const { return constraints.end(); }

protected:
};

#endif /* _CONSTRAINT_H */
