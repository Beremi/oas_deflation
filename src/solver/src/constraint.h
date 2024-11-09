#ifndef _CONSTRAINT_H
#define _CONSTRAINT_H

#include "linalg.h"
#include "function.h"


class Node; //forward declaration
class NodeContainer; //forward declaration
class ElementContainer;  // forward declaration
class ConstraintContainer;  // forward declaration
class BCContainer;  // forward declaration
class Solver;  // forward declaration

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
    virtual ~JointDoF() {};
    JointDoF(Node *s, const unsigned &dir, const std :: vector< Node * > &m, const std :: vector< unsigned > &dirs, const std :: vector< double > &mult, const std :: vector< Function * > &fns = {}, const std :: vector< double > &time_mult = {});
    void readFromLine(std :: istringstream &iss, NodeContainer *nodes);
    void print();
    virtual void init(Solver *solver);
    unsigned giveSlaveDoF() const;
    Node *giveSlaveNode() const { return slaveNode; };
    unsigned giveSlaveDir() const { return direction; };
    std :: vector< Node * >giveMasterNodes() { return masters; };
    virtual unsigned giveNumOfDoFMasters() const { return masters.size(); }; //TODO: warning C4267: 'return': conversion from 'size_t' to 'unsigned int', possible loss of data
    virtual unsigned giveNumOfConjugateMasters() const { return 0; };
    Node *giveMasterNode(unsigned k) { return masters [ k ]; };
    unsigned giveMasterDoF(unsigned k) const;
    std :: vector< unsigned >giveMasterDirs() { return directions; };
    unsigned giveMasterDir(unsigned k) const { return directions [ k ]; };
    std :: vector< double >giveMasterMultipliers() { return multipliers; };
    double giveMasterMultiplier(unsigned k, const double time_now = 0.0) const { ( void ) time_now; return multipliers [ k ]; };
    double giveFnDepPart(unsigned k, const double time_now = 0.0) const {
        return ( time_fns [ k ] == nullptr ) ? 0.0 : additional_term [ k ] * time_fns [ k ]->giveY(time_now);
    };
    std :: vector< Function * >giveTimeFns() { return time_fns; };
    Function *giveTimeFn(unsigned k) const { return time_fns [ k ]; };
    bool isTimeDependent() { return !time_fns.empty(); };
};

//////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////
class VolumetricAverage : public JointDoF
{
protected:
    ElementContainer *elems;
    ConstraintContainer *constraints;
    std :: vector< Node * >nodes;
    std :: vector< unsigned >dirs;
    Node *masternode;
    unsigned masterdir;
public:
    VolumetricAverage(std :: vector< Node * > &n, std :: vector< unsigned > &d, Node *mn, unsigned md, ElementContainer *ec, ConstraintContainer *cc);
    virtual ~VolumetricAverage() {};
    virtual void init(Solver *solver);
};


//////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////
// UPDATE PRESSURE BASED ON TOTAL LOAD
class DoFDependentOnConjugates : public JointDoF
{
protected:
public:
    DoFDependentOnConjugates() {};
    DoFDependentOnConjugates(Node *s, const unsigned &dir, const std :: vector< Node * > &m, const std :: vector< unsigned > &dirs, const std :: vector< double > &mult, const std :: vector< Function * > &fns = {}, const std :: vector< double > &time_mult = {});
    ;
    ~DoFDependentOnConjugates() {};
    virtual unsigned giveNumOfDoFMasters() const { return 0; }; //TODO: warning C4267: 'return': conversion from 'size_t' to 'unsigned int', possible loss of data
    virtual unsigned giveNumOfConjugateMasters() const { return masters.size(); };
    virtual void init(Solver *solver);
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
    CoordinateIndexedSparseMatrix X_full;  // for export
    bool time_dependent = false;

public:
    ConstraintContainer() = default;
    ~ConstraintContainer();
    void readFromFile(const std :: string filename, const unsigned ndim, NodeContainer *nodes);
    // void calculateSlaveDoFfield(NodeContainer *nodes);
    void init(NodeContainer *nodes, BCContainer *bconds, Solver *solver); // here matrix X will be created
    void initFull(NodeContainer *nodes, BCContainer *bconds, Solver *solver); // here matrix X will be created
    CoordinateIndexedSparseMatrix giveFullMatrixX(NodeContainer *nodecont, BCContainer *bccont, Solver *solver);
    void clear();
    void transformToConstraintSpace(CoordinateIndexedSparseMatrix &K, const double time_now = 0);
    void calculateDependentDoFs(Vector &fullDoFs, const double time_now = 0.0, const bool all = false) const;
    void calculateDoFsDependentOnConjugates(Vector &full_ddr, const Vector &fullDoFs, const Vector &fullFExt) const;
    void calculateMasterForces(Vector &fullForces);
    JointDoF *giveConstraint(const unsigned &i) { return constraints [ i ]; };
    void addConstraint(JointDoF *jd) { constraints.push_back(jd); };
    size_t giveSize() { return constraints.size(); };
    bool isActive() const { return !constraints.empty(); }
    void removeConstraint(unsigned i);

    std :: vector< JointDoF * > :: iterator begin() { return constraints.begin(); }
    std :: vector< JointDoF * > :: iterator end() { return constraints.end(); }
    std :: vector< JointDoF * > :: const_iterator begin() const { return constraints.begin(); }
    std :: vector< JointDoF * > :: const_iterator end() const { return constraints.end(); }

protected:
};

#endif /* _CONSTRAINT_H */
