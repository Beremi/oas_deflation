#ifndef _CONSTRAINT_H
#define _CONSTRAINT_H

#include "linear_algebra.h"

class Node; //forward declaration
class NodeContainer; //forward declaration

class ConstraintContainer;  // forward declaration

class JointDoF
{
private:
  Node * slaveNode;
  unsigned direction;
  std::vector< Node * > masters;
  std::vector< unsigned > directions;
  std::vector< double > multipliers;
public:
  JointDoF() {};
  ~JointDoF() {};
  JointDoF(Node * s, unsigned &dir, std::vector< Node * > &m, std::vector< unsigned > &dirs, std::vector< double > &mult);
  void readFromLine(istringstream &iss, NodeContainer *nodes);
  void print();
  // void init();
  unsigned giveSlaveDoF() const ;
  Node * giveSlaveNode() const { return slaveNode; };
  unsigned giveSlaveDir() const { return direction; };
  std::vector< Node * > giveMasterDoFs(){ return masters; };
  std::vector< unsigned > giveDirs(){ return directions; };
  std::vector< double > giveMultipliers(){ return multipliers; };
protected:

};

class ConstraintContainer
{
private:
  std::vector< JointDoF * > constraints;
  CoordinateIndexedSparseMatrix X;

public:
  ConstraintContainer() {};
  ~ConstraintContainer() {};
  [[deprecated("rigid plates moved to preprocessing_block, specify them in any of PBlockFiles instead of ConstrFiles")]]
  void readRigidPlate(istringstream &iss, const unsigned ndim, NodeContainer *nodes);
  [[deprecated("rigid plates moved to preprocessing_block, specify them in any of PBlockFiles instead of ConstrFiles")]]
  void readCoordRigidPlate(istringstream &iss, const unsigned ndim, NodeContainer *nodes);
  void readFromFile(const string filename, const unsigned ndim, NodeContainer *nodes);
  // void calculateSlaveDoFfield(NodeContainer *nodes);
  void init(NodeContainer *nodes);  // here matrix X will be created
  void transformToConstraintSpace(CoordinateIndexedSparseMatrix &K);
  void calculateDependentDoFs(Vector &fullDoFs);
  void calculateMasterForces(Vector &fullForces);
  JointDoF* giveConstraint(const unsigned &i){ return constraints[ i ]; };
  void addConstraint(JointDoF *jd){ constraints.push_back(jd); };
  void connectSlaveMaster(Node *slave, Node *master, unsigned const &ndim);
  size_t giveSize(){ return constraints.size(); };
  bool isActive() const { return !constraints.empty(); }

  std :: vector< JointDoF * > :: iterator begin(){return constraints.begin();}
  std :: vector< JointDoF * > :: iterator end(){return constraints.end();}
  std :: vector< JointDoF * > :: const_iterator begin() const {return constraints.begin();}
  std :: vector< JointDoF * > :: const_iterator end() const {return constraints.end();}

protected:

};

bool isInBlock(const Point &P, const Point &leftBottom, const Point &rightTop);

#endif /* _CONSTRAINT_H */
