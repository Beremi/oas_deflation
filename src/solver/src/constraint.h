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
  std::vector< JointDoF * > Constraints;
  CoordinateIndexedSparseMatrix X;

  void connectSlaveMaster(Node *slave, Node *master, unsigned const &ndim);
public:
  ConstraintContainer() {};
  ~ConstraintContainer() {};
  void readRigidPlate(istringstream &iss, const unsigned ndim, NodeContainer *nodes);
  void readCoordRigidPlate(istringstream &iss, const unsigned ndim, NodeContainer *nodes);
  void readFromFile(const string filename, const unsigned ndim, NodeContainer *nodes);
  // void calculateSlaveDoFfield(NodeContainer *nodes);
  void init(NodeContainer *nodes);  // here matrix X will be created
  void transformToConstraintSpace(CoordinateIndexedSparseMatrix &K);
  void calculateDependentDoFs(Vector &fullDoFs);
  void calculateMasterForces(Vector &fullForces);
  JointDoF* giveConstraint(const unsigned &i){ return Constraints[ i ]; };


  unsigned size() const { return Constraints.size(); }
  bool isActive() const { return !Constraints.empty(); }

  std :: vector< JointDoF * > :: iterator begin(){return Constraints.begin();}
  std :: vector< JointDoF * > :: iterator end(){return Constraints.end();}
  std :: vector< JointDoF * > :: const_iterator begin() const {return Constraints.begin();}
  std :: vector< JointDoF * > :: const_iterator end() const {return Constraints.end();}

protected:

};


#endif /* _CONSTRAINT_H */
