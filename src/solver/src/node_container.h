#ifndef _NODE_CONTAINER_H
#define _NODE_CONTAINER_H

#include "node.h"
#include <vector>
#include <iostream>
#include <fstream>
#include "boundary_condition.h"
#include "constraint.h"

class Solver; //forward declaration;

class NodeContainer
{
private:
    std :: vector< Node * >nodes;
    std :: vector< unsigned >DoFid;         //mapping from full DoF order to free DoF order for all DoFs
    std :: vector< unsigned >invDoFid;      //mapping from free DoF order to full DoF order for free DoFs
    std :: vector< unsigned >DoF2nodes;     //mapping from full DoF order to nodes
    std :: vector< unsigned >loadedDoFs;
    std :: vector< unsigned >bodyForceDoFs;
    std :: vector< unsigned >blockedDoFid;
    unsigned totalDoFs, freeDoFs, constrDoFs;
    BCContainer *BC = nullptr;

    std :: vector< unsigned >physicalFieldsDoF;  //tells which physical field is involved
    std :: vector< unsigned >constrainedDoFid;  //mapping from particle order to DoF order for constrained DoFs
    ConstraintContainer *constr;

    void establishDoFArray();

public:
    NodeContainer() { totalDoFs = 0; };
    ~NodeContainer();

    void setContainers(BCContainer *bc, ConstraintContainer *c) { constr = c; BC = bc; }
    void readFromFile(const std :: string filename, const int dim);
    void saveToFile(const std :: string &filepath, std :: vector< unsigned > &nodes_to_save) const;
    unsigned giveNodeId(const Node *node) const;
    Node *giveNode(unsigned const num) const;
    size_t giveSize() { return nodes.size(); };
    unsigned giveDoFid(unsigned i) const { return DoFid [ i ]; }
    unsigned giveInvDoFid(unsigned i) const { return invDoFid [ i ]; }
    unsigned giveNodeNumberOfDoFID(unsigned i) const { return DoF2nodes [ i ]; };
    Node *giveNodePointerOfDoFID(unsigned i) const { return nodes [ DoF2nodes [ i ] ]; };
    void establishDoFArray(BCContainer *BC);
    unsigned giveTotalNumDoFs() const { return totalDoFs; };
    unsigned giveNumFreeDoFs() const { return freeDoFs; };
    unsigned giveNumConstrDoFs() const { return constrDoFs; };
    std :: vector< unsigned >givePhysicalFieldsOfDoFs() const { return physicalFieldsDoF; };
    unsigned givePhysicalFieldOfDoF(unsigned i) const { return physicalFieldsDoF[i]; };    
    void init();
    void clear();
    void addRHS_nodalLoad(Vector &f, double time) const;
    void updateDirrichletBC(Vector &r, double time) const;
    void giveFullDoFArray(const Vector &fDoFs, Vector &fullDoFs) const;
    void updateFullDoFsByDependenciesOnConjugates(Vector &ddr, const Vector &trial_r, const Vector &f_ext) const;  // accounts also for constraints between master and conjugate variables
    void giveReducedForceArray(Vector &fullDoFs, Vector &fDoFs) const;
    void giveReducedDoFArray(const Vector &fullDoFs, Vector &fDoFs) const;
    void updateExternalForcesByReactions(Vector &f_int, Vector &load, Vector &f_dam, Vector &f_acc, Vector &f_ext, const Vector &full_r) const;
    Node *findClosestMechanicalNode(const Point A, double *distance) const;
    Node *findClosestAuxiliaryNode(const Point A, double *distance) const;
    Node *findClosestTransportNode(const Point A, double *distance) const;
    unsigned giveNodeNumber(const Node *n) const;
    void addNode(Node *n);
    void initSimplices();
    void updateSimplexVolumetricStrains(const Vector &fullDoFs);

    ConstraintContainer *giveConstraints() const { return constr; };
    Vector readInitialConditions(std :: string initfile) const;

    unsigned giveNumOfPhysicalFieldNodes(unsigned i) const;

    std :: vector< Node * > :: iterator begin() { return nodes.begin(); }
    std :: vector< Node * > :: iterator end() { return nodes.end(); }
    std :: vector< Node * > :: const_iterator begin() const { return nodes.begin(); }
    std :: vector< Node * > :: const_iterator end() const { return nodes.end(); }

protected:
};

#endif /* _NODE_CONTAINER_H */
