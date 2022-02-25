#ifndef _NODE_C_H
#define _NODE_C_H

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
    vector< Node * >nodes;
    vector< unsigned >DoFid;         //mapping from particle order to DoF order for all DoFs
    vector< unsigned >loadedDoFs;
    vector< unsigned >bodyForceDoFs;
    vector< unsigned >blockedDoFid;
    unsigned totalDoFs, freeDoFs, constrDoFs;
    BCContainer *BC;

    vector< bool >mechDoFs, transpDoFs; //tells if the DoF is mechanical or Transport

    // #constraint
    vector< unsigned >constrainedDoFid;  //mapping from particle order to DoF order for constrained DoFs
    ConstraintContainer *constr;

    void establishDoFArray();

public:
    NodeContainer() { totalDoFs = 0; };
    ~NodeContainer();

    void setContainers(BCContainer *bc, ConstraintContainer *c) { constr = c; BC = bc; }
    void readFromFile(const string filename, const int dim);
    void saveToFile(const std :: string &filepath, std :: vector< unsigned > &nodes_to_save) const;
    unsigned giveNodeId(const Node *node) const;
    Node *giveNode(unsigned const num) const;
    size_t giveSize() { return nodes.size(); };
    unsigned giveDoFid(unsigned i) const { return DoFid [ i ]; }
    void establishDoFArray(BCContainer *BC);
    unsigned giveTotalNumDoFs() const { return totalDoFs; };
    unsigned giveNumFreeDoFs() const { return freeDoFs; };
    unsigned giveNumConstrDoFs() const { return constrDoFs; };
    void init();
    void clear();
    void addRHS_nodalLoad(Vector &f, double time) const;
    void updateDirrichletBC(Vector &r, double time) const;
    void giveFullDoFArray(const Vector &fDoFs, Vector &fullDoFs) const;
    void updateFullDoFsByDependenciesOnConjugates(Vector &ddr, const Vector &trial_r, const Vector &f_ext) const; // accounts also for constraints between master and conjugate variables
    void giveReducedForceArray(Vector &fullDoFs, Vector &fDoFs) const;
    void giveReducedDoFArray(const Vector &fullDoFs, Vector &fDoFs) const;
    void updateExternalForcesByReactions(Vector &f_int, const Vector &load, Vector &f_dam, Vector &f_acc, Vector &f_ext) const;
    Node *findClosestMechanicalNode(const Point A, double *distance) const;
    Node *findClosestAuxiliaryNode(const Point A, double *distance) const;
    Node *findClosestTransportNode(const Point A, double *distance) const;
    unsigned giveNodeNumber(const Node *n) const;
    void addNode(Node *n) { n->setID( nodes.size() ); nodes.push_back(n); };
    vector< bool >giveMechDoFsIndicator() { return mechDoFs; }
    vector< bool >giveTranspDoFsIndicator() { return transpDoFs; }
    void initSimplices();
    void updateSimplexVolumetricStrains(const Vector &fullDoFs);

    ConstraintContainer *giveConstraints() const { return constr; };
    Vector readInitialConditions(string initfile) const;

    vector< Node * > :: iterator begin() { return nodes.begin(); }
    vector< Node * > :: iterator end() { return nodes.end(); }
    vector< Node * > :: const_iterator begin() const { return nodes.begin(); }
    vector< Node * > :: const_iterator end() const { return nodes.end(); }

protected:
};

#endif /* _NODE_C_H */
