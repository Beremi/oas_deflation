#ifndef _NODE_C_H
#define _NODE_C_H

#include "node.h"
#include <vector>
#include <iostream>
#include <fstream>
#include "boundary_condition.h"
#include "constraint.h"

class NodeContainer
{
private:
    vector< Node * >nodes;
    vector< unsigned >DoFid;         //mapping from particle order to DoF order for all DoFs
    vector< unsigned >loadedDoFs;
    vector< unsigned >bodyForceDoFs;
    vector< unsigned >blockedDoFid;
    unsigned totalDoFs, freeDoFs;
    BCContainer *BC;

    vector< bool >mechDoFs, transpDoFs; //tells if the DoF is mechanical or Transport

    // #constraint
    size_t constrDoFs;
    vector< unsigned >constrainedDoFid;  //mapping from particle order to DoF order for constrained DoFs
    ConstraintContainer *constr;

    void establishDoFArray();

public:
    NodeContainer() { totalDoFs = 0; };
    ~NodeContainer();

    void setContainers(BCContainer *bc, ConstraintContainer *c) { constr = c; BC = bc; }
    void readFromFile(const string filename, const int dim);
    void saveToFile(const std :: string &filepath, std :: vector< unsigned > &nodes_to_save) const;
    unsigned giveNodeId(const Node * node) const;
    Node *giveNode(unsigned const num) const { return nodes [ num ]; }
    size_t giveSize() { return nodes.size(); };
    unsigned giveDoFid(unsigned i) const { return DoFid [ i ]; }
    void establishDoFArray(BCContainer *BC);
    unsigned giveTotalNumDoFs() const { return totalDoFs; };
    unsigned giveNumFreeDoFs() const { return freeDoFs; };
    void init();
    void addRHS_nodalLoad(Vector &RHS, double time) const;
    void updateDirrichletBC(Vector &bc, double time) const;
    void giveFullDoFArray(const Vector &fDoFs, Vector &fullDoFs) const;
    void giveReducedForceArray(Vector &fullDoFs, Vector &fDoFs) const;
    void giveReducedDoFArray(const Vector &fullDoFs, Vector &fDoFs) const;
    void updateExternalForcesByReactions(Vector &f_int, const Vector &load, Vector &f_dam, Vector &f_acc, Vector &f_ext) const;
    Node *findClosestMechanicalNode(Point A) const;
    void addNode(Node *n) { nodes.push_back(n); };
    vector< bool >giveMechDoFsIndicator() { return mechDoFs; }
    vector< bool >giveTranspDoFsIndicator() { return transpDoFs; }

    size_t giveNumConstrDoFs() const { return constrDoFs; };
    ConstraintContainer *giveConstraints() const { return constr; };

    vector< Node * > :: iterator begin() { return nodes.begin(); }
    vector< Node * > :: iterator end() { return nodes.end(); }
    vector< Node * > :: const_iterator begin() const { return nodes.begin(); }
    vector< Node * > :: const_iterator end() const { return nodes.end(); }

protected:
};

#endif /* _NODE_C_H */
