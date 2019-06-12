#ifndef _NODE_C_H
#define _NODE_C_H

#include "node.h"
#include <vector>
#include <iostream>
#include <fstream>
#include "boundary_condition.h"

class NodeContainer
{
private:
    vector< Node * >nodes;
    vector< unsigned >DoFid;         //mapping from particle order to DoF order for all DoFs
    vector< unsigned >loadedDoFid;     //mapping from particle order to DoF order for loaded DoFs
    vector< unsigned >blockedDoFid;  //mapping from particle order to DoF order for blocked DoFs
    unsigned totalDoFs, freeDoFs;
    BCContainer *BC;

    void establishDoFArray();

public:
    NodeContainer(BCContainer *bc) { BC = bc; totalDoFs = freeDoFs = 0; };
    ~NodeContainer();

    void readFromFile(const string filename, const int dim);
    Node *giveNode(unsigned const num) { return nodes [ num ]; }
    unsigned giveSize() { return nodes.size(); };
    unsigned giveDoFid(unsigned i) const { return DoFid [ i ]; }
    void establishDoFArray(BCContainer *BC);
    unsigned giveTotalNumDoFs() const { return totalDoFs; };
    unsigned giveNumFreeDoFs() const { return freeDoFs; };
    void init();
    void addRHS_nodalLoad(Vector &RHS, double time) const;
    void updateDirrichletBC(Vector &bc, double time) const;
    void giveFullDoFArray(const Vector &fDoFs, const Vector &bDoFs, Vector &fullDoFs) const;
    void giveFullDoFArray(const Vector &fDoFs, Vector &fullDoFs) const;
    void giveReducedDoFArray(const Vector &fullDoFs, Vector &fDoFs) const;
    void updateExteranlForcesByReactions(const Vector &f_int, const Vector &load, Vector &f_ext) const;
    Node *findClosestMechanicalNode(Point A) const;

protected:
};

#endif /* _NODE_C_H */
