#ifndef _SIMPLEX_H
#define _SIMPLEX_H

#include "linalg.h"
#include <vector>
#include <iostream>
#include <unordered_set>


class Node; //forward declaration
class Particle; //forward declaration
class RigidBodyContact; //forward declaration

class Simplex
{
protected:
    Node *center;
    std :: vector< Particle * >nodes; ///list of nodes
    bool valid, transport;
    double volume, volstrain;
    double pressure;
    unsigned pressureDoF;
    std :: vector< unsigned >DoFs;  ///corresponding degrees of freedom
    std :: vector< double >DoFweights;  ///corresponding weights
    std :: vector< RigidBodyContact * >elems;
    std :: unordered_set< Simplex * >neighbors;  //neighboring simplexes used to steal volumetric strain

public:
    Simplex(Node *c) { center = c; valid = false; transport = false; volstrain = 0; };
    virtual ~Simplex() {};
    void addElement(RigidBodyContact *rbc);
    void init(unsigned ndim);
    void findNeighbors();
    double giveVolumetricStrain() const;
    double givePressure() const;
    void computeVolumetricStrain(const Vector &DoFs);
    void stealVolumetricStrain(const Vector &DoFs);
    bool isValid() const { return valid; }
    bool hasPressure() const { return transport; }
    std :: vector< RigidBodyContact * >giveElements() { return elems; }
    bool isUsable() const {return (valid || neighbors.size()>0);};
};


#endif /* _SIMPLEX_H */
