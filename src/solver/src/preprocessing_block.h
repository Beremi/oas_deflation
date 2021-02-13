#ifndef _PBLOCK_H
#define _PBLOCK_H

#include "linear_algebra.h"
#include <vector>
#include <iostream>
#include <fstream>
#include <algorithm>
#include <typeinfo>

// #include "boundary_condition.h"
// #include "node_container.h"
// #include "element_container.h"
#include "data_exporter.h"

//////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////
// MASTER CLASS
class PBlock
{
private:
protected:
    unsigned dim;
public:
    PBlock() {};
    virtual ~PBlock() {};
    virtual void apply(NodeContainer *n, ElementContainer *e, BCContainer *b, ConstraintContainer *c, FunctionContainer *funcs, ExporterContainer *ex)  = 0;
    virtual void readFromLine(istringstream &iss, unsigned d) = 0;
};

//////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////
// Mechanical Periodic Boundary Condition on prism
class MechanicalPeriodicBC : public PBlock
{
protected:
    string name;
    vector< double >PUCsize;
    vector< unsigned >masters;
    vector< unsigned >slaves;
    vector< int >strainFunc;
    vector< int >stressFunc;
    double volume;
    bool use_half_gammas;
    unsigned initalNodeNum;

    int volumetricAverageRigidBC; ///< new boundary condition prescribing average value of pressure

    virtual void generateNewDoFs(NodeContainer *nodes);
    virtual void generateConstraints(NodeContainer *nodes, ConstraintContainer *constrs);
    virtual void generateExporters(NodeContainer *nodes, ExporterContainer *ex);
    virtual void readLoading(istringstream &iss);
    virtual void generateRigidBodyBC(NodeContainer *nodes, ElementContainer *elems, BCContainer *bcs, ConstraintContainer *constrs, FunctionContainer *funcs);
public:
    MechanicalPeriodicBC() { name = "MechanicalPeriodicBC"; };
    virtual ~MechanicalPeriodicBC() {};
    virtual void apply(NodeContainer *nodes, ElementContainer *e, BCContainer *bcs, ConstraintContainer *constrs, FunctionContainer *funcs, ExporterContainer *ex);
    virtual void readFromLine(istringstream &iss, unsigned d);
};


//////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////
// Transport Periodic Boundary Condition on prism
class TransportPeriodicBC : public MechanicalPeriodicBC
{
protected:
    vector< int >microscaleSources; ///< sources at nodes due to mohogenized macroscale pressure gradient
    virtual void generateNewDoFs(NodeContainer *nodes);
    virtual void generateConstraints(NodeContainer *nodes, ConstraintContainer *constrs);
    virtual void generateExporters(NodeContainer *nodes, ExporterContainer *ex);
    virtual void readLoading(istringstream &iss);
    virtual void generateRigidBodyBC(NodeContainer *nodes, ElementContainer *elems, BCContainer *bcs, ConstraintContainer *constrs, FunctionContainer *funcs);
public:
    TransportPeriodicBC() { name = "TransportPeriodicBC"; };
    virtual ~TransportPeriodicBC() {};
};

//////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////
// Rigid Plate
class RigidPlate : public PBlock
{
private:
    std :: vector< unsigned >slave_ids;
public:
    RigidPlate() {};
    virtual ~RigidPlate() {};
    virtual void apply(NodeContainer *nodes, ElementContainer *e, BCContainer *bcs, ConstraintContainer *constrs, FunctionContainer *funcs, ExporterContainer *ex);
    virtual void readFromLine(istringstream &iss, unsigned d);
protected:
    bool transport = false;
    unsigned master_id, ndim;
    std :: string which;  ///< which direction to fix (e.g. to leave expansion in perpendicualr direction)
    void checkMechTransport(Node *master);
};

class CoordRigidPlate : public RigidPlate
{
private:
    Point leftBottom, rightTop;
public:
    CoordRigidPlate() {};
    virtual ~CoordRigidPlate() {};
    virtual void apply(NodeContainer *nodes, ElementContainer *e, BCContainer *bcs, ConstraintContainer *constrs, FunctionContainer *funcs, ExporterContainer *ex);
    virtual void readFromLine(istringstream &iss, unsigned d);
protected:
};


// rigid plate constraining nodes in holow cylindric
class RingRigidPlate : public RigidPlate
{
private:
    Point center, axis;
    double r_inner, r_outer;
    unsigned direction;
public:
    RingRigidPlate() {};
    virtual ~RingRigidPlate() {};
    virtual void apply(NodeContainer *nodes, ElementContainer *e, BCContainer *bcs, ConstraintContainer *constrs, FunctionContainer *funcs, ExporterContainer *ex);
    virtual void readFromLine(istringstream &iss, unsigned d);
protected:
};


//////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////
// CONTAINER FOR PREPROCESSOR BLOCKS
class PBlockContainer
{
private:
    vector< PBlock * >blocks;
    NodeContainer *nodes;
    ElementContainer *elems;
    BCContainer *bcs;
    ConstraintContainer *constrs;
    FunctionContainer *funcs;
    ExporterContainer *exporters;
public:
    PBlockContainer() {};
    virtual ~PBlockContainer();
    void readFromFile(const string filename, unsigned dim);
    void setContainers(NodeContainer *n, ElementContainer *e, BCContainer *b, ConstraintContainer *c, FunctionContainer *f, ExporterContainer *ex);
    void apply();
protected:
};


#endif /* _PBLOCK_H */
