#ifndef _PBLOCK_H
#define _PBLOCK_H

#include "linear_algebra.h"
#include <vector>
#include <iostream>
#include <fstream>
#include <algorithm>
#include <typeinfo>

#include "boundary_condition.h"
#include "node_container.h"
#include "element_container.h"
#include "constraint.h"
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
    unsigned intialNodeNum;

    virtual void genereteNewDoFs(NodeContainer *nodes);
    virtual void genereteConstraints(NodeContainer *nodes, ConstraintContainer *constrs);
    virtual void genereteExporters(NodeContainer *nodes, ExporterContainer *ex);
    virtual void readLoading(istringstream &iss);
    virtual void genereteRigidBodyBC(NodeContainer *nodes, ElementContainer *elems, BCContainer *bcs, ConstraintContainer *constrs, FunctionContainer *funcs);
public:
    MechanicalPeriodicBC() {name = "MechanicalPeriodicBC";};
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
    int volumetricAverageRigidBC; ///< new boundary condition prescribing average value of pressure
    virtual void genereteNewDoFs(NodeContainer *nodes);
    virtual void genereteConstraints(NodeContainer *nodes, ConstraintContainer *constrs);
    virtual void genereteExporters(NodeContainer *nodes, ExporterContainer *ex);
    virtual void readLoading(istringstream &iss);
    virtual void genereteRigidBodyBC(NodeContainer *nodes, ElementContainer *elems, BCContainer *bcs, ConstraintContainer *constrs, FunctionContainer *funcs);
public:
    TransportPeriodicBC() {name = "TransportPeriodicBC";};
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
    unsigned master_id, ndim;
    std :: string which;  ///< which direction to fix (e.g. to leave expansion in perpendicualr direction)
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
