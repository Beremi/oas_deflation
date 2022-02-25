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
#include "geometry.h"

void connectSlaveMasterRigid(ConstraintContainer *constrs, Node *slave, Node *master, unsigned const &ndim, const vector< bool > &activeDirs, const bool trsp = false);

void connectSlaveMasterExpansion(ConstraintContainer *constrs, Node *slave, Node *master, unsigned const &ndim, const bool trsp = false, Function *fn = nullptr);

void connectSlaveMasterExpansionFLoad(ConstraintContainer *constrs, Node *slave, Node *master, Node *expMaster, unsigned const &ndim);

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
    virtual void apply(NodeContainer *n, ElementContainer *e, BCContainer *b, ConstraintContainer *c, FunctionContainer *funcs, ExporterContainer *ex, MaterialContainer *mats, Solver *solver)  = 0;
    virtual void readFromLine(istringstream &iss, unsigned d) = 0;
};

// TODO JK: regions with separate material / elastic regions
//////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////
// MaterialRegion
class MaterialRegion : public PBlock
{
private:
    Region *reg;
    // Block block;
    bool transport = false;
    unsigned material_id;
public:
    MaterialRegion() { };
    virtual ~MaterialRegion() {};
    virtual void apply(NodeContainer *nodes, ElementContainer *e, BCContainer *bcs, ConstraintContainer *constrs, FunctionContainer *funcs, ExporterContainer *ex, MaterialContainer *mats, Solver *solver);
    virtual void readFromLine(istringstream &iss, unsigned d);
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
    bool nonsymmetric_shear;
    unsigned initalNodeNum;

    int volumetricAverageRigidBC; ///< new boundary condition prescribing average value of pressure

    virtual void generateNewDoFs(NodeContainer *nodes);
    virtual void generateConstraints(NodeContainer *nodes, ConstraintContainer *constrs);
    virtual void generateExporters(NodeContainer *nodes, ExporterContainer *ex);
    virtual void readLoading(istringstream &iss);
    virtual void generateRigidBodyBC(NodeContainer *nodes, ElementContainer *elems, BCContainer *bcs, ConstraintContainer *constrs, FunctionContainer *funcs);
public:
    MechanicalPeriodicBC() { name = "MechanicalPeriodicBC"; nonsymmetric_shear = false; };
    virtual ~MechanicalPeriodicBC() {};
    virtual void apply(NodeContainer *nodes, ElementContainer *e, BCContainer *bcs, ConstraintContainer *constrs, FunctionContainer *funcs, ExporterContainer *ex, MaterialContainer *mats, Solver *solver);
    virtual void readFromLine(istringstream &iss, unsigned d);
    vector< double >giveDimensions() const { return PUCsize; };
    double giveVolume() const;
};

//////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////
// Mechanical Periodic Boundary Condition on prism with Voigt's constraint
class MechanicalPeriodicBCwithVoigtConstraint : public MechanicalPeriodicBC
{
protected:
    virtual void generateConstraints(NodeContainer *nodes, ConstraintContainer *constrs);
    virtual void generateRigidBodyBC(NodeContainer *nodes, ElementContainer *elems, BCContainer *bcs, ConstraintContainer *constrs, FunctionContainer *funcs);
public:
    MechanicalPeriodicBCwithVoigtConstraint() { name = "MechanicalPeriodicBCwithVoigtConstraint"; };
    virtual ~MechanicalPeriodicBCwithVoigtConstraint() {};
};


//////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////
// Mechanical Periodic Boundary Condition on prism with Elastic constraint
class MechanicalPeriodicBCwithElasticConstraint : public MechanicalPeriodicBC
{
protected:
public:
    MechanicalPeriodicBCwithElasticConstraint() { name = "MechanicalPeriodicBCwithElasticConstraint"; };
    virtual ~MechanicalPeriodicBCwithElasticConstraint() {};
    void apply(NodeContainer *nodes, ElementContainer *e, BCContainer *bcs, ConstraintContainer *constrs, FunctionContainer *funcs, ExporterContainer *ex, MaterialContainer *mats, Solver *solver);
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
// Voigt's constraint
class VoigtConstraint : public PBlock
{
    double volume;
    vector< int >strainFunc, stressFunc;
public:
    VoigtConstraint();
    virtual ~VoigtConstraint();
    virtual void apply(NodeContainer *n, ElementContainer *e, BCContainer *b, ConstraintContainer *c, FunctionContainer *funcs, ExporterContainer *ex, MaterialContainer *mats, Solver *solver);
    virtual void readFromLine(istringstream &iss, unsigned d);
};

//////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////
class PressureFromMechanicalLoad : public PBlock
{
    double multiplier;
    unsigned master;
    unsigned direction, masterdirection;
    unsigned materialnum;
    vector< unsigned >trsprtnodes;

public:
    PressureFromMechanicalLoad();
    virtual ~PressureFromMechanicalLoad();
    virtual void apply(NodeContainer *n, ElementContainer *e, BCContainer *b, ConstraintContainer *c, FunctionContainer *funcs, ExporterContainer *ex, MaterialContainer *mats, Solver *solver);
    virtual void readFromLine(istringstream &iss, unsigned d);
};


//////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////
// Rigid Plate
class RigidPlate : public PBlock
{
private:
    std :: vector< unsigned >slave_ids;
public:
    RigidPlate() { };
    virtual ~RigidPlate() {};
    virtual void apply(NodeContainer *nodes, ElementContainer *e, BCContainer *bcs, ConstraintContainer *constrs, FunctionContainer *funcs, ExporterContainer *ex, MaterialContainer *mats, Solver *solver);
    virtual void readFromLine(istringstream &iss, unsigned d);
protected:
    bool transport = false;
    unsigned master_id, ndim;
    std :: string which;  ///< which direction to fix (e.g. to leave expansion in perpendicualr direction)
    vector< bool >activeDirs;
    void checkMechTransport(Node *master);
    void setDirectionToFix(istringstream &iss);
};

//////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////
class CoordRigidPlate : public RigidPlate
{
private:
    Point leftBottom, rightTop;
public:
    CoordRigidPlate() {};
    virtual ~CoordRigidPlate() {};
    virtual void apply(NodeContainer *nodes, ElementContainer *e, BCContainer *bcs, ConstraintContainer *constrs, FunctionContainer *funcs, ExporterContainer *ex, MaterialContainer *mats, Solver *solver);
    virtual void readFromLine(istringstream &iss, unsigned d);
protected:
};

//////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////
// rigid plate constraining nodes in holow cylindric
class RingRigidPlate : public RigidPlate
{
private:
public:
    RingRigidPlate() {};
    virtual ~RingRigidPlate() {};
    virtual void apply(NodeContainer *nodes, ElementContainer *e, BCContainer *bcs, ConstraintContainer *constrs, FunctionContainer *funcs, ExporterContainer *ex, MaterialContainer *mats, Solver *solver);
    virtual void readFromLine(istringstream &iss, unsigned d);
protected:
    Point center, axis;
    double r_inner, r_outer;
    unsigned direction;
};

//////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////
// rigid plate constraining nodes in holow cylindric
class ExpansionRing : public RingRigidPlate
{
private:
    unsigned fn_id;
public:
    ExpansionRing() {};
    virtual ~ExpansionRing() {};
    virtual void apply(NodeContainer *nodes, ElementContainer *e, BCContainer *bcs, ConstraintContainer *constrs, FunctionContainer *funcs, ExporterContainer *ex, MaterialContainer *mats, Solver *solver);
    virtual void readFromLine(istringstream &iss, unsigned d);
protected:
};

//////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////
// rigid plate constraining nodes in holow cylindric
class ExpansionRingDoFLoad : public RingRigidPlate
{
public:
    ExpansionRingDoFLoad() {};
    virtual ~ExpansionRingDoFLoad() {};
    virtual void apply(NodeContainer *nodes, ElementContainer *e, BCContainer *bcs, ConstraintContainer *constrs, FunctionContainer *funcs, ExporterContainer *ex, MaterialContainer *mats, Solver *solver);
    virtual void readFromLine(istringstream &iss, unsigned d);
protected:
    unsigned expansion_master_id;
};

//////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////
// rigid plate constraining nodes in holow cylindric
class ExpansionRingSingleDoFLoad : public ExpansionRingDoFLoad
{
public:
    ExpansionRingSingleDoFLoad() {};
    virtual ~ExpansionRingSingleDoFLoad() {};
    virtual void apply(NodeContainer *nodes, ElementContainer *e, BCContainer *bcs, ConstraintContainer *constrs, FunctionContainer *funcs, ExporterContainer *ex, MaterialContainer *mats, Solver *solv);
    // virtual void readFromLine(istringstream &iss, unsigned d);
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
    MaterialContainer *materials;
    Solver *solver;
public:
    PBlockContainer() {};
    virtual ~PBlockContainer();
    void readFromFile(const string filename, unsigned dim);
    void setContainers(NodeContainer *n, ElementContainer *e, BCContainer *b, ConstraintContainer *c, FunctionContainer *f, ExporterContainer *ex, MaterialContainer *mats, Solver *solver);
    void apply();
    void clear();
    unsigned giveSize() const { return blocks.size(); };
    PBlock *givePBlock(unsigned i) { return blocks [ i ]; };
protected:
};


#endif /* _PBLOCK_H */
