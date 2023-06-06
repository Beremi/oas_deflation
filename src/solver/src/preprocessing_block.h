#ifndef _PREPROCESSING_BLOCK_H
#define _PREPROCESSING_BLOCK_H

#include "linalg.h"
#include <vector>
#include <iostream>
#include <fstream>
#include <algorithm>
#include <typeinfo>

// #include "boundary_condition.h"
// #include "node_container.h"
// #include "element_container.h"
#include "data_exporter.h"
#include "element_discrete.h"
#include "geometry.h"

void connectSlaveMasterRigid(ConstraintContainer *constrs, Node *slave, Node *master, unsigned const &ndim, const std :: vector< bool > &activeDirs, bool includeRigidBodyRotation = true);

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
    std :: string name;
    std :: vector< unsigned >insideRegions;
    std :: vector< unsigned >outsideRegions;
public:
    PBlock() {};
    virtual ~PBlock() {};
    virtual void apply(NodeContainer *nodes, ElementContainer *elems, BCContainer *bcs, ConstraintContainer *constrs, FunctionContainer *funcs, ExporterContainer *ex, MaterialContainer *mats, RegionContainer *regions, Solver *solv)  = 0;
    virtual void readFromLine(std :: istringstream &iss, unsigned d) = 0;
    std :: string giveName() const { return name; }
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
    bool all_nodes = false;  ///> only elements with all nodes inside the region are included if true
    unsigned material_id;
public:
    MaterialRegion() { };
    virtual ~MaterialRegion() {};
    virtual void apply(NodeContainer *nodes, ElementContainer *elems, BCContainer *bcs, ConstraintContainer *constrs, FunctionContainer *funcs, ExporterContainer *ex, MaterialContainer *mats, RegionContainer *regions, Solver *solv);
    virtual void readFromLine(std :: istringstream &iss, unsigned d);
};
//////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////
// Mechanical Periodic Boundary Condition on prism
class MechanicalPeriodicBC : public PBlock
{
protected:
    std :: vector< double >PUCsize;
    std :: vector< unsigned >masters;
    std :: vector< unsigned >slaves;
    std :: vector< int >strainFunc;
    std :: vector< int >stressFunc;
    double volume;
    bool nonsymmetric_shear;
    unsigned initalNodeNum;

    int volumetricAverageRigidBC; ///< new boundary condition prescribing average value of pressure

    virtual void generateNewDoFs(NodeContainer *nodes);
    virtual void generateConstraints(NodeContainer *nodes, ConstraintContainer *constrs);
    virtual void generateExporters(NodeContainer *nodes, ExporterContainer *ex);
    virtual void readLoading(std :: istringstream &iss);
    virtual void generateRigidBodyBC(NodeContainer *nodes, ElementContainer *elems, BCContainer *bcs, ConstraintContainer *constrs, FunctionContainer *funcs);
public:
    MechanicalPeriodicBC() { name = "MechanicalPeriodicBC"; nonsymmetric_shear = false; };
    virtual ~MechanicalPeriodicBC() {};
    virtual void apply(NodeContainer *nodes, ElementContainer *elems, BCContainer *bcs, ConstraintContainer *constrs, FunctionContainer *funcs, ExporterContainer *ex, MaterialContainer *mats, RegionContainer *regions, Solver *solv);
    virtual void readFromLine(std :: istringstream &iss, unsigned d);
    std :: vector< double >giveDimensions() const { return PUCsize; };
    double giveVolume() const;
    virtual void calculateVolume();
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
    void apply(NodeContainer *nodes, ElementContainer *elems, BCContainer *bcs, ConstraintContainer *constrs, FunctionContainer *funcs, ExporterContainer *ex, MaterialContainer *mats, RegionContainer *regions, Solver *solv);
};

//////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////
// Transport Periodic Boundary Condition on prism
class TransportPeriodicBC : public MechanicalPeriodicBC
{
protected:
    std :: vector< int >microscaleSources; ///< sources at nodes due to mohogenized macroscale pressure gradient
    virtual void generateNewDoFs(NodeContainer *nodes);
    virtual void generateConstraints(NodeContainer *nodes, ConstraintContainer *constrs);
    virtual void generateExporters(NodeContainer *nodes, ExporterContainer *ex);
    virtual void readLoading(std :: istringstream &iss);
    virtual void generateRigidBodyBC(NodeContainer *nodes, ElementContainer *elems, BCContainer *bcs, ConstraintContainer *constrs, FunctionContainer *funcs);
public:
    TransportPeriodicBC() { name = "TransportPeriodicBC"; };
    virtual ~TransportPeriodicBC() {};
};

//////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////
// Transport Periodic Boundary Condition on sphere
class MechanicalSphericalPeriodicBC : public MechanicalPeriodicBC
{
protected:
    virtual void generateConstraints(NodeContainer *nodes, ConstraintContainer *constrs);
public:
    MechanicalSphericalPeriodicBC() { name = "MechanicalSphericalPeriodicBC"; nonsymmetric_shear = false; };
    virtual ~MechanicalSphericalPeriodicBC() {};
    virtual void calculateVolume();
};

//////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////
// Voigt's constraint
class VoigtConstraint : public PBlock
{
    double volume;
    std :: vector< int >strainFunc, stressFunc;
public:
    VoigtConstraint();
    virtual ~VoigtConstraint();
    virtual void apply(NodeContainer *nodes, ElementContainer *elems, BCContainer *bcs, ConstraintContainer *constrs, FunctionContainer *funcs, ExporterContainer *ex, MaterialContainer *mats, RegionContainer *regions, Solver *solv);
    virtual void readFromLine(std :: istringstream &iss, unsigned d);
};

//////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////
class PressureFromMechanicalLoad : public PBlock
{
    double multiplier;
    unsigned master;
    unsigned direction, masterdirection;
    unsigned materialnum;
    std :: vector< unsigned >trsprtnodes;

public:
    PressureFromMechanicalLoad();
    virtual ~PressureFromMechanicalLoad();
    virtual void apply(NodeContainer *nodes, ElementContainer *elems, BCContainer *bcs, ConstraintContainer *constrs, FunctionContainer *funcs, ExporterContainer *ex, MaterialContainer *mats, RegionContainer *regions, Solver *solv);
    virtual void readFromLine(std :: istringstream &iss, unsigned d);
};


//////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////
// Rigid Plate
class RigidPlate : public PBlock
{
private:
    std :: vector< unsigned >slave_ids;
    std :: vector< unsigned >insideRegions, outsideRegions;
public:
    RigidPlate() { };
    virtual ~RigidPlate() {};
    virtual void apply(NodeContainer *nodes, ElementContainer *elems, BCContainer *bcs, ConstraintContainer *constrs, FunctionContainer *funcs, ExporterContainer *ex, MaterialContainer *mats, RegionContainer *regions, Solver *solv);
    virtual void readFromLine(std :: istringstream &iss, unsigned d);
protected:
    bool transport = false;
    unsigned master_id, ndim;
    std :: string which;  ///< which direction to fix (e.g. to leave expansion in perpendicualr direction)
    std :: vector< bool >activeDirs;
    void checkPhysicalField(Node *master);
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
    virtual void apply(NodeContainer *nodes, ElementContainer *elems, BCContainer *bcs, ConstraintContainer *constrs, FunctionContainer *funcs, ExporterContainer *ex, MaterialContainer *mats, RegionContainer *regions, Solver *solv);
    virtual void readFromLine(std :: istringstream &iss, unsigned d);
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
    virtual void apply(NodeContainer *nodes, ElementContainer *elems, BCContainer *bcs, ConstraintContainer *constrs, FunctionContainer *funcs, ExporterContainer *ex, MaterialContainer *mats, RegionContainer *regions, Solver *solv);
    virtual void readFromLine(std :: istringstream &iss, unsigned d);
protected:
    Point center, axis;
    double r_inner, r_outer, w0, w1;
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
    virtual void apply(NodeContainer *nodes, ElementContainer *elems, BCContainer *bcs, ConstraintContainer *constrs, FunctionContainer *funcs, ExporterContainer *ex, MaterialContainer *mats, RegionContainer *regions, Solver *solv);
    virtual void readFromLine(std :: istringstream &iss, unsigned d);
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
    virtual void apply(NodeContainer *nodes, ElementContainer *elems, BCContainer *bcs, ConstraintContainer *constrs, FunctionContainer *funcs, ExporterContainer *ex, MaterialContainer *mats, RegionContainer *regions, Solver *solv);
    virtual void readFromLine(std :: istringstream &iss, unsigned d);
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
    virtual void apply(NodeContainer *nodes, ElementContainer *elems, BCContainer *bcs, ConstraintContainer *constrs, FunctionContainer *funcs, ExporterContainer *ex, MaterialContainer *mats, RegionContainer *regions, Solver *solv);
    // virtual void readFromLine(istringstream &iss, unsigned d);
};

//////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////
// Surface load, pressure
class NormalSurfaceLoad : public PBlock
{
protected:
    unsigned fnID;
public:
    NormalSurfaceLoad() { name = "NormalSurfaceLoad"; };
    virtual ~NormalSurfaceLoad() {};
    virtual void apply(NodeContainer *nodes, ElementContainer *elems, BCContainer *bcs, ConstraintContainer *constrs, FunctionContainer *funcs, ExporterContainer *ex, MaterialContainer *mats, RegionContainer *regions, Solver *solv);
    virtual void readFromLine(std :: istringstream &iss, unsigned d);
};

//////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////
// Surface load, pressure
class MechHangingNode : public PBlock
{
protected:
    unsigned nodeid;
    unsigned elemid;
public:
    MechHangingNode() { name = "MechHangingNode"; };
    virtual ~MechHangingNode() {};
    virtual void apply(NodeContainer *nodes, ElementContainer *elems, BCContainer *bcs, ConstraintContainer *constrs, FunctionContainer *funcs, ExporterContainer *ex, MaterialContainer *mats, RegionContainer *regions, Solver *solv);
    virtual void readFromLine(std :: istringstream &iss, unsigned d);
};

//////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////
// CONTAINER FOR PREPROCESSOR BLOCKS
class PBlockContainer
{
private:
    std :: vector< PBlock * >blocks;
    NodeContainer *nodes;
    ElementContainer *elems;
    BCContainer *bcs;
    ConstraintContainer *constrs;
    FunctionContainer *funcs;
    ExporterContainer *exporters;
    MaterialContainer *materials;
    RegionContainer *regions;
    Solver *solver;
public:
    PBlockContainer() { nodes = nullptr; elems = nullptr; bcs = nullptr; constrs = nullptr; funcs = nullptr; exporters = nullptr; materials = nullptr; regions = nullptr; solver = nullptr; };
    virtual ~PBlockContainer();
    void readFromFile(const std :: string filename, unsigned dim);
    void setContainers(NodeContainer *n, ElementContainer *e, BCContainer *b, ConstraintContainer *c, FunctionContainer *f, ExporterContainer *ex, MaterialContainer *mats, RegionContainer *r, Solver *solver);
    void init();
    void clear();
    unsigned giveSize() const { return blocks.size(); };
    PBlock *givePBlock(unsigned i) { return blocks [ i ]; };
protected:
};


#endif /* _PREPROCESSING_BLOCK_H */
