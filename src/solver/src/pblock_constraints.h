#ifndef _PREPROCESSING_BLOCK_COSTRAINTS_H
#define _PREPROCESSING_BLOCK_COSTRAINTS_H

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
    virtual void apply(Model *model)  = 0;
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
    virtual void apply(Model *model);
    virtual void readFromLine(std :: istringstream &iss, unsigned d);
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
    virtual void apply(Model *model);
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
    virtual void apply(Model *model);
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
    virtual void apply(Model *model);
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
    virtual void apply(Model *model);
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
    virtual void apply(Model *model);
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
    virtual void apply(Model *model);
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
    virtual void apply(Model *model);
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
    virtual void apply(Model *model);
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
    virtual void apply(Model *model);
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
    virtual void apply(Model *model);
    virtual void readFromLine(std :: istringstream &iss, unsigned d);
};


#endif /* _PREPROCESSING_BLOCK_COSTRAINTS_H */
