#ifndef _PERIODIC_BC_H
#define _PERIODIC_BC_H

#include "pblock_constraints.h"

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

    Vector crack_normal;
    bool crack_active;
    
    int volumetricAverageRigidBC; ///< new boundary condition prescribing average value of pressure

    virtual void generateNewDoFs(NodeContainer *nodes);
    virtual void generateConstraints(NodeContainer *nodes, ConstraintContainer *constrs);
    virtual void generateExporters(NodeContainer *nodes, ExporterContainer *ex);
    virtual void readLoading(std :: istringstream &iss);
    virtual void generateRigidBodyBC(NodeContainer *nodes, ElementContainer *elems, BCContainer *bcs, ConstraintContainer *constrs, FunctionContainer *funcs);
public:
    MechanicalPeriodicBC() { name = "MechanicalPeriodicBC"; nonsymmetric_shear = false; };
    virtual ~MechanicalPeriodicBC() {};
    virtual void apply(Model *model);
    virtual void readFromLine(std :: istringstream &iss, unsigned d);
    std :: vector< double >giveDimensions() const { return PUCsize; };
    double giveVolume() const;
    virtual void calculateVolume(ElementContainer *elems);
    std :: vector< unsigned >giveMasters() const { return masters; };
    std :: vector< unsigned >giveSlaves() const { return slaves; };
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
    void apply(Model *model);
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
// Mechanical Periodic Boundary Condition on sphere
class MechanicalSphericalPeriodicBC : public MechanicalPeriodicBC
{
protected:
    virtual void generateConstraints(NodeContainer *nodes, ConstraintContainer *constrs);
public:
    MechanicalSphericalPeriodicBC() { name = "MechanicalSphericalPeriodicBC"; nonsymmetric_shear = false; };
    virtual ~MechanicalSphericalPeriodicBC() {};
    virtual void calculateVolume(ElementContainer *elems);
};

//////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////
// Mechanical Periodic Boundary Condition on sphere experimental
class MechanicalSphericalPeriodicBCExperimental : public MechanicalPeriodicBC
{
protected:
    virtual void generateConstraints(NodeContainer *nodes, ConstraintContainer *constrs);
    virtual void generateRigidBodyBC(NodeContainer *nodes, ElementContainer *elems, BCContainer *bcs, ConstraintContainer *constrs, FunctionContainer *funcs);
    void constrainRegular(NodeContainer *nodes, ConstraintContainer *constrs, Node *m, Node *s, Point n, Point t);
    void constrainRotation(NodeContainer *nodes, ConstraintContainer *constrs, Node *m, Node *s, Point n, Point t);
public:
    MechanicalSphericalPeriodicBCExperimental() { name = "MechanicalSphericalPeriodicBCExperimental"; nonsymmetric_shear = false; };
    virtual ~MechanicalSphericalPeriodicBCExperimental() {};
    virtual void calculateVolume(ElementContainer *elems);
};

//////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////
// Cosserat Mechanical Periodic Boundary Condition on prism
class CosseratMechanicalPeriodicBC : public MechanicalPeriodicBC
{
protected:
    virtual void generateRigidBodyBC(NodeContainer *nodes, ElementContainer *elems, BCContainer *bcs, ConstraintContainer *constrs, FunctionContainer *funcs);
    virtual void generateNewDoFs(NodeContainer *nodes);
    virtual void generateConstraints(NodeContainer *nodes, ConstraintContainer *constrs);
    virtual void generateExporters(NodeContainer *nodes, ExporterContainer *ex);
    virtual void readLoading(std :: istringstream &iss);
public:
    CosseratMechanicalPeriodicBC() { name = "CosseratMechanicalPeriodicBC"; };
    virtual ~CosseratMechanicalPeriodicBC() {};
};
#endif /* _PERIODIC_BC_H */
