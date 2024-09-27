#ifndef _NODE_H
#define _NODE_H

#include "linalg.h"
#include "simplex.h"
#include <vector>
#include <iostream>
#include <fstream>

class BoundaryCondition; //forward declaration
class RigidBodyContact; //forward declaration
class Transp1D; //forward declaration
class Solver; //forward declaration;

/**
 * \brief A more elaborate Node class description.
 *
 *  A more elaborate Node class description. DETAILED
 */

//physicalFields Mechanics, Transport, Thermal, Humidity

class Node
{
private:

protected:
    Point point;  /// center of voronoi cell
    unsigned dim; /// number of dimensions
    unsigned nDoFs; //
    unsigned firstDoF;
    std :: string name;
    std :: vector< bool >physicalFields;
    std :: vector< unsigned >physicalFieldsDoFNum;
    unsigned id;
    Simplex *simplex; //simplex used to determine volumetric strain
public:
    Node(unsigned dimension);
    virtual ~Node() { if ( simplex != nullptr ) { delete simplex; } };
    virtual void readFromLine(std :: istringstream &iss);
    virtual void init();
    virtual std :: string giveLineToSave() const;
    Point givePoint() const { return point; };
    Point *givePointPointer() { return & point; };
    void setPoint(const Point &P) { point = P; };
    void setNumberOfDoFs(const int num) { nDoFs = num; };
    unsigned giveNumberOfDoFs() const { return nDoFs; };
    void setStartingDoF(unsigned num) { firstDoF = num; };
    unsigned giveStartingDoF() const { return firstDoF; };
    std :: string giveName() const { return name; }
    void setName(const std :: string &newName) { this->name = newName; };
    bool doesMechanics() const { return physicalFields [ 0 ]; };
    bool doesTransport() const { return physicalFields [ 1 ]; };
    bool doesTemperature() const { return physicalFields [ 2 ]; };
    virtual bool giveDoFBasedValues(std :: string code, const Solver *solver, Vector &result) const;
    unsigned giveRelativeDoFPhysicalFieldNum(unsigned k) const;
    unsigned giveAbsoluteDoFPhysicalFieldNum(unsigned k) const;
    std :: vector< unsigned >givePhysicalFieldsDoFNum() const { return physicalFieldsDoFNum; };
    std :: vector< unsigned >givePhysicalFieldNumForAllDoFs() const;
    void subtructFromPoint(Point *p) { point -= ( * p ); };
    virtual unsigned giveOrderOfEnergyConjugateCode(std :: string code) const;
    unsigned giveID() const { return id; };
    void setID(unsigned p) { id = p; };
    Simplex *addElementToSimplex(RigidBodyContact *rbc);
    void initSimplex();
    void updateSimplexVolumetricStrain(const Vector &fullDoFs);
    bool stealSimplexVolumetricStrain();
    Simplex *giveSimplex() { return simplex; }
    bool hasSimplex()const { return simplex != nullptr; }
    bool hasValidSimplex() const { if ( !simplex ) { return false; } return simplex->isValid(); }
    unsigned giveDimension() const { return dim; };
};

//////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////
// Auxiliary NODE - only geometrical, without DoFs
class AuxNode : public Node
{
private:

protected:
public:
    AuxNode(unsigned dimension) : Node(dimension) { name = "AuxNode"; };
    virtual ~AuxNode() {};
};

//////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////
// MECHANICAL NODE - translational DoFs
class MechNode : public Node
{
private:
protected:
public:
    MechNode(unsigned dimension) : Node(dimension) { name = "MechNode"; physicalFields [ 0 ] = true; physicalFieldsDoFNum [ 0 ] = dim; };
    virtual ~MechNode() {};
};

//////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////
// TRANSPORT NODE - transport DoF
class TrsNode : public Node
{
private:

protected:
public:
    TrsNode(unsigned dimension) : Node(dimension) { name = "TransportNode"; physicalFields [ 1 ] = true; physicalFieldsDoFNum [ 1 ] = 1; physicalFields [ 0 ] = false; physicalFieldsDoFNum [ 0 ] = 0;};
    virtual ~TrsNode() {};
};

//////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////
// TEMPERATURE NODE
class TempNode : public Node
{
private:

protected:
public:
    TempNode(unsigned dimension) : Node(dimension) { name = "TempNode"; physicalFields [ 2 ] = true; physicalFieldsDoFNum [ 2 ] = 1;  physicalFields [ 0 ] = false; physicalFieldsDoFNum [ 0 ] = 0;};
    virtual ~TempNode() {};
};


//////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////
// Master Dof - additional DoF governing other dependent DoFs
class FreeDoF : public Node
{
private:

protected:
public:
    FreeDoF(unsigned dimension, unsigned numDoFs);
    virtual ~FreeDoF() {};
    void readFromLine(std :: istringstream &iss);
    virtual void init();
};

//////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////
class MechDoF : public FreeDoF
{
private:

protected:
public:
    MechDoF(unsigned dimension, unsigned numDoFs);
    virtual ~MechDoF() {};
    virtual void init();
};

//////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////
class TrsDoF : public FreeDoF
{
private:
protected:
public:
    TrsDoF(unsigned dimension, unsigned numDoFs);
    virtual ~TrsDoF() {};
    virtual void init();
};

//////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////
class TempDoF : public FreeDoF
{
private:
protected:
public:
    TempDoF(unsigned dimension, unsigned numDoFs);
    virtual ~TempDoF() {};
    virtual void init();
};

//////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////
class HumidityDoF : public FreeDoF
{
private:
protected:
public:
    HumidityDoF(unsigned dimension, unsigned numDoFs);
    virtual ~HumidityDoF() {};
    virtual void init();
};

//////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////
// MECHANICAL NODE - translational and rotational DoFs
class Particle : public MechNode
{
private:
protected:
    double r;  // radius in case of power tessellation
public:
    Particle(unsigned dimension) : MechNode(dimension) { name = "Particle";  physicalFieldsDoFNum [ 0 ] = 3 * ( dim - 1 ); };
    virtual ~Particle() {};
    virtual void readFromLine(std :: istringstream &iss);
    double giveRadius() const { return r; };
    void setRadius(const double num) { r = num; };
    virtual bool giveDoFBasedValues(std :: string code, const Solver *solver, Vector &result) const;
    virtual unsigned giveOrderOfEnergyConjugateCode(std :: string code) const;
    virtual std :: string giveLineToSave() const;
    Vector calculateRigidBodyMotionVector(const Point *x, const Vector &DoFs) const;
    Point calculateRigidBodyMotionPoint(const Point *x, const Vector &DoFs) const;
    Matrix giveRigidBodyMotionMatrix(const Point *x) const;
};

//////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////
// COUPLED PARTICLE - translational and rotational and transport DoFs
class ParticleWithTransport : public Particle
{
private:
protected:
public:
    ParticleWithTransport(unsigned dimension) : Particle(dimension) { name = "ParticleWithTransport";  physicalFields [ 1 ] = true; physicalFieldsDoFNum [ 1 ] = 1; };
    virtual ~ParticleWithTransport() {};
};

//////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////
// PARTICLE COUPLED WITH TEMPARATURE
class ParticleWithTemperature : public Particle
{
private:
protected:
public:
    ParticleWithTemperature(unsigned dimension) : Particle(dimension) { nDoFs = 3 * ( dim - 1 ) + 1; name = "ParticleWithTemperature";  physicalFields [ 2 ] = true; physicalFieldsDoFNum [ 2 ] = 1; };
    virtual ~ParticleWithTemperature() {};
};

//////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////
// TRANSPORT + TEMPERATURE NODE
class TempHumidityNode : public TempNode
{
private:

protected:
public:
    TempHumidityNode(unsigned dimension) : TempNode(dimension) { name = "TempHumidityNode"; physicalFields [ 3 ] = true; physicalFieldsDoFNum [ 3 ] = 1; };
    virtual ~TempHumidityNode() {};
};

#endif /* _NODE_H */
