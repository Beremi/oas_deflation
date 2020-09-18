#ifndef _NODE_H
#define _NODE_H

#include "linear_algebra.h"
#include <vector>
#include <iostream>
#include <fstream>


class BoundaryCondition; //forward declaration, needed to have pointer here, but defined in boundary_condition.h


/**
 * \brief A more elaborate Node class description.
 *
 *  A more elaborate Node class description. DETAILED
 */
class Node
{
private:

protected:
    Point point;  /// center of voronoi cell
    unsigned dim; /// number of dimensions
    unsigned nDoFs; //
    unsigned firstDoF;
    string name;
    bool isMechanical;
    bool isTransport;
public:
    Node() {name = "basic node"; isMechanical = false; isTransport = false; };
    virtual ~Node() {};

    /// A pure virtual member.
    /**
     * \sa testMe()
     * \param c1 the first argument.
     * \param c2 the second argument.
     */
    virtual void readFromLine(istringstream &iss);
    Point givePoint() const { return point; };
    Point *givePointPointer() { return & point; };
    void setPoint(const Point &P) { point = P; };
    void setNumberOfDoFs(const int num) { nDoFs = num; };
    unsigned giveNumberOfDoFs() const { return nDoFs; };
    void setStartingDoF(unsigned num) { firstDoF = num; };
    unsigned giveStartingDoF() const { return firstDoF; };
    string giveName() const { return name; }
    bool doesMechanics() const { return isMechanical; };
    bool doesTransport() const { return isTransport; };
    virtual double giveDoFBasedValue(string code, const Vector &DoFs) const { ( void ) code; ( void ) DoFs; return 0; };
    bool isDoFMechanical(unsigned k) { ( void ) k; return isMechanical; }; //in future we might have node with both fields
    bool isDoFTransport(unsigned k) { ( void ) k; return isTransport; };   //in future we might have node with both fields
};

//////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////
// Auxiliary NODE - only geometrical, without DoFs
class AuxNode : public Node
{
private:

protected:
public:
    AuxNode(unsigned dimension) { dim = dimension; nDoFs = 0; name = "auxiliary node"; };
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
    MechNode(unsigned dimension) { dim = dimension; nDoFs = dim; name = "mechanical node"; isMechanical = true; };
    virtual ~MechNode() {};
    virtual double giveDoFBasedValue(string code, const Vector &DoFs) const;
};

//////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////
// Master Dof - additional DoF governing other dependent DoFs
class MechDoF : public MechNode
{
private:

protected:
public:
    MechDoF(unsigned dimension);
    virtual ~MechDoF() {};
    void readFromLine(istringstream &iss);
};

//////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////
// TRANSPORT NODE - transport DoF
class TrsNode : public Node
{
private:

protected:
public:
    TrsNode(unsigned dimension) { dim = dimension; nDoFs = 1; name = "transport node"; isTransport = true; };
    virtual ~TrsNode() {};
    virtual double giveDoFBasedValue(string code, const Vector &DoFs) const;
};



//////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////
// Master Dof - additional DoF governing other dependent DoFs
class TrsDoF : public TrsNode
{
private:
protected:
public:
    TrsDoF(unsigned dimension);
    virtual ~TrsDoF() {};
    void readFromLine(istringstream &iss);
};


//////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////
// MECHANICAL NODE - translational and rotational DoFs
class Particle : public MechNode
{
private:
    double r;  // radius in case of power tessellation
protected:
public:
    Particle(unsigned dimension) : MechNode(dimension) { nDoFs = 3 * ( dim - 1 ); name = "particle"; };
    virtual ~Particle() {};

    virtual void readFromLine(istringstream &iss);
    double giveRadius() const { return r; };
    void setRadius(const double num) { r = num; };
    virtual double giveDoFBasedValue(string code, const Vector &DoFs) const;
};

#endif /* _NODE_H */
