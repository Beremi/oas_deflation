#ifndef _NODE_H
#define _NODE_H

#include "linear_algebra.h"
#include <vector>
#include <iostream>
#include <fstream>


class BoundaryCondition; //forward declaration, needed to have pointer here, but defined in boundary_condition.h


//////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////
// BASIC NODE - ONLY MASTER CLASS
class Node {
private:

protected:
    Point point;  // center of voronoi cell
    int nDoFs; //
    BoundaryCondition *bc;
    unsigned firstDoF;
    string name;
    bool isMechanical;
    bool isTransport;
public:
    Node(){bc = NULL; name = "basic node"; isMechanical = false; isTransport = false;};
    virtual ~Node(){};

    virtual void readFromLine(istringstream &iss, int dim);
    Point givePoint() const {return point;};
    void setPoint(const Point &P){point = P;};
    void setNumberOfDoFs(const int num){nDoFs = num;};
    int giveNumberOfDoFs() const {return nDoFs;};
    int giveNumberOfFreeDoFs() const;
    void setBC(BoundaryCondition *nbc){ bc = nbc;}
    void setStartingDoF(unsigned num) {firstDoF = num;};
    unsigned giveStartingDoF() const {return firstDoF;};
    string giveName()const{return name;}
    bool doesMechanics()const{return isMechanical;};
    bool doesTransport()const{return isTransport;};
    virtual double giveDoFBasedValue(string code, const Vector &DoFs) const {return 0;};
};

//////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////
// Auxiliary NODE - only geometrical, without DoFs
class AuxNode: public Node{
private:

protected:
public:
    AuxNode(int dim){nDoFs = 0; name = "auxiliary node";};
    virtual ~AuxNode(){};
};

//////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////
// MECHANICAL NODE - translational DoFs
class MechNode: public Node{
private:
protected:
    MechNode(){isMechanical = true;};
public:
    MechNode(int dim){MechNode(); nDoFs = dim; name = "mechanical node";};
    virtual ~MechNode(){};
    virtual double giveDoFBasedValue(string code, const Vector &DoFs) const;
};

//////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////
// TRANSPORT NODE - transport DoF
class TrsNode: public Node{
private:

protected:
public:
    TrsNode(int dim){nDoFs = 1; name = "transport node"; isTransport=true;};
    virtual ~TrsNode(){};
    virtual double giveDoFBasedValue(string code, const Vector &DoFs) const;
};

//////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////
// MECHANICAL NODE - translational and rotational DoFs
class Particle: public MechNode{
private:
    double r;  // radius in case of power tessellation
protected:
public:
    Particle(int dim){nDoFs = 3*(dim-1); name = "particle";};
    virtual ~Particle(){};    

    virtual void readFromLine(istringstream &iss, int dim);
    double giveRadius() const {return r;};
    void setRadius(const double num) {r = num;};
    virtual double giveDoFBasedValue(string code, const Vector &DoFs) const;
};

#endif /* _NODE_H */
