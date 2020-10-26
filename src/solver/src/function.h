#ifndef _FUNCTION_H
#define _FUNCTION_H

#include "linear_algebra.h"
#include <vector>
#include <iostream>
#include <fstream>
#include <algorithm>
#include <typeinfo>


//////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////
// BASIC FUNCTION - MASTER CLASS
class Function
{
private:
    bool active = false;
public:
    Function() {};
    virtual ~Function() {};
    virtual double giveY(double t)  = 0;
    virtual void readFromLine(istringstream &iss) = 0;
    virtual double giveNextEtreme(const double &t) const = 0;
    virtual void setActive() { active = true; };
    virtual bool isActive() const { return active; };
protected:
};

//////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////
// PIECE-WISE LINEAR FUNCTION
class PieceWiseLinearFunction : public Function
{
private:
    vector< double >x;
    vector< double >y;
public:
    PieceWiseLinearFunction() {};
    PieceWiseLinearFunction(vector< double >nx, vector< double >ny) { x = nx; y = ny; };
    virtual ~PieceWiseLinearFunction() {};
    void readFromLine(istringstream &iss);
    double giveY(double t);
    virtual double giveNextEtreme(const double &t) const;
    void setYValue(double yv, unsigned i) { y [ i ] = yv; };
protected:
};

//////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////
// SAW TOOTH FUNCTION WITH CONSTANT MAX VALUE
class ConstSawToothFunction : public Function
{
private:
    double upper, lower, period, time_shift;
    int multip = 1;
public:
    ConstSawToothFunction() {};
    virtual ~ConstSawToothFunction() {};
    void readFromLine(istringstream &iss);
    double giveY(double t);
    virtual double giveNextEtreme(const double &t) const;
protected:
};


//////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////
// SAW TOOTH FUNCTION WITH linearly increasing MAX VALUE
class LinSawToothFunction : public ConstSawToothFunction
{
private:
    double time_multiplier;
public:
    LinSawToothFunction() {};
    virtual ~LinSawToothFunction() {};
    void readFromLine(istringstream &iss);
    double giveY(double t);
protected:
};

//////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////
// SAW TOOTH FUNCTION WITH VARYING MAX VALUE
// INHERITED FROM TWO CLASSES
// BE CAREFUL FURTHER INHERITING THIS CLASS!!!
class VaryingSawToothFunction : public ConstSawToothFunction, public PieceWiseLinearFunction
{
private:
    double shift = 0;
public:
    VaryingSawToothFunction() {};
    virtual ~VaryingSawToothFunction() {};
    void readFromLine(istringstream &iss);
    double giveY(double t);
protected:
};

//////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////
// JM: Rotation function using angle multiplier from a ConstSawToothFunction
class ConstSawToothRotationFunction : public ConstSawToothFunction
{
private:
    Point initNodePosition;
    Point rotationAngles;
    unsigned int displacementType;
    double currentTime;
    double previousTime;
public:
    ConstSawToothRotationFunction() {};
    virtual ~ConstSawToothRotationFunction() {};
    void readFromLine(istringstream &iss);
    double giveY(double t);
    void setCurrentTime(double t);
protected:
};

//////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////
// JM: Shear function using multiplier from a ConstSawToothFunction
class ConstSawToothShearFunction : public ConstSawToothFunction
{
private:
    Point initNodePosition;
    Point governingPoint;
    unsigned int displacementType;
    double currentTime;
    double previousTime;
public:
    ConstSawToothShearFunction() {};
    virtual ~ConstSawToothShearFunction() {};
    void readFromLine(istringstream &iss);
    double giveY(double t);
    void setCurrentTime(double t);
protected:
};

//////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////
// SINUS FUNCTION
// TODO make parent functin "PeriodicFunction" for both SinusFn and SawToothFn
class SinusFunction : public Function
{
private:
    double period, amplitude, shift;
public:
    SinusFunction() {};
    virtual ~SinusFunction() {};
    void readFromLine(istringstream &iss);
    double giveY(double t);
    virtual double giveNextEtreme(const double &t) const;
protected:
};

//////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////
// CONTAINER FOR FUNCTIONS
class FunctionContainer
{
private:
    vector< Function * >functions;
public:
    FunctionContainer() {};
    virtual ~FunctionContainer();
    void readFromFile(const string filename);
    double giveY(unsigned f, double t);
    double giveTimeOfNextExtreme(const double &t) const;
    void setActive(const unsigned &fid) { functions [ fid ]->setActive(); };
    bool isActive(const unsigned &fid) const { return functions [ fid ]->isActive(); };
    size_t giveSize() const { return functions.size(); };
    void addFunction(Function *f) { functions.push_back(f); };
    Function *giveFunction(unsigned k) { return functions [ k ]; };
protected:
};


#endif /* _FUNCTION_H */
