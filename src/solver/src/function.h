#ifndef _FUNCTION_H
#define _FUNCTION_H

#include "linalg.h"
#include <vector>
#include <iostream>
#include <fstream>
#include <algorithm>
#include <typeinfo>

//parser for mathematic expressions
#ifdef __EXPRTK_MODULE
 #include "exprtk.hpp"

typedef exprtk :: symbol_table< double >symbol_table_t;
typedef exprtk :: expression< double >expression_t;
typedef exprtk :: parser< double >parser_t;
#endif

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
    virtual double giveY(double t) const { ( void ) t; return 0; };
    virtual double giveY(const Point *t) { ( void ) t; return 0; }
    virtual void readFromLine(std :: istringstream &iss) = 0;
    virtual double giveNextEtreme(const double &t) const = 0;
    virtual void setActive() { active = true; };
    virtual bool isActive() const { return active; };
protected:
};


#ifdef __EXPRTK_MODULE
//////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////
// GENERAL SPATIAL FUNCTION
class GeneralSpatialFunction : public Function
{
private:
    string expression_string;
    double x, y, z;
    expression_t expression;
public:
    GeneralSpatialFunction() {};
    virtual ~GeneralSpatialFunction() {};
    void readFromLine(istringstream &iss);
    double giveY(const Point *xyz);
    virtual double giveNextEtreme(const double &t) const;
protected:
};
#endif

//////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////
// PIECE-WISE LINEAR FUNCTION
class PieceWiseLinearFunction : public Function
{
protected:
    std :: vector< double >x;
    std :: vector< double >y;
public:
    PieceWiseLinearFunction() {};
    PieceWiseLinearFunction(std :: vector< double >nx, std :: vector< double >ny) { x = nx; y = ny; };
    virtual ~PieceWiseLinearFunction() {};
    void readFromLine(std :: istringstream &iss);
    double giveY(double t) const;
    virtual double giveNextEtreme(const double &t) const;
    void setYValue(double yv, unsigned i) { y [ i ] = yv; };
};

//////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////
// PIECE-WISE LINEAR FUNCTION WITH TIME STEP CONTROLL
class PieceWiseLinearFunctionWithExtremes : public PieceWiseLinearFunction
{
public:
    PieceWiseLinearFunctionWithExtremes() {};
    PieceWiseLinearFunctionWithExtremes(std :: vector< double >nx, std :: vector< double >ny) : PieceWiseLinearFunction(nx, ny) {};
    virtual double giveNextEtreme(const double &t) const;
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
    void readFromLine(std :: istringstream &iss);
    double giveY(double t) const;
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
    void readFromLine(std :: istringstream &iss);
    double giveY(double t) const;
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
    void readFromLine(std :: istringstream &iss);
    double giveY(double t) const;
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
    void readFromLine(std :: istringstream &iss);
    double giveY(double t) const;
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
    void readFromLine(std :: istringstream &iss);
    double giveY(double t) const;
    virtual double giveNextEtreme(const double &t) const;
protected:
};

//////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////
// CONTAINER FOR FUNCTIONS
class FunctionContainer
{
private:
    std :: vector< Function * >functions;
public:
    FunctionContainer() {};
    virtual ~FunctionContainer();
    void readFromFile(const std :: string filename);
    double giveTimeOfNextExtreme(const double &t) const;
    void setActive(const unsigned &fid) { functions [ fid ]->setActive(); };
    bool isActive(const unsigned &fid) const { return functions [ fid ]->isActive(); };
    size_t giveSize() const { return functions.size(); };
    void addFunction(Function *f) { functions.push_back(f); };
    Function *giveFunction(unsigned k);
    void removeFunction(unsigned i);
protected:
};


#endif /* _FUNCTION_H */
