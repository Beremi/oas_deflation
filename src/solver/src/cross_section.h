#ifndef _CROSS_SECTION_H
#define _CROSS_SECTION_H

#include "globals.h"

//////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////
// CROSS SECTION OF BEAM
class CrossSection
{
protected: 
    double area, Iz, Iy, J, kappaY, kappaZ;
public:
    CrossSection();
    virtual ~CrossSection() { };  
    virtual void readFromLine(std :: istringstream &iss, const  unsigned ndim);
    virtual void init();
    double giveArea() const {return area;};
    double giveIy() const {return Iy;};
    double giveIz() const {return Iz;};
    double giveJ() const {return J;};
    double giveKappaY() const {return kappaY;};  
    double giveKappaZ() const {return kappaZ;};  
};

//////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////
// CIRCULAR CROSS SECTION
class CircularCrossSection : public CrossSection
{
protected: 
    double radius;
    double nu; //Poisson number for shear reduction coeff. kappa
public:
    CircularCrossSection();
    CircularCrossSection(double radius, double nu);
    virtual ~CircularCrossSection() { };  
    virtual void readFromLine(std :: istringstream &iss, const  unsigned ndim);
    virtual void init();
    double giveDiameter()const{return 2*radius;};
};

//////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////
// SQUARE CROSS SECTION
class SquareCrossSection : public CrossSection
{
protected: 
    double a;
    double nu; //Poisson number for shear reduction coeff. kappa
public:
    SquareCrossSection();
    virtual ~SquareCrossSection() { };  
    virtual void readFromLine(std :: istringstream &iss, const  unsigned ndim);
    virtual void init();
};

//////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////
// RECTANGULAR CROSS SECTION
class RectangularCrossSection : public CrossSection
{
protected: 
    double b,h;
    double nu; //Poisson number for shear reduction coeff. kappa
public:
    RectangularCrossSection();
    virtual ~RectangularCrossSection() { };  
    virtual void readFromLine(std :: istringstream &iss, const  unsigned ndim);
    virtual void init();
};

//////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////
// CROSS SECTION CONTAINER
class CrossSectionContainer
{
private:
    std :: vector< CrossSection * >css;
public:
    CrossSectionContainer();
    virtual ~CrossSectionContainer();
    void readFromFile(const std :: string filename, const unsigned ndim);
    void init();
    void clear();
    size_t giveSize() const { return css.size(); }
    CrossSection* giveCrossSection(unsigned i){return css[i];};
    void addCrossSection(CrossSection *cs){css.push_back(cs);};

    std :: vector< CrossSection * > :: iterator begin() { return css.begin(); }
    std :: vector< CrossSection * > :: iterator end() { return css.end(); }
    std :: vector< CrossSection * > :: const_iterator begin() const { return css.begin(); }
    std :: vector< CrossSection * > :: const_iterator end() const { return css.end(); }
};



#endif  /* _CROSS_SECTION_H */
