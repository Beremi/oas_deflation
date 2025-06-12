#include "cross_section.h"
#include <math.h>

using namespace std;


//////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////
// CROSS SECTION OF BEAM
//////////////////////////////////////////////////////////
CrossSection::CrossSection(){
}

/////////////////////////////////////////////////////////
void CrossSection::readFromLine(std :: istringstream &iss, const  unsigned ndim){
    (void) ndim;
    string param;
    bool barea, bIz, bIy, bJ, bkappaY, bkappaZ;
    barea = bIz = bIy = bJ = bkappaY = bkappaZ = false;
    while ( iss >> param ) {
        if ( param.compare("area") == 0 ) {
            iss >> area;
            barea = true;
        } else if ( param.compare("Ix") == 0 ) {
            iss >> Iz;
            bIz = true;
        } else if ( param.compare("Iy") == 0 ) {
            iss >> Iy;
            bIy = true;
        } else if ( param.compare("J") == 0 ) {
            iss >> J;
            bJ = true;
        } else if ( param.compare("kappaY") == 0 ) {
            iss >> kappaY;
            bkappaY = true;
        } else if ( param.compare("kappaZ") == 0 ) {
            iss >> kappaZ;
            bkappaZ = true;
        }
    }
    if (! barea){ 
        cerr << "Error: Cross section parameter 'area' was not specified" << endl;
        exit(1);
    }
    if (! bIz){ 
        cerr << "Error: Cross section parameter 'Iz' was not specified" << endl;
        exit(1);
    }
    if (! bIy){ 
        cerr << "Error: Cross section parameter 'Iy' was not specified" << endl;
        exit(1);
    }
    if (! bJ){ 
        cerr << "Error: Cross section parameter 'J' was not specified" << endl;
        exit(1);
    }
    if (! bkappaY){ 
        cerr << "Error: Cross section parameter 'kappaY' was not specified" << endl;
        exit(1);
    }
    if (! bkappaZ){ 
        cerr << "Error: Cross section parameter 'kappaZ' was not specified" << endl;
        exit(1);
    }

}

/////////////////////////////////////////////////////////
void CrossSection::init(){
}


//////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////
// CIRCULAR CROSS SECTION
//////////////////////////////////////////////////////////
CircularCrossSection::CircularCrossSection(){
}

/////////////////////////////////////////////////////////
void CircularCrossSection::readFromLine(std :: istringstream &iss, const  unsigned ndim){
    (void) ndim;
    string param;
    bool bradius = false;
    bool bnu = false;
    while ( iss >> param ) {
        if ( param.compare("radius") == 0 ) {
            iss >> radius;
            bradius = true;
        } else if ( param.compare("nu") == 0 ) {
            iss >> nu;
            bnu = true;
        }
    }
    if (! bradius){ 
        cerr << "Error: Cross section parameter 'radius' was not specified" << endl;
        exit(1);
    }
    if (! bnu){ 
        cerr << "Error: Cross section parameter 'nu' (Poisson's ratio to compute shear reduction parameter) was not specified" << endl;
        exit(1);
    }
}

/////////////////////////////////////////////////////////
void CircularCrossSection::init(){
    area = M_PI*pow(radius,2);
    Iz = M_PI*pow(radius,4)/4;
    Iy = Iz;
    J = Iy+Iz;
    kappaY = 6*(1+nu)/(7+6*nu); // COWPER, G R. The Shear Coefficient in Timoshenko's Beam Theory
    kappaZ = kappaY;
}
//////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////
// SQUARE CROSS SECTION
//////////////////////////////////////////////////////////
SquareCrossSection::SquareCrossSection(){
}

/////////////////////////////////////////////////////////
void SquareCrossSection::readFromLine(std :: istringstream &iss, const  unsigned ndim){
    (void) ndim;
    string param;
    bool ba = false;
    bool bnu = false;
    while ( iss >> param ) {
        if ( param.compare("a") == 0 ) {
            iss >> a;
            ba = true;
        } else if ( param.compare("nu") == 0 ) {
            iss >> nu;
            bnu = true;
        }
    }
    if (! ba){ 
        cerr << "Error: Cross section parameter 'a' was not specified" << endl;
        exit(1);
    }
    if (! bnu){ 
        cerr << "Error: Cross section parameter 'nu' (Poisson's ratio to compute shear reduction parameter) was not specified" << endl;
        exit(1);
    }
}

/////////////////////////////////////////////////////////
void SquareCrossSection::init(){
    area = a*a;
    Iz = pow(a,4)/12;
    Iy = Iz;
    J = Iz+Iy;
    kappaY = 10*(1+nu)/(12+11*nu); // COWPER, G R. The Shear Coefficient in Timoshenko's Beam Theory
    kappaZ = kappaY;
}

//////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////
// RECTANGULAR CROSS SECTION
//////////////////////////////////////////////////////////
RectangularCrossSection::RectangularCrossSection(){
}

/////////////////////////////////////////////////////////
void RectangularCrossSection::readFromLine(std :: istringstream &iss, const  unsigned ndim){
    (void) ndim;
    string param;
    bool bb = false;
    bool bh = false;
    bool bnu = false;
    while ( iss >> param ) {
        if ( param.compare("b") == 0 ) {
            iss >> b;
            bb = true;
        } else if ( param.compare("h") == 0 ) {
            iss >> h;
            bh = true;
        } else if ( param.compare("nu") == 0 ) {
            iss >> nu;
            bnu = true;
        }
    }
    if (! bb){ 
        cerr << "Error: Cross section parameter 'b' was not specified" << endl;
        exit(1);
    }
    if (! bh){ 
        cerr << "Error: Cross section parameter 'h' was not specified" << endl;
        exit(1);
    }
    if (! bnu){ 
        cerr << "Error: Cross section parameter 'nu' (Poisson's ratio to compute shear reduction parameter) was not specified" << endl;
        exit(1);
    }
}

/////////////////////////////////////////////////////////
void RectangularCrossSection::init(){
    area = b*h;
    Iz = b*pow(h,3)/12;
    Iy = h*pow(b,3)/12;
    J = Iz + Iy;
    kappaY = 10*(1+nu)/(12+11*nu); // COWPER, G R. The Shear Coefficient in Timoshenko's Beam Theory
    kappaZ = kappaY;
}

//////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////
// CROSS SECTION CONTAINER
CrossSectionContainer:: CrossSectionContainer() {
};

/////////////////////////////////////////////////////////
CrossSectionContainer::~CrossSectionContainer(){
    this->clear();
}

/////////////////////////////////////////////////////////
void CrossSectionContainer::readFromFile(const std :: string filename, const unsigned ndim){
    cout << "Input file '" <<  filename;
    size_t origsize = css.size();
    string line, csType;
    ifstream inputfile( filename.c_str() );
    if ( inputfile.is_open() ) {
        while ( getline(inputfile >> std :: ws, line) ) {
            if ( line.empty() || ( line.at(0) == '#' ) ) {
                continue;
            }
            istringstream iss(line);
            iss >> csType;
            if ( !( csType.rfind("#", 0) == 0 ) ) {
                if ( csType.compare("CrossSection") == 0) {
                    CrossSection *newcs = new CrossSection();
                    newcs->readFromLine(iss, ndim);
                    css.push_back(newcs);
                } else {
                    cerr << "Error: cross section '" <<  csType <<  "' does not exists" << endl;
                    exit(EXIT_FAILURE);
                }
            }
        }
        inputfile.close();
        cout << "' succesfully loaded; " << css.size() - origsize << " cross sections found" << endl;
    } else {
        cerr << "Error: unable to open input file '" <<  filename <<  "'" << endl;
        exit(EXIT_FAILURE);
    }
}

/////////////////////////////////////////////////////////
void CrossSectionContainer::init(){
    for ( vector< CrossSection * > :: iterator c = css.begin(); c != css.end(); ++c ) {
        (*c)->init();
    }
}

/////////////////////////////////////////////////////////
void CrossSectionContainer::clear(){
    for ( vector< CrossSection * > :: iterator c = css.begin(); c != css.end(); ++c ) {
        if ( * c != nullptr ) {
            delete * c;
        }
    }
}



