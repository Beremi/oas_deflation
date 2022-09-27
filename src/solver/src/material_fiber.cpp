#include "material_fiber.h"

using namespace std;

//////////////////////////////////////////////////////////
// FIBER MATERIAL STATUS

FiberMaterialStatus :: FiberMaterialStatus(FiberMaterial *m, Element *e, unsigned ipnum) : MaterialStatus(m, e, ipnum) {
    name = "fiber mat. status";
    crack_opening = 0;
}

//////////////////////////////////////////////////////////
void FiberMaterialStatus :: giveValues(string code, Vector &result) const {
    MaterialStatus :: giveValues(code, result);
}

//////////////////////////////////////////////////////////
void FiberMaterialStatus :: init() {}

//////////////////////////////////////////////////////////
void FiberMaterialStatus :: update() {
    MaterialStatus :: update();
}

//////////////////////////////////////////////////////////
void FiberMaterialStatus :: resetTemporaryVariables() {
    MaterialStatus :: resetTemporaryVariables();
}



//////////////////////////////////////////////////////////
Matrix FiberMaterialStatus :: giveStiffnessTensor(string type, unsigned dim) const {
    //TODO
    ( void ) type;
    ( void ) dim;
    cout << "FiberMaterialStatus::giveStiffnessTensor" << endl;
    cout.flush();

    return Matrix :: Zero(dim, dim);
}

//////////////////////////////////////////////////////////
Vector FiberMaterialStatus :: giveStress(const Vector &strain, double timeStep) {
    //TODO
    ( void ) strain;
    ( void ) timeStep;
    cout << "FiberMaterialStatus::giveStress" << endl;
    cout.flush();
    return Vector :: Zero( strain.size() );
}

//////////////////////////////////////////////////////////
Vector FiberMaterialStatus :: giveStressWithFrozenIntVars(const Vector &strain, double timeStep) {
    //TODO
    ( void ) strain;
    ( void ) timeStep;
    cout << "FiberMaterialStatus::giveStressWithFrozenIntVars" << endl;
    cout.flush();
    return Vector :: Zero( strain.size() );
}

//////////////////////////////////////////////////////////
void FiberMaterialStatus :: setParameterValue(string code, double value) {
    if ( code.compare("crack_opening") == 0 ) {
        crack_opening = value;
    } else {
        MaterialStatus :: setParameterValue(code, value);
    }
}


//////////////////////////////////////////////////////////
// FIBER MATERIAL

//////////////////////////////////////////////////////////
void FiberMaterial :: readFromLine(istringstream &iss) {
    Material :: readFromLine(iss); //read elastic parameters

    iss.clear(); // clear string stream
    iss.seekg(0, iss.beg); //reset position in string stream

    // initialize all values to zero (NOTE probably no ned in linux, but in windows necessary)
    Ef = Gd = tau0 = betaf = ft = 0;

    string param;
    bool bft, bGd, btau0, bbetaf, bEf;
    bft = bGd = btau0 = bbetaf = bEf =  false;

    while ( !iss.eof() ) {
        iss >> param;
        if ( param.compare("Ef") == 0 ) {
            bEf = true;
            iss >> Ef;
        } else if ( param.compare("Gd") == 0 ) {
            bGd = true;
            iss >> Gd;
        } else if ( param.compare("tau0") == 0 ) {
            btau0 = true;
            iss >> tau0;
        } else if ( param.compare("ft") == 0 ) {
            bft = true;
            iss >> ft;
        } else if ( param.compare("betaf") == 0 ) {
            bbetaf = true;
            iss >> betaf;
        }
    }
    if ( !bEf ) {
        cerr << name << ": material parameter 'Ef' was not specified" << endl;
        exit(EXIT_FAILURE);
    }
    ;
    if ( !bGd ) {
        cerr << name << ": material parameter 'Gd' was not specified" << endl;
        exit(EXIT_FAILURE);
    }
    ;
    if ( !btau0 ) {
        cerr << name << ": material parameter 'tau0' was not specified" << endl;
        exit(EXIT_FAILURE);
    }
    ;
    if ( !bft ) {
        cerr << name << ": material parameter 'ft' was not specified" << endl;
        exit(EXIT_FAILURE);
    }
    ;
    if ( !bbetaf ) {
        cerr << name << ": material parameter 'betaf' was not specified" << endl;
        exit(EXIT_FAILURE);
    }
    ;
};

//////////////////////////////////////////////////////////
MaterialStatus *FiberMaterial :: giveNewMaterialStatus(Element *e, unsigned ipnum) {
    FiberMaterialStatus *newStatus = new FiberMaterialStatus(this, e, ipnum); //needs to be deleted manually
    return newStatus;
};


//////////////////////////////////////////////////////////
void FiberMaterial :: init() {
    // if variables not specified on the input, use default multipliers
    Ef = Gd = tau0 = betaf = ft = 0;
};
