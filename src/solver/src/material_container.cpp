#include "material_container.h"

#include "material_rve.h"
#include "material_csl.h"
#include "material_ldpm.h"
#include "material_fatigue.h"
#include "material_misc.h"
#include "material_slide_3_2.h"
#include "material_htc.h"
#include "material_coulomb_friction.h"
#include "material_fiber.h"

using namespace std;

//////////////////////////////////////////////////////////
MaterialContainer :: ~MaterialContainer() {
    for ( vector< Material * > :: iterator m = matrs.begin(); m != matrs.end(); ++m ) {
        if ( * m != nullptr ) {
            delete * m;
        }
    }
}

//////////////////////////////////////////////////////////
void MaterialContainer :: init() {
    for ( vector< Material * > :: iterator m = matrs.begin(); m != matrs.end(); ++m ) {
        ( * m )->init();
    }
}

//////////////////////////////////////////////////////////
Material *MaterialContainer :: giveMaterial(unsigned const mat) {
    if ( mat >= matrs.size() ) {
        cerr << "MaterialContainer Error: material " << mat << " requested, but only " << matrs.size() << " materials exist" << endl;
        exit(1);
    }
    return matrs [ mat ];
}

//////////////////////////////////////////////////////////
void MaterialContainer :: readFromFile(const string filename) {
    size_t origsize = matrs.size();
    string line, matType;
    ifstream inputfile(filename.c_str() );
    unsigned id = 0;
    if ( inputfile.is_open() ) {
        while ( getline(inputfile >> std :: ws, line) ) {
            if ( line.empty() || ( line.at(0) == '#' ) ) {
                continue;
            }
            istringstream iss(line);
            iss >> matType;
            if ( !( matType.rfind("#", 0) == 0 ) ) {
                if ( matType.compare("DisMechMaterial") == 0 ) {
                    DisMechMaterial *newmat = new DisMechMaterial();
                    newmat->readFromLine(iss);
                    matrs.push_back(newmat);
                } else if ( matType.compare("ElasticMechMaterial") == 0 ) {
                    ElasticMechMaterial *newmat = new ElasticMechMaterial();
                    newmat->readFromLine(iss);
                    matrs.push_back(newmat);
                } else if ( matType.compare("CosseratMechMaterial") == 0 ) {
                    CosseratMechMaterial *newmat = new CosseratMechMaterial();
                    newmat->readFromLine(iss);
                    matrs.push_back(newmat);
                } else if ( matType.compare("TrsprtMaterial") == 0 ) {
                    TrsprtMaterial *newmat = new TrsprtMaterial();
                    newmat->readFromLine(iss);
                    matrs.push_back(newmat);
                } else if ( matType.compare("DiscreteTrsprtMaterial") == 0 ) {
                    DiscreteTrsprtMaterial *newmat = new DiscreteTrsprtMaterial();
                    newmat->readFromLine(iss);
                    matrs.push_back(newmat);
                } else if ( matType.compare("DiscreteTransportRVEMaterial") == 0 ) {
                    DiscreteTransportRVEMaterial *newmat = new DiscreteTransportRVEMaterial();
                    newmat->readFromLine(iss);
                    matrs.push_back(newmat);
                } else if ( matType.compare("DiscreteMechanicalRVEMaterial") == 0 ) {
                    DiscreteMechanicalRVEMaterial *newmat = new DiscreteMechanicalRVEMaterial();
                    newmat->readFromLine(iss);
                    matrs.push_back(newmat);
                } else if ( matType.compare("HTCMaterial") == 0 ) {
                    HTCMaterial *newmat = new HTCMaterial();
                    newmat->readFromLine(iss);
                    matrs.push_back(newmat);
                } else if ( matType.compare("DiscreteCoupledRVEMaterial") == 0 ) {
                    DiscreteCoupledRVEMaterial *newmat = new DiscreteCoupledRVEMaterial();
                    newmat->readFromLine(iss);
                    matrs.push_back(newmat);
                } else if ( matType.compare("DiscreteTrsprtCoupledMaterial") == 0 ) {
                    DiscreteTrsprtCoupledMaterial *newmat = new DiscreteTrsprtCoupledMaterial();
                    newmat->readFromLine(iss);
                    matrs.push_back(newmat);
                } else if ( matType.compare("CSLMaterial") == 0 ) {
                    CSLMaterial *newmat = new CSLMaterial();
                    newmat->readFromLine(iss);
                    matrs.push_back(newmat);
                } else if ( matType.compare("LDPMMaterial") == 0 ) {
                    LDPMMaterial *newmat = new LDPMMaterial();
                    newmat->readFromLine(iss);
                    matrs.push_back(newmat);
                } else if ( matType.compare("CoupledCSLMaterial") == 0 ) {
                    CoupledCSLMaterial *newmat = new CoupledCSLMaterial();
                    newmat->readFromLine(iss);
                    matrs.push_back(newmat);
                } else if ( matType.compare("FatigueShearMaterial") == 0 ) {
                    FatigueShearMaterial *newmat = new FatigueShearMaterial();
                    newmat->readFromLine(iss);
                    matrs.push_back(newmat);
                } else if ( matType.compare("DamagePlasticMaterial") == 0 ) {
                    DamagePlasticMaterial *newmat = new DamagePlasticMaterial();
                    newmat->readFromLine(iss);
                    matrs.push_back(newmat);
                } else if ( matType.compare("FatigueMaterial") == 0 ) {
                    FatigueMaterial *newmat = new FatigueMaterial();
                    newmat->readFromLine(iss);
                    matrs.push_back( ( FatigueShearMaterial * ) newmat );
                } else if ( matType.compare("Slide32Material") == 0 ) {
                    Slide32Material *newmat = new Slide32Material();
                    newmat->readFromLine(iss);
                    matrs.push_back(newmat);
                } else if ( matType.compare("AllicheMaterial") == 0 ) {
                    AllicheMaterial *newmat = new AllicheMaterial();
                    newmat->readFromLine(iss);
                    matrs.push_back(newmat);
                } else if ( matType.compare("DesmoratMaterial") == 0 ) {
                    DesmoratMaterial *newmat = new DesmoratMaterial();
                    newmat->readFromLine(iss);
                    matrs.push_back(newmat);
                } else if ( matType.compare("BrittleMaterial") == 0 ) {
                    BrittleMaterial *newmat = new BrittleMaterial();
                    newmat->readFromLine(iss);
                    matrs.push_back(newmat);
                } else if ( matType.compare("ContactMaterial") == 0 ) {
                    ContactMaterial *newmat = new ContactMaterial();
                    newmat->readFromLine(iss);
                    matrs.push_back(newmat);
                } else if ( matType.compare("CoulombFrictionMaterial") == 0 ) {
                    CoulombFrictionMaterial *newmat = new CoulombFrictionMaterial();
                    newmat->readFromLine(iss);
                    matrs.push_back(newmat);
                } else if ( matType.compare("FiberMaterial") == 0 ) {
                    FiberMaterial *newmat = new FiberMaterial();
                    newmat->readFromLine(iss);
                    matrs.push_back(newmat);
                } else {
                    cerr << "Error: material '" <<  matType <<  "' does not exists" << endl;
                    exit(EXIT_FAILURE);
                }
            }
            matrs.back()->setId(id);
            id++;
        }
        inputfile.close();
        cout << "Input file '" <<  filename << "' succesfully loaded; " << matrs.size() - origsize << " materials found" << endl;
    } else {
        cerr << "Error: unable to open input file '" <<  filename <<  "'" << endl;
        exit(EXIT_FAILURE);
    }
}
