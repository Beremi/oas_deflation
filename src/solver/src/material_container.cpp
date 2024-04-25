#include "material_container.h"
#include "material_tensorial.h"
#include "material_vectorial.h"
#include "material_rve.h"
#include "material_csl.h"
#include "material_ldpm.h"
#include "material_fatigue.h"
#include "material_misc.h"
#include "material_slide_3_2.h"
#include "material_htc.h"
#include "material_coulomb_friction.h"
#include "material_fiber.h"
#include "material_thermomechanical.h"
#include "material_plasticity.h"

#include "element_container.h"

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
        ( * m )->init(this);
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
void MaterialContainer :: readFromFile(const string filename, unsigned dim) {
    cout << "Input file '" <<  filename;
    size_t origsize = matrs.size();
    string line, matType;
    ifstream inputfile(filename.c_str() );
    unsigned id = 0;
    CSLMaterialWithTensorialStressUpdate* CSLMaterialWithTensorialStressUpdateMaster = nullptr;       
    if ( inputfile.is_open() ) {
        while ( getline(inputfile >> std :: ws, line) ) {
            if ( line.empty() || ( line.at(0) == '#' ) ) {
                continue;
            }
            istringstream iss(line);
            iss >> matType;     
            if ( !( matType.rfind("#", 0) == 0 ) ) {
                if ( matType.compare("VectMechMaterial") == 0 ) {
                    VectMechMaterial *newmat = new VectMechMaterial(dim);
                    newmat->readFromLine(iss);
                    matrs.push_back(newmat);
                } else if ( matType.compare("TensMechMaterial") == 0 ) {
                    TensMechMaterial *newmat = new TensMechMaterial(dim);
                    newmat->readFromLine(iss);
                    matrs.push_back(newmat);
                } else if ( matType.compare("TensCosseratMechMaterial") == 0 ) {
                    TensCosseratMechMaterial *newmat = new TensCosseratMechMaterial(dim);
                    newmat->readFromLine(iss);
                    matrs.push_back(newmat);
                } else if ( matType.compare("TensTrsprtMaterial") == 0 ) {
                    TensTrsprtMaterial *newmat = new TensTrsprtMaterial(dim);
                    newmat->readFromLine(iss);
                    matrs.push_back(newmat);
                } else if ( matType.compare("TensHeatConductionMaterial") == 0 ) {
                    TensHeatConductionMaterial *newmat = new TensHeatConductionMaterial(dim);
                    newmat->readFromLine(iss);
                    matrs.push_back(newmat);
                } else if ( matType.compare("VectHeatConductionMaterial") == 0 ) {
                    VectHeatConductionMaterial *newmat = new VectHeatConductionMaterial(dim);
                    newmat->readFromLine(iss);
                    matrs.push_back(newmat);
                } else if ( matType.compare("VectTrsprtMaterial") == 0 ) {
                    VectTrsprtMaterial *newmat = new VectTrsprtMaterial(dim);
                    newmat->readFromLine(iss);
                    matrs.push_back(newmat);
                } else if ( matType.compare("DiscreteTransportRVEMaterial") == 0 ) {
                    DiscreteTransportRVEMaterial *newmat = new DiscreteTransportRVEMaterial(dim);
                    newmat->readFromLine(iss);
                    matrs.push_back(newmat);
                } else if ( matType.compare("DiscreteMechanicalRVEMaterial") == 0 ) {
                    DiscreteMechanicalRVEMaterial *newmat = new DiscreteMechanicalRVEMaterial(dim);
                    newmat->readFromLine(iss);
                    matrs.push_back(newmat);
                } else if ( matType.compare("HTCMaterial") == 0 ) {
                    HTCMaterial *newmat = new HTCMaterial(dim);
                    newmat->readFromLine(iss);
                    matrs.push_back(newmat);
                } else if ( matType.compare("DiscreteCoupledRVEMaterial") == 0 ) {
                    DiscreteCoupledRVEMaterial *newmat = new DiscreteCoupledRVEMaterial(dim);
                    newmat->readFromLine(iss);
                    matrs.push_back(newmat);
                } else if ( matType.compare("VectTrsprtCoupledMaterial") == 0 ) {
                    VectTrsprtCoupledMaterial *newmat = new VectTrsprtCoupledMaterial(dim);
                    newmat->readFromLine(iss);
                    matrs.push_back(newmat);
                } else if ( matType.compare("CSLMaterial") == 0 ) {
                    CSLMaterial *newmat = new CSLMaterial(dim);
                    newmat->readFromLine(iss);
                    matrs.push_back(newmat);
                } else if ( matType.compare("CSLMaterialWithTensorialStressUpdate") == 0 ) {
                    CSLMaterialWithTensorialStressUpdate *newmat = new CSLMaterialWithTensorialStressUpdate(dim,CSLMaterialWithTensorialStressUpdateMaster);
                    newmat->readFromLine(iss);
                    matrs.push_back(newmat);
                    if (CSLMaterialWithTensorialStressUpdateMaster==nullptr) CSLMaterialWithTensorialStressUpdateMaster = newmat;
                } else if ( matType.compare("LDPMMaterial") == 0 ) {
                    LDPMMaterial *newmat = new LDPMMaterial(dim);
                    newmat->readFromLine(iss);
                    matrs.push_back(newmat);
                } else if ( matType.compare("CoupledCSLMaterial") == 0 ) {
                    CoupledCSLMaterial *newmat = new CoupledCSLMaterial(dim);
                    newmat->readFromLine(iss);
                    matrs.push_back(newmat);
                } else if ( matType.compare("ThermoMechanicalMaterial") == 0 ) {
                    ThermoMechanicalMaterial *newmat = new ThermoMechanicalMaterial(dim);
                    newmat->readFromLine(iss);
                    matrs.push_back(newmat);
                } else if ( matType.compare("FatigueShearMaterial") == 0 ) {
                    FatigueShearMaterial *newmat = new FatigueShearMaterial(dim);
                    newmat->readFromLine(iss);
                    matrs.push_back(newmat);
                } else if ( matType.compare("DamagePlasticMaterial") == 0 ) {
                    DamagePlasticMaterial *newmat = new DamagePlasticMaterial(dim);
                    newmat->readFromLine(iss);
                    matrs.push_back(newmat);
                } else if ( matType.compare("FatigueMaterial") == 0 ) {
                    FatigueMaterial *newmat = new FatigueMaterial(dim);
                    newmat->readFromLine(iss);
                    matrs.push_back( ( FatigueShearMaterial * ) newmat );
                } else if ( matType.compare("Slide32Material") == 0 ) {
                    Slide32Material *newmat = new Slide32Material(dim);
                    newmat->readFromLine(iss);
                    matrs.push_back(newmat);
                } else if ( matType.compare("AllicheMaterial") == 0 ) {
                    AllicheMaterial *newmat = new AllicheMaterial(dim);
                    newmat->readFromLine(iss);
                    matrs.push_back(newmat);
                } else if ( matType.compare("DesmoratMaterial") == 0 ) {
                    DesmoratMaterial *newmat = new DesmoratMaterial(dim);
                    newmat->readFromLine(iss);
                    matrs.push_back(newmat);
                } else if ( matType.compare("BrittleMaterial") == 0 ) {
                    BrittleMaterial *newmat = new BrittleMaterial(dim);
                    newmat->readFromLine(iss);
                    matrs.push_back(newmat);
                } else if ( matType.compare("ContactMaterial") == 0 ) {
                    ContactMaterial *newmat = new ContactMaterial(dim);
                    newmat->readFromLine(iss);
                    matrs.push_back(newmat);
                } else if ( matType.compare("CoulombFrictionMaterial") == 0 ) {
                    CoulombFrictionMaterial *newmat = new CoulombFrictionMaterial(dim);
                    newmat->readFromLine(iss);
                    matrs.push_back(newmat);
                } else if ( matType.compare("FiberMaterial") == 0 ) {
                    FiberMaterial *newmat = new FiberMaterial(dim);
                    newmat->readFromLine(iss);
                    matrs.push_back(newmat);
                } else if ( matType.compare("VectMechMaterialWithRotationalStiffness") == 0 ) {
                    VectMechMaterialWithRotationalStiffness *newmat = new VectMechMaterialWithRotationalStiffness(dim);
                    newmat->readFromLine(iss);
                    matrs.push_back(newmat);
                } else if ( matType.compare("VonMisesPlasticMaterial") == 0 ) {
                    VonMisesPlasticMaterial *newmat = new VonMisesPlasticMaterial(dim);
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
        cout << "' succesfully loaded; " << matrs.size() - origsize << " materials found" << endl;
    } else {
        cerr << "Error: unable to open input file '" <<  filename <<  "'" << endl;
        exit(EXIT_FAILURE);
    }
}

//////////////////////////////////////////////////////////
void MaterialContainer :: runPreparationForStressEvaluation(ElementContainer *elems) {
    for ( auto &m : matrs ) {
        m->prepareForStressEvaluation(elems);
    }
}
