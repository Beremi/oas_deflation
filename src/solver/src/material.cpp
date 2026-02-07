#include "material.h"
#include "element.h"
#include "element_discrete.h"
#include "material_container.h"

using namespace std;

//////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////
// BASIC MATERIAL
//////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////
MaterialStatus :: MaterialStatus(Material *m, Element *e, unsigned ipnum) {
    name = "basic mat. status";
    mat = m;
    element = e;
    idx = ipnum;
    totalEnergyDensity = 0;
    strainEnergyDensity = 0;
    dissipEnergyDensity = 0;
    updt_dissip_energy = 0;
}

//////////////////////////////////////////////////////////
MaterialStatus :: ~MaterialStatus() {
    for ( vector< MaterialStatus * > :: iterator n = matStatComponents.begin(); n != matStatComponents.end(); ++n ) {
        if ( * n != nullptr ) {
            delete * n;
        }
    }
}

//////////////////////////////////////////////////////////
Vector MaterialStatus :: addEigenStrain(const Vector &totalStrain) const {
    if ( eigenstrain.size() > 0 ) {
        if ( eigenstrain.size() != totalStrain.size() ) {
            cerr << "Material status error: cannot apply eigenstrain of size " << eigenstrain.size() << " to  total strain of size " << totalStrain.size() << endl;
            exit(1);
        }
        Vector activeStrain = totalStrain - eigenstrain;
        return activeStrain;
    } else {
        return totalStrain;
    }
}

//////////////////////////////////////////////////////////
double MaterialStatus :: giveEnergyDissipationIncrement() const {
    return 0;
}

//////////////////////////////////////////////////////////
void MaterialStatus :: computeEnergyDensities() {
    totalEnergyDensity += ( ( temp_stress + updt_stress ).dot(temp_strain - updt_strain) ) / 2.;
    dissipEnergyDensity = updt_dissip_energy + giveEnergyDissipationIncrement();
    strainEnergyDensity = totalEnergyDensity - dissipEnergyDensity;
    //unsigned ndim = element->giveDimension();
    //totalEnergyDensity *= ndim;     //TODO: fix, this works only for discrete material
    //strainEnergyDensity *= ndim;    //TODO: fix, this works only for discrete material
}


//////////////////////////////////////////////////////////
void MaterialStatus :: update() {
    computeEnergyDensities();
    updt_strain = temp_strain;
    updt_strain_total = temp_strain_total;
    updt_stress = temp_stress;
    updt_dissip_energy = dissipEnergyDensity;
}

//////////////////////////////////////////////////////////
void MaterialStatus :: resetTemporaryVariables() {
    temp_strain = updt_strain;
    temp_stress = updt_stress;
}

//////////////////////////////////////////////////////////
void MaterialStatus :: addToEigenStrain(const Vector &x) {
    eigenstrain += x;
}

//////////////////////////////////////////////////////////
bool MaterialStatus :: isElastic(const bool &now) const {
    ( void ) now;
    if ( this->name != "basic mat. status" ) {
        std :: cout << "using elastic check for base class MaterialStatus, if this is not a desire, you need to implement method \'isElastic\' for " << this->name << '\n';
    }
    return true;
}

//////////////////////////////////////////////////////////
void MaterialStatus :: readFromLine(istringstream &iss) {
    std :: string param;
    unsigned num;
    while (  iss >> param ) {
        if ( param.compare("eigenstrain") == 0 ) {
            iss >> num;
            if ( num != mat->giveStrainSize() ) {
                cerr << "Error: strain size is " << mat->giveStrainSize() << " but eigenstrain of size " << num << " provided" << endl;
                exit(1);
            }
            Vector eigs = Vector :: Zero(num);
            for ( unsigned i = 0; i < num; i++ ) {
                iss >> eigs [ i ];
            }
            addToEigenStrain(eigs);
        }
    }
}

//////////////////////////////////////////////////////////
bool MaterialStatus :: giveValues(std :: string code, Vector &result) const {
    if ( code.compare("stress") == 0 || code.compare("stresses") == 0 ) {
        result.resize(temp_strain.size() );
        for ( unsigned i = 0; i < result.size(); i++ ) {
            result [ i ] = temp_stress [ i ];
        }
        return true;
    } else if ( code.compare("strain") == 0  || code.compare("strains") == 0 ) {
        result.resize(temp_strain.size() );
        for ( unsigned i = 0; i < result.size(); i++ ) {
            result [ i ] = temp_strain [ i ];
        }
        return true;
    } else if ( code.compare("materialID") == 0 || code.compare("materialId") == 0 ) {
        result.resize(1);
        result [ 0 ] = mat->giveId();
        return true;
    } else if ( code.compare("total_energy_density") == 0 ) {
        result.resize(1);
        result [ 0 ] = totalEnergyDensity;
        return true;
    } else if ( code.compare("strain_energy_density") == 0 ) {
        result.resize(1);
        result [ 0 ] = strainEnergyDensity;
        return true;
    } else if ( code.compare("dissipated_energy_density") == 0 ) {
        result.resize(1);
        result [ 0 ] = dissipEnergyDensity;
        return true;
    } else {
        result.resize(0);
        return false;
    }
}

//////////////////////////////////////////////////////////
void MaterialStatus :: initializeStressAndStrainVector(unsigned num) {
    if ( num != giveMaterial()->giveStrainSize() ) {
        cerr << "Error in " << giveName() << ": strain size does not match" << endl;
        exit(1);
    }
    temp_stress = temp_strain = updt_stress = updt_strain = Vector :: Zero(num);
}

//////////////////////////////////////////////////////////
Material :: ~Material() {
    for ( vector< Material * > :: iterator n = matComponents.begin(); n != matComponents.end(); ++n ) {
        if ( * n != nullptr ) {
            delete * n;
        }
    }
}

//////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////
// COUPLED MATERIAL
//////////////////////////////////////////////////////////

CoupledMaterialStatus :: CoupledMaterialStatus(Material *m, Element *e, unsigned ipnum) : MaterialStatus(m, e, ipnum) {
    CoupledMaterial *cm = dynamic_cast< CoupledMaterial * >( m );
    if ( !cm ) {
        cerr << name << " " << "Error: material is not derived from  CoupledMaterialStatus" << endl;
        exit(1);
    }

    name = "generic coupled material status";


    vector< Material * >mats = cm->giveMaterials();
    stats.resize( mats.size() );
    for ( unsigned i = 0; i < mats.size(); i++ ) {
        stats [ i ] = mats [ i ]->giveNewMaterialStatus(e, ipnum);
    }
}

//////////////////////////////////////////////////////////
CoupledMaterialStatus :: ~CoupledMaterialStatus() {
    for ( vector< MaterialStatus * > :: iterator m = stats.begin(); m != stats.end(); ++m ) {
        if ( * m != nullptr ) {
            delete * m;
        }
    }
}


//////////////////////////////////////////////////////////
void CoupledMaterialStatus :: init() {
    for ( auto &s:stats ) {
        s->init();
    }
}

//////////////////////////////////////////////////////////
Matrix CoupledMaterialStatus :: giveStiffnessTensor(string type) const {
    ( void ) type;
    unsigned ss = mat->giveStrainSize();
    Matrix D = Matrix :: Zero(ss, ss);
    unsigned k = 0;
    Matrix mD;
    unsigned mk = 0;
    for ( auto &s:stats ) {
        mD = s->giveStiffnessTensor(type);
        mk = s->giveMaterial()->giveStrainSize();
        for ( unsigned i = 0; i < mk; i++ ) {
            for ( unsigned j = 0; j < mk; j++ ) {
                D(k + i, k + j) += mD(i, j);
            }
        }
        k += mk;
    }
    return D;
};

//////////////////////////////////////////////////////////
void CoupledMaterialStatus :: setParameterValue(std :: string code, double value) {
    MaterialStatus :: setParameterValue(code, value);
    for ( auto &s:stats ) {
        s->setParameterValue(code, value);
    }
}

//////////////////////////////////////////////////////////
void CoupledMaterialStatus :: computeStress(double timeStep) {
    Vector mstress;
    unsigned k = 0;
    unsigned h, i;
    for ( auto &s:stats ) {
        s->computeStress(timeStep);
        mstress = s->giveTempStress();
        h = mstress.size();
        for ( i = 0; i < h; i++ ) {
            temp_stress [ k + i ] = mstress [ i ];
        }
        k += h;
    }
}


//////////////////////////////////////////////////////////
void CoupledMaterialStatus :: computeStressWithFrozenIntVars(double timeStep) {
    Vector mstress;
    unsigned k = 0;
    unsigned h, i;
    for ( auto &s:stats ) {
        s->computeStressWithFrozenIntVars(timeStep);
        mstress = s->giveTempStress();
        h = mstress.size();
        for ( i = 0; i < h; i++ ) {
            temp_stress [ k + i ] = mstress [ i ];
        }
        k += h;
    }
}


//////////////////////////////////////////////////////////
void CoupledMaterialStatus :: setTotalTempStrain(Vector str) {
    unsigned k = 0;
    unsigned h, i;
    Vector mstrain;
    for ( auto &s:stats ) {
        h = s->giveMaterial()->giveStrainSize();
        mstrain.resize(h);
        for ( i = 0; i < h; i++ ) {
            mstrain [ i ] = str [ k + i ];
        }
        s->setTotalTempStrain(mstrain);
        k += h;
    }
    temp_strain_total = str;
}

//////////////////////////////////////////////////////////
bool CoupledMaterialStatus :: giveValues(std :: string code, Vector &result) const {
    bool found = 0;
    for ( auto &s:stats ) {
        found = s->giveValues(code, result);
        if ( found ) {
            return true;
        }
    }
    return MaterialStatus :: giveValues(code, result);
}

//////////////////////////////////////////////////////////
void CoupledMaterialStatus :: update() {
    for ( auto &s:stats ) {
        s->update();
    }
    MaterialStatus :: update();
}

//////////////////////////////////////////////////////////
MaterialStatus * CoupledMaterialStatus :: giveMechanicalMaterialStatus() {
    for ( auto &s:stats ) {
        if ( s->giveMechanicalMaterialStatus() ) {
            return s;
        }
    }
    return nullptr;
}

//////////////////////////////////////////////////////////
MaterialStatus * CoupledMaterialStatus :: giveTransportMaterialStatus() {
    for ( auto &s:stats ) {
        if ( s->giveTransportMaterialStatus() ) {
            return s;
        }
    }
    return nullptr;
}

//////////////////////////////////////////////////////////
MaterialStatus * CoupledMaterialStatus :: giveHeatConductionMaterialStatus() {
    for ( auto &s:stats ) {
        if ( s->giveHeatConductionMaterialStatus() ) {
            return s;
        }
    }
    return nullptr;
}

//////////////////////////////////////////////////////////
void CoupledMaterialStatus :: initializeStressAndStrainVector(unsigned num) {
    MaterialStatus :: initializeStressAndStrainVector(num);
    for ( auto &s:stats ) {
        s->initializeStressAndStrainVector( s->giveMaterial()->giveStrainSize() );
    }
}

//////////////////////////////////////////////////////////
void CoupledMaterialStatus :: addToEigenVolumetricStrain(double x) {
    MaterialStatus :: addEigenVolumetricStrain(x);
    for ( auto &s:stats ) {
        if ( s->giveMechanicalMaterialStatus() ) {
            s->addToEigenVolumetricStrain(x);
        }
    }
}

//////////////////////////////////////////////////////////
void CoupledMaterialStatus :: addToEigenStrain(const Vector &x) {
    MaterialStatus :: addToEigenStrain(x);
    unsigned k = 0;
    unsigned h, i;
    Vector mstrain;
    for ( auto &s:stats ) {
        h = s->giveMaterial()->giveStrainSize();
        mstrain.resize(h);
        for ( i = 0; i < h; i++ ) {
            mstrain [ i ] = x [ k + i ];
        }
        s->addToEigenStrain(mstrain);
        k += h;
    }
}

//////////////////////////////////////////////////////////
void CoupledMaterial :: init(MaterialContainer *matcont) {
    mats.resize(nmats);
    for ( unsigned i = 0; i < nmats; i++ ) {
        mats [ i ] = matcont->giveMaterial(matnums [ i ]);
    }

    strainsize = 0;
    for ( auto &m: mats ) {
        strainsize += m->giveStrainSize();
    }

    for ( auto &i : matnums ) {
        if ( i >= idx ) {
            cerr << name << " Error: individual materials inside the coupled material must be specified in advance" << endl;
            exit(1);
        }
    }
}

//////////////////////////////////////////////////////////
CoupledMaterial :: ~CoupledMaterial() {
    //does not apply, the material is deleted by material container automatically
    //for ( vector< Material * > :: iterator m = mats.begin(); m != mats.end(); ++m ) {
    //    if ( * m != nullptr ) {
    //        delete * m;
    //    }
    //}
}

//////////////////////////////////////////////////////////
void CoupledMaterial :: readFromLine(std :: istringstream &iss) {
    iss >> nmats;
    matnums.resize(nmats);
    for ( unsigned i = 0; i < nmats; i++ ) {
        iss >> matnums [ i ];
    }
}

//////////////////////////////////////////////////////////
MaterialStatus *CoupledMaterial :: giveNewMaterialStatus(Element *e, unsigned ipnum) {
    CoupledMaterialStatus *newStatus = new CoupledMaterialStatus(this, e, ipnum);
    return newStatus;
};

//////////////////////////////////////////////////////////
Material *CoupledMaterial :: giveMaterial(unsigned i) const {
    if ( i >= mats.size() ) {
        cerr << name << " Error: requested material number " << i << " but only " << mats.size() << " materials exist." << endl;
        exit(1);
    }
    return mats [ i ];
};

//////////////////////////////////////////////////////////
Material * CoupledMaterial :: giveMechanicalMaterial() {
    for ( auto &s:mats ) {
        if ( s->giveMechanicalMaterial() ) {
            return s;
        }
    }
    return nullptr;
}

//////////////////////////////////////////////////////////
Material * CoupledMaterial :: giveTransportMaterial() {
    for ( auto &s:mats ) {
        if ( s->giveTransportMaterial() ) {
            return s;
        }
    }
    return nullptr;
}

//////////////////////////////////////////////////////////
Material * CoupledMaterial :: giveHeatConductionMaterial() {
    for ( auto &s:mats ) {
        if ( s->giveHeatConductionMaterial() ) {
            return s;
        }
    }
    return nullptr;
}
