#include "material_container.h"

//////////////////////////////////////////////////////////
MaterialContainer :: ~MaterialContainer() {
    for ( vector< Material * > :: iterator m = matrs.begin(); m != matrs.end(); ++m ) {
        delete * m;
    }
}

//////////////////////////////////////////////////////////
void MaterialContainer :: init() {
    for ( vector< Material * > :: iterator m = matrs.begin(); m != matrs.end(); ++m ) {
        ( * m )->init();
    }
}

//////////////////////////////////////////////////////////
void MaterialContainer :: readFromFile(const string filename) {
    size_t origsize = matrs.size();
    string line, matType;
    ifstream inputfile(filename.c_str() );
    if ( inputfile.is_open() ) {
        while ( getline(inputfile >> std :: ws, line) ) {
            if ( line.empty() ) {
                continue;
            }
            if ( line.at(0) == '#' ) {
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
                } else if ( matType.compare("TrsprtRVEMaterial") == 0 ) {
                    TrsprtRVEMaterial *newmat = new TrsprtRVEMaterial();
                    newmat->readFromLine(iss);
                    matrs.push_back(newmat);
                } else if ( matType.compare("TrsprtCoupledMaterial") == 0 ) {
                    TrsprtCoupledMaterial *newmat = new TrsprtCoupledMaterial();
                    newmat->readFromLine(iss);
                    matrs.push_back(newmat);
                } else if ( matType.compare("MarsMaterial") == 0 ) {
                    MarsMaterial *newmat = new MarsMaterial();
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
                } else if ( matType.compare("AllicheMaterial") == 0 ) {
                    AllicheMaterial *newmat = new AllicheMaterial();
                    newmat->readFromLine(iss);
                    matrs.push_back( ( AllicheMaterial * ) newmat );
                } else if ( matType.compare("DesmoratMaterial") == 0 ) {
                    DesmoratMaterial *newmat = new DesmoratMaterial();
                    newmat->readFromLine(iss);
                    matrs.push_back( ( DesmoratMaterial * ) newmat );
                } else if ( matType.compare("BrittleMaterial") == 0 ) {
                    BrittleMaterial *newmat = new BrittleMaterial();
                    newmat->readFromLine(iss);
                    matrs.push_back( ( BrittleMaterial * ) newmat );
                } else if ( matType.compare("ContactMaterial") == 0 ) {
                    ContactMaterial *newmat = new ContactMaterial();
                    newmat->readFromLine(iss);
                    matrs.push_back( ( ContactMaterial * ) newmat );
                } else {
                    cerr << "Error: material '" <<  matType <<  "' does not exists" << endl;
                    exit(EXIT_FAILURE);
                }
            }
        }
        inputfile.close();
        cout << "Input file '" <<  filename << "' succesfully loaded; " << matrs.size() - origsize << " materials found" << endl;
    } else {
        cerr << "Error: unable to open input file '" <<  filename <<  "'" << endl;
        exit(EXIT_FAILURE);
    }
}
