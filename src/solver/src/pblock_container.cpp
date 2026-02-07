#include "pblock_container.h"
#include "model.h"
#include "material_vectorial.h"
#include "element_ldpm.h"
#include "solver_implicit.h"
#include "pblock_constraints.h"
#include "pblock_periodic_bc.h"
#include "pblock_rebar.h"

using namespace std;


//////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////
// CONTAINER FOR PBLOCKS
PBlockContainer :: ~PBlockContainer() {
    for ( vector< PBlock * > :: iterator f = blocks.begin(); f != blocks.end(); ++f ) {
        if ( * f != nullptr ) {
            delete * f;
        }
    }
}

//////////////////////////////////////////////////////////
void PBlockContainer :: clear() {
    for ( vector< PBlock * > :: iterator f = blocks.begin(); f != blocks.end(); ++f ) {
        if ( * f != nullptr ) {
            delete * f;
        }
    }
}

//////////////////////////////////////////////////////////
void PBlockContainer :: setModel(Model *m) {
    model = m;
}

//////////////////////////////////////////////////////////
void PBlockContainer :: init() {
    for ( auto b: blocks ) {
        b->apply(model);
    }
}

//////////////////////////////////////////////////////////
void PBlockContainer :: readFromFile(const string filename, unsigned dim) {
    unsigned origsize = blocks.size(); //todo: warning C4267: 'initializing': conversion from 'size_t' to 'unsigned int', possible loss of dat
    string line, ftype;
    ifstream inputfile( filename.c_str() );
    if ( inputfile.is_open() ) {
        while ( getline(inputfile >> std :: ws, line) ) {
            if ( line.empty() || ( line.at(0) == '#' ) ) {
                continue;
            }
            istringstream iss(line);
            iss >> ftype;
            if ( !( ftype.rfind("#", 0) == 0 ) ) {
                if ( ftype.compare("MechanicalPeriodicBC") == 0 ) {
                    MechanicalPeriodicBC *newblock = new MechanicalPeriodicBC();
                    newblock->readFromLine(iss, dim);
                    blocks.push_back(newblock);
                } else if ( ftype.compare("MechanicalSphericalPeriodicBCExperimental") == 0 ) {
                    MechanicalSphericalPeriodicBCExperimental *newblock = new MechanicalSphericalPeriodicBCExperimental();
                    newblock->readFromLine(iss, dim);
                    blocks.push_back(newblock);
                } else if ( ftype.compare("MechanicalSphericalPeriodicBC") == 0 ) {
                    MechanicalSphericalPeriodicBC *newblock = new MechanicalSphericalPeriodicBC();
                    newblock->readFromLine(iss, dim);
                    blocks.push_back(newblock);
                } else if ( ftype.compare("MechanicalPeriodicBCwithVoigtConstraint") == 0 ) {
                    MechanicalPeriodicBCwithVoigtConstraint *newblock = new MechanicalPeriodicBCwithVoigtConstraint();
                    newblock->readFromLine(iss, dim);
                    blocks.push_back(newblock);
                } else if ( ftype.compare("MechanicalPeriodicBCwithElasticConstraint") == 0 ) {
                    MechanicalPeriodicBCwithElasticConstraint *newblock = new MechanicalPeriodicBCwithElasticConstraint();
                    newblock->readFromLine(iss, dim);
                    blocks.push_back(newblock);
                } else if ( ftype.compare("TransportPeriodicBC") == 0 ) {
                    MechanicalPeriodicBC *newblock = new TransportPeriodicBC();
                    newblock->readFromLine(iss, dim);
                    blocks.push_back(newblock);
                } else if ( ftype.compare("RigidPlate") == 0 ) {
                    RigidPlate *newblock = new RigidPlate();
                    newblock->readFromLine(iss, dim);
                    blocks.push_back(newblock);
                } else if ( ftype.compare("CoordRigidPlate") == 0 ) {
                    CoordRigidPlate *newblock = new CoordRigidPlate();
                    newblock->readFromLine(iss, dim);
                    blocks.push_back(newblock);
                } else if ( ftype.compare("RingRigidPlate") == 0 ) {
                    RingRigidPlate *newblock = new RingRigidPlate();
                    newblock->readFromLine(iss, dim);
                    blocks.push_back(newblock);
                } else if ( ftype.compare("ExpansionRing") == 0 ) {
                    ExpansionRing *newblock = new ExpansionRing();
                    newblock->readFromLine(iss, dim);
                    blocks.push_back(newblock);
                } else if ( ftype.compare("ExpansionRingDoFLoad") == 0 ) {
                    ExpansionRingDoFLoad *newblock = new ExpansionRingDoFLoad();
                    newblock->readFromLine(iss, dim);
                    blocks.push_back(newblock);
                } else if ( ftype.compare("ExpansionRingSingleDoFLoad") == 0 ) {
                    ExpansionRingSingleDoFLoad *newblock = new ExpansionRingSingleDoFLoad();
                    newblock->readFromLine(iss, dim);
                    blocks.push_back(newblock);
                } else if ( ftype.compare("VoigtConstraint") == 0 ) {
                    VoigtConstraint *newblock = new VoigtConstraint();
                    newblock->readFromLine(iss, dim);
                    blocks.push_back(newblock);
                } else if ( ftype.compare("PressureFromMechanicalLoad") == 0 ) {
                    PressureFromMechanicalLoad *newblock = new PressureFromMechanicalLoad();
                    newblock->readFromLine(iss, dim);
                    blocks.push_back(newblock);
                } else if ( ftype.compare("MaterialRegion") == 0 ) {
                    MaterialRegion *newblock = new MaterialRegion();
                    newblock->readFromLine(iss, dim);
                    blocks.push_back(newblock);
                } else if ( ftype.compare("NormalSurfaceLoad") == 0 ) {
                    NormalSurfaceLoad *newblock = new NormalSurfaceLoad();
                    newblock->readFromLine(iss, dim);
                    blocks.push_back(newblock);
                } else if ( ftype.compare("MechHangingNode") == 0 ) {
                    MechHangingNode *newblock = new MechHangingNode();
                    newblock->readFromLine(iss, dim);
                    blocks.push_back(newblock);
                } else if ( ftype.compare("CosseratMechanicalPeriodicBC") == 0 ) {
                    CosseratMechanicalPeriodicBC *newblock = new CosseratMechanicalPeriodicBC();
                    newblock->readFromLine(iss, dim);
                    blocks.push_back(newblock);
                } else if ( ftype.compare("Rebar") == 0 ) {
                    Rebar *newblock = new Rebar();
                    newblock->readFromLine(iss, dim);
                    blocks.push_back(newblock);
                } else {
                    cerr << "Error: preprocessor block '" <<  ftype <<  "' is not implemented yet." << endl;
                    exit(EXIT_FAILURE);
                }
            }
        }
        inputfile.close();
        cout << "Input file '" <<  filename << "' succesfully loaded; " << blocks.size() - origsize << " preprocessing blocks found" << endl;
    } else {
        cerr << "Error: unable to open input file '" <<  filename <<  "'" << endl;
        exit(EXIT_FAILURE);
    }
}
