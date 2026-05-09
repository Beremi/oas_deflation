#include "element_container.h"
#include "element_discrete.h"
#include "element_continuous.h"
#include "element_fiber.h"
#include "element_polyhedral.h"
#include "element_ldpm.h"
#ifdef ML_TORCH_FOUND
 #include "element_superelem.h"
#endif // TORCH_FOUND
#include <algorithm>
#include <chrono>
#include "model.h"
#include "pblock_periodic_bc.h"
#include "constraint.h"
#include "cross_section.h"
#include "element_beam.h"
#include "openmp_utils.h"

using namespace std;

namespace {
void recordElementRuntimePhase(Model *model, const std :: string &phase, double durationSeconds, const std :: string &detail = "", long long sampleCount = 1) {
    if ( model && model->giveSolver() ) {
        model->giveSolver()->recordExternalRuntimePhaseSamples(phase, durationSeconds, sampleCount, detail, "ElementContainer");
    }
}
}

//////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////
// ELEMENT CONATINER
ElementContainer :: ~ElementContainer() {
    for ( vector< Element * > :: iterator e = elems.begin(); e != elems.end(); ++e ) {
        if ( * e != nullptr ) {
            delete * e;
        }
    }
}

//////////////////////////////////////////////////////////
void ElementContainer :: clear() {
    clearStiffnessEntryMap();
    for ( vector< Element * > :: iterator e = elems.begin(); e != elems.end(); ++e ) {
        if ( * e != nullptr ) {
            delete * e;
        }
    }
}


//////////////////////////////////////////////////////////
void ElementContainer :: setModel(Model *mod) {
    model = mod;
    nodes = model->giveNodes();
    bconds = model->giveBoundaryConditions();
    crosssects = model->giveCrossSections();
};

//////////////////////////////////////////////////////////
void ElementContainer :: readFromFile(const string filename, const unsigned ndim, MaterialContainer *matrs) {
    cout << "Input file '" <<  filename;
    this->materials = matrs;
    size_t origsize = elems.size();
    string line, elemType;
    ifstream inputfile( filename.c_str() );
    if ( inputfile.is_open() ) {
        while ( getline(inputfile >> std :: ws, line) ) {
            if ( line.empty() || ( line.at(0) == '#' ) ) {
                continue;
            }
            istringstream iss(line);
            iss >> elemType;
            if ( !( elemType.rfind("#", 0) == 0 ) ) {
                if ( elemType.compare("CSLElem") == 0 || elemType.compare("LTCBEAM") == 0 ) {
                    RigidBodyContact *newelem = new RigidBodyContact(ndim);
                    newelem->readFromLine(iss, nodes, matrs);
                    elems.push_back(newelem);
                } else if ( elemType.compare("CSLCoupledElem") == 0 || elemType.compare("LTCBEAMCoupled") == 0 ) {
                    RigidBodyContactCoupled *newelem = new RigidBodyContactCoupled(ndim);
                    newelem->readFromLine(iss, nodes, matrs);
                    elems.push_back(newelem);
                } else if ( elemType.compare("LDPMTetra") == 0 ) {
                    LDPMTetra *newelem = new LDPMTetra(ndim);
                    newelem->readFromLine(iss, nodes, matrs);
                    elems.push_back(newelem);
                } else if ( elemType.compare("LDPMTetraWithTransportAndHeatConduction") == 0 ) {
                    LDPMTetraWithTransportAndHeatConduction *newelem = new LDPMTetraWithTransportAndHeatConduction(ndim);
                    newelem->readFromLine(iss, nodes, matrs);
                    elems.push_back(newelem);
                } else if ( elemType.compare("LDPMTetraWithTransport") == 0 ) {
                    LDPMTetraWithTransport *newelem = new LDPMTetraWithTransport(ndim);
                    newelem->readFromLine(iss, nodes, matrs);
                    elems.push_back(newelem);
                } else if ( elemType.compare("LDPMEdgeTransport") == 0 ) {
                    LDPMEdgeTransport *newelem = new LDPMEdgeTransport(this);
                    newelem->readFromLine(iss, nodes, matrs);
                    elems.push_back(newelem);
                } else if ( elemType.compare("LDPMEdgeTransportBoundary") == 0 ) {
                    LDPMEdgeTransportBoundary *newelem = new LDPMEdgeTransportBoundary(this);
                    newelem->readFromLine(iss, nodes, matrs);
                    elems.push_back(newelem);
                } else if ( elemType.compare("CLSBoundaryElem") == 0 || elemType.compare("LTCBoundary") == 0 ) {
                    RigidBodyBoundary *newelem = new RigidBodyBoundary(ndim);
                    newelem->readFromLine(iss, nodes, matrs);
                    elems.push_back(newelem);
                } else if ( elemType.compare("CLSBoundaryCoupledElem") == 0 || elemType.compare("LTCBoundaryCoupled") == 0 ) {
                    RigidBodyBoundaryCoupled *newelem = new RigidBodyBoundaryCoupled(ndim);
                    newelem->readFromLine(iss, nodes, matrs);
                    elems.push_back(newelem);
                } else if ( elemType.compare("Truss") == 0 ) {
                    Truss *newelem = new Truss(ndim);
                    newelem->readFromLine(iss, nodes, matrs);
                    elems.push_back(newelem);
                } else if ( elemType.compare("DiscreteTransportElem") == 0 || elemType.compare("LTCTRSP") == 0 ) {
                    DiscreteTrsprtElem *newelem = new DiscreteTrsprtElem(ndim);
                    newelem->readFromLine(iss, nodes, matrs);
                    elems.push_back(newelem);
                } else if ( elemType.compare("DiscreteTrsprtCoupledElem") == 0 || elemType.compare("LTCTRSPCoupled") == 0 ) {
                    DiscreteTrsprtCoupledElem *newelem = new DiscreteTrsprtCoupledElem(ndim);
                    newelem->readFromLine(iss, nodes, matrs);
                    elems.push_back(newelem);
                } else if ( elemType.compare("TrsprtTri") == 0 ) {
                    TrsprtTriangle *newelem = new TrsprtTriangle();
                    newelem->readFromLine(iss, nodes, matrs);
                    elems.push_back(newelem);
                } else if ( elemType.compare("TrsprtQuad") == 0 ) {
                    TrsprtQuad *newelem = new TrsprtQuad();
                    newelem->readFromLine(iss, nodes, matrs);
                    elems.push_back(newelem);
                } else if ( elemType.compare("TrsprtTet") == 0 ) {
                    TrsprtTetra *newelem = new TrsprtTetra();
                    newelem->readFromLine(iss, nodes, matrs);
                    elems.push_back(newelem);
                } else if ( elemType.compare("TrsprtBrick") == 0 ) {
                    TrsprtBrick *newelem = new TrsprtBrick();
                    newelem->readFromLine(iss, nodes, matrs);
                    elems.push_back(newelem);
                } else if ( elemType.compare("TrsprtTemprtrCoupledBrick") == 0 ) {
                    TrsprtTemprtrCoupledBrick *newelem = new TrsprtTemprtrCoupledBrick();
                    newelem->readFromLine(iss, nodes, matrs);
                    elems.push_back(newelem);
                } else if ( elemType.compare("MechanicalTri") == 0 ) {
                    MechanicalTriangle *newelem = new MechanicalTriangle();
                    newelem->readFromLine(iss, nodes, matrs);
                    elems.push_back(newelem);
                } else if ( elemType.compare("MechanicalQuad") == 0 ) {
                    MechanicalQuad *newelem = new MechanicalQuad();
                    newelem->readFromLine(iss, nodes, matrs);
                    elems.push_back(newelem);
                } else if ( elemType.compare("MechanicalTet") == 0 ) {
                    MechanicalTetra *newelem = new MechanicalTetra();
                    newelem->readFromLine(iss, nodes, matrs);
                    elems.push_back(newelem);
                } else if ( elemType.compare("MechanicalBrick") == 0 ) {
                    MechanicalBrick *newelem = new MechanicalBrick();
                    newelem->readFromLine(iss, nodes, matrs);
                    elems.push_back(newelem);
                } else if ( elemType.compare("CosseratQuad") == 0 ) {
                    CosseratQuad *newelem = new CosseratQuad();
                    newelem->readFromLine(iss, nodes, matrs);
                    elems.push_back(newelem);
                } else if ( elemType.compare("CosseratBrick") == 0 ) {
                    CosseratBrick *newelem = new CosseratBrick();
                    newelem->readFromLine(iss, nodes, matrs);
                    elems.push_back(newelem);
                } else if ( elemType.compare("CoupledCosseratTransportQuad") == 0 ) {
                    CoupledCosseratTransportQuad *newelem = new CoupledCosseratTransportQuad();
                    newelem->readFromLine(iss, nodes, matrs);
                    elems.push_back(newelem);
                } else if ( elemType.compare("CoupledCosseratTransportBrick") == 0 ) {
                    CoupledCosseratTransportBrick *newelem = new CoupledCosseratTransportBrick();
                    newelem->readFromLine(iss, nodes, matrs);
                    elems.push_back(newelem);
                } else if ( elemType.compare("Fiber") == 0 ) {
                    Fiber *newelem = new Fiber(ndim);
                    newelem->readFromLine(iss, nodes, matrs);
                    elems.push_back(newelem);
                } else if ( elemType.compare("CoupledCosseratBrickWithDependentUpperZLayer") == 0 ) {
                    CoupledCosseratBrickWithDependentUpperZLayer *newelem = new CoupledCosseratBrickWithDependentUpperZLayer();
                    newelem->readFromLine(iss, nodes, matrs);
                    elems.push_back(newelem);
                } else if ( elemType.compare("TranspPolygonal") == 0 ) {
                    TranspPolygonal *newelem = new TranspPolygonal(ndim);
                    newelem->readFromLine(iss, nodes, matrs);
                    elems.push_back(newelem);
                } else if ( elemType.compare("TranspVirtPolygonal") == 0 ) {
                    TranspVirtPolygonal *newelem = new TranspVirtPolygonal(ndim);
                    newelem->readFromLine(iss, nodes, matrs);
                    elems.push_back(newelem);
                } else if ( elemType.compare("TranspCondensedPolygonal") == 0 ) {
                    TranspCondensedPolygonal *newelem = new TranspCondensedPolygonal(ndim);
                    newelem->readFromLine(iss, nodes, matrs);
                    elems.push_back(newelem);
                } else if ( elemType.compare("DiscreteHeatConductionElem") == 0 ) {
                    DiscreteHeatConductionElem *newelem = new DiscreteHeatConductionElem(ndim);
                    newelem->readFromLine(iss, nodes, matrs);
                    elems.push_back(newelem);
                } else if ( elemType.compare("RigidBodyContactWithHeatConduction") == 0 ) {
                    RigidBodyContactWithHeatConduction *newelem = new RigidBodyContactWithHeatConduction(ndim);
                    newelem->readFromLine(iss, nodes, matrs);
                    elems.push_back(newelem);
                    //} else if ( elemType.compare("EulerBernoulliBeam3D") == 0 ) {
                    //    TimoshenkoBeam3D *newelem = new TimoshenkoBeam3D(ndim);
                    //    newelem->readFromLine(iss, nodes, matrs, crosssects);
                    //    elems.push_back(newelem);
                } else if ( elemType.compare("TimoshenkoBeam3D") == 0 ) {
                    TimoshenkoBeam3D *newelem = new TimoshenkoBeam3D(ndim);
                    newelem->readFromLine(iss, nodes, matrs, crosssects);
                    elems.push_back(newelem);
                } else if ( elemType.compare("MLMechElement") == 0 ) {
#ifdef ML_TORCH_FOUND
                    MLMechElement *newelem = new MLMechElement(ndim);
                    newelem->readFromLine(iss, nodes, matrs);
                    elems.push_back(newelem);
#else
                    cerr << "Error: This OAS executable compiled without MLMechElement support." << endl;
                    exit(EXIT_FAILURE);
#endif                 // TORCH_FOUND
                } else if ( elemType.compare("MaterialTestElement") == 0 ) {
                    MaterialTestElement *newelem = new MaterialTestElement(ndim);
                    newelem->readFromLine(iss, nodes, matrs);
                    elems.push_back(newelem);
                    /*} else if ( elemType.compare("PolyhedralFace") == 0 ) {
                     * PolyhedralFace *newelem = new PolyhedralFace(ndim);
                     * newelem->readFromLine(iss, nodes, matrs);
                     * elems.push_back(newelem);
                     *  } else if ( elemType.compare("TranspPolyhedral") == 0 ) {
                     *  TranspPolyhedral *newelem = new TranspPolyhedral(ndim);
                     *  newelem->readFromLine(iss, nodes, matrs);
                     *  elems.push_back(newelem);
                     *
                     * } else if ( elemType.compare("TranspVirtPolyhedral") == 0 ) {
                     * TranspVirtPolyhedral *newelem = new TranspVirtPolyhedral(ndim);
                     * newelem->readFromLine(iss, nodes, matrs);
                     * elems.push_back(newelem);
                     * } else if ( elemType.compare("TranspCondensedPolyhedral") == 0 ) {
                     * TranspCondensedPolyhedral *newelem = new TranspCondensedPolyhedral(ndim);
                     * newelem->readFromLine(iss, nodes, matrs);
                     * elems.push_back(newelem);*/
                } else if ( elemType.compare("RigidBodyContactWithRotationalStiffness") == 0 ) {
                    RigidBodyContactWithRotationalStiffness *newelem = new RigidBodyContactWithRotationalStiffness(ndim);
                    newelem->readFromLine(iss, nodes, matrs);
                    elems.push_back(newelem);
                } else if ( elemType.compare("RigidBodyContactWithDecoupledRotationsAndTranslations") == 0 ) {
                    RigidBodyContactWithDecoupledRotationsAndTranslations *newelem = new RigidBodyContactWithDecoupledRotationsAndTranslations(ndim);
                    newelem->readFromLine(iss, nodes, matrs);
                    elems.push_back(newelem);
                } else {
                    cerr << "Error: element '" <<  elemType <<  "' does not exists" << endl;
                    exit(EXIT_FAILURE);
                }
            }
        }
        inputfile.close();
        cout << "' succesfully loaded; " << elems.size() - origsize << " elements found" << endl;
    } else {
        cerr << "Error: unable to open input file '" <<  filename <<  "'" << endl;
        exit(EXIT_FAILURE);
    }
}

//////////////////////////////////////////////////////////
Element *ElementContainer :: giveElement(unsigned const num) const {
    if ( num < elems.size() ) {
        return elems [ num ];
    }
    cerr << "ElementContainer Error: requested element no. " << num << " but only " << elems.size() << " exist" << endl;
    exit(1);
}

//////////////////////////////////////////////////////////
void ElementContainer :: addElement(Element *newelem) {
    clearStiffnessEntryMap();
    newelem->init();
    newelem->setID(elems.size() );
    newelem->initMaterialStatuses();
    elems.push_back(newelem);
}


//////////////////////////////////////////////////////////
// void ElementContainer ::  saveToFile(const string &filepath, std :: vector< unsigned > &elems_to_save) const {
//     std :: ofstream outputfile( filepath );
//     if ( outputfile.is_open() ) {
//         outputfile << "#elements saved from calculation";
//         for ( auto const &elem_id : elems_to_save) {
//           outputfile << this->giveElement(elem_id)->giveLineToSave(this->nodes) << '\n';
//         }
//         outputfile.close();
//     }
// }

unsigned ElementContainer :: giveElemId(const Element *elem) const {
    // do not use this method for node that is not a part of this (nodeContainer)
    auto res = std :: find(std :: begin(this->elems), std :: end(this->elems), elem);
    if ( res == this->elems.end() ) {
        // if node is not in container, return zero (but zero can be also for the first node)
        // just to prevent errors here
        return 0;
    }
    return std :: distance(std :: begin(this->elems), res);
}

//////////////////////////////////////////////////////////
void ElementContainer :: saveElemStatsToFile(const string &filepath, const std :: vector< unsigned > &elems_to_save, const double time_now, const unsigned step, const bool saveNodeIds, const double idc_time, const double time_step) const {
    std :: ofstream outputfile(filepath);
    unsigned stat_id = 0;
    unsigned num;
    if ( outputfile.is_open() ) {
        outputfile << "#elem_id stat_id mat_id ";
        if ( saveNodeIds ) {
            outputfile << "num_nodes node ids ... ";
        }
        outputfile << "internal_variables ...";
        if ( time_now != 0 && step != 0 ) {
            outputfile <<  scientific;
            outputfile << "\ntime " << time_now << " " << step << " " << time_step << " " << idc_time;
        }
        // if ( time_step != 0 ) {
        //     outputfile << "\ntime_step " << time_step;
        // }
        for ( auto const &elem_id : elems_to_save ) {
            stat_id = 0;
            for ( auto const &mat_stat : this->giveElement(elem_id)->giveMaterialStats() ) {
                // elem_id - stat_id -  mat_id - internal_variables
                outputfile << "\nmatStat\t" << elem_id << '\t' << stat_id++ << '\t' << mat_stat->giveMaterial()->giveId() << '\t';
                if ( saveNodeIds ) {
                    num = this->giveElement(elem_id)->giveNodes().size();
                    outputfile << '\t' << num;
                    for ( auto const &nod : this->giveElement(elem_id)->giveNodes() ) {
                        // std::cout << "node iD = " << nodes->giveNodeId( nod ) << '\n';
                        outputfile << '\t' << nodes->giveNodeId(nod);
                    }
                }
                outputfile << '\t' << mat_stat->giveLineToSave();
                outputfile << '\t' << this->giveElement(elem_id)->giveName();  // JK not necessary, but for quick check usefull
            }
        }
        outputfile.close();
    }
}

//////////////////////////////////////////////////////////
// these methods must be separate, because at first, I need to set time of calculation in solver
// but at that time, mat_stats are still before init() so all the matStats vould be reset after that
void ElementContainer :: setFileToLoadStatsFrom(const std :: string &str) {
    this->file_to_load_from.push_back(str);
    //the commented part belonged to adaptivity feature
    /*
     * // TODO JK: make the following more universally, now it works only for LD, but the file can be named by any other name
     * if ( this->file_to_load_from.size() == 1 ) {
     *  // if LD file exists, it will be deleted, so rename it
     *  unsigned LD_num = 0;
     *  std :: string fnm = "LD";
     *  // NOTE here GlobPaths :: RESULTDIR would be usefull
     *  std :: string fnm_ini = ( masterModel->resultDir / ( fnm + ".out" ) ).string();
     *  std :: string fnm_fin;
     *  while ( LD_num < 1000 ) {
     *      fnm_fin = ( masterModel->resultDir / ( fnm + std :: to_string(LD_num) + ".out" ) ).string();
     *      if ( !fs :: exists(fnm_fin) ) {
     *          if ( fs :: exists(fnm_ini) ) {
     *              std :: rename( fnm_ini.c_str(), fnm_fin.c_str() );
     *              std :: cout << "file \'" << fnm_ini << "\' from previous calculation succesfully renamed to \'" << fnm_fin << '\'' << '\n';
     *              break;
     *          } else {
     *              std :: cerr << "could not rename \'" << fnm_ini << "\', file does not exists" << '\n';
     *          }
     *      }
     *      LD_num++;
     *  }
     * }
     */
};

//////////////////////////////////////////////////////////
void ElementContainer :: readMatStatsFromFile(double &ini_time, unsigned &ini_step, double &ini_time_step, double &ini_idc_time, const bool get_time_from_file) {
    ( void ) ini_time;
    ( void ) ini_step;
    ( void ) ini_time_step;
    ( void ) ini_idc_time;
    ( void ) get_time_from_file;
    if ( this->file_to_load_from.size() != 0 ) {
        string line, param;
        unsigned elem_id, stat_id;
        for ( auto const &file_with_stats : this->file_to_load_from ) {
            ifstream inputfile( file_with_stats.c_str() );
            if ( inputfile.is_open() ) {
                while ( getline(inputfile >> std :: ws, line) ) {
                    if ( line.empty() || ( line.at(0) == '#' ) ) {
                        continue;
                    }
                    istringstream iss(line);
                    iss >> param;
                    if ( param.compare("matStat") == 0 ) {
                        iss >> elem_id >> stat_id;
                        this->giveElement(elem_id)->giveMatStatus(stat_id)->readFromLine(iss);
                    }
                }
                inputfile.close();
            } else {
                std :: cerr << "there is no such file " << file_with_stats << '\n';
                exit(EXIT_FAILURE);
            }
        }
    }
}


//////////////////////////////////////////////////////////
void ElementContainer :: init() {
    max_sol_order = 0;
    unsigned num = 0;
    for ( vector< Element * > :: iterator e = elems.begin(); e != elems.end(); ++e, num++ ) {
        ( * e )->setID(num);
        ( * e )->init();
        ( * e )->initMaterialStatuses();
        max_sol_order = max( max_sol_order, ( * e )->giveSolutionOrder() );
    }

    //update neighborhood information
    for ( vector< Element * > :: iterator e = elems.begin(); e != elems.end(); ++e, num++ ) {
        ( * e )->collectInformationsFromNeigborhood();
    }
}

//////////////////////////////////////////////////////////
void ElementContainer :: updateMaterialStatuses() {
    if ( OmpUtils :: shouldParallelize(elems.size() ) ) {
#pragma omp parallel for schedule(static)
        for ( ptrdiff_t i = 0; i < static_cast< ptrdiff_t >( elems.size() ); i++ ) {
            elems [ i ]->updateMaterialStatuses();
        }
    } else {
        for ( vector< Element * > :: iterator e = elems.begin(); e != elems.end(); ++e ) {
            ( * e )->updateMaterialStatuses();
        }
    }
}

//////////////////////////////////////////////////////////
void ElementContainer :: resetMaterialStatuses() {
    for ( vector< Element * > :: iterator e = elems.begin(); e != elems.end(); ++e ) {
        ( * e )->resetMaterialStatuses();
    }
}


//////////////////////////////////////////////////////////
void ElementContainer :: prepareStructuralMatrix(CoordinateIndexedSparseMatrix &K, unsigned diffType, bool lumped, bool BC_applied) const {
    auto phaseStart = std :: chrono :: steady_clock :: now();
    clearStiffnessEntryMap();
    recordElementRuntimePhase(model, "matrix.prepare_pattern.clear_entry_map", std :: chrono :: duration< double >(std :: chrono :: steady_clock :: now() - phaseStart).count(), "all");

    std :: vector< Ttripletd >tripletList;

    ( void ) diffType; //not needed, matrix size is the same

    unsigned nfreeDoFs;
    if ( BC_applied ) {
        nfreeDoFs = nodes->giveTotalNumDoFs() - bconds->giveNumBlockedDoFs();
    } else {
        nfreeDoFs = nodes->giveTotalNumDoFs();
    }
    // unsigned nfreeDoFs = nodes->giveTotalNumDoFs() - bconds->giveNumBlockedDoFs();
    // cout << "nfreeDoFs " <<  nfreeDoFs << endl;

    unsigned DoFi, DoFj;
    if ( diffType == 0 ) {
        phaseStart = std :: chrono :: steady_clock :: now();
        for ( unsigned i = 0; i < constcont->giveLagrangeMultsSize(); i++ ) {
            LagrangeMultiplier *lm = constcont->giveLagrangeMultiplier(i);
            DoFi = nodes->giveDoFid( lm->giveSlaveDoF() );
            for ( unsigned j = 0; j < lm->giveNumOfDoFMasters(); j++ ) {
                DoFj = nodes->giveDoFid( lm->giveMasterDoF(j) );
                if ( DoFi < nfreeDoFs && DoFj < nfreeDoFs ) {
                    tripletList.push_back( Ttripletd(DoFi, DoFj, 0.0) );
                    tripletList.push_back( Ttripletd(DoFj, DoFi, 0.0) );
                }
            }
        }
        recordElementRuntimePhase(model, "matrix.prepare_pattern.lagrange_triplets", std :: chrono :: duration< double >(std :: chrono :: steady_clock :: now() - phaseStart).count(), "stiffness", constcont->giveLagrangeMultsSize() );
    }

    vector< unsigned >elDoFs;
    phaseStart = std :: chrono :: steady_clock :: now();
    for ( vector< Element * > :: const_iterator e = elems.begin(); e != elems.end(); ++e ) {
        elDoFs = ( * e )->giveDoFs();
        for ( unsigned i = 0; i < elDoFs.size(); i++ ) {
            for ( unsigned j = i; j < elDoFs.size(); j++ ) {
                DoFi = nodes->giveDoFid(elDoFs [ i ]);
                DoFj = nodes->giveDoFid(elDoFs [ j ]);
                //diagonal
                if ( DoFi == DoFj ) {
                    if ( DoFi < nfreeDoFs ) {
                        tripletList.push_back( Ttripletd(DoFi, DoFi, 0.0) );
                    }
                } else if ( !lumped ) {
                    //remaining items
                    if ( DoFi < nfreeDoFs && DoFj < nfreeDoFs ) {
                        tripletList.push_back( Ttripletd(DoFi, DoFj, 0.0) );
                        tripletList.push_back( Ttripletd(DoFj, DoFi, 0.0) );
                    }
                }
            }
        }
    }
    recordElementRuntimePhase(model, "matrix.prepare_pattern.element_triplets", std :: chrono :: duration< double >(std :: chrono :: steady_clock :: now() - phaseStart).count(), "all", elems.size() );

    if ( nfreeDoFs > 0 ) {
        phaseStart = std :: chrono :: steady_clock :: now();
        K.resize(nfreeDoFs, nfreeDoFs);
        K.setFromTriplets( tripletList.begin(), tripletList.end() );
        K.makeCompressed();
        recordElementRuntimePhase(model, "matrix.prepare_pattern.compress", std :: chrono :: duration< double >(std :: chrono :: steady_clock :: now() - phaseStart).count(), "all");
    }
}

//////////////////////////////////////////////////////////
void ElementContainer :: prepareStiffnessMatrix(CoordinateIndexedSparseMatrix &K) const {
    prepareStructuralMatrix(K, 0, false);
}

//////////////////////////////////////////////////////////
void ElementContainer :: prepareDampingMatrix(CoordinateIndexedSparseMatrix &C) const {
    prepareStructuralMatrix(C, 1, false);
}

//////////////////////////////////////////////////////////
void ElementContainer :: prepareMassMatrix(CoordinateIndexedSparseMatrix &M, bool lumped) const {
    prepareStructuralMatrix(M, 2, lumped);
}

//////////////////////////////////////////////////////////
void ElementContainer :: updateStructuralMatrix(CoordinateIndexedSparseMatrix &K, unsigned diffType, string matrixType, bool lumped, bool BC_applied, bool solver_numbering) const { // last two parameters only for export
    if ( diffType == 0 && !lumped && updateStiffnessMatrixParallel(K, matrixType, BC_applied, solver_numbering) ) {
        return;
    }
    updateStructuralMatrixSerial(K, diffType, matrixType, lumped, BC_applied, solver_numbering);
}

//////////////////////////////////////////////////////////
void ElementContainer :: updateStructuralMatrixSerial(CoordinateIndexedSparseMatrix &K, unsigned diffType, string matrixType, bool lumped, bool BC_applied, bool solver_numbering) const {
    if ( K.rows() == 0 ) {
        return;
    }
    auto phaseStart = std :: chrono :: steady_clock :: now();
    K = K * 0; //set everything to zero
    recordElementRuntimePhase(model, "matrix.assembly.zero_values", std :: chrono :: duration< double >(std :: chrono :: steady_clock :: now() - phaseStart).count(), matrixType);

    unsigned nfreeDoFs;
    if ( BC_applied ) {
        nfreeDoFs = nodes->giveTotalNumDoFs() - bconds->giveNumBlockedDoFs();
    } else {
        nfreeDoFs = nodes->giveTotalNumDoFs();
    }
    // unsigned nfreeDoFs = nodes->giveTotalNumDoFs() - bconds->giveNumBlockedDoFs();
    // cout << "nfreeDoFs " <<  nfreeDoFs << endl;


    unsigned DoFi, DoFj;
    vector< unsigned >elDoFs;
    Vector elDoFValues;
    Matrix k;


    if ( diffType == 0 ) {
        phaseStart = std :: chrono :: steady_clock :: now();
        for ( unsigned i = 0; i < constcont->giveLagrangeMultsSize(); i++ ) {
            LagrangeMultiplier *lm = constcont->giveLagrangeMultiplier(i);
            DoFi = nodes->giveDoFid( lm->giveSlaveDoF() );
            for ( unsigned j = 0; j < lm->giveNumOfDoFMasters(); j++ ) {
                DoFj = nodes->giveDoFid( lm->giveMasterDoF(j) );
                if ( DoFi < nfreeDoFs && DoFj < nfreeDoFs ) {
                    K.coeffRef(DoFi, DoFj) += lm->giveMasterMultiplier(j);
                    K.coeffRef(DoFj, DoFi) += lm->giveMasterMultiplier(j);
                }
            }
        }
        recordElementRuntimePhase(model, "matrix.assembly.lagrange_entries", std :: chrono :: duration< double >(std :: chrono :: steady_clock :: now() - phaseStart).count(), matrixType, constcont->giveLagrangeMultsSize() );
    }

    double localMatrixSeconds = 0.;
    double scatterSeconds = 0.;
    long long localMatrixCount = 0;
    long long scatterCount = 0;

    for ( vector< Element * > :: const_iterator e = elems.begin(); e != elems.end(); ++e ) {
        phaseStart = std :: chrono :: steady_clock :: now();
        if      ( diffType == 0 ) {
            k = ( * e )->giveStiffnessMatrix(matrixType);                    //stiffness or conductivity
        } else if ( diffType == 1 ) {
            k = ( * e )->giveDampingMatrix();                    //damping or capacity
        } else if ( diffType == 2 && lumped ) {
            k = ( * e )->giveLumpedMassMatrix();                    //mass 
        } else if ( diffType == 2 ) {
            k = ( * e )->giveMassMatrix();
        } else {
            cerr << "ElementContainer Error: time derivative matrix type " << matrixType << " unknown" << endl;
            exit(1);
        }
        localMatrixSeconds += std :: chrono :: duration< double >(std :: chrono :: steady_clock :: now() - phaseStart).count();
        localMatrixCount++;

        phaseStart = std :: chrono :: steady_clock :: now();
        elDoFs = ( * e )->giveDoFs();
        // cout << "\n  elDoFs: " ;
        for ( unsigned i = 0; i < elDoFs.size(); i++ ) {
            if ( solver_numbering ) {
                DoFi = nodes->giveDoFid(elDoFs [ i ]);
            } else {
                DoFi = elDoFs [ i ];
            }
            // cout << DoFi << " ";
            // DoFi = nodes->giveDoFid(elDoFs [ i ]);

            for ( unsigned j = i; j < elDoFs.size(); j++ ) {
                if ( solver_numbering ) {
                    DoFj = nodes->giveDoFid(elDoFs [ j ]);
                } else {
                    DoFj = elDoFs [ j ];
                }
                // DoFj = nodes->giveDoFid(elDoFs [ j ]);

                //diagonal
                if ( DoFi == DoFj ) {
                    if ( DoFi < nfreeDoFs ) {
                        K.coeffRef(DoFi, DoFi) += k(i, j);
                    }
                } else if ( !lumped ) {
                    //remaining items
                    if ( DoFi < nfreeDoFs && DoFj < nfreeDoFs ) {
                        K.coeffRef(DoFi, DoFj) += k(i, j);
                        K.coeffRef(DoFj, DoFi) += k(j, i);
                    }
                }
            }
        }
        // cout << "\n ";
        scatterSeconds += std :: chrono :: duration< double >(std :: chrono :: steady_clock :: now() - phaseStart).count();
        scatterCount++;
    }
    recordElementRuntimePhase(model, "matrix.assembly.local_matrix", localMatrixSeconds, matrixType, localMatrixCount);
    recordElementRuntimePhase(model, "matrix.assembly.scatter", scatterSeconds, matrixType, scatterCount);
    /*
     * for(size_t i=0; i<K.RowCount; i++){
     *  if (abs(K[i][i])<1E-30){         //JE:test matrix singularity
     *      cerr<< "Error in ElementContainer: stiffness matrix has zero on diagonal " << endl;
     *      exit(1);
     *  }
     * }
     */
}

/*
 * //////////////////////////////////////////////////////////
 * void ElementContainer :: updateLumpedMassMatrix(Vector &M) const {
 *  unsigned nfreeDoFs = nodes->giveTotalNumDoFs() - bconds->giveNumBlockedDoFs();
 *  M = Vector :: Zero(nfreeDoFs);
 *  vector< unsigned >elDoFs;
 *  Vector m;
 *  unsigned DoFi;
 *  for ( vector< Element * > :: const_iterator e = elems.begin(); e != elems.end(); ++e ) {
 *      m = ( * e )->giveLumpedMassMatrix();
 *      elDoFs = ( * e )->giveDoFs();
 *      for ( unsigned i = 0; i < elDoFs.size(); i++ ) {
 *          DoFi = nodes->giveDoFid(elDoFs [ i ]);
 *          if ( DoFi < nfreeDoFs ) {
 *              M [ DoFi ] += m [ i ];
 *          }
 *      }
 *  }
 * }
 */

//////////////////////////////////////////////////////////
void ElementContainer :: updateStiffnessMatrix(CoordinateIndexedSparseMatrix &K, string param) const {
    updateStructuralMatrix(K, 0, param, 0);
}

//////////////////////////////////////////////////////////
void ElementContainer :: updateDampingMatrix(CoordinateIndexedSparseMatrix &C) const {
    updateStructuralMatrix(C, 1, "", 0);
}

//////////////////////////////////////////////////////////
void ElementContainer :: updateMassMatrix(CoordinateIndexedSparseMatrix &M, bool lumped) const {
    updateStructuralMatrix(M, 2, "", lumped);
}

//////////////////////////////////////////////////////////
ptrdiff_t ElementContainer :: findSparseValueIndex(const CoordinateIndexedSparseMatrix &K, unsigned row, unsigned col) const {
    if ( row >= static_cast< unsigned >( K.rows() ) || col >= static_cast< unsigned >( K.cols() ) || !K.isCompressed() ) {
        return -1;
    }
    const auto *outer = K.outerIndexPtr();
    const auto *inner = K.innerIndexPtr();
    for ( ptrdiff_t k = outer [ col ]; k < outer [ col + 1 ]; k++ ) {
        if ( inner [ k ] == static_cast< CoordinateIndexedSparseMatrix :: StorageIndex >( row ) ) {
            return k;
        }
    }
    return -1;
}

//////////////////////////////////////////////////////////
void ElementContainer :: clearStiffnessEntryMap() const {
    stiffnessEntryMapMatrix = nullptr;
    stiffnessEntryMapRows = 0;
    stiffnessEntryMapCols = 0;
    stiffnessEntryMapNonzeros = 0;
    stiffnessEntryMap.clear();
}

//////////////////////////////////////////////////////////
bool ElementContainer :: addLagrangeMultiplierEntries(CoordinateIndexedSparseMatrix &K, unsigned nfreeDoFs) const {
    double *values = K.valuePtr();
    for ( unsigned i = 0; i < constcont->giveLagrangeMultsSize(); i++ ) {
        LagrangeMultiplier *lm = constcont->giveLagrangeMultiplier(i);
        const unsigned DoFi = nodes->giveDoFid( lm->giveSlaveDoF() );
        for ( unsigned j = 0; j < lm->giveNumOfDoFMasters(); j++ ) {
            const unsigned DoFj = nodes->giveDoFid( lm->giveMasterDoF(j) );
            if ( DoFi < nfreeDoFs && DoFj < nfreeDoFs ) {
                const ptrdiff_t ij = findSparseValueIndex(K, DoFi, DoFj);
                const ptrdiff_t ji = findSparseValueIndex(K, DoFj, DoFi);
                if ( ij < 0 || ji < 0 ) {
                    return false;
                }
                values [ ij ] += lm->giveMasterMultiplier(j);
                values [ ji ] += lm->giveMasterMultiplier(j);
            }
        }
    }
    return true;
}

//////////////////////////////////////////////////////////
bool ElementContainer :: prepareStiffnessEntryMap(const CoordinateIndexedSparseMatrix &K, bool BC_applied, bool solver_numbering) const {
    if ( stiffnessEntryMapMatrix == &K
         && stiffnessEntryMapRows == K.rows()
         && stiffnessEntryMapCols == K.cols()
         && stiffnessEntryMapNonzeros == K.nonZeros()
         && stiffnessEntryMap.size() == elems.size() ) {
        return true;
    }
    if ( !K.isCompressed() || K.rows() == 0 || K.cols() == 0 ) {
        return false;
    }

    unsigned nfreeDoFs;
    if ( BC_applied ) {
        nfreeDoFs = nodes->giveTotalNumDoFs() - bconds->giveNumBlockedDoFs();
    } else {
        nfreeDoFs = nodes->giveTotalNumDoFs();
    }

    std :: vector< std :: vector< StructuralMatrixEntry > >newMap(elems.size() );
    for ( unsigned elemIndex = 0; elemIndex < elems.size(); elemIndex++ ) {
        std :: vector< unsigned >elDoFs = elems [ elemIndex ]->giveDoFs();
        for ( unsigned i = 0; i < elDoFs.size(); i++ ) {
            const unsigned DoFi = solver_numbering ? nodes->giveDoFid(elDoFs [ i ]) : elDoFs [ i ];
            for ( unsigned j = i; j < elDoFs.size(); j++ ) {
                const unsigned DoFj = solver_numbering ? nodes->giveDoFid(elDoFs [ j ]) : elDoFs [ j ];
                if ( DoFi == DoFj ) {
                    if ( DoFi < nfreeDoFs ) {
                        const ptrdiff_t index = findSparseValueIndex(K, DoFi, DoFi);
                        if ( index < 0 ) {
                            return false;
                        }
                        newMap [ elemIndex ].push_back( {i, j, index} );
                    }
                } else if ( DoFi < nfreeDoFs && DoFj < nfreeDoFs ) {
                    const ptrdiff_t ij = findSparseValueIndex(K, DoFi, DoFj);
                    const ptrdiff_t ji = findSparseValueIndex(K, DoFj, DoFi);
                    if ( ij < 0 || ji < 0 ) {
                        return false;
                    }
                    newMap [ elemIndex ].push_back( {i, j, ij} );
                    newMap [ elemIndex ].push_back( {j, i, ji} );
                }
            }
        }
    }

    stiffnessEntryMapMatrix = &K;
    stiffnessEntryMapRows = K.rows();
    stiffnessEntryMapCols = K.cols();
    stiffnessEntryMapNonzeros = K.nonZeros();
    stiffnessEntryMap.swap(newMap);
    return true;
}

//////////////////////////////////////////////////////////
bool ElementContainer :: updateStiffnessMatrixParallel(CoordinateIndexedSparseMatrix &K, string matrixType, bool BC_applied, bool solver_numbering) const {
    if ( !OmpUtils :: shouldParallelize(elems.size() ) || !BC_applied || !solver_numbering || K.rows() == 0 || !K.isCompressed() ) {
        return false;
    }
    auto phaseStart = std :: chrono :: steady_clock :: now();
    if ( !prepareStiffnessEntryMap(K, BC_applied, solver_numbering) ) {
        return false;
    }
    recordElementRuntimePhase(model, "matrix.assembly.entry_map", std :: chrono :: duration< double >(std :: chrono :: steady_clock :: now() - phaseStart).count(), matrixType, elems.size() );

    phaseStart = std :: chrono :: steady_clock :: now();
    std :: fill(K.valuePtr(), K.valuePtr() + K.nonZeros(), 0.);
    recordElementRuntimePhase(model, "matrix.assembly.zero_values", std :: chrono :: duration< double >(std :: chrono :: steady_clock :: now() - phaseStart).count(), matrixType);

    unsigned nfreeDoFs = nodes->giveTotalNumDoFs() - bconds->giveNumBlockedDoFs();
    phaseStart = std :: chrono :: steady_clock :: now();
    if ( !addLagrangeMultiplierEntries(K, nfreeDoFs) ) {
        return false;
    }
    recordElementRuntimePhase(model, "matrix.assembly.lagrange_entries", std :: chrono :: duration< double >(std :: chrono :: steady_clock :: now() - phaseStart).count(), matrixType, constcont->giveLagrangeMultsSize() );

    double *values = K.valuePtr();
    const ptrdiff_t numElems = static_cast< ptrdiff_t >( elems.size() );
    const ptrdiff_t blockSize = 2048;
    double localMatrixSeconds = 0.;
    double scatterSeconds = 0.;
    for ( ptrdiff_t blockBegin = 0; blockBegin < numElems; blockBegin += blockSize ) {
        const ptrdiff_t blockEnd = std :: min(blockBegin + blockSize, numElems);
        std :: vector< std :: vector< double > > elementValues(blockEnd - blockBegin);

        phaseStart = std :: chrono :: steady_clock :: now();
#pragma omp parallel for schedule(dynamic)
        for ( ptrdiff_t elemIndex = blockBegin; elemIndex < blockEnd; elemIndex++ ) {
            Matrix k = elems [ elemIndex ]->giveStiffnessMatrix(matrixType);
            const std :: vector< StructuralMatrixEntry > &entries = stiffnessEntryMap [ elemIndex ];
            std :: vector< double > &localValues = elementValues [ elemIndex - blockBegin ];
            localValues.resize(entries.size() );
            for ( unsigned i = 0; i < entries.size(); i++ ) {
                localValues [ i ] = k(entries [ i ].localRow, entries [ i ].localCol);
            }
        }
        localMatrixSeconds += std :: chrono :: duration< double >(std :: chrono :: steady_clock :: now() - phaseStart).count();

        phaseStart = std :: chrono :: steady_clock :: now();
        for ( ptrdiff_t elemIndex = blockBegin; elemIndex < blockEnd; elemIndex++ ) {
            const std :: vector< StructuralMatrixEntry > &entries = stiffnessEntryMap [ elemIndex ];
            const std :: vector< double > &localValues = elementValues [ elemIndex - blockBegin ];
            for ( unsigned i = 0; i < entries.size(); i++ ) {
                values [ entries [ i ].valueIndex ] += localValues [ i ];
            }
        }
        scatterSeconds += std :: chrono :: duration< double >(std :: chrono :: steady_clock :: now() - phaseStart).count();
    }
    recordElementRuntimePhase(model, "matrix.assembly.local_matrix", localMatrixSeconds, matrixType, elems.size() );
    recordElementRuntimePhase(model, "matrix.assembly.scatter", scatterSeconds, matrixType, elems.size() );
    return true;
}

//////////////////////////////////////////////////////////
CoordinateIndexedSparseMatrix ElementContainer :: prepareOutputStiffnessMatrix(bool BC_applied) const {
    CoordinateIndexedSparseMatrix K_out;
    prepareStructuralMatrix(K_out, 0, 0, BC_applied);
    return K_out;
}

//////////////////////////////////////////////////////////
CoordinateIndexedSparseMatrix ElementContainer :: updateOutputStiffnessMatrix(CoordinateIndexedSparseMatrix K_out, string param, bool BC_applied, bool solver_numbering) const {
    updateStructuralMatrix(K_out, 0, param, 0, BC_applied, solver_numbering);
    return K_out;
}

//////////////////////////////////////////////////////////
void ElementContainer :: resetEigenStrain(double time) {
    if ( OmpUtils :: shouldParallelize(elems.size() ) ) {
#pragma omp parallel for schedule(static)
        for ( ptrdiff_t elemIndex = 0; elemIndex < static_cast< ptrdiff_t >( elems.size() ); elemIndex++ ) {
            elems [ elemIndex ]->removeEigenStrain();
        }
    } else {
        for ( auto &e : elems ) {
            e->removeEigenStrain();
        }
    }
    bconds->applyEigenStrainLoads(time);   
}

//////////////////////////////////////////////////////////
void ElementContainer :: integrateInternalForces(const Vector &full_r, Vector &full_f, bool frozen, double time, double timeStep) {
    const std :: string detail = frozen ? "frozen" : "active";
    auto phaseStart = std :: chrono :: steady_clock :: now();
    full_f.setZero();  // clear array
    recordElementRuntimePhase(model, "forces.integrate.zero_output", std :: chrono :: duration< double >(std :: chrono :: steady_clock :: now() - phaseStart).count(), detail);

    phaseStart = std :: chrono :: steady_clock :: now();
    resetEigenStrain(time);
    recordElementRuntimePhase(model, "forces.integrate.reset_eigenstrain", std :: chrono :: duration< double >(std :: chrono :: steady_clock :: now() - phaseStart).count(), detail);

    phaseStart = std :: chrono :: steady_clock :: now();
    materials->runPreparationForStressEvaluation(this);
    recordElementRuntimePhase(model, "forces.integrate.material_preparation", std :: chrono :: duration< double >(std :: chrono :: steady_clock :: now() - phaseStart).count(), detail);

    if ( OmpUtils :: shouldParallelize(elems.size() ) ) {
        phaseStart = std :: chrono :: steady_clock :: now();
        for ( unsigned so = 0; so <= max_sol_order; so++ ) {
#pragma omp parallel for schedule(dynamic)
            for ( ptrdiff_t elemIndex = 0; elemIndex < static_cast< ptrdiff_t >( elems.size() ); elemIndex++ ) {
                Element *e = elems [ elemIndex ];
                if ( e->giveSolutionOrder() != so ) {
                    continue;                                  //correct order must be used;
                }
                std :: vector< unsigned >elDoFs = e->giveDoFs();
                Vector elDoFvalues = Vector :: Zero( elDoFs.size() );
                for ( unsigned i = 0; i < elDoFs.size(); i++ ) {
                    elDoFvalues [ i ] = full_r [ elDoFs [ i ] ];
                }
                e->evaluateStrains(elDoFvalues);
                e->evaluateStresses(frozen, timeStep);
            }
        }
        recordElementRuntimePhase(model, "forces.integrate.strain_stress", std :: chrono :: duration< double >(std :: chrono :: steady_clock :: now() - phaseStart).count(), detail, elems.size() );

        const int threads = OmpUtils :: usableThreads(elems.size() );
        std :: vector< Vector >threadForces;
        threadForces.reserve(threads);
        for ( int i = 0; i < threads; i++ ) {
            threadForces.push_back( Vector :: Zero(full_f.size() ) );
        }

        phaseStart = std :: chrono :: steady_clock :: now();
#pragma omp parallel num_threads(threads)
        {
            const int thread = OmpUtils :: threadNum();
            Vector &localFullF = threadForces [ thread ];
#pragma omp for schedule(static)
            for ( ptrdiff_t elemIndex = 0; elemIndex < static_cast< ptrdiff_t >( elems.size() ); elemIndex++ ) {
                Element *e = elems [ elemIndex ];
                Vector elForces = e->giveInternalForces();
                std :: vector< unsigned >elDoFs = e->giveDoFs();
                for ( unsigned i = 0; i < elDoFs.size(); i++ ) {
                    localFullF [ elDoFs [ i ] ] += elForces [ i ];
                }
            }
        }
        recordElementRuntimePhase(model, "forces.integrate.element_force_scatter", std :: chrono :: duration< double >(std :: chrono :: steady_clock :: now() - phaseStart).count(), detail, elems.size() );

        phaseStart = std :: chrono :: steady_clock :: now();
        for ( int i = 0; i < threads; i++ ) {
            full_f += threadForces [ i ];
        }
        recordElementRuntimePhase(model, "forces.integrate.thread_force_reduce", std :: chrono :: duration< double >(std :: chrono :: steady_clock :: now() - phaseStart).count(), detail, threads);
    } else {
        Vector elDoFvalues, elForces;
        vector< unsigned >elDoFs;
        phaseStart = std :: chrono :: steady_clock :: now();
        for ( unsigned so = 0; so <= max_sol_order; so++ ) {
            for ( vector< Element * > :: iterator e = elems.begin(); e != elems.end(); ++e ) {
                if ( ( * e )->giveSolutionOrder() != so ) {
                    continue;                                  //correct order must be used;
                }
                elDoFs = ( * e )->giveDoFs();
                elDoFvalues.resize( elDoFs.size() );
                elDoFvalues.setZero();
                for ( unsigned i = 0; i < elDoFs.size(); i++ ) {
                    elDoFvalues [ i ] = full_r [ elDoFs [ i ] ];
                }
                ( * e )->evaluateStrains(elDoFvalues);
                ( * e )->evaluateStresses(frozen, timeStep);
            }
        }
        recordElementRuntimePhase(model, "forces.integrate.strain_stress", std :: chrono :: duration< double >(std :: chrono :: steady_clock :: now() - phaseStart).count(), detail, elems.size() );

        phaseStart = std :: chrono :: steady_clock :: now();
        for ( vector< Element * > :: iterator e = elems.begin(); e != elems.end(); ++e ) {
            elForces = ( * e )->giveInternalForces();
            elDoFs = ( * e )->giveDoFs();
            for ( unsigned i = 0; i < elDoFs.size(); i++ ) {
                full_f [ elDoFs [ i ] ] += elForces [ i ];
            }
        }
        recordElementRuntimePhase(model, "forces.integrate.element_force_scatter", std :: chrono :: duration< double >(std :: chrono :: steady_clock :: now() - phaseStart).count(), detail, elems.size() );
    }
}

//////////////////////////////////////////////////////////
double ElementContainer :: integrateKineticEnergy(const Vector &velocity) const {
    Vector elVelocities;
    vector< unsigned >elDoFs;
    double W_kin = 0;

    for ( vector< Element * > :: const_iterator e = elems.begin(); e != elems.end(); ++e ) {
        elDoFs = ( * e )->giveDoFs();
        elVelocities.resize( elDoFs.size() );
        for ( unsigned i = 0; i < elDoFs.size(); i++ ) {
            elVelocities [ i ] = velocity [ elDoFs [ i ] ];
        }
        W_kin += ( * e )->giveKineticEnergy(elVelocities);
    }
    return W_kin;
}

//////////////////////////////////////////////////////////
void ElementContainer :: integrateDampingOrInertiaForces(const Vector &full_v, Vector &full_f, unsigned diffType) const {
    Vector elDoFvalues, elForces;
    vector< unsigned >elDoFs;
    full_f.setZero(); // *= 0; //clear array

    for ( vector< Element * > :: const_iterator e = elems.begin(); e != elems.end(); ++e ) {
        elDoFs = ( * e )->giveDoFs();
        elDoFvalues.resize( elDoFs.size() );
        elDoFvalues.setZero();
        for ( unsigned i = 0; i < elDoFs.size(); i++ ) {
            elDoFvalues [ i ] = full_v [ elDoFs [ i ] ];
        }
        if      ( diffType == 1 ) {
            elForces = ( * e )->giveDampingMatrix() * elDoFvalues;                 //damping or conductivity
        } else if ( diffType == 2 ) {
            elForces = ( * e )->giveMassMatrix() * elDoFvalues;                    //inertia
        } else {
            cerr << "ElementContainer Error: time derivative matrix type " << diffType << " unknown" << endl;
            exit(1);
        }
        for ( unsigned i = 0; i < elDoFs.size(); i++ ) {
            full_f [ elDoFs [ i ] ] += elForces [ i ];
        }
    }
}

//////////////////////////////////////////////////////////
void ElementContainer :: integrateDampingForces(const Vector &full_v, Vector &full_f) const {
    integrateDampingOrInertiaForces(full_v, full_f, 1);
}

//////////////////////////////////////////////////////////
void ElementContainer :: integrateInertiaForces(const Vector &full_a, Vector &full_f) const {
    integrateDampingOrInertiaForces(full_a, full_f, 2);
}

//////////////////////////////////////////////////////////
void ElementContainer :: integrateInternalForces(Vector &full_r, Vector &full_f, double time, double timeStep) {
    integrateInternalForces(full_r, full_f, false, time, timeStep);
}

//////////////////////////////////////////////////////////
void ElementContainer :: integrateInternalForcesWithFrozenIntVariables(Vector &full_r, Vector &full_f, double time, double timeStep) {
    integrateInternalForces(full_r, full_f, true, time, timeStep);
}

//////////////////////////////////////////////////////////
void ElementContainer :: findElementFriends() {
    cout << "building background tetra mesh for volumetric strain calculation" << endl;
    for ( vector< Element * > :: iterator e = elems.begin(); e != elems.end(); ++e ) {
        ( * e )->findElementFriends(this);
    }
}

//////////////////////////////////////////////////////////
Element *ElementContainer :: giveElementConnectingNodes(std :: vector< unsigned > &node_ids) const {
    std :: sort( node_ids.begin(), node_ids.end() );
    // std::cout << "this elem should connect nodes";
    // for ( auto const &nid : node_ids ) {
    //   std::cout << " " << nid;
    // }
    // std::cout << '\n';

    std :: vector< unsigned >elem_node_ids;
    for ( auto const &el : this->elems ) {
        for ( auto const &nod : el->giveNodes() ) {
            for ( auto const &nID : node_ids ) {
                if ( this->nodes->giveNodeId(nod) == nID ) {
                    // std::cout << "this elem connects nodes";
                    for ( auto const &n : el->giveNodes() ) {
                        // std::cout << " " << this->nodes->giveNodeId(n);
                        elem_node_ids.push_back( this->nodes->giveNodeId(n) );
                    }
                    // std::cout << '\n';
                    if ( elem_node_ids.size() == node_ids.size() ) { ///< for other than rbc elems
                        // std::cout << "and what about here?" << '\n';
                        std :: sort( elem_node_ids.begin(), elem_node_ids.end() );
                        for ( unsigned i = 0; i < node_ids.size(); i++ ) {
                            if ( elem_node_ids [ i ] != node_ids [ i ] ) {
                                break;
                            } else {
                                if ( i == node_ids.size() - 1 ) {
                                    return el;
                                }
                            }
                        }
                    }
                    elem_node_ids.clear();
                }
            }
        }
    }
    std :: cerr << "did not find any element connecting nodes";
    for ( auto const &nid : node_ids ) {
        std :: cerr << " " << nid;
    }
    // std::cerr << '\n';
    return nullptr;
}

//////////////////////////////////////////////////////////
bool ElementContainer :: findElementOwningPoint(Element **elem, Point *xn, const Point *x) const {
    ( void ) elem;
    bool found = false;
    for ( auto &e:elems ) {
        found = e->isPointInside(xn, x);
        if ( found ) {
            * elem = e;
            return true;
        }
    }
    return false;
}

//////////////////////////////////////////////////////////
Element *ElementContainer :: findClosestElement(const Point *x) const {
    double mindist = 1e100;
    double dist;
    Element *mine = nullptr;
    Point C;
    for ( auto &e:elems ) {
        C = e->giveApproxCenter();
        dist = ( C - ( * x ) ).norm();
        if ( dist < mindist ) {
            mindist = dist;
            mine = e;
        }
    }
    return mine;
}

//////////////////////////////////////////////////////////
void ElementContainer :: extrapolateValuesFromIntegrationPointsToNodes(string code, vector< Vector > &result) const {
    //delete everythink inside
    size_t p;
    result.clear(); // result.resize(0);
    result.resize( nodes->giveSize() );
    Vector weights = Vector :: Zero( nodes->giveSize() );

    //fill with data
    vector< Vector >res;
    Vector wei;
    unsigned nodeid;
    for ( vector< Element * > :: const_iterator ee = elems.begin(); ee != elems.end(); ++ee ) {
        ( * ee )->extrapolateIPValuesToNodes(code, res, wei);
        size_t reslen = res.size();
        for ( p = 0; p < ( * ee )->giveNumOfNodes(); p++ ) {
            nodeid = ( * ee )->giveNode(p)->giveID();
            weights [ nodeid ] += wei [ p ];
            if ( reslen > ( size_t ) result [ nodeid ].size() ) {
                result [ nodeid ].resize(reslen);
                result [ nodeid ].setZero();
            }
            for ( size_t m = 0; m < min< size_t >( reslen, result [ nodeid ].size() ); m++ ) {
                result [ nodeid ] [ m ] += res [ m ] [ p ];
            }
        }
    }

    //extract pairs
    set< pair< unsigned, unsigned > >periodicPairs;
    for ( unsigned i = 0; i < model->givePreprocessingBlocks()->giveSize(); i++ ) {
        MechanicalPeriodicBC *pb = dynamic_cast< MechanicalPeriodicBC * >( model->givePreprocessingBlocks()->givePBlock(i) );
        if ( pb ) {
            vector< unsigned >masters = pb->giveMasters();
            vector< unsigned >slaves = pb->giveSlaves();
            for ( unsigned k = 0; k < masters.size(); k++ ) {
                periodicPairs.insert( make_pair(masters [ k ], slaves [ k ]) );
            }
        }
    }
    //add slaves to masters
    for ( auto q = periodicPairs.begin(); q != periodicPairs.end(); ++q ) {
        weights [ q->first ] += weights [ q->second ];
        result [ q->first ] += result [ q->second ];
    }
    //copy masters to slaves
    for ( auto q = periodicPairs.begin(); q != periodicPairs.end(); ++q ) {
        weights [ q->second ] = weights [ q->first ];
        result [ q->second ] = result [ q->first ];
    }


    //normalize by number of attached elements
    for ( p = 0; p < nodes->giveSize(); p++ ) {
        result [ p ] /= weights [ p ];
    }
}

//////////////////////////////////////////////////////////
void ElementContainer :: sumFromElements(std :: string code, Vector &result) const {
    Vector help;
    result.resize(0);
    for ( auto &e: elems ) {
        e->giveValues(code, help);
        if ( help.size() > result.size() ) {
            size_t oldsize = result.size();
            result.resize(help.size() );
            for ( size_t i = oldsize; i < ( size_t ) result.size(); i++ ) {
                result [ i ] = 0.;
            }
        }
        for ( unsigned i = 0; i < help.size(); i++ ) {
            result [ i ] += help [ i ];
        }
    }
}

//////////////////////////////////////////////////////////
void ElementContainer :: replaceTrueMassMatricesByLumpedOnes() {
    for ( auto &e: elems ) {
        e->setMassMatrix(e->giveLumpedMassMatrix() );
    }
}


//////////////////////////////////////////////////////////
void ElementContainer :: giveValues(std :: string code, Vector &result) const {
    if ( code.compare("strain_energy") == 0 || code.compare("elastic_energy") == 0 ) {
        sumFromElements("strain_energy", result);
    } else if ( code.compare("total_energy") == 0 ) {
        sumFromElements("total_energy", result);
    } else if ( code.compare("dissipated_energy") == 0 ) {
        sumFromElements("dissipated_energy", result);
    } else if ( code.compare("kinetic_energy") == 0 ) {
        sumFromElements("kinetic_energy", result);
    } else {
        result.resize(0);
    }
}


//////////////////////////////////////////////////////////
vector< Vector >ElementContainer :: computePrincipalStresses() const {
    vector< Vector >tensstress;
    string code = "stress";
    extrapolateValuesFromIntegrationPointsToNodes(code, tensstress);
    return tensstress;
}

//////////////////////////////////////////////////////////
double ElementContainer :: giveDissipatedEnergy() const {
    double E = 0;
    for ( auto &e: elems ) {
        E += e->giveDissipatedEnergy();
    }
    return E;
}


//////////////////////////////////////////////////////////
void ElementContainer :: assignFibersToElems() {
    vector< Fiber * >fibers;
    RigidBodyContact *rbc;
    Fiber *fib;
    for ( auto &e: elems ) {
        fib = dynamic_cast< Fiber * >( e );
        if ( fib ) {
            fibers.push_back(fib);
        }
    }
    if ( fibers.size() == 0 ) {
        return;
    }

    unsigned ndim = fibers [ 0 ]->giveDimension();

    Vector bbox;
    Point normal, dirvec, intersec;
    Point loc;
    Point *a, *b, *r, *s;
    Point auxA, auxB;
    double d, t, length;
    bool bintersect;
    double vol1, vol2, vol3;
    std :: vector< Node * >verts;
    for ( auto &ee:elems ) {
        rbc = dynamic_cast< RigidBodyContact * >( ee );
        if ( rbc ) {
            bbox = rbc->giveFacetBoundingBox();
            normal = rbc->giveNormal();
            loc = rbc->giveIPLoc(0);
            d = -( loc [ 0 ] * normal [ 0 ] + loc [ 1 ] * normal [ 1 ] + loc [ 2 ] * normal [ 2 ] );
            verts = rbc->giveVertices();
            for ( auto &f: fibers ) {
                a = f->giveNode(0)->givePointPointer();
                b = f->giveNode(1)->givePointPointer();

                //detect bbox intersection
                bintersect = true;
                for ( unsigned i = 0; i < ndim; i++ ) {
                    if ( min( ( * a ) [ i ], ( * b ) [ i ]) > bbox [ 2 * i + 1 ] || max( ( * a ) [ i ], ( * b ) [ i ]) < bbox [ 2 * i ] ) {
                        bintersect = false;
                        break;
                    }
                }
                if ( !bintersect ) {
                    continue;
                }

                //compute intersection with facet plane
                dirvec = f->giveDirVector();
                length = f->giveLength();
                t = -( normal [ 0 ] * ( * a ) [ 0 ] + normal [ 1 ] * ( * a ) [ 1 ] + normal [ 2 ] * ( * a ) [ 2 ] + d ) / ( normal [ 0 ] * dirvec [ 0 ] + normal [ 1 ] * dirvec [ 1 ] + normal [ 2 ] * dirvec [ 2 ] );
                if ( t <= 0 || t >= length ) {
                    continue;
                }
                intersec = ( * a ) + t * dirvec;

                //check that it is inside the facet
                //according to https://stackoverflow.com/questions/42740765/intersection-between-line-and-triangle-in-3d
                bintersect = false;
                auxA = intersec + dirvec * ( 5 * sqrt(rbc->giveArea() ) );
                auxB = intersec - dirvec * ( 5 * sqrt(rbc->giveArea() ) );
                for ( unsigned i = 0; i < verts.size() && !bintersect; i++ ) {
                    s = verts [ i ]->givePointPointer();
                    if ( i == 0 ) {
                        r = verts [ verts.size() - 1 ]->givePointPointer();
                    } else {
                        r = verts [ i - 1 ]->givePointPointer();
                    }
                    vol1 = tetraVolumeSigned(& auxA, & auxB, & loc, r);
                    vol2 = tetraVolumeSigned(& auxA, & auxB, r, s);
                    vol3 = tetraVolumeSigned(& auxA, & auxB, s, & loc);
                    if ( ( vol1 < 0 && vol2 < 0 && vol3 < 0 ) || ( vol1 > 0 && vol2 > 0 && vol3 > 0 ) ) {
                        bintersect = true;
                    }
                }
                if ( !bintersect ) {
                    continue;
                }

                f->createNewCrossing(intersec, rbc);
            }
        }
    }
    for ( auto &f:fibers ) {
        f->setUpCrossings();
    }
}
