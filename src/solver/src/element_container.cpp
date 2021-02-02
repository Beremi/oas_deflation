#include "element_container.h"

//////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////
// ELEMENT CONATINER
ElementContainer :: ~ElementContainer() {
    for ( vector< Element * > :: iterator e = elems.begin(); e != elems.end(); ++e ) {
        delete * e;
    }
}

//////////////////////////////////////////////////////////
void ElementContainer :: readFromFile(const string filename, const unsigned ndim, MaterialContainer *matrs) {
    size_t origsize = elems.size();
    string line, elemType;
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
            iss >> elemType;
            if ( !( elemType.rfind("#", 0) == 0 ) ) {
                if ( elemType.compare("LTCBEAM") == 0 ) {
                    RigidBodyContact *newelem = new RigidBodyContact(ndim);
                    newelem->readFromLine(iss, nodes, matrs);
                    elems.push_back(newelem);
                } else if ( elemType.compare("Truss") == 0 ) {
                    Truss *newelem = new Truss(ndim);
                    newelem->readFromLine(iss, nodes, matrs);
                    elems.push_back(newelem);
                } else if ( elemType.compare("LTCTRSP") == 0 ) {
                    Transp1D *newelem = new Transp1D(ndim);
                    newelem->readFromLine(iss, nodes, matrs);
                    elems.push_back(newelem);
                } else if ( elemType.compare("LTCTRSPCoupled") == 0 ) {
                    Transp1DCoupled *newelem = new Transp1DCoupled(ndim);
                    newelem->readFromLine(iss, nodes, matrs);
                    elems.push_back(newelem);
                } else if ( elemType.compare("TrsprtQuad") == 0 ) {
                    TranspQuad *newelem = new TranspQuad();
                    newelem->readFromLine(iss, nodes, matrs);
                    elems.push_back(newelem);
                } else if ( elemType.compare("MechanicalQuad") == 0 ) {
                    MechanicalQuad *newelem = new MechanicalQuad();
                    newelem->readFromLine(iss, nodes, matrs);
                    elems.push_back(newelem);
                } else if ( elemType.compare("CosseratQuad") == 0 ) {
                    CosseratQuad *newelem = new CosseratQuad();
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
                    elems.push_back(newelem);/*
                                              * } else if ( elemType.compare("PolyhedralFace") == 0 ) {
                                              * PolyhedralFace *newelem = new PolyhedralFace(ndim);
                                              * newelem->readFromLine(iss, nodes, matrs);
                                              * elems.push_back(newelem);*/
                    /*} else if ( elemType.compare("TranspPolyhedral") == 0 ) {
                     *  TranspPolyhedral *newelem = new TranspPolyhedral(ndim);
                     *  newelem->readFromLine(iss, nodes, matrs);
                     *  elems.push_back(newelem);
                     */
                    /*
                     * } else if ( elemType.compare("TranspVirtPolyhedral") == 0 ) {
                     * TranspVirtPolyhedral *newelem = new TranspVirtPolyhedral(ndim);
                     * newelem->readFromLine(iss, nodes, matrs);
                     * elems.push_back(newelem);
                     * } else if ( elemType.compare("TranspCondensedPolyhedral") == 0 ) {
                     * TranspCondensedPolyhedral *newelem = new TranspCondensedPolyhedral(ndim);
                     * newelem->readFromLine(iss, nodes, matrs);
                     * elems.push_back(newelem);*/
                } else {
                    cerr << "Error: element '" <<  elemType <<  "' does not exists" << endl;
                    exit(EXIT_FAILURE);
                }
            }
        }
        inputfile.close();
        cout << "Input file '" <<  filename << "' succesfully loaded; " << elems.size() - origsize << " elements found" << endl;
    } else {
        cerr << "Error: unable to open input file '" <<  filename <<  "'" << endl;
        exit(EXIT_FAILURE);
    }
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

//////////////////////////////////////////////////////////
void ElementContainer :: saveElemStatsToFile(const string &filepath, const std :: vector< unsigned > &elems_to_save, const double &time_now) const {
    std :: ofstream outputfile( filepath );
    unsigned stat_id = 0;
    if ( outputfile.is_open() ) {
        outputfile << "#elem_id stat_id internal_variables...";
        outputfile << "\ntime " << time_now;
        for ( auto const &elem_id : elems_to_save) {
          stat_id = 0;
          for ( auto const &mat_stat : this->giveElement(elem_id)->giveMaterialStats()){
            // elem_id - stat_id -  mat_id - internal_variables
            outputfile << '\n' << elem_id << '\t' << stat_id++ << '\t' << mat_stat->giveMaterial()->giveId() << '\t' << mat_stat->giveLineToSave();
          }
        }
        outputfile.close();
    }
}

void ElementContainer :: readElemStatsFromFile(const string filename, const unsigned ndim, MaterialContainer *matrs) {
    string line, param;
    double time_now;
    unsigned elem_id, stat_id, mat_id;
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
            iss >> param;
            if ( param.compare("time") == 0 ){
              iss >> time_now;
            } else {
              elem_id = stoi(param);
              iss >> stat_id >> mat_id;
              this->giveElement(elem_id)->giveMatStatus(stat_id)->readFromLine(iss, matrs->giveMaterial(mat_id) );
            }

        }
        inputfile.close();
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
        max_sol_order = max(max_sol_order, ( * e )->giveSolutionOrder() );
    }
}

//////////////////////////////////////////////////////////
void ElementContainer :: updateMaterialStatuses() {
    for ( vector< Element * > :: iterator e = elems.begin(); e != elems.end(); ++e ) {
        ( * e )->updateMaterialStatuses();
    }
}


//////////////////////////////////////////////////////////
void ElementContainer :: prepareStructuralMatrix(CoordinateIndexedSparseMatrix &K, unsigned diffType) const {
    map< pair< size_t, size_t >, double >indices11;

    (void) diffType; //not needed, matrix size is the same

    unsigned nfreeDoFs = nodes->giveTotalNumDoFs() - bconds->giveNumBlockedDoFs();
    unsigned DoFi, DoFj;
    vector< unsigned >elDoFs;
    for ( vector< Element * > :: const_iterator e = elems.begin(); e != elems.end(); ++e ) {
        elDoFs = ( * e )->giveDoFs();
        for ( unsigned i = 0; i < elDoFs.size(); i++ ) {
            for ( unsigned j = i; j < elDoFs.size(); j++ ) {
                DoFi = nodes->giveDoFid(elDoFs [ i ]);
                DoFj = nodes->giveDoFid(elDoFs [ j ]);
                //diagonal
                if ( DoFi == DoFj ) {
                    if ( DoFi < nfreeDoFs ) {
                        indices11.insert(pair< pair< size_t, size_t >, double >(pair< size_t, size_t >(DoFi, DoFi), 0.0) );
                    }
                } else {
                    //remaining items
                    if ( DoFi < nfreeDoFs && DoFj < nfreeDoFs ) {
                        indices11.insert(pair< pair< size_t, size_t >, double >(pair< size_t, size_t >(DoFi, DoFj), 0.0) );
                        indices11.insert(pair< pair< size_t, size_t >, double >(pair< size_t, size_t >(DoFj, DoFi), 0.0) );
                    }
                }
            }
        }
    }
    if ( nfreeDoFs > 0 ) {
        K = CoordinateIndexedSparseMatrix(indices11, nfreeDoFs, nfreeDoFs);
    }
}

//////////////////////////////////////////////////////////
void ElementContainer :: prepareStiffnessMatrix(CoordinateIndexedSparseMatrix &K) const {
    prepareStructuralMatrix(K, 0);
}

//////////////////////////////////////////////////////////
void ElementContainer :: prepareDampingMatrix(CoordinateIndexedSparseMatrix &C) const {
    prepareStructuralMatrix(C, 1);
}

//////////////////////////////////////////////////////////
void ElementContainer :: prepareMassMatrix(CoordinateIndexedSparseMatrix &M) const {
    prepareStructuralMatrix(M, 2);
}

//////////////////////////////////////////////////////////
void ElementContainer :: updateStructuralMatrix(CoordinateIndexedSparseMatrix &K, unsigned diffType, string matrixType) const {

    unsigned nfreeDoFs = nodes->giveTotalNumDoFs() - bconds->giveNumBlockedDoFs();
    unsigned DoFi, DoFj;
    vector< unsigned >elDoFs;
    Vector elDoFValues;
    Matrix k;

    for ( vector< Element * > :: const_iterator e = elems.begin(); e != elems.end(); ++e ) {
        if      ( diffType == 0 ) k = (*e)->giveStiffnessMatrix(matrixType); //stiffness or conductivity
        else if ( diffType == 1 ) k = (*e)->giveDampingMatrix(); //damping or capacity
        else if ( diffType == 2 ) k = (*e)->giveMassMatrix(); //mass
        else {
            cerr << "ElementContainer Error: time derivative matrix type " << matrixType << " unknown" << endl;
            exit(1);
        }
        elDoFs = ( * e )->giveDoFs();
        for ( unsigned i = 0; i < elDoFs.size(); i++ ) {
            for ( unsigned j = i; j < elDoFs.size(); j++ ) {
                DoFi = nodes->giveDoFid(elDoFs [ i ]);
                DoFj = nodes->giveDoFid(elDoFs [ j ]);

                //diagonal
                if ( DoFi == DoFj ) {
                    if ( DoFi < nfreeDoFs ) {
                        K [ DoFi ] [ DoFi ] += k [ i ] [ i ];
                    }
                } else {
                    //remaining items
                    if ( DoFi < nfreeDoFs && DoFj < nfreeDoFs ) {
                        K [ DoFi ] [ DoFj ] += k [ i ] [ j ];
                        K [ DoFj ] [ DoFi ] += k [ j ] [ i ];
                    }
                }
            }
        }
    }
    if ( nodes->giveConstraints()->isActive() ) {
        nodes->giveConstraints()->transformToConstraintSpace(K);
    }

    /*
     * for(size_t i=0; i<K.RowCount; i++){
     *  if (abs(K[i][i])<1E-30){         //JE:test matrix singularity
     *      cerr<< "Error in ElementContainer: stiffness matrix has zero on diagonal " << endl;
     *      exit(1);
     *  }
     * }
     */
}

//////////////////////////////////////////////////////////
void ElementContainer :: updateStiffnessMatrix(CoordinateIndexedSparseMatrix &K, string param) const {
    updateStructuralMatrix(K, 0, param);
}

//////////////////////////////////////////////////////////
void ElementContainer :: updateDampingMatrix(CoordinateIndexedSparseMatrix &C) const {
    updateStructuralMatrix(C, 1, "");
}

//////////////////////////////////////////////////////////
void ElementContainer :: updateMassMatrix(CoordinateIndexedSparseMatrix &M) const {
    updateStructuralMatrix(M, 2, "");
}

//////////////////////////////////////////////////////////
void ElementContainer :: integrateInternalForces(const Vector &full_r, Vector &full_f, bool frozen) {
    Vector elDoFvalues, elForces;
    vector< unsigned >elDoFs;
    full_f *= 0;  // clear array

    for ( unsigned so = 0; so <= max_sol_order; so++ ) {
        for ( vector< Element * > :: iterator e = elems.begin(); e != elems.end(); ++e ) {
            if ( ( * e )->giveSolutionOrder() != so ) {
                continue;                                  //correct order must be used;
            }
            elDoFs = ( * e )->giveDoFs();
            elDoFvalues.resize(elDoFs.size() );
            for ( unsigned i = 0; i < elDoFs.size(); i++ ) {
                elDoFvalues [ i ] = full_r [ elDoFs [ i ] ];
            }
            elForces = ( * e )->giveInternalForces(elDoFvalues, frozen);
            for ( unsigned i = 0; i < elDoFs.size(); i++ ) {
                full_f [ elDoFs [ i ] ] += elForces [ i ];
            }
        }
    }
}

//////////////////////////////////////////////////////////
void ElementContainer :: integrateDampingOrInertiaForces(const Vector &full_v, Vector &full_f, unsigned diffType) const{
    Vector elDoFvalues, elForces;
    vector< unsigned >elDoFs;
    full_f *= 0; //clear array

    for ( vector< Element * > :: const_iterator e = elems.begin(); e != elems.end(); ++e ) {
        elDoFs = ( * e )->giveDoFs();
        elDoFvalues.resize(elDoFs.size() );
        for ( unsigned i = 0; i < elDoFs.size(); i++ ) {
            elDoFvalues [ i ] = full_v [ elDoFs [ i ] ];
        }
        if      (diffType==1) elForces = ( * e )->giveDampingMatrix() * elDoFvalues;  //damping or conductivity
        else if (diffType==2) elForces = ( * e )->giveMassMatrix() * elDoFvalues;  //inertia
        else {
            cerr << "ElementContainer Error: time derivative matrix type " << diffType << " unknown" << endl;
            exit(1);
        }
        for ( unsigned i = 0; i < elDoFs.size(); i++ ) {
            full_f [ elDoFs [ i ] ] += elForces [ i ];
        }
    }
}

//////////////////////////////////////////////////////////
void ElementContainer :: integrateDampingForces(const Vector &full_v, Vector &full_f) const{
    integrateDampingOrInertiaForces(full_v, full_f, 1);
}

//////////////////////////////////////////////////////////
void ElementContainer :: integrateInertiaForces(const Vector &full_a, Vector &full_f) const {
    integrateDampingOrInertiaForces(full_a, full_f, 2);
}

//////////////////////////////////////////////////////////
void ElementContainer :: integrateInternalForces(Vector &full_r, Vector &full_f) {
    integrateInternalForces(full_r, full_f, false);
}

//////////////////////////////////////////////////////////
void ElementContainer :: integrateInternalForcesWithFrozenIntVariables(Vector &full_r, Vector &full_f) {
    integrateInternalForces(full_r, full_f, true);
}

//////////////////////////////////////////////////////////
void ElementContainer :: findElementFriends() {
    for ( vector< Element * > :: iterator e = elems.begin(); e != elems.end(); ++e ) {
        ( * e )->findElementFriends(this);
    }
}
