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
        while ( getline(inputfile >> std::ws, line) ) {
            if ( line.empty() ){
                continue;
            }
            if ( line.at(0) == '#' ) {
                continue;
            }
            istringstream iss(line);
            iss >> elemType;
            if ( !elemType.rfind("#", 0) == 0 ) {
                if ( elemType.compare("LTCBEAM") == 0 ) {
                    RigidBodyContact *newelem = new RigidBodyContact(ndim);
                    newelem->readFromLine(iss, nodes, matrs);
                    elems.push_back(newelem);
                } else if ( elemType.compare("Truss") == 0 )    {
                    Truss *newelem = new Truss(ndim);
                    newelem->readFromLine(iss, nodes, matrs);
                    elems.push_back(newelem);
                } else if ( elemType.compare("LTCTRSP") == 0 )    {
                    Transp1D *newelem = new Transp1D(ndim);
                    newelem->readFromLine(iss, nodes, matrs);
                    elems.push_back(newelem);
                } else if ( elemType.compare("TranspPolyhedral") == 0 )    {
                    TranspPolyhedral *newelem = new TranspPolyhedral(ndim);
                    newelem->readFromLine(iss, nodes, matrs);
                    elems.push_back(newelem);
                } else if ( elemType.compare("TranspVirtPolyhedral") == 0 )    {
                    TranspVirtPolyhedral *newelem = new TranspVirtPolyhedral(ndim);
                    newelem->readFromLine(iss, nodes, matrs);
                    elems.push_back(newelem);
                } else if ( elemType.compare("TranspCondensedPolyhedral") == 0 )    {
                    TranspCondensedPolyhedral *newelem = new TranspCondensedPolyhedral(ndim);
                    newelem->readFromLine(iss, nodes, matrs);
                    elems.push_back(newelem);
                } else  {
                    cerr << "Error: element '" <<  elemType <<  "' does not exists" << endl;
                    exit(0);
                }
            }
        }
        inputfile.close();
        cout << "Input file '" <<  filename << "' succesfully loaded; " << elems.size() - origsize << " elements found" << endl;
    } else {
        cerr << "Error: unable to open input file '" <<  filename <<  "'" << endl;
        exit(0);
    }
}

//////////////////////////////////////////////////////////
void ElementContainer :: init() {
    for ( vector< Element * > :: iterator e = elems.begin(); e != elems.end(); ++e ) {
        ( * e )->init();
        ( * e )->initMaterialStatuses();
    }
}

//////////////////////////////////////////////////////////
void ElementContainer :: updateMaterialStatuses() {
    for ( vector< Element * > :: iterator e = elems.begin(); e != elems.end(); ++e ) {
        ( * e )->updateMaterialStatuses();
    }
}


//////////////////////////////////////////////////////////
void ElementContainer :: prepareSteadyStateMatrices(CoordinateIndexedSparseMatrix &K) const {

    map< pair< size_t, size_t >, double >indices11;
    //map<pair<size_t, size_t>, double> indices12;

    unsigned nfreeDoFs = nodes->giveNumFreeDoFs();
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
                } else  {
                    //remaining items
                    if ( DoFi < nfreeDoFs && DoFj < nfreeDoFs ) {
                        indices11.insert(pair< pair< size_t, size_t >, double >(pair< size_t, size_t >(DoFi, DoFj), 0.0) );
                        indices11.insert(pair< pair< size_t, size_t >, double >(pair< size_t, size_t >(DoFj, DoFi), 0.0) );
                    }//else if (DoFi < nfreeDoFs && DoFj >= nfreeDoFs) indices12.insert(pair<pair<size_t, size_t>, double>(pair<size_t, size_t > (DoFi, DoFj-nfreeDoFs), 0.0));
                     //else if (DoFj < nfreeDoFs && DoFi >= nfreeDoFs) indices12.insert(pair<pair<size_t, size_t>, double>(pair<size_t, size_t > (DoFj, DoFi-nfreeDoFs), 0.0));
                }
            }
        }
    }

    K = CoordinateIndexedSparseMatrix(indices11, nfreeDoFs, nfreeDoFs);
    //K12 = CoordinateIndexedSparseMatrix(indices12, nfreeDoFs, nDoFs-nfreeDoFs);
}

//////////////////////////////////////////////////////////
void ElementContainer :: updateSteadyStateMatrices(CoordinateIndexedSparseMatrix &K, string matrixType) const {
    // K = K * 0.;
    //K12 = K12*0.;

    unsigned nfreeDoFs = nodes->giveNumFreeDoFs();
    unsigned DoFi, DoFj;
    vector< unsigned >elDoFs;
    Vector elDoFValues;
    Matrix k;
    for ( vector< Element * > :: const_iterator e = elems.begin(); e != elems.end(); ++e ) {
        elDoFs = ( * e )->giveDoFs();
        k = ( * e )->giveSteadyStateMatrix(matrixType);
        for ( unsigned i = 0; i < elDoFs.size(); i++ ) {
            for ( unsigned j = i; j < elDoFs.size(); j++ ) {
                DoFi = nodes->giveDoFid(elDoFs [ i ]);
                DoFj = nodes->giveDoFid(elDoFs [ j ]);
                //diagonal
                if ( DoFi == DoFj ) {
                    if ( DoFi < nfreeDoFs ) {
                        K [ DoFi ] [ DoFi ] += k [ i ] [ i ];
                    }
                } else  {
                    //remaining items
                    if ( DoFi < nfreeDoFs && DoFj < nfreeDoFs ) {
                        K [ DoFi ] [ DoFj ] += k [ i ] [ j ];
                        K [ DoFj ] [ DoFi ] += k [ j ] [ i ];
                    }
                    //else if (DoFi < nfreeDoFs && DoFj >= nfreeDoFs) K12[DoFi][DoFj-nfreeDoFs] += k[i][j];
                    //else if (DoFj < nfreeDoFs && DoFi >= nfreeDoFs) K12[DoFj][DoFi-nfreeDoFs] += k[j][i];
                }
            }
        }
    }
    // cout << "Steady state matrices updated" << endl;
    if (nodes->giveConstraints()->isActive()){
      nodes->giveConstraints()->transformToConstraintSpace(K);
    }

    for(size_t i=0; i<K.RowCount; i++){
        if (abs(K[i][i])<1E-30){         //JE:test matrix singularity
            cerr<< "Error in ElementContainer: stiffness matrix has zero on diagonal " << endl; 
            exit(1);
        }
    }
}

//////////////////////////////////////////////////////////
void ElementContainer :: giveInternalForces(Vector &full_r, Vector &full_f) {
    Vector elDoFvalues, elForces;
    vector< unsigned >elDoFs;
    full_f *= 0;  // clear array

    if (nodes->giveConstraints()->isActive()){
      nodes->giveConstraints()->calculateDependentDoFs(full_r);
    }

    for ( vector< Element * > :: const_iterator e = elems.begin(); e != elems.end(); ++e ) {
        elDoFs = ( * e )->giveDoFs();
        elDoFvalues.resize(elDoFs.size() );
        for ( unsigned i = 0; i < elDoFs.size(); i++ ) {
            elDoFvalues [ i ] = full_r [ elDoFs [ i ] ];
        }
        elForces = ( * e )->giveInternalForces(elDoFvalues);
        for ( unsigned i = 0; i < elDoFs.size(); i++ ) {
            full_f [ elDoFs [ i ] ] += elForces [ i ];
        }
    }

    if (nodes->giveConstraints()->isActive()){
      nodes->giveConstraints()->calculateMasterForces(full_f);
    }
}

//////////////////////////////////////////////////////////
void ElementContainer :: addBodyForces(Vector &R, double time) const {
    //here comes distributed load, self weight
    //TO BE DONE
}
