#include "node_container.h"

//////////////////////////////////////////////////////////
NodeContainer :: ~NodeContainer() {
    for ( vector< Node * > :: iterator n = nodes.begin(); n != nodes.end(); ++n ) {
        delete * n;
    }
}

//////////////////////////////////////////////////////////
void NodeContainer :: readFromFile(const string filename, const int dim) {
    size_t origsize = nodes.size();
    string line, nodeType;
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
            iss >> std :: ws >> nodeType;
            if ( !nodeType.rfind("#", 0) == 0 ) {
                if ( nodeType.compare("TrsprtNode") == 0 ) {
                    TrsNode *newnode = new TrsNode(dim);
                    newnode->readFromLine(iss);
                    nodes.push_back(newnode);
                } else if ( nodeType.compare("MechNode") == 0 ) {
                    MechNode *newnode = new MechNode(dim);
                    newnode->readFromLine(iss);
                    nodes.push_back(newnode);
                } else if ( nodeType.compare("Particle") == 0 ) {
                    Particle *newnode = new Particle(dim);
                    newnode->readFromLine(iss);
                    nodes.push_back(newnode);
                } else if ( nodeType.compare("AuxNode") == 0 ) {
                    AuxNode *newnode = new AuxNode(dim);
                    newnode->readFromLine(iss);
                    nodes.push_back(newnode);
                } else if ( nodeType.compare("MasterDoF") == 0 ) {
                    MasterDoF *newnode = new MasterDoF(dim);
                    newnode->readFromLine(iss);
                    nodes.push_back(newnode);
                } else if ( nodeType.compare("MasterNode") == 0 ) {
                    MasterNode *newnode = new MasterNode(dim);
                    newnode->readFromLine(iss);
                    nodes.push_back(newnode);
                } else {
                    cerr << "Error: node type '" <<  nodeType <<  "' does not exists" << endl;
                    exit(EXIT_FAILURE);
                }
            }
        }
        inputfile.close();
        cout << "Input file '" <<  filename << "' succesfully loaded; " << nodes.size() - origsize << " nodes found" << endl;
    } else {
        cerr << "Error: unable to open input file '" <<  filename <<  "'" << endl;
        exit(EXIT_FAILURE);
    }
}

//////////////////////////////////////////////////////////
void NodeContainer :: init() {
    establishDoFArray();
}

//////////////////////////////////////////////////////////
void NodeContainer :: establishDoFArray() {
    totalDoFs = 0;
    freeDoFs = 0;

    for ( vector< Node * > :: iterator n = nodes.begin(); n != nodes.end(); ++n ) {
        ( * n )->setStartingDoF(totalDoFs);
        totalDoFs += ( * n )->giveNumberOfDoFs();
        freeDoFs += ( * n )->giveNumberOfFreeDoFs();
    }

    BC->calculateDoFfields();
    DoFid.resize(totalDoFs);
    vector< unsigned >blocked = BC->giveArrayOfBlockedDoFs();
    vector< unsigned >loaded = BC->giveArrayOfLoadedDoFs();
    blockedDoFid.resize(blocked.size() );
    loadedDoFid.resize(loaded.size() );

    /////////////////////////////////////////////////////////////////
    // #constraint
    constrDoFs = constr->giveSize();
    constrainedDoFid.resize(constrDoFs);
    //sort DoFs, keep track of indices
    vector< pair< unsigned, unsigned > >cstr;
    cstr.resize(constr->giveSize() );
    for ( unsigned j = 0; j < constr->giveSize(); j++ ) {
        cstr [ j ].first = constr->giveConstraint(j)->giveSlaveDoF();
        cstr [ j ].second = j;
    }
    sort(cstr.begin(), cstr.end() );

    /////////////////////////////////////////////////////////////////

    //sort DoFs, keep track of indices
    vector< pair< unsigned, unsigned > >a;
    a.resize(blocked.size() );
    for ( unsigned i = 0; i < blocked.size(); i++ ) {
        a [ i ].first = blocked [ i ];
        a [ i ].second = i;
    }
    sort(a.begin(), a.end() );

    unsigned cs = 0;
    unsigned k = 0;
    unsigned id = 0;
    for ( vector< unsigned > :: iterator d = DoFid.begin(); d != DoFid.end(); ++d, id++ ) {
        if ( k < a.size() && id == a [ k ].first ) {
            * d = freeDoFs + k;
            blockedDoFid [ a [ k ].second ] = id;
            k++;
            if ( cs < cstr.size() && id == cstr [ cs ].first ) {
                std :: cerr << "Error in establishDoFArray: cannot assign Dirichlet BC to slave node" << '\n';
                exit(1);
            }
        } else if ( cs < cstr.size() && id == cstr [ cs ].first ) {
            // #constraint
            * d = freeDoFs - constrDoFs + cs;
            constrainedDoFid [ cstr [ cs ].second ] = id;
            cs++;
        } else {
            * d = id - k - cs;
        }
    }
    for ( unsigned i = 0; i < loaded.size(); i++ ) {
        loadedDoFid [ i ] = loaded [ i ];
    }
    cout << "Loaded problem contains " << freeDoFs - constrDoFs << " DoF; additional " << constrDoFs << " DoF are dictated by constraint and "  << totalDoFs - freeDoFs << " DoF are directly prescribed" << endl;
}

//////////////////////////////////////////////////////////
void NodeContainer :: addRHS_nodalLoad(Vector &RHS, double time) const {
    vector< double >loaded = BC->giveLoadedDoFValues(time);
    for ( unsigned k = 0; k < loaded.size(); k++ ) {
        // JK QUESTION why is here += when after every step load *= 0. Can loadedDoFid[k] indeces repeat? or can any value appear in place of other than loaded DoF?
        //  in case not, I would remove load *= 0 after every step and here replace += with just =
        // JK ANSWER (based on discussion with JE): in future, more stuff will be necessary to add to load and it will be nonzero
        RHS [ loadedDoFid [ k ] ] += loaded [ k ];
    }
}

//////////////////////////////////////////////////////////
void NodeContainer :: updateDirrichletBC(Vector &r, double time) const {
    vector< double >blocked = BC->giveBlockedDoFValues(time);
    for ( unsigned k = 0; k < blocked.size(); k++ ) {
        r [ blockedDoFid [ k ] ] = blocked [ k ];
    }

    if ( this->giveConstraints()->isActive() ) {
        this->giveConstraints()->calculateDependentDoFs(r);
    }
}


//////////////////////////////////////////////////////////
void NodeContainer :: giveFullDoFArray(const Vector &fDoFs, Vector &fullDoFs) const {
    for ( unsigned i = 0; i < totalDoFs; i++ ) {
        if ( DoFid [ i ] < freeDoFs - constrDoFs ) {            
            fullDoFs [ i ] = fDoFs [ DoFid [ i ] ];
        }
    }
    // #constr_new
    if ( this->giveConstraints()->isActive() ) {
        this->giveConstraints()->calculateDependentDoFs(fullDoFs);
    }
}

//////////////////////////////////////////////////////////
void NodeContainer :: giveReducedDoFArray(const Vector &fullDoFs, Vector &fDoFs) const {
    for ( unsigned i = 0; i < totalDoFs; i++ ) {
        if ( DoFid [ i ] < freeDoFs - constrDoFs ) {
            fDoFs [ DoFid [ i ] ] = fullDoFs [ i ];
        }
    }
}

//////////////////////////////////////////////////////////
void NodeContainer :: updateExteranlForcesByReactions(Vector &f_int, const Vector &load, Vector &f_ext) const {
    // #constr_new
    if ( this->giveConstraints()->isActive() ) {
        this->giveConstraints()->calculateMasterForces(f_int);
    }

    for ( unsigned k = 0; k < totalDoFs; k++ ) {
        f_ext [ k ] = load [ k ];
        if ( DoFid [ k ] >= freeDoFs - constrDoFs ) {
            f_ext [ k ] += f_int [ k ];
        }
    }
}

//////////////////////////////////////////////////////////
Node *NodeContainer :: findClosestMechanicalNode(Point A) const {
    Node *closest = nullptr;
    double minDist = 1e20;
    double distance2;
    for ( vector< Node * > :: const_iterator n = nodes.begin(); n != nodes.end(); ++n ) {
        if ( ( * n )->doesMechanics() ) {
            distance2  = ( ( * n )->givePoint() - A ).sqNorm();
            if ( distance2 < minDist ) {
                minDist = distance2;
                closest = ( * n );
            }
        }
    }
    return closest;
}
