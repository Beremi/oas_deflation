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
                } else if ( nodeType.compare("CoupledParticle") == 0 ) {
                    CoupledParticle *newnode = new CoupledParticle(dim);
                    newnode->readFromLine(iss);
                    nodes.push_back(newnode);
                } else if ( nodeType.compare("AuxNode") == 0 ) {
                    AuxNode *newnode = new AuxNode(dim);
                    newnode->readFromLine(iss);
                    nodes.push_back(newnode);
                } else if ( nodeType.compare("MechDoF") == 0 ) {
                    MechDoF *newnode = new MechDoF(dim);
                    newnode->readFromLine(iss);
                    nodes.push_back(newnode);
                } else if ( nodeType.compare("TrsDoF") == 0 ) {
                    TrsDoF *newnode = new TrsDoF(dim);
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
void NodeContainer :: saveToFile(const std :: string &filepath, std :: vector< unsigned > &nodes_to_save) const {
    std :: ofstream outputfile(filepath);
    if ( outputfile.is_open() ) {
        outputfile << "#nodes saved from calculation\n";
        for ( auto const &node_id : nodes_to_save ) {
            outputfile << this->giveNode(node_id)->giveLineToSave() << '\n';
        }
        outputfile.close();
    }
}

unsigned NodeContainer :: giveNodeId(const Node *node) const {
    // do not use this method for node that is not a part of this (nodeContainer)
    auto res = std :: find(std :: begin(this->nodes), std :: end(this->nodes), node);
    if ( res == this->nodes.end() ) {
        // if node is not in container, return zero (but zero can be also for the first node)
        // just to prevent errors here
        return 0;
    }
    return std :: distance(std :: begin(this->nodes), res);
}


//////////////////////////////////////////////////////////
void NodeContainer :: init() {
    establishDoFArray();
}

//////////////////////////////////////////////////////////
void NodeContainer :: establishDoFArray() {
    totalDoFs = 0;

    for ( vector< Node * > :: iterator n = nodes.begin(); n != nodes.end(); ++n ) {
        ( * n )->setStartingDoF(totalDoFs);
        totalDoFs += ( * n )->giveNumberOfDoFs();
    }

    BC->calculateDoFfields();
    DoFid.resize(totalDoFs);
    vector< unsigned >blocked = BC->giveArrayOfBlockedDoFs();
    loadedDoFs = BC->giveArrayOfLoadedDoFs();
    bodyForceDoFs = BC->giveArrayOfBodyForceDoFs();
    blockedDoFid.resize(blocked.size() );
    freeDoFs = totalDoFs - blocked.size();


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

    //check that there are no two Dirichlet BC assigned to one DoF
    if ( a.size() > 0 ) {
        vector< pair< unsigned, unsigned > > :: const_iterator prev = a.begin();
        for ( vector< pair< unsigned, unsigned > > :: const_iterator cur = prev + 1; cur != a.end(); ++cur ) {
            if ( prev->first == cur->first ) {
                cerr << "Node Container Error: two Dirichlet BC assigned to the same DoF number " << cur->first << endl;
                exit(1);
            }
            prev = cur;
        }
    } else  {
        cerr << "WARNING: no Dirichlet BC, the model does not prevent rigid-body motion" << endl;
    }

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
            * d = freeDoFs - constrDoFs + cs; //todo:  warning C4267: '=': conversion from 'size_t' to 'unsigned int', possible loss of data
            constrainedDoFid [ cstr [ cs ].second ] = id;
            cs++;
        } else {
            * d = id - k - cs;
        }
    }

    //identify whether the DoF is mechanical or Transport
    mechDoFs.resize(totalDoFs);
    transpDoFs.resize(totalDoFs);
    unsigned i = 0;
    unsigned ndofs;
    for ( vector< Node * > :: iterator n = nodes.begin(); n != nodes.end(); ++n ) {
        ndofs = ( * n )->giveNumberOfDoFs();
        for ( unsigned k = 0; k < ndofs; k++, i++ ) { //todo:  warning C4456: declaration of 'k' hides previous local declaration
            mechDoFs [ i ]   = ( * n )->isDoFMechanical(k);
            transpDoFs [ i ] = ( * n )->isDoFTransport(k);
        }
    }


    cout << "Loaded problem contains " << freeDoFs - constrDoFs << " DoF; additional " << constrDoFs << " DoF are dictated by constraint and "  << totalDoFs - freeDoFs << " DoF are directly prescribed" << endl;
}

//////////////////////////////////////////////////////////
void NodeContainer :: addRHS_nodalLoad(Vector &RHS, double time) const {
    vector< double >nodalLoad = BC->giveLoadedDoFValues(time);
    for ( unsigned k = 0; k < nodalLoad.size(); k++ ) {
        RHS [ loadedDoFs [ k ] ] += nodalLoad [ k ];
    }

    vector< double >bodyLoad = BC->giveBodyForceDoFValues(time);
    for ( unsigned k = 0; k < bodyLoad.size(); k++ ) {
        RHS [ bodyForceDoFs [ k ] ] += bodyLoad [ k ];
    }
}

//////////////////////////////////////////////////////////
void NodeContainer :: updateDirrichletBC(Vector &r, double time) const {
    vector< double >blocked = BC->giveBlockedDoFValues(time);
    for ( unsigned k = 0; k < blocked.size(); k++ ) {
        r [ blockedDoFid [ k ] ] = blocked [ k ];
    }
        this->giveConstraints()->calculateDependentDoFs(r, time, true);
}


//////////////////////////////////////////////////////////
void NodeContainer :: giveFullDoFArray(const Vector &fDoFs, Vector &fullDoFs) const {
    for ( unsigned i = 0; i < totalDoFs; i++ ) {
        if ( DoFid [ i ] < freeDoFs - constrDoFs ) {
            fullDoFs [ i ] = fDoFs [ DoFid [ i ] ];
        }
    }

    // #constr_new
    this->giveConstraints()->calculateDependentDoFs(fullDoFs);
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
void NodeContainer :: giveReducedForceArray(Vector &fullf, Vector &f) const {
    this->giveConstraints()->calculateMasterForces(fullf);

    for ( unsigned i = 0; i < totalDoFs; i++ ) {
        if ( DoFid [ i ] < freeDoFs - constrDoFs  ) {
            f [ DoFid [ i ] ] = fullf [ i ];
        }
    }
}

//////////////////////////////////////////////////////////
void NodeContainer :: updateExternalForcesByReactions(Vector &f_int, const Vector &load, Vector &f_dam, Vector &f_acc, Vector &f_ext) const {
    // #constr_new
    this->giveConstraints()->calculateMasterForces(f_int);
    this->giveConstraints()->calculateMasterForces(f_dam);
    this->giveConstraints()->calculateMasterForces(f_acc);


    for ( unsigned k = 0; k < totalDoFs; k++ ) {
        f_ext [ k ] = load [ k ];
        if ( DoFid [ k ] >= freeDoFs - constrDoFs ) {
            f_ext [ k ] += f_int [ k ] + f_dam [ k ] + f_acc [ k ];
        }
    }
}

//////////////////////////////////////////////////////////
Node *NodeContainer :: findClosestMechanicalNode(const Point A, double *distance) const {
    Node *closest = nullptr;
    double minDist = 1e20;
    double distance2=0;
    for ( vector< Node * > :: const_iterator n = nodes.begin(); n != nodes.end(); ++n ) {
        if ( ( * n )->doesMechanics() ) {
            distance2  = ( ( * n )->givePoint() - A ).sqNorm();
            if ( distance2 < minDist ) {
                minDist = distance2;
                closest = ( * n );
            }
        }
    }
    *distance = sqrt(minDist);
    return closest;
}

//////////////////////////////////////////////////////////
Node *NodeContainer :: findClosestAuxiliaryNode(const Point A, double *distance) const {
    Node *closest = nullptr;
    double minDist = 1e20;
    double distance2=0;
    for ( vector< Node * > :: const_iterator n = nodes.begin(); n != nodes.end(); ++n ) {
        if ( !( * n )->doesMechanics() && !( * n )->doesTransport() ) {
            distance2  = ( ( * n )->givePoint() - A ).sqNorm();
            if ( distance2 < minDist ) {
                minDist = distance2;
                closest = ( * n );
            }
        }
    }
    *distance = sqrt(minDist);
    return closest;
}

//////////////////////////////////////////////////////////
Node *NodeContainer :: findClosestTransportNode(const Point A, double *distance) const {
    Node *closest = nullptr;
    double minDist = 1e20;
    double distance2=0;
    for ( vector< Node * > :: const_iterator n = nodes.begin(); n != nodes.end(); ++n ) {
        if ( ( * n )->doesTransport() ) {
            distance2  = ( ( * n )->givePoint() - A ).sqNorm();
            if ( distance2 < minDist ) {
                minDist = distance2;
                closest = ( * n );
            }
        }
    }
    *distance = sqrt(minDist);
    return closest;
}

//////////////////////////////////////////////////////////
Node *NodeContainer :: giveNode(unsigned const num) const {
    if ( num >= nodes.size() ) {
        cerr << "NodeContainer Error: node " << num << " requested, but only " << nodes.size() << " nodes exist" << endl;
        exit(1);
    }
    return nodes [ num ];
}


//////////////////////////////////////////////////////////
unsigned NodeContainer :: giveNodeNumber(const Node* n) const {
    auto it = std::find (nodes.begin(), nodes.end(), n);
    return it-nodes.begin();
}
