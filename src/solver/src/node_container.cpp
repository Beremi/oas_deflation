#include "node_container.h"

//////////////////////////////////////////////////////////
NodeContainer :: ~NodeContainer() {
    for ( vector< Node * > :: iterator n = nodes.begin(); n != nodes.end(); ++n ) {
        delete * n;
    }
}

//////////////////////////////////////////////////////////
void NodeContainer :: readFromFile(const string filename, const int dim) {
    int origsize = nodes.size();
    string line, nodeType;
    ifstream inputfile(filename.c_str() );
    if ( inputfile.is_open() ) {
        while ( getline(inputfile >> std::ws, line) ) {
            if ( line.empty() ){
                cout << "EMPTY LINE" << endl;
                continue;
            }
            if ( line.at(0) == '#' ) {
                continue;
            }
            istringstream iss(line);
            iss >> std::ws >> nodeType;
            if ( !nodeType.rfind("#", 0) == 0 ) {
                if ( nodeType.compare("TrsprtNode") == 0 ) {
                    TrsNode *newnode = new TrsNode(dim);
                    newnode->readFromLine(iss, dim);
                    nodes.push_back(newnode);
                } else if ( nodeType.compare("Particle") == 0 )    {
                    Particle *newnode = new Particle(dim);
                    newnode->readFromLine(iss, dim);
                    nodes.push_back(newnode);
                } else if ( nodeType.compare("AuxNode") == 0 )    {
                    AuxNode *newnode = new AuxNode(dim);
                    newnode->readFromLine(iss, dim);
                    nodes.push_back(newnode);
                } else if ( nodeType.compare("MasterDoF") == 0 )    {
                    MasterDoF *newnode = new MasterDoF(dim);
                    newnode->readFromLine(iss, dim);
                    nodes.push_back(newnode);
                } else if ( nodeType.compare("MasterNode") == 0 )    {
                    MasterNode *newnode = new MasterNode(dim);
                    newnode->readFromLine(iss, dim);
                    nodes.push_back(newnode);
                } else  {
                    cerr << "Error: node type '" <<  nodeType <<  "' does not exists" << endl;
                    exit(0);
                }
            }
        }
        inputfile.close();
        cout << "Input file '" <<  filename << "' succesfully loaded; " << nodes.size() - origsize << " nodes found" << endl;
    } else {
        cerr << "Error: unable to open input file '" <<  filename <<  "'" << endl;
        exit(0);
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
    constrDoFs = constr->size();
    constrainedDoFid.resize(constrDoFs);
    //sort DoFs, keep track of indices
    vector< pair< unsigned, unsigned > >cstr;

    for (unsigned j = 0; j < constr->size(); j++){
      cstr.push_back(make_pair(constr->giveConstraint( j )->giveSlaveDoF(), j) );
    }
    sort(cstr.begin(), cstr.end() );
    unsigned cs = 0;  // constrained
    /////////////////////////////////////////////////////////////////

    //sort DoFs, keep track of indices
    vector< pair< unsigned, unsigned > >a;
    for ( unsigned i = 0; i < blocked.size(); i++ ) {
        a.push_back(make_pair(blocked [ i ], i) );
    }
    sort(a.begin(), a.end() );
    unsigned k = 0;
    unsigned id = 0;
    for ( vector< unsigned > :: iterator d = DoFid.begin(); d != DoFid.end(); ++d, id++ ) {
        if ( id == a [ k ].first && k < a.size()) {
          // condition < a.size() is necessary because otherwise it continues to evaluate the values from following memory (if constraint present, from constraint)
            * d = freeDoFs + k;
            blockedDoFid [ a [ k ].second ] = id;
            k++;
        } else if ( constrDoFs > 0  &&  // the later alone would not work with no constraint
                                        id == cstr [ cs ].first && cs < cstr.size()) {
            // #constraint
            * d = freeDoFs - constrDoFs + cs;
            constrainedDoFid [ cstr [ cs ].second ] = id;
            cs++;
        } else   {
            * d = id - k - cs;
        }
    }
    for ( unsigned i = 0; i < loaded.size(); i++ ) {
        loadedDoFid [ i ] = loaded [ i ];
    }
    cout << "Loaded problem contains " << freeDoFs << " degrees of freedom; additional " << totalDoFs - freeDoFs << " degrees of freedom were prescribed" << endl;
}

//////////////////////////////////////////////////////////
void NodeContainer :: addRHS_nodalLoad(Vector &RHS, double time) const {
    vector< double >loaded = BC->giveLoadedDoFValues(time);
    for ( unsigned k = 0; k < loaded.size(); k++ ) {
        // QUESTION why is here += when after every step load *= 0. Can loadedDoFid[k] indeces repeat? or can any value appear in place of other than loaded DoF?
        //  in case not, I would remove load *= 0 after every step and here replace += with just =
        RHS [ loadedDoFid [ k ] ] += loaded [ k ];
    }
}

//////////////////////////////////////////////////////////
void NodeContainer :: updateDirrichletBC(Vector &r, double time) const {
    vector< double >blocked = BC->giveBlockedDoFValues(time);
    for ( unsigned k = 0; k < blocked.size(); k++ ) {
        r [ blockedDoFid [ k ] ] = blocked [ k ];
    }
}


//////////////////////////////////////////////////////////
void NodeContainer :: giveFullDoFArray(const Vector &fDoFs, const Vector &bDoFs, Vector &fullDoFs) const {
    for ( unsigned i = 0; i < totalDoFs; i++ ) {
        if ( DoFid [ i ] < freeDoFs ) {
            fullDoFs [ i ] = fDoFs [ DoFid [ i ] ];
        } else {
            fullDoFs [ i ] = bDoFs [ DoFid [ i ] - freeDoFs ];
        }
    }
}

//////////////////////////////////////////////////////////
void NodeContainer :: giveFullDoFArray(const Vector &fDoFs, Vector &fullDoFs) const {
    for ( unsigned i = 0; i < totalDoFs; i++ ) {
        if ( DoFid [ i ] < freeDoFs ) {
            fullDoFs [ i ] = fDoFs [ DoFid [ i ] ];
        }
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
void NodeContainer :: updateExteranlForcesByReactions(const Vector &f_int, const Vector &load, Vector &f_ext) const {
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
