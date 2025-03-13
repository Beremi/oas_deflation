#include "node_container.h"
#include "solver.h"

using namespace std;

//////////////////////////////////////////////////////////
NodeContainer :: ~NodeContainer() {
    for ( vector< Node * > :: iterator n = nodes.begin(); n != nodes.end(); ++n ) {
        if ( * n != nullptr ) {
            delete * n;
        }
    }
}

//////////////////////////////////////////////////////////
void NodeContainer :: clear() {
    for ( vector< Node * > :: iterator n = nodes.begin(); n != nodes.end(); ++n ) {
        if ( * n != nullptr ) {
            delete * n;
        }
    }
    totalDoFs = 0;
}

//////////////////////////////////////////////////////////
void NodeContainer :: addNode(Node *n) {
    n->init();
    n->setID(nodes.size() );
    n->setStartingDoF(totalDoFs);
    nodes.push_back(n);
    totalDoFs += n->giveNumberOfDoFs();
}

//////////////////////////////////////////////////////////
void NodeContainer :: readFromFile(const string filename, const int dim) {
    cout << "Input file '" <<  filename;
    size_t origsize = nodes.size();
    string line, nodeType;
    ifstream inputfile( filename.c_str() );
    if ( inputfile.is_open() ) {
        while ( getline(inputfile >> std :: ws, line) ) {
            if ( line.empty() || ( line.at(0) == '#' ) ) {
                continue;
            }
            istringstream iss(line);
            iss >> std :: ws >> nodeType;
            if ( !( nodeType.rfind("#", 0) == 0 ) ) {
                if ( nodeType.compare("TrsprtNode") == 0 ) {
                    TrsNode *newnode = new TrsNode(dim);
                    newnode->readFromLine(iss);
                    addNode(newnode);
                } else if ( nodeType.compare("MechNode") == 0 ) {
                    MechNode *newnode = new MechNode(dim);
                    newnode->readFromLine(iss);
                    addNode(newnode);
                } else if ( nodeType.compare("Particle") == 0 ) {
                    Particle *newnode = new Particle(dim);
                    newnode->readFromLine(iss);
                    addNode(newnode);
                } else if ( nodeType.compare("ParticleWithTransport") == 0 ) {
                    ParticleWithTransport *newnode = new ParticleWithTransport(dim);
                    newnode->readFromLine(iss);
                    addNode(newnode);
                } else if ( nodeType.compare("AuxNode") == 0 ) {
                    AuxNode *newnode = new AuxNode(dim);
                    newnode->readFromLine(iss);
                    addNode(newnode);
                } else if ( nodeType.compare("MechDoF") == 0 ) {
                    MechDoF *newnode = new MechDoF(dim, 0);
                    newnode->readFromLine(iss);
                    addNode(newnode);
                } else if ( nodeType.compare("TrsDoF") == 0 ) {
                    TrsDoF *newnode = new TrsDoF(dim, 0);
                    newnode->readFromLine(iss);
                    addNode(newnode);
                } else if ( nodeType.compare("TempDoF") == 0 ) {
                    TempDoF *newnode = new TempDoF(dim, 0);
                    newnode->readFromLine(iss);
                    addNode(newnode);
                } else if ( nodeType.compare("HumidityDoF") == 0 ) {
                    HumidityDoF *newnode = new HumidityDoF(dim, 0);
                    newnode->readFromLine(iss);
                    addNode(newnode);
                } else if ( nodeType.compare("TempNode") == 0 ) {
                    TempNode *newnode = new TempNode(dim);
                    newnode->readFromLine(iss);
                    addNode(newnode);
                } else if ( nodeType.compare("ParticleWithTemperature") == 0 ) {
                    ParticleWithTemperature *newnode = new ParticleWithTemperature(dim);
                    newnode->readFromLine(iss);
                    addNode(newnode);
                } else {
                    cerr << "Error: node type '" <<  nodeType <<  "' does not exists" << endl;
                    exit(EXIT_FAILURE);
                }
            }
        }
        inputfile.close();
        cout << "' succesfully loaded; " << nodes.size() - origsize << " nodes found" << endl;
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

//////////////////////////////////////////////////////////
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
void NodeContainer :: initSimplices() {
    for ( auto &n:nodes ) {
        n->initSimplex();
    }
    Simplex *s;
    for ( auto &n:nodes ) {
        s = n->giveSimplex();
        if ( s && !s->isValid() ) {
            s->findNeighbors(this);
        }
    }
}

//////////////////////////////////////////////////////////
void NodeContainer :: updateSimplexVolumetricStrains(const Vector &fullDoFs) {
    //update valid simplices
    for ( auto &n:nodes ) {
        if ( n->hasSimplex() ) {
            n->updateSimplexVolumetricStrain(fullDoFs);
        }
    }
    //update from neighbours
    for ( auto &n:nodes ) {
        if ( n->hasSimplex() ) {
            n->stealSimplexVolumetricStrain();
        }
    }
}

//////////////////////////////////////////////////////////
void NodeContainer :: establishDoFArray() {
    BC->calculateDoFfields();
    DoFid.resize(totalDoFs);
    DoF2nodes.resize(totalDoFs);
    vector< unsigned >blocked = BC->giveArrayOfBlockedDoFs();
    loadedDoFs = BC->giveArrayOfLoadedDoFs();
    bodyForceDoFs = BC->giveArrayOfBodyForceDoFs();
    blockedDoFid.resize( blocked.size() );

    /////////////////////////////////////////////////////////////////
    // #constraint
    constrDoFs = constr->giveConstraintsSize();
    constrainedDoFid.resize(constrDoFs);
    //sort DoFs, keep track of indices
    vector< pair< unsigned, unsigned > >cstr;
    cstr.resize( constr->giveConstraintsSize() );
    for ( unsigned j = 0; j < constr->giveConstraintsSize(); j++ ) {
        cstr [ j ].first = constr->giveConstraint(j)->giveSlaveDoF();
        cstr [ j ].second = j;
    }
    sort( cstr.begin(), cstr.end() );

    /////////////////////////////////////////////////////////////////
    freeDoFs = totalDoFs - constrDoFs - blocked.size();

    invDoFid.resize(freeDoFs);

    //sort DoFs, keep track of indices
    vector< pair< unsigned, unsigned > >a;
    a.resize( blocked.size() );
    for ( unsigned i = 0; i < blocked.size(); i++ ) {
        a [ i ].first = blocked [ i ];
        a [ i ].second = i;
    }
    sort( a.begin(), a.end() );

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
    } else {
        cerr << "WARNING: no Dirichlet BC, the model does not prevent rigid-body motion" << endl;
    }

    unsigned cs = 0;
    unsigned k = 0;
    unsigned id = 0;
    for ( vector< unsigned > :: iterator d = DoFid.begin(); d != DoFid.end(); ++d, id++ ) {
        if ( k < a.size() && id == a [ k ].first ) {
            * d = freeDoFs + constrDoFs + k;
            blockedDoFid [ a [ k ].second ] = id;
            k++;
            if ( cs < cstr.size() && id == cstr [ cs ].first ) {
                std :: cerr << "Error in establishDoFArray: cannot assign Dirichlet BC to slave node" << '\n';
                exit(1);
            }
        } else if ( cs < cstr.size() && id == cstr [ cs ].first ) {
            // #constraint
            * d = freeDoFs + cs; //todo:  warning C4267: '=': conversion from 'size_t' to 'unsigned int', possible loss of data
            constrainedDoFid [ cstr [ cs ].second ] = id;
            cs++;
        } else {
            * d = id - k - cs;
            invDoFid [ * d ] = id;
        }
    }

    //identify physical fields of DoFs
    physicalFieldsDoF.resize(totalDoFs);
    unsigned i = 0;
    vector< unsigned >nodePhysFields;
    for ( auto &n:nodes ) {
        nodePhysFields = n->givePhysicalFieldNumForAllDoFs();
        for ( auto &q:nodePhysFields ) {
            physicalFieldsDoF [ i ]   = q;
            DoF2nodes [ i ] = n->giveID();
            i++;
        }
    }

    cout << "Loaded problem contains " << freeDoFs << " DoF; additional " << constrDoFs << " DoF are dictated by constraint and "  << totalDoFs - freeDoFs - constrDoFs << " DoF are directly prescribed" << endl;
}

//////////////////////////////////////////////////////////
void NodeContainer :: addRHS_nodalLoad(Vector &f, double time) const {
    vector< double >bodyLoad = BC->giveBodyForceDoFValues(time);
    for ( unsigned k = 0; k < bodyLoad.size(); k++ ) {
        f [ bodyForceDoFs [ k ] ] += bodyLoad [ k ];
    }

    vector< double >nodalLoad = BC->giveLoadedDoFValues(time);
    for ( unsigned k = 0; k < nodalLoad.size(); k++ ) {
        f [ loadedDoFs [ k ] ] += nodalLoad [ k ];
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
        if ( DoFid [ i ] < freeDoFs ) {
            fullDoFs [ i ] = fDoFs [ DoFid [ i ] ];
        }
    }
    this->giveConstraints()->calculateDependentDoFs(fullDoFs);
}

//////////////////////////////////////////////////////////
void NodeContainer :: updateFullDoFsByDependenciesOnConjugates(Vector &ddr, const Vector &trial_r, const Vector &f_ext) const {  // accounts also for constraints between master and conjugate variables
    this->giveConstraints()->calculateDoFsDependentOnConjugates(ddr, trial_r, f_ext);
}


//////////////////////////////////////////////////////////
void NodeContainer :: giveReducedDoFArray(const Vector &fullDoFs, Vector &fDoFs) const {
    for ( unsigned i = 0; i < totalDoFs; i++ ) {
        if ( DoFid [ i ] < freeDoFs ) {
            fDoFs [ DoFid [ i ] ] = fullDoFs [ i ];
        }
    }
}

//////////////////////////////////////////////////////////
void NodeContainer :: giveReducedForceArray(Vector &fullf, Vector &f) const {
    this->giveConstraints()->calculateMasterForces(fullf);

    for ( unsigned i = 0; i < totalDoFs; i++ ) {
        if ( DoFid [ i ] < freeDoFs ) {
            f [ DoFid [ i ] ] = fullf [ i ];
        }
    }
}

//////////////////////////////////////////////////////////
void NodeContainer :: updateExternalForcesByReactions(Vector &f_int, Vector &load, Vector &f_dam, Vector &f_acc, Vector &f_ext, const Vector &full_r) const {
    // #constr_new
    this->giveConstraints()->calculateMasterForces(f_int);
    this->giveConstraints()->calculateMasterForces(f_dam);
    this->giveConstraints()->calculateMasterForces(f_acc);
    this->giveConstraints()->calculateMasterForces(load);

    for ( unsigned k = 0; k < totalDoFs; k++ ) {
        f_ext [ k ] = load [ k ];
        if ( DoFid [ k ] >= freeDoFs + constrDoFs ) {
            f_ext [ k ] = f_int [ k ] + f_dam [ k ] + f_acc [ k ];
        }
    }

    double lag;
    for ( unsigned i = 0; i < this->giveConstraints()->giveLagrangeMultsSize(); i++ ) {
        LagrangeMultiplier *lm = this->giveConstraints()->giveLagrangeMultiplier(i);
        lag = full_r [ lm->giveSlaveDoF() ];
        for ( unsigned j = 0; j < lm->giveNumOfDoFMasters(); j++ ) {
            f_ext [ lm->giveMasterDoF(j) ] -= lag * lm->giveMasterMultiplier(j);
        }
    }
}

//////////////////////////////////////////////////////////
Node *NodeContainer :: findClosestMechanicalNode(const Point A, double *distance) const {
    Node *closest = nullptr;
    double minDist = 1e20;
    double distance2 = 0;
    for ( vector< Node * > :: const_iterator n = nodes.begin(); n != nodes.end(); ++n ) {
        if ( ( * n )->doesMechanics() ) {
            distance2  = ( ( * n )->givePoint() - A ).squaredNorm();
            if ( distance2 < minDist ) {
                minDist = distance2;
                closest = ( * n );
            }
        }
    }
    * distance = sqrt(minDist);
    return closest;
}

//////////////////////////////////////////////////////////
Node *NodeContainer :: findClosestAuxiliaryNode(const Point A, double *distance) const {
    Node *closest = nullptr;
    double minDist = 1e20;
    double distance2 = 0;
    for ( vector< Node * > :: const_iterator n = nodes.begin(); n != nodes.end(); ++n ) {
        if ( !( * n )->doesMechanics() && !( * n )->doesTransport() ) {
            distance2  = ( ( * n )->givePoint() - A ).squaredNorm();
            if ( distance2 < minDist ) {
                minDist = distance2;
                closest = ( * n );
            }
        }
    }
    * distance = sqrt(minDist);
    return closest;
}

//////////////////////////////////////////////////////////
Node *NodeContainer :: findClosestTransportNode(const Point A, double *distance) const {
    Node *closest = nullptr;
    double minDist = 1e20;
    double distance2 = 0;
    for ( vector< Node * > :: const_iterator n = nodes.begin(); n != nodes.end(); ++n ) {
        if ( ( * n )->doesTransport() ) {
            distance2  = ( ( * n )->givePoint() - A ).squaredNorm();
            if ( distance2 < minDist ) {
                minDist = distance2;
                closest = ( * n );
            }
        }
    }
    * distance = sqrt(minDist);
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
unsigned NodeContainer :: giveNodeNumber(const Node *n) const {
    auto it = std :: find(nodes.begin(), nodes.end(), n);
    return it - nodes.begin();
}

//////////////////////////////////////////////////////////
Vector NodeContainer :: readInitialConditions(string initfile) const {
    string line;
    unsigned numi, startDoF;
    double numd;
    Vector initvalues = Vector :: Zero(totalDoFs);
    ifstream inputfile( initfile.c_str() );
    if ( inputfile.is_open() ) {
        while ( getline(inputfile >> std :: ws, line) ) {
            istringstream iss(line);
            iss >> numi;
            startDoF = nodes [ numi ]->giveStartingDoF();
            for ( unsigned v = 0; v < nodes [ numi ]->giveNumberOfDoFs(); v++ ) {
                iss >> numd;
                initvalues [ startDoF + v ] = numd;
            }
        }
        inputfile.close();

        Vector initreduced = Vector :: Zero(freeDoFs);
        giveReducedDoFArray(initvalues, initreduced);// to propagate intial master field through constraints
        giveFullDoFArray(initreduced, initvalues);

        cout << "Input file '" <<  initfile << "' succesfully loaded" << endl;
    } else {
        cerr << "Error: unable to open input file '" <<  initfile <<  "'" << endl;
        exit(EXIT_FAILURE);
    }
    return initvalues;
}
