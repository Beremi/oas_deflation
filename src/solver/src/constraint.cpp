#include "constraint.h"
#include "node_container.h"
#include "element_container.h"
#include "element_discrete.h"
#include "boundary_condition.h"
#include "solver.h"
#include "solver_implicit.h"
#include <fstream>

using namespace std;

//////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////
// Joint Degree of Freedom
JointDoF :: JointDoF(Node *s, const unsigned &dir, const std :: vector< Node * > &m, const std :: vector< unsigned > &dirs, const std :: vector< double > &mult, const std :: vector< Function * > &fns, const std :: vector< double > &time_mult) {
    slaveNode = s;
    direction = dir;
    masters = m;
    directions = dirs;
    multipliers = mult;
    if ( fns.empty() ) {
        this->time_fns = std :: vector< Function * >(multipliers.size(), nullptr);
    } else {
        this->time_fns = fns;
        this->additional_term = time_mult;
    }
}

//////////////////////////////////////////////////////////
void JointDoF :: readFromLine(std :: istringstream &iss, NodeContainer *nodes) {
    // # type      node    direction   numMasters master1 multiplier1 master2...
    unsigned intn, intmas, dir;
    double mult;
    iss >> intn;
    slaveNode = nodes->giveNode(intn);
    iss >> direction;
    iss >> intn;
    for ( unsigned i = 0; i < intn; i++ ) {
        iss >> intmas >> dir >> mult;
        masters.push_back( nodes->giveNode(intmas) );
        directions.push_back(dir);
        time_fns.push_back(nullptr);
        multipliers.push_back(mult);
    }
}

//////////////////////////////////////////////////////////
unsigned JointDoF :: giveSlaveDoF() const {
    return slaveNode->giveStartingDoF() + direction;
}

//////////////////////////////////////////////////////////
unsigned JointDoF :: giveMasterDoF(unsigned k) const {
    return masters [ k ]->giveStartingDoF() + directions [ k ];
}
//////////////////////////////////////////////////////////
void JointDoF :: print() {
    std :: cout << "slave DoF = " << giveSlaveDoF() << '\n';
    for ( unsigned i = 0; i < masters.size(); i++ ) {
        std :: cout << "master DoF[ " << i << " ] = " << masters [ i ]->giveStartingDoF() + directions [ i ] << " with multiplier = " << multipliers [ i ];
        if ( this->time_fns [ i ] != nullptr ) {
            std :: cout << " and time_dep_mult = " << additional_term [ i ];
        }
        std :: cout << '\n';
    }
}

//////////////////////////////////////////////////////////
void JointDoF :: init(Solver *solver) {
    ( void ) solver;
    //consolidate (do not allow repetition of masters)
    for ( unsigned k = masters.size() - 1; k > 0; k-- ) {
        for ( unsigned l = 0; l < k; l++ ) {
            if ( masters [ l ] == masters [ k ] && directions [ l ] == directions [ k ] ) {
                multipliers [ l ] += multipliers [ k ];
                masters.erase(masters.begin() + k);
                directions.erase(directions.begin() + k);
                multipliers.erase(multipliers.begin() + k);
                break;
            }
        }
    }
}

//////////////////////////////////////////////////////////
bool JointDoF :: replaceDependentMasters(vector<unsigned> &depms, vector<JointDoF*> &depmsJDs){

    bool replaced = false;
    //todo: might also involve dependency on some function ...    
    
    //collect masters
    unsigned masterID;
    unsigned pos;
    for ( unsigned k = masters.size(); k>0; k-- ) {
        masterID = giveMasterDoF(k-1);
        auto posx = find(depms.begin(), depms.end(), masterID);
        if(posx != depms.end()){
            pos = posx - depms.begin();

            vector<Node*> mas = depmsJDs[pos]->giveMasterNodes();
            masters.erase(masters.begin() + k-1);
            masters.insert(masters.end(), mas.begin(), mas.end());

            vector<unsigned> dir = depmsJDs[pos]->giveMasterDirs();
            directions.erase(directions.begin() + k-1);
            directions.insert(directions.end(), dir.begin(), dir.end());
        
            vector<double> mul = depmsJDs[pos]->giveMasterMultipliers();
            for(vector<double>::const_iterator m=mul.begin(); m!=mul.end(); ++m){
                multipliers.push_back( (*m) * multipliers[k-1]);
            }
            multipliers.erase(multipliers.begin() + k-1);

            replaced = true;
        }
    }
    
    return replaced;
}

//////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////
// Lagrange multiplier
LagrangeMultiplier :: LagrangeMultiplier(const std :: vector< Node * > &m, const std :: vector< unsigned > &dirs, const std :: vector< double > &mult, const std :: vector< Function * > &fns, const std :: vector< double > &time_mult) {
    slaveNode = nullptr;
    direction = -1;
    masters = m;
    directions = dirs;
    multipliers = mult;
    if ( fns.empty() ) {
        this->time_fns = std :: vector< Function * >(multipliers.size(), nullptr);
    } else {
        this->time_fns = fns;
        this->additional_term = time_mult;
    }
}

//////////////////////////////////////////////////////////
void LagrangeMultiplier :: init(Solver *solver) {
    JointDoF :: init(solver);
}

//////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////
// Volumetric Average
VolumetricAverage :: VolumetricAverage(std :: vector< Node * > &n, std :: vector< unsigned > &d, Node *mn, unsigned md, ElementContainer *ec, ConstraintContainer *cc) {
    nodes = n;
    dirs = d;
    masternode = mn;
    masterdir = md;
    elems = ec;
    constraints = cc;

    //collect all slaves from other joint DoFs and also all involved DOFs
    std :: vector< Node * >excludedNodes;
    std :: vector< unsigned >excudedDirs;
    std :: vector< Node * >jdm;
    std :: vector< unsigned >jdd;
    JointDoF *jd;
    VolumetricAverage *va;
    for ( unsigned j = 0; j < constraints->giveConstraintsSize(); j++ ) {
        jd = constraints->giveConstraint(j);
        va = dynamic_cast< VolumetricAverage * >( jd );
        if ( !va ) {
            excludedNodes.push_back( jd->giveSlaveNode() );
            jdm = jd->giveMasterNodes();
            excludedNodes.insert( excludedNodes.end(), jdm.begin(), jdm.end() );
            excudedDirs.push_back( jd->giveSlaveDir() );
            jdd = jd->giveMasterDirs();
            excudedDirs.insert( excudedDirs.end(), jdd.begin(), jdd.end() );
        }
    }


    // find DoF to become a slave (cannot be slave or master of any other constraint)
    unsigned k;
    bool found;
    for ( k = 0; k < nodes.size(); k++ ) {
        found = false;
        for ( unsigned l = 0; l < excludedNodes.size(); l++ ) {
            if ( excludedNodes [ l ] == nodes [ k ] && excudedDirs [ l ] == dirs [ k ] ) {
                found = true;
                break;
            }
        }
        if ( !found ) {
            break;
        }
    }
    if ( found ) { //Todo potentially uninitialized variable 'found' used!
        std :: cerr << "volumetric average didn't find any suitable slave DoF" << std :: endl;
        exit(1);
    }
    slaveNode = nodes [ k ];
    direction = dirs [ k ];
}

//////////////////////////////////////////////////////////
void VolumetricAverage :: init(Solver *solver) {
    //collect all slaves from other joint DoFs
    std :: vector< unsigned >otherslaves;
    otherslaves.resize( constraints->giveConstraintsSize() );
    JointDoF *jd;
    for ( unsigned j = 0; j < constraints->giveConstraintsSize(); j++ ) {
        jd = constraints->giveConstraint(j);
        otherslaves [ j ] = jd->giveSlaveDoF();
    }

    //find slave position in array
    unsigned slaveid = find(nodes.begin(), nodes.end(), slaveNode) - nodes.begin();
    if ( slaveid == nodes.size() ) {
        std :: cerr << "Volumetric average Error: slave node not included in node list" << std :: endl;
        exit(1);
    }

    // calculate volumes associated with nodes
    unsigned r;
    unsigned s = nodes.size();
    std :: vector< double >m;
    m.resize( nodes.size() );
    for ( auto e = elems->begin(); e != elems->end(); ++e ) {
        for ( unsigned p = 0; p < ( * e )->giveNumOfNodes(); p++ ) {
            r = ( std :: find( nodes.begin(), nodes.end(), ( * e )->giveNode(p) ) - nodes.begin() );
            if ( r < s ) {
                m [ r ] += ( * e )->giveVolumeAssociatedWithNode(p);
            }
        }
    }

    //rearrange everything and update weights
    double factor = m [ slaveid ];
    double fullVolume = 0;
    for ( auto &ss: m ) {
        fullVolume += ss;
        ss /= -factor;
    }
    std :: cout << "volumetric Average: check of volume " << fullVolume << std :: endl;
    m [ slaveid ] = fullVolume / factor;
    nodes [ slaveid ] = masternode;
    dirs [ slaveid ] = masterdir;

    std :: vector< Node * >jdmasters;
    std :: vector< unsigned >jddirs;
    std :: vector< double >jdmults;
    for ( unsigned i = 0; i < nodes.size(); i++ ) {
        r = ( std :: find(otherslaves.begin(), otherslaves.end(), nodes [ i ]->giveStartingDoF() + dirs [ i ]) - otherslaves.begin() );
        if ( r == otherslaves.size() ) { //true masters
            masters.push_back(nodes [ i ]);
            directions.push_back(dirs [ i ]);
            multipliers.push_back(m [ i ]);
        } else {
            jd = constraints->giveConstraint(r);
            jdmasters = jd->giveMasterNodes();
            jddirs = jd->giveMasterDirs();
            jdmults = jd->giveMasterMultipliers();
            for ( auto &k:jdmults ) {
                k *= m [ i ];
            }
            masters.insert( masters.end(), jdmasters.begin(), jdmasters.end() );
            directions.insert( directions.end(), jddirs.begin(), jddirs.end() );
            multipliers.insert( multipliers.end(), jdmults.begin(), jdmults.end() );
        }
    }
    time_fns = std :: vector< Function * >(multipliers.size(), nullptr);

    JointDoF :: init(solver);
}


//////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////
// BC DEPENDENT ON OTHER DOFS OR CONJUGATE VARIABLES
DoFDependentOnConjugates :: DoFDependentOnConjugates(Node *s, const unsigned &dir, const std :: vector< Node * > &m, const std :: vector< unsigned > &dirs, const std :: vector< double > &mult, const std :: vector< Function * > &fns, const std :: vector< double > &time_mult) : JointDoF(s, dir, m, dirs, mult, fns, time_mult) {}


//////////////////////////////////////////////////////////
void DoFDependentOnConjugates :: init(Solver *solver) {
    JointDoF :: init(solver);
    TransientLinearTransportSolver *ts = dynamic_cast< TransientLinearTransportSolver * >( solver );
    if ( ts ) {
        std :: cerr << "*******************************" << std :: endl;
        std :: cerr << "Warning: using DoFDependentOnConjugates is strongly discouraged as the time integration might not take place at the end of time step so the update of primary variables at the end will be incorrect" << std :: endl;
        std :: cerr << "*******************************" << std :: endl;
    }
}

//////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////
// Container
//////////////////////////////////////////////////////////
void ConstraintContainer :: readFromFile(const std :: string filename, const unsigned ndim, NodeContainer *nodecont) {
    ( void ) ndim;
    std :: cout << "Input file '" <<  filename;
    unsigned origsize = constraints.size();
    unsigned lagorigsize = lagmults.size();
    std :: string line, ConstrType;
    std :: ifstream inputfile( filename.c_str() );
    if ( inputfile.is_open() ) {
        while ( getline(inputfile >> std :: ws, line) ) {
            if ( line.empty() || ( line.at(0) == '#' ) ) {
                continue;
            }
            std :: istringstream iss(line);
            iss >> std :: ws >> ConstrType;
            if ( !( ConstrType.rfind("#", 0) == 0 ) ) {
                if ( ConstrType.compare("jointDoF") == 0 || ConstrType.compare("JointDoF") == 0 ) {
                    JointDoF *newJD = new JointDoF();
                    newJD->readFromLine(iss, nodecont);
                    constraints.push_back(newJD);
                } else if ( ConstrType.compare("LagrangeMultiplier") == 0 || ConstrType.compare("lagrangemultiplier") == 0 ) {
                    LagrangeMultiplier *newLM = new LagrangeMultiplier();
                    newLM->readFromLine(iss, nodecont);
                    lagmults.push_back(newLM);
                } else {
                    std :: cerr << "Error: constraint '" <<  ConstrType <<  "' is not implemented yet." << std :: endl;
                    exit(EXIT_FAILURE);
                }
            }
        }
        std :: cout << "' succesfully loaded; " << constraints.size() - origsize << " dependent DoFs and " << lagmults.size() - lagorigsize << " Lagrangian multipliers found" << std :: endl;
        inputfile.close();
    } else {
        std :: cerr << "Error: unable to open input file '" <<  filename <<  "'" << std :: endl;
        exit(1);
    }
}

// void ConstraintContainer :: calculateSlaveDoFfield(NodeContainer *nodes){
//     for (auto &jD : constraints){
//
//     }
// }

//////////////////////////////////////////////////////////
ConstraintContainer :: ~ConstraintContainer() {
    for ( auto &c: constraints ) {
        if ( c != nullptr ) {
            delete c;
        }
    }
    for ( auto &c: lagmults ) {
        if ( c != nullptr ) {
            delete c;
        }
    }
}


void ConstraintContainer :: clear() {
    for ( auto &c: constraints ) {
        if ( c != nullptr ) {
            delete c;
        }
    }
    for ( auto &c: lagmults ) {
        if ( c != nullptr ) {
            delete c;
        }
    }
}


//////////////////////////////////////////////////////////
void ConstraintContainer :: checkInternalDependencies() {

    unsigned numM;
    vector<unsigned> allDependents(constraints.size());
    set<unsigned> allMasters;
    vector<unsigned>::const_iterator prev;
    unsigned i =0;
    vector<unsigned> depms;
    vector<JointDoF*> depmsJD;
    JointDoF *repJD;
    unsigned repJDnum;
    for ( auto const &jD : constraints ) {
        prev = find(allDependents.begin(), allDependents.begin()+i+1, jD->giveSlaveDoF());
        if (prev-allDependents.begin()<i){
            cerr << "Error in constraints, DoF " << jD->giveSlaveDoF() << " (node "<< jD->giveSlaveNode()->giveID() << " direction " << jD->giveSlaveDir() << " ) is contrained twice" << endl;
            exit(1);
        }
        allDependents[i] = jD->giveSlaveDoF();
        
        numM = jD->giveNumOfDoFMasters();
        for ( unsigned ind = 0; ind < numM; ind++ ) {
            allMasters.insert(jD->giveMasterDoF(ind));
        }
        i++;
    }
    for (set<unsigned>::iterator depm = allMasters.begin(); depm != allMasters.end(); depm++){
        prev = find(allDependents.begin(), allDependents.end(), *depm);
        if (prev!=allDependents.end()){
            repJDnum = prev  - allDependents.begin();
            repJD = constraints[repJDnum];
            depms.push_back(*depm);
            depmsJD.push_back(repJD);
        }
    }

    if (depms.size()>0){
        //first replace those that provides the dependencies, do it repetitively until no dependencies are inside hidden
        for ( auto &jD : depmsJD ) {
            while (jD->replaceDependentMasters(depms, depmsJD)) {};        
        }
        //apply replacement to all joint dofs
        for ( auto &jD : constraints ) {
            jD->replaceDependentMasters(depms, depmsJD);        
        }
    }
}


//////////////////////////////////////////////////////////
void ConstraintContainer :: init(NodeContainer *nodecont, BCContainer *bccont, Solver *solver) {
    //initiate volumetric averages

    unsigned numFreeDoFs = nodecont->giveTotalNumDoFs() - bccont->giveNumBlockedDoFs();

    nodes = nodecont;
    bconds = bccont;

    for ( auto const &jD : constraints ) {
        jD->init(solver);
    }
    VolumetricAverage *va;
    for ( int k = constraints.size() - 1; k >= 0; k-- ) {
        va = dynamic_cast< VolumetricAverage * >( constraints [ k ] );
        if ( va ) {
            constraints.push_back(va);
            constraints.erase(constraints.begin() + k);
        }
    }
    for ( auto const &LM : lagmults ) {
        LM->init(solver);
    }

    std :: vector< Ttripletd >tripletList;
    // // this should remain empty, unless you apply BC on Constrained DoF (Do not do it!)
    // // map<pair<size_t, size_t>, double> indeces21;
    // map<pair<size_t, size_t>, double> indeces22;
    if ( !this->isActive() ) {
        return;
    }

    ///////////////////////////////////////////////////
    // fill the matrix with corresponding multipliers for slaveDoFs
    unsigned i; // row index =  slave DoF
    unsigned j, numM; // column index = master DoF
    bool update;
    for ( auto const &jD : constraints ) {
        // jD->print();
        i = nodes->giveDoFid( jD->giveSlaveDoF() );

        //std::cout << jD->giveSlaveDoF() << " " << nodes->giveTotalNumDoFs() << " i = " << i << ", numFreeDoFs = " << numFreeDoFs << '\n';
        // auto res = std::find(nodes->begin(), nodes->end(), jD->giveSlaveNode());
        // std::cout << "node ID = " << std::distance(nodes->begin(), res) << '\n';
        // std::cout << "DoF num = " << jD->giveSlaveDoF() << '\n';

        if ( i < numFreeDoFs - constraints.size() ) {
            std :: cerr << "CONSTRAINT error: should never come here, constraint application unsuccesfull (hint: you are applying bondary conditions on constrained DoF) " << std :: endl;
            std :: cout << i << " " << jD->giveSlaveDoF() << " " <<  numFreeDoFs << " " << numFreeDoFs - constraints.size() << std :: endl;
            exit(1);
        } else if ( i >= numFreeDoFs ) {
            // Point A = jD->giveSlaveNode()->givePoint();
            // std::cout << "node name = " << jD->giveSlaveNode()->giveName() << '\n';
            // std::cout << "Point(" << A.getX() << ", " << A.getY() << ", " << A.getZ() << ")" << '\n';
            std :: cerr << "constraint applied simultaneously with boundary conditions" << '\n';
            exit(1);
        }

        numM = jD->giveNumOfDoFMasters();            
        for ( unsigned ind = 0; ind < numM; ind++ ) {
            j = nodes->giveDoFid( jD->giveMasterDoF(ind) );            
            if ( j < numFreeDoFs - constraints.size() ) {
                // master DoF is free
                tripletList.push_back( Ttripletd( i, j, jD->giveMasterMultiplier(ind) ) );
            }
        }
    }

    // here fill in value 1 for all other DoFs
    ///////////////////////////////////////////////////
    for ( i = 0; i < numFreeDoFs - constraints.size(); i++ ) {
        // fill the matrix with 1 for each unrestrained DoF (diagonal)
        tripletList.push_back( Ttripletd(i, i, 1) );
    }

    X.resize( numFreeDoFs, numFreeDoFs - constraints.size() );
    X.setFromTriplets( tripletList.begin(), tripletList.end() );
    X.makeCompressed();
    // // JK - left for testing:
    // for ( auto const &cn : this->constraints ){
    //     cn->print();
    // }
}

//////////////////////////////////////////////////////////
void ConstraintContainer :: initFull(NodeContainer *nodecont, BCContainer *bccont, Solver *solver) {
    //initiate volumetric averages
    unsigned numFreeDoFs = nodecont->giveTotalNumDoFs();

    nodes = nodecont;
    bconds = bccont;

    for ( auto const &jD : constraints ) {
        jD->init(solver);
    }
    VolumetricAverage *va;
    for ( int k = constraints.size() - 1; k >= 0; k-- ) {
        va = dynamic_cast< VolumetricAverage * >( constraints [ k ] );
        if ( va ) {
            constraints.push_back(va);
            constraints.erase(constraints.begin() + k);
        }
    }

    for ( auto const &LM : lagmults ) {
        LM->init(solver);
    }

    std :: vector< Ttripletd >tripletList;
    // // this should remain empty, unless you apply BC on Constrained DoF (Do not do it!)
    // // map<pair<size_t, size_t>, double> indeces21;
    // map<pair<size_t, size_t>, double> indeces22;
    if ( !this->isActive() ) {
        return;
    }

    ///////////////////////////////////////////////////
    // fill the matrix with corresponding multipliers for slaveDoFs
    unsigned i; // row index =  slave DoF
    unsigned j, numM; // column index = master DoF

    std :: vector< unsigned >NonSlaveIDs;
    for (i = 0; i < numFreeDoFs; ++i) {
        NonSlaveIDs.push_back(i);
    }
    std :: vector< unsigned >masterIDs;

    for ( auto const &jD : constraints ) {
        // jD->print();
        i = jD->giveSlaveDoF();
        NonSlaveIDs.erase( std :: remove(NonSlaveIDs.begin(), NonSlaveIDs.end(), i), NonSlaveIDs.end() );
    }

    for ( auto const &jD : constraints ) {
        i = jD->giveSlaveDoF();

        numM = jD->giveNumOfDoFMasters();
        for ( unsigned ind = 0; ind < numM; ind++ ) {
            j = jD->giveMasterDoF(ind);
            masterIDs.push_back(j);

            int n = std :: distance( NonSlaveIDs.begin(), std :: find(NonSlaveIDs.begin(), NonSlaveIDs.end(), j) ); //fins index of value "j"
            tripletList.push_back( Ttripletd( i, n, jD->giveMasterMultiplier(ind) ) );
            // }
        }
    }

    // here fill in value 1 for all other DoFs
    ///////////////////////////////////////////////////
    for ( i = 0; i < numFreeDoFs - constraints.size(); i++ ) {
        // fill the matrix with 1 for each unrestrained DoF (diagonal)
        tripletList.push_back( Ttripletd(NonSlaveIDs [ i ], i, 1) );
    }

    X_full.resize( numFreeDoFs, numFreeDoFs - constraints.size() );
    X_full.setFromTriplets( tripletList.begin(), tripletList.end() );
    X_full.makeCompressed();

    fullMasterIDs = NonSlaveIDs;
}

//////////////////////////////////////////////////////////
CoordinateIndexedSparseMatrix ConstraintContainer :: giveMatrixX(NodeContainer *nodecont, BCContainer *bccont, Solver *solver, bool BC_applied) {
    if ( BC_applied == true ) { // X matrix with BC applied, solver numbering
        return X;
    } else {
        initFull(nodecont, bccont, solver);  // X matrix without BC applied, DoFs numbering
        return X_full;
    }
}

//////////////////////////////////////////////////////////
std :: vector< unsigned int >ConstraintContainer :: giveFullMasterIDs() {
    return fullMasterIDs;
}

//////////////////////////////////////////////////////////
void ConstraintContainer :: transformToConstraintSpace(CoordinateIndexedSparseMatrix &K, const double time_now) {
    // if ( this->isTimeDependent() ) {
    //   // std::cout << "transformToConstraintSpace, time = " << time_now << '\n';
    //   this->init(this->nodes, this->bconds);
    // }
    ( void ) time_now;
    // TODO: its possible that it will work
    CoordinateIndexedSparseMatrix Knew;
    Knew = X.transpose() * K * X;
    K = Knew;
    /*
     * if ( X.cols() > 0 ) {
     *  CoordinateIndexedSparseMatrix Knew;
     *  Knew = X.transpose() * K * X;
     *  K = Knew;
     * } else {
     *  //map< pair< size_t, size_t >, double >indices11;
     *  //K = CoordinateIndexedSparseMatrix(indices11, 0, 0);
     *  //Ttripletd tripletList;
     *  //K.setFromTriplets(tripletList.begin(), tripletList.end() );
     *  //K.makeCompressed();
     * }*/
}

//////////////////////////////////////////////////////////
void ConstraintContainer :: calculateDependentDoFs(Vector &fullDoFs, const double time_now, const bool all) const {
    // bool all explanation: if false (by default), only multipliers are taken into account, if true,  also timeFunction-dependent parts
    if ( this->isActive() ) {
        for ( auto const &jD : constraints ) {
            fullDoFs [ jD->giveSlaveDoF() ] = 0; // to be sure that there is zero value
            //standard constraint, imposed relations between dofs
            for ( unsigned i = 0; i < jD->giveNumOfDoFMasters(); i++ ) {
                fullDoFs [ jD->giveSlaveDoF() ] += fullDoFs [ jD->giveMasterDoF(i) ] * jD->giveMasterMultiplier(i) + ( all ? jD->giveFnDepPart(i, time_now) : 0.0 );
            }
        }
    }
}

//////////////////////////////////////////////////////////
void ConstraintContainer :: calculateDoFsDependentOnConjugates(Vector &full_ddr, const Vector &trial_r, const Vector &fullFExt) const {
    // bool all explanation: if false (by default), only multipliers are taken into account, if true,  also timeFunction-dependent parts
    if ( this->isActive() ) {
        for ( auto const &jD : constraints ) {
            //full_ddr [ jD->giveSlaveDoF() ] = 0; // to be sure that there is zero value
            //constraint relating DoFs to conjugate variables
            for ( unsigned i = 0; i < jD->giveNumOfConjugateMasters(); i++ ) {
                full_ddr [ jD->giveSlaveDoF() ] += fullFExt [ jD->giveMasterDoF(i) ] * jD->giveMasterMultiplier(i);
            }
            if ( jD->giveNumOfConjugateMasters() > 0 ) {
                full_ddr [ jD->giveSlaveDoF() ] -= trial_r [ jD->giveSlaveDoF() ];   //this is needed because ddr will be later summed with trial_r, but the calculated pressures are already absolute (after the summation)
            }
        }
    }
}

//////////////////////////////////////////////////////////
void ConstraintContainer :: calculateMasterForces(Vector &fullForces) {
    if ( this->isActive() ) {
        // std::cout << "calculateMasterForces, time = " << time_now << '\n';
        for ( auto const &jD : constraints ) {
            for ( unsigned i = 0; i < jD->giveNumOfDoFMasters(); i++ ) {
                fullForces [ jD->giveMasterDoF(i) ] += fullForces [ jD->giveSlaveDoF() ] * ( jD->giveMasterMultiplier(i) );
            }
            fullForces [ jD->giveSlaveDoF() ] = 0;
        }
    }
}


//////////////////////////////////////////////////////////
void ConstraintContainer :: removeConstraint(unsigned i) {
    if ( i > constraints.size() - 1 ) {
        std :: cerr << "ConstraintContainer Error: requester constraint number " << i << " out of " << constraints.size() << std :: endl;
        exit(1);
    }
    delete constraints [ i ];
    constraints [ i ] = nullptr;
    constraints.erase(constraints.begin() + i);
}

//////////////////////////////////////////////////////////
void ConstraintContainer :: addRigidArmConstraint(unsigned dim, Node* dependent, Node* primary, bool includeRotations){
    Point diff = dependent->givePoint() - primary->givePoint();
    vector< Node * >vm;
    vector< unsigned >dirs;
    vector< double >mults;    
    JointDoF* jd;

    if ( dim == 3 ) {
        vm.resize(3,primary);
        mults.resize(3);
        dirs.resize(3);
        //direction X
        mults[0] = 1.;
        mults[1] = diff[2];
        mults[2] = -diff[1];
        dirs[0] = 0;
        dirs[1] = 4;
        dirs[2] = 5;
        jd = new JointDoF(dependent, dirs [ 0 ], vm, dirs, mults);
        addConstraint(jd);

        //direction Y
        mults[1] = -diff[2];
        mults[2] = diff[0];
        dirs[0] = 1;
        dirs[1] = 3;
        dirs[2] = 5;
        jd = new JointDoF(dependent, dirs [ 0 ], vm, dirs, mults);
        addConstraint(jd);

        //direction Z
        mults[1] = diff[1];
        mults[2] = -diff[0];
        dirs[0] = 2;
        dirs[1] = 3;
        dirs[2] = 4;
        jd = new JointDoF(dependent, dirs [ 0 ], vm, dirs, mults);
        addConstraint(jd);

        if (includeRotations){
            vm.resize(1,primary);
            mults.resize(1,1);
            dirs.resize(1,3);
            jd = new JointDoF(dependent, dirs [ 0 ], vm, dirs, mults);
            addConstraint(jd);
            dirs.resize(1,4);
            jd = new JointDoF(dependent, dirs [ 0 ], vm, dirs, mults);
            addConstraint(jd);
            dirs.resize(1,5);
            jd = new JointDoF(dependent, dirs [ 0 ], vm, dirs, mults);
            addConstraint(jd);
        }

    } else if( dim == 2 ) {
        vm.resize(2,primary);
        mults.resize(2);
        dirs.resize(2);
        //direction X
        mults[0] = 1.;
        mults[1] = -diff[1];
        dirs[0] = 0;
        dirs[1] = 2;
        jd = new JointDoF(dependent, dirs [ 0 ], vm, dirs, mults);
        addConstraint(jd);

        //direction Y
        mults[0] = 1.;
        mults[1] = diff[0];
        dirs[0] = 1;
        jd = new JointDoF(dependent, dirs [ 0 ], vm, dirs, mults);
        addConstraint(jd);

        if (includeRotations){
            vm.resize(1,primary);
            mults.resize(1,1);
            dirs.resize(1,2);
            jd = new JointDoF(dependent, dirs [ 0 ], vm, dirs, mults);
            addConstraint(jd);
        }
    }
}
