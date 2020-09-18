#include "boundary_condition.h"
#include "node_container.h"
#include "linear_algebra.h"

//TODO one currently cannot have two boundary conditions assigned to one node

//////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////
// DIRICHLET AND NEUMANN BOUNDARY CONDITION
void BoundaryCondition :: init(FunctionContainer *funcs) {

    blockedDoFNum = 0;
    loadedDoFNum = 0;

    dirichF.resize(dirichBC.size());
    for ( unsigned i = 0; i < dirichBC.size(); i++ ) {
        if ( dirichBC[i] >= 0 ) {
            blockedDoFNum++;
            dirichF[i] = funcs->giveFunction(dirichBC[i]);            
        }else{
            dirichF[i] = nullptr;
        }
    }

    neumannF.resize(neumannBC.size());
    for ( unsigned i = 0; i < neumannBC.size(); i++ ) {
        if ( neumannBC[i] >= 0 ) {
            loadedDoFNum++;
            neumannF[i] = funcs->giveFunction(neumannBC[i]);            
        }else{
            neumannF[i] = nullptr;
        }
    }
}

//////////////////////////////////////////////////////////
vector< unsigned >BoundaryCondition :: giveBlockedDoFs() const {
    vector< unsigned >blocked;
    blocked.resize(blockedDoFNum);
    unsigned k = 0;
    unsigned s = 0;
    for ( vector< int > :: const_iterator i = dirichBC.begin(); i != dirichBC.end(); ++i, k++ ) {
        if ( * i >= 0 ) {
            blocked [ s ] = node->giveStartingDoF() + k;
            s++;
        }
    }
    ;
    return blocked;
}

//////////////////////////////////////////////////////////
vector< unsigned >BoundaryCondition :: giveLoadedDoFs() const {
    vector< unsigned >loaded;
    loaded.resize(loadedDoFNum);
    unsigned k = 0;
    unsigned s = 0;
    for ( vector< int > :: const_iterator i = neumannBC.begin(); i != neumannBC.end(); ++i, k++ ) {
        if ( * i >= 0 ) {
            loaded [ s ] = node->giveStartingDoF() + k;
            s++;
        }
    }
    ;
    return loaded;
}

//////////////////////////////////////////////////////////
vector< double >BoundaryCondition :: giveBlockedDoFValues(double t) const {
    vector< double >blocked;
    blocked.resize(blockedDoFNum);
    unsigned s = 0;
    unsigned i = 0;
    for ( vector<Function*>::const_iterator f = dirichF.begin(); f!=dirichF.end(); ++f, i++) {
        if ( *f ) {
            blocked [ s ] = (*f)->giveY(t) * multipliers [ i ];
            s++;
        }
    }
    return blocked;
}

//////////////////////////////////////////////////////////
vector< double >BoundaryCondition :: giveLoadedDoFValues(double t) const {
    vector< double >loaded;
    loaded.resize(loadedDoFNum);
    unsigned s = 0;
    unsigned i = 0;
    for ( vector<Function*>::const_iterator f = neumannF.begin(); f!=neumannF.end(); ++f, i++) {
        if ( *f ) {
            loaded [ s ] = (*f)->giveY(t) * multipliers [ i ];
            s++;
        }
    }
    return loaded;
}


//////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////
// CONTAINER FOR BOUNDARY CONDITIONS
BCContainer :: ~BCContainer() {
    for ( vector< BoundaryCondition * > :: iterator bc = BC.begin(); bc != BC.end(); ++bc ) {
        delete * bc;
    }
}

//////////////////////////////////////////////////////////
void BCContainer :: readFromFile(const string filename, NodeContainer *nodes) {
    size_t origsize = BC.size();
    string line, aux;
    unsigned intnum, nDoFs;
    vector< int >dirichBC, neumannBC;
    Node *node;
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
            iss >> intnum;
            node = nodes->giveNode(intnum);
            nDoFs = node->giveNumberOfDoFs();
            dirichBC.resize(nDoFs);
            neumannBC.resize(nDoFs);
            for ( unsigned i = 0; i < nDoFs; i++ ) {
                iss >> dirichBC [ i ];
            }
            for ( unsigned i = 0; i < nDoFs; i++ ) {
                iss >> neumannBC [ i ];
            }
            for ( unsigned i = 0; i < nDoFs; i++ ) {
                if ( neumannBC [ i ] >= 0 && dirichBC [ i ] >= 0 ) {
                    cerr << "Error: Dirichlet and Neumann boundary conditions assigned simulatneuosly" << endl;
                    cerr << line << endl;
                    exit(EXIT_FAILURE);
                }
            }
            ;
            BoundaryCondition *newBC = new BoundaryCondition(node, dirichBC, neumannBC);
            BC.push_back(newBC);
        }
        inputfile.close();
        cout << "Input file '" <<  filename << "' succesfully loaded; " << BC.size() - origsize << " boundary conditions found" << endl;
    } else {
        cerr << "Error: unable to open input file '" <<  filename <<  "'" << endl;
        exit(EXIT_FAILURE);
    }
}

//////////////////////////////////////////////////////////
void BCContainer :: init() {
    dirichDoFs.resize(0);
    neumannDoFs.resize(0);

    for ( vector< BoundaryCondition * > :: iterator bc = BC.begin(); bc != BC.end(); ++bc ) {
        ( * bc )->init(functions);
    }
}

//////////////////////////////////////////////////////////
void BCContainer :: calculateDoFfields() {
    vector< unsigned >help;
    vector< double >help2;
    for ( vector< BoundaryCondition * > :: iterator bc = BC.begin(); bc != BC.end(); ++bc ) {
        help = ( * bc )->giveBlockedDoFs();
        dirichDoFs.insert(dirichDoFs.end(), help.begin(), help.end() );
        help = ( * bc )->giveLoadedDoFs();
        neumannDoFs.insert(neumannDoFs.end(), help.begin(), help.end() );
    }

    // NOTE know which fns are actually used //JE WHY? What is this information for? Nobody cares. // JK to prevent restricting time step in extreme points of unused fns, it is probably not necessary and if anyone does not comment (or remove fn from fn file) fn that is not used it is his problem, can be removed then
    // for (auto const &f_id : dirichF ){
    //   if ( !functions->isActive(f_id) ){
    //     functions->setActive(f_id);
    //   }
    // }
    // for (auto const &f_id : neumannF ){
    //   if ( !functions->isActive(f_id) ){
    //     functions->setActive(f_id);
    //   }
    // }
}

//////////////////////////////////////////////////////////
vector< double >BCContainer :: giveBlockedDoFValues(double t) const {
    vector< double >blocked(dirichDoFs.size() );
    unsigned i, s = 0;
    vector < double > b;
    for (auto &bc: BC){
        b = bc->giveBlockedDoFValues(t);
        for(i=0; i<b.size(); i++, s++)
        blocked [ s ] = b[i];
    }
    return blocked;
}

//////////////////////////////////////////////////////////
vector< double >BCContainer :: giveLoadedDoFValues(double t) const {
    vector< double >loaded(neumannDoFs.size() );
    unsigned i, s = 0;
    vector < double > b;
    for (auto &bc: BC){
        b = bc->giveLoadedDoFValues(t);
        for(i=0; i<b.size(); i++, s++)
        loaded [ s ] = b[i];
    }
    return loaded;
}
