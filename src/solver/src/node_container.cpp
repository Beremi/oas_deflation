#include "node_container.h"

//////////////////////////////////////////////////////////
NodeContainer::~NodeContainer(){
    for(vector<Node*>::iterator n=nodes.begin(); n!=nodes.end(); ++n) delete *n;
}

//////////////////////////////////////////////////////////
void NodeContainer::readFromFile(const string filename, const int dim){
    int origsize = nodes.size();
    string line, nodeType;
    ifstream inputfile (filename.c_str());
    if (inputfile.is_open()) {
        while (getline(inputfile, line)){
            if (line.at(0) == '#') continue;
            istringstream iss(line); 
            iss >> nodeType;            
            if (not nodeType.rfind("#", 0) == 0){
                if (nodeType.compare("TrsprtNode") == 0){
                    TrsNode* newnode = new TrsNode(dim);
                    newnode->readFromLine(iss, dim);     
                    nodes.push_back(newnode);
                }else if (nodeType.compare("Particle") == 0){ 
                    Particle* newnode=new Particle(dim); 
                    newnode->readFromLine(iss, dim);              
                    nodes.push_back(newnode);
                }else if (nodeType.compare("AuxNode") == 0){ 
                    AuxNode* newnode=new AuxNode(dim); 
                    newnode->readFromLine(iss, dim);              
                    nodes.push_back(newnode);
                }else{
                    cerr << "Error: node type '" <<  nodeType <<  "' does not exists" << endl;
                    exit(0);
                }
            }      
        }   
        inputfile.close();  
        cout << "Input file '" <<  filename << "' succesfully loaded; "<< nodes.size()-origsize << " nodes found" << endl;
    } else {
        cerr << "Error: unable to open input file '" <<  filename <<  "'" << endl;
        exit(0);
    }
}  

//////////////////////////////////////////////////////////
void NodeContainer::init(){
    establishDoFArray();
}

//////////////////////////////////////////////////////////
void NodeContainer::establishDoFArray(){
    totalDoFs = 0;
    freeDoFs = 0;

    for(vector<Node*>::iterator n = nodes.begin(); n!=nodes.end(); ++n){
        (*n)->setStartingDoF(totalDoFs);
        totalDoFs += (*n)->giveNumberOfDoFs();
        freeDoFs += (*n)->giveNumberOfFreeDoFs();        
    }   
    BC->calculateDoFfields();
    DoFid.resize(totalDoFs);
    vector<unsigned> blocked = BC->giveArrayOfBlockedDoFs();
    vector<unsigned> loaded = BC->giveArrayOfLoadedDoFs();
    blockedDoFid.resize(blocked.size());
    loadedDoFid.resize(loaded.size());

    //sort DoFs, keep track of indices
    vector<pair<int,int> >a;
    for (int i = 0 ;i < blocked.size() ; i++) a.push_back (make_pair (blocked[i],i));
    sort (a.begin(),a.end());
     
    unsigned k=0;
    unsigned id=0;
    for(vector<unsigned>::iterator d=DoFid.begin(); d!=DoFid.end(); ++d,id++){
        if (id==a[k].first){
            *d = freeDoFs+k;
            blockedDoFid[a[k].second]=id;
            k++;
        }
        else *d = id-k;
    }
    for(int i=0; i<loaded.size(); i++)  loadedDoFid[i]=loaded[i];
    cout << "Loaded problem contains " << freeDoFs << " degrees of freedom; additional " << totalDoFs-freeDoFs << " degrees of freedom were prescribed"<< endl;
}

//////////////////////////////////////////////////////////
void NodeContainer::addRHS_nodalLoad(Vector &RHS, double time) const{
    vector<double> loaded = BC->giveLoadedDoFValues(time);      
    for(int k=0; k<loaded.size(); k++) {
        RHS[loadedDoFid[k]] += loaded[k];  
    }
}

//////////////////////////////////////////////////////////
void NodeContainer::updateDirrichletBC(Vector &r, double time) const {
    vector<double> blocked = BC->giveBlockedDoFValues(time);
    for(int k=0; k<blocked.size(); k++) {r[blockedDoFid[k]] = blocked[k];
    }
}


//////////////////////////////////////////////////////////
void NodeContainer::giveFullDoFArray(const Vector &fDoFs, const Vector &bDoFs, Vector &fullDoFs) const{
    for(int i=0; i<totalDoFs; i++){       
        if(DoFid[i]<freeDoFs) fullDoFs[i] = fDoFs[DoFid[i]];
        else fullDoFs[i] = bDoFs[DoFid[i]-freeDoFs];
    }
}

//////////////////////////////////////////////////////////
void NodeContainer::giveFullDoFArray(const Vector &fDoFs, Vector &fullDoFs) const{
    for(int i=0; i<totalDoFs; i++){       
        if(DoFid[i]<freeDoFs) fullDoFs[i] = fDoFs[DoFid[i]];
    }
}

//////////////////////////////////////////////////////////
void NodeContainer::giveReducedDoFArray(const Vector &fullDoFs, Vector &fDoFs) const{
    for(int i=0; i<totalDoFs; i++){       
        if(DoFid[i]<freeDoFs) fDoFs[DoFid[i]] = fullDoFs[i];
        //else bDoFs[DoFid[i]-freeDoFs] = fullDoFs[i];
    }
}

//////////////////////////////////////////////////////////
void NodeContainer::updateExteranlForcesByReactions(const Vector &f_int, const Vector &load, Vector &f_ext) const{
    for(int k=0; k<totalDoFs; k++) {
        f_ext[k] = load[k];
        if(DoFid[k]>=freeDoFs) f_ext[k] += f_int[k];
    }
}

//////////////////////////////////////////////////////////
Node* NodeContainer::findClosestMechanicalNode(Point A) const{
    
    Node* closest;
    double minDist = 1e20;
    double distance2;
    for(vector<Node*>::const_iterator n = nodes.begin(); n!=nodes.end(); ++n){
        if ((*n)->doesMechanics()){
            distance2  = ((*n)->givePoint()-A).sqNorm();
            if(distance2<minDist){
                minDist = distance2;
                closest = (*n);
            }
        }
    }
    return closest;
    
}

