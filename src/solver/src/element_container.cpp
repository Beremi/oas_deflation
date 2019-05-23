#include "element_container.h"

//////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////
// ELEMENT CONATINER
ElementContainer::~ElementContainer(){
    for(vector<Element*>::iterator e=elems.begin(); e!=elems.end(); ++e) delete *e;    
}

//////////////////////////////////////////////////////////
void ElementContainer::readFromFile(const string filename, const unsigned ndim, MaterialContainer *matrs){
    int origsize = elems.size();
    string line, elemType;
    ifstream inputfile (filename.c_str());
    if (inputfile.is_open()) {
        while (getline(inputfile, line)){
            if (line.at(0) == '#') continue;
            istringstream iss(line);
            iss >> elemType;            
            if (not elemType.rfind("#", 0) == 0){
                if (elemType.compare("LTCBEAM") == 0){
                    RigidBodyContact* newelem = new RigidBodyContact(ndim);
                    newelem->readFromLine(iss, nodes, matrs);     
                    elems.push_back(newelem);
                }else if (elemType.compare("LTCTRSP") == 0){ 
                    Transp1D* newelem=new Transp1D(ndim); 
                    newelem->readFromLine(iss, nodes, matrs);              
                    elems.push_back(newelem);
                }else{
                    cerr << "Error: element '" <<  elemType <<  "' does not exists" << endl;
                    exit(0);
                }
            }
        }
        inputfile.close();    
        cout << "Input file '" <<  filename << "' succesfully loaded; "<< elems.size()-origsize << " elements found" << endl;
    } else {
        cerr << "Error: unable to open input file '" <<  filename <<  "'" << endl;
        exit(0);
    }
}  

//////////////////////////////////////////////////////////
void ElementContainer::init(){
    for(vector<Element*>::iterator e=elems.begin(); e!=elems.end(); ++e) (*e)->init();
}

//////////////////////////////////////////////////////////
void ElementContainer::prepareSteadyStateMatrices(CoordinateIndexedSparseMatrix &K11, CoordinateIndexedSparseMatrix &K12)const {
    
    map<pair<size_t, size_t>, double> indices11;
    map<pair<size_t, size_t>, double> indices12;

    unsigned nDoFs = nodes->giveTotalNumDoFs();
    unsigned nfreeDoFs = nodes->giveNumFreeDoFs();
    unsigned DoFi, DoFj;
    vector<unsigned> DoFs;
    for(vector<Element*>::const_iterator e=elems.begin(); e!=elems.end(); ++e){
        DoFs = (*e)->giveDoFs();
        for(int i=0; i<DoFs.size(); i++){
            for(int j=i; j<DoFs.size(); j++){                
                DoFi = nodes->giveDoFid(DoFs[i]);
                DoFj = nodes->giveDoFid(DoFs[j]);
                //diagonal
                if (DoFi==DoFj){
                    if (DoFi < nfreeDoFs) indices11.insert(pair<pair<size_t, size_t>, double>(pair<size_t, size_t > (DoFi, DoFi), 0.0));
                }else{
                    //remaining items
                    if (DoFi < nfreeDoFs && DoFj < nfreeDoFs) {
                        indices11.insert(pair<pair<size_t, size_t>, double>(pair<size_t, size_t > (DoFi, DoFj), 0.0));
                        indices11.insert(pair<pair<size_t, size_t>, double>(pair<size_t, size_t > (DoFj, DoFi), 0.0));
                    }else if (DoFi < nfreeDoFs && DoFj >= nfreeDoFs) indices12.insert(pair<pair<size_t, size_t>, double>(pair<size_t, size_t > (DoFi, DoFj-nfreeDoFs), 0.0));
                    else if (DoFj < nfreeDoFs && DoFi >= nfreeDoFs) indices12.insert(pair<pair<size_t, size_t>, double>(pair<size_t, size_t > (DoFj, DoFi-nfreeDoFs), 0.0));
                }
            }
        }
    }

    K11 = CoordinateIndexedSparseMatrix(indices11, nfreeDoFs, nfreeDoFs);
    K12 = CoordinateIndexedSparseMatrix(indices12, nfreeDoFs, nDoFs-nfreeDoFs);
}

//////////////////////////////////////////////////////////
void ElementContainer::updateSteadyStateMatrices(CoordinateIndexedSparseMatrix &K11, CoordinateIndexedSparseMatrix &K12)const {
    
    K11 = K11*0.;
    K12 = K12*0.;

    unsigned nDoFs = nodes->giveTotalNumDoFs();
    unsigned nfreeDoFs = nodes->giveNumFreeDoFs();
    unsigned DoFi, DoFj;
    vector<unsigned> DoFs;
    Matrix K;
    for(vector<Element*>::const_iterator e=elems.begin(); e!=elems.end(); ++e){
        DoFs = (*e)->giveDoFs();
        K = (*e)->giveSteadyStateMatrix();
        for(int i=0; i<DoFs.size(); i++){
            for(int j=i; j<DoFs.size(); j++){    
                DoFi = nodes->giveDoFid(DoFs[i]);
                DoFj = nodes->giveDoFid(DoFs[j]);                
                //diagonal
                if (DoFi==DoFj){
                    if (DoFi < nfreeDoFs) K11[DoFi][DoFi] += K[i][i];
                }else{
                    //remaining items
                    if (DoFi < nfreeDoFs && DoFj < nfreeDoFs) {K11[DoFi][DoFj] += K[i][j]; K11[DoFj][DoFi] += K[j][i];}                    
                    else if (DoFi < nfreeDoFs && DoFj >= nfreeDoFs) K12[DoFi][DoFj-nfreeDoFs] += K[i][j];
                    else if (DoFj < nfreeDoFs && DoFi >= nfreeDoFs) K12[DoFj][DoFi-nfreeDoFs] += K[j][i];
                }
            }
        }
    } 
    cout << "Steady state matrices updated" << endl;
}

//////////////////////////////////////////////////////////
void ElementContainer::addBodyForces(Vector &R, double time) const{
    //here comes distributed load, self weight 
    //TO BE DONE
}

