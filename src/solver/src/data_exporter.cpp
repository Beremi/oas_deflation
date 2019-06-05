#include "data_exporter.h"

//////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////
// BASIC DATA EXPORTER - master class
void DataExporter::giveFileName(int step, char* buffer) const {
    sprintf (buffer, "%s_%05d.out", filename.c_str(),step); 
}

//////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////
// EXPORT FROM NODES TO TXT
void TXTNodalExporter::readFromLine(istringstream &iss, int dimension){
    iss >> filename;
    unsigned num;
    iss >> num;
    codes.resize(num);
    for(int i =0; i<num; i++){
        iss >> codes[i];
    }
}

//////////////////////////////////////////////////////////
void TXTNodalExporter::exportData(int step, const Vector &DoFs, const Vector &reactions) const{

    char buffer [100];
    Node *nn;
    double value;
    giveFileName(step, buffer);
    ofstream outputfile (buffer);
    if (outputfile.is_open()) {
        outputfile << std::scientific;
        for(int n = 0; n<nodes->giveSize(); n++){            
            nn = nodes->giveNode(n);
            for(vector<string>::const_iterator c=codes.begin(); c!=codes.end(); ++c){
                value = nn->giveDoFBasedValue(*c,DoFs);
                outputfile << value;                            
                if (c!=codes.end()-1) outputfile << "\t";            
            }   
            outputfile << endl;         
        }
        outputfile.close();
    }       
}

//////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////
// EXPORT FROM ELEMENTS TO TXT
void TXTElementExporter::readFromLine(istringstream &iss, int dimension){
    iss >> filename;
    unsigned num;
    iss >> num;
    codes.resize(num);
    for(int i =0; i<num; i++){
        iss >> codes[i];
    }
}

//////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////
// EXPORT FROM GAUSS POINTS TO TXT
void TXTGaussPointExporter::readFromLine(istringstream &iss, int dimension){
    iss >> filename;
    unsigned num;
    iss >> num;
    codes.resize(num);
    for(int i =0; i<num; i++){
        iss >> codes[i];
    }
}

//////////////////////////////////////////////////////////
void TXTGaussPointExporter::exportData(int step, const Vector &DoFs, const Vector &reactions) const{

    char buffer [100];
    Element *ee;
    double value;
    int nIP;
    giveFileName(step, buffer);
    ofstream outputfile (buffer);
    
    if (outputfile.is_open()) {
        outputfile << std::scientific;
        for(int e = 0; e<elems->giveSize(); e++){            
            ee = elems->giveElement(e);
            nIP = ee->giveIPNum();
            for(unsigned k=0; k<nIP; k++){
                for(vector<string>::const_iterator c=codes.begin(); c!=codes.end(); ++c){
                    value = ee->giveIPValue(*c,k);
                    outputfile << value;                            
                    if (c!=codes.end()-1) outputfile << "\t";            
                }   
                outputfile << endl;         
            }
        }
        outputfile.close();
    }       
}

//////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////
// GAUGE EXPORTERS
void Gauge::giveFileName(int step, char* buffer) const {
    sprintf (buffer, "%s.out", filename.c_str()); 
}


//////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////
// EXPORT OF FORCES
void ForceGauge::readFromLine(istringstream &iss, int dimension){
    iss >> filename;
    iss >> name;
    codes.resize(1);
    iss >> codes[0];
    unsigned num;
    iss >> num;
    n.resize(num);
    for(int i =0; i<num; i++) iss >> n[i];  
}
//////////////////////////////////////////////////////////
void ForceGauge::init(){
    unsigned DoFpos = 0;
    if (codes[0].compare("fx")==0) DoFpos=0;
    else if (codes[0].compare("fy")==0) DoFpos=1;
    else if (codes[0].compare("fz")==0) DoFpos=2;
    else{
        cerr<<"Error: only 'fx', 'fy' or 'fz' can be exported by ForceGauge" << endl;
        exit(0);
    }    

    DoFs.resize(n.size());
    for(int i =0; i<n.size(); i++) DoFs[i] = nodes->giveNode(n[i])->giveStartingDoF()+DoFpos;
}


//////////////////////////////////////////////////////////
void ForceGauge::exportData(int step, const Vector &full_f, const Vector &reactions) const{

    char buffer [100];
    double value = 0;
    giveFileName(step, buffer);
    ofstream outputfile;
    outputfile.open(buffer,ios::app);
    if (outputfile.good()) {
        outputfile << std::scientific;
        for(int i = 0; i<DoFs.size(); i++) value += reactions[DoFs[i]];
        outputfile <<  "\t" << value;
    }       
    outputfile.close();
}

//////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////
// EXPORT OF DISPLACEMENTS
void DisplacementGauge::readFromLine(istringstream &iss, int dim){
    iss >> filename;
    iss >> name;
    codes.resize(1);
    iss >> codes[0];
    double x,y,z;
    if (dim==2){
        iss >> x >> y;
        pointA = Point(x,y);
        iss >> x >> y;
        pointB = Point(x,y);
    }else if (dim==3){
        iss >> x >> y >> z;
        pointA = Point(x,y,z);
        iss >> x >> y >> z;
        pointB = Point(x,y,z);
    }
}
//////////////////////////////////////////////////////////
void DisplacementGauge::init(){
    //find closest point
    //TODO: better to use direct displacements at given points A and B
    nodeA = nodes->findClosestMechanicalNode(pointA);
    nodeB = nodes->findClosestMechanicalNode(pointB);    
}

//////////////////////////////////////////////////////////
void DisplacementGauge::exportData(int step, const Vector &DoFs, const Vector &reactions) const{

    char buffer [100];
    double value = 0;
    giveFileName(step, buffer);
    ofstream outputfile;
    outputfile.open(buffer,ios::app);
    if (outputfile.good()) {
        outputfile << std::scientific;
        value = nodeB->giveDoFBasedValue(codes[0],DoFs)-nodeA->giveDoFBasedValue(codes[0],DoFs); 
        outputfile << "\t" << value;
    }
    outputfile.close();
}

//////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////
// CONTAINER FOR EXPORTERS
ExporterContainer::~ExporterContainer(){
    for(vector<DataExporter*>::iterator e=exporters.begin(); e!=exporters.end(); ++e) delete *e;    
}

//////////////////////////////////////////////////////////
void ExporterContainer::readFromFile(const string filename,NodeContainer *n, ElementContainer *e, unsigned dimension){
    int origsize = exporters.size();
    string line, exptype;
    ifstream inputfile (filename.c_str());
    if (inputfile.is_open()) {
        while (getline(inputfile, line)){
            if (line.at(0) == '#') continue;
            istringstream iss(line); 
            iss >> exptype;            
            if (not exptype.rfind("#", 0) == 0){
                if (exptype.compare("TXTNodalExporter") == 0){
                    TXTNodalExporter* newexp = new TXTNodalExporter(n);
                    newexp->readFromLine(iss, dimension);     
                    exporters.push_back(newexp);
                }else if (exptype.compare("TXTElementExporter") == 0){
                    TXTElementExporter* newexp = new TXTElementExporter(e);
                    newexp->readFromLine(iss, dimension);     
                    exporters.push_back(newexp);
                }else if (exptype.compare("ForceGauge") == 0){
                    ForceGauge* newexp = new ForceGauge(n);
                    newexp->readFromLine(iss, dimension);     
                    exporters.push_back(newexp);
                }else if (exptype.compare("DisplacementGauge") == 0){
                    DisplacementGauge* newexp = new DisplacementGauge(n,e);
                    newexp->readFromLine(iss, dimension);     
                    exporters.push_back(newexp);
                }else if (exptype.compare("TXTGaussPointExporter") == 0){
                    TXTGaussPointExporter* newexp = new TXTGaussPointExporter(e);
                    newexp->readFromLine(iss, dimension);     
                    exporters.push_back(newexp);
                }else{


                    cerr << "Error: Data exporter '" <<  exptype <<  "' is not implemented yet." << endl;
                    exit(0);
                }
            }                                
        }
        inputfile.close();   
        cout << "Input file '" <<  filename << "' succesfully loaded; "<< exporters.size()-origsize << " exporters found" << endl;
    } else {
        cerr << "Error ExporterContainer: unable to open input file '" <<  filename <<  "'" << endl;
        exit(0);
    }
}  

//////////////////////////////////////////////////////////
void ExporterContainer::init(){
    bool newname;
    for(vector<DataExporter*>::const_iterator d=exporters.begin(); d!=exporters.end(); ++d){
        (*d)->init();
        Gauge *g = dynamic_cast<Gauge*> (*d);
        if(g){
            newname = true;
            for(vector<DataExporter*>::const_iterator unique=unique_file_exporters.begin(); unique!=unique_file_exporters.end(); ++unique){
                if((*unique)->giveFileName().compare((*d)->giveFileName())==0) newname=false;
            }
            if(newname) unique_file_exporters.push_back(*d);
        }
    }

    //gauge files header
    char buffer [100];
    for(vector<DataExporter*>::const_iterator unique=unique_file_exporters.begin(); unique!=unique_file_exporters.end(); ++unique){
        (*unique)->giveFileName(0,buffer);
        ofstream outputfile;
        outputfile.open(buffer);
        if (outputfile.good()) outputfile << "#step";
        outputfile.close();
    }  

    for(vector<DataExporter*>::const_iterator d=exporters.begin(); d!=exporters.end(); ++d){
        Gauge *g = dynamic_cast<Gauge*> (*d);
        if(g){
            (*d)->giveFileName(0,buffer);
            ofstream outputfile;
            outputfile.open(buffer,ios::app);
            if (outputfile.good()) outputfile << "\t" << g->giveName();
            outputfile.close();
        }
    }
    
    for(vector<DataExporter*>::const_iterator unique=unique_file_exporters.begin(); unique!=unique_file_exporters.end(); ++unique){
        (*unique)->giveFileName(0,buffer);
        ofstream outputfile;
        outputfile.open(buffer,ios::app);
        if (outputfile.good()) outputfile << endl;
        outputfile.close();
    }  
};

//////////////////////////////////////////////////////////
void ExporterContainer::exportData(int step, const Vector &DoFs, const Vector &reactions)const{
    //add step number to gauge exporter files
    char buffer [100];
    for(vector<DataExporter*>::const_iterator unique=unique_file_exporters.begin(); unique!=unique_file_exporters.end(); ++unique){
        (*unique)->giveFileName(0,buffer);
        ofstream outputfile;
        outputfile.open(buffer,ios::app);
        if (outputfile.good()) outputfile << step;
        outputfile.close();
    }    

    //export
    for(vector<DataExporter*>::const_iterator d=exporters.begin(); d!=exporters.end(); ++d) (*d)->exportData(step, DoFs, reactions);

    //add end line to gauge exporter files
    for(vector<DataExporter*>::const_iterator unique=unique_file_exporters.begin(); unique!=unique_file_exporters.end(); ++unique){
        (*unique)->giveFileName(0,buffer);
        ofstream outputfile;
        outputfile.open(buffer,ios::app);
        if (outputfile.good()) outputfile << endl;
        outputfile.close();
    }   
};
