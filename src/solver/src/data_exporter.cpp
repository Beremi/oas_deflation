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
void TXTNodalExporter::readFromLine(istringstream &iss){
    iss >> filename;
    unsigned num;
    iss >> num;
    codes.resize(num);
    for(int i =0; i<num; i++){
        iss >> codes[i];
    }
}

//////////////////////////////////////////////////////////
void TXTNodalExporter::exportData(int step, const Vector &DoFs) const{

    char buffer [100];
    Node *nn;
    double value;
    DataExporter:: giveFileName(step, buffer);
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
void TXTElementExporter::readFromLine(istringstream &iss){
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
// CONTAINER FOR EXPORTERS
ExporterContainer::~ExporterContainer(){
    for(vector<DataExporter*>::iterator e=exporters.begin(); e!=exporters.end(); ++e) delete *e;    
}

//////////////////////////////////////////////////////////
void ExporterContainer::readFromFile(const string filename,NodeContainer *n, ElementContainer *e){
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
                    newexp->readFromLine(iss);     
                    exporters.push_back(newexp);
                }else if (exptype.compare("TXTElementExporter") == 0){
                    TXTElementExporter* newexp = new TXTElementExporter(e);
                    newexp->readFromLine(iss);     
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
void ExporterContainer::exportData(int step, const Vector &DoFs)const{
    for(vector<DataExporter*>::const_iterator d=exporters.begin(); d!=exporters.end(); ++d) (*d)->exportData(step, DoFs);
};
