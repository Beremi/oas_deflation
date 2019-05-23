#ifndef _EXPORTER_C_H
#define _EXPORTER_C_H

#include "node_container.h"
#include "element_container.h"
#include <vector>
#include <iostream>
#include <fstream>

//////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////
// BASIC DATA EXPORTER - master class
class DataExporter {
private:
   
public:
    DataExporter(){};
    ~DataExporter(){};
    virtual void readFromLine(istringstream &iss){};
    virtual void exportData(int step, const Vector &DoFs)const{};
    void giveFileName(int step, char* buffer) const;           
protected:
    string filename;
    vector<string> codes;
};

//////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////
// EXPORT FROM NODES TO TXT
class TXTNodalExporter: public DataExporter  {
private:
    NodeContainer *nodes;
public:
    TXTNodalExporter(NodeContainer *n){nodes=n;};
    ~TXTNodalExporter();
    void readFromLine(istringstream &iss);
    virtual void exportData(int step, const Vector &DoFs) const;
protected:

};

//////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////
// EXPORT FROM ELEMENTS TO TXT
class TXTElementExporter: public DataExporter {
private:
    ElementContainer *elems;   
public:
    TXTElementExporter(ElementContainer *e){elems=e;};
    ~TXTElementExporter(){};
    void readFromLine(istringstream &iss);
    virtual void exportData(int step, Vector &DoFs) const{};
protected:

};

//////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////
// CONTAINER FOR EXPORTERS
class ExporterContainer {
private:
    vector <DataExporter*> exporters;
public:
    ExporterContainer(){};
    ~ExporterContainer();
    void readFromFile(const string filename,NodeContainer *n, ElementContainer *e);
    void exportData(int step,const Vector &DoFs)const;
protected:

};



#endif /* _EXPORTER_C_H */
