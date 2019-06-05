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
    virtual void readFromLine(istringstream &iss, int dimension){};
    virtual void exportData(int step, const Vector &DoFs, const Vector &reactions)const{};
    virtual void giveFileName(int step, char* buffer) const; 
    string giveFileName() const{return filename;};
    virtual void init(){};             
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
    void readFromLine(istringstream &iss, int dimension);
    virtual void exportData(int step, const Vector &DoFs, const Vector &reactions) const;
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
    void readFromLine(istringstream &iss, int dimension);
    virtual void exportData(int step, const Vector &DoFs, const Vector &reactions) const{};
protected:

};

//////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////
// EXPORT FROM GAUSS POINTS TO TXT
class TXTGaussPointExporter: public DataExporter {
private:
    ElementContainer *elems;   
public:
    TXTGaussPointExporter(ElementContainer *e){elems=e;};
    ~TXTGaussPointExporter(){};
    void readFromLine(istringstream &iss, int dimension);
    virtual void exportData(int step, const Vector &DoFs, const Vector &reactions) const;
protected:

};


//////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////
// GAUGE EXPORTER
class Gauge: public DataExporter {
protected:
    string name;
public:
    Gauge(){};
    ~Gauge(){};
    virtual void giveFileName(int step, char* buffer) const; 
    string giveName(){return name;};
};

//////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////
// EXPORT OF FORCES
class ForceGauge: public Gauge {
private:
    NodeContainer *nodes;  
    vector<unsigned> DoFs;
    vector<unsigned> n;
public:
    ForceGauge(NodeContainer *n){nodes=n;};
    ~ForceGauge(){};
    void readFromLine(istringstream &iss, int dimension);
    virtual void exportData(int step, const Vector &DoFs, const Vector &reactions) const;
    virtual void init();
protected:

};

//////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////
// EXPORT OF DISPLACEMENTS
class DisplacementGauge: public Gauge {
private:
    NodeContainer *nodes;
    ElementContainer *elems;
    Node *nodeA, *nodeB;
    Point pointA, pointB;
public:
    DisplacementGauge(NodeContainer *n, ElementContainer *e){nodes=n; elems=e;};
    ~DisplacementGauge(){};
    void readFromLine(istringstream &iss, int dimension);
    virtual void exportData(int step, const Vector &DoFs, const Vector &reactions) const;
    virtual void init();
protected:
};

//////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////
// CONTAINER FOR EXPORTERS
class ExporterContainer {
private:
    vector <DataExporter*> exporters;
    vector <DataExporter*> unique_file_exporters;
public:
    ExporterContainer(){};
    ~ExporterContainer();
    void readFromFile(const string filename,NodeContainer *n, ElementContainer *e, unsigned dimension);
    void exportData(int step,const Vector &DoFs, const Vector &reactions)const;
    void init();
protected:

};



#endif /* _EXPORTER_C_H */
