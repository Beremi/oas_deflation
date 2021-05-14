#ifndef _EXPORTER_C_H
#define _EXPORTER_C_H

#include "globals.h"
#include "node_container.h"
#include "element_container.h"
#include <vector>
#include <iostream>
#include <fstream>


//////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////
// BASIC DATA EXPORTER - master class
class DataExporter
{
private:
public:
    DataExporter(unsigned dimension) { dim = dimension; precision = 6; };
    virtual ~DataExporter() {};
    virtual void readFromLine(istringstream &iss);
    virtual bool doExportNow(const double &time);
    virtual void exportData(unsigned step, const Vector &DoFs, const Vector &reactions, fs :: path resultDir) const = 0;
    virtual void giveFileName(unsigned step, char *buffer) const;
    string giveFileName() const { return filename; };
    virtual void init() {};
    void appendToName(string app) { filename = filename + app; };
protected:
    unsigned dim;
    string filename;
    unsigned precision;
    vector< string >codes;
    double time_each, time_last;
};

//////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////
// EXPORT FROM NODES TO TXT
class TXTNodalExporter : public DataExporter
{
private:
    NodeContainer *nodes;
public:
    TXTNodalExporter(NodeContainer *n, unsigned dimension) : DataExporter(dimension) { nodes = n; };
    ~TXTNodalExporter() {};
    void readFromLine(istringstream &iss);
    virtual void exportData(unsigned step, const Vector &DoFs, const Vector &reactions, fs :: path resultDir) const;
protected:
};

//////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////
// EXPORT FROM ELEMENTS TO TXT
class TXTElementExporter : public DataExporter
{
private:
    ElementContainer *elems;
public:
    TXTElementExporter(ElementContainer *e, unsigned dimension) : DataExporter(dimension) { elems = e; };
    ~TXTElementExporter() {};
    void readFromLine(istringstream &iss);
    virtual void exportData(unsigned step, const Vector &DoFs, const Vector &reactions, fs :: path resultDir) const { ( void ) step; ( void ) DoFs; ( void ) reactions; ( void ) resultDir; };
protected:
};

//////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////
// EXPORT FROM GAUSS POINTS TO TXT
class TXTGaussPointExporter : public DataExporter
{
private:
    ElementContainer *elems;
public:
    TXTGaussPointExporter(ElementContainer *e, unsigned dimension) : DataExporter(dimension) { elems = e; };
    ~TXTGaussPointExporter() {};
    void readFromLine(istringstream &iss);
    virtual void exportData(unsigned step, const Vector &DoFs, const Vector &reactions, fs :: path resultDir) const;
protected:
};

//////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////
// GAUGE EXPORTER
class Gauge : public DataExporter
{
protected:
    string name;
    double multiplier;
public:
    Gauge(unsigned dimension) : DataExporter(dimension) {};
    ~Gauge() {};
    virtual void giveFileName(unsigned step, char *buffer) const;
    string giveName() { return name; };
};

//////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////
// EXPORT OF FORCES
class ForceGauge : public Gauge
{
protected:
    NodeContainer *nodes;
    vector< unsigned >DoFs;
    vector< unsigned >n;
public:
    ForceGauge(NodeContainer *nc, unsigned dimension) : Gauge(dimension) { nodes = nc; multiplier = 1; };
    ForceGauge(string &f, string &gname, string &c, vector< unsigned > &nn, NodeContainer *nc, double m, unsigned dimension);
    ~ForceGauge() {};
    void readFromLine(istringstream &iss);
    virtual void exportData(unsigned step, const Vector &DoFs, const Vector &reactions, fs :: path resultDir) const;
    virtual void init();
};


//////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////
// EXPORT OF DoFs
class DoFGauge : public ForceGauge
{
public:
    DoFGauge(NodeContainer *nc, unsigned dimension) : ForceGauge(nc, dimension) {};
    DoFGauge(string &f, string &gname, string &c, vector< unsigned > &nn, NodeContainer *nc, double m, unsigned dimension);
    ~DoFGauge() {};
    virtual void exportData(unsigned step, const Vector &DoFs, const Vector &reactions, fs :: path resultDir) const;
    virtual void init();
};

//////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////
// EXPORT OF IP VALUES
class IPGauge : public Gauge
{
protected:
    vector< unsigned > elems;
    vector< unsigned > ipnums;
    ElementContainer * elemcont;
public:
    IPGauge(ElementContainer *ec, unsigned dimension) : Gauge(dimension) { elemcont = ec; multiplier = 1; };
    ~IPGauge() {};
    void readFromLine(istringstream &iss);
    virtual void exportData(unsigned step, const Vector &DoFs, const Vector &reactions, fs :: path resultDir) const;
    virtual void init();
};



//////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////
// EXPORT OF DISPLACEMENTS DIFFERENCES
class DisplacementGauge : public Gauge
{
private:
    NodeContainer *nodes;
    ElementContainer *elems;
    Node *nodeA, *nodeB;
    Point pointA, pointB;
public:
    DisplacementGauge(NodeContainer *n, ElementContainer *e, unsigned dimension) : Gauge(dimension) { nodes = n; elems = e;  multiplier = 1; };
    ~DisplacementGauge() {};
    void readFromLine(istringstream &iss);
    virtual void exportData(unsigned step, const Vector &DoFs, const Vector &reactions, fs :: path resultDir) const;
    virtual void init();
protected:
};

//////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////
// EXPORT OF SUM OF VALUES
class StructuralExporter : public Gauge
{
private:
    NodeContainer *nodes;
    ElementContainer *elems;
    double calcValue() const;
public:
    StructuralExporter(NodeContainer *n, ElementContainer *e, unsigned dimension) : Gauge(dimension) { nodes = n; elems = e;  multiplier = 1; };
    ~StructuralExporter() {};
    void readFromLine(istringstream &iss);
    virtual void exportData(unsigned step, const Vector &DoFs, const Vector &reactions, fs :: path resultDir) const;
    virtual void init();
protected:
};

//////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////
// CONTAINER FOR EXPORTERS
class ExporterContainer
{
private:
    fs :: path resultDir;
    vector< DataExporter * >exporters;
    vector< DataExporter * >unique_file_exporters;
public:
    ExporterContainer() {};
    ~ExporterContainer();
    void readFromFile(const string filename, NodeContainer *n, ElementContainer *e, unsigned dimension);
    void exportData(unsigned step, double time, const Vector &DoFs, const Vector &reactions, const bool &exportAll) const;
    void addExporter(DataExporter *de) { exporters.push_back(de); };
    size_t giveSize() { return exporters.size(); }
    void init(const bool &initial = true);
    void setResultDirectory(fs :: path directory) { resultDir = directory; }
    fs :: path giveDirectoryPath() { return resultDir; }
    void appendToAllNames(string app);
protected:
};

void ExportAllElementsNodalStress(std :: vector< Matrix > &stress, const Vector &DoFs, const Vector &reactions, const NodeContainer *nodes, const ElementContainer *elems, const unsigned &dim);

void saveNodes(const NodeContainer &nodes, const std :: vector< std :: string > &NodeTypes, fs :: path resultDir);

void saveElems(const ElementContainer &elems, const std :: vector< std :: string > &ElemTypes, fs :: path resultDir);

#endif /* _EXPORTER_C_H */
