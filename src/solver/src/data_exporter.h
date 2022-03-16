#ifndef _DATA_EXPORTER_H
#define _DATA_EXPORTER_H

#include "globals.h"
#include "node_container.h"
#include "element_container.h"
#include <vector>
#include <iostream>
#include <fstream>

class Solver; //forward declaration

//////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////
// BASIC DATA EXPORTER - master class
class DataExporter
{
private:
public:
    DataExporter(unsigned dimension) { dim = dimension; precision = 6; };
    virtual ~DataExporter() {};
    virtual void readFromLine(std :: istringstream &iss);
    virtual bool doExportNow(const double &time);
    virtual void exportData(unsigned step, const Vector &DoFs, const Vector &reactions, fs :: path resultDir) const = 0;
    virtual void giveFileName(unsigned step, char *buffer) const;
    std :: string giveFileName() const { return filename; };
    virtual void init() {};
    size_t giveMaxSize(unsigned c) const { return maxsize [ c ]; }
    void appendToName(std :: string app) { filename = filename + app; };
protected:
    unsigned dim;
    std :: string filename;
    unsigned precision;
    std :: vector< std :: string >codes;
    std :: vector< size_t >maxsize;
    double time_each, time_last;
};

//////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////
// EXPORT FROM NODES TO TXT
class TXTNodalExporter : public DataExporter
{
private:
    NodeContainer *nodes;
    ElementContainer *elems;
public:
    TXTNodalExporter(NodeContainer *n, ElementContainer *e, unsigned dimension) : DataExporter(dimension) { nodes = n; elems = e; };
    ~TXTNodalExporter() {};
    virtual void init();
    void readFromLine(std :: istringstream &iss);
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
    virtual void init();
    void readFromLine(std :: istringstream &iss);
    virtual void exportData(unsigned step, const Vector &DoFs, const Vector &reactions, fs :: path resultDir) const;
protected:
};

//////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////
// EXPORT FROM GAUSS POINTS TO TXT
class TXTIntegrationPointExporter : public DataExporter
{
private:
    ElementContainer *elems;
public:
    TXTIntegrationPointExporter(ElementContainer *e, unsigned dimension) : DataExporter(dimension) { elems = e; };
    ~TXTIntegrationPointExporter() {};
    virtual void init();
    void readFromLine(std :: istringstream &iss);
    virtual void exportData(unsigned step, const Vector &DoFs, const Vector &reactions, fs :: path resultDir) const;
protected:
};

//////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////
// GAUGE EXPORTER
class Gauge : public DataExporter
{
protected:
    std :: string name;
    double multiplier;
public:
    Gauge(unsigned dimension) : DataExporter(dimension) {};
    ~Gauge() {};
    virtual void giveFileName(unsigned step, char *buffer) const;
    std :: string giveName() { return name; };
};

//////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////
// EXPORT OF FORCES
class ForceGauge : public Gauge
{
protected:
    NodeContainer *nodes;
    std :: vector< unsigned >DoFs;
    std :: vector< unsigned >n;
public:
    ForceGauge(NodeContainer *nc, unsigned dimension) : Gauge(dimension) { nodes = nc; multiplier = 1; };
    ForceGauge(std :: string &f, std :: string &gname, std :: string &c, std :: vector< unsigned > &nn, NodeContainer *nc, double m, unsigned dimension);
    ~ForceGauge() {};
    void readFromLine(std :: istringstream &iss);
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
    DoFGauge(std :: string &f, std :: string &gname, std :: string &c, std :: vector< unsigned > &nn, NodeContainer *nc, double m, unsigned dimension);
    ~DoFGauge() {};
    virtual void exportData(unsigned step, const Vector &DoFs, const Vector &reactions, fs :: path resultDir) const;
    virtual void init();
};

//////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////
// EXPORT OF IP VALUES
class IntegrationPointGauge : public Gauge
{
protected:
    std :: vector< unsigned >elems;
    std :: vector< unsigned >ipnums;
    ElementContainer *elemcont;
public:
    IntegrationPointGauge(ElementContainer *ec, unsigned dimension) : Gauge(dimension) { elemcont = ec; multiplier = 1; };
    ~IntegrationPointGauge() {};
    void readFromLine(std :: istringstream &iss);
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
    Point natCoordsA, natCoordsB;
    Element *elemA, *elemB;
public:
    DisplacementGauge(NodeContainer *n, ElementContainer *e, unsigned dimension) : Gauge(dimension) { nodes = n; elems = e;  multiplier = 1; elemA = nullptr; elemB = nullptr; };
    ~DisplacementGauge() {};
    void readFromLine(std :: istringstream &iss);
    virtual void exportData(unsigned step, const Vector &DoFs, const Vector &reactions, fs :: path resultDir) const;
    virtual void init();
protected:
};

//////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////
// EXPORT OF SOLVER VALUES
class SolverGauge : public Gauge
{
private:
    Solver *solver;
public:
    SolverGauge(unsigned dimension) : Gauge(dimension) {};
    ~SolverGauge() {};
    void readFromLine(std :: istringstream &iss);
    virtual void exportData(unsigned step, const Vector &DoFs, const Vector &reactions, fs :: path resultDir) const;
    virtual void init();
    void setSolverPointer(Solver *s);
protected:
};


/*
 * //////////////////////////////////////////////////////////
 * //////////////////////////////////////////////////////////
 * // EXPORT OF SUM OF VALUES
 * class StructuralExporter : public Gauge
 * {
 * private:
 *  NodeContainer *nodes;
 *  ElementContainer *elems;
 *  double calcValue() const;
 * public:
 *  StructuralExporter(NodeContainer *n, ElementContainer *e, unsigned dimension) : Gauge(dimension) { nodes = n; elems = e;  multiplier = 1; };
 *  ~StructuralExporter() {};
 *  void readFromLine(istringstream &iss);
 *  virtual void exportData(unsigned step, const MyVector &DoFs, const MyVector &reactions, fs :: path resultDir) const;
 *  virtual void init();
 * protected:
 * };
 */

//////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////
// CONTAINER FOR EXPORTERS
class ExporterContainer
{
private:
    fs :: path resultDir;
    std :: vector< DataExporter * >exporters;
    std :: vector< DataExporter * >unique_file_exporters;
public:
    ExporterContainer() {};
    ~ExporterContainer();
    void readFromFile(const std :: string filename, NodeContainer *n, ElementContainer *e, unsigned dimension);
    void exportData(unsigned step, double time, const Vector &DoFs, const Vector &reactions, const bool &exportAll) const;
    void addExporter(DataExporter *de) { exporters.push_back(de); };
    size_t giveSize() { return exporters.size(); }
    void init(const bool &initial = true);
    void clear();
    void setResultDirectory(fs :: path directory) { resultDir = directory; }
    fs :: path giveDirectoryPath() { return resultDir; }
    void appendToAllNames(std :: string app);
    void setSolver(Solver *s);
protected:
};

void ExportAllElementsNodalStress(std :: vector< Matrix > &stress, const Vector &DoFs, const Vector &reactions, const NodeContainer *nodes, const ElementContainer *elems, const unsigned &dim);

void saveNodes(const NodeContainer &nodes, const std :: vector< std :: string > &NodeTypes, fs :: path resultDir);

void saveElems(const ElementContainer &elems, const std :: vector< std :: string > &ElemTypes, fs :: path resultDir);

#endif /* _DATA_EXPORTER_H */
