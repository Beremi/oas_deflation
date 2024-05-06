#ifndef _DATA_EXPORTER_H
#define _DATA_EXPORTER_H

#include "globals.h"
#include "node_container.h"
#include "element_container.h"
#include <vector>
#include <iostream>
#include <fstream>
//#include <bits/stdc++.h> //lower case string

class Solver; //forward declaration

//////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////
// BASIC DATA EXPORTER - master class
class DataExporter
{
private:
public:
    DataExporter(unsigned dimension) { dim = dimension; precision = 6; multiplier = 1.; };
    virtual ~DataExporter() {};
    virtual void readFromLine(std :: istringstream &iss);
    virtual bool doExportNow(const double &time, const unsigned &step);
    virtual void updateNextTimeToSave(const double &time);
    virtual void updateNextStepToSave(const unsigned &step);
    virtual void exportData(unsigned step, fs :: path resultDir) const = 0;
    virtual void giveFileName(unsigned step, char *buffer) const;
    std :: string giveFileName() const { return filename; };
    virtual void init();
    size_t giveMaxSize(unsigned c) const { return maxsize [ c ]; }
    void appendToName(std :: string app) { filename = filename + app; };
    void setSolverPointer(Solver *s);
    std :: string giveName() { return name; };
protected:
    unsigned dim;
    std :: string name;
    std :: string filename;
    unsigned precision;
    std :: vector< std :: string >codes;
    std :: vector< size_t >maxsize;
    double multiplier;
    Solver *solver;

    std :: vector< double >times_to_save;  // vector to store times for export
    double saveTime_each = std :: numeric_limits< double > :: max(); // time period
    double saveTime_last = 0; // last saved time using time period
    double time_last = 0; // time of the last export
    double next_time_to_save = std :: numeric_limits< double > :: max(); // time to perform the next export
    unsigned saveTimes_idx = 0; // possition in vector of times_to_save

    std :: vector< unsigned >steps_to_save;  // vector to store steps for export
    unsigned saveStep_each = std :: numeric_limits< unsigned > :: max(); // step period
    unsigned saveStep_last = 0; // last saved step using step period
    unsigned step_last = 0; // step of the last export
    unsigned next_step_to_save = std :: numeric_limits< unsigned > :: max(); // step to perform the next export
    unsigned saveSteps_idx = 0; // possition in vector of steps_to_save
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
    TXTNodalExporter(NodeContainer *n, ElementContainer *e, unsigned dimension) : DataExporter(dimension) { nodes = n; elems = e;  name = "TXTNodalExporter";};
    ~TXTNodalExporter() {};
    virtual void init();
    void readFromLine(std :: istringstream &iss);
    virtual void exportData(unsigned step, fs :: path resultDir) const;
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
    TXTElementExporter(ElementContainer *e, unsigned dimension) : DataExporter(dimension) { elems = e; name = "TXTElementExporter";};
    ~TXTElementExporter() {};
    virtual void init();
    void readFromLine(std :: istringstream &iss);
    virtual void exportData(unsigned step, fs :: path resultDir) const;
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
    TXTIntegrationPointExporter(ElementContainer *e, unsigned dimension) : DataExporter(dimension) { elems = e;   name = "TXTIntegrationPointExporter";};
    ~TXTIntegrationPointExporter() {};
    virtual void init();
    void readFromLine(std :: istringstream &iss);
    virtual void exportData(unsigned step, fs :: path resultDir) const;
protected:
};

//////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////
// GAUGE EXPORTER
class Gauge : public DataExporter
{
protected:
    std :: string gname;
public:
    Gauge(unsigned dimension) : DataExporter(dimension) {   name = "Gauge"; };
    ~Gauge() {};
    virtual void giveFileName(unsigned step, char *buffer) const;
    std :: string giveGaugeName() const { return gname;}
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
    ForceGauge(NodeContainer *nc, unsigned dimension) : Gauge(dimension) { nodes = nc;   name = "ForceGauge";};
    ForceGauge(std :: string &f, std :: string &gaugename, std :: string &c, std :: vector< unsigned > &nn, NodeContainer *nc, double m, unsigned dimension);
    ~ForceGauge() {};
    void readFromLine(std :: istringstream &iss);
    virtual void exportData(unsigned step, fs :: path resultDir) const;
    virtual void init();
};


//////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////
// EXPORT OF DoFs
class DoFGauge : public ForceGauge
{
public:
    DoFGauge(NodeContainer *nc, unsigned dimension) : ForceGauge(nc, dimension) {  name = "DoFGauge";};
    DoFGauge(std :: string &f, std :: string &gaugename, std :: string &c, std :: vector< unsigned > &nn, NodeContainer *nc, double m, unsigned dimension);
    ~DoFGauge() {};
    virtual void exportData(unsigned step, fs :: path resultDir) const;
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
    IntegrationPointGauge(ElementContainer *ec, unsigned dimension) : Gauge(dimension) { elemcont = ec; multiplier = 1;   name = "IntegrationPointGauge";};
    ~IntegrationPointGauge() {};
    void readFromLine(std :: istringstream &iss);
    virtual void exportData(unsigned step, fs :: path resultDir) const;
    virtual void init();
};

//////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////
// EXPORT OF SUMMATIONS AND AVERAGES FROM ELEMENT CONTAINER
class ElementContainerGauge : public Gauge
{
protected:
    ElementContainer *elemcont;
public:
    ElementContainerGauge(ElementContainer *ec, unsigned dimension) : Gauge(dimension) { elemcont = ec; multiplier = 1;   name = "ElementContainerGauge";};
    ~ElementContainerGauge() {};
    void readFromLine(std :: istringstream &iss);
    virtual void exportData(unsigned step, fs :: path resultDir) const;
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
    DisplacementGauge(NodeContainer *n, ElementContainer *e, unsigned dimension) : Gauge(dimension) { nodes = n; elems = e;  multiplier = 1; elemA = nullptr; elemB = nullptr;   name = "DisplacementGauge";};
    ~DisplacementGauge() {};
    void readFromLine(std :: istringstream &iss);
    virtual void exportData(unsigned step, fs :: path resultDir) const;
    virtual void init();
protected:
};

//////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////
// EXPORT OF SOLVER VALUES
class SolverGauge : public Gauge
{
private:
public:
    SolverGauge(unsigned dimension) : Gauge(dimension) {  name = "SolverGauge";};
    ~SolverGauge() {};
    void readFromLine(std :: istringstream &iss);
    virtual void exportData(unsigned step, fs :: path resultDir) const;
    virtual void init();
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
    void exportData(unsigned step, double time, const bool &exportAll) const;
    void addExporter(DataExporter *de) { exporters.push_back(de); };
    size_t giveSize() { return exporters.size(); }
    void init(const bool &initial = true);
    void clear();
    void setResultDirectory(fs :: path directory) { resultDir = directory; }
    fs :: path giveDirectoryPath() { return resultDir; }
    void appendToAllNames(std :: string app);
    void setSolver(Solver *s);
    void updateAllTimeAndStepToSave(unsigned step, double time) const;
protected:
};

void ExportAllElementsNodalStress(std :: vector< Matrix > &stress, const Vector &DoFs, const Vector &reactions, const NodeContainer *nodes, const ElementContainer *elems, const unsigned &dim);

void saveNodes(const NodeContainer &nodes, const std :: vector< std :: string > &NodeTypes, fs :: path resultDir);

void saveElems(const ElementContainer &elems, const std :: vector< std :: string > &ElemTypes, fs :: path resultDir);

#endif /* _DATA_EXPORTER_H */
