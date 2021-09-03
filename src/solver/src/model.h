#ifndef MODEL_H
#define MODEL_H

#include <iostream>
#include <fstream>
#include <cmath>
#include <stdio.h>
#include <string>
#include <vector>
#include <valarray>
//#include <io.h>

// time management:
#include <chrono>
#include <ctime>

#include "material_container.h"
#include "node_container.h"
#include "element_container.h"
#include "boundary_condition.h"
#include "data_exporter.h"
#include "preprocessing_block.h"
#include "solver.h"

class Model
{
public:
    Model(bool pT);
    ~Model();
    // JK: different read for reloading the geometry (e.g. after adaptive remesh)
    void readFromFile(const string filename, const bool &initial = true);
    void init(const bool &initial = true);
    void solve();

    void clear();

    ElementContainer *giveElements() { return & elems; };
    NodeContainer *giveNodes() { return & nodes; };
    Solver *giveSolver() { return solver; };
    void setSolver(Solver *s) { solver = s; };
    FunctionContainer *giveFunctions() { return & funcs; };
    BCContainer *giveBC() { return & bconds; };
    ConstraintContainer *giveConstraints() { return & constr; };
    ExporterContainer *giveExporters() { return & exporters; };
    MaterialContainer *giveMaterials() { return & matrs; };
    fs :: path giveResultDirectory() const { return resultDir; };
    PBlockContainer *givePBlockContainer() { return & pblocks; };

    unsigned giveDimension() const { return ndim; };
    void resetTime() { solver->setTime(0); solver->setStep(0); }

    void print_res_dir() { std :: cout << "resultdir: " << resultDir.string() << '\n'; }

    fs :: path baseDir;
    fs :: path resultDir;

    string initialFieldFile, initialTimeDerFieldFile; //files with initial conditions

protected:
    unsigned ndim;
    bool printTime;
    FunctionContainer funcs;
    BCContainer bconds;
    ConstraintContainer constr;
    NodeContainer nodes;
    MaterialContainer matrs;
    ElementContainer elems;
    ExporterContainer exporters;
    PBlockContainer pblocks;
    Solver *solver;
};


#endif
