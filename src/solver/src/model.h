#ifndef _MODEL_H
#define _MODEL_H

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
#include "cross_section.h"
#include "solver.h"
#include "geometry.h"

class Model
{
public:
    Model(bool pT);
    ~Model();
    // JK: different read for reloading the geometry (e.g. after adaptive remesh)
    void readFromFile(const std :: string filename, const bool &initial = true);
    void init(const bool &initial = true);
    void solve();
    void jumpToNextStage();
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
    RegionContainer *giveRegionContainer() { return & regions; };
    CrossSectionContainer *giveCSContainer() { return & crosssects; };

    unsigned giveDimension() const { return ndim; };
    void resetTime() { solver->setTime(0); solver->setStep(0); }

    void print_res_dir() { std :: cout << "resultdir: " << resultDir.string() << '\n'; }
    std :: chrono :: time_point< std :: chrono :: system_clock >giveStartTime() const { return time_of_sim_start; }

    fs :: path baseDir;
    fs :: path resultDir;

    std :: string initialFieldFile, initialTimeDerFieldFile; //files with initial conditions

protected:
    unsigned ndim;
    bool printTime;
    std :: chrono :: time_point< std :: chrono :: system_clock >time_of_sim_start;

    FunctionContainer funcs;
    BCContainer bconds;
    ConstraintContainer constr;
    NodeContainer nodes;
    MaterialContainer matrs;
    ElementContainer elems;
    ExporterContainer exporters;
    PBlockContainer pblocks;
    RegionContainer regions;
    CrossSectionContainer crosssects;
    Solver *solver;
};


#endif /* _MODEL_H */
