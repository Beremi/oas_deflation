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

#include "globals.h"
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
    ~Model(){ delete solver;};
    void readFromFile(const string filename);
    void init();
    void solve();

    ElementContainer* giveElements() {return &elems;};
    NodeContainer* giveNodes() {return &nodes;};
    Solver* giveSolver() {return solver;};
    FunctionContainer* giveFunctions() {return &funcs;};
    BCContainer* giveBC() {return &bconds;};

protected:
    fs :: path baseDir;

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
