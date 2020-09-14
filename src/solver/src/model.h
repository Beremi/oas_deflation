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

using namespace std;

/*
//////////////////////////////////////////////////////////////////////////////////////
string convertTimeToString(std :: chrono :: duration< double >time_interval) {
    int hours = time_interval.count() / 3600;
    int minutes = time_interval.count() / 60 - hours * 60;
    int seconds = time_interval.count() - hours * 3600 - minutes * 60;
    int miliseconds = ( time_interval.count() - hours * 3600 - minutes * 60 - seconds ) * 100;
    stringstream ss;
    ss << std :: setw(2) << std :: setfill('0') << hours << ":"
       << std :: setw(2) << std :: setfill('0') << minutes << ":"
       << std :: setw(2) << std :: setfill('0') << seconds << "."
       << std :: setw(3) << std :: setfill('0') << miliseconds;
    return ss.str();
}
*/

//////////////////////////////////////////////////////////////////////////////////////

class Model
{
public:
    Model(bool pT);
    ~Model(){ delete solver;};
    void readFromFile(const string filename);
    void init();
    void solve();

    fs :: path baseDir;
protected:
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
