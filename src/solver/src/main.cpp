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
#include "solver.h"
using namespace std;

#define PRINT_TIME true

string convertTimeToString(std :: chrono :: duration< double >time_interval) {
    int hours = time_interval.count() / 3600;
    int minutes = time_interval.count() / 60 - hours * 60;
    int seconds = time_interval.count() - hours * 3600 - minutes * 60;
    stringstream ss;
    ss << hours << " : " << minutes << " : " << seconds;
    return ss.str();
}


Solver *readMasterFile(const string filename, NodeContainer *nodes, MaterialContainer *matrs, ElementContainer *elems, FunctionContainer *funcs, BCContainer *bcconds, ExporterContainer *exporters) {
    string istr, line;
    int iint, dimension;
    Solver *newsolver;
    ifstream inputfile(filename.c_str() );
    if ( inputfile.is_open() ) {
        while ( getline(inputfile, line) ) {
            if ( line.at(0) == '#' ) {
                continue;
            }
            istringstream iss(line);
            iss >> istr;
            if ( istr.compare("Dimension") == 0 ) {
                iss >> dimension;
            } else if ( istr.compare("NodeFiles") == 0 )    {
                iss >> iint;
                for ( int i = 0; i < iint; i++ ) {
                    iss >> istr;
                    nodes->readFromFile(istr, dimension);
                }
            } else if ( istr.compare("MatFiles") == 0 )    {
                iss >> iint;
                for ( int i = 0; i < iint; i++ ) {
                    iss >> istr;
                    matrs->readFromFile(istr);
                }
            } else if ( istr.compare("ElemFiles") == 0 )    {
                iss >> iint;
                for ( int i = 0; i < iint; i++ ) {
                    iss >> istr;
                    elems->readFromFile(istr, dimension, matrs);
                }
            } else if ( istr.compare("BCFiles") == 0 )    {
                iss >> iint;
                for ( int i = 0; i < iint; i++ ) {
                    iss >> istr;
                    bcconds->readFromFile(istr, nodes);
                }
            } else if ( istr.compare("FunctionFiles") == 0 )    {
                iss >> iint;
                for ( int i = 0; i < iint; i++ ) {
                    iss >> istr;
                    funcs->readFromFile(istr);
                }
            } else if ( istr.compare("ExporterFiles") == 0 )    {
                iss >> iint;
                for ( int i = 0; i < iint; i++ ) {
                    iss >> istr;
                    exporters->readFromFile(istr, nodes, elems, dimension);
                }
            } else if ( istr.compare("Solver") == 0 )    {
                Solver auxs;
                newsolver = auxs.readFromLine(iss);
                newsolver->setElementContainer(elems);
                newsolver->setNodeContainer(nodes);
            }
        }
        inputfile.close();
        cout << "Master input file '" <<  filename << "' succesfully loaded" << endl;
        return newsolver;
    } else {
        cerr << "Error: unable to open input file '" <<  filename <<  "'" << endl;
        exit(0);
    }
}


//////////////////////////////////////////////////////////////////////////////////////

int main(int argc, char **argv) {
    if ( argc == 1 ) {
        cerr << "Error: no input file specified, please try again with input file" << endl;
        exit(0);
    }

    auto start = std :: chrono :: system_clock :: now();
    auto now = start;
    auto start_part = start;
    std :: chrono :: duration< double >elapsed_seconds;
    if ( PRINT_TIME ) {
        std :: time_t time_now = std :: chrono :: system_clock :: to_time_t(start);
        string nowstring = ctime(& time_now);
        std :: cout << "######### start of calculation on: " << nowstring.substr(0, nowstring.length() - 1) << " #########" << endl;
    }

    //read all files
    FunctionContainer funcs;
    BCContainer bcconds(& funcs);
    NodeContainer nodes(& bcconds);
    MaterialContainer matrs;
    ElementContainer elems;
    elems.setNodeContainer(& nodes);
    ExporterContainer exporters;
    Solver *solver = readMasterFile(argv [ 1 ], & nodes, & matrs, & elems, & funcs, & bcconds, & exporters);

    //initialization
    bcconds.init();
    nodes.init();
    matrs.init();
    elems.init();
    exporters.init();
    solver->init();


    //solution
    while ( not solver->isTerminated() ) {
        start_part = std :: chrono :: system_clock :: now();
        solver->solveStep();
        exporters.exportData(solver->giveStepNumber(), solver->giveDoFValues(), solver->giveNodalForces() );
        if ( PRINT_TIME ) {
            now = std :: chrono :: system_clock :: now();

            elapsed_seconds = now - start_part;
            std :: cout << "step duration: " << convertTimeToString(elapsed_seconds) << endl;
        }
    }

    if ( PRINT_TIME ) {
        now = std :: chrono :: system_clock :: now();
        elapsed_seconds = now - start;
        std :: time_t time_now = std :: chrono :: system_clock :: to_time_t(now);
        string nowstring = ctime(& time_now);
        std :: cout << "######### end of calculation on: " << nowstring.substr(0, nowstring.length() - 1) << " #########" << endl;
        std :: cout << "######### total duration: " << convertTimeToString(elapsed_seconds) << " #########" << endl;
    }

    delete solver;
    return 0;
}
