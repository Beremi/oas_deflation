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

#include "version.h"
#include "globals.h"
#include "material_container.h"
#include "node_container.h"
#include "element_container.h"
#include "boundary_condition.h"
#include "data_exporter.h"
#include "preprocessing_block.h"
#include "solver.h"
using namespace std;

#define PRINT_TIME true

fs :: path GlobPaths :: INPUT;
fs :: path GlobPaths :: INPUTFILENAME;
fs :: path GlobPaths :: BASEDIR;
fs :: path GlobPaths :: RESULTDIR;

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


Solver *readMasterFile(const string filename, NodeContainer *nodes, MaterialContainer *matrs, ElementContainer *elems, FunctionContainer *funcs, BCContainer *bconds, ConstraintContainer *constr, ExporterContainer *exporters, PBlockContainer *pblocks) {
    string istr, line;
    int iint, dimension;
    Solver *newsolver = nullptr;
    ifstream inputfile( ( GlobPaths :: BASEDIR / filename ).string() );
    if ( inputfile.is_open() ) {
        while ( getline(inputfile >> std :: ws, line) ) {
            if ( line.empty() ) {
                continue;
            }
            if ( line.at(0) == '#' ) {
                continue;
            }
            istringstream iss(line);
            iss >> std :: skipws >> istr;
            if ( istr.compare("Dimension") == 0 ) {
                iss >> dimension;
            } else if ( istr.compare("NodeFiles") == 0 ) {
                iss >> iint;
                for ( int i = 0; i < iint; i++ ) {
                    iss >> istr;
                    nodes->readFromFile( ( GlobPaths :: BASEDIR / istr ).string(), dimension );
                }
            } else if ( istr.compare("MatFiles") == 0 ) {
                iss >> iint;
                for ( int i = 0; i < iint; i++ ) {
                    iss >> istr;
                    matrs->readFromFile( ( GlobPaths :: BASEDIR / istr ).string() );
                }
            } else if ( istr.compare("ElemFiles") == 0 ) {
                iss >> iint;
                for ( int i = 0; i < iint; i++ ) {
                    iss >> istr;
                    elems->readFromFile( ( GlobPaths :: BASEDIR / istr ).string(), dimension, matrs );
                }
            } else if ( istr.compare("ConstrFiles") == 0 ) {
                // read constraint files
                iss >> iint;
                for ( int i = 0; i < iint; i++ ) {
                    iss >> istr;
                    constr->readFromFile( ( GlobPaths :: BASEDIR / istr ).string(), dimension, nodes );
                }
            } else if ( istr.compare("BCFiles") == 0 ) {
                iss >> iint;
                for ( int i = 0; i < iint; i++ ) {
                    iss >> istr;
                    bconds->readFromFile( ( GlobPaths :: BASEDIR / istr ).string(), nodes );
                }
            } else if ( istr.compare("FunctionFiles") == 0 ) {
                iss >> std :: skipws >> iint;
                for ( int i = 0; i < iint; i++ ) {
                    iss >> std :: skipws >> istr;
                    funcs->readFromFile( ( GlobPaths :: BASEDIR / istr ).string() );
                }
            } else if ( istr.compare("ExporterFiles") == 0 ) {
                iss >> iint;
                for ( int i = 0; i < iint; i++ ) {
                    iss >> istr;
                    exporters->readFromFile( ( GlobPaths :: BASEDIR / istr ).string(), nodes, elems, dimension );
                }
            } else if ( istr.compare("PBlockFiles") == 0 )  {
                iss >> iint;
                for ( int i = 0; i < iint; i++ ) {
                    iss >> istr;
                    pblocks->readFromFile( ( GlobPaths :: BASEDIR / istr ).string(), dimension );
                }
            } else if ( istr.compare("Solver") == 0 ) {
                Solver auxs;
                newsolver = auxs.readFromLine(iss);
                newsolver->setElementContainer(elems);
                newsolver->setNodeContainer(nodes);
                newsolver->setFunctionContainer(funcs);
            }
        }
        inputfile.close();
        cout << "Master input file '" <<  filename << "' succesfully loaded" << endl;
        return newsolver;
    } else {
        cerr << "Error: unable to open input file '" <<  filename <<  "'" << endl;
        exit(EXIT_FAILURE);
    }
}

//////////////////////////////////////////////////////////////////////////////////////

string version_info() {
    string s = "This code has been built from version " + GIT_HASH + " : " + TIME_STRING + "\n";
    s += "OS name: " + OS_NAME + "\n";
    s += "OS sub-type: " + OS_RELEASE + "\n";
    s += "OS build ID: " + OS_VERSION + "\n";
    s += "OS platform: " + OS_PLATFORM + "\n";
    s += "Build type: " + BUILD_TYPE;
    return s;
}

//////////////////////////////////////////////////////////////////////////////////////

int main(int argc, char **argv) {
    if ( argc == 1 ) {
        cerr << "Error: no input file specified, please try again with input file" << endl;
        cerr << "Usage: DiscreteModel [path/to/master.inp]" << endl;
        //cerr << std::setfill('=') << setw(50) << "" << endl;
        cerr << string(80, '=') << endl;
        cerr << version_info() << endl;
        exit(EXIT_FAILURE);
    }

    GlobPaths :: INPUT = fs :: absolute(argv [ 1 ]);
    GlobPaths :: INPUTFILENAME = GlobPaths :: INPUT.filename();
    GlobPaths :: BASEDIR = GlobPaths :: INPUT.parent_path();
    GlobPaths :: RESULTDIR = GlobPaths :: BASEDIR / "results";

    /*cout << GlobPaths::INPUTFILE << endl
     *   << GlobPaths::BASEDIR << endl
     *   << GlobPaths::RESULTDIR << endl;*/
    // create directory for results
    fs :: create_directories(GlobPaths :: RESULTDIR);

    ofstream outputfile( ( GlobPaths :: RESULTDIR / "version.txt" ).string() );
    if ( outputfile.is_open() ) {
        outputfile << std :: scientific;
        outputfile << version_info();
        outputfile << endl;
    }
    outputfile.close();

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
    BCContainer bconds(& funcs);
    ConstraintContainer constr;
    NodeContainer nodes(& bconds);
    nodes.setConstraintContainer(& constr);
    MaterialContainer matrs;
    ElementContainer elems;
    elems.setNodeContainer(& nodes);
    ExporterContainer exporters;
    PBlockContainer pblocks;
    pblocks.setContainers(& nodes, & elems, & bconds, & constr, & funcs, & exporters);
    Solver *solver = readMasterFile(GlobPaths :: INPUTFILENAME.string(), & nodes, & matrs, & elems, & funcs, & bconds, & constr, & exporters, & pblocks);

    //initialization
    pblocks.apply();
    bconds.init();
    nodes.init();
    constr.init(& nodes);
    matrs.init();
    elems.init();
    exporters.init();
    solver->init();
    //solution
    while ( !solver->isTerminated() ) {
        start_part = std :: chrono :: system_clock :: now();
        solver->solveStep();
        exporters.exportData(solver->giveStepNumber(), solver->giveTime(), solver->giveDoFValues(), solver->giveNodalForces(), solver->isTerminated() );
        if ( PRINT_TIME ) {
            now = std :: chrono :: system_clock :: now();

            elapsed_seconds = now - start_part;
            std :: cout << "step duration: " << convertTimeToString(elapsed_seconds) << endl;
            cout.flush();
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

    int terminationStatus = solver->giveTerminationStatus();
    std :: cout << "termination status = " << terminationStatus << '\n';
    delete solver;
    return terminationStatus;
    // JK: why is this deleted? does it cause memory leak? after that, the calculation is teminated, it should not cause any problems, without deleting it, main coul just return solver terminationStatus as following (instead of three previous rows)
    // return solver->giveTerminationStatus();
}
