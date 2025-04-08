#include "model.h"
#include <signal.h>

using namespace std;

volatile sig_atomic_t TERMINATED = 0;

void my_handler(int s) {
    ( void ) s;
    if ( TERMINATED ) {
        exit(EXIT_FAILURE);
    }
    TERMINATED = 1;
}

//////////////////////////////////////////////////////////
Model :: ~Model() {
    if ( solver != nullptr ) {
        delete solver;
    }
};

//////////////////////////////////////////////////////////
Model :: Model(bool pT) {
    printTime = pT;
    nodes.setContainers(& bconds, & constr);
    bconds.setContainers(& funcs);
    elems.setModel(this);
    elems.setContainers(& constr);
    solver = nullptr;
    initialFieldFile = "";
    initialTimeDerFieldFile = "";
    time_of_sim_start = std :: chrono :: system_clock :: now();
}

//////////////////////////////////////////////////////////
void Model :: init(const bool &initial) {     //initialization
    if ( initial ) {
        cout << "initialization of materials" << endl;
        cout.flush();
        matrs.init();
    }
    cout << "initialization of elements" << endl;
    cout.flush();
    elems.init();
    cout << "initialization of preprocessing blocks" << endl;
    cout.flush();
    pblocks.init();
    cout << "initialization of boundary conditions" << endl;
    cout.flush();
    bconds.init(solver->giveTime() );
    cout << "initialization of nodes" << endl;
    cout.flush();
    nodes.init();
    nodes.initSimplices();
    cout << "initialization of constraints" << endl;
    cout.flush();
    constr.init(& nodes, & bconds, solver);
    elems.assignFibersToElems();
    elems.findElementFriends();
    if ( initialFieldFile.compare("") != 0 ) {
        initialFieldFile = ( baseDir / initialFieldFile ).string();
    }
    if ( initialTimeDerFieldFile.compare("") != 0 ) {
        initialTimeDerFieldFile = ( baseDir / initialTimeDerFieldFile ).string();
    }
    cout << "initialization of solver" << endl;
    cout.flush();
    solver->init(initialFieldFile, initialTimeDerFieldFile, initial);
    cout << "initialization of exporters" << endl;
    cout.flush();
    exporters.setResultDirectory(resultDir);
    exporters.setSolver(solver);
    exporters.init(initial);
    bconds.setInitialDoFFields(solver);
    cout << "Model succesfully initialized" << endl;
    // exporters.updateAllTimeAndStepToSave(solver->giveTime(), solver->giveStepNumber());  // needed especially in adaptivity
}

//////////////////////////////////////////////////////////
void Model :: jumpToNextStage() {
    cout << "** updating model to the next computational stage ** " << endl;
    bconds.init(solver->giveTime() + 1e-12);
    bconds.setInitialDoFFields(solver);
    nodes.init();
    constr.init(& nodes, & bconds, solver);
    solver->rebuild();
}

//////////////////////////////////////////////////////////
void Model :: solve() {
    //solution
    signal(SIGINT, my_handler);
    exporters.exportData( solver->giveStepNumber(), -1, solver->giveTime(), solver->isTerminated() );
    while ( !solver->isTerminated() && TERMINATED == 0 ) {
        auto start_part = std :: chrono :: system_clock :: now();
        solver->solveStep();
        exporters.exportData( solver->giveStepNumber(), -1, solver->giveTime(), solver->isTerminated() );
        if ( printTime && solver->showStepTime() ) {
            auto now = std :: chrono :: system_clock :: now();
            auto elapsed_seconds = now - start_part;
            std :: cout << "step duration: " << convertTimeToString(elapsed_seconds) << endl;
            cout.flush();
        }
    }
}

//////////////////////////////////////////////////////////
void Model :: readFromFile(const string filename, const bool &initial) {
    fs :: path fullPath = fs :: absolute(filename);
    baseDir = fullPath.parent_path();
    std :: string result_dir_name = "results";

    string istr, line;
    int iint;
    ifstream inputfile( fullPath.string() );
    if ( inputfile.is_open() ) {
        while ( getline(inputfile >> std :: ws, line) ) {
            if ( line.empty() || ( line.at(0) == '#' ) ) {
                continue;
            }
            istringstream iss(line);
            iss >> std :: skipws >> istr;
            if ( istr.compare("Dimension") == 0 ) {
                iss >> ndim;
            } else if ( istr.compare("NodeFiles") == 0 ) {
                iss >> iint;
                for ( int i = 0; i < iint; i++ ) {
                    iss >> istr;
                    nodes.readFromFile( ( baseDir / istr ).string(), ndim);
                }
            } else if ( initial && istr.compare("MatFiles") == 0 ) {
                iss >> iint;
                for ( int i = 0; i < iint; i++ ) {
                    iss >> istr;
                    matrs.readFromFile( ( baseDir / istr ).string(), ndim);
                }
            } else if ( istr.compare("ElemFiles") == 0 ) {
                iss >> iint;
                for ( int i = 0; i < iint; i++ ) {
                    iss >> istr;
                    elems.readFromFile( ( baseDir / istr ).string(), ndim, & matrs);
                }
            } else if ( istr.compare("MatStatFiles") == 0 ) {
                iss >> iint;
                for ( int i = 0; i < iint; i++ ) {
                    iss >> istr;
                    elems.setFileToLoadStatsFrom( ( baseDir / istr ).string() );
                }
            } else if ( istr.compare("ConstrFiles") == 0 ) {
                // read constraint files
                iss >> iint;
                for ( int i = 0; i < iint; i++ ) {
                    iss >> istr;
                    constr.readFromFile( ( baseDir / istr ).string(), ndim, & nodes);
                }
            } else if ( istr.compare("BCFiles") == 0 ) {
                iss >> iint;
                for ( int i = 0; i < iint; i++ ) {
                    iss >> istr;
                    bconds.readFromFile( ( baseDir / istr ).string(), & nodes, & elems);
                }
            } else if ( initial && istr.compare("FunctionFiles") == 0 ) {  // functions are constant during whole calculation, even in adaptive case
                iss >> std :: skipws >> iint;
                for ( int i = 0; i < iint; i++ ) {
                    iss >> std :: skipws >> istr;
                    funcs.readFromFile( ( baseDir / istr ).string() );
                }
            } else if ( istr.compare("ExporterFiles") == 0 ) {
                iss >> iint;
                for ( int i = 0; i < iint; i++ ) {
                    iss >> istr;
                    exporters.readFromFile( ( baseDir / istr ).string(), & nodes, & elems, & constr, & bconds, ndim);
                }
            } else if ( istr.compare("PBlockFiles") == 0 ) {
                iss >> iint;
                for ( int i = 0; i < iint; i++ ) {
                    iss >> istr;
                    pblocks.readFromFile( ( baseDir / istr ).string(), ndim);
                }
            } else if ( istr.compare("RegionFiles") == 0 ) {
                iss >> iint;
                for ( int i = 0; i < iint; i++ ) {
                    iss >> istr;
                    regions.readFromFile( ( baseDir / istr ).string(), ndim);
                }
            } else if ( initial && istr.compare("Solver") == 0 ) {
                iss >> istr;
                //solver = new Solver;
                //Solver* ptr = solver;
                //solver = solver->readFromFile(( baseDir / istr ).string() );
                //delete ptr;
                solver = Solver().readFromFile( ( baseDir / istr ).string() );
                // QUESTION JK: why is this here and not in the constructor? together with new Solver() ?
                solver->setContainers(& elems, & nodes, & funcs, & bconds, &exporters);
            } else if ( initial && istr.compare("initial_master_field") == 0 ) {
                iss >> initialFieldFile;
            } else if ( initial && istr.compare("initial_master_time_derivative_field") == 0 ) {
                iss >> initialTimeDerFieldFile;
            }  else if ( initial && istr.compare("result_dir") == 0 ) {
                iss >> result_dir_name;
            }
        }
        inputfile.close();
        cout << "Master input file '" <<  fullPath.string() << "' succesfully loaded" << endl;
    } else {
        cerr << "Error: unable to open input file '" <<  fullPath.string() <<  "'" << endl;
        exit(EXIT_FAILURE);
    }

    if ( initial ) {
        resultDir = baseDir / result_dir_name;
    }

    pblocks.setContainers(& nodes, & elems, & bconds, & constr, & funcs, & exporters, & matrs, & regions, solver);
}


void Model :: clear() {
    // TODO JK: check containers for memory leaks
    // initialize new model with clear geometry, only solver remains
    bconds.clear();
    constr.clear();
    nodes.clear();
    elems.clear();
    exporters.clear();
    pblocks.clear();

    // funcs = FunctionContainer();  // functions remain too
    bconds = BCContainer();
    constr = ConstraintContainer();
    nodes = NodeContainer();
    // matrs = MaterialContainer();  // materials remain too
    elems = ElementContainer();
    exporters = ExporterContainer();
    pblocks = PBlockContainer();

    nodes.setContainers(& bconds, & constr);
    bconds.setContainers(& funcs);
    elems.setModel(this);
    pblocks.setContainers(& nodes, & elems, & bconds, & constr, & funcs, & exporters, & matrs, & regions, solver);
    // std :: cout << "step: " << solver->giveStepNumber() << ", time: " << solver->giveTime() << '\n';

    solver->setContainers(& elems, & nodes, & funcs, & bconds, & exporters);
    // std :: cout << "step: " << solver->giveStepNumber() << ", time: " << solver->giveTime() << '\n';
}
