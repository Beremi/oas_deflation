#include "model.h"

using namespace std;

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
    elems.setContainers(& nodes, & bconds);
    solver = nullptr;
    initialFieldFile = "";
    initialTimeDerFieldFile = "";
    time_of_sim_start = std :: chrono :: system_clock :: now();
}

//////////////////////////////////////////////////////////
void Model :: init(const bool &initial) {     //initialization
    //pblocks.apply(); //moved to reader
    bconds.init();
    nodes.init();
    if ( initial ) {
        matrs.init();
    }
    elems.init();
    nodes.initSimplices();
    constr.init(& nodes, & bconds, solver);
    elems.findElementFriends();
    if ( initialFieldFile.compare("") != 0 ) {
        initialFieldFile = ( baseDir / initialFieldFile ).string();
    }
    if ( initialTimeDerFieldFile.compare("") != 0 ) {
        initialTimeDerFieldFile = ( baseDir / initialTimeDerFieldFile ).string();
    }
    solver->init(initialFieldFile, initialTimeDerFieldFile, initial);
    exporters.init(initial);
}

//////////////////////////////////////////////////////////
void Model :: solve() {
    //solution
    exporters.exportData( solver->giveStepNumber(), solver->giveTime(), solver->giveDoFValues(), solver->giveNodalForces(), solver->isTerminated() );
    while ( !solver->isTerminated() ) {
        auto start_part = std :: chrono :: system_clock :: now();
        solver->solveStep();
        exporters.exportData( solver->giveStepNumber(), solver->giveTime(), solver->giveDoFValues(), solver->giveNodalForces(), solver->isTerminated() );
        if ( printTime ) {
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
            if ( line.empty() || (line.at(0) == '#') ) {
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
                    matrs.readFromFile( ( baseDir / istr ).string() );
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
                    exporters.readFromFile( ( baseDir / istr ).string(), & nodes, & elems, ndim);
                }
            } else if ( istr.compare("PBlockFiles") == 0 ) {
                iss >> iint;
                for ( int i = 0; i < iint; i++ ) {
                    iss >> istr;
                    pblocks.readFromFile( ( baseDir / istr ).string(), ndim);
                }
            } else if ( initial && istr.compare("Solver") == 0 ) {
                iss >> istr;
                //solver = new Solver;
                //Solver* ptr = solver;
                //solver = solver->readFromFile(( baseDir / istr ).string() );
                //delete ptr;
                solver = Solver().readFromFile( ( baseDir / istr ).string() );
                // QUESTION JK: why is this here and not in the constructor? together with new Solver() ?
                solver->setContainers(& elems, & nodes, & funcs);
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

    //here we apply periodic blocks to generate all the necessary objects
    //it was removed from Model initialization, because it had to be called in advance for RVE materials
    pblocks.setContainers(& nodes, & elems, & bconds, & constr, & funcs, & exporters, & matrs, solver);
    pblocks.apply();


    exporters.setResultDirectory(resultDir);
    exporters.setSolver(solver);
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
    elems.setContainers(& nodes, & bconds);
    pblocks.setContainers(& nodes, & elems, & bconds, & constr, & funcs, & exporters, & matrs, solver);
    // std :: cout << "step: " << solver->giveStepNumber() << ", time: " << solver->giveTime() << '\n';

    solver->setContainers(& elems, & nodes, & funcs);
    // std :: cout << "step: " << solver->giveStepNumber() << ", time: " << solver->giveTime() << '\n';
}
