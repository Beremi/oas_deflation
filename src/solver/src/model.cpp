#include "model.h"

//////////////////////////////////////////////////////////
Model :: Model(bool pT){
    printTime = pT;
    nodes.setContainers(&bconds, &constr);
    bconds.setContainers(&funcs);
    elems.setContainers(& nodes);
    pblocks.setContainers(& nodes, & elems, & bconds, & constr, & funcs, & exporters);
}

//////////////////////////////////////////////////////////
void Model::init() {    //initialization
    pblocks.apply();
    bconds.init();
    nodes.init();
    matrs.init();
    elems.init();
    constr.init(& nodes);
    elems.findElementFriends();
    exporters.init();
    solver->init();
}

//////////////////////////////////////////////////////////
void Model:: solve(){
    //solution
    while ( !solver->isTerminated() ) {
        auto start_part = std :: chrono :: system_clock :: now();
        solver->solveStep();
        exporters.exportData(solver->giveStepNumber(), solver->giveTime(), solver->giveDoFValues(), solver->giveNodalForces(), solver->isTerminated() );
        if ( printTime ) {
            auto now = std :: chrono :: system_clock :: now();
            auto elapsed_seconds = now - start_part;
            std :: cout << "step duration: " << convertTimeToString(elapsed_seconds) << endl;
            cout.flush();
        }
    }
}

//////////////////////////////////////////////////////////
void Model::readFromFile(const string filename) {

    fs :: path fullPath = fs :: absolute(filename);
    baseDir = fullPath.parent_path();
    fs :: path resultDir = baseDir / "results";
    
    exporters.setResultDirectory(resultDir);

    string istr, line;
    int iint, dimension;
    ifstream inputfile( fullPath.string() );
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
                    nodes.readFromFile( ( baseDir / istr ).string(), dimension );
                }
            } else if ( istr.compare("MatFiles") == 0 ) {
                iss >> iint;
                for ( int i = 0; i < iint; i++ ) {
                    iss >> istr;
                    matrs.readFromFile( ( baseDir / istr ).string() );
                }
            } else if ( istr.compare("ElemFiles") == 0 ) {
                iss >> iint;
                for ( int i = 0; i < iint; i++ ) {
                    iss >> istr;
                    elems.readFromFile( ( baseDir / istr ).string(), dimension, &matrs );
                }
            } else if ( istr.compare("ConstrFiles") == 0 ) {
                // read constraint files
                iss >> iint;
                for ( int i = 0; i < iint; i++ ) {
                    iss >> istr;
                    constr.readFromFile( ( baseDir / istr ).string(), dimension, &nodes );
                }
            } else if ( istr.compare("BCFiles") == 0 ) {
                iss >> iint;
                for ( int i = 0; i < iint; i++ ) {
                    iss >> istr;
                    bconds.readFromFile( ( baseDir / istr ).string(), &nodes );
                }
            } else if ( istr.compare("FunctionFiles") == 0 ) {
                iss >> std :: skipws >> iint;
                for ( int i = 0; i < iint; i++ ) {
                    iss >> std :: skipws >> istr;
                    funcs.readFromFile( ( baseDir / istr ).string() );
                }
            } else if ( istr.compare("ExporterFiles") == 0 ) {
                iss >> iint;
                for ( int i = 0; i < iint; i++ ) {
                    iss >> istr;
                    exporters.readFromFile( ( baseDir / istr ).string(), &nodes, &elems, dimension );
                }
            } else if ( istr.compare("PBlockFiles") == 0 ) {
                iss >> iint;
                for ( int i = 0; i < iint; i++ ) {
                    iss >> istr;
                    pblocks.readFromFile( ( baseDir / istr ).string(), dimension );
                }
            } else if ( istr.compare("Solver") == 0 ) {
                iss >> istr;
                solver = new Solver();
                solver = solver->readFromFile( ( baseDir / istr ).string() );
                solver->setContainers(&elems, &nodes, &funcs);
            }
        }
        inputfile.close();
        cout << "Master input file '" <<  fullPath.string() << "' succesfully loaded" << endl;
    } else {
        cerr << "Error: unable to open input file '" <<  fullPath.string() <<  "'" << endl;
        exit(EXIT_FAILURE);
    }
    
}
