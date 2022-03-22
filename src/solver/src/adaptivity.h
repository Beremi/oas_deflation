#ifndef _ADAPTIVITY_H
#define _ADAPTIVITY_H

#include <vector>
#include <iostream>
#include <fstream>
#include <iomanip>      // std::setw
#include <iterator>
#include <memory>
// #include <sanitizer/lsan_interface.h>
#include "solver.h"
#include "data_exporter.h"
#include "geometry.h"
#include "model.h"

#define PRINT_TEST false


///////////////////////////////////////////////////////////////////////////////
//// on conditional inheritance:
//// https://stackoverflow.com/questions/16358804/c-is-conditional-inheritance-possible
///////////////////////////////////////////////////////////////////////////////

/*
 * NOTE JK:
 * adaptivity is based on maximum principal fabric_stress calculated in nodes
 */



template< typename BaseSolver >
class AdaptiveSolver : public BaseSolver
{
private:
    std :: string prepInput; ///> modelName used in preprocessor
    std :: string remeshDir;
    fs :: path preprocessorDir;

    fs :: path nodeFile; // path to initial node file - particles must be in a single file

    NodeContainer *nodesFine = nullptr;
    fs :: path pathToFineNodes;

    // split string:
    // http://www.martinbroadhurst.com/how-to-split-a-string-in-c.html
    fs :: path elemStatuses;

    unsigned dim;
    double adaptThreshold; // user input value of stress
    double radius; // user input value of radius to remesh around node over threshold
    double radius2; // user input value of radius of transitional area between remesher and inital lmin
    double remesherLmin = 0;

    double time_before_step; // if remesh, step will be restatrted

    std :: vector< Point >nodeCentersToRmesh;   // clear after remesh

    std :: string regionsToSkip;
    std :: vector< std :: unique_ptr< Region > >regionsNotToRemesh;  // can be specified on the input
    std :: vector< std :: unique_ptr< Region > >fineRegions;   // regions previously remeshed
    unsigned remeshMaterialId;

    std :: vector< unsigned >nodesToKeep;
    unsigned remesherSeed = 1;
    bool reseted = true;
    std :: string additional_procedures = "none";

    //////////////////////////////////////////////////////////////////////////////
    void saveCenters(const std :: string &centersFName, const std :: vector< Point > &centersPoints) {
        std :: ofstream outputfile(fs :: path(this->remeshDir) / centersFName);
        if ( outputfile.is_open() ) {
            outputfile << std :: scientific;
            outputfile.precision(6);
            outputfile << "#nodeCentersToRemesh";
            for ( auto const &c : centersPoints ) {
                outputfile << "\ncenter\t" << c.x() << '\t' << c.y() << '\t' << c.z();
            }
            outputfile.close();
        }
    }

    void saveCenters() {
        // overloading previous method to give it default arguments (to be able to pass them by reference)
        saveCenters("centersToRemesh.out", this->nodeCentersToRmesh);
    }

    void saveNodesToKeep() {
        Node *n;
        Point p;
        // create vector of centers to remove nodes from
        std :: vector< std :: unique_ptr< Region > >regionsToRemove;
        // regionsToRemove.resize( this->nodeCentersToRmesh.size() );
        // unsigned rr = 0;
        // for ( auto const &cent : this->nodeCentersToRmesh ) {
        //     regionsToRemove [ rr++ ] = new Sphere(cent, this->radius2);
        // }
        for ( auto const &cent : this->nodeCentersToRmesh ) {
            regionsToRemove.push_back( std :: make_unique< Sphere >(cent, this->radius2) );
        }
        // save nodes that are going to be kept
        // maybe here can be nodes.out to distinguish between old and the new ones
        std :: string node_file = "nodes.out";
        for ( unsigned i = 0; i < BaseSolver :: nodes->giveSize(); i++ ) { // foreach loop does not work here
            n = BaseSolver :: nodes->giveNode(i);
            if ( n->giveName().compare("particle") == 0 || n->giveName().compare("Particle") == 0 ) {
                p = n->givePoint();
                // check if node is in region to remesh
                if ( isInsideRegions(regionsToRemove, p) && !isInsideRegions(this->regionsNotToRemesh, p) && !isInsideRegions(this->fineRegions, p) ) {
                    continue;
                }
                this->nodesToKeep.push_back(i);
            }
        }

        BaseSolver :: nodes->saveToFile(
            ( fs :: path(this->remeshDir) / node_file ).string(),
            this->nodesToKeep
            );
        // for ( auto delr : regionsToRemove ) {
        //     delete delr;
        // }
    };

    void saveNodesFine() {
        Node *n;
        Point p;
        std :: vector< unsigned >nodeIdsToSave;
        // create vector of centers to remove nodes from
        // std :: vector< Region * >regionsToRemove;
        // regionsToRemove.resize( this->nodeCentersToRmesh.size() );
        // unsigned rr = 0;
        // for ( auto const &cent : this->nodeCentersToRmesh ) {
        //     regionsToRemove [ rr++ ] = new Sphere(cent, this->radius);
        // }
        std :: vector< std :: unique_ptr< Region > >regionsToRemove;
        for ( auto const &cent : this->nodeCentersToRmesh ) {
            regionsToRemove.push_back( std :: make_unique< Sphere >(cent, this->radius2) );
        }
        // save nodes that are going to be kept
        // maybe here can be nodes.out to distinguish between old and the new ones
        std :: string node_file = "nodesFine.out";
        for ( unsigned i = 0; i < this->nodesFine->giveSize(); i++ ) { // foreach loop does not work here
            n = this->nodesFine->giveNode(i);
            if ( n->giveName().compare("particle") == 0 || n->giveName().compare("Particle") == 0 ) {
                p = n->givePoint();
                // check if node is in region to remesh
                if ( isInsideRegions(regionsToRemove, p) && !isInsideRegions(this->regionsNotToRemesh, p) && !isInsideRegions(this->fineRegions, p) ) {
                    nodeIdsToSave.push_back(i);
                }
            }
        }

        this->nodesFine->saveToFile(
            ( fs :: path(this->remeshDir) / node_file ).string(),
            nodeIdsToSave
            );
        // for ( auto delr : regionsToRemove ) {
        //     delete delr;
        // }
    };

    //////////////////////////////////////////////////////////////////////////////
    void saveElemStatuses() {
        std :: vector< unsigned >elems_to_save;
        Element *el;
        Sphere sph;
        for ( unsigned i = 0; i < BaseSolver :: elems->giveSize(); i++ ) {
            el = BaseSolver :: elems->giveElement(i);
            if ( isInsideRegions(this->fineRegions, el) ) {
                for ( auto const &mstat : el->giveMaterialStats() ) {
                    if ( !mstat->isElastic(false) ) { // save elems that were already damaged (in past) now=false checks damage (true checks temp_damage)
                        elems_to_save.push_back(el->giveID() );
                        break;
                    }
                }
            }
        }
        this->elemStatuses = fs :: path(this->remeshDir) / "elemStats.out";
        BaseSolver :: elems->saveElemStatsToFile( ( this->elemStatuses ).string(),
                                                  elems_to_save
                                                   );
    }
    //////////////////////////////////////////////////////////////////////////////
    void saveRemeshData() {
        if ( PRINT_TEST ) { std :: cout << "adaptivity remesh II d - saveRemeshData" << '\n'; }
        this->saveCenters(); // save centersToRemesh
        std :: vector< Point >fine_centers;
        for ( auto const &reg : this->fineRegions ) {
            fine_centers.push_back(reg->giveMainPoint() );
        }
        this->saveCenters("centersFine.out", fine_centers);   // save any specified vector of points

        this->saveNodesToKeep();
        if ( this->nodesFine ) {
            this->saveNodesFine();
        }

        this->saveElemStatuses();
    }

    //////////////////////////////////////////////////////////////////////////////
    void updateGeometry() {
        if ( PRINT_TEST ) { std :: cout << "adaptivity remesh II e - updateGeoemtry" << '\n'; }
        if ( system(NULL) ) {
            std :: cout << "preparing to run preprocessor to remesh geometry" << '\n';
            ;
        } else {
            std :: cerr << "cannot run system command" << '\n';
            exit(EXIT_FAILURE);
        }
        std :: string remeshCmd = "python " + ( this->preprocessorDir / "Remesh.py" ).string() +
                                  " " + this->prepInput + " " + this->remeshDir + " " +
                                  GlobPaths :: BASEDIR.string() + " " +
                                  std :: to_string(this->radius) + " " +
                                  std :: to_string(this->radius2)
                                  + " " +
                                  std :: to_string(int( this->nodesFine != nullptr ) );
        ;
        remeshCmd = remeshCmd + " " + std :: to_string(this->remesherSeed);
        // regionsNotToRemesh
        remeshCmd = remeshCmd + " " + regionsToSkip;
        // cout << remeshCmd << endl;
        if ( this->remesherLmin != 0 ) {
            remeshCmd = remeshCmd + " " + std :: to_string(this->remesherLmin);
        }



        std :: cout << "system cmd " << remeshCmd << '\n';

        if ( system(remeshCmd.c_str() ) != 0 ) {
            std :: cerr << "something went wrong during remesher run" << '\n';
            exit(EXIT_FAILURE);
        }

        if ( this->additional_procedures != "none" ) {
            remeshCmd = "python " + this->additional_procedures +
            " " + this->remeshDir;

            std :: cout << "additional_python_script cmd " << remeshCmd << '\n';

            if ( system(remeshCmd.c_str() ) != 0 ) {
                std :: cerr << "something went wrong during remesher additional procedures" << '\n';
                exit(EXIT_FAILURE);
            }
            std::cout << "\n\n" << '\n';
        }

    }

    //////////////////////////////////////////////////////////////////////////////
    unsigned giveNewNodeId(const unsigned &oldNodeId) {
        auto res = std :: find(std :: begin(this->nodesToKeep), std :: end(this->nodesToKeep), oldNodeId);
        if ( res == this->nodesToKeep.end() ) {
            // if old id is not in vector, return given oldNodeId
            return oldNodeId;
        }
        return std :: distance(std :: begin(this->nodesToKeep), res);
    }


    void setMaterialInFineRegions() {
        // only the elements in fine regions have nonlinear material
        Element *el;
        unsigned change = 0;
        for ( unsigned i = 0; i < BaseSolver :: elems->giveSize(); i++ ) {
            // change = 0;
            el = BaseSolver :: elems->giveElement(i);
            if ( el->giveNode(0)->doesMechanics() && // NOTE JK: adaptivity is based on mechanical stress only
                 isInsideRegions(this->fineRegions, el) ) {
                if ( PRINT_TEST ) { std :: cout << "adaptivity remesh II g - setMaterialInFineRegions " << change++ << ", " << el->giveName() << '\n'; }
                el->changeMaterial(masterModel->giveMaterials()->giveMaterial(this->remeshMaterialId) );
            }
        }
    }

    void readElemStatuses() {
        std :: string line, param;
        unsigned elem_id, stat_id, mat_id, num, node_id;
        std :: vector< unsigned >node_ids;
        Element *el;
        std :: string node_ids_string;

        std :: ifstream inputfile(this->elemStatuses.string().c_str() );
        if ( inputfile.is_open() ) {
            while ( getline(inputfile >> std :: ws, line) ) {
                if ( line.at(0) == '#' || line.empty() ) {
                    continue;
                }
                std :: istringstream iss(line);
                iss >> param;
                if ( param.compare("matStat") == 0 ) {
                    iss >> elem_id >> stat_id >> mat_id >> num;
                    node_ids_string = "";
                    for ( unsigned i = 0; i < num; i++ ) {
                        iss >> node_id;
                        node_ids_string += "\t" + std :: to_string(node_id);
                        // map old nodes to new
                        node_ids.push_back(this->giveNewNodeId(node_id) );
                    }
                    // find element connecting these nodes
                    el = BaseSolver :: elems->giveElementConnectingNodes(node_ids);
                    if ( el ) {
                        // read material status from file
                        // std::cout << "elem found" << '\n';
                        el->giveMatStatus(stat_id)->readFromLine(iss);
                    } else {
                        // if could not find elem connecting nodes, el = nullptr
                        std :: cerr << " old node ids: " << node_ids_string << '\n';
                    }
                    // std::cout << "------------------------------------------------" << '\n';
                }
                node_ids.clear();
            }
            inputfile.close();
        } else {
            std :: cerr << "there is no such file " << this->elemStatuses.string() << '\n';
            exit(EXIT_FAILURE);
        }
    }

    //////////////////////////////////////////////////////////////////////////////
    void loadRemeshData() {
        if ( PRINT_TEST ) { std :: cout << "adaptivity remesh II f - loadRemeshData" << '\n'; }
        masterModel->clear();

        masterModel->readFromFile( ( fs :: path(this->remeshDir) / "master.inp" ).string(), false );

        // update the dof fields etc
        std :: cout << "updated model initialization ..." << '\n';
        masterModel->init(false);

        this->setMaterialInFineRegions();
        this->readElemStatuses();
        BaseSolver :: elems->updateStiffnessMatrix(BaseSolver :: K, "secant");
    }


    //////////////////////////////////////////////////////////////////////////////
    bool checkNodes() {
        if ( PRINT_TEST ) { std :: cout << "adaptivity check nodes II b" << '\n'; }
        std :: vector< Matrix >nodal_stress;
        nodal_stress.resize(BaseSolver :: nodes->giveSize(), Matrix :: Zero(this->dim, this->dim) );
        // calculate nodal stresses
        ExportAllElementsNodalStress(nodal_stress, BaseSolver :: giveDoFValues(), BaseSolver :: giveNodalForces(), BaseSolver :: nodes, BaseSolver :: elems, this->dim);

        std :: vector< Vector >tensorial_stress;
        BaseSolver :: elems->extrapolateValuesFromIntegrationPointsToNodes("solid_stress", tensorial_stress);
        // calculate principal stresses
        Node *n;
        std :: vector< Vector >eigvecs;
        Vector eignums;

        for ( unsigned i = 0; i < BaseSolver :: nodes->giveSize(); i++ ) { // foreach loop does not work here
            n = BaseSolver :: nodes->giveNode(i);
            if ( ( n->giveName().compare("particle") == 0 || n->giveName().compare("Particle") == 0 ) && (tensorial_stress [ i ].size()>0) ) {
                if ( (
                         // JK check point from regions not to remesh, only do not remesh nodes inside of it, that't why the following line commented
                         // !isInsideRegions( this->regionsNotToRemesh, n->givePoint() ) &&
                         !isInsideRegions(this->fineRegions, n->givePoint() )
                         )
                      ) {
                    LinalgEigenSolver(tensorial_stress [ i ], eignums, eigvecs);
                    if ( eignums.maxCoeff() > this->adaptThreshold ) {
                        nodeCentersToRmesh.push_back(n->givePoint() );
                    }
                }
            }
        }
        // if vector empty check returns false, it means not to remesh
        return !this->nodeCentersToRmesh.empty();
    }

    //////////////////////////////////////////////////////////////////////////////
    void remeshGeometry() {
        if ( PRINT_TEST ) { std :: cout << "adaptivity remesh II a" << '\n'; }
        if ( checkNodes() ) {
            if ( PRINT_TEST ) { std :: cout << "adaptivity remesh II c" << '\n'; }
            std :: ostringstream stringStream;
            stringStream << "remesh_" << BaseSolver :: step;
            this->remeshDir = ( GlobPaths :: BASEDIR / stringStream.str() ).string();
            // TODO do the routines to remesh
            std :: cout << "remeshDir: " << this->remeshDir << '\n';
            if ( !fs :: exists(this->remeshDir) ) {
                fs :: create_directories(this->remeshDir);
            }

            this->saveRemeshData(); // save nodes, regions to remesh and element statuses

            this->updateGeometry(); // run python preprocessor

            for ( auto const &p : nodeCentersToRmesh ) {
                this->fineRegions.push_back(std :: make_unique< Sphere >(p, this->radius) );
            }

            this->loadRemeshData(); // load updated geometry


            this->nodeCentersToRmesh.clear();
            this->nodesToKeep.clear();

            if ( PRINT_TEST ) { std :: cout << "-------------------->>>>>>>>>>> solve after remesh" << '\n'; }

            BaseSolver :: dt = BaseSolver :: time - this->time_before_step;
            BaseSolver :: step--;
            BaseSolver :: time = this->time_before_step;
            this->reseted = true;
            BaseSolver :: reset();  // zakázat další adaptivitu v resetovaném kroku            
            BaseSolver :: runBeforeEachStep();




            // BaseSolver :: load *= 0.0;
            // BaseSolver :: ddr *= 0.0;

            BaseSolver :: solve();
            if ( PRINT_TEST ) { std :: cout << "solution after remesh  done" << '\n'; }

            // return true;
        }
        // return false;
    };



    void readAdaptivity(const std :: string filename) {
        std :: string param, path, line;
        bool bat, br, br2, bptp, bmn, brl;
        bat = br = br2 = bptp = bmn = brl = false;
        std :: ifstream inputfile(filename.c_str() );
        if ( inputfile.is_open() ) {
            while ( getline(inputfile >> std :: ws, line) ) {
                if ( line.empty() ) {
                    continue;
                }
                if ( line.at(0) == '#' ) {
                    continue;
                }
                std :: istringstream iss(line);
                iss >> param;
                if ( param.compare("adaptThreshold") == 0 ) {
                    iss >> this->adaptThreshold;
                    bat = true;
                } else if ( param.compare("remesherLmin") == 0 ) {
                    iss >> this->remesherLmin;
                    brl = true;
                } else if ( param.compare("radius") == 0 ) {
                    iss >> this->radius;
                    br = true;
                } else if ( param.compare("radius2") == 0 ) {
                    iss >> this->radius2;
                    br2 = true;
                } else if ( param.compare("preprocessorDir") == 0 ) {
                    iss >> path;
                    this->preprocessorDir = fs :: absolute(path);
                    bptp = true;
                } else if ( param.compare("elemStatuses") == 0 ) {
                    iss >> path;
                    this->elemStatuses = GlobPaths :: BASEDIR / path;
                } else if ( param.compare("prepInput") == 0 ) {
                    iss >> this->prepInput;
                    bmn = true;
                } else if ( param.compare("remesherSeed") == 0 ) {
                    iss >> this->remesherSeed;
                } else if ( param.compare("regionsToSkip") == 0 ) {
                    // std::cout << "reading regions >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>" << '\n';
                    iss >> path;
                    regionsToSkip = ( GlobPaths :: BASEDIR / path ).string();
                    readRegions(regionsToSkip, this->regionsNotToRemesh);
                } else if ( param.compare("remeshMaterialId") == 0 ) {
                    // std::cout << "reading regions >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>" << '\n';
                    iss >> this->remeshMaterialId;
                } else if ( param.compare("additional_python_script") == 0 ) {
                    iss >> this->additional_procedures;
                } else if ( param.compare("pathToFineNodes") == 0 ) {
                    iss >> path;
                    this->pathToFineNodes = GlobPaths :: BASEDIR / path;
                    if ( fs :: exists(this->pathToFineNodes) ) {
                        this->nodesFine = new NodeContainer();
                    } else {
                        std :: cerr << "file with specified fine geoemtry ' " << this->pathToFineNodes.string() << "' does not exist, will use randomly generated geoemtry for refined regions" << '\n';
                    }
                }
            }
            inputfile.close();
        }
        if ( !bat ) {
            std :: cerr << BaseSolver :: name << " parameter 'adaptThreshold' not specified" << '\n';
            exit(EXIT_FAILURE);
        }
        if ( !br ) {
            std :: cerr << BaseSolver :: name << " parameter 'radius' not specified" << '\n';
            exit(EXIT_FAILURE);
        }
        if ( !br2 ) {
            radius2 = 2.0 * radius;
            std :: cerr << BaseSolver :: name << " parameter 'radius2' not specified, taking radius2 = " << radius2 << '\n';
        }
        if ( !bptp ) {
            std :: cerr << BaseSolver :: name << " 'pathToPreprocessor' not specified" << '\n';
            exit(EXIT_FAILURE);
        }
        if ( !bmn ) {
            std :: cerr << BaseSolver :: name << " parameter  'modelName' not specified" << '\n';
            exit(EXIT_FAILURE);
        }
        if ( !brl ) {
            std :: cerr << BaseSolver :: name << " parameter  'remesherLmin' not specified, will take 1/3 of initial Lmin" << '\n';
            // exit(EXIT_FAILURE);
        }
    };

protected:


public:
    AdaptiveSolver() { BaseSolver :: name.append("-Adaptive"); };
    virtual ~AdaptiveSolver() {};

    virtual void init(std :: string init_r_file, std :: string init_v_file, const bool initial = true) {
        unsigned numParticles = 0;
        Node *n;
        // NOTE JK this is intended for field of nodes to keep
        for ( unsigned i = 0; i < BaseSolver :: nodes->giveSize(); i++ ) {
            n = BaseSolver :: nodes->giveNode(i);
            if ( n->giveName().compare("particle") == 0 || n->giveName().compare("Particle") == 0 ) {
                numParticles++;
            }
        }
        // dimension is needed for calculation of nodal stresses (TODO calulate particle volume somewhere at the beginning)
        this->dim = BaseSolver :: elems->giveElement(0)->giveDimension();
        // std::cout << "elem 0: " << BaseSolver :: elems->giveElement(0)->giveName() << '\n';

        if ( initial && this->nodesFine ) {
            std :: cout << "Adaptivity: loading fine geometry ..." << '\n';
            this->nodesFine->readFromFile( ( this->pathToFineNodes ).string(), this->dim );
        }

        BaseSolver :: init(init_r_file, init_v_file, initial);
    };


    virtual Solver *readFromFile(const std :: string filename) {
        BaseSolver :: readFromFile(filename);

        std :: string param, path, line;
        bool bfa = false;
        std :: ifstream inputfile(filename.c_str() );
        if ( inputfile.is_open() ) {
            while ( getline(inputfile >> std :: ws, line) ) {
                if ( line.empty() ) {
                    continue;
                }
                if ( line.at(0) == '#' ) {
                    continue;
                }
                std :: istringstream iss(line);
                iss >> param;
                if ( param.compare("Adaptivity") == 0 || param.compare("adaptivity") == 0 ) {
                    iss >> path;
                    this->readAdaptivity( ( GlobPaths :: BASEDIR / path ).string() );
                    bfa = true;
                }
            }
            inputfile.close();
        }
        if ( !bfa ) {
            std :: cerr << BaseSolver :: name << " no file containing adaptivity data specified" << '\n';
            exit(EXIT_FAILURE);
        }
        return this;
    };


    virtual void runBeforeEachStep() {
        if ( PRINT_TEST ) { std :: cout << "adaptivity before each step I" << '\n'; }
        this->time_before_step = BaseSolver :: time;
        BaseSolver :: runBeforeEachStep();
    };


    virtual void runAfterEachStep() {
        if ( PRINT_TEST ) { std :: cout << "adaptivity after each step II" << '\n'; }
        if ( this->reseted ) {
            this->reseted = false;
        } else {
            remeshGeometry();
        }

        BaseSolver :: runAfterEachStep();
        // __lsan_do_leak_check();
    };

    virtual void solve() {
        if ( PRINT_TEST ) { std :: cout << "adaptivity solving each step II" << '\n'; }
        // TODO JK: check nodes before solving - export fabric stress with frozen variables
        // std :: cout << "solving with " << BaseSolver :: name << '\n';
        BaseSolver :: solve();
    };
};
////////////////////////////////////////////////////////////////
// Foo<Advanced>  fooWithAdvanced;
// Foo<Basic>     fooWithBasic;
// Foo<OtherBase> fooWithOtherSolver;

#endif /* _ADAPTIVITY_H */
