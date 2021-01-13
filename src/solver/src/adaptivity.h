#ifndef SOLVER_ADAPTIVE_H
#define SOLVER_ADAPTIVE_H

#include <vector>
#include <iostream>
#include <fstream>
#include <iomanip>      // std::setw
#include <iterator>
#include "solver.h"
#include "data_exporter.h"
#include "geometry.h"
#include "model.h"


Vector calcPrincipalStress(const Matrix &stress) {
  Vector principalStress;
  principalStress.resize(stress.numCols());
  if ( principalStress.size() == 2 ){
    // 2d case
    principalStress[0] = 0.5 * (stress[0][0] + stress[1][1]) +
        sqrt(
          pow( 0.5 * (stress[0][0] + stress[1][1]), 2) +
          pow( 0.5 * (stress[0][1] + stress[1][0]), 2)
        );
    principalStress[1] = 0.5 * (stress[0][0] + stress[1][1]) -
        sqrt(
          pow( 0.5 * (stress[0][0] + stress[1][1]), 2) +
          pow( 0.5 * (stress[0][1] + stress[1][0]), 2)
        );
  } else if ( principalStress.size() == 3 ) {
    // http://www.continuummechanics.org/principalstress.html
    double I1, I2, I3;  // invariants
    double Q, R, theta;
    // calculate invariants
    I1 = stress[0][0] + stress[1][1] + stress[2][2];
    I2 = stress[0][0] * stress[1][1] +
         stress[1][1] * stress[2][2] +
         stress[0][0] * stress[2][2] +
         stress[0][1] * stress[1][0] +   // for symetric matrix = stress12^2
         stress[1][2] * stress[2][1] +
         stress[0][2] * stress[2][0];
    I3 = stress[0][0] * stress[1][1] * stress[2][2] +
         stress[0][0] * stress[1][2] * stress[2][1] +
         stress[1][1] * stress[0][2] * stress[2][0] +
         stress[2][2] * stress[0][2] * stress[2][0] +
         stress[0][1] * stress[1][2] * stress[0][2] +
         stress[1][0] * stress[2][1] * stress[2][0];   // last two could be done once and multiplied by 2
     // calculate intermediate quantities
     Q = (3 * I2 - pow(I1, 2) ) / 9;
     R = ( 2 * pow(I1, 3) - 9 * I1 * I2 + 27 * I3) / 54;
     if ( Q >= 0 ){
       theta = 0;
     } else {
       if ( R / sqrt( - pow(Q, 3) ) > 1.0 ){
          theta = acos( 1.0 );
       } else {
         theta = acos(R / sqrt( - pow(Q, 3) ) );
       }
     }
     // calculate eigen values
     principalStress[0] = 2 * sqrt(-Q) * cos(theta / 3) + I1 / 3;
     principalStress[1] = 2 * sqrt(-Q) * cos((theta + 2 * M_PI) / 3) + I1 / 3;
     principalStress[2] = 2 * sqrt(-Q) * cos((theta + 4 * M_PI) / 3) + I1 / 3;
     // stress.print();
     // std::cout << "principalStress ( " << principalStress[0] << ", " << principalStress[1] << ", " << principalStress[2] << " ), Q = " << Q << ", R = " << R << ", theta = " << theta << ", acos argument = " << R / sqrt( - pow(Q, 3) ) << '\n';
  }
  return principalStress;
}

double calcMaxPrincipalStress(const Matrix &stress) {
    return calcPrincipalStress(stress).max();
}


///////////////////////////////////////////////////////////////////////////////
//// on conditional inheritance:
//// https://stackoverflow.com/questions/16358804/c-is-conditional-inheritance-possible
///////////////////////////////////////////////////////////////////////////////

/*
 * NOTE JK:
 * adaptivity is based on maximum principal fabric_stress calculated in nodes - nodal_stress in VTK exporters
 */



template< typename BaseSolver >
class AdaptiveSolver : public BaseSolver
{
private:
    std :: string prepInput; ///> modelName used in preprocessor
    std :: string remeshDir;
    fs :: path preprocessorDir;

    fs :: path nodeFile; // path to initial node file - particles must be in a single file

    // split string:
    // http://www.martinbroadhurst.com/how-to-split-a-string-in-c.html
    fs :: path elemStatuses;
    fs :: path nodeStatuses;


    unsigned dim;
    double adaptThreshold; // user input value of stress
    double radius; // user input value of radius to remesh around node over threshold
    double radius2; // user input value of radius of transitional area between remesher and inital lmin
    double remesherLmin = 0;

    double time_before_step; // if remesh, step will be restatrted

    // true = keep forever, false = remesh if needed
    std :: vector< bool >keepThisNodeForever;  // set on init() and set it only for particles

    std :: vector< Point >nodeCentersToRmesh;  // clear after remesh
    // std :: vector< Point > centersPreviouslyRmeshed;   // clear after remesh

    std :: vector< Region * >regionsNotToRemesh; // can be specified on the input, as well as by centers of already remeshed circles
    // furthermore will be stored in memory

    //////////////////////////////////////////////////////////////////////////////
    void saveCentersToRemesh() {
        std :: string centers = "centersToRemesh.out";
        ofstream outputfile(fs :: path(this->remeshDir) / centers);
        if ( outputfile.is_open() ) {
            outputfile << std :: scientific;
            outputfile.precision(6);
            outputfile << "#nodeCentersToRemesh";
            for ( auto const &c : nodeCentersToRmesh ) {
                outputfile << "\ncenter\t" << c.getX() << '\t' << c.getY() << '\t' << c.getZ();
            }
            outputfile.close();
        }
    }

    void saveNodesToKeep() {
        Node *n;
        Point p;
        bool to_be_exported;
        // create vector of centers to remove nodes from
        std :: vector< Sphere >regionsToRemove;
        for ( auto const &cent : nodeCentersToRmesh ) {
            regionsToRemove.push_back(Sphere(cent, this->radius2) );
        }

        // save nodes that are going to be kept
        // maybe here can be nodes.out to distinguish between old and the new ones
        std :: string node_file = "nodes.inp";
        ofstream outputfile(fs :: path(this->remeshDir) / node_file);
        if ( outputfile.is_open() ) {
            outputfile << std :: scientific;
            outputfile.precision(6);
            outputfile << "#nodesToLoad";
            for ( unsigned i = 0; i < BaseSolver :: nodes->giveSize(); i++ ) { // foreach loop does not work here
                n = BaseSolver :: nodes->giveNode(i);
                if ( n->giveName().compare("particle") == 0 ) {
                    to_be_exported = true;
                    Particle *part = static_cast< Particle * >( n );
                    p = n->givePoint();
                    // check if node is in region to remesh
                    for ( auto const &regRemove : regionsToRemove ) {
                        if ( regRemove.isInside(p) ) {
                            to_be_exported = false;
                            break;
                        }
                    }
                    // also check if it is not in region where remesh is restricted
                    // TODO this should be done in better way
                    for ( auto const &regKeep : regionsNotToRemesh ) {
                        if ( regKeep->isInside(p) ) {
                            to_be_exported = true;
                            break;
                        }
                    }
                    if ( to_be_exported ) {
                        outputfile << "\nparticle\t" << p.getX() << '\t' << p.getY() << '\t' << p.getZ() << '\t' << part->giveRadius();
                    }
                }
            }
            outputfile.close();
        }
    };

    //////////////////////////////////////////////////////////////////////////////
    // void saveNodeStatuses() {
    //   //TODO  save this only for Particles, otherwise it will be loading and saving a lot of useless data
    //   // NOTE - need to distinguish between structural particle and master nodes that are purely virtual
    //   std :: string centers = "nodalStatuses.out";
    //   ofstream outputfile( fs :: path(this->remeshDir) / centers );
    //   if ( outputfile.is_open() ) {
    //     outputfile << std :: scientific;
    //     outputfile.precision(6);
    //     outputfile << "#nodalStatuses - keepThisNodeForever or not?";
    //     for ( auto const & ke : keepThisNodeForever ) {
    //       outputfile << "\n" << int(ke);
    //     }
    //     outputfile.close();
    //   }
    // }

    //////////////////////////////////////////////////////////////////////////////
    void saveElemStatuses() {
        // TODO
    }
    //////////////////////////////////////////////////////////////////////////////
    void saveRemeshData() {
        saveCentersToRemesh();
        // saveNodeStatuses();

        saveNodesToKeep();

        // NOTE maybe there is no need to save elems, maybe I can just save field relating elems and nodes and after remsh just get the information about which elem it was previously and copy it into some new container
        // TODO prepare saveElemStatuses for the future anyway
        saveElemStatuses();
    }

    //////////////////////////////////////////////////////////////////////////////
    void updateGeometry() {
        if ( system(NULL) ) {
            std :: cout << "preparing to run preprocessor to remesh geometry" << '\n';
            ;
        } else   {
            std :: cerr << "cannot run system command" << '\n';
            exit(EXIT_FAILURE);
        }
        std :: string remeshCmd = "python " + ( this->preprocessorDir / "Remesh.py" ).string() +
                                  " " + this->prepInput + " " + this->remeshDir + " " +
                                  GlobPaths :: BASEDIR.string() + " " +
                                  std :: to_string(this->radius) + " " +
                                  std :: to_string(this->radius2);
        // cout << remeshCmd << endl;
        if ( this->remesherLmin != 0 ) {
            remeshCmd = remeshCmd + " " + std :: to_string(remesherLmin);
        }

        // std::cout << "system return " << system ( remeshCmd.c_str() ) << '\n';

        if ( system(remeshCmd.c_str() ) != 0 ) {
            std :: cerr << "something went wrong during remesher run" << '\n';
            exit(EXIT_FAILURE);
        }
    }

    //////////////////////////////////////////////////////////////////////////////
    // void readNodeStatuses() {
    //   string line;
    //   bool remeshAllowed;
    //   ifstream inputfile(this->nodeStatuses.c_str() );
    //   if ( inputfile.is_open() ) {
    //     while ( getline(inputfile >> std :: ws, line) ) {
    //       if ( line.empty() ) {
    //           continue;
    //       }
    //       if ( line.at(0) == '#' ) {
    //           continue;
    //       }
    //       istringstream iss(line);
    //       iss >> remeshAllowed;
    //       this->keepThisNodeForever.push_back(remeshAllowed);
    //     }
    //     inputfile.close();
    //     cout << "Input file '" <<  this->nodeStatuses.string() << "' succesfully loaded; " << keepThisNodeForever.size() << " nodes found" << '\n';
    //   } else {
    //       cerr << "Error: unable to open input file '" <<  this->nodeStatuses.string() <<  "'" << endl;
    //       exit(EXIT_FAILURE);
    //   }
    // }

    //////////////////////////////////////////////////////////////////////////////
    void readElemStatuses() {
        // TODO
    }

    //////////////////////////////////////////////////////////////////////////////
    void loadRemeshData() {
        std :: cout << "old nodes size: " << BaseSolver :: nodes->giveSize() << std :: endl;
        std :: cout << "old elem size: " << BaseSolver :: elems->giveSize() << std :: endl;

        masterModel->clear();

        std :: cout << "step: " << BaseSolver :: giveStepNumber() << ", time: " << BaseSolver :: giveTime() << '\n';

        masterModel->readFromFile( ( fs :: path(this->remeshDir) / "master.inp" ).string(), false );
        // masterModel->print_res_dir();


        // update the dof fields etc
        std :: cout << "updated model initialization ..." << '\n';
        masterModel->init(false);

        std :: cout << "step: " << BaseSolver :: giveStepNumber() << ", time: " << BaseSolver :: giveTime() << '\n';


        std :: cout << "updated nodes size: " << BaseSolver :: nodes->giveSize() << std :: endl;
        std :: cout << "updated elem size: " << BaseSolver :: elems->giveSize() << std :: endl;

        // TODO JK how to keep previous elem container just to keep statuses corresponding to not remeshed elems?
    }


    //////////////////////////////////////////////////////////////////////////////
    bool checkNodes() {
        std :: vector< Matrix >nodal_stress;
        nodal_stress.resize(BaseSolver :: nodes->giveSize(), Matrix(this->dim, this->dim) );
        // calculate nodal stresses
        ExportAllElementsNodalStress(nodal_stress, BaseSolver :: giveDoFValues(), BaseSolver :: giveNodalForces(), BaseSolver :: nodes, BaseSolver :: elems, this->dim);
        // calcuulate principal stresses
        Node *n;
        bool checkStress;
        for ( unsigned i = 0; i < BaseSolver :: nodes->giveSize(); i++ ) { // foreach loop does not work here
            checkStress = true;
            n = BaseSolver :: nodes->giveNode(i);
            if ( n->giveName().compare("particle") == 0 ) {
                for ( auto const &reg : regionsNotToRemesh ) {
                    if ( reg->isInside(n->givePoint() ) ) {
                        checkStress = false;
                        break;
                    }
                }
                if ( checkStress ) {
                    if ( calcMaxPrincipalStress(nodal_stress [ i ]) > adaptThreshold ) {
                        nodeCentersToRmesh.push_back(n->givePoint() );
                    }
                }
            }
        }

        // if vector empty check returns false means not remesh
        return !nodeCentersToRmesh.empty();
    }

    //////////////////////////////////////////////////////////////////////////////
    bool remeshGeometry() {
        if ( checkNodes() ) {
            std :: ostringstream stringStream;
            stringStream << "remesh_" << BaseSolver :: step;
            this->remeshDir = stringStream.str();
            // TODO do the routines to remesh
            std :: cout << "remeshDir: " << this->remeshDir << '\n';
            if ( !fs :: exists(GlobPaths :: BASEDIR / this->remeshDir) ) {
                fs :: create_directories(GlobPaths :: BASEDIR / this->remeshDir);
            }

            saveRemeshData(); // save nodes, regions to remesh and element statuses

            updateGeometry(); // run python preprocessor

            loadRemeshData(); // load updated geometry

            for ( auto const &p : nodeCentersToRmesh ) {
                regionsNotToRemesh.push_back(new Sphere(p, this->radius) );
            }

            nodeCentersToRmesh.clear();

            BaseSolver :: solve();

            // std::cout << "regionsNotToRemesh.size(): " << regionsNotToRemesh.size() << '\n';
            // std::cout << "<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<" << '\n';
            // exit(0);
            return true;
        }
        return false;
    };



    void readAdaptivity(const string filename) {
        string param, path, line;
        bool bat, br, br2, bptp, bmn, brl;
        bat = br = br2 = bptp = bmn = brl = false;
        ifstream inputfile(filename.c_str() );
        if ( inputfile.is_open() ) {
            while ( getline(inputfile >> std :: ws, line) ) {
                if ( line.empty() ) {
                    continue;
                }
                if ( line.at(0) == '#' ) {
                    continue;
                }
                istringstream iss(line);
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
                } else if ( param.compare("nodeStatuses") == 0 ) {
                    iss >> path;
                    this->nodeStatuses = GlobPaths :: BASEDIR / path;
                } else if ( param.compare("prepInput") == 0 ) {
                    iss >> this->prepInput;
                    bmn = true;
                } else if ( param.compare("regionsToSkip") == 0 ) {
                    // std::cout << "reading regions >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>" << '\n';
                    iss >> path;
                    readRegions(path, this->regionsNotToRemesh);
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

    virtual void init(const bool &initial = true) {
        unsigned numParticles = 0;
        Node *n;
        // NOTE JK this is intended for field of nodes to keep
        for ( unsigned i = 0; i < BaseSolver :: nodes->giveSize(); i++ ) {
            n = BaseSolver :: nodes->giveNode(i);
            if ( n->giveName().compare("particle") == 0 ) {
                numParticles++;
            }
        }
        // initially nodestatuses does not have to be specified (usually at the beginning of calculation)
        // NOTE this was replaced with regions not to remesh!
        // if ( fs :: exists( this->nodeStatuses ) ){
        //   // readNodeStatuses();
        //   if ( numParticles != keepThisNodeForever.size() ) {
        //     std::cerr << "number of nodalStatuses = " << keepThisNodeForever.size() << " does not correspond to number of particles in the model (" << numParticles << ")" << '\n';
        //     exit(EXIT_FAILURE);
        //   }
        // } else {
        //   keepThisNodeForever.resize(numParticles);
        // }

        // dimension is needed for calculation of nodal stresses (TODO calulate particle volume somewhere at the beginning)
        this->dim = BaseSolver :: elems->giveElement(0)->giveDimension();

        BaseSolver :: init(initial);
    };


    virtual Solver *readFromFile(const string filename) {
        BaseSolver :: readFromFile(filename);

        string param, path, line;
        bool bfa = false;
        ifstream inputfile(filename.c_str() );
        if ( inputfile.is_open() ) {
            while ( getline(inputfile >> std :: ws, line) ) {
                if ( line.empty() ) {
                    continue;
                }
                if ( line.at(0) == '#' ) {
                    continue;
                }
                istringstream iss(line);
                iss >> param;
                if ( param.compare("Adaptivity") == 0 || param.compare("adaptivity") == 0 ) {
                    iss >> path;
                    this->readAdaptivity(path);
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
        this->time_before_step = BaseSolver :: time;
        BaseSolver :: runBeforeEachStep();
    };



    virtual void runAfterEachStep() {
        // TODO change if to while
        // if ( remeshGeometry() ) {
        // BaseSolver :: time = this->time_before_step;  //
        // BaseSolver :: step--;
        // }

        while ( remeshGeometry() ) {
            // instead of this, put circles into regions not to remesh
            // centersPreviouslyRmeshed.insert(centersPreviouslyRmeshed.end(),
            //                                 nodeCentersToRmesh.begin(),
            //                                 nodeCentersToRmesh.end() );
        }
        BaseSolver :: runAfterEachStep();
    };

    virtual void solve() {
        std :: cout << "solving with " << BaseSolver :: name << '\n';
        BaseSolver :: solve();
    };
};
////////////////////////////////////////////////////////////////
// Foo<Advanced>  fooWithAdvanced;
// Foo<Basic>     fooWithBasic;
// Foo<OtherBase> fooWithOtherSolver;

#endif
