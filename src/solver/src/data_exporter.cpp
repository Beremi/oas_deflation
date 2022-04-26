#include <algorithm>
#include <iterator>
#include "data_exporter.h"
#include "vtk_exporter.h"
#include "exporter_model.h"
#include "geometry.h"
#include "element_discrete.h"
#include "solver.h"

using namespace std;

//////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////
// BASIC DATA EXPORTER - master class
void DataExporter :: giveFileName(unsigned step, char *buffer) const {
    sprintf(buffer, "%s_%05d.out", filename.c_str(), step);
}

//////////////////////////////////////////////////////////
/*!
 *  DataExporter optional parameters.
 *  Keywords:
 *  - saveEveryTime [float] - save each nth time (TODO: better description)
 *  - saveEveryStep [int] - save each nth step (TODO: better description)
 *  - saveTimes [count] [floats] - save in specified times
 *  - saveSteps [count] [ints] - save in specified steps
 *  - timeShift [float] - time shift (TODO: better description)
 *  - precision [float] - precision of stored values
 */
void DataExporter :: readFromLine(istringstream &iss) {
    iss.clear(); // clear string stream
    iss.seekg(0, iss.beg); //reset position in string stream
    string param;
    // initiate variables in case they are not specified
    saveTime_each = numeric_limits<double>::max();
    saveTime_last = 0;
    saveStep_each = numeric_limits<unsigned>::max();
    saveStep_last = 0;
    saveSteps_idx = 0;
    saveTimes_idx = 0;
    next_time_to_save = 0;
    next_step_to_save = 0;
    int num = 0;
    bool saveTimeStepWasConfigured = false;
    while ( !iss.eof() ) {
        iss >> param;
        if ( param.compare("saveEveryTime") == 0 || param.compare("saveEvery") == 0 || param.compare("timeEach") == 0 ) {
            iss >> saveTime_each;
            saveTimeStepWasConfigured = true;
        } else if ( param.compare("saveEveryStep") == 0 ) {
            iss >> saveStep_each;
            saveTimeStepWasConfigured = true;
        } else if ( param.compare("saveTimes") == 0 ) {
            iss >> num;
            double val = 0;
            for (int i = 0; i < num; i++ ) {
                iss >> val;
                times_to_save.push_back(val);
            }
            sort(times_to_save.begin(), times_to_save.end()); // sort times
            times_to_save.erase(std::unique(times_to_save.begin(), times_to_save.end()), times_to_save.end()); // store only unique values
            saveTimeStepWasConfigured = true;
        } else if ( param.compare("saveSteps") == 0 ) {
            iss >> num;
            unsigned val = 0;
            for (int i = 0; i < num; i++ ) {
                iss >> val;
                steps_to_save.push_back(val);
            }
            sort(steps_to_save.begin(), steps_to_save.end()); // sort steps
            steps_to_save.erase(std::unique(steps_to_save.begin(), steps_to_save.end()), steps_to_save.end()); // store only unique values
            saveTimeStepWasConfigured = true;
        } else if ( param.compare("timeShift") == 0 ) {
            iss >> saveTime_last;
        } else if ( param.compare("precision") == 0 ) {
            iss >> precision;
        } else if ( param.compare("multiplier") == 0 ) {
            iss >> multiplier;
        }
    }
    if (!saveTimeStepWasConfigured) saveStep_each = 1; // save in each step because no export frequency was set
    updateNextTimeToSave(0);
    updateNextStepToSave(0);
}

//////////////////////////////////////////////////////////
bool DataExporter :: doExportNow(const double &time, const unsigned &step) {
    cout << "Time a step: " << time << " : " << next_time_to_save << " : " << step<< " : " << next_step_to_save  << endl;
    if ( (time > next_time_to_save) || (step == next_step_to_save) ) {
        updateNextTimeToSave(time);
        updateNextStepToSave(step);
        return true;
    } else {
        return false;
    }
}

void DataExporter::updateNextTimeToSave(const double &time)
{
    if (time > next_time_to_save) {
        double t = saveTime_last + saveTime_each;
        if (saveTimes_idx < times_to_save.size()){
            if (t > times_to_save[saveTimes_idx]){
                t = times_to_save[saveTimes_idx];
                saveTimes_idx++;
            } else {
                saveTime_last = t;
            }
        } else {
            saveTime_last = t;
        }
        next_time_to_save = t - 1e-12;
        time_last = time;
    }
}


void DataExporter::updateNextStepToSave(const unsigned &step)
{
    if (step == next_step_to_save) {
        unsigned s = saveStep_last + saveStep_each;
        if (saveSteps_idx < steps_to_save.size()){
            if (s > steps_to_save[saveSteps_idx]){
                s = steps_to_save[saveSteps_idx];
                saveSteps_idx++;
            } else {
                saveStep_last = s;
            }
        } else {
            saveStep_last = s;
        }
        next_step_to_save = s;
        step_last = step;
    }
}

//////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////
// EXPORT FROM NODES TO TXT
/*!
 *  Export from nodes to txt.
 *  Parameters:
 *  - filename [string] - file to store results (e.g. LD -> results/LD.out)
 *  - num [int] - number of components to read
 *  - codes [string] - component name (ux, uy, uz, pressure, rotx, roty, rotz) (TODO: links to component names)
 *
 *  These parameters can be followed by optional keywords:
 *  - see DataExporter::readFromLine
 */
void TXTNodalExporter :: readFromLine(istringstream &iss) {
    iss >> filename;
    unsigned num;
    iss >> num;
    codes.resize(num);
    for ( unsigned i = 0; i < num; i++ ) {
        iss >> codes [ i ];
    }
    DataExporter :: readFromLine(iss);
}

//////////////////////////////////////////////////////////
void TXTNodalExporter :: init() {
    unsigned ncod = codes.size();
    maxsize.resize(ncod);
    Vector res;
    Vector fakeDoFs = Vector :: Zero( nodes->giveTotalNumDoFs() );
    for ( unsigned i = 0; i < ncod; i++ ) {
        maxsize [ i ] = 0;
        for ( unsigned n = 0; n < nodes->giveSize(); n++ ) {
            nodes->giveNode(n)->giveDoFBasedValues(codes [ i ], fakeDoFs, res);
            maxsize [ i ] = max< size_t >( maxsize [ i ], res.size() ); // (maxsize[i] < res.size()) ? res.size() : maxsize[i];
        }
    }
}

//////////////////////////////////////////////////////////
void TXTNodalExporter :: exportData(unsigned step, const Vector &DoFs, const Vector &reactions, fs :: path resultDir) const {
    ( void ) reactions;
    char buffer [ 100 ];
    Node *nn;
    giveFileName(step, buffer);
    ofstream outputfile( ( resultDir / buffer ).string() );

    unsigned p;
    if ( outputfile.is_open() ) {
        outputfile << "#nodeID";
        for ( unsigned c = 0; c < codes.size(); c++ ) {
            if ( maxsize [ c ] == 1 ) {
                outputfile << "\t" << codes [ c ];
            } else {
                for ( p = 0; p < maxsize [ c ]; p++ ) {
                    outputfile << "\t" << codes [ c ] << "_" << p;
                }
            }
        }
        outputfile << "\n";

        outputfile << std :: scientific;
        outputfile.precision(precision);
        Vector res;
        for ( unsigned n = 0; n < nodes->giveSize(); n++ ) {
            nn = nodes->giveNode(n);
            outputfile << nn->giveID();
            for ( unsigned c = 0; c < codes.size(); c++ ) {
                nn->giveDoFBasedValues(codes [ c ], DoFs, res);
                for ( p = 0; p < min< size_t >( maxsize [ c ], res.size() ); p++ ) {
                    outputfile << "\t" << res [ p ]*multiplier;
                }
                for ( ; p < maxsize [ c ]; p++ ) {
                    outputfile <<  "\t" << 0;
                }
            }
            outputfile << endl;
        }
        outputfile.close();
    }
}

//////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////
// EXPORT FROM ELEMENTS TO TXT
/*!
 *  Export from elements to txt.
 *  Parameters:
 *  - filename [string] - file to store results (e.g. elem -> results/elem.out)
 *  - num [int] - number of components to read
 *  - codes [string] - component labels
 *
 *  These parameters can be followed by optional keywords:
 *  - see DataExporter::readFromLine
 */
void TXTElementExporter :: readFromLine(istringstream &iss) {
    iss >> filename;
    unsigned num;
    iss >> num;
    codes.resize(num);
    for ( unsigned i = 0; i < num; i++ ) {
        iss >> codes [ i ];
    }
    DataExporter :: readFromLine(iss);
}

//////////////////////////////////////////////////////////
void TXTElementExporter :: init() {
    unsigned ncod = codes.size();
    maxsize.resize(ncod);
    Vector res;
    for ( unsigned i = 0; i < ncod; i++ ) {
        maxsize [ i ] = 0;
        for ( unsigned e = 0; e < elems->giveSize(); e++ ) {
            elems->giveElement(e)->giveValues(codes [ i ], res);
            maxsize [ i ] = max< size_t >( maxsize [ i ], res.size() );
        }
    }
}

//////////////////////////////////////////////////////////
void TXTElementExporter :: exportData(unsigned step, const Vector &DoFs, const Vector &reactions, fs :: path resultDir) const {
    ( void ) DoFs;
    ( void ) reactions;
    char buffer [ 100 ];
    Element *ee;
    giveFileName(step, buffer);
    ofstream outputfile( ( resultDir / buffer ).string() );
    Vector res;
    unsigned p;
    if ( outputfile.is_open() ) {
        outputfile << "#elementID";
        for ( unsigned c = 0; c < codes.size(); c++ ) {
            if ( maxsize [ c ] == 1 ) {
                outputfile << "\t" << codes [ c ];
            } else {
                for ( p = 0; p < maxsize [ c ]; p++ ) {
                    outputfile << "\t" << codes [ c ] << "_" << p;
                }
            }
        }
        outputfile << "\n";

        outputfile << std :: scientific;
        outputfile.precision(precision);
        for ( unsigned e = 0; e < elems->giveSize(); e++ ) {
            ee = elems->giveElement(e);
            outputfile << ee->giveID();
            for ( unsigned c = 0; c < codes.size(); c++ ) {
                ee->giveValues(codes [ c ], res);
                for ( p = 0; p < min< size_t >( maxsize [ c ], res.size() ); p++ ) {
                    outputfile << "\t" << res [ p ]*multiplier;
                }
                for ( ; p < maxsize [ c ]; p++ ) {
                    outputfile << "\t" << 0;
                }
            }
            outputfile << endl;
        }
        outputfile.close();
    }
}

//////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////
// EXPORT FROM INTEGRATION POINTS TO TXT
/*!
 *  Export from Gauss points to txt.
 *  Parameters:
 *  - filename [string] - file to store results (e.g. gauss -> results/gauss.out)
 *  - num [int] - number of components to read
 *  - codes [string] - component labels (TODO: links to possible values)
 *
 *  These parameters can be followed by optional keywords:
 *  - see DataExporter::readFromLine
 */
void TXTIntegrationPointExporter :: readFromLine(istringstream &iss) {
    iss >> filename;
    unsigned num;
    iss >> num;
    codes.resize(num);
    for ( unsigned i = 0; i < num; i++ ) {
        iss >> codes [ i ];
    }
    DataExporter :: readFromLine(iss);
}

//////////////////////////////////////////////////////////
void TXTIntegrationPointExporter :: init() {
    unsigned ncod = codes.size();
    maxsize.resize(ncod);
    Vector res;
    Element *ee;
    size_t nIP;
    for ( unsigned i = 0; i < ncod; i++ ) {
        maxsize [ i ] = 0;
        for ( unsigned e = 0; e < elems->giveSize(); e++ ) {
            ee = elems->giveElement(e);
            nIP = ee->giveNumIP();
            for ( unsigned k = 0; k < nIP; k++ ) {
                ee->giveIPValues(codes [ i ], k, res);
                maxsize [ i ] = max< size_t >( maxsize [ i ], res.size() );
            }
        }
    }
}

//////////////////////////////////////////////////////////
void TXTIntegrationPointExporter :: exportData(unsigned step, const Vector &DoFs, const Vector &reactions, fs :: path resultDir) const {
    ( void ) DoFs;
    ( void ) reactions;
    char buffer [ 100 ];
    Element *ee;
    size_t nIP;
    Vector res;
    giveFileName(step, buffer);
    ofstream outputfile( ( resultDir / buffer ).string() );

    unsigned p;
    if ( outputfile.is_open() ) {
        outputfile << "#elementID\tintpointID";
        for ( unsigned c = 0; c < codes.size(); c++ ) {
            if ( maxsize [ c ] == 1 ) {
                outputfile << "\t" << codes [ c ];
            } else {
                for ( p = 0; p < maxsize [ c ]; p++ ) {
                    outputfile << "\t" << codes [ c ] << "_" << p;
                }
            }
        }
        outputfile << "\n";

        outputfile << std :: scientific;
        outputfile.precision(precision);
        for ( unsigned e = 0; e < elems->giveSize(); e++ ) {
            ee = elems->giveElement(e);
            nIP = ee->giveNumIP();
            for ( unsigned k = 0; k < nIP; k++ ) {
                outputfile << ee->giveID() << "\t" << k;
                for ( unsigned c = 0; c < codes.size(); c++ ) {
                    ee->giveIPValues(codes [ c ], k, res);
                    for ( p = 0; p < min< size_t >( maxsize [ c ], res.size() ); p++ ) {
                        outputfile << "\t" << res [ p ]*multiplier;
                    }
                    for ( ; p < maxsize [ c ]; p++ ) {
                        outputfile << "\t" << 0;
                    }
                }
                outputfile << endl;
            }
        }
        outputfile.close();
    }
}

//////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////
// GAUGE EXPORTERS
void Gauge :: giveFileName(unsigned step, char *buffer) const {
    ( void ) step;
    sprintf( buffer, "%s.out", filename.c_str() );
}


//////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////
// EXPORT OF FORCES
/*!
 *  Export Forces to txt.
 *  Parameters:
 *  - filename [string] - file to store results (e.g. LD -> results/LD.out)
 *  - name - label placed in table header
 *  - codes [string] - force component name (2D - fx, fy, mz) and (3D fx, fy, fz, mx, my, mz)
 *  - num [int] - number of nodes to read
 *  - n [int] - node numbers
 *
 *  These parameters can be followed by optional keywords:
 *  - see DataExporter::readFromLine
 */
void ForceGauge :: readFromLine(istringstream &iss) {
    iss >> this->filename;
    iss >> this->name;
    this->codes.resize(1);
    iss >> this->codes [ 0 ];
    unsigned num;
    std :: string param;
    iss >> param;
    if ( param.compare("block") == 0 || param.compare("coords") == 0 ) {
        std :: string param2;
        bool mech = true;
        iss >> param2;
        if ( param2.compare("mech") == 0 ) {
            mech = true;
        } else if ( param2.compare("trsp") == 0 ) {
            mech = false;
        } else {
            std :: cout << "type of force 'mech' or 'trsp' for ForceGauge not determined, by default, 'mech' is considered" << '\n';
        }
        Block bl;
        bl.readFromLine(iss);
        for ( auto const &nod : * nodes ) {
            if ( bl.isInside( nod->givePoint() ) ) {
                if ( ( nod->doesMechanics() && mech ) || ( nod->doesTransport() && !mech ) ) {
                    this->n.push_back( nodes->giveNodeId(nod) );
                }
            }
        }
    } else {
        num = std :: stoul( param.c_str() );
        this->n.resize(num);
        for ( unsigned i = 0; i < num; i++ ) {
            iss >> this->n [ i ];
        }
    }
    DataExporter :: readFromLine(iss);
}

//////////////////////////////////////////////////////////
ForceGauge :: ForceGauge(string &f, string &gname, string &c, vector< unsigned > &nn, NodeContainer *nc, double m, unsigned dimension) : Gauge(dimension) {
    nodes = nc;
    filename = f;
    name = gname;
    n = nn;
    codes.resize(1);
    codes [ 0 ] = c;
    multiplier = m;
}

//////////////////////////////////////////////////////////
void ForceGauge :: init() {
    saveTime_each = 0;
    saveTime_last = 0;
    DoFs.resize( n.size() );
    for ( unsigned i = 0; i < n.size(); i++ ) {
        DoFs [ i ] = nodes->giveNode(n [ i ])->giveStartingDoF() + nodes->giveNode(n [ i ])->giveOrderOfForceCode(codes [ 0 ]);
    }
    maxsize.resize(1);
    maxsize [ 0 ] = 1;
}


//////////////////////////////////////////////////////////
void ForceGauge :: exportData(unsigned step, const Vector &full_f, const Vector &reactions, fs :: path resultDir) const {
    ( void ) full_f;
    char buffer [ 100 ];
    double value = 0;
    giveFileName(step, buffer);
    ofstream outputfile;
    outputfile.open( ( resultDir / buffer ).string(), ios :: app);
    if ( outputfile.good() ) {
        outputfile << std :: scientific;
        outputfile.precision(precision);
        for ( unsigned i = 0; i < DoFs.size(); i++ ) {
            value += reactions [ DoFs [ i ] ];
        }
        outputfile <<  "\t" << value * multiplier;
    }
    outputfile.close();
}

//////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////
// EXPORT OF DEGREES OF FREEDOM
/*!
 *  Export DoF from nodes to txt.
 *  Parameters:
 *  - filename [string] - file to store results (e.g. LD -> results/LD.out)
 *  - name - label placed in table header
 *  - codes [string] - component name (2D - ux, uy, rz) and (3D ux, uy, uz, rx, ry, rz)
 *  - num [int] - number of nodes to read
 *  - n [int] - node numbers
 *
 *  These parameters can be followed by optional keywords:
 *  - see DataExporter::readFromLine
 */

//////////////////////////////////////////////////////////
DoFGauge :: DoFGauge(string &f, string &gname, string &c, vector< unsigned > &nn, NodeContainer *nc, double m, unsigned dimension) : ForceGauge(f, gname, c, nn, nc, m, dimension) {}

//////////////////////////////////////////////////////////
void DoFGauge :: init() {
    saveTime_each = 0;
    saveTime_last = 0;
    unsigned DoFpos = 0;
    if ( codes [ 0 ].compare("ux") == 0 ) {
        DoFpos = 0;
    } else if ( codes [ 0 ].compare("uy") == 0 ) {
        DoFpos = 1;
    } else if ( codes [ 0 ].compare("uz") == 0 && dim > 2 ) {
        DoFpos = 2;
    } else if ( codes [ 0 ].compare("rx") == 0 && dim > 2 ) {
        DoFpos = 3;
    } else if ( codes [ 0 ].compare("ry") == 0 && dim > 2 ) {
        DoFpos = 4;
    } else if ( codes [ 0 ].compare("rz") == 0 ) {
        DoFpos = 5;
        if ( dim == 2 ) {
            DoFpos = 2;
        }
    } else if ( all_of(codes [ 0 ].begin(), codes [ 0 ].end(), :: isdigit) ) {
        DoFpos = atoi(codes [ 0 ].c_str() );
    } else {
        if ( dim == 3 ) {
            cerr << "Error in DoFGauge: only 'ux', 'uy', 'uz', 'rx', 'ry' or 'rz' can be exported by DoFGauge in 3D model" << endl;
            exit(EXIT_FAILURE);
        } else if ( dim == 2 ) {
            cerr << "Error in DoFGauge: only 'ux', 'uy' or 'rz' can be exported by DoFGauge in 2D model" << endl;
            exit(EXIT_FAILURE);
        }
    }

    DoFs.resize( n.size() );
    for ( unsigned i = 0; i < n.size(); i++ ) {
        DoFs [ i ] = nodes->giveNode(n [ i ])->giveStartingDoF() + DoFpos;
    }
    maxsize.resize(1);
    maxsize [ 0 ] = 1;
}


//////////////////////////////////////////////////////////
void DoFGauge :: exportData(unsigned step, const Vector &full_f, const Vector &reactions, fs :: path resultDir) const {
    ( void ) reactions;
    char buffer [ 100 ];
    double value = 0;
    giveFileName(step, buffer);
    ofstream outputfile;
    outputfile.open( ( resultDir / buffer ).string(), ios :: app);
    if ( outputfile.good() ) {
        outputfile << std :: scientific;
        outputfile.precision(precision);
        for ( unsigned i = 0; i < DoFs.size(); i++ ) {
            value += full_f [ DoFs [ i ] ];
        }
        outputfile <<  "\t" << value * multiplier;
    }
    outputfile.close();
}

//////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////
// EXPORT OF IP VALUES
void IntegrationPointGauge :: readFromLine(istringstream &iss) {
    iss >> this->filename;
    iss >> this->name;
    this->codes.resize(1);
    iss >> this->codes [ 0 ];
    unsigned num;
    iss >> num;
    elems.resize(num);
    ipnums.resize(num);
    for ( unsigned i = 0; i < num; i++ ) {
        iss >> elems [ i ];
        iss >> ipnums [ i ];
    }
    DataExporter :: readFromLine(iss);
}

//////////////////////////////////////////////////////////
void IntegrationPointGauge :: init() {
    saveTime_each = 0;
    saveTime_last = 0;

    maxsize.resize(1);
    Vector res;
    Element *ee;
    maxsize [ 0 ] = 0;
    for ( unsigned e = 0; e < elems.size(); e++ ) {
        ee = elemcont->giveElement(e);
        ee->giveIPValues(codes [ 0 ], e, res);
        maxsize [ 0 ] = max< size_t >( maxsize [ 0 ], res.size() );
    }
}


//////////////////////////////////////////////////////////
void IntegrationPointGauge :: exportData(unsigned step, const Vector &full_f, const Vector &reactions, fs :: path resultDir) const {
    ( void ) full_f;
    ( void ) reactions;
    char buffer [ 100 ];
    Vector values;
    giveFileName(step, buffer);
    ofstream outputfile;
    outputfile.open( ( resultDir / buffer ).string(), ios :: app);
    Vector res, sum;
    sum.resize(maxsize [ 0 ]);

    Element *e;
    if ( outputfile.good() ) {
        outputfile << std :: scientific;
        outputfile.precision(precision);
        for ( unsigned i = 0; i < elems.size(); i++ ) {
            e = elemcont->giveElement(elems [ i ]);
            e->giveIPValues(codes [ 0 ], i, res);
            for ( unsigned p = 0; p < min< size_t >( maxsize [ 0 ], res.size() ); p++ ) {
                sum [ p ] += res [ p ];
            }
        }
#if EIGEN_VERSION_AT_LEAST(3, 4, 0)
        for ( auto &p: sum ) {
            outputfile <<  "\t" << p * multiplier;
        }
#else
        for ( long i = 0; i < sum.size(); i++ ) {
            outputfile <<  "\t" << sum [ i ] * multiplier;
        }
#endif
    }
    outputfile.close();
}


//////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////
// EXPORT OF DISPLACEMENTS
/*!
 *  Export displacements (nodeB-nodeA) to txt.
 *  Parameters:
 *  - filename [string] - file to store results (e.g. LD -> results/LD.out)
 *  - name [string] - label placed in table header
 *  - codes [string] - displacement code (ux, uy, uz)
 *  - coordinatesA [float] - coordinates of gauge point A (2D x y) (2D x y z)
 *  - coordinatesB [float] - coordinates of gauge point B (2D x y) (2D x y z)
 *
 *  These parameters can be followed by optional keywords:
 *  - see DataExporter::readFromLine
 */
void DisplacementGauge :: readFromLine(istringstream &iss) {
    iss >> filename;
    iss >> name;
    codes.resize(1);
    iss >> codes [ 0 ];
    double x, y, z;
    if ( dim == 2 ) {
        iss >> x >> y;
        pointA = Point(x, y, 0);
        iss >> x >> y;
        pointB = Point(x, y, 0);
    } else if ( dim == 3 ) {
        iss >> x >> y >> z;
        pointA = Point(x, y, z);
        iss >> x >> y >> z;
        pointB = Point(x, y, z);
    }
    DataExporter :: readFromLine(iss);
}
//////////////////////////////////////////////////////////
void DisplacementGauge :: init() {
    saveTime_each = 0;
    saveTime_last = 0;
    //find element or closest point
    double dist;
    bool foundA = elems->findElementOwningPoint(& elemA, & natCoordsA, & pointA);
    if ( !foundA ) {
        elemA = nullptr;
        nodeA = nodes->findClosestMechanicalNode(pointA, & dist);
    }
    bool foundB = elems->findElementOwningPoint(& elemB, & natCoordsB, & pointB);
    if ( !foundB ) {
        elemB = nullptr;
        nodeB = nodes->findClosestMechanicalNode(pointB, & dist);
    }
    maxsize.resize(1);
    maxsize [ 0 ] = 1;
}

//////////////////////////////////////////////////////////
void DisplacementGauge :: exportData(unsigned step, const Vector &DoFs, const Vector &reactions, fs :: path resultDir) const {
    ( void ) reactions;
    char buffer [ 100 ];
    Vector res;
    double valueA = 0;
    double valueB = 0;
    giveFileName(step, buffer);
    ofstream outputfile;
    outputfile.open( ( resultDir / buffer ).string(), ios :: app);
    if ( outputfile.good() ) {
        outputfile << std :: scientific;
        outputfile.precision(precision);
        if ( elemA ) {
            Vector mv = elemA->giveMasterVariables( & natCoordsA, elemA->giveElemDoFsFromFullDoFs(DoFs) );
            valueA = 0;
            if ( codes [ 0 ].compare("ux") == 0 ) {
                valueA = mv [ 0 ];
            } else if ( dim > 1 && codes [ 0 ].compare("uy") == 0 ) {
                valueA = mv [ 1 ];
            } else if ( dim > 2 && codes [ 0 ].compare("uz") == 0 ) {
                valueA = mv [ 2 ];
            }
        } else {
            nodeA->giveDoFBasedValues(codes [ 0 ], DoFs, res);
            if ( res.size() > 0 ) {
                valueA = res [ 0 ];
            } else {
                valueA = 0;
            }
        }
        if ( elemB ) {
            Vector mv = elemB->giveMasterVariables( & natCoordsB, elemB->giveElemDoFsFromFullDoFs(DoFs) );
            valueB = 0;
            if ( codes [ 0 ].compare("ux") == 0 ) {
                valueB = mv [ 0 ];
            } else if ( dim > 1 && codes [ 0 ].compare("uy") == 0 ) {
                valueB = mv [ 1 ];
            } else if ( dim > 2 && codes [ 0 ].compare("uz") == 0 ) {
                valueB = mv [ 2 ];
            }
        } else {
            nodeB->giveDoFBasedValues(codes [ 0 ], DoFs, res);
            if ( res.size() > 0 ) {
                valueB = res [ 0 ];
            } else {
                valueB = 0;
            }
        }
        outputfile << "\t" << ( valueB - valueA ) * multiplier;
    }
    outputfile.close();
}

//////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////
// EXPORT OF SOLVER VALUES
void SolverGauge :: readFromLine(istringstream &iss) {
    iss >> filename;
    iss >> name;
    codes.resize(1);
    iss >> codes [ 0 ];
    DataExporter :: readFromLine(iss);
}
//////////////////////////////////////////////////////////
void SolverGauge :: init() {
    saveTime_each = 0;
    saveTime_last = 0;

    maxsize.resize(1);
    Vector res;
    solver->giveValues(codes [ 0 ], res);
    maxsize [ 0 ] = res.size();
}

//////////////////////////////////////////////////////////
void SolverGauge :: setSolverPointer(Solver *s) {
    solver = s;
}


//////////////////////////////////////////////////////////
void SolverGauge :: exportData(unsigned step, const Vector &DoFs, const Vector &reactions, fs :: path resultDir) const {
    ( void ) reactions;
    ( void ) DoFs;
    Vector res;
    size_t p;
    char buffer [ 100 ];
    giveFileName(step, buffer);
    ofstream outputfile;
    outputfile.open( ( resultDir / buffer ).string(), ios :: app);
    if ( outputfile.good() ) {
        outputfile << std :: scientific;
        outputfile.precision(precision);
        solver->giveValues(codes [ 0 ], res);
        for ( p = 0; p < min< size_t >( maxsize [ 0 ], res.size() ); p++ ) {
            outputfile << "\t" << res [ p ]*multiplier;
        }
        for ( ; p < maxsize [ 0 ]; p++ ) {
            outputfile << "\t" << 0;
        }
    }
    outputfile.close();
}

//////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////
/*
 * void StructuralExporter :: readFromLine(istringstream &iss) {
 *  iss >> filename;
 *  iss >> name;
 *  codes.resize(1);
 *  iss >> codes [ 0 ];
 *  DataExporter :: readFromLine(iss);
 * }
 * //////////////////////////////////////////////////////////
 * void StructuralExporter :: init() {
 *  time_each = 0;
 *  time_last = 0;
 * }
 *
 * //////////////////////////////////////////////////////////
 * void StructuralExporter :: exportData(unsigned step, const MyVector &DoFs, const MyVector &reactions, fs :: path resultDir) const {
 *  ( void ) DoFs;
 *  ( void ) reactions;
 *  char buffer [ 100 ];
 *  double value = calcValue();
 *  giveFileName(step, buffer);
 *  ofstream outputfile;
 *  outputfile.open( ( resultDir / buffer ).string(), ios :: app );
 *  if ( outputfile.good() ) {
 *      outputfile << std :: scientific;
 *      outputfile.precision(precision);
 *      outputfile << "\t" << value * multiplier;
 *  }
 *  outputfile.close();
 * }
 *
 * double StructuralExporter :: calcValue() const {
 *  double value = 0;
 *  for ( auto const &e : * elems ) {
 *      for ( unsigned i = 0; i < e->giveNumIP(); i++ ) {
 *          value += e->giveIPValue(codes [ 0 ], i);
 *      }
 *  }
 *  return value;
 * }
 */

//////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////
// CONTAINER FOR EXPORTERS
ExporterContainer :: ~ExporterContainer() {
    for ( vector< DataExporter * > :: iterator e = exporters.begin(); e != exporters.end(); ++e ) {
        if ( * e != nullptr ) {
            delete * e;
        }
    }
}


//////////////////////////////////////////////////////////
//!  Container class for Exporters.
/*!
 * Input file keywords:
 * - %TXTNodalExporter - for other parameters see TXTNodalExporter::readFromLine
 * - %TXTElementExporter - for other parameters see TXTElementExporter::readFromLine
 * - %ForceGauge - for other parameters see ForceGauge::readFromLine
 * - %DisplacementGauge - for other parameters see DisplacementGauge::readFromLine
 * - %ValueGauge || %StructuralExporter - for other parameters see StructuralExporter::readFromLine
 * - %DoFGauge - for other parameters see DoFGauge::readFromLine
 * - %TXTGaussPointExporter - for other parameters see TXTGaussPointExporter::readFromLine
 * - %VTKElementExporter - for other parameters see VTKElementExporter::readFromLine
 * - %VTKRBExporter - for other parameters see VTKRB2DExporter::readFromLine
 * - %VTKRCExporter - for other parameters see VTKRCExporter::readFromLine
 */
void ExporterContainer :: readFromFile(const string filename, NodeContainer *n, ElementContainer *e, unsigned dimension) {
    size_t origsize = exporters.size();
    string line, exptype;
    ifstream inputfile( filename.c_str() );
    if ( inputfile.is_open() ) {
        while ( getline(inputfile >> std :: ws, line) ) {
            if ( line.empty() || (line.at(0) == '#') ) {
                continue;
            }
            istringstream iss(line);
            iss >> exptype;

            if ( !( exptype.rfind("#", 0) == 0 ) ) {
                if ( exptype.compare("TXTNodalExporter") == 0 ) {
                    TXTNodalExporter *newexp = new TXTNodalExporter(n, e, dimension);
                    newexp->readFromLine(iss);
                    exporters.push_back(newexp);
                } else if ( exptype.compare("TXTElementExporter") == 0 ) {
                    TXTElementExporter *newexp = new TXTElementExporter(e, dimension);
                    newexp->readFromLine(iss);
                    exporters.push_back(newexp);
                } else if ( exptype.compare("ForceGauge") == 0 ) {
                    ForceGauge *newexp = new ForceGauge(n, dimension);
                    newexp->readFromLine(iss);
                    exporters.push_back(newexp);
                } else if ( exptype.compare("DisplacementGauge") == 0 ) {
                    DisplacementGauge *newexp = new DisplacementGauge(n, e, dimension);
                    newexp->readFromLine(iss);
                    exporters.push_back(newexp);
                    /*} else if ( exptype.compare("ValueGauge") == 0 ||
                     *          exptype.compare("StructuralExporter") == 0 ) {
                     *  StructuralExporter *newexp = new StructuralExporter(n, e, dimension);
                     *  newexp->readFromLine(iss);
                     *  exporters.push_back(newexp);
                     */
                } else if ( exptype.compare("DoFGauge") == 0 ) {
                    DoFGauge *newexp = new DoFGauge(n, dimension);
                    newexp->readFromLine(iss);
                    exporters.push_back(newexp);
                } else if ( exptype.compare("IntegrationPointGauge") == 0 ) {
                    IntegrationPointGauge *newexp = new IntegrationPointGauge(e, dimension);
                    newexp->readFromLine(iss);
                    exporters.push_back(newexp);
                } else if ( exptype.compare("SolverGauge") == 0 ) {
                    SolverGauge *newexp = new SolverGauge(dimension);
                    newexp->readFromLine(iss);
                    exporters.push_back(newexp);
                } else if ( exptype.compare("TXTIntegrationPointExporter") == 0 ) {
                    TXTIntegrationPointExporter *newexp = new TXTIntegrationPointExporter(e, dimension);
                    newexp->readFromLine(iss);
                    exporters.push_back(newexp);
                } else if ( exptype.compare("VTKElementExporter") == 0 ) {
                    VTKElementExporter *newexp = new VTKElementExporter(e, n, dimension);
                    newexp->readFromLine(iss);
                    exporters.push_back(newexp);
                } else if ( exptype.compare("VTKRBExporter") == 0 ) {
                    if ( dimension == 2 ) {
                        VTKRB2DExporter *newexp = new VTKRB2DExporter(e, n, dimension);
                        newexp->readFromLine(iss);
                        exporters.push_back(newexp);
                    } else {
                        std :: cout << "no rigid body exporter for dimension " << dimension << '\n';
                    }
                } else if ( exptype.compare("VTKRCExporter") == 0 ) {
                    if ( dimension == 2 || dimension == 3 ) {
                        VTKRCExporter *newexp = new VTKRCExporter(e, n, dimension);
                        newexp->readFromLine(iss);
                        exporters.push_back(newexp);
                    } else {
                        std :: cout << "no rigid body contacts exporter for dimension " << dimension << '\n';
                    }
                } else if ( exptype.compare("ElementStatsExporter") == 0 ) {
                    ElementStatsExporter *newexp = new ElementStatsExporter(e, dimension);
                    newexp->readFromLine(iss);
                    exporters.push_back(newexp);
                } else {
                    cerr << "Error: Data exporter '" <<  exptype <<  "' is not implemented yet." << endl;
                    exit(EXIT_FAILURE);
                }
            }
        }
        inputfile.close();
        cout << "Input file '" <<  filename << "' succesfully loaded; " << exporters.size() - origsize << " exporters found" << endl;
    } else {
        cerr << "Error ExporterContainer: unable to open input file '" <<  filename <<  "'" << endl;
        exit(EXIT_FAILURE);
    }
}

//////////////////////////////////////////////////////////
void ExporterContainer :: setSolver(Solver *s) {
    for ( auto &exp : exporters ) {
        SolverGauge *sg = dynamic_cast< SolverGauge * >( exp );
        if ( sg ) {
            sg->setSolverPointer(s);
        }
    }
}

//////////////////////////////////////////////////////////
void ExporterContainer :: init(const bool &initial) {
    fs :: create_directories(resultDir);

    bool newname;
    for ( vector< DataExporter * > :: const_iterator d = exporters.begin(); d != exporters.end(); ++d ) {
        ( * d )->init();
        Gauge *g = dynamic_cast< Gauge * >( * d );
        if ( g ) {
            newname = true;
            for ( vector< DataExporter * > :: const_iterator unique = unique_file_exporters.begin(); unique != unique_file_exporters.end(); ++unique ) {
                if ( ( * unique )->giveFileName().compare( ( * d )->giveFileName() ) == 0 ) {
                    newname = false;
                }
            }
            if ( newname ) {
                unique_file_exporters.push_back(* d);
            }
        }
    }

    if ( initial ) {
        //gauge files header
        char buffer [ 100 ];
        for ( vector< DataExporter * > :: const_iterator unique = unique_file_exporters.begin(); unique != unique_file_exporters.end(); ++unique ) {
            ( * unique )->giveFileName(0, buffer);
            ofstream outputfile;
            outputfile.open( ( resultDir / buffer ).string() );
            if ( outputfile.good() ) {
                outputfile << "#step" << "\t" << "time";
            }
            outputfile.close();
        }

        size_t p, maxsize;
        for ( vector< DataExporter * > :: const_iterator d = exporters.begin(); d != exporters.end(); ++d ) {
            Gauge *g = dynamic_cast< Gauge * >( * d );
            if ( g ) {
                g->giveFileName(0, buffer);
                ofstream outputfile;
                outputfile.open( ( resultDir / buffer ).string(), ios :: app);
                if ( outputfile.good() ) {
                    maxsize = g->giveMaxSize(0);
                    if ( maxsize == 1 ) {
                        outputfile << "\t" << g->giveName();
                    } else {
                        for ( p = 0; p < maxsize; p++ ) {
                            outputfile << "\t" << g->giveName() << "_" << p;
                        }
                    }
                }
                outputfile.close();
            }
        }

        for ( vector< DataExporter * > :: const_iterator unique = unique_file_exporters.begin(); unique != unique_file_exporters.end(); ++unique ) {
            ( * unique )->giveFileName(0, buffer);
            ofstream outputfile;
            outputfile.open( ( resultDir / buffer ).string(), ios :: app);
            if ( outputfile.good() ) {
                outputfile << endl;
            }
            outputfile.close();
        }
    }
};


void ExporterContainer :: clear() {
    for ( vector< DataExporter * > :: iterator e = exporters.begin(); e != exporters.end(); ++e ) {
        if ( * e != nullptr ) {
            delete * e;
        }
    }
}

//////////////////////////////////////////////////////////


//////////////////////////////////////////////////////////
void ExporterContainer :: exportData(unsigned step, double time, const Vector &DoFs, const Vector &reactions, const bool &exportAll) const {
    //add step number to gauge exporter files
    char buffer [ 100 ];
    for ( vector< DataExporter * > :: const_iterator unique = unique_file_exporters.begin(); unique != unique_file_exporters.end(); ++unique ) {
        ( * unique )->giveFileName(0, buffer);
        ofstream outputfile;
        outputfile.open( ( resultDir / buffer ).string(), ios :: app);
        if ( outputfile.good() ) {
            outputfile << std :: scientific;
            outputfile << step << "\t" << time;
        }
        outputfile.close();
    }

    // export
    for ( vector< DataExporter * > :: const_iterator d = exporters.begin(); d != exporters.end(); ++d ) {
        if ( ( * d )->doExportNow(time, step) || exportAll ) {
            ( * d )->exportData(step, DoFs, reactions, resultDir);
        }
    }

    // add end line to gauge exporter files
    for ( vector< DataExporter * > :: const_iterator unique = unique_file_exporters.begin(); unique != unique_file_exporters.end(); ++unique ) {
        ( * unique )->giveFileName(0, buffer);
        ofstream outputfile;
        outputfile.open( ( resultDir / buffer ).string(), ios :: app);
        if ( outputfile.good() ) {
            outputfile << endl;
        }
        outputfile.close();
    }
};

//////////////////////////////////////////////////////////
void ExporterContainer :: appendToAllNames(string app) {
    for ( auto &d : exporters ) {
        d->appendToName(app);
    }
}

//////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////
void ExportAllElementsNodalStress(std :: vector< Matrix > &stress, const Vector &DoFs, const Vector &reactions, const NodeContainer *nodes, const ElementContainer *elems, const unsigned &dim) {
    // MyVector stressXYZ, stress_zero;
    // stress_zero = MyVector((double)0, dim);
    ( void ) reactions;

    unsigned node_id, ni;
    double first;
    Vector intF0 = Vector :: Zero(2 * dim);
    Vector intF1 = Vector :: Zero(dim);
    Vector intF2 = Vector :: Zero(dim);
    Vector elDoFvalues, strainNT;
    vector< unsigned >elDoFs;

    vector< double >Volume(stress.size(), 0);
    double single_volume;
    RigidBodyContact *rbc;

    for ( auto const &el : * elems ) {
        // use only elements that are derived from LTCBEAM
        // but dynamic cast costs a lot and this is performed every step (in adaptivity), or on every model vtk save
        // rfind overload to compare start of a string
        if ( el->giveName().rfind("LTCB", 0) == 0 ) {
            rbc = static_cast< RigidBodyContact * >( el );
            elDoFs = el->giveDoFs();
            elDoFvalues.resize( elDoFs.size() );
            for ( unsigned i = 0; i < elDoFs.size(); i++ ) {
                elDoFvalues [ i ] = DoFs [ elDoFs [ i ] ];
            }

            // to each node correspond 0.5 of volume
            // TODO must be repaired for power tessellation
            single_volume = 0.5 * rbc->giveLength() * rbc->giveArea() / dim;
            first = 1;
            ni = 0;
            for ( auto const &n : el->giveNodes() ) {
                auto res = std :: find(begin(* nodes), end(* nodes), n);
                node_id = std :: distance(begin(* nodes), res);
                Volume [ node_id ] += single_volume;

                stress [ node_id ] += dyadicProduct(
                    (
                        rbc->giveContactStressXYZ()
                        * rbc->giveArea() // vyhodit
                    )
                    * first
                    , rbc->giveVectorToNode(ni, 0) ); // tady může být jen poloha  IP (HonzaE článek 2020)
                // TODO here is probably missing some value corresponding to internal moments
                // for node corresponding to end of element, traction needs to be reversed
                first = -1;
                ni++;
            }
        } else {
            //
        }
    }
    for ( unsigned si = 0; si < stress.size(); si++ ) {
        if ( Volume [ si ] <= 0 ) {
            stress [ si ].setZero();
        } else {
            stress [ si ] /= Volume [ si ];
        }
    }
}


void saveNodes(const NodeContainer &nodes, const std :: vector< std :: string > &NodeTypes, fs :: path resultDir) {
    // if NodeTypes.empty() then save all nodes
    // TODO finish this, now (for adaptivity) just save path to file with particles
    ( void ) nodes;
    ( void ) NodeTypes;
    ( void ) resultDir;
}

void saveElems(const ElementContainer &elems, const std :: vector< std :: string > &ElemTypes, fs :: path resultDir) {
    // if ElemTypes.empty() then save all elems
    ( void ) elems;
    ( void ) ElemTypes;
    ( void ) resultDir;
}
