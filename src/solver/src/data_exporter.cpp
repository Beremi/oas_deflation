#include "data_exporter.h"
#include "vtk_exporter.h"
#include "exporter_model.h"
#include "geometry.h"
#include "element_discrete.h"

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
 *  - timeEach [float] - time each (TODO: better description)
 *  - time_last [float] - time shift (TODO: better description)
 *  - precision [float] - precision of stored values
 */
void DataExporter :: readFromLine(istringstream &iss) {
    iss.clear(); // clear string stream
    iss.seekg(0, iss.beg); //reset position in string stream
    string param;
    // initiate variables in case they are not specified
    time_each = 0;
    time_last = 0;
    while ( !iss.eof() ) {
        iss >> param;
        if ( param.compare("saveEvery") == 0 || param.compare("timeEach") == 0 ) {
            iss >> time_each;
        } else if ( param.compare("timeShift") == 0 ) {
            iss >> time_last;
        } else if ( param.compare("precision") == 0 ) {
            iss >> precision;
        }
    }
}

//////////////////////////////////////////////////////////
bool DataExporter :: doExportNow(const double &time) {
    if ( time < time_last + time_each - 1e-12 ) {
        return false;
    } else {
        time_last = time;
        return true;
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
void TXTNodalExporter :: exportData(unsigned step, const Vector &DoFs, const Vector &reactions, fs :: path resultDir) const {
    ( void ) reactions;
    char buffer [ 100 ];
    Node *nn;
    double value;
    giveFileName(step, buffer);
    ofstream outputfile( ( resultDir / buffer ).string() );

    // first export nodal_stress if wanted
    string nds = "nodal_stress";
    vector< Matrix >nodal_stress;
    if ( std :: find(codes.begin(), codes.end(), nds) != codes.end() ) {
        // reserve space only if nodal stresses should be exported
        nodal_stress.resize( nodes->giveSize(), Matrix(this->dim, this->dim) );
        // export nodal stresses:
        ExportAllElementsNodalStress(nodal_stress, DoFs, reactions, this->nodes, this->elems, this->dim);
    }

    if ( outputfile.is_open() ) {
        outputfile << std :: scientific;
        outputfile.precision(precision);
        unsigned i = 0;
        for ( unsigned n = 0; n < nodes->giveSize(); n++, i++ ) {
            nn = nodes->giveNode(n);
            for ( vector< string > :: const_iterator c = codes.begin(); c != codes.end(); ++c ) {
                if ( c->compare(nds) == 0 ) {
                    for ( auto const &d : MatrixToStdVectForParaview(nodal_stress [ i ], dim) ) {
                        outputfile << d << '\t';
                    }
                } else {
                    value = nn->giveDoFBasedValue(* c, DoFs);
                    outputfile << value;
                    if ( c != codes.end() - 1 ) {
                        outputfile << "\t";
                    }
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
void TXTElementExporter :: exportData(unsigned step, const Vector &DoFs, const Vector &reactions, fs :: path resultDir) const {
    ( void ) DoFs;
    ( void ) reactions;
    char buffer [ 100 ];
    Element *ee;
    double value;
    size_t nIP;
    giveFileName(step, buffer);
    ofstream outputfile( ( resultDir / buffer ).string() );

    if ( outputfile.is_open() ) {
        outputfile << std :: scientific;
        outputfile.precision(precision);
        for ( unsigned e = 0; e < elems->giveSize(); e++ ) {
            ee = elems->giveElement(e);
            nIP = ee->giveNumIP();
            for ( unsigned k = 0; k < nIP; k++ ) {
                for ( vector< string > :: const_iterator c = codes.begin(); c != codes.end(); ++c ) {
                    value = ee->giveIPValue(* c, k);
                    outputfile << value;
                    if ( c != codes.end() - 1 ) {
                        outputfile << "\t";
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
// EXPORT FROM GAUSS POINTS TO TXT
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
void TXTGaussPointExporter :: readFromLine(istringstream &iss) {
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
void TXTGaussPointExporter :: exportData(unsigned step, const Vector &DoFs, const Vector &reactions, fs :: path resultDir) const {
    ( void ) DoFs;
    ( void ) reactions;
    char buffer [ 100 ];
    Element *ee;
    double value;
    size_t nIP;
    giveFileName(step, buffer);
    ofstream outputfile( ( resultDir / buffer ).string() );

    if ( outputfile.is_open() ) {
        outputfile << std :: scientific;
        outputfile.precision(precision);
        for ( unsigned e = 0; e < elems->giveSize(); e++ ) {
            ee = elems->giveElement(e);
            nIP = ee->giveNumIP();
            for ( unsigned k = 0; k < nIP; k++ ) {
                for ( vector< string > :: const_iterator c = codes.begin(); c != codes.end(); ++c ) {
                    value = ee->giveIPValue(* c, k);
                    outputfile << value;
                    if ( c != codes.end() - 1 ) {
                        outputfile << "\t";
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
    time_each = 0;
    time_last = 0;
    DoFs.resize( n.size() );
    for ( unsigned i = 0; i < n.size(); i++ ) {
        DoFs [ i ] = nodes->giveNode(n [ i ])->giveStartingDoF() + nodes->giveNode(n [ i ])->giveOrderOfForceCode(codes [ 0 ]);
    }
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
    time_each = 0;
    time_last = 0;
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
    } else if ( all_of(codes [ 0 ].begin(), codes [ 0 ].end(), ::isdigit) ) {
        DoFpos = atoi(codes [ 0 ].c_str()); 
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
void IPGauge :: readFromLine(istringstream &iss) {
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
void IPGauge :: init() {
    time_each = 0;
    time_last = 0;
}


//////////////////////////////////////////////////////////
void IPGauge :: exportData(unsigned step, const Vector &full_f, const Vector &reactions, fs :: path resultDir) const {
    ( void ) full_f;
    ( void ) reactions;
    char buffer [ 100 ];
    double value = 0;
    giveFileName(step, buffer);
    ofstream outputfile;
    outputfile.open( ( resultDir / buffer ).string(), ios :: app);
    if ( outputfile.good() ) {
        outputfile << std :: scientific;
        outputfile.precision(precision);
        for ( unsigned i = 0; i < elems.size(); i++ ) {
            value += elemcont->giveElement(elems [ i ])->giveMatStatus(ipnums [ i ])->giveValue(codes [ 0 ]);
        }
        outputfile <<  "\t" << value * multiplier;
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
        pointA = Point(x, y);
        iss >> x >> y;
        pointB = Point(x, y);
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
    time_each = 0;
    time_last = 0;
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
}

//////////////////////////////////////////////////////////
void DisplacementGauge :: exportData(unsigned step, const Vector &DoFs, const Vector &reactions, fs :: path resultDir) const {
    ( void ) reactions;
    char buffer [ 100 ];
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
            } else if ( dim > 2 && codes [ 0 ].compare("uz") == 0 )                                                                      {
                valueA = mv [ 2 ];
            }
        } else  {
            valueA = nodeA->giveDoFBasedValue(codes [ 0 ], DoFs);
        }
        if ( elemB ) {
            Vector mv = elemB->giveMasterVariables( & natCoordsB, elemB->giveElemDoFsFromFullDoFs(DoFs) );
            valueB = 0;
            if ( codes [ 0 ].compare("ux") == 0 ) {
                valueB = mv [ 0 ];
            } else if ( dim > 1 && codes [ 0 ].compare("uy") == 0 ) {
                valueB = mv [ 1 ];
            } else if ( dim > 2 && codes [ 0 ].compare("uz") == 0 )                                                                      {
                valueB = mv [ 2 ];
            }
        } else  {
            valueB = nodeB->giveDoFBasedValue(codes [ 0 ], DoFs);
        }
        outputfile << "\t" << ( valueB - valueA ) * multiplier;
    }
    outputfile.close();
}

//////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////
// EXPORT OF DISPLACEMENTS
/*!
 *  Export displacements of the structure to txt.
 *  Parameters:
 *  - filename [string] - file to store results (e.g. LD -> results/LD.out)
 *  - name [string]- label placed in table header
 *  - codes [string] - displacement code (one of ux, uy, uz)
 *
 *  These parameters can be followed by optional keywords:
 *  - see DataExporter::readFromLine
 */
void StructuralExporter :: readFromLine(istringstream &iss) {
    iss >> filename;
    iss >> name;
    codes.resize(1);
    iss >> codes [ 0 ];
    DataExporter :: readFromLine(iss);
}
//////////////////////////////////////////////////////////
void StructuralExporter :: init() {
    time_each = 0;
    time_last = 0;
}

//////////////////////////////////////////////////////////
void StructuralExporter :: exportData(unsigned step, const Vector &DoFs, const Vector &reactions, fs :: path resultDir) const {
    ( void ) DoFs;
    ( void ) reactions;
    char buffer [ 100 ];
    double value = calcValue();
    giveFileName(step, buffer);
    ofstream outputfile;
    outputfile.open( ( resultDir / buffer ).string(), ios :: app);
    if ( outputfile.good() ) {
        outputfile << std :: scientific;
        outputfile.precision(precision);
        outputfile << "\t" << value * multiplier;
    }
    outputfile.close();
}

double StructuralExporter :: calcValue() const {
    double value = 0;
    for ( auto const &e : * elems ) {
        for ( unsigned i = 0; i < e->giveNumIP(); i++ ) {
            value += e->giveIPValue(codes [ 0 ], i);
        }
    }
    return value;
}

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
        while ( getline(inputfile, line) ) {
            if ( line.empty() ) {
                continue;
            }
            if ( line.at(0) == '#' ) {
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
                } else if ( exptype.compare("ValueGauge") == 0 ||
                            exptype.compare("StructuralExporter") == 0 ) {
                    StructuralExporter *newexp = new StructuralExporter(n, e, dimension);
                    newexp->readFromLine(iss);
                    exporters.push_back(newexp);
                } else if ( exptype.compare("DoFGauge") == 0 ) {
                    DoFGauge *newexp = new DoFGauge(n, dimension);
                    newexp->readFromLine(iss);
                    exporters.push_back(newexp);
                } else if ( exptype.compare("IPGauge") == 0 ) {
                    IPGauge *newexp = new IPGauge(e, dimension);
                    newexp->readFromLine(iss);
                    exporters.push_back(newexp);
                } else if ( exptype.compare("TXTGaussPointExporter") == 0 ) {
                    TXTGaussPointExporter *newexp = new TXTGaussPointExporter(e, dimension);
                    newexp->readFromLine(iss);
                    exporters.push_back(newexp);
                } else if ( exptype.compare("VTKElementExporter") == 0 ) {
                    VTKElementExporter *newexp = new VTKElementExporter(e, n, dimension);
                    newexp->readFromLine(iss);
                    exporters.push_back(newexp);
                } else if ( exptype.compare("VTKElement2Exporter") == 0 ) {
                    VTKElement2Exporter *newexp = new VTKElement2Exporter(e, n, dimension);
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

        for ( vector< DataExporter * > :: const_iterator d = exporters.begin(); d != exporters.end(); ++d ) {
            Gauge *g = dynamic_cast< Gauge * >( * d );
            if ( g ) {
                ( * d )->giveFileName(0, buffer);
                ofstream outputfile;
                outputfile.open( ( resultDir / buffer ).string(), ios :: app);
                if ( outputfile.good() ) {
                    outputfile << "\t" << g->giveName();
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
            outputfile << step << "\t" << time;
        }
        outputfile.close();
    }

    // export
    for ( vector< DataExporter * > :: const_iterator d = exporters.begin(); d != exporters.end(); ++d ) {
        if ( ( * d )->doExportNow(time) || exportAll ) {
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
    // Vector stressXYZ, stress_zero;
    // stress_zero = Vector((double)0, dim);
    ( void ) reactions;

    unsigned node_id, ni;
    double first;
    Vector intF0(2 * dim);
    Vector intF1(dim);
    Vector intF2(dim);
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
            stress [ si ] *= 0;
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
