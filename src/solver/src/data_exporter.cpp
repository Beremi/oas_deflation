#include "data_exporter.h"
#include "vtk_exporter.h"

//////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////
// BASIC DATA EXPORTER - master class
void DataExporter :: giveFileName(unsigned step, char *buffer) const {
    sprintf(buffer, "%s_%05d.out", filename.c_str(), step);
}

//////////////////////////////////////////////////////////
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
    if ( outputfile.is_open() ) {
        outputfile << std :: scientific;
        outputfile.precision(precision);
        for ( unsigned n = 0; n < nodes->giveSize(); n++ ) {
            nn = nodes->giveNode(n);
            for ( vector< string > :: const_iterator c = codes.begin(); c != codes.end(); ++c ) {
                value = nn->giveDoFBasedValue(* c, DoFs);
                outputfile << value;
                if ( c != codes.end() - 1 ) {
                    outputfile << "\t";
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
//////////////////////////////////////////////////////////
// EXPORT FROM GAUSS POINTS TO TXT
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
            nIP = ee->giveIPNum();
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
    sprintf(buffer, "%s.out", filename.c_str() );
}


//////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////
// EXPORT OF FORCES
void ForceGauge :: readFromLine(istringstream &iss) {
    iss >> filename;
    iss >> name;
    codes.resize(1);
    iss >> codes [ 0 ];
    unsigned num;
    iss >> num;
    n.resize(num);
    for ( unsigned i = 0; i < num; i++ ) {
        iss >> n [ i ];
    }
    DataExporter :: readFromLine(iss);
}

//////////////////////////////////////////////////////////
ForceGauge :: ForceGauge(string &f, string &gname, vector< string > &c, vector< unsigned > &nn, NodeContainer *nc, double m, unsigned dimension) : Gauge(dimension) {
    nodes = nc;
    filename = f;
    name = gname;
    n = nn;
    codes = c;
    multiplier = m;
}

//////////////////////////////////////////////////////////
void ForceGauge :: init() {
    time_each = 0;
    time_last = 0;
    unsigned DoFpos = 0;
    if ( codes [ 0 ].compare("fx") == 0 ) {
        DoFpos = 0;
    } else if ( codes [ 0 ].compare("fy") == 0 ) {
        DoFpos = 1;
    } else if ( codes [ 0 ].compare("fz") == 0 && dim > 2 ) {
        DoFpos = 2;
    } else if ( codes [ 0 ].compare("mx") == 0 && dim > 2 ) {
        DoFpos = 3;
    } else if ( codes [ 0 ].compare("my") == 0 && dim > 2 ) {
        DoFpos = 4;
    } else if ( codes [ 0 ].compare("mz") == 0 ) {
        DoFpos = 5;
        if ( dim == 2 ) {
            DoFpos = 2;
        }
    } else {
        if ( dim == 3 ) {
            cerr << "Error in ForceGauge: only 'fx', 'fy', 'fz', 'mx', 'my' or 'mz' can be exported by ForceGauge in 3D model" << endl;
            exit(EXIT_FAILURE);
        } else if ( dim == 2 ) {
            cerr << "Error in ForceGauge: only 'fx', 'fy' or 'mz' can be exported by ForceGauge in 2D model" << endl;
            exit(EXIT_FAILURE);
        }
    }

    DoFs.resize(n.size() );
    for ( unsigned i = 0; i < n.size(); i++ ) {
        DoFs [ i ] = nodes->giveNode(n [ i ])->giveStartingDoF() + DoFpos;
    }
}


//////////////////////////////////////////////////////////
void ForceGauge :: exportData(unsigned step, const Vector &full_f, const Vector &reactions, fs :: path resultDir) const {
    ( void ) full_f;
    char buffer [ 100 ];
    double value = 0;
    giveFileName(step, buffer);
    ofstream outputfile;
    outputfile.open( ( resultDir / buffer ).string(), ios :: app );
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
void DoFGauge :: readFromLine(istringstream &iss) {
    iss >> filename;
    iss >> name;
    iss >> nodenum;
    iss >> direction;
    DataExporter :: readFromLine(iss);
}

//////////////////////////////////////////////////////////
DoFGauge :: DoFGauge(string &f, string &gname, unsigned n, unsigned dir, NodeContainer *nn, double m, unsigned dimension) : Gauge(dimension) {
    filename = f;
    name = gname;
    direction = dir;
    nodes = nn;
    nodenum = n;
    multiplier = m;
}

//////////////////////////////////////////////////////////
void DoFGauge :: init() {
    n = nodes->giveNode(nodenum);
    DoF = n->giveStartingDoF() + direction;
    time_each = 0;
    time_last = 0;
}


//////////////////////////////////////////////////////////
void DoFGauge :: exportData(unsigned step, const Vector &full_f, const Vector &reactions, fs :: path resultDir) const {
    ( void ) step;
    ( void ) reactions;
    char buffer [ 100 ];
    giveFileName(step, buffer);
    ofstream outputfile;
    outputfile.open( ( resultDir / buffer ).string(), ios :: app );
    if ( outputfile.good() ) {
        outputfile << std :: scientific;
        outputfile.precision(precision);
        outputfile << "\t" <<  full_f [ DoF ] * multiplier;
    }
    outputfile.close();
}

//////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////
// EXPORT OF DISPLACEMENTS
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
    //find closest point
    //TODO: better to use direct displacements at given points A and B
    nodeA = nodes->findClosestMechanicalNode(pointA);
    nodeB = nodes->findClosestMechanicalNode(pointB);
}

//////////////////////////////////////////////////////////
void DisplacementGauge :: exportData(unsigned step, const Vector &DoFs, const Vector &reactions, fs :: path resultDir) const {
    ( void ) reactions;
    char buffer [ 100 ];
    double value = 0;
    giveFileName(step, buffer);
    ofstream outputfile;
    outputfile.open( ( resultDir / buffer ).string(), ios :: app );
    if ( outputfile.good() ) {
        outputfile << std :: scientific;
        outputfile.precision(precision);
        value = nodeB->giveDoFBasedValue(codes [ 0 ], DoFs) - nodeA->giveDoFBasedValue(codes [ 0 ], DoFs);
        outputfile << "\t" << value * multiplier;
    }
    outputfile.close();
}

//////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////
// EXPORT OF DISPLACEMENTS
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
    outputfile.open( ( resultDir / buffer ).string(), ios :: app );
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
        for ( unsigned i = 0; i < e->giveIPNum(); i++ ) {
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
        delete * e;
    }
}

//////////////////////////////////////////////////////////
void ExporterContainer :: readFromFile(const string filename, NodeContainer *n, ElementContainer *e, unsigned dimension) {
    size_t origsize = exporters.size();
    string line, exptype;
    ifstream inputfile(filename.c_str() );
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
            if ( !exptype.rfind("#", 0) == 0 ) {
                if ( exptype.compare("TXTNodalExporter") == 0 ) {
                    TXTNodalExporter *newexp = new TXTNodalExporter(n, dimension);
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
                } else if ( exptype.compare("TXTGaussPointExporter") == 0 ) {
                    TXTGaussPointExporter *newexp = new TXTGaussPointExporter(e, dimension);
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
void ExporterContainer :: init() {

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
            outputfile.open( ( resultDir / buffer ).string(), ios :: app );
            if ( outputfile.good() ) {
                outputfile << "\t" << g->giveName();
            }
            outputfile.close();
        }
    }

    for ( vector< DataExporter * > :: const_iterator unique = unique_file_exporters.begin(); unique != unique_file_exporters.end(); ++unique ) {
        ( * unique )->giveFileName(0, buffer);
        ofstream outputfile;
        outputfile.open( ( resultDir / buffer ).string(), ios :: app );
        if ( outputfile.good() ) {
            outputfile << endl;
        }
        outputfile.close();
    }
};

//////////////////////////////////////////////////////////
void ExporterContainer :: exportData(unsigned step, double time, const Vector &DoFs, const Vector &reactions, const bool &exportAll) const {
    //add step number to gauge exporter files
    char buffer [ 100 ];
    for ( vector< DataExporter * > :: const_iterator unique = unique_file_exporters.begin(); unique != unique_file_exporters.end(); ++unique ) {
        ( * unique )->giveFileName(0, buffer);
        ofstream outputfile;
        outputfile.open( ( resultDir / buffer ).string(), ios :: app );
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
        outputfile.open( ( resultDir / buffer ).string(), ios :: app );
        if ( outputfile.good() ) {
            outputfile << endl;
        }
        outputfile.close();
    }
};

void ExportAllElementsNodalStress(std::vector< Matrix > &stress, const Vector &DoFs, const NodeContainer *nodes, const ElementContainer *elems, const unsigned &dim) {

  // Vector stressXYZ, stress_zero;
  // stress_zero = Vector((double)0, dim);

  unsigned node_id, ni;
  int first;
  Vector elDoFvalues, strainNT;
  vector< unsigned >elDoFs;

  vector < double > Volume( stress.size(), 0);
  double single_volume;
  RigidBodyContact *rbc;

  for ( auto const & el : *elems) {
    if ( el->giveName().compare("RigidBodyContact") == 0 ) {
      rbc = static_cast< RigidBodyContact * >( el );
      elDoFs = el->giveDoFs();
      elDoFvalues.resize(elDoFs.size() );
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

          stress [ node_id ] += dyadicProduct(el->giveInternalForces(elDoFvalues, false) * rbc->giveArea() * first, rbc->giveDistanceToNode(ni, 0));
          // for node corresponding to end of element, traction needs to be reversed
          first = -1;
          ni++;

      }


    } else {
        //
    }
  }
  for ( unsigned si = 0; si < stress.size(); si++){
    if ( Volume [ si ] <= 0 ) stress [ si ] *= 0;
    else stress [ si ] /= Volume [ si ];
  }
}
