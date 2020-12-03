#include "vtk_exporter.h"
#include "misc.h"

//////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////
// VTK EXPORTERS

void VTKExporter :: giveFileName(unsigned step, char *buffer) const {
    sprintf(buffer, "%s_%05d.vtu", filename.c_str(), step);
}

//////////////////////////////////////////////////////////
void VTKExporter :: readFromLine(istringstream &iss) {
    string param;
    unsigned num;
    vector< string >cell_data, point_data;
    cell_data.resize(0);
    point_data.resize(0);
    iss >> filename;
    while ( !iss.eof() ) {
        iss >> param;
        if ( param.compare("cellData") == 0 ) {
            iss >> num;
            cell_data.resize(num);
            for ( unsigned i = 0; i < num; i++ ) {
                iss >> cell_data [ i ];
            }
        } else if ( param.compare("pointData") == 0 ) {
            iss >> num;
            point_data.resize(num);
            for ( unsigned i = 0; i < num; i++ ) {
                iss >> point_data [ i ];
            }
        }
    }
    codes.resize(cell_data.size() + point_data.size() );
    num = 0;
    for ( auto const &cel : cell_data ) {
        codes [ num++ ] = cel;
    }
    cell_data_size = cell_data.size();
    for ( auto const &po : point_data ) {
        codes [ num++ ] = po;
    }
    DataExporter :: readFromLine(iss);
}

bool isAddonCellScalarData(const string &var) {
    std :: vector< string >addon_list = { "normal_strain", "strainN",
                                          "sliding_strain", "strainT",
                                          "strainTY", "strainTZ",
                                          "strainT1", "strainT2",
                                          "crack_opening", "crack_sliding",
                                          "stressN", "stressT",
                                          "stressTY", "stressTZ"
                                          "stressT2", "stressT2"
                                          // "strainXYZ", "stressXYZ"
    };

    return std :: find(addon_list.begin(), addon_list.end(), var) != addon_list.end();
}

bool isAddonCellVectorialData(const string &var) {
    std :: vector< string >addon_list = { "strainXYZ", "stressXYZ" };

    return std :: find(addon_list.begin(), addon_list.end(), var) != addon_list.end();
}

bool isAddonPointVectorialData(const string &var) {
    std :: vector< string >addon_list = { "nodal_stress" };

    return std :: find(addon_list.begin(), addon_list.end(), var) != addon_list.end();
}

void exportAddonVectorialCellData(const unsigned &dim, const ElementContainer *elems, const Vector &DoFs, const vector< string > &codes, vector< unsigned > &indeces, vector< vector< Vector > > &cell_vect_data, bool doubled = false, bool rbcOnly = false) {
    Vector elDoFvalues, strainNT;
    vector< unsigned >elDoFs;

    Vector vect_ini(dim);
    // set all elements of a vector to zero
    for ( auto &a : vect_ini ) {
        a = 0;
    }
    Vector vect_data = vect_ini;
    RigidBodyContact *rbc;

    for ( auto const &e : * elems ) {
        rbc = nullptr;
        // NOTE do not use this for transport elements
        if ( e->giveName().compare("RigidBodyContact") == 0 ) {
            elDoFs = e->giveDoFs();
            elDoFvalues.resize(elDoFs.size() );
            for ( unsigned i = 0; i < elDoFs.size(); i++ ) {
                elDoFvalues [ i ] = DoFs [ elDoFs [ i ] ];
            }
            rbc = static_cast< RigidBodyContact * >( e );

            strainNT = rbc->giveContactStrainNT(elDoFvalues);
            for ( unsigned i = 0; i < indeces.size(); i++ ) {
                if ( codes [ indeces [ i ] ].compare("strainXYZ") == 0 ) {
                    vect_data = rbc->giveContactStrainXYZ(elDoFvalues);
                    // vect_data = matrix_vector_multiply(rbc->giveRMatrix().transpose(),strainNT);
                } else if ( codes [ indeces [ i ] ].compare("stressXYZ") == 0 ) {
                    vect_data = rbc->giveInternalForces(elDoFvalues, false) / rbc->giveArea();
                } else {
                    vect_data = vect_ini;
                }
                cell_vect_data [ i ].push_back(vect_data);
                if ( doubled ) {
                    cell_vect_data [ i ].push_back(vect_data);
                }
            }
        } else {
            // for any other element, only one number will be stored
            if ( !rbcOnly ) {
                for ( unsigned i = 0; i < indeces.size(); i++ ) {
                    cell_vect_data [ i ].push_back(vect_ini);
                }
            }
        }
    }
}


void exportAddonScalarCellData(const ElementContainer *elems, const Vector &DoFs,
                               const vector< bool > &codes_positions, const vector< string > &codes,
                               vector< vector< double > > &cell_data, bool doubled = false, bool rbcOnly = false) {
    Vector elDoFvalues, strainNT;
    vector< unsigned >elDoFs;
    double data;
    RigidBodyContact *rbc;

    for ( auto const &e : * elems ) {
        rbc = nullptr;
        // NOTE do not use this for transport elements
        if ( e->giveName().compare("RigidBodyContact") == 0 ) {
            elDoFs = e->giveDoFs();
            elDoFvalues.resize(elDoFs.size() );
            for ( unsigned i = 0; i < elDoFs.size(); i++ ) {
                elDoFvalues [ i ] = DoFs [ elDoFs [ i ] ];
            }
            rbc = static_cast< RigidBodyContact * >( e );
            strainNT = rbc->giveContactStrainNT(elDoFvalues);
            for ( unsigned i = 0; i < codes.size(); i++ ) {
                if ( codes_positions [ i ] ) {
                    if ( codes [ i ].compare("strainN") == 0 ) {
                        data = strainNT [ 0 ];
                    } else if ( codes [ i ].compare("strainT") == 0 ) {
                        double strT = 0;
                        for ( unsigned j = 1; j < strainNT.size(); j++ ) {
                            strT += pow(strainNT [ j ], 2);
                        }
                        data = sqrt(strT);
                    } else if ( codes [ i ].compare("strainTY") == 0 ||
                                codes [ i ].compare("strainT1") == 0 ) {
                        data = strainNT [ 1 ];
                    } else if ( codes [ i ].compare("strainTZ") == 0 ||
                                codes [ i ].compare("strainT2") == 0 ) {
                        if ( strainNT.size() > 2 ) {
                            data = strainNT [ 2 ];
                        } else {
                            data = 0;
                        }
                    } else if ( codes [ i ].rfind("stress", 0) == 0 ) {
                        MaterialStatus *stats = static_cast< MaterialStatus * >( rbc->giveMaterialStats() [ 0 ] );
                        Vector stressNT = stats->giveStress(strainNT);
                        if ( codes [ i ].compare("stressN") == 0 ) {
                            data = stressNT [ 0 ];
                        } else if ( codes [ i ].compare("stressT") == 0 ) {
                            double strT = 0;
                            for ( unsigned j = 1; j < stressNT.size(); j++ ) {
                                strT += pow(stressNT [ j ], 2);
                            }
                            data = sqrt(strT);
                        } else if ( codes [ i ].compare("stressTY") == 0 ||
                                    codes [ i ].compare("stressT1") == 0 ) {
                            data = stressNT [ 1 ];
                        } else if ( codes [ i ].compare("stressTZ") == 0 ||
                                    codes [ i ].compare("stressT2") == 0 ) {
                            if ( stressNT.size() > 2 ) {
                                data = stressNT [ 2 ];
                            } else {
                                data = 0;
                            }
                        }
                    } else if ( codes [ i ].compare("crack_opening") == 0 ) {
                        data = strainNT [ 0 ] * rbc->giveLength() * rbc->giveIPValue("damageN", 0);
                    } else if ( codes [ i ].compare("crack_sliding") == 0 ) {
                        double strT = 0;
                        for ( unsigned j = 1; j < strainNT.size(); j++ ) {
                            strT += pow(strainNT [ j ], 2);
                        }
                        data = sqrt(strT) * rbc->giveLength() * rbc->giveIPValue("damageT", 0);
                    } else {
                        data = 0;
                    }
                    cell_data [ i ].push_back(data);
                    if ( doubled ) {
                        cell_data [ i ].push_back(data);
                    }
                }
            }
        } else {
            // for any other element, only one number will be stored
            if ( !rbcOnly ) {
                for ( unsigned i = 0; i < codes.size(); i++ ) {
                    if ( codes_positions [ i ] ) {
                        cell_data [ i ].push_back(0);
                    }
                }
            }
        }
    }
}


vector< double >MatrixToStdVectForParaview(const Matrix &s, const unsigned &dim) {
    vector< double >data;
    if ( dim == 2 ) {
        data.resize(3);
    } else {
        data.resize(6);  // NOTE other case than 2D or 3D not considered
    }
    for ( unsigned i = 0; i < s.numRows(); i++ ) {
        data [ i ] = s [ i ] [ i ];
    }
    data [ ( dim - 1 ) * 3 - 1 ] = 0.5 * ( s [ dim - 1 ] [ 0 ] + s [ 0 ] [ dim - 1 ] );
    if ( dim > 2 ) {
        data [ 3 ] = 0.5 * ( s [ 0 ] [ 1 ] + s [ 1 ] [ 0 ] );
        data [ 4 ] = 0.5 * ( s [ 2 ] [ 1 ] + s [ 1 ] [ 2 ] );
    }
    return data;
}


//////////////////////////////////////////////////////////
// ELEMENTS TO VTU FILE
//////////////////////////////////////////////////////////
void VTKElementExporter :: exportData(unsigned step, const Vector &DoFs, const Vector &reactions, fs :: path resultDir) const {
    // Export of elements into vtu xml file format (vtu = vtk for unstructured grid)
    // NOTE this is messy construction of xml file, will be remade using some of xml libraries for cpp
    char buffer [ 100 ];
    // Point P;
    // Element *ee;

    vector< int >points_id;
    vector< vector< int > >all_points_id;

    // vector of nodal stresses - Matrices (tensors) dim x dim
    // TODO make this only for particles, now vector nodal_stress has length of all nodes, vertices and auxnodes (then node_id is needed for each matrix of nodal stresss - either pair <unsigned, Matrix> or two vectors - vector<unsigned> + vector<Matrix>) combine it at the beginning in init() - then it will apply also for adaptivity
    vector< Matrix >nodal_stress;
    bool export_nodal_stress = isStringInVect("nodal_stress", codes);

    vector< int >cell_types;
    vector< int >offsets;
    vector< Point >displ;

    vector< string >materials;
    unsigned matI;

    vector< bool >codes_positions(cell_data_size); // position indeces of data that cannot be exported from points or elements directly (they are not stored there)
    for ( unsigned i = 0; i < cell_data_size; i++ ) {
        codes_positions [ i ] = isAddonCellScalarData(codes [ i ]);
    }

    vector< vector< Vector > >cell_vect_data;
    vector< unsigned >vector_data_code_indeces;  // indeces of data to be exported as vectors

    vector< vector< double > >cell_data;
    cell_data.resize(cell_data_size);

    vector< vector< double > >point_data;
    point_data.resize(codes.size() - cell_data_size);

    size_t offset = 0;
    for ( auto const &el : * elems ) {
        for ( auto const &n : el->giveNodes() ) {
            auto res = std :: find(begin(* nodes), end(* nodes), n);
            points_id.push_back(std :: distance(begin(* nodes), res) );
            // points_id.push_back(n - *nodes->begin());
        }
        all_points_id.push_back(points_id);
        cell_types.push_back(points_id.size() * 2 - 1); // NOTE this works for line (type 3), triangle (type 5), be careful with quad (type 9), but closed polygon is type 7, needs to be enhanced for bricks etc...
        offset += points_id.size();
        offsets.push_back(offset);
        for ( unsigned i = 0; i < cell_data_size; i++ ) {
            if ( isAddonCellScalarData(codes [ i ]) || isAddonPointVectorialData(codes [ i ]) ) {
                // TODO this needs to be improved, it is duplicated in these two functions
                continue;
            } else if ( isAddonCellVectorialData(codes [ i ]) ) {
                vector_data_code_indeces.push_back(i);
                continue;
            } else if ( codes [ i ].compare("material") == 0 ) {
                matI = 0;
                while ( true ) {
                    if ( matI >= materials.size() ) {
                        materials.push_back(el->giveMaterial()->giveName() );
                        break;
                    } else if ( materials.size() > 0 && materials [ matI ].compare(el->giveMaterial()->giveName() ) == 0 ) {
                        break;
                    }
                    matI++;
                }
                cell_data [ i ].push_back(matI);
            } else {
                cell_data [ i ].push_back(el->giveIPValue(codes [ i ], 0) ); // so far for single IP point
            }
        }
        points_id.clear();
    }
    cell_vect_data.resize(vector_data_code_indeces.size() );

    if ( !std :: none_of(codes_positions.begin(), codes_positions.end(), [ ](bool i) {
        // TODO toto je nutný udělat jinak, teď, když jsou tady i jiný sady (vektorový data), není to prostě buď jedno nebo druhý
        return i == true;
    }) ) {
        exportAddonScalarCellData(elems, DoFs, codes_positions, codes, cell_data);
    }
    if ( vector_data_code_indeces.size() > 0 ) {
        exportAddonVectorialCellData(this->dim, elems, DoFs, codes, vector_data_code_indeces, cell_vect_data);
    }
    if ( export_nodal_stress ) {
        // reserve space only if nodal stresses should be exported
        nodal_stress.resize(nodes->giveSize(), Matrix(this->dim, this->dim) );
        // export nodal stresses:
        ExportAllElementsNodalStress(nodal_stress, DoFs, reactions, nodes, elems, this->dim);
    }

    giveFileName(step, buffer);
    ofstream outputfile( ( resultDir / buffer ).string() );

    if ( outputfile.is_open() ) {
        outputfile << std :: scientific;
        outputfile.precision(precision);
        outputfile << "<?xml version=\"1.0\"?>\n";
        outputfile << "<VTKFile type=\"UnstructuredGrid\" version=\"0.1\" byte_order=\"LittleEndian\">" << '\n';
        outputfile << "<UnstructuredGrid>" << '\n';
        outputfile << "<Piece NumberOfPoints=\"" << nodes->giveSize() << "\" NumberOfCells=\"" << elems->giveSize() << "\">" << '\n';


        outputfile << "<Points>" << '\n';
        outputfile << "<DataArray type=\"Float32\" NumberOfComponents=\"3\" format=\"ascii\">" << '\n';
        for ( auto const &n : * nodes ) {
            outputfile << n->givePoint().getX() << "\t" << n->givePoint().getY() << "\t" << n->givePoint().getZ() << '\n';

            displ.push_back(Point(n->giveDoFBasedValue("ux", DoFs),
                                  n->giveDoFBasedValue("uy", DoFs),
                                  dim == 3 ? n->giveDoFBasedValue("uz", DoFs) : 0
                                  )
                            );

            for ( unsigned i = cell_data_size; i < codes.size(); i++ ) {
                point_data [ i - cell_data_size ].push_back(n->giveDoFBasedValue(codes [ i ], DoFs) );
            }
        }
        // for (unsigned n; n < nodes->giveSize(); n++) {
        //   P = nodes->giveNode(n)->givePoint();
        // }
        outputfile << "</DataArray>" << "\n";
        outputfile << "</Points>" << "\n";
        // /*
        outputfile << "<Cells>" << '\n';
        outputfile << "<DataArray type=\"Int32\" Name=\"connectivity\" format=\"ascii\">" << '\n';
        for ( auto const &value : all_points_id ) {
            for ( auto const &id : value ) {
                outputfile << "\t" << id;
            }
            outputfile << '\n';
        }
        outputfile << "</DataArray>" << '\n';
        outputfile << "<DataArray type=\"Int32\" Name=\"offsets\" format=\"ascii\">" << '\n';
        for ( auto const &value : offsets ) {
            outputfile << value << '\n';
        }
        outputfile << "</DataArray>" << '\n';
        outputfile << "<DataArray type=\"UInt8\" Name=\"types\" format=\"ascii\">" << '\n';
        for ( auto const &value : cell_types ) {
            outputfile << value << '\n';
        }
        outputfile << "</DataArray>" << '\n';
        outputfile << "</Cells>" << '\n';
        // */
        outputfile << "<PointData Scalars=\"scalars\">" << "\n";
        outputfile << "<DataArray type=\"Float32\" Name=\"displacement\" NumberOfComponents=\"3\" format=\"ascii\">" << '\n';
        for ( auto const &p : displ ) {
            outputfile << p.getX() << '\t' << p.getY() << '\t' << p.getZ() << '\n';
        }
        outputfile << "</DataArray>" << '\n';
        //////////////////////////////////////////////////////////////////////////
        if ( export_nodal_stress ) {
            outputfile << "<DataArray type=\"Float32\" Name=\"nodal_stress\" NumberOfComponents=\"" << ( dim - 1 ) * 3 << "\" format=\"ascii\">" << '\n';
            vector< double >data;

            for ( auto const &s : nodal_stress ) {
                data = MatrixToStdVectForParaview(s, dim);

                for ( auto const &d : data ) {
                    outputfile << d << '\t';
                }
                outputfile << '\n';
                // data.clear();
            }
            outputfile << "</DataArray>" << '\n';
        }
        //////////////////////////////////////////////////////////////////////////
        for ( unsigned i = 0; i < point_data.size(); i++ ) {
            if ( codes [ i + cell_data_size ].compare("nodal_stress") == 0 ) {
                continue;
            }
            outputfile << "<DataArray type=\"Float32\" Name=\" " << codes [ i + cell_data_size ] << "\" format=\"ascii\">" << '\n';
            for ( auto const &p : point_data [ i ] ) {
                outputfile << p << '\n';
            }
            outputfile << "</DataArray>" << '\n';
        }
        //////////////////////////////////////////////////////////////////////////
        outputfile << "</PointData>" << '\n';
        // /*
        outputfile << "<CellData Scalars=\"scalars\">" << "\n";
        unsigned num_vectors = 0;
        for ( unsigned i = 0; i < cell_data.size(); i++ ) {
            if ( isInVect(i, vector_data_code_indeces) ) {
                outputfile << "<DataArray type=\"Float32\" Name=\"" << codes [ i ] <<
                    "\" NumberOfComponents=\"" << dim << "\"  format=\"ascii\">" << '\n';
                for ( auto const &value : cell_vect_data [ num_vectors ] ) {
                    for ( auto const &a : value ) {
                        outputfile << a << '\t';
                    }
                    outputfile << '\n';
                }
                outputfile << "</DataArray>" << '\n';
                num_vectors += 1;
            } else {
                outputfile << "<DataArray type=\"Float32\" Name=\"" << codes [ i ] << "\" format=\"ascii\">" << '\n';
                for ( auto const &value : cell_data [ i ] ) {
                    outputfile << value << '\n';
                }
                outputfile << "</DataArray>" << '\n';
            }
        }
        outputfile << "</CellData>" << '\n';
        // */
        outputfile << "</Piece>" << '\n';
        outputfile << "</UnstructuredGrid>" << '\n';
        outputfile << "</VTKFile>" << '\n';
        outputfile.close();
    }
}

//////////////////////////////////////////////////////////
// function tahat calculates displacement of any point of rigid body from its rotations and
Point calculateVertexDisplacement(const RigidBodyContact &rbc, const Node *v, const Node *a, const Vector &DoFs, const unsigned &dim) {
    Matrix A = rbc.giveAMatrix(a->givePoint(), v->givePoint() );

    unsigned DofsPerNode = ( dim - 1 ) * 3;
    Matrix U(DofsPerNode, 1);
    for ( unsigned i = 0; i < DofsPerNode; i++ ) {
        U [ i ] [ 0 ] = DoFs [ a->giveStartingDoF() + i ];
    }

    Matrix P = A * U;

    return Point(P [ 0 ] [ 0 ], P [ 1 ] [ 0 ], P.numRows() > 2 ? P [ 2 ] [ 0 ] : 0);
}
//////////////////////////////////////////////////////////

// RIGID POLYGONS TO VTU FILE
//////////////////////////////////////////////////////////
void VTKRB2DExporter :: exportData(unsigned step, const Vector &DoFs, const Vector &reactions, fs :: path resultDir) const {
    // Export of elements into vtu xml file format (vtu = vtk for unstructured grid)
    // NOTE this is messy construction of xml file, will be remade using some of xml libraries for cpp
    char buffer [ 100 ];
    // Point P;
    // Element *ee;

    vector< Point >all_vertices_twice; // because the vertices are needed twice (on each side of the contact)
    vector< Point >vertices_displ; // displaceemnt of the contact
    vector< int >points_id;
    vector< int >node_id;
    vector< vector< int > >all_points_id;
    vector< int >cell_types;
    vector< int >offsets;
    vector< Point >displ;
    vector< double >damage; // test version, this and more will be specified on the exporter input
    int offset = 0;
    for ( auto const &el : * elems ) {
        RigidBodyContact *rbc = static_cast< RigidBodyContact * >( el );
        for ( auto const &n : rbc->giveNodes() ) {
            // Particle *part = static_cast< Particle * >( n );
            auto nod_id_ptr = std :: find(begin(* nodes), end(* nodes), n);
            for ( auto const &v : rbc->giveVertices() ) {
                all_vertices_twice.push_back(v->givePoint() );
                points_id.push_back(nodes->giveSize() + all_vertices_twice.size() - 1);
                vertices_displ.push_back(calculateVertexDisplacement(* rbc, v, n, DoFs, this->dim) );
            }
            points_id.push_back(std :: distance(begin(* nodes), nod_id_ptr) );
            node_id.push_back(std :: distance(begin(* nodes), nod_id_ptr) );
            all_points_id.push_back(points_id);
            cell_types.push_back(points_id.size() * 2 - 1); // NOTE this works for line (type 3), triangle (type 5), be careful with quad (type 9), but closed polygon is type 7, needs to be enhanced for bricks etc...
            offset += points_id.size();
            offsets.push_back(offset);
            damage.push_back(el->giveValue("damage") );
            points_id.clear();
        }
    }
    // for (unsigned e; e < elems->giveSize(); e++){
    // }

    giveFileName(step, buffer);
    ofstream outputfile( ( resultDir / buffer ).string() );

    if ( outputfile.is_open() ) {
        outputfile << std :: scientific;
        outputfile.precision(precision);
        outputfile << "<?xml version=\"1.0\"?>\n";
        outputfile << "<VTKFile type=\"UnstructuredGrid\" version=\"0.1\" byte_order=\"LittleEndian\">" << '\n';
        outputfile << "<UnstructuredGrid>" << '\n';
        outputfile << "<Piece NumberOfPoints=\"" << nodes->giveSize() + all_vertices_twice.size() << "\" NumberOfCells=\"" << elems->giveSize() * 2 << "\">" << '\n';


        outputfile << "<Points>" << '\n';
        outputfile << "<DataArray type=\"Float32\" NumberOfComponents=\"3\" format=\"ascii\">" << '\n';
        for ( auto const &n : * nodes ) {
            displ.push_back(Point(n->giveDoFBasedValue("ux", DoFs),
                                  n->giveDoFBasedValue("uy", DoFs),
                                  dim == 3 ? n->giveDoFBasedValue("uz", DoFs) : 0
                                  ) ); // for 2D this DOF stays for rotation
            outputfile << n->givePoint().getX() << "\t" << n->givePoint().getY() << "\t" << n->givePoint().getZ() << '\n';
        }
        for ( auto const &p : all_vertices_twice ) {
            outputfile << p.getX() << "\t" << p.getY() << "\t" << p.getZ() << '\n';
        }
        // for (unsigned n; n < nodes->giveSize(); n++) {
        //   P = nodes->giveNode(n)->givePoint();
        // }
        outputfile << "</DataArray>" << "\n";
        outputfile << "</Points>" << "\n";
        // /*
        outputfile << "<Cells>" << '\n';
        outputfile << "<DataArray type=\"Int32\" Name=\"connectivity\" format=\"ascii\">" << '\n';
        for ( auto const &value : all_points_id ) {
            for ( auto const &id : value ) {
                outputfile << "\t" << id;
            }
            outputfile << '\n';
        }
        outputfile << "</DataArray>" << '\n';
        outputfile << "<DataArray type=\"Int32\" Name=\"offsets\" format=\"ascii\">" << '\n';
        for ( auto const &value : offsets ) {
            outputfile << value << '\n';
        }
        outputfile << "</DataArray>" << '\n';
        outputfile << "<DataArray type=\"UInt8\" Name=\"types\" format=\"ascii\">" << '\n';
        for ( auto const &value : cell_types ) {
            outputfile << value << '\n';
        }
        outputfile << "</DataArray>" << '\n';
        outputfile << "</Cells>" << '\n';
        // */
        outputfile << "<PointData Scalars=\"scalars\">" << "\n";
        outputfile << "<DataArray type=\"Float32\" Name=\"displacement\" NumberOfComponents=\"3\" format=\"ascii\">" << '\n';
        for ( auto const &p : displ ) {
            outputfile << p.getX() << '\t' << p.getY() << '\t' << p.getZ() << '\n';
        }
        for ( auto const &p : vertices_displ ) {
            outputfile << p.getX() << '\t' << p.getY() << '\t' << p.getZ() << '\n';
        }
        outputfile << "</DataArray>" << '\n';
        outputfile << "</PointData>" << '\n';
        // /*
        outputfile << "<CellData Scalars=\"scalars\">" << "\n";
        outputfile << "<DataArray type=\"Float32\" Name=\"damage\" format=\"ascii\">" << '\n';
        for ( auto const &value : damage ) {
            outputfile << value << '\n';
        }
        outputfile << "</DataArray>" << '\n';

        outputfile << "<DataArray type=\"Float32\" Name=\"node_id\" format=\"ascii\">" << '\n';
        for ( auto const &value : node_id ) {
            outputfile << value << '\n';
        }
        outputfile << "</DataArray>" << '\n';

        outputfile << "</CellData>" << '\n';
        // */
        outputfile << "</Piece>" << '\n';
        outputfile << "</UnstructuredGrid>" << '\n';
        outputfile << "</VTKFile>" << '\n';
        outputfile.close();
    }
}


// RIGID contacts TO VTU FILE
//////////////////////////////////////////////////////////
void VTKRCExporter :: exportData(unsigned step, const Vector &DoFs, const Vector &reactions, fs :: path resultDir) const {
    // Export of elements into vtu xml file format (vtu = vtk for unstructured grid)
    // NOTE this is messy construction of xml file, will be remade using some of xml libraries for cpp
    char buffer [ 100 ];
    // Point P;
    // Element *ee;

    vector< Point >all_vertices_twice; // because the vertices are needed twice (on each side of the contact)
    vector< Point >vertices_displ; // displaceemnt of the contact
    vector< int >points_id;
    vector< int >node_id;
    vector< vector< int > >all_points_id;
    vector< int >cell_types;
    vector< int >offsets;
    vector< Point >displ;

    // vector of nodal stresses - Matrices (tensors) dim x dim
    vector< Matrix >nodal_stress;
    vector< unsigned >vertex_correponding_node_id;
    bool export_nodal_stress = isStringInVect("nodal_stress", codes);

    vector< bool >codes_positions(cell_data_size); // position indeces of data that cannot be exported from points or elements directly (they are not stored there)
    for ( unsigned i = 0; i < cell_data_size; i++ ) {
        codes_positions [ i ] = isAddonCellScalarData(codes [ i ]);
    }

    vector< vector< Vector > >cell_vect_data;
    vector< unsigned >vector_data_code_indeces;  // indeces of data to

    vector< vector< double > >cell_data; // test version, this and more will be specified on the exporter input
    cell_data.resize(cell_data_size);

    // TODO be able to export also point data
    // vector< vector< double > > point_data;
    // point_data.resize(codes.size() - cell_data_size);
    unsigned num_rbcs = 0;

    int offset = 0;
    unsigned node_id_i;
    for ( auto const &el : * elems ) {
        if ( el->giveName().compare("RigidBodyContact") != 0 ) {
            continue;
        }
        RigidBodyContact *rbc = static_cast< RigidBodyContact * >( el );
        num_rbcs++;
        for ( auto const &n : rbc->giveNodes() ) {
            // Particle *part = static_cast< Particle * >( n );
            auto nod_id_ptr = std :: find(begin(* nodes), end(* nodes), n);
            node_id_i = std :: distance(begin(* nodes), nod_id_ptr);
            for ( auto const &v : rbc->giveVertices() ) {
                all_vertices_twice.push_back(v->givePoint() );
                points_id.push_back(nodes->giveSize() + all_vertices_twice.size() - 1);
                vertices_displ.push_back(calculateVertexDisplacement(* rbc, v, n, DoFs, this->dim) );
                if ( export_nodal_stress ) {
                    vertex_correponding_node_id.push_back(node_id_i);
                }
            }
            // points_id.push_back(std::distance(begin(*nodes), nod_id_ptr));
            node_id.push_back(node_id_i);
            all_points_id.push_back(points_id);
            cell_types.push_back( ( points_id.size() > 2 ) ? 7 : points_id.size() * 2 - 1 );
            // cell_types.push_back(points_id.size()*2 - 1);  // NOTE this works for line (type 3), triangle (type 5), be careful with quad (type 9), but closed polygon is type 7, needs to be enhanced for bricks etc...
            offset += points_id.size();
            offsets.push_back(offset);

            for ( unsigned i = 0; i < cell_data_size; i++ ) {
                if ( isAddonCellScalarData(codes [ i ]) || isAddonPointVectorialData(codes [ i ]) ) {
                    continue;
                } else if ( isAddonCellVectorialData(codes [ i ]) ) {
                    vector_data_code_indeces.push_back(i);
                    continue;
                } else {
                    cell_data [ i ].push_back(el->giveIPValue(codes [ i ], 0) ); // so far for single IP point
                }
            }
            points_id.clear();
        }
    }
    cell_vect_data.resize(vector_data_code_indeces.size() );

    if ( !std :: none_of(codes_positions.begin(), codes_positions.end(), [ ](bool i) {
        return i == true;
    }) ) {
        exportAddonScalarCellData(elems, DoFs, codes_positions, codes, cell_data, true, true);
    }
    if ( vector_data_code_indeces.size() > 0 ) {
        exportAddonVectorialCellData(this->dim, elems, DoFs, codes, vector_data_code_indeces, cell_vect_data, true, true);
    }
    if ( export_nodal_stress ) {
        // reserve space only if nodal stresses should be exported
        nodal_stress.resize(nodes->giveSize(), Matrix(this->dim, this->dim) );
        // export nodal stresses:
        ExportAllElementsNodalStress(nodal_stress, DoFs, reactions, nodes, elems, this->dim);
    }

    // unsigned iii = 0;
    // for (auto const &da : cell_data){
    //   std::cout << "cell_data[" << iii << "] = " << da.size() << '\n';
    // }
    // std::cout << "cell_types.size() = " << cell_types.size() << '\n';
    // std::cout << "offsets.size() = " << offsets.size() << '\n';
    // std::cout << "all_points_id.size() = " << all_points_id.size() << '\n';

    giveFileName(step, buffer);
    ofstream outputfile( ( resultDir / buffer ).string() );

    if ( outputfile.is_open() ) {
        outputfile << std :: scientific;
        outputfile.precision(precision);
        outputfile << "<?xml version=\"1.0\"?>\n";
        outputfile << "<VTKFile type=\"UnstructuredGrid\" version=\"0.1\" byte_order=\"LittleEndian\">" << '\n';
        outputfile << "<UnstructuredGrid>" << '\n';
        outputfile << "<Piece NumberOfPoints=\"" << nodes->giveSize() + all_vertices_twice.size() << "\" NumberOfCells=\"" << num_rbcs * 2 << "\">" << '\n';


        outputfile << "<Points>" << '\n';
        outputfile << "<DataArray type=\"Float32\" NumberOfComponents=\"3\" format=\"ascii\">" << '\n';
        for ( auto const &n : * nodes ) {
            displ.push_back(Point(n->giveDoFBasedValue("ux", DoFs),
                                  n->giveDoFBasedValue("uy", DoFs),
                                  0) ); // for 2D this DOF stays for rotation
            outputfile << n->givePoint().getX() << "\t" << n->givePoint().getY() << "\t" << n->givePoint().getZ() << '\n';
        }
        for ( auto const &p : all_vertices_twice ) {
            outputfile << p.getX() << "\t" << p.getY() << "\t" << p.getZ() << '\n';
        }
        // for (unsigned n; n < nodes->giveSize(); n++) {
        //   P = nodes->giveNode(n)->givePoint();
        // }
        outputfile << "</DataArray>" << "\n";
        outputfile << "</Points>" << "\n";
        // /*
        outputfile << "<Cells>" << '\n';
        outputfile << "<DataArray type=\"Int32\" Name=\"connectivity\" format=\"ascii\">" << '\n';
        for ( auto const &value : all_points_id ) {
            for ( auto const &id : value ) {
                outputfile << "\t" << id;
            }
            outputfile << '\n';
        }
        outputfile << "</DataArray>" << '\n';
        outputfile << "<DataArray type=\"Int32\" Name=\"offsets\" format=\"ascii\">" << '\n';
        for ( auto const &value : offsets ) {
            outputfile << value << '\n';
        }
        outputfile << "</DataArray>" << '\n';
        outputfile << "<DataArray type=\"UInt8\" Name=\"types\" format=\"ascii\">" << '\n';
        for ( auto const &value : cell_types ) {
            outputfile << value << '\n';
        }
        outputfile << "</DataArray>" << '\n';
        outputfile << "</Cells>" << '\n';
        // */
        outputfile << "<PointData Scalars=\"scalars\">" << "\n";
        outputfile << "<DataArray type=\"Float32\" Name=\"displacement\" NumberOfComponents=\"3\" format=\"ascii\">" << '\n';
        for ( auto const &p : displ ) {
            outputfile << p.getX() << '\t' << p.getY() << '\t' << p.getZ() << '\n';
        }
        for ( auto const &p : vertices_displ ) {
            outputfile << p.getX() << '\t' << p.getY() << '\t' << p.getZ() << '\n';
        }
        outputfile << "</DataArray>" << '\n';

        ///////////////////////////////////////////////////////////////////////
        if ( export_nodal_stress ) {
            outputfile << "<DataArray type=\"Float32\" Name=\"nodal_stress\" NumberOfComponents=\"" << ( dim - 1 ) * 3 << "\" format=\"ascii\">" << '\n';
            vector< double >data;

            for ( auto const &s : nodal_stress ) {
                data = MatrixToStdVectForParaview(s, dim);

                for ( auto const &d : data ) {
                    outputfile << d << '\t';
                }
                outputfile << '\n';
                // data.clear();
            }
            for ( auto const &is : vertex_correponding_node_id ) {
                data = MatrixToStdVectForParaview(nodal_stress [ is ], dim);

                for ( auto const &d : data ) {
                    outputfile << d << '\t';
                }
                outputfile << '\n';
                // data.clear();
            }


            outputfile << "</DataArray>" << '\n';
        }
        ///////////////////////////////////////////////////////////////////////

        outputfile << "</PointData>" << '\n';
        // /*
        unsigned num_vectors = 0;
        outputfile << "<CellData Scalars=\"scalars\">" << "\n";
        for ( unsigned i = 0; i < cell_data.size(); i++ ) {
            if ( isInVect(i, vector_data_code_indeces) ) {
                outputfile << "<DataArray type=\"Float32\" Name=\"" << codes [ i ] <<
                    "\" NumberOfComponents=\"" << dim << "\"  format=\"ascii\">" << '\n';
                for ( auto const &value : cell_vect_data [ num_vectors ] ) {
                    for ( auto const &a : value ) {
                        outputfile << a << '\t';
                    }
                    outputfile << '\n';
                }
                outputfile << "</DataArray>" << '\n';
                num_vectors += 1;
            } else {
                outputfile << "<DataArray type=\"Float32\" Name=\"" << codes [ i ] << "\" format=\"ascii\">" << '\n';
                for ( auto const &value : cell_data [ i ] ) {
                    outputfile << value << '\n';
                }
                outputfile << "</DataArray>" << '\n';
            }
        }

        outputfile << "<DataArray type=\"Float32\" Name=\"node_id\" format=\"ascii\">" << '\n';
        for ( auto const &value : node_id ) {
            outputfile << value << '\n';
        }
        outputfile << "</DataArray>" << '\n';

        outputfile << "</CellData>" << '\n';
        // */
        outputfile << "</Piece>" << '\n';
        outputfile << "</UnstructuredGrid>" << '\n';
        outputfile << "</VTKFile>" << '\n';
        outputfile.close();
    }
}
