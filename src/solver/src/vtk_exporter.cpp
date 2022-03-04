#include "vtk_exporter.h"
#include "element_discrete.h"
#include "misc.h"

#ifdef __VTK_MODULE
    #include <vtkCellArray.h>
    #include <vtkNew.h>
    #include <vtkPointData.h>
    #include <vtkProperty.h>
    #include <vtkPlane.h>
    #include <vtkBitArray.h>
    #include <vtkUnstructuredGrid.h>
    #include <vtkXMLUnstructuredGridWriter.h>
    #include <vtkPoints.h>
    #include <vtkIntArray.h>
    #include <vtkDoubleArray.h>
    #include <vtkCellData.h>
    #include <vtkXMLPUnstructuredGridWriter.h>
    #include <vtkSmartPointer.h>
    #include <vtkCutter.h>
    #include <vtkXMLPolyDataWriter.h>

    //#include <vtkActor.h>
    //#include <vtkDataSetMapper.h>
    //#include <vtkNamedColors.h>
    //#include <vtkRenderWindow.h>
    //#include <vtkRenderWindowInteractor.h>
    //#include <vtkRenderer.h>
    //#include <vtkTetra.h>
    //#include <vtkVertexGlyphFilter.h>
    //#include <vtkXMLUnstructuredGridReader.h>
    //#include <vtkPolygon.h>
    //#include <vtkCamera.h>
    //#include <vtkPolyData.h>
    //#include <vtkPolyDataMapper.h>
    //#include <vtkXMLPolyDataReader.h>
    //#include <vtkMassProperties.h>
#endif


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
    vector< string >cellData, pointData, extPointData;
    iss >> filename;
    while ( !iss.eof() ) {
        iss >> param;
        if ( param.compare("cellData") == 0 ) {
            iss >> num;
            cellData.resize(num);
            for ( unsigned i = 0; i < num; i++ ) {
                iss >> cellData [ i ];
            }
        } else if ( param.compare("pointData") == 0 || param.compare("nodeData") == 0 ) {
            iss >> num;
            pointData.resize(num);
            for ( unsigned i = 0; i < num; i++ ) {
                iss >> pointData [ i ];
            }
        } else if ( param.compare("extrapolatedNodeData") == 0 ) {
            iss >> num;
            extPointData.resize(num);
            for ( unsigned i = 0; i < num; i++ ) {
                iss >> extPointData [ i ];
            }
        } else if ( param.compare("binary") == 0 ) {
            binaryswitch = true;
        }
    }
    codes.resize(cellData.size() + pointData.size() + extPointData.size() );
    num = 0;
    for ( auto const &cel : cellData ) {
        codes [ num++ ] = cel;
    }
    cell_data_size = cellData.size();
    for ( auto const &po : pointData ) {
        codes [ num++ ] = po;
    }
    node_data_size = pointData.size();
    for ( auto const &ex : extPointData ) {
        codes [ num++ ] = ex;
    }

    DataExporter :: readFromLine(iss);
}


//////////////////////////////////////////////////////////
// ELEMENTS TO VTU FILE
//////////////////////////////////////////////////////////
void VTKElementExporter :: exportData(unsigned step, const Vector &DoFs, const Vector &reactions, fs :: path resultDir) const {
    // Export of elements into vtu xml file format (vtu = vtk for unstructured grid)
    char buffer [ 100 ];
    giveFileName(step, buffer);
    // Point P;
    // Element *ee;

    vector< int >points_id;
    vector< vector< int > >all_points_id;

    #ifdef __VTK_MODULE
        vtkSmartPointer<vtkUnstructuredGrid> unstructuredGrid = vtkSmartPointer<vtkUnstructuredGrid>::New();

        vtkSmartPointer<vtkPoints> points = vtkSmartPointer<vtkPoints>::New();

        Point *pp;       
        for ( unsigned n = 0; n < nodes->giveSize(); n++ ) {
            pp = nodes->giveNode(n)->givePointPointer();
            points->InsertNextPoint(pp->getX(), pp->getY(), pp->getZ());
        }
        unstructuredGrid->SetPoints(points);

        vector <Node*> elnodes;
        Element *el;  
        for ( unsigned e = 0; e < elems->giveSize(); e++ ) {
            el = elems->giveElement(e);
            elnodes = el->giveNodes();
            vtkSmartPointer<vtkIdList> elindices = vtkSmartPointer<vtkIdList>::New();
            for(unsigned p=0; p<elnodes.size(); p++) elindices->InsertNextId(elnodes[p]->giveID());
            unstructuredGrid->InsertNextCell(el->giveVTKCellType(), elindices );
        }

        unsigned i,j;
        size_t msize;
        vector< Vector > data;
        unsigned p;
        // ****************** cell data
        data.resize(elems->giveSize());
        for (p = 0; p < cell_data_size; p++){
            msize = 1;
            i = 0;
            for (vector< Element * > :: const_iterator ee = elems->begin(); ee != elems->end(); ++ee , i++){
                (*ee)->giveValues(codes[p].c_str(),data[i]);
                msize = max(msize,data[i].size());
            }
            vtkSmartPointer< vtkDoubleArray > cellDataArray = vtkSmartPointer< vtkDoubleArray > :: New();
            cellDataArray->SetName(codes[p].c_str());
            cellDataArray->SetNumberOfComponents(msize);
            cellDataArray->SetNumberOfValues(elems->giveSize()*msize);
            i = 0;
            for(vector< Vector > :: const_iterator d = data.begin(); d != data.end(); ++d , i++){
                for(j=0; j<min(msize,d->size()); j++){
                    cellDataArray->SetValue(msize*i+j, (*d)[j]);
                }
                for(; j<msize; j++){
                    cellDataArray->SetValue(msize*i+j, 0);
                }
            }
            unstructuredGrid->GetCellData()->AddArray(cellDataArray);
        }


        // ****************** node data
        data.resize(nodes->giveSize());
        for (; p < node_data_size + cell_data_size; p++){
            msize = 1;
            i = 0;
            for (vector< Node * > :: const_iterator nn = nodes->begin(); nn != nodes->end(); ++nn , i++){
                (*nn)->giveDoFBasedValues(codes[p].c_str(),DoFs,data[i]);
                msize = max(msize,data[i].size());
            }
            vtkSmartPointer< vtkDoubleArray > pointDataArray = vtkSmartPointer< vtkDoubleArray > :: New();
            pointDataArray->SetName(codes[p].c_str());
            pointDataArray->SetNumberOfComponents(msize);
            pointDataArray->SetNumberOfValues(nodes->giveSize()*msize);
            i = 0;
            for(vector< Vector > :: const_iterator d = data.begin(); d != data.end(); ++d , i++){
                for(j=0; j<min(msize,d->size()); j++){
                    pointDataArray->SetValue(msize*i+j, (*d)[j]);
                }
                for(; j<msize; j++){
                    pointDataArray->SetValue(msize*i+j, 0);
                }
            }
            unstructuredGrid->GetPointData()->AddArray(pointDataArray);
        }

        // ****************** extrapolated node data        
        for (; p < codes.size() ; p++){
            msize = 1;

            elems->extrapolateValuesFromIntegrationPointsToNodes(codes[p], data);
            for (auto &v: data) msize = max(msize,v.size());

            vtkSmartPointer< vtkDoubleArray > pointDataArray = vtkSmartPointer< vtkDoubleArray > :: New();
            pointDataArray->SetName(codes[p].c_str());
            pointDataArray->SetNumberOfComponents(msize);
            pointDataArray->SetNumberOfValues(nodes->giveSize()*msize);
            i = 0;
            for(vector< Vector > :: const_iterator d = data.begin(); d != data.end(); ++d , i++){
                for(j=0; j<min(msize,d->size()); j++){
                    pointDataArray->SetValue(msize*i+j, (*d)[j]);
                }
                for(; j<msize; j++){
                    pointDataArray->SetValue(msize*i+j, 0);
                }
            }
            unstructuredGrid->GetPointData()->AddArray(pointDataArray);
        }
        
        //vtkNew<vtkXMLUnstructuredGridWriter> writer;
        vtkSmartPointer<vtkXMLUnstructuredGridWriter> writer = vtkSmartPointer<vtkXMLUnstructuredGridWriter>::New();
        writer->SetFileName(  ( resultDir / buffer ).string().c_str() );
        writer->SetInputData(unstructuredGrid);
        if (binaryswitch) writer->SetDataModeToBinary();
        else writer->SetDataModeToAscii();
        //writer->SetCompressorType();
        writer->Write();
    #endif
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

    (void) step; (void) DoFs; (void) reactions; (void) resultDir;
    /*
    // Export of elements into vtu xml file format (vtu = vtk for unstructured grid)
    // NOTE this is messy construction of xml file, will be remade using some of xml libraries for cpp
    char buffer [ 100 ];
    // Point P;
    // Element *ee;
    ( void ) reactions;

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
                                  ) );  // for 2D this DOF stays for rotation
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
        outputfile << "</Piece>" << '\n';
        outputfile << "</UnstructuredGrid>" << '\n';
        outputfile << "</VTKFile>" << '\n';
        outputfile.close();
    }
    */
}


// RIGID contacts TO VTU FILE
//////////////////////////////////////////////////////////
void VTKRCExporter :: exportData(unsigned step, const Vector &DoFs, const Vector &reactions, fs :: path resultDir) const {
    (void) step; (void) DoFs; (void) reactions; (void) resultDir;
    /*
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
        if ( el->giveName().compare("LTCBEAM") != 0 && el->giveName().compare("LTCBEAMCoupled") != 0 ) {
            continue;
        }
        RigidBodyContact *rbc = static_cast< RigidBodyContact * >( el );
        num_rbcs++;
        for ( auto const &n : rbc->giveNodes() ) {
            // Particle *part = static_cast< Particle * >( n );
            auto nod_id_ptr = std :: find(begin(* nodes), end(* nodes), n);
            node_id_i = std :: distance(begin(* nodes), nod_id_ptr);// todo warning C4244: '=': conversion from '__int64' to 'unsigned int', possible loss of data
            for ( auto const &v : rbc->giveVertices() ) {
                all_vertices_twice.push_back(v->givePoint() );
                points_id.push_back(nodes->giveSize() + all_vertices_twice.size() - 1); //tofo warning C4267: 'argument': conversion from 'size_t' to '_Ty', possible loss of data
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
                    cell_data [ i ].push_back(el->giveIPValue(codes [ i ], 0) );   // so far for single IP point
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
                                  0) );  // for 2D this DOF stays for rotation
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
        outputfile << "</Piece>" << '\n';
        outputfile << "</UnstructuredGrid>" << '\n';
        outputfile << "</VTKFile>" << '\n';
        outputfile.close();
    }
    */
}
