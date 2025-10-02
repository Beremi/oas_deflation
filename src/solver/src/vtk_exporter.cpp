#include "vtk_exporter.h"
#include "element_discrete.h"
#include "element_beam.h"
#include "element_ldpm.h"
#include "solver.h"
#include "misc.h"

using namespace std;

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

std::string VTKExporter :: giveFileName(unsigned step, int iteration) const {
    std::ostringstream oss;
    if (iteration < 0) {
        oss << filename << "_" << std::setfill('0') << std::setw(5) << step << ".vtu";
    } else {
        oss << filename << "_" << std::setfill('0') << std::setw(5) << step 
            << "_iter_" << std::setfill('0') << std::setw(5) << iteration << ".vtu";
    }
    return oss.str();
}

//////////////////////////////////////////////////////////
void VTKExporter :: readFromLine(istringstream &iss) {
    string param;
    unsigned num;
    vector< string >cellData, pointData, extPointData;
    iss >> filename;
    while (  iss >> param ) {
        if ( param.compare("cellData") == 0 ||  param.compare("elementData") == 0 ) {
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
        } else if ( param.compare("ascii") == 0 ) {
            binaryswitch = false;
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
void VTKElementExporter :: exportData(unsigned step, int iteration, fs :: path resultDir) const {
    // Export of elements into vtu xml file format (vtu = vtk for unstructured grid)
    std::string fname = giveFileName(step, iteration);
    // Point P;
    // Element *ee;

#ifdef __VTK_MODULE

    Vector DoFs = solver->giveTrialDoFValues();

    vtkSmartPointer< vtkUnstructuredGrid >unstructuredGrid = vtkSmartPointer< vtkUnstructuredGrid > :: New();

    vtkSmartPointer< vtkPoints >points = vtkSmartPointer< vtkPoints > :: New();

    Point *pp;
    for ( unsigned n = 0; n < nodes->giveSize(); n++ ) {
        pp = nodes->giveNode(n)->givePointPointer();
        points->InsertNextPoint(pp->x(), pp->y(), pp->z() );
    }
    unstructuredGrid->SetPoints(points);

    vector< Node * >elnodes;
    Element *el;
    for ( unsigned e = 0; e < elems->giveSize(); e++ ) {
        el = elems->giveElement(e);
        elnodes = el->giveNodes();
        vtkSmartPointer< vtkIdList >elindices = vtkSmartPointer< vtkIdList > :: New();
        for ( unsigned p = 0; p < elnodes.size(); p++ ) {
            elindices->InsertNextId(elnodes [ p ]->giveID() );
        }
        unstructuredGrid->InsertNextCell(el->giveVTKCellType(), elindices);
    }

    unsigned i, j;
    size_t msize;
    vector< Vector >data;
    unsigned p;
    // ****************** cell data
    data.resize(elems->giveSize() );
    for ( p = 0; p < cell_data_size; p++ ) {
        msize = 1;
        i = 0;
        for ( vector< Element * > :: const_iterator ee = elems->begin(); ee != elems->end(); ++ee, i++ ) {
            ( * ee )->giveValues(codes [ p ].c_str(), data [ i ]);
            msize = max< size_t >(msize, data [ i ].size() );
        }
        vtkSmartPointer< vtkDoubleArray >cellDataArray = vtkSmartPointer< vtkDoubleArray > :: New();
        cellDataArray->SetName(codes [ p ].c_str() );
        cellDataArray->SetNumberOfComponents(msize);
        cellDataArray->SetNumberOfValues(elems->giveSize() * msize);
        i = 0;
        for ( vector< Vector > :: const_iterator d = data.begin(); d != data.end(); ++d, i++ ) {
            for ( j = 0; j < min< size_t >(msize, d->size() ); j++ ) {
                cellDataArray->SetValue(msize * i + j, ( * d ) [ j ]);
            }
            for ( ; j < msize; j++ ) {
                cellDataArray->SetValue(msize * i + j, 0);
            }
        }
        unstructuredGrid->GetCellData()->AddArray(cellDataArray);
    }


    // ****************** node data
    data.resize(nodes->giveSize() );
    for ( ; p < node_data_size + cell_data_size; p++ ) {
        msize = 1;
        i = 0;
        if ( codes [ p ].compare("residuals") == 0 ) {
            msize = 0;
            unsigned stDoF = 0;
            unsigned numDoF = 0;
            Vector res = solver->giveResiduals();
            for ( vector< Node * > :: const_iterator nn = nodes->begin(); nn != nodes->end(); ++nn, i++ ) {
                stDoF = ( * nn )->giveStartingDoF();
                numDoF = ( * nn )->giveNumberOfDoFs();
                data [ i ] = Vector(numDoF);
                msize = max< size_t >(msize, numDoF);
                for ( unsigned t = 0; t < numDoF; t++ ) {
                    data [ i ] [ t ] = res [ stDoF + t ];
                }
            }
        } else {
            for ( vector< Node * > :: const_iterator nn = nodes->begin(); nn != nodes->end(); ++nn, i++ ) {
                ( * nn )->giveValues(codes [ p ].c_str(), solver, data [ i ]);
                msize = max< size_t >(msize, data [ i ].size() );
            }
        }
        vtkSmartPointer< vtkDoubleArray >pointDataArray = vtkSmartPointer< vtkDoubleArray > :: New();
        pointDataArray->SetName(codes [ p ].c_str() );
        pointDataArray->SetNumberOfComponents(msize);
        pointDataArray->SetNumberOfValues(nodes->giveSize() * msize);
        i = 0;
        for ( vector< Vector > :: const_iterator d = data.begin(); d != data.end(); ++d, i++ ) {
            for ( j = 0; j < min< size_t >(msize, d->size() ); j++ ) {
                pointDataArray->SetValue(msize * i + j, ( * d ) [ j ]);
            }
            for ( ; j < msize; j++ ) {
                pointDataArray->SetValue(msize * i + j, 0);
            }
        }
        unstructuredGrid->GetPointData()->AddArray(pointDataArray);
    }

    // ****************** extrapolated node data
    for ( ; p < codes.size(); p++ ) {
        msize = 1;

        elems->extrapolateValuesFromIntegrationPointsToNodes(codes [ p ], data);
        for ( auto &v: data ) {
            msize = max< size_t >(msize, v.size() );
        }

        vtkSmartPointer< vtkDoubleArray >pointDataArray = vtkSmartPointer< vtkDoubleArray > :: New();
        pointDataArray->SetName(codes [ p ].c_str() );
        pointDataArray->SetNumberOfComponents(msize);
        pointDataArray->SetNumberOfValues(nodes->giveSize() * msize);
        i = 0;
        for ( vector< Vector > :: const_iterator d = data.begin(); d != data.end(); ++d, i++ ) {
            for ( j = 0; j < min< size_t >(msize, d->size() ); j++ ) {
                pointDataArray->SetValue(msize * i + j, ( * d ) [ j ]);
            }
            for ( ; j < msize; j++ ) {
                pointDataArray->SetValue(msize * i + j, 0);
            }
        }
        unstructuredGrid->GetPointData()->AddArray(pointDataArray);
    }

    //vtkNew<vtkXMLUnstructuredGridWriter> writer;
    vtkSmartPointer< vtkXMLUnstructuredGridWriter >writer = vtkSmartPointer< vtkXMLUnstructuredGridWriter > :: New();
    
    // Use safe path handling for Windows/MXE compatibility
    fs::path outputPath = resultDir / fname;
    writer->SetFileName( outputPath.string().c_str() );
    
    writer->SetInputData(unstructuredGrid);
    if ( binaryswitch ) {
        writer->SetDataModeToBinary();
    } else {
        writer->SetDataModeToAscii();
    }
    //writer->SetCompressorType();
    writer->Write();
#else
    cout << "VTK library not install, export of Elements skipped" << endl;
#endif
}

//////////////////////////////////////////////////////////

// RIGID POLYGONS TO VTU FILE
//////////////////////////////////////////////////////////
void VTKRB2DExporter :: exportData(unsigned step, int iteration, fs :: path resultDir) const {
    ( void ) step;
    ( void ) resultDir;
    ( void ) iteration;
    /*
     * // Export of elements into vtu xml file format (vtu = vtk for unstructured grid)
     * // NOTE this is messy construction of xml file, will be remade using some of xml libraries for cpp
     * char buffer [ 100 ];
     * // Point P;
     * // Element *ee;
     * ( void ) reactions;
     *
     * vector< Point >all_vertices_twice; // because the vertices are needed twice (on each side of the contact)
     * vector< Point >vertices_displ; // displaceemnt of the contact
     * vector< int >points_id;
     * vector< int >node_id;
     * vector< vector< int > >all_points_id;
     * vector< int >cell_types;
     * vector< int >offsets;
     * vector< Point >displ;
     * vector< double >damage; // test version, this and more will be specified on the exporter input
     * int offset = 0;
     * for ( auto const &el : * elems ) {
     *  RigidBodyContact *rbc = static_cast< RigidBodyContact * >( el );
     *  for ( auto const &n : rbc->giveNodes() ) {
     *      // Particle *part = static_cast< Particle * >( n );
     *      auto nod_id_ptr = std :: find(begin(* nodes), end(* nodes), n);
     *      for ( auto const &v : rbc->giveVertices() ) {
     *          all_vertices_twice.push_back(v->givePoint() );
     *          points_id.push_back(nodes->giveSize() + all_vertices_twice.size() - 1);
     *          vertices_displ.push_back(calculateVertexDisplacement(* rbc, v, n, DoFs, this->dim) );
     *      }
     *      points_id.push_back(std :: distance(begin(* nodes), nod_id_ptr) );
     *      node_id.push_back(std :: distance(begin(* nodes), nod_id_ptr) );
     *      all_points_id.push_back(points_id);
     *      cell_types.push_back(points_id.size() * 2 - 1); // NOTE this works for line (type 3), triangle (type 5), be careful with quad (type 9), but closed polygon is type 7, needs to be enhanced for bricks etc...
     *      offset += points_id.size();
     *      offsets.push_back(offset);
     *      damage.push_back(el->giveValue("damage") );
     *      points_id.clear();
     *  }
     * }
     * // for (unsigned e; e < elems->giveSize(); e++){
     * // }
     *
     * giveFileName(step, buffer);
     * ofstream outputfile( ( resultDir / buffer ).string() );
     *
     * if ( outputfile.is_open() ) {
     *  outputfile << std :: scientific;
     *  outputfile.precision(precision);
     *  outputfile << "<?xml version=\"1.0\"?>\n";
     *  outputfile << "<VTKFile type=\"UnstructuredGrid\" version=\"0.1\" byte_order=\"LittleEndian\">" << '\n';
     *  outputfile << "<UnstructuredGrid>" << '\n';
     *  outputfile << "<Piece NumberOfPoints=\"" << nodes->giveSize() + all_vertices_twice.size() << "\" NumberOfCells=\"" << elems->giveSize() * 2 << "\">" << '\n';
     *
     *
     *  outputfile << "<Points>" << '\n';
     *  outputfile << "<DataArray type=\"Float32\" NumberOfComponents=\"3\" format=\"ascii\">" << '\n';
     *  for ( auto const &n : * nodes ) {
     *      displ.push_back(Point(n->giveValue("ux", DoFs),
     *                            n->giveValue("uy", DoFs),
     *                            dim == 3 ? n->giveValue("uz", DoFs) : 0
     *                            ) );  // for 2D this DOF stays for rotation
     *      outputfile << n->givePoint().getX() << "\t" << n->givePoint().getY() << "\t" << n->givePoint().getZ() << '\n';
     *  }
     *  for ( auto const &p : all_vertices_twice ) {
     *      outputfile << p.getX() << "\t" << p.getY() << "\t" << p.getZ() << '\n';
     *  }
     *  // for (unsigned n; n < nodes->giveSize(); n++) {
     *  //   P = nodes->giveNode(n)->givePoint();
     *  // }
     *  outputfile << "</DataArray>" << "\n";
     *  outputfile << "</Points>" << "\n";
     *  outputfile << "<Cells>" << '\n';
     *  outputfile << "<DataArray type=\"Int32\" Name=\"connectivity\" format=\"ascii\">" << '\n';
     *  for ( auto const &value : all_points_id ) {
     *      for ( auto const &id : value ) {
     *          outputfile << "\t" << id;
     *      }
     *      outputfile << '\n';
     *  }
     *  outputfile << "</DataArray>" << '\n';
     *  outputfile << "<DataArray type=\"Int32\" Name=\"offsets\" format=\"ascii\">" << '\n';
     *  for ( auto const &value : offsets ) {
     *      outputfile << value << '\n';
     *  }
     *  outputfile << "</DataArray>" << '\n';
     *  outputfile << "<DataArray type=\"UInt8\" Name=\"types\" format=\"ascii\">" << '\n';
     *  for ( auto const &value : cell_types ) {
     *      outputfile << value << '\n';
     *  }
     *  outputfile << "</DataArray>" << '\n';
     *  outputfile << "</Cells>" << '\n';
     *  outputfile << "<PointData Scalars=\"scalars\">" << "\n";
     *  outputfile << "<DataArray type=\"Float32\" Name=\"displacement\" NumberOfComponents=\"3\" format=\"ascii\">" << '\n';
     *  for ( auto const &p : displ ) {
     *      outputfile << p.getX() << '\t' << p.getY() << '\t' << p.getZ() << '\n';
     *  }
     *  for ( auto const &p : vertices_displ ) {
     *      outputfile << p.getX() << '\t' << p.getY() << '\t' << p.getZ() << '\n';
     *  }
     *  outputfile << "</DataArray>" << '\n';
     *  outputfile << "</PointData>" << '\n';
     *  outputfile << "<CellData Scalars=\"scalars\">" << "\n";
     *  outputfile << "<DataArray type=\"Float32\" Name=\"damage\" format=\"ascii\">" << '\n';
     *  for ( auto const &value : damage ) {
     *      outputfile << value << '\n';
     *  }
     *  outputfile << "</DataArray>" << '\n';
     *
     *  outputfile << "<DataArray type=\"Float32\" Name=\"node_id\" format=\"ascii\">" << '\n';
     *  for ( auto const &value : node_id ) {
     *      outputfile << value << '\n';
     *  }
     *  outputfile << "</DataArray>" << '\n';
     *
     *  outputfile << "</CellData>" << '\n';
     *  outputfile << "</Piece>" << '\n';
     *  outputfile << "</UnstructuredGrid>" << '\n';
     *  outputfile << "</VTKFile>" << '\n';
     *  outputfile.close();
     * }
     */
}


// RIGID contacts TO VTU FILE
//////////////////////////////////////////////////////////
void VTKRCExporter :: exportData(unsigned step, int iteration, fs :: path resultDir) const {
    std::string fname = giveFileName(step, iteration);
    Vector DoFs = solver->giveTrialDoFValues();

#ifdef __VTK_MODULE
    vtkSmartPointer< vtkUnstructuredGrid >unstructuredGrid = vtkSmartPointer< vtkUnstructuredGrid > :: New();

    vtkSmartPointer< vtkPoints >points = vtkSmartPointer< vtkPoints > :: New();

    vector< RigidBodyContact * >exportedElemsRBC;
    vector< LDPMTetra * >exportedElemsTET;
    vector< Node * >vertices, elemnodes;
    Point *pp;
    unsigned pointID = 0;
    RigidBodyContact *rbc;
    LDPMTetra *tetra;
    unsigned celtype = 3;     //line
    vector< unsigned >fcodes;
    vector< unsigned >ncodes;


    unsigned nfaces = 0;
    if ( dim == 3 ) {
        celtype = 7;            //polygon
    }
    for ( vector< Element * > :: const_iterator ee = elems->begin(); ee != elems->end(); ++ee ) {
        rbc = dynamic_cast< RigidBodyContact * >( * ee );
        if ( rbc ) {
            exportedElemsRBC.push_back(rbc);
            vertices = rbc->giveVertices();
            vtkSmartPointer< vtkIdList >elindicesA = vtkSmartPointer< vtkIdList > :: New();
            vtkSmartPointer< vtkIdList >elindicesB = vtkSmartPointer< vtkIdList > :: New();
            for ( auto &p:vertices ) {
                pp = p->givePointPointer();
                elindicesA->InsertNextId(pointID);
                elindicesB->InsertNextId(pointID + 1);
                points->InsertNextPoint(pp->x(), pp->y(), pp->z() );          //every node twice
                points->InsertNextPoint(pp->x(), pp->y(), pp->z() );          //every node twice
                pointID++;
                pointID++;
            }
            unstructuredGrid->InsertNextCell(celtype, elindicesA);
            unstructuredGrid->InsertNextCell(celtype, elindicesB);
            nfaces++;
        }
    }
    for ( vector< Element * > :: const_iterator ff = elems->begin(); ff != elems->end(); ++ff ) {
        tetra = dynamic_cast< LDPMTetra * >( * ff );
        if ( tetra ) {
            exportedElemsTET.push_back(tetra);
            vertices = tetra->giveVertices();
            for ( unsigned i = 0; i < tetra->giveNumOfFacets(); i++ ) {
                fcodes = tetra->giveFacetVertCodes(i);
                vtkSmartPointer< vtkIdList >elindicesA = vtkSmartPointer< vtkIdList > :: New();
                vtkSmartPointer< vtkIdList >elindicesB = vtkSmartPointer< vtkIdList > :: New();
                for ( auto &p:fcodes ) {
                    pp = vertices [ p ]->givePointPointer();
                    elindicesA->InsertNextId(pointID);
                    elindicesB->InsertNextId(pointID + 1);
                    points->InsertNextPoint(pp->x(), pp->y(), pp->z() );              //every node twice
                    points->InsertNextPoint(pp->x(), pp->y(), pp->z() );              //every node twice
                    pointID++;
                    pointID++;
                }
                unstructuredGrid->InsertNextCell(celtype, elindicesA);
                unstructuredGrid->InsertNextCell(celtype, elindicesB);
                nfaces++;
            }
        }
    }
    unstructuredGrid->SetPoints(points);
    unsigned numOfPoints = pointID;

    unsigned i, j;
    size_t msize;
    vector< Vector >data;
    unsigned p;
    // ****************** DISPLACEMENTS
    Vector dataA, dataB;
    msize = 3;
    vtkSmartPointer< vtkDoubleArray >pointDataArray = vtkSmartPointer< vtkDoubleArray > :: New();
    pointDataArray->SetName("displacements");
    pointDataArray->SetNumberOfComponents(msize);
    pointDataArray->SetNumberOfValues(numOfPoints * msize);
    pointID = 0;
    Point displ;
    Particle *part;
    for ( vector< RigidBodyContact * > :: const_iterator ee = exportedElemsRBC.begin(); ee != exportedElemsRBC.end(); ++ee ) {
        vertices = ( * ee )->giveVertices();
        for ( auto &q:vertices ) {
            pp = q->givePointPointer();
            for ( unsigned k = 0; k < 2; k++ ) {
                part = static_cast< Particle * >( ( * ee )->giveNode(k) );
                displ = part->calculateRigidBodyMotionPoint(pp, DoFs);
                for ( p = 0; p < msize; p++ ) {
                    pointDataArray->SetValue(msize * pointID + p,   displ(p) );
                }
                pointID++;
            }
        }
    }


    for ( const auto &tet: exportedElemsTET ) {
        vertices = tet->giveVertices();
        elemnodes = tet->giveNodes();
        for ( unsigned h = 0; h < tet->giveNumOfFacets(); h++ ) {
            fcodes = tet->giveFacetVertCodes(h);
            ncodes = tet->giveFacetNodeCodes(h);
            for ( auto &q:fcodes ) {
                pp = vertices [ q ]->givePointPointer();
                for ( unsigned k = 0; k < 2; k++ ) {
                    part = static_cast< Particle * >( elemnodes [ ncodes [ k ] ] );
                    displ = part->calculateRigidBodyMotionPoint(pp, DoFs);
                    for ( p = 0; p < msize; p++ ) {
                        pointDataArray->SetValue(msize * pointID + p,   displ(p) );
                    }
                    pointID++;
                }
            }
        }
    }

    unstructuredGrid->GetPointData()->AddArray(pointDataArray);

    // ****************** cell data
    data.resize(nfaces);
    for ( p = 0; p < cell_data_size; p++ ) {
        msize = 1;
        i = 0;
        for ( vector< RigidBodyContact * > :: const_iterator ee = exportedElemsRBC.begin(); ee != exportedElemsRBC.end(); ++ee ) {
            ( * ee )->giveValues(codes [ p ].c_str(), data [ i ]);
            msize = max< size_t >(msize, data [ i ].size() );
            i++;
        }
        for ( const auto &tet: exportedElemsTET ) {
            for ( unsigned k = 0; k < tet->giveNumOfFacets(); k++ ) {
                tet->giveIPValues(codes [ p ].c_str(), k, data [ i ]);
                msize = max< size_t >(msize, data [ i ].size() );
                i++;
            }
        }
        vtkSmartPointer< vtkDoubleArray >cellDataArray = vtkSmartPointer< vtkDoubleArray > :: New();
        cellDataArray->SetName(codes [ p ].c_str() );
        cellDataArray->SetNumberOfComponents(msize);
        cellDataArray->SetNumberOfValues(2 * nfaces * msize);
        i = 0;
        for ( vector< Vector > :: const_iterator d = data.begin(); d != data.end(); ++d, i++ ) {
            for ( j = 0; j < min< size_t >(msize, d->size() ); j++ ) {
                cellDataArray->SetValue(msize * ( i * 2 ) + j, ( * d ) [ j ]);
                cellDataArray->SetValue(msize * ( i * 2 + 1 ) + j, ( * d ) [ j ]);
            }
            for ( ; j < msize; j++ ) {
                cellDataArray->SetValue(msize * ( i * 2 ) + j, 0);
                cellDataArray->SetValue(msize * ( i * 2 + 1 ) + j, 0);
            }
        }
        unstructuredGrid->GetCellData()->AddArray(cellDataArray);
    }


    // ****************** node data
    data.resize(nodes->giveSize() );
    for ( ; p < node_data_size + cell_data_size; p++ ) {
        if ( codes [ 0 ].compare("displacements") == 0 ) {
            continue;
        }
        msize = 1;
        i = 0;
        for ( vector< Node * > :: const_iterator nn = nodes->begin(); nn != nodes->end(); ++nn, i++ ) {
            if ( static_cast< Particle * >( * nn ) ) {
                ( * nn )->giveValues(codes [ p ].c_str(), solver, data [ i ]);
                msize = max< size_t >(msize, data [ i ].size() );
            } else {
                data [ i ].resize(0);
            }
        }
        vtkSmartPointer< vtkDoubleArray >cellDataArray = vtkSmartPointer< vtkDoubleArray > :: New();
        cellDataArray->SetName(codes [ p ].c_str() );
        cellDataArray->SetNumberOfComponents(msize);
        cellDataArray->SetNumberOfValues(2 * elems->giveSize() * msize);
        i = 0;
        for ( vector< RigidBodyContact * > :: const_iterator ee = exportedElemsRBC.begin(); ee != exportedElemsRBC.end(); ++ee, i++ ) {
            dataA = data [ ( * ee )->giveNode(0)->giveID() ];
            dataB = data [ ( * ee )->giveNode(1)->giveID() ];
            for ( j = 0; j < min< size_t >(msize, dataA.size() ); j++ ) {
                cellDataArray->SetValue(msize * ( i * 2 ) + j, dataA [ j ]);
            }
            for ( ; j < msize; j++ ) {
                cellDataArray->SetValue(msize * ( i * 2 ) + j, 0);
            }
            for ( j = 0; j < min< size_t >(msize, dataB.size() ); j++ ) {
                cellDataArray->SetValue(msize * ( i * 2 + 1 ) + j, dataB [ j ]);
            }
            for ( ; j < msize; j++ ) {
                cellDataArray->SetValue(msize * ( i * 2 + 1 ) + j, 0);
            }
        }
        unstructuredGrid->GetCellData()->AddArray(cellDataArray);
    }

    // ****************** extrapolated node data
    for ( ; p < codes.size(); p++ ) {
        cout << "Extrapolated node data not implemented in RC exporter" << endl;
    }


    //vtkNew<vtkXMLUnstructuredGridWriter> writer;
    vtkSmartPointer< vtkXMLUnstructuredGridWriter >writer = vtkSmartPointer< vtkXMLUnstructuredGridWriter > :: New();
    
    // Use safe path handling for Windows/MXE compatibility  
    fs::path outputPath = resultDir / fname;
    writer->SetFileName( outputPath.string().c_str() );
    
    writer->SetInputData(unstructuredGrid);
    if ( binaryswitch ) {
        writer->SetDataModeToBinary();
    } else {
        writer->SetDataModeToAscii();
    }
    //writer->SetCompressorType();
    writer->Write();
#else
    cout << "VTK library not install, export of Rigid Contacts skipped" << endl;
#endif
}


//////////////////////////////////////////////////////////
// REBARS TO VTU FILE
//////////////////////////////////////////////////////////
void VTKRebarExporter :: exportData(unsigned step, int iteration, fs :: path resultDir) const {
    // Export of elements into vtu xml file format (vtu = vtk for unstructured grid)
    std::string fname = giveFileName(step, iteration);
    // Point P;
    // Element *ee;

#ifdef __VTK_MODULE

    //SELECT ELEMENTS
    vector< Element * >selected_elements;
    Element *el;
    for ( unsigned e = 0; e < elems->giveSize(); e++ ) {
        el = elems->giveElement(e);
        TimoshenkoBeam3D *tb = dynamic_cast< TimoshenkoBeam3D * >( el );
        if ( tb ) {
            selected_elements.push_back(el);
        }
    }

    Vector DoFs = solver->giveTrialDoFValues();

    vtkSmartPointer< vtkUnstructuredGrid >unstructuredGrid = vtkSmartPointer< vtkUnstructuredGrid > :: New();

    vtkSmartPointer< vtkPoints >points = vtkSmartPointer< vtkPoints > :: New();

    Point *pp;
    for ( unsigned n = 0; n < nodes->giveSize(); n++ ) {
        pp = nodes->giveNode(n)->givePointPointer();
        points->InsertNextPoint(pp->x(), pp->y(), pp->z() );
    }
    unstructuredGrid->SetPoints(points);

    vector< Node * >elnodes;
    for ( unsigned e = 0; e < selected_elements.size(); e++ ) {
        el = selected_elements [ e ];
        elnodes = el->giveNodes();
        vtkSmartPointer< vtkIdList >elindices = vtkSmartPointer< vtkIdList > :: New();
        for ( unsigned p = 0; p < elnodes.size(); p++ ) {
            elindices->InsertNextId(elnodes [ p ]->giveID() );
        }
        unstructuredGrid->InsertNextCell(el->giveVTKCellType(), elindices);
    }

    unsigned i, j;
    size_t msize;
    vector< Vector >data;
    unsigned p;
    // ****************** cell data
    data.resize(selected_elements.size() );
    for ( p = 0; p < cell_data_size; p++ ) {
        msize = 1;
        i = 0;
        for ( vector< Element * > :: const_iterator ee = selected_elements.begin(); ee != selected_elements.end(); ++ee, i++ ) {
            ( * ee )->giveValues(codes [ p ].c_str(), data [ i ]);
            msize = max< size_t >(msize, data [ i ].size() );
        }
        vtkSmartPointer< vtkDoubleArray >cellDataArray = vtkSmartPointer< vtkDoubleArray > :: New();
        cellDataArray->SetName(codes [ p ].c_str() );
        cellDataArray->SetNumberOfComponents(msize);
        cellDataArray->SetNumberOfValues(selected_elements.size() * msize);
        i = 0;
        for ( vector< Vector > :: const_iterator d = data.begin(); d != data.end(); ++d, i++ ) {
            for ( j = 0; j < min< size_t >(msize, d->size() ); j++ ) {
                cellDataArray->SetValue(msize * i + j, ( * d ) [ j ]);
            }
            for ( ; j < msize; j++ ) {
                cellDataArray->SetValue(msize * i + j, 0);
            }
        }
        unstructuredGrid->GetCellData()->AddArray(cellDataArray);
    }


    // ****************** node data
    data.resize(nodes->giveSize() );
    for ( ; p < node_data_size + cell_data_size; p++ ) {
        msize = 1;
        i = 0;
        if ( codes [ p ].compare("residuals") == 0 ) {
            msize = 0;
            unsigned stDoF = 0;
            unsigned numDoF = 0;
            Vector res = solver->giveResiduals();
            for ( vector< Node * > :: const_iterator nn = nodes->begin(); nn != nodes->end(); ++nn, i++ ) {
                stDoF = ( * nn )->giveStartingDoF();
                numDoF = ( * nn )->giveNumberOfDoFs();
                data [ i ] = Vector(numDoF);
                msize = max< size_t >(msize, numDoF);
                for ( unsigned t = 0; t < numDoF; t++ ) {
                    data [ i ] [ t ] = res [ stDoF + t ];
                }
            }
        } else {
            for ( vector< Node * > :: const_iterator nn = nodes->begin(); nn != nodes->end(); ++nn, i++ ) {
                ( * nn )->giveValues(codes [ p ].c_str(), solver, data [ i ]);
                msize = max< size_t >(msize, data [ i ].size() );
            }
        }
        vtkSmartPointer< vtkDoubleArray >pointDataArray = vtkSmartPointer< vtkDoubleArray > :: New();
        pointDataArray->SetName(codes [ p ].c_str() );
        pointDataArray->SetNumberOfComponents(msize);
        pointDataArray->SetNumberOfValues(nodes->giveSize() * msize);
        i = 0;
        for ( vector< Vector > :: const_iterator d = data.begin(); d != data.end(); ++d, i++ ) {
            for ( j = 0; j < min< size_t >(msize, d->size() ); j++ ) {
                pointDataArray->SetValue(msize * i + j, ( * d ) [ j ]);
            }
            for ( ; j < msize; j++ ) {
                pointDataArray->SetValue(msize * i + j, 0);
            }
        }
        unstructuredGrid->GetPointData()->AddArray(pointDataArray);
    }

    // ****************** extrapolated node data
    for ( ; p < codes.size(); p++ ) {
        msize = 1;

        elems->extrapolateValuesFromIntegrationPointsToNodes(codes [ p ], data);
        for ( auto &v: data ) {
            msize = max< size_t >(msize, v.size() );
        }

        vtkSmartPointer< vtkDoubleArray >pointDataArray = vtkSmartPointer< vtkDoubleArray > :: New();
        pointDataArray->SetName(codes [ p ].c_str() );
        pointDataArray->SetNumberOfComponents(msize);
        pointDataArray->SetNumberOfValues(nodes->giveSize() * msize);
        i = 0;
        for ( vector< Vector > :: const_iterator d = data.begin(); d != data.end(); ++d, i++ ) {
            for ( j = 0; j < min< size_t >(msize, d->size() ); j++ ) {
                pointDataArray->SetValue(msize * i + j, ( * d ) [ j ]);
            }
            for ( ; j < msize; j++ ) {
                pointDataArray->SetValue(msize * i + j, 0);
            }
        }
        unstructuredGrid->GetPointData()->AddArray(pointDataArray);
    }

    //vtkNew<vtkXMLUnstructuredGridWriter> writer;
    vtkSmartPointer< vtkXMLUnstructuredGridWriter >writer = vtkSmartPointer< vtkXMLUnstructuredGridWriter > :: New();
    
    fs::path outputPath = resultDir / fname;
    writer->SetFileName( outputPath.string().c_str() );
    writer->SetInputData(unstructuredGrid);
    if ( binaryswitch ) {
        writer->SetDataModeToBinary();
    } else {
        writer->SetDataModeToAscii();
    }
    //writer->SetCompressorType();
    writer->Write();
#else
    cout << "VTK library not install, export of Rebars skipped" << endl;
#endif
}
