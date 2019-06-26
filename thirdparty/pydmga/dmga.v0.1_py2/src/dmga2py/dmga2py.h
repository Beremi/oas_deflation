/*
 * dmgapython.h
 *
 *  Created on: 23-04-2013
 *      Author: robson
 */

#ifndef DMGAPYTHON_H_
#define DMGAPYTHON_H_

#include <dmga/exceptions/exceptions.hpp>
#include <dmga/geometry/geometry.hpp>
#include <dmga/model/model.hpp>
#include <dmga/utils/utils.hpp>
#include <dmga/container/container.hpp>
#include <dmga/trajectory/trajectory.hpp>
#include <dmga/diagram/diagram.hpp>
#include <dmga/diagram/kinetic_diagram.hpp>
#include <dmga/shape/shape.hpp>
#include <dmga/graph/computation_graph.hpp>

typedef dmga::geometry::VoroppGeometry PythonGeometry;
typedef dmga::geometry::GeometryPreset PythonGeometryPreset;
typedef dmga::model::WeightedParticle<int> PythonParticle;
typedef dmga::container::VectorContainer<PythonGeometry, PythonParticle> PythonContainer;
typedef dmga::model::CellEx PythonCell;
typedef typename dmga::diagram::VoroppDiagram<PythonContainer, PythonCell> PythonDiagram;
typedef typename dmga::shape::VoroppAlphaShape<PythonDiagram> PythonAlphaShape;
typedef typename dmga::shape::AlphaShapeCell<PythonCell> PythonAlphaShapeCell;
typedef typename dmga::shape::AlphaComplex PythonAlphaComplex;
typedef typename dmga::shape::VoroppSASAShape<PythonDiagram> PythonSASAShape;
typedef typename dmga::shape::SASAShapeCell<PythonCell> PythonSASAShapeCell;
typedef typename dmga::shape::SASAContours PythonSASAContours;
typedef voro::CellIterator PythonCellIterator;
typedef PythonSASAContours::SASAContourIterator PythonSASAContourIterator;


#endif /* DMGAPYTHON_H_ */
