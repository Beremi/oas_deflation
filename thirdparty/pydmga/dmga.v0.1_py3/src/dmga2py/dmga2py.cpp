/*
 * dmga.cpp
 *
 * Python bindings C-file
 *
 *  Created on: 22-04-2013
 *      Author: robson
 */
#include <Python.h>

// #define DEBUG 0

#include <dmga/base.h>
#include "dmga2py.h"
#include <iostream>
#include <stdio.h>
#include <string>
#include <vector>

#define WITH_FIRST_HANDLE(s)

extern "C"{

static PyObject *
dmga2py_free_object(PyObject *self, PyObject *args)
{
    Py_ssize_t handle;
    if (!PyArg_ParseTuple(args, "n", &handle))
        return NULL;
    dmga::DMGAObject* object = reinterpret_cast<dmga::DMGAObject*>(handle);
    dmga::utils::safeDelete(object);
    return Py_BuildValue("i", 1);
}

static PyObject *
dmga2py_new_orthogonal_geometry(PyObject *self, PyObject *args)
{
    double x1;
    double y1;
    double z1;
    double x2;
	double y2;
	double z2;
    int x_per;
	int y_per;
	int z_per;
    if (!PyArg_ParseTuple(args, "ddddddiii",
    						&x1, &y1, &z1,
    						&x2, &y2, &z2,
    						&x_per,
    						&y_per,
    						&z_per))
        return NULL;

    dmga::geometry::OrthogonalGeometry* geom;
    geom = new dmga::geometry::OrthogonalGeometry(x1, y1, z1,
    											  x2, y2, z2,
    											  (bool)x_per,
    											  (bool)y_per,
    											  (bool)z_per);
    Py_ssize_t handle = reinterpret_cast<Py_ssize_t>(geom);
    return Py_BuildValue("n", handle);
}

static PyObject *
dmga2py_new_geometry_preset_sphere(PyObject *self, PyObject *args)
{
	double x, y, z, r;
    if (!PyArg_ParseTuple(args, "dddd", &x, &y, &z, &r))
        return NULL;

    PythonGeometryPreset* preset;
    preset = new dmga::geometry::SphereGeometryPreset(x, y, z, r);
    Py_ssize_t handle = reinterpret_cast<Py_ssize_t>(preset);
    return Py_BuildValue("n", handle);
}

static PyObject *
dmga2py_new_geometry_preset_one_cut_cylinder(PyObject *self, PyObject *args)
{
	double vx, vy, vz, x, y, z, r;
    if (!PyArg_ParseTuple(args, "ddddddd", &r, &x, &y, &z, &vx, &vy, &vz))
        return NULL;

    PythonGeometryPreset* preset;
    preset = new dmga::geometry::OneCutCylinderGeometryPreset(r, x, y, z, vx, vy, vz);
    Py_ssize_t handle = reinterpret_cast<Py_ssize_t>(preset);
    return Py_BuildValue("n", handle);
}

static PyObject *
dmga2py_new_geometry_preset_approximate_cylinder(PyObject *self, PyObject *args)
{
	double r, H;
	int n;
    if (!PyArg_ParseTuple(args, "ddi", &r, &H, &n))
        return NULL;

    //dmga::geometry::ApproximateCylinderGeometryPreset* preset;
    PythonGeometryPreset* preset;
    preset = new dmga::geometry::ApproximateCylinderGeometryPreset(r, H, n);
    Py_ssize_t handle = reinterpret_cast<Py_ssize_t>(preset);
    return Py_BuildValue("n", handle);
}

static PyObject *
dmga2py_add_geometry_preset(PyObject *self, PyObject *args)
{
	Py_ssize_t geo_handle, pres_handle;
    if (!PyArg_ParseTuple(args, "nn", &geo_handle, &pres_handle))
        return NULL;

    bool res = reinterpret_cast<PythonGeometry*>(geo_handle)->addPreset(reinterpret_cast<PythonGeometryPreset*>(pres_handle));
    return Py_BuildValue("i", (res ? 1 : 0));
}

static PyObject *
dmga2py_new_basic_container(PyObject *self, PyObject *args)
{
	Py_ssize_t geom_handle;
    if (!PyArg_ParseTuple(args, "n", &geom_handle))
        return NULL;

    PythonContainer* cont = new PythonContainer(*reinterpret_cast<PythonGeometry*>(geom_handle));
    Py_ssize_t handle = reinterpret_cast<Py_ssize_t>(cont);
    return Py_BuildValue("n", handle);
}

static PyObject *
dmga2py_new_diagram(PyObject *self, PyObject *args)
{
	Py_ssize_t container_handle;
    int cache_policy;
    if (!PyArg_ParseTuple(args, "ni", &container_handle, &cache_policy))
        return NULL;

    PythonDiagram* diag = new PythonDiagram(*reinterpret_cast<PythonContainer*>(container_handle), cache_policy);
    Py_ssize_t handle = reinterpret_cast<Py_ssize_t>(diag);
    return Py_BuildValue("n", handle);
}

static PyObject *
dmga2py_basic_container_add(PyObject *self, PyObject *args)
{
	Py_ssize_t container_handle;
    int id;
    double x, y, z, r;
    if (!PyArg_ParseTuple(args, "nidddd", &container_handle, &id, &x, &y, &z, &r))
        return NULL;

    PythonContainer* c = reinterpret_cast<PythonContainer*>(container_handle);
    c->add(PythonParticle(id, x, y, z, r));
    return Py_BuildValue("i", 1);
}

static PyObject *
dmga2py_basic_container_get_coords(PyObject *self, PyObject *args)
{
	Py_ssize_t container_handle;
    int num;
    if (!PyArg_ParseTuple(args, "ni", &container_handle, &num))
        return NULL;

    PythonContainer* c = reinterpret_cast<PythonContainer*>(container_handle);
    double* raw = c->getRawCoords(num);
    return Py_BuildValue("(dddd)", raw[0], raw[1], raw[2], raw[3]);
}

static PyObject *
dmga2py_basic_container_find(PyObject *self, PyObject *args)
{
	Py_ssize_t container_handle;
    int id;
    if (!PyArg_ParseTuple(args, "ni", &container_handle, &id))
        return NULL;

    PythonContainer* c = reinterpret_cast<PythonContainer*>(container_handle);
    PythonParticle data = c->find(id);
    return Py_BuildValue("(idddd)", data.key, data.getX(), data.getY(), data.getZ(), data.getR());
}

static PyObject *
dmga2py_basic_container_get(PyObject *self, PyObject *args)
{
	Py_ssize_t container_handle;
    int num;
    if (!PyArg_ParseTuple(args, "ni", &container_handle, &num))
        return NULL;

    PythonContainer* c = reinterpret_cast<PythonContainer*>(container_handle);
    PythonParticle data = c->get(num);
    return Py_BuildValue("(idddd)", data.key, data.getX(), data.getY(), data.getZ(), data.getR());
}

static PyObject *
dmga2py_diagram_get_cell(PyObject *self, PyObject *args)
{
	Py_ssize_t diagram_handle;
    int num;

    if (!PyArg_ParseTuple(args, "ni", &diagram_handle, &num))
        return NULL;

    PythonDiagram* d = reinterpret_cast<PythonDiagram*>(diagram_handle);
    PythonCell* ce = d->getCell(num);
    //TODO: tutaj trzeba sprawdzać, czy juz nie zwrócilismy i czy jest cached typ ustawiony w diagramie... :/
    Py_ssize_t cell_handle = reinterpret_cast<Py_ssize_t>(ce);

    return Py_BuildValue("n", cell_handle);
}

static PyObject *
dmga2py_cell_get_volume(PyObject *self, PyObject *args)
{
	Py_ssize_t cell_handle; if (!PyArg_ParseTuple(args, "n", &cell_handle)) return NULL;
    double vol = reinterpret_cast<PythonCell*>(cell_handle)->volume(); return Py_BuildValue("d", vol);
}

static PyObject *
dmga2py_cell_get_area(PyObject *self, PyObject *args)
{
	Py_ssize_t cell_handle; if (!PyArg_ParseTuple(args, "n", &cell_handle)) return NULL;
	double area = reinterpret_cast<PythonCell*>(cell_handle)->area(); return Py_BuildValue("d", area);
}

static PyObject *
dmga2py_cell_get_vertex_count(PyObject *self, PyObject *args)
{
	Py_ssize_t cell_handle; if (!PyArg_ParseTuple(args, "n", &cell_handle)) return NULL;
    int vertex_count = reinterpret_cast<PythonCell*>(cell_handle)->getVerticesCount(); return Py_BuildValue("i", vertex_count);
}

static PyObject *
dmga2py_cell_get_sides_count(PyObject *self, PyObject *args)
{
	Py_ssize_t cell_handle; if (!PyArg_ParseTuple(args, "n", &cell_handle)) return NULL;
    int sides_count = reinterpret_cast<PythonCell*>(cell_handle)->getSidesCount(); return Py_BuildValue("i", sides_count);
}

static PyObject *
dmga2py_cell_get_edges_count(PyObject *self, PyObject *args)
{
	Py_ssize_t cell_handle; if (!PyArg_ParseTuple(args, "n", &cell_handle)) return NULL;
    int sides_count = reinterpret_cast<PythonCell*>(cell_handle)->getEdgesCount(); return Py_BuildValue("i", sides_count);
}

static PyObject *
dmga2py_cell_get_vertex_coords(PyObject *self, PyObject *args) {
	Py_ssize_t cell_handle; int v; if (!PyArg_ParseTuple(args, "ni", &cell_handle, &v)) return NULL;
    double* coords = reinterpret_cast<PythonCell*>(cell_handle)->getVertexRaw(v); return Py_BuildValue("(ddd)", coords[0], coords[1], coords[2]);
}

static PyObject *
dmga2py_cell_edge_get_inverse(PyObject *self, PyObject *args) {
	Py_ssize_t cell_handle; int v, j; if (!PyArg_ParseTuple(args, "nii", &cell_handle, &v, &j)) return NULL;
	PythonCell* c = reinterpret_cast<PythonCell*>(cell_handle);
    int u, k; voro::edges::inv(c->voropp_cell, v, j, u, k); return Py_BuildValue("(ii)", u, k);
}

static PyObject *
dmga2py_cell_edge_get_next(PyObject *self, PyObject *args) {
	Py_ssize_t cell_handle; int v, j; if (!PyArg_ParseTuple(args, "nii", &cell_handle, &v, &j)) return NULL;
	PythonCell* c = reinterpret_cast<PythonCell*>(cell_handle);
    int u, k; voro::edges::next(c->voropp_cell, v, j, u, k); return Py_BuildValue("(ii)", u, k);
}

static PyObject *
dmga2py_cell_edge_get_previous(PyObject *self, PyObject *args) {
	Py_ssize_t cell_handle; int v, j; if (!PyArg_ParseTuple(args, "nii", &cell_handle, &v, &j)) return NULL;
	PythonCell* c = reinterpret_cast<PythonCell*>(cell_handle);
    int u, k; voro::edges::prev(c->voropp_cell, v, j, u, k); return Py_BuildValue("(ii)", u, k);
}

static PyObject *
dmga2py_cell_edge_get_cw(PyObject *self, PyObject *args) {
	Py_ssize_t cell_handle; int v, j; if (!PyArg_ParseTuple(args, "nii", &cell_handle, &v, &j)) return NULL;
	PythonCell* c = reinterpret_cast<PythonCell*>(cell_handle);
    int u, k; voro::edges::cw(c->voropp_cell, v, j, u, k); return Py_BuildValue("(ii)", u, k);
}

static PyObject *
dmga2py_cell_edge_get_ccw(PyObject *self, PyObject *args) {
	Py_ssize_t cell_handle; int v, j; if (!PyArg_ParseTuple(args, "nii", &cell_handle, &v, &j)) return NULL;
	PythonCell* c = reinterpret_cast<PythonCell*>(cell_handle);
    int u, k; voro::edges::ccw(c->voropp_cell, v, j, u, k); return Py_BuildValue("(ii)", u, k);
}

static PyObject *
dmga2py_cell_edge_get_neighbour_id(PyObject *self, PyObject *args) {
	Py_ssize_t cell_handle; int v, j; if (!PyArg_ParseTuple(args, "nii", &cell_handle, &v, &j)) return NULL;
	PythonCell* c = reinterpret_cast<PythonCell*>(cell_handle);
	int neighbour_id = c->getNeighbourId(v, j); return Py_BuildValue("i", neighbour_id);
}

static PyObject *
dmga2py_cell_iterator(PyObject *self, PyObject *args)
{
	Py_ssize_t cell_handle; if (!PyArg_ParseTuple(args, "n", &cell_handle)) return NULL;

    voro::CellIterator* cell_iterator_obj = reinterpret_cast<PythonCell*>(cell_handle)->newCellIterator();
    Py_ssize_t cell_it_handle = reinterpret_cast<Py_ssize_t>(cell_iterator_obj);
    return Py_BuildValue("n", cell_it_handle);
}

static PyObject *
dmga2py_cell_iterator_current_vertex(PyObject *self, PyObject *args) {
	Py_ssize_t cell_it_handle; if (!PyArg_ParseTuple(args, "n", &cell_it_handle)) return NULL;
    int v = reinterpret_cast<voro::CellIterator*>(cell_it_handle)->current().v; return Py_BuildValue("i", v);
}

static PyObject *
dmga2py_cell_iterator_next_vertex(PyObject *self, PyObject *args) {
	Py_ssize_t cell_it_handle; if (!PyArg_ParseTuple(args, "n", &cell_it_handle)) return NULL;
    int v = reinterpret_cast<voro::CellIterator*>(cell_it_handle)->next().v; return Py_BuildValue("i", v);
}

static PyObject *
dmga2py_cell_iterator_prev_vertex(PyObject *self, PyObject *args) {
	Py_ssize_t cell_it_handle; if (!PyArg_ParseTuple(args, "n", &cell_it_handle)) return NULL;
    int v = reinterpret_cast<voro::CellIterator*>(cell_it_handle)->prev().v; return Py_BuildValue("i", v);
}

static PyObject *
dmga2py_cell_iterator_current_edge(PyObject *self, PyObject *args) {
	Py_ssize_t cell_it_handle; if (!PyArg_ParseTuple(args, "n", &cell_it_handle)) return NULL;
    int j = reinterpret_cast<voro::CellIterator*>(cell_it_handle)->current().j; return Py_BuildValue("i", j);
}

static PyObject *
dmga2py_cell_iterator_inverse_edge(PyObject *self, PyObject *args) {
	Py_ssize_t cell_it_handle; if (!PyArg_ParseTuple(args, "n", &cell_it_handle)) return NULL;
    int j = reinterpret_cast<voro::CellIterator*>(cell_it_handle)->inv().j; return Py_BuildValue("i", j);
}

static PyObject *
dmga2py_cell_iterator_current(PyObject *self, PyObject *args) {
	Py_ssize_t cell_it_handle; if (!PyArg_ParseTuple(args, "n", &cell_it_handle)) return NULL;
    voro::EdgeDesc a = reinterpret_cast<voro::CellIterator*>(cell_it_handle)->current(); return Py_BuildValue("(ii)", a.v, a.j);
}

static PyObject *
dmga2py_cell_iterator_inverse(PyObject *self, PyObject *args) {
	Py_ssize_t cell_it_handle; if (!PyArg_ParseTuple(args, "n", &cell_it_handle)) return NULL;
    voro::EdgeDesc a = reinterpret_cast<voro::CellIterator*>(cell_it_handle)->inv(); return Py_BuildValue("(ii)", a.v, a.j);
}

static PyObject *
dmga2py_cell_iterator_vertex_coords(PyObject *self, PyObject *args) {
	Py_ssize_t cell_it_handle; int v; if (!PyArg_ParseTuple(args, "ni", &cell_it_handle, &v)) return NULL;
    double* coords = reinterpret_cast<voro::CellIterator*>(cell_it_handle)->vertexRaw(v); return Py_BuildValue("(ddd)", coords[0], coords[1], coords[2]);
}

static PyObject *
dmga2py_cell_iterator_current_vertex_coords(PyObject *self, PyObject *args) {
	Py_ssize_t cell_it_handle; if (!PyArg_ParseTuple(args, "n", &cell_it_handle)) return NULL;
    double* coords = reinterpret_cast<voro::CellIterator*>(cell_it_handle)->vertexRaw(); return Py_BuildValue("(ddd)", coords[0], coords[1], coords[2]);
}

static PyObject *
dmga2py_cell_iterator_current_side_neighbour_id(PyObject *self, PyObject *args) {
	Py_ssize_t cell_handle, cell_it_handle; if (!PyArg_ParseTuple(args, "nn", &cell_handle, &cell_it_handle)) return NULL;
    int j = reinterpret_cast<voro::CellIterator*>(cell_it_handle)->current().j;
    int v = reinterpret_cast<voro::CellIterator*>(cell_it_handle)->current().v;
    int nid = reinterpret_cast<PythonCell*>(cell_handle)->getNeighbourId(v, j);
    return Py_BuildValue("i", nid);
}

static PyObject *
dmga2py_cell_iterator_reset(PyObject *self, PyObject *args) {
	Py_ssize_t cell_it_handle; if (!PyArg_ParseTuple(args, "n", &cell_it_handle)) return NULL;
	reinterpret_cast<voro::CellIterator*>(cell_it_handle)->reset(); return Py_BuildValue("i", 1);
}

static PyObject *
dmga2py_cell_iterator_is_finished(PyObject *self, PyObject *args) {
	Py_ssize_t cell_it_handle; if (!PyArg_ParseTuple(args, "n", &cell_it_handle)) return NULL;
    int state = reinterpret_cast<voro::CellIterator*>(cell_it_handle)->isFinished() ? 1 : 0; return Py_BuildValue("i", state);
}

static PyObject *
dmga2py_cell_iterator_jump(PyObject *self, PyObject *args) {
	Py_ssize_t cell_it_handle; if (!PyArg_ParseTuple(args, "n", &cell_it_handle)) return NULL;
    int state = reinterpret_cast<voro::CellIterator*>(cell_it_handle)->jump() ? 1 : 0; return Py_BuildValue("i", state);
}

static PyObject *
dmga2py_cell_iterator_forward(PyObject *self, PyObject *args) {
	Py_ssize_t cell_it_handle; if (!PyArg_ParseTuple(args, "n", &cell_it_handle)) return NULL;
	reinterpret_cast<voro::CellIterator*>(cell_it_handle)->forward(); return Py_BuildValue("i", 1);
}

static PyObject *
dmga2py_cell_iterator_backward(PyObject *self, PyObject *args) {
	Py_ssize_t cell_it_handle; if (!PyArg_ParseTuple(args, "n", &cell_it_handle)) return NULL;
	reinterpret_cast<voro::CellIterator*>(cell_it_handle)->backward(); return Py_BuildValue("i", 1);
}

static PyObject *
dmga2py_cell_iterator_reverse(PyObject *self, PyObject *args) {
	Py_ssize_t cell_it_handle; if (!PyArg_ParseTuple(args, "n", &cell_it_handle)) return NULL;
	reinterpret_cast<voro::CellIterator*>(cell_it_handle)->reverse(); return Py_BuildValue("i", 1);
}

static PyObject *
dmga2py_cell_iterator_is_marked(PyObject *self, PyObject *args) {
	Py_ssize_t cell_it_handle; if (!PyArg_ParseTuple(args, "n", &cell_it_handle)) return NULL;
    int state = reinterpret_cast<voro::CellIterator*>(cell_it_handle)->isMarked() ? 1 : 0; return Py_BuildValue("i", state);
}

static PyObject *
dmga2py_cell_iterator_mark(PyObject *self, PyObject *args) {
	Py_ssize_t cell_it_handle; if (!PyArg_ParseTuple(args, "n", &cell_it_handle)) return NULL;
	reinterpret_cast<voro::CellIterator*>(cell_it_handle)->mark(); return Py_BuildValue("i", 1);
}

static PyObject *
dmga2py_cell_iterator_unmark(PyObject *self, PyObject *args) {
	Py_ssize_t cell_it_handle; if (!PyArg_ParseTuple(args, "n", &cell_it_handle)) return NULL;
	reinterpret_cast<voro::CellIterator*>(cell_it_handle)->unmark(); return Py_BuildValue("i", 1);
}

static PyObject *
dmga2py_new_alpha_shape(PyObject *self, PyObject *args)
{
	Py_ssize_t diagram_handle; if (!PyArg_ParseTuple(args, "n", &diagram_handle)) return NULL;

    PythonAlphaShape* shape;
    shape = new PythonAlphaShape(*reinterpret_cast<PythonDiagram*>(diagram_handle));
    Py_ssize_t handle = reinterpret_cast<Py_ssize_t>(shape);
    return Py_BuildValue("n", handle);
}

static PyObject *
dmga2py_alpha_shape_get_shape_cell(PyObject *self, PyObject *args)
{
	Py_ssize_t shape_handle;
    int num;
    if (!PyArg_ParseTuple(args, "ni", &shape_handle, &num))
        return NULL;

    PythonAlphaShape* s = reinterpret_cast<PythonAlphaShape*>(shape_handle);
    PythonAlphaShapeCell* ce = s->getCellShape(num);
    //TODO: tutaj trzeba sprawdzać, czy juz nie zwrócilismy i czy jest cached typ ustawiony w diagramie... :/
    Py_ssize_t shape_cell_handle = reinterpret_cast<Py_ssize_t>(ce);
    Py_ssize_t cell_handle = reinterpret_cast<Py_ssize_t>(ce->cell);
    Py_ssize_t complex_handle = reinterpret_cast<Py_ssize_t>(ce->complex);

    return Py_BuildValue("(nnn)", shape_cell_handle, cell_handle, complex_handle);
}

static PyObject *
dmga2py_alpha_shape_cell_as_tuple(PyObject *self, PyObject *args)
{
	Py_ssize_t alpha_shape_cell_handle;
    if (!PyArg_ParseTuple(args, "n", &alpha_shape_cell_handle))
        return NULL;

    PythonAlphaShapeCell* c = reinterpret_cast<PythonAlphaShapeCell*>(alpha_shape_cell_handle);
    Py_ssize_t cell_handle = reinterpret_cast<Py_ssize_t>(c->cell);
    Py_ssize_t complex_handle = reinterpret_cast<Py_ssize_t>(c->complex);

    return Py_BuildValue("(nn)", cell_handle, complex_handle);
}

static PyObject *
dmga2py_alpha_shape_get_complex(PyObject *self, PyObject *args)
{
	Py_ssize_t shape_handle;
    int num;
    if (!PyArg_ParseTuple(args, "ni", &shape_handle, &num))
        return NULL;

    PythonAlphaShape* s = reinterpret_cast<PythonAlphaShape*>(shape_handle);
    PythonAlphaComplex* cmplx = s->getAlphaComplex(num);
    //TODO: tutaj trzeba sprawdzać, czy juz nie zwrócilismy i czy jest cached typ ustawiony w diagramie... :/
    Py_ssize_t complex_handle = reinterpret_cast<Py_ssize_t>(cmplx);

    return Py_BuildValue("n", complex_handle);
}

static PyObject *
dmga2py_alpha_shape_get_max_alpha_threshold(PyObject *self, PyObject *args)
{
	Py_ssize_t shape_handle;

    if (!PyArg_ParseTuple(args, "n", &shape_handle))
        return NULL;

    PythonAlphaShape* s = reinterpret_cast<PythonAlphaShape*>(shape_handle);

    return Py_BuildValue("d", s->max_alpha_threshold);
}

static PyObject *
dmga2py_alpha_shape_cell_get_max_alpha_threshold(PyObject *self, PyObject *args)
{
	Py_ssize_t alpha_shape_cell_handle;

    if (!PyArg_ParseTuple(args, "n", &alpha_shape_cell_handle))
        return NULL;

    PythonAlphaShapeCell* cell = reinterpret_cast<PythonAlphaShapeCell*>(alpha_shape_cell_handle);

    return Py_BuildValue("d", cell->complex->max_alpha_threshold);
}

static PyObject *
dmga2py_alpha_shape_complex_get_max_alpha_threshold(PyObject *self, PyObject *args)
{
	Py_ssize_t alpha_shape_complex_handle;

    if (!PyArg_ParseTuple(args, "n", &alpha_shape_complex_handle))
        return NULL;

    PythonAlphaComplex* cmplx = reinterpret_cast<PythonAlphaComplex*>(alpha_shape_complex_handle);

    return Py_BuildValue("d", cmplx->max_alpha_threshold);
}

static PyObject *
dmga2py_alpha_shape_complex_print(PyObject *self, PyObject *args)
{
	Py_ssize_t alpha_shape_complex_handle;

    if (!PyArg_ParseTuple(args, "n", &alpha_shape_complex_handle))
        return NULL;

    PythonAlphaComplex* cmplx = reinterpret_cast<PythonAlphaComplex*>(alpha_shape_complex_handle);
    std::cout << alpha_shape_complex_handle << " complex print: \n";
    std::cout << (*cmplx) << "\n";
    std::cout << alpha_shape_complex_handle << " complex print finished.\n";

    return Py_BuildValue("i", 1);
}

static PyObject *
dmga2py_alpha_shape_complex_upper_bound(PyObject *self, PyObject *args)
{
	Py_ssize_t alpha_shape_complex_handle;
	double alpha;
    if (!PyArg_ParseTuple(args, "nd", &alpha_shape_complex_handle, &alpha))
        return NULL;

    PythonAlphaComplex* cmplx = reinterpret_cast<PythonAlphaComplex*>(alpha_shape_complex_handle);
    return Py_BuildValue("i", cmplx->end(alpha) - cmplx->begin());
}

static PyObject *
dmga2py_alpha_shape_complex_get_simplex(PyObject *self, PyObject *args)
{
	Py_ssize_t alpha_shape_complex_handle;
	int number;
    if (!PyArg_ParseTuple(args, "ni", &alpha_shape_complex_handle, &number))
        return NULL;

    PythonAlphaComplex* cmplx = reinterpret_cast<PythonAlphaComplex*>(alpha_shape_complex_handle);
    dmga::shape::AlphaSimplex& s = cmplx->item(number);
    //we return dim,v,  j,  alpha, id, is_killed
    //          int,int,int,double,int,bool(as int)
    return Py_BuildValue("(iiidii)", s.dim, s.v, s.j, s.alpha_threshold, s.id, s.is_killed);
}

static PyObject *
dmga2py_new_sasa_shape(PyObject *self, PyObject *args)
{
	Py_ssize_t diagram_handle; if (!PyArg_ParseTuple(args, "n", &diagram_handle)) return NULL;

    PythonSASAShape* shape;
    shape = new PythonSASAShape(*reinterpret_cast<PythonDiagram*>(diagram_handle));
    Py_ssize_t handle = reinterpret_cast<Py_ssize_t>(shape);
    return Py_BuildValue("n", handle);
}

static PyObject *
dmga2py_sasa_shape_get_shape_cell(PyObject *self, PyObject *args)
{
	Py_ssize_t shape_handle;
    int num;
    if (!PyArg_ParseTuple(args, "ni", &shape_handle, &num))
        return NULL;

    PythonSASAShape* s = reinterpret_cast<PythonSASAShape*>(shape_handle);
    PythonSASAShapeCell* ce = s->getCellShape(num);
    //TODO: tutaj trzeba sprawdzać, czy juz nie zwrócilismy i czy jest cached typ ustawiony w diagramie... :/
    Py_ssize_t shape_cell_handle = reinterpret_cast<Py_ssize_t>(ce);
    Py_ssize_t cell_handle = reinterpret_cast<Py_ssize_t>(ce->cell);
    Py_ssize_t contours_handle = reinterpret_cast<Py_ssize_t>(ce->contours);

    return Py_BuildValue("(nnn)", shape_cell_handle, cell_handle, contours_handle);
}

static PyObject *
dmga2py_sasa_shape_cell_as_tuple(PyObject *self, PyObject *args)
{
	Py_ssize_t sasa_cell_handle;
    if (!PyArg_ParseTuple(args, "n", &sasa_cell_handle))
        return NULL;

    PythonSASAShapeCell* ce = reinterpret_cast<PythonSASAShapeCell*>(sasa_cell_handle);
    Py_ssize_t cell_handle = reinterpret_cast<Py_ssize_t>(ce->cell);
    Py_ssize_t contours_handle = reinterpret_cast<Py_ssize_t>(ce->contours);

    return Py_BuildValue("(nn)", cell_handle, contours_handle);
}

static PyObject *
dmga2py_sasa_shape_get_area(PyObject *self, PyObject *args)
{
	Py_ssize_t shape_handle; if (!PyArg_ParseTuple(args, "n", &shape_handle)) return NULL;
    PythonSASAShape* s = reinterpret_cast<PythonSASAShape*>(shape_handle);
    return Py_BuildValue("d", s->area());
}

static PyObject *
dmga2py_sasa_shape_cell_get_area(PyObject *self, PyObject *args)
{
	Py_ssize_t shape_cell_handle; if (!PyArg_ParseTuple(args, "n", &shape_cell_handle)) return NULL;
	PythonSASAShapeCell* c = reinterpret_cast<PythonSASAShapeCell*>(shape_cell_handle);
    return Py_BuildValue("d", c->area());
}

static PyObject *
dmga2py_sasa_shape_cell_get_excluded_area(PyObject *self, PyObject *args)
{
	Py_ssize_t shape_cell_handle; if (!PyArg_ParseTuple(args, "n", &shape_cell_handle)) return NULL;
	PythonSASAShapeCell* c = reinterpret_cast<PythonSASAShapeCell*>(shape_cell_handle);
    return Py_BuildValue("d", c->contours->excluded_area);
}

static PyObject *
dmga2py_sasa_shape_cell_get_volume(PyObject *self, PyObject *args)
{
	Py_ssize_t shape_cell_handle; if (!PyArg_ParseTuple(args, "n", &shape_cell_handle)) return NULL;
	PythonSASAShapeCell* c = reinterpret_cast<PythonSASAShapeCell*>(shape_cell_handle);
    return Py_BuildValue("d", c->volume());
}

static PyObject *
dmga2py_sasa_shape_cell_get_excluded_volume(PyObject *self, PyObject *args)
{
	Py_ssize_t shape_cell_handle; if (!PyArg_ParseTuple(args, "n", &shape_cell_handle)) return NULL;
	PythonSASAShapeCell* c = reinterpret_cast<PythonSASAShapeCell*>(shape_cell_handle);
    return Py_BuildValue("d", c->contours->excluded_volume);
}


/////SASAContourIterator
static PyObject *
dmga2py_sasa_contours_new_border_iterator(PyObject *self, PyObject *args)
{
	Py_ssize_t contour_handle; if (!PyArg_ParseTuple(args, "n", &contour_handle)) return NULL;
	PythonSASAContours* c = reinterpret_cast<PythonSASAContours*>(contour_handle);
    return Py_BuildValue("n", reinterpret_cast<Py_ssize_t>(new PythonSASAContours::SASABorderIterator(c)));
}

static PyObject *
dmga2py_sasa_contours_new_polygon_iterator(PyObject *self, PyObject *args)
{
	Py_ssize_t contour_handle; if (!PyArg_ParseTuple(args, "n", &contour_handle)) return NULL;
	PythonSASAContours* c = reinterpret_cast<PythonSASAContours*>(contour_handle);
	return Py_BuildValue("n", reinterpret_cast<Py_ssize_t>(new PythonSASAContours::SASAPolygonIterator(c)));
}

static PyObject *
dmga2py_sasa_contours_iterator_next_vertex(PyObject *self, PyObject *args)
{
	Py_ssize_t iterator_handle; int vertex; if (!PyArg_ParseTuple(args, "ni", &iterator_handle, &vertex)) return NULL;
	PythonSASAContourIterator* cit = reinterpret_cast<PythonSASAContourIterator*>(iterator_handle);
	return Py_BuildValue("i", cit->nextVertex(vertex));
}

static PyObject *
dmga2py_sasa_contours_iterator_next_center(PyObject *self, PyObject *args)
{
	Py_ssize_t iterator_handle; int vertex; if (!PyArg_ParseTuple(args, "ni", &iterator_handle, &vertex)) return NULL;
	PythonSASAContourIterator* cit = reinterpret_cast<PythonSASAContourIterator*>(iterator_handle);
	return Py_BuildValue("i", cit->nextCenter(vertex));
}

static PyObject *
dmga2py_sasa_contours_iterator_is_marked(PyObject *self, PyObject *args)
{
	Py_ssize_t iterator_handle; int vertex; if (!PyArg_ParseTuple(args, "ni", &iterator_handle, &vertex)) return NULL;
	PythonSASAContourIterator* cit = reinterpret_cast<PythonSASAContourIterator*>(iterator_handle);
	return Py_BuildValue("i", cit->isMarked(vertex) ? 1 : 0);
}

static PyObject *
dmga2py_sasa_contours_iterator_is_marked_current(PyObject *self, PyObject *args)
{
	Py_ssize_t iterator_handle; if (!PyArg_ParseTuple(args, "n", &iterator_handle)) return NULL;
	PythonSASAContourIterator* cit = reinterpret_cast<PythonSASAContourIterator*>(iterator_handle);
	return Py_BuildValue("i", cit->isMarked() ? 1 : 0);
}

static PyObject *
dmga2py_sasa_contours_iterator_mark(PyObject *self, PyObject *args)
{
	Py_ssize_t iterator_handle; int vertex; if (!PyArg_ParseTuple(args, "ni", &iterator_handle, &vertex)) return NULL;
	PythonSASAContourIterator* cit = reinterpret_cast<PythonSASAContourIterator*>(iterator_handle);
	cit->mark(vertex);
	return Py_BuildValue("i", 1);
}

static PyObject *
dmga2py_sasa_contours_iterator_mark_current(PyObject *self, PyObject *args)
{
	Py_ssize_t iterator_handle; if (!PyArg_ParseTuple(args, "n", &iterator_handle)) return NULL;
	PythonSASAContourIterator* cit = reinterpret_cast<PythonSASAContourIterator*>(iterator_handle);
	cit->mark();
	return Py_BuildValue("i", 1);
}

static PyObject *
dmga2py_sasa_contours_iterator_unmark(PyObject *self, PyObject *args)
{
	Py_ssize_t iterator_handle; int vertex; if (!PyArg_ParseTuple(args, "ni", &iterator_handle, &vertex)) return NULL;
	PythonSASAContourIterator* cit = reinterpret_cast<PythonSASAContourIterator*>(iterator_handle);
	cit->unmark(vertex);
	return Py_BuildValue("i", 1);
}

static PyObject *
dmga2py_sasa_contours_iterator_unmark_current(PyObject *self, PyObject *args)
{
	Py_ssize_t iterator_handle; if (!PyArg_ParseTuple(args, "n", &iterator_handle)) return NULL;
	PythonSASAContourIterator* cit = reinterpret_cast<PythonSASAContourIterator*>(iterator_handle);
	cit->unmark();
	return Py_BuildValue("i", 1);
}

static PyObject *
dmga2py_sasa_contours_iterator_current_vertex_coords(PyObject *self, PyObject *args)
{
	Py_ssize_t iterator_handle; if (!PyArg_ParseTuple(args, "n", &iterator_handle)) return NULL;
	PythonSASAContourIterator* cit = reinterpret_cast<PythonSASAContourIterator*>(iterator_handle);
	const dmga::model::Point3D& p = cit->currentPoint();
	return Py_BuildValue("(ddd)", p.getX(), p.getY(), p.getZ());
}

static PyObject *
dmga2py_sasa_contours_iterator_next_vertex_coords(PyObject *self, PyObject *args)
{
	Py_ssize_t iterator_handle; if (!PyArg_ParseTuple(args, "n", &iterator_handle)) return NULL;
	PythonSASAContourIterator* cit = reinterpret_cast<PythonSASAContourIterator*>(iterator_handle);
	const dmga::model::Point3D& p = cit->nextPoint();
	return Py_BuildValue("(ddd)", p.getX(), p.getY(), p.getZ());
}

static PyObject *
dmga2py_sasa_contours_iterator_center_on_plane_coords(PyObject *self, PyObject *args)
{
	Py_ssize_t iterator_handle; if (!PyArg_ParseTuple(args, "n", &iterator_handle)) return NULL;
	PythonSASAContourIterator* cit = reinterpret_cast<PythonSASAContourIterator*>(iterator_handle);
	const dmga::model::Point3D& p = cit->centerPointOnPlane();
	return Py_BuildValue("(ddd)", p.getX(), p.getY(), p.getZ());
}

static PyObject *
dmga2py_sasa_contours_iterator_center_on_sphere_coords(PyObject *self, PyObject *args)
{
	Py_ssize_t iterator_handle; if (!PyArg_ParseTuple(args, "n", &iterator_handle)) return NULL;
	PythonSASAContourIterator* cit = reinterpret_cast<PythonSASAContourIterator*>(iterator_handle);
	const dmga::model::Point3D& p = cit->centerPointOnSphere();
	return Py_BuildValue("(ddd)", p.getX(), p.getY(), p.getZ());
}

static PyObject *
dmga2py_sasa_contours_iterator_ready(PyObject *self, PyObject *args)
{
	Py_ssize_t iterator_handle; if (!PyArg_ParseTuple(args, "n", &iterator_handle)) return NULL;
	PythonSASAContourIterator* cit = reinterpret_cast<PythonSASAContourIterator*>(iterator_handle);
	return Py_BuildValue("i", cit->ready() ? 1 : 0);
}

static PyObject *
dmga2py_sasa_contours_iterator_reset(PyObject *self, PyObject *args)
{
	Py_ssize_t iterator_handle; if (!PyArg_ParseTuple(args, "n", &iterator_handle)) return NULL;
	PythonSASAContourIterator* cit = reinterpret_cast<PythonSASAContourIterator*>(iterator_handle);
	return Py_BuildValue("i", cit->reset() ? 1 : 0);
}

static PyObject *
dmga2py_sasa_contours_iterator_forward(PyObject *self, PyObject *args)
{
	Py_ssize_t iterator_handle; if (!PyArg_ParseTuple(args, "n", &iterator_handle)) return NULL;
	PythonSASAContourIterator* cit = reinterpret_cast<PythonSASAContourIterator*>(iterator_handle);
	return Py_BuildValue("i", cit->forward());
}

static PyObject *
dmga2py_sasa_contours_iterator_jump(PyObject *self, PyObject *args)
{
	Py_ssize_t iterator_handle; if (!PyArg_ParseTuple(args, "n", &iterator_handle)) return NULL;
	PythonSASAContourIterator* cit = reinterpret_cast<PythonSASAContourIterator*>(iterator_handle);
	return Py_BuildValue("i", cit->jump());
}

static PyObject *
dmga2py_sasa_contours_iterator_iterate(PyObject *self, PyObject *args)
{
	Py_ssize_t iterator_handle; if (!PyArg_ParseTuple(args, "n", &iterator_handle)) return NULL;
	PythonSASAContourIterator* cit = reinterpret_cast<PythonSASAContourIterator*>(iterator_handle);
	return Py_BuildValue("i", cit->iterate());
}

static PyObject *
dmga2py_sasa_contours_iterator_as_tuple(PyObject *self, PyObject *args)
{
	Py_ssize_t iterator_handle; if (!PyArg_ParseTuple(args, "n", &iterator_handle)) return NULL;
	PythonSASAContourIterator* cit = reinterpret_cast<PythonSASAContourIterator*>(iterator_handle);
	return Py_BuildValue("(iii)", cit->current, cit->center, cit->next);
}

static PyObject *
dmga2py_sasa_contours_iterator_as_tuple_coords(PyObject *self, PyObject *args)
{
	Py_ssize_t iterator_handle; if (!PyArg_ParseTuple(args, "n", &iterator_handle)) return NULL;
	PythonSASAContourIterator* cit = reinterpret_cast<PythonSASAContourIterator*>(iterator_handle);
	const dmga::model::Point3D& first = cit->currentPoint();
	const dmga::model::Point3D& on_plane = cit->centerPointOnPlane();
	const dmga::model::Point3D& on_sphere = cit->centerPointOnSphere();
	const dmga::model::Point3D& second = cit->nextPoint();
	return Py_BuildValue("((ddd)(ddd)(ddd)(ddd))",
							first.getX(), first.getY(), first.getZ(),
							on_plane.getX(), on_plane.getY(), on_plane.getZ(),
							on_sphere.getX(), on_sphere.getY(), on_sphere.getZ(),
							second.getX(), second.getY(), second.getZ()
						);
}



static PyMethodDef Dmga2pyMethods[] = {
    {
		"free_object",
		dmga2py_free_object,
		METH_VARARGS,
		"Given a handle (to any kind of object of type DMGAObject) it frees the memory."
    },
    {
		"new_orthogonal_geometry",
		dmga2py_new_orthogonal_geometry,
		METH_VARARGS,
		"Params: x_size, y_size, z_size, x_periodic, y_periodic, z_periodic. Creates a new OrthogonalGeometry object and returns its handle."
    },
    {
		"new_geometry_preset_sphere",
		dmga2py_new_geometry_preset_sphere,
		METH_VARARGS,
		"TODO: doc"
    },
    {
		"new_geometry_preset_one_cut_cylinder",
		dmga2py_new_geometry_preset_one_cut_cylinder,
		METH_VARARGS,
		"TODO: doc"
	},
	{
		"new_geometry_preset_approximate_cylinder",
		dmga2py_new_geometry_preset_approximate_cylinder,
		METH_VARARGS,
		"TODO: doc"
	},
    {
		"add_geometry_preset",
		dmga2py_add_geometry_preset,
		METH_VARARGS,
		"TODO: doc"
	},
    {
		"new_basic_container",
		dmga2py_new_basic_container,
		METH_VARARGS,
		"Params: geometry_handle. Creates a new BasicContainer object and returns its handle."
    },
    {
		"new_diagram",
		dmga2py_new_diagram,
		METH_VARARGS,
		"Params: container_handle. Creates a new Diagram object and returns its handle."
    },
    {
		"basic_container_get",
		dmga2py_basic_container_get,
		METH_VARARGS,
		"TODO: doc"
    },
    {
		"basic_container_find",
		dmga2py_basic_container_find,
		METH_VARARGS,
		"TODO: doc"
    },
    {
		"basic_container_get_coords",
		dmga2py_basic_container_get_coords,
		METH_VARARGS,
		"TODO: doc"
    },
    {
		"basic_container_add",
		dmga2py_basic_container_add,
		METH_VARARGS,
		"TODO: doc"
    },
    {
		"diagram_get_cell",
		dmga2py_diagram_get_cell,
		METH_VARARGS,
		"TODO: doc"
    },
    {
		"cell_get_volume",
		dmga2py_cell_get_volume,
		METH_VARARGS,
		"TODO: doc"
    },
    {
		"cell_get_area",
		dmga2py_cell_get_area,
		METH_VARARGS,
		"TODO: doc"
    },
    {
		"cell_get_vertex_count",
		dmga2py_cell_get_vertex_count,
		METH_VARARGS,
		"TODO: doc"
	},
	{
		"cell_get_sides_count",
		dmga2py_cell_get_sides_count,
		METH_VARARGS,
		"TODO: doc"
	},
	{
		"cell_get_edges_count",
		dmga2py_cell_get_edges_count,
		METH_VARARGS,
		"TODO: doc"
	},
	{
		"cell_get_vertex_coords",
		dmga2py_cell_get_vertex_coords,
		METH_VARARGS,
		"TODO: doc"
	},
	{
		"cell_edge_get_inverse",
		dmga2py_cell_edge_get_inverse,
		METH_VARARGS,
		"TODO: doc"
	},
	{
		"cell_edge_get_next",
		dmga2py_cell_edge_get_next,
		METH_VARARGS,
		"TODO: doc"
	},
	{
		"cell_edge_get_previous",
		dmga2py_cell_edge_get_previous,
		METH_VARARGS,
		"TODO: doc"
	},
	{
		"cell_edge_get_cw",
		dmga2py_cell_edge_get_cw,
		METH_VARARGS,
		"TODO: doc"
	},
	{
		"cell_edge_get_ccw",
		dmga2py_cell_edge_get_ccw,
		METH_VARARGS,
		"TODO: doc"
	},
	{
		"cell_edge_get_neighbour_id",
		dmga2py_cell_edge_get_neighbour_id,
		METH_VARARGS,
		"TODO: doc"
	},
    {
    	"cell_iterator",
    	dmga2py_cell_iterator,
    	METH_VARARGS,
    	"TODO: doc"
    },
    {
    	"cell_iterator_current_vertex",
    	dmga2py_cell_iterator_current_vertex,
    	METH_VARARGS,
    	"TODO: doc"
    },
    {
    	"cell_iterator_next_vertex",
    	dmga2py_cell_iterator_next_vertex,
    	METH_VARARGS,
    	"TODO: doc"
    },
    {
		"cell_iterator_prev_vertex",
		dmga2py_cell_iterator_prev_vertex,
		METH_VARARGS,
		"TODO: doc"
    },
    {
		"cell_iterator_current_edge",
		dmga2py_cell_iterator_current_edge,
		METH_VARARGS,
		"TODO: doc"
    },
    {
		"cell_iterator_inverse_edge",
		dmga2py_cell_iterator_inverse_edge,
		METH_VARARGS,
		"TODO: doc"
    },
    {
		"cell_iterator_inverse",
		dmga2py_cell_iterator_inverse,
		METH_VARARGS,
		"TODO: doc"
    },
    {
		"cell_iterator_current",
		dmga2py_cell_iterator_current,
		METH_VARARGS,
		"TODO: doc"
    },
    {
    	"cell_iterator_vertex_coords",
    	dmga2py_cell_iterator_vertex_coords,
    	METH_VARARGS,
    	"TODO: doc"
    },
    {
    	"cell_iterator_current_vertex_coords",
    	dmga2py_cell_iterator_current_vertex_coords,
    	METH_VARARGS,
    	"TODO: doc"
    },
    {
		"cell_iterator_current_side_neighbour_id",
		dmga2py_cell_iterator_current_side_neighbour_id,
		METH_VARARGS,
		"TODO: doc"
	},
    {
    	"cell_iterator_reset",
    	dmga2py_cell_iterator_reset,
    	METH_VARARGS,
    	"TODO: doc"
    },
    {
    	"cell_iterator_is_finished",
    	dmga2py_cell_iterator_is_finished,
    	METH_VARARGS,
    	"TODO: doc"
    },
    {
    	"cell_iterator_jump",
    	dmga2py_cell_iterator_jump,
    	METH_VARARGS,
    	"TODO: doc"
    },
    {
    	"cell_iterator_forward",
    	dmga2py_cell_iterator_forward,
    	METH_VARARGS,
    	"TODO: doc"
    },
    {
    	"cell_iterator_backward",
    	dmga2py_cell_iterator_backward,
    	METH_VARARGS,
    	"TODO: doc"
    },
    {
    	"cell_iterator_reverse",
    	dmga2py_cell_iterator_reverse,
    	METH_VARARGS,
    	"TODO: doc"
    },
    {
    	"cell_iterator_is_marked",
    	dmga2py_cell_iterator_is_marked,
    	METH_VARARGS,
    	"TODO: doc"
    },
    {
    	"cell_iterator_mark",
    	dmga2py_cell_iterator_mark,
    	METH_VARARGS,
    	"TODO: doc"
    },
    {
    	"cell_iterator_unmark",
    	dmga2py_cell_iterator_unmark,
    	METH_VARARGS,
    	"TODO: doc"
    },
    {
		"new_alpha_shape",
		dmga2py_new_alpha_shape,
		METH_VARARGS,
		"TODO: doc"
    },
    {
		"alpha_shape_get_shape_cell",
		dmga2py_alpha_shape_get_shape_cell,
		METH_VARARGS,
		"TODO: doc"
	},
	{
		"alpha_shape_cell_as_tuple",
		dmga2py_alpha_shape_cell_as_tuple,
		METH_VARARGS,
		"TODO: doc"
	},
	{
		"alpha_shape_get_complex",
		dmga2py_alpha_shape_get_complex,
		METH_VARARGS,
		"TODO: doc"
	},
    {
		"alpha_shape_get_max_alpha_threshold",
		dmga2py_alpha_shape_get_max_alpha_threshold,
		METH_VARARGS,
		"TODO: doc"
    },
    {
		"alpha_shape_cell_get_max_alpha_threshold",
		dmga2py_alpha_shape_cell_get_max_alpha_threshold,
		METH_VARARGS,
		"TODO: doc"
    },
    {
		"alpha_shape_complex_get_max_alpha_threshold",
		dmga2py_alpha_shape_complex_get_max_alpha_threshold,
		METH_VARARGS,
		"TODO: doc"
    },
    {
		"alpha_shape_complex_print",
		dmga2py_alpha_shape_complex_print,
		METH_VARARGS,
		"TODO: doc"
    },
    {
		"alpha_shape_complex_get_simplex",
		dmga2py_alpha_shape_complex_get_simplex,
		METH_VARARGS,
		"TODO: doc"
	},
	{
		"alpha_shape_complex_upper_bound",
		dmga2py_alpha_shape_complex_upper_bound,
		METH_VARARGS,
		"TODO: doc"
	},
	{
		"alpha_shape_complex_upper_bound",
		dmga2py_alpha_shape_complex_upper_bound,
		METH_VARARGS,
		"TODO: doc"
	},
	{
		"new_sasa_shape",
		dmga2py_new_sasa_shape,
		METH_VARARGS,
		"TODO: doc"
	},
	{
		"sasa_shape_get_shape_cell",
		dmga2py_sasa_shape_get_shape_cell,
		METH_VARARGS,
		"TODO: doc"
	},
	{
		"sasa_shape_cell_as_tuple",
		dmga2py_sasa_shape_cell_as_tuple,
		METH_VARARGS,
		"TODO: doc"
	},
	{
		"sasa_shape_get_area",
		dmga2py_sasa_shape_get_area,
		METH_VARARGS,
		"TODO: doc"
	},
	{
		"sasa_shape_cell_get_area",
		dmga2py_sasa_shape_cell_get_area,
		METH_VARARGS,
		"TODO: doc"
	},
	{
		"sasa_shape_cell_get_excluded_area",
		dmga2py_sasa_shape_cell_get_excluded_area,
		METH_VARARGS,
		"TODO: doc"
	},
	{
		"sasa_shape_cell_get_volume",
		dmga2py_sasa_shape_cell_get_volume,
		METH_VARARGS,
		"TODO: doc"
	},
	{
		"sasa_shape_cell_get_excluded_volume",
		dmga2py_sasa_shape_cell_get_excluded_volume,
		METH_VARARGS,
		"TODO: doc"
	},
	{
		"sasa_contours_new_border_iterator",
		dmga2py_sasa_contours_new_border_iterator,
		METH_VARARGS,
		"TODO: doc"
	},
	{
		"sasa_contours_new_polygon_iterator",
		dmga2py_sasa_contours_new_polygon_iterator,
		METH_VARARGS,
		"TODO: doc"
	},
	{
		"sasa_contours_iterator_next_vertex",
		dmga2py_sasa_contours_iterator_next_vertex,
		METH_VARARGS,
		"TODO: doc"
	},
	{
		"sasa_contours_iterator_next_center",
		dmga2py_sasa_contours_iterator_next_center,
		METH_VARARGS,
		"TODO: doc"
	},
	{
		"sasa_contours_iterator_is_marked",
		dmga2py_sasa_contours_iterator_is_marked,
		METH_VARARGS,
		"TODO: doc"
	},
	{
		"sasa_contours_iterator_is_marked_current",
		dmga2py_sasa_contours_iterator_is_marked_current,
		METH_VARARGS,
		"TODO: doc"
	},
	{
		"sasa_contours_iterator_mark",
		dmga2py_sasa_contours_iterator_mark,
		METH_VARARGS,
		"TODO: doc"
	},
	{
		"sasa_contours_iterator_mark_current",
		dmga2py_sasa_contours_iterator_mark_current,
		METH_VARARGS,
		"TODO: doc"
	},
	{
		"sasa_contours_iterator_unmark",
		dmga2py_sasa_contours_iterator_unmark,
		METH_VARARGS,
		"TODO: doc"
	},
	{
		"sasa_contours_iterator_unmark_current",
		dmga2py_sasa_contours_iterator_unmark_current,
		METH_VARARGS,
		"TODO: doc"
	},
	{
		"sasa_contours_iterator_current_vertex_coords",
		dmga2py_sasa_contours_iterator_current_vertex_coords,
		METH_VARARGS,
		"TODO: doc"
	},
	{
		"sasa_contours_iterator_next_vertex_coords",
		dmga2py_sasa_contours_iterator_next_vertex_coords,
		METH_VARARGS,
		"TODO: doc"
	},
	{
		"sasa_contours_iterator_center_on_plane_coords",
		dmga2py_sasa_contours_iterator_center_on_plane_coords,
		METH_VARARGS,
		"TODO: doc"
	},
	{
		"sasa_contours_iterator_center_on_sphere_coords",
		dmga2py_sasa_contours_iterator_center_on_sphere_coords,
		METH_VARARGS,
		"TODO: doc"
	},
	{
		"sasa_contours_iterator_ready",
		dmga2py_sasa_contours_iterator_ready,
		METH_VARARGS,
		"TODO: doc"
	},
	{
		"sasa_contours_iterator_reset",
		dmga2py_sasa_contours_iterator_reset,
		METH_VARARGS,
		"TODO: doc"
	},
	{
		"sasa_contours_iterator_forward",
		dmga2py_sasa_contours_iterator_forward,
		METH_VARARGS,
		"TODO: doc"
	},
	{
		"sasa_contours_iterator_jump",
		dmga2py_sasa_contours_iterator_jump,
		METH_VARARGS,
		"TODO: doc"
	},
	{
		"sasa_contours_iterator_iterate",
		dmga2py_sasa_contours_iterator_iterate,
		METH_VARARGS,
		"TODO: doc"
	},
	{
		"sasa_contours_iterator_as_tuple",
		dmga2py_sasa_contours_iterator_as_tuple,
		METH_VARARGS,
		"TODO: doc"
	},
	{
		"sasa_contours_iterator_as_tuple_coords",
		dmga2py_sasa_contours_iterator_as_tuple_coords,
		METH_VARARGS,
		"TODO: doc"
	},

    {
		NULL,
		NULL,
		0,
		NULL
    }
};



static struct PyModuleDef dmga2py = {
  PyModuleDef_HEAD_INIT,
  "dmga2py",
  NULL,
  -1,
  Dmga2pyMethods,
  NULL,
  NULL,
  NULL,
  NULL
};

PyMODINIT_FUNC PyInit_dmga2py(void) {
  return PyModule_Create(&dmga2py);
}

int
main(int argc, char *argv[])
{
    wchar_t *program = Py_DecodeLocale(argv[0], NULL);
    if (program == NULL) {
        fprintf(stderr, "Fatal error: cannot decode argv[0]\n");
        exit(1);
    }
    /* Add a built-in module, before Py_Initialize */
    PyImport_AppendInittab("dmga2py", PyInit_dmga2py);

    /* Pass argv[0] to the Python interpreter */
    Py_SetProgramName(program);

    /* Initialize the Python interpreter.  Required. */
    Py_Initialize();

    /* Add a static module */
    PyInit_dmga2py();

    PyMem_RawFree(program);
	return 0;
}



}
