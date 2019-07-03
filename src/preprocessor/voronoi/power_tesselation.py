import logging
import time
from collections import defaultdict
from pydmga.geometry import OrthogonalGeometry
from pydmga.container import Container
from pydmga.diagram import Diagram
import numpy as np

logging.basicConfig(level=logging.INFO,
                    format='%(asctime)s %(levelname)s - %(message)s')


class PowerTesselation(object):
    """
    PowerTesselation(points, weights=None, limits=None)
    Power diagrams in N dimensions using pydmga.
    Parameters
    ----------
    points : ndarray of floats, shape (npoints, ndim)
        Coordinates of points to construct a convex hull from
    weights : ndarray of floats, shape (npoints), optional
        Weights of points to construct power diagram. Default: None
    limits: ndarray of floats, shape (npoints), optional
        Dimensions of box - (dx, dy [, dz])
    Attributes
    ----------
    points : ndarray of double, shape (npoints, ndim)
        Coordinates of input points.
    vertices : ndarray of double, shape (nvertices, ndim)
        Coordinates of the Voronoi vertices.
    ridge_points : ndarray of ints, shape ``(nridges, 2)``
        Indices of the points between which each Voronoi ridge lies.
    ridge_vertices : list of list of ints, shape ``(nridges, *)``
        Indices of the Voronoi vertices forming each Voronoi ridge.
    regions : list of list of ints, shape ``(nregions, *)``
        Indices of the Voronoi vertices forming each Voronoi region.
        -1 indicates vertex outside the Voronoi diagram.
    point_region : ndarray of ints, shape (npoints)
        Index of the Voronoi region for each input point.
    Raises
    ------
    ValueError
        Raised if an incompatible array is given as input.
    Notes
    -----
    The power diagram is computed using the
    `DMG-alpha <http://37.187.139.172/~robson/dmga_doc/>`__.
    Examples
    --------
    Voronoi diagram for a set of point:
    >>> points = np.array([[0.25, 0.5], [0.75, 0.5]])
    >>> radii = np.array([.1, .1])
    >>> from power_tesselation import PowerTesselation
    >>> vor = PowerTesselation(points, weights=radii)

    Plot it:
    >>> import matplotlib.pyplot as plt
    >>> from scipy.spatial import Voronoi, voronoi_plot_2d
    >>> fig = voronoi_plot_2d(vor)
    >>> ax = fig.gca()
    >>> _ = [ax.text(v[0], v[1], '%d' % vi) for vi, v in enumerate(vor.vertices)]
    >>> _ = ax.set_xlim(0, 1)
    >>> _ = ax.set_ylim(0, 1)
    >>> plt.show()

    The Voronoi vertices:
    >>> vor.vertices
    array([[ 0. ,  0. ],
           [ 0. ,  1. ],
           [ 0.5,  0. ],
           [ 0.5,  1. ],
           [ 1. ,  0. ],
           [ 1. ,  1. ]])

    There is a single finite Voronoi region, and four finite Voronoi
    ridges:
    >>> vor.regions
    [[], [0, 2, 3, 1], [2, 4, 5, 3]]
    >>> vor.ridge_vertices
    [[0, 2], [0, 1], [2, 3], [1, 3], [2, 4], [4, 5], [3, 5]]

    The ridges are perpendicular between lines drawn between the following
    input points:
    >>> vor.ridge_points
    array([[-3,  0],
           [-1,  0],
           [ 0,  1],
           [-4,  0],
           [-3,  1],
           [-2,  1],
           [-4,  1]])
    """

    def __init__(self, points, weights=None, limits=None):
        self._points = np.ascontiguousarray(points, dtype=np.double)
        self._npoints, self._ndim = points.shape
        if self._ndim not in [2, 3]:
            raise ValueError('Shape of array has to be (npoints, ndim) for ndim = [2, 3] but ndim = {}.'.format(self._ndim))
        self._weights = weights
        # convert points+weights to array for pydgma
        self._points_pydgma = self._get_points_pydgma()
        if limits is None:
            dx, dy, dz = 1, 1, 1
        else:
            dx = limits[0]
            dy = limits[1]
            if self._ndim == 2:
                dz = 1
            else:
                dz = limits[2]

        start = time.time()
        self.geometry = OrthogonalGeometry(dx, dy, dz, False, False, False)
        self.container = Container(self.geometry)
        self.container.add(self._points_pydgma)
        self.diagram = Diagram(self.container, False)
        logging.info('Time to generate diagram = {} s'.format(time.time() - start))

        start = time.time()
        if self._ndim == 2:
            self._generate_voronoi_parts_2d()
        if self._ndim == 3:
            self._generate_voronoi_parts_3d()
        logging.info('Time to generate voronoi parts = {} s'.format(time.time() - start))

    def _get_points_pydgma(self):
        points_pydgma = []
        if self._weights is None:
            self._weights = np.zeros(self._npoints)
        for idx, row in enumerate(self._points):
            row_tmp = [idx]
            row_tmp.extend(row.tolist())
            if self._ndim == 2:
                row_tmp.append(0)
            row_tmp.append(self._weights[idx])
            points_pydgma.append(tuple(row_tmp))
        return points_pydgma

    @property
    def points(self):
        return self._points

    @property
    def weights(self):
        return self._weights

    @property
    def vertices(self):
        return np.array(self._vertices)[:, 1:(self._ndim + 1)]

    @property
    def ridge_points(self):
        return np.array(self._ridge_points)

    @property
    def ridge_vertices(self):
        return self._ridge_vertices

    @property
    def regions(self):
        return self._regions

    @property
    def point_region(self):
        return np.array(self._point_region)

    def _generate_voronoi_parts_3d(self):
        self._total_volume = 0.0
        connection_list = []
        ridge_vertices = []
        vertices = []
        vertex_start_num = 0
        regions = [[], ]
        point_region = []
        for i, cell in enumerate(self.diagram):
            region = []
            point_region.append(i + 1)
            (id_, x, y, z, r) = self.container.get(i)
            self._total_volume += cell.volume()
            for side in cell.sides:
                point_con = tuple(sorted((i, side.neighbour)))
                if point_con not in connection_list:
                    connection_list.append(point_con)
                    ridge_vertices.append(np.array(side.as_list()) + vertex_start_num)

            for vertex in cell.vertices:
                id_, x, y, z = vertex.as_tuple()
                vertices.append((vertex_start_num + id_, x, y, z))
                region.append(vertex_start_num + id_)
            regions.append(region)

            vertex_start_num += cell.vertices.size()
        self._ridge_vertices = ridge_vertices
        self._vertices = vertices
        self._regions = regions
        self._ridge_points = connection_list
        self._point_region = point_region
        self._merge_duplicate_vertices()

    def _generate_voronoi_parts_2d(self):
        self._total_volume = 0.0
        connection_list = []
        ridge_vertices = []
        vertices = []
        vertex_start_num = 0
        regions = [[], ]
        point_region = []
        for i, cell in enumerate(self.diagram):
            point_region.append(i+1)
            (id_, x, y, z, r) = self.container.get(i)
            self._total_volume += cell.volume()
            for side in cell.sides:
                if not np.all(np.array(side.as_coords())[:, -1] == 0):
                    continue
                else:
                    verts = np.array(side.as_list()) + vertex_start_num
                    for id_, (x, y, z) in zip(verts, side.as_coords()):
                        vertices.append((id_, x, y, z))
                    regions.append(verts)

            for side in cell.sides:
                if (np.allclose(np.array(side.as_coords())[:, 2], 0)
                        or np.allclose(np.array(side.as_coords())[:, 2], 1)):
                    continue

                point_con = tuple(sorted((i, side.neighbour)))
                if point_con not in connection_list:
                    connection_list.append(point_con)
                    ridge_vertices.append(np.array(side.as_list()) + vertex_start_num)

            vertex_start_num += cell.vertices.size()
        self._ridge_vertices = ridge_vertices
        self._vertices = vertices
        self._regions = regions
        self._ridge_points = connection_list
        self._point_region = point_region
        self._merge_duplicate_vertices()

    def _merge_duplicate_vertices(self):
        def replace(val, dict_):
            return dict_[val] if val in dict_ else -1
        replace = np.vectorize(replace)
        unique_dict = defaultdict(list)
        for id_, x, y, z in self._vertices:
            if (self._ndim == 2) and (z != 0):
                continue
            key = tuple(round(i, 13) for i in (x, y, z))
            unique_dict[key].append(id_)
        replace_dict = {}
        vertices_unique = []
        for idx, (key, val) in enumerate(sorted(unique_dict.items())):
            for v in val:
                replace_dict[v] = idx
            vertices_unique.append([idx] + list(key))

        self._vertices = vertices_unique

        for ri, r in enumerate(self._regions):
            if ri == 0:
                continue
            r_replaced = replace(r, replace_dict)
            self._regions[ri] = list(r_replaced)

        for rvi, rv in enumerate(self._ridge_vertices):
            rv_replaced = replace(rv, replace_dict)
            if self._ndim == 2:
                rv_replaced = [i for i in rv_replaced if i != -1]
            self._ridge_vertices[rvi] = rv_replaced


if __name__ == '__main__':
    import doctest
    doctest.testmod()
