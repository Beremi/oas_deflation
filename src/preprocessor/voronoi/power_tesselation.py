import logging
import datetime
import time
from collections import defaultdict
from pydmga.geometry import BoxGeometry
from pydmga.container import Container
from pydmga.diagram import Diagram
import numpy as np
from scipy.spatial.distance import cdist

class TimeFilter(logging.Filter):

    def filter(self, record):
        try:
          last = self.last
        except AttributeError:
          last = record.relativeCreated

        delta = datetime.datetime.fromtimestamp(record.relativeCreated/1000.0) - datetime.datetime.fromtimestamp(last/1000.0)

        record.relative = '{0:.6f}'.format(delta.seconds + delta.microseconds/1000000.0)

        self.last = record.relativeCreated
        return True

logging.basicConfig(level=logging.INFO)
fmt = logging.Formatter(fmt="%(asctime)s %(levelname)s (%(relative)ss) - %(message)s")
log = logging.getLogger()
[hndl.addFilter(TimeFilter()) for hndl in log.handlers]
[hndl.setFormatter(fmt) for hndl in log.handlers]


class PowerTesselation(object):
    """
    PowerTesselation(points, weights=None, maxLim=None)
    Power diagrams in N dimensions using pydmga.
    Parameters
    ----------
    points : ndarray of floats, shape (npoints, ndim)
        Coordinates of points to construct a convex hull from
    weights : ndarray of floats, shape (npoints), optional
        Weights of points to construct power diagram. Default: None
    maxLim: ndarray of floats, shape (npoints), optional
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

    def __init__(self, points, weights=None, limits='unit'):
        self._points = np.asarray(points, dtype=np.float64)
        self._npoints, self._ndim = points.shape
        if weights is None:
            self._weights = np.zeros(self._npoints)
        else:
            self._weights = np.asarray(weights, dtype=np.float64)
        if self._ndim not in [2, 3]:
            raise ValueError('Shape of array has to be (npoints, ndim) for ndim = [2, 3] but ndim = {}.'.format(self._ndim))
        if not self._points.shape[0] == self._weights.shape[0]:
            raise ValueError('Number of points is not eqaul to number of weights {} != {}.'.format(points.shape[0], weights.shape[0]))
        # convert points+weights to array for pydgma
        self._points_pydgma = self._get_points_pydgma()
        if limits == 'auto':
            if self._ndim == 2:
                xmin, ymin = points.min(axis=0) - 1
                xmax, ymax = points.max(axis=0) + 1
                zmin, zmax = 0, 1
            else:
                xmin, ymin, zmin = points.min(axis=0) - 1 # 0, 0, 0
                xmax, ymax, zmax = points.max(axis=0) + 1 # 1, 1, 1
        elif limits == 'unit':
            xmin, ymin, zmin = 0, 0, 0
            xmax, ymax, zmax = 1, 1, 1
        else:
            if self._ndim == 2:
                xmin, ymin, xmax, ymax = limits
                zmin, zmax = 0, 1
            else:
                xmin, ymin, zmin, xmax, ymax, zmax = limits

        point_distances = cdist(points, points)
        point_distances = point_distances[np.triu_indices(point_distances.shape[0], k=1)]
        if np.nanmin(point_distances) == 0:
            raise ValueError('There are identical points (distance is zero).')

        start = time.time()
        #dx = xmax - xmin
        #dy = ymax - ymin
        #dz = zmax - zmin
        self.geometry = BoxGeometry(xmin, ymin, zmin, xmax, ymax, zmax, False, False, False)
        self.container = Container(self.geometry)
        self.container.add(self._points_pydgma)
        self.diagram = Diagram(self.container, False)
        logging.info('Time to generate pydgma diagram = {} s'.format(time.time() - start))

        start = time.time()
        if self._ndim == 2:
            self._generate_voronoi_parts_2d()
        if self._ndim == 3:
            self._generate_voronoi_parts_3d()
        logging.info('Time to generate voronoi parts = {} s'.format(time.time() - start))
        #print('#'*50)
        #print('points:', self.points)
        #print('vertices:', self.vertices)
        #print('ridge_points:', self.ridge_points)
        #print('ridge_vertices:', self.ridge_vertices)
        #print('regions:', self.regions)
        #print('point_region:', self.point_region)

    def _get_points_pydgma(self):
        points_pydgma = []
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
        return self._vertices

    @property
    def ridge_points(self):
        return self._ridge_points

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
        connection_set = set()
        ridge_vertices = []
        vertices = []
        vertex_start_num = 0
        regions = [[], ]
        point_region = []
        start = time.time()

        for i, cell in enumerate(self.diagram):
            region = []
            point_region.append(i + 1)
            (id_, x, y, z, r) = self.container.get(i)
            self._total_volume += cell.volume()
            for side in cell.sides:
                point_con = tuple(sorted((i, side.neighbour)))
                if point_con not in connection_set:
                    connection_list.append(point_con)
                    connection_set.add(point_con)
                    ridge_vertices.append(np.array(side.as_list()) + vertex_start_num)

            for vertex in cell.vertices:
                id_, x, y, z = vertex.as_tuple()
                vertices.append((vertex_start_num + id_, x, y, z))
                region.append(vertex_start_num + id_)
            regions.append(region)

            vertex_start_num += cell.vertices.size()
        logging.debug('Time to generate voronoi parts - part = {} s'.format(time.time() - start))
        self._vertices = vertices

        self._ridge_vertices = ridge_vertices
        self._regions = regions
        #self._ridge_points = connection_list
        self._point_region = point_region
        start = time.time()
        self._merge_duplicate_vertices()
        logging.debug('Time to generate voronoi parts - merge = {} s'.format(time.time() - start))

        ridge_points = np.array(connection_list)
        neg = ridge_points < 0
        ridge_points[neg] = np.max(ridge_points) - ridge_points[neg]
        self._ridge_points = ridge_points
        self._vertices = np.array(self._vertices)[:, 1:(self._ndim + 1)]

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
                # TODO: why not - if not np.allclose(np.array(side.as_coords())[:, -1], 0):
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
        #self._ridge_points = connection_list
        self._point_region = point_region
        self._merge_duplicate_vertices()

        ridge_points = np.array(connection_list)
        neg = ridge_points < 0
        ridge_points[neg] = np.max(ridge_points) - ridge_points[neg]
        self._ridge_points = ridge_points
        self._vertices = np.array(self._vertices)[:, 1:(self._ndim + 1)]

    def _merge_duplicate_vertices(self):
        def replace(val, dict_):
            return dict_[val] if val in dict_ else -1
        replace = np.vectorize(replace)
        unique_dict = defaultdict(list)
        orig_dict = {}
        for id_, x, y, z in self._vertices:
            if (self._ndim == 2) and (z != 0):
                continue
            key = tuple(round(i, 10) for i in (x, y, z))
            unique_dict[key].append(id_)
            orig_dict[key] = (x, y, z)
        replace_dict = {}
        vertices_unique = []
        for idx, (key, val) in enumerate(sorted(unique_dict.items())):
            for v in val:
                replace_dict[v] = idx
            vertices_unique.append([idx] + list(orig_dict[key]))

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
    #points = np.array([[0.25, 0.5], [0.75, 0.5]])
    #radii = np.array([.1, .1])
    #vor = PowerTesselation(points, weights=radii)
    import doctest
    doctest.testmod()
