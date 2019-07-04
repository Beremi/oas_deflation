import sys
import pathlib
import time
import logging
from scipy.spatial import Voronoi, voronoi_plot_2d
import matplotlib.pyplot as plt
import numpy as np
print(pathlib.Path(__file__))
PROJECT_DIR = pathlib.Path(__file__).resolve().parents[1]
print(PROJECT_DIR)

sys.path.append(str(PROJECT_DIR / 'src' / 'preprocessor' / 'voronoi'))

from point_generators import generateTesC
from voronoi import mirror_dataBeam
from voronoi_viewer import voronoi_plot_3d_vtk
from power_tesselation import PowerTesselation

logging.basicConfig(level=logging.INFO,
                    format='%(asctime)s %(levelname)s - %(message)s')

def test_power_voronoi_2d():
    points, radii = generateTesC(10, 10, seed=1234)
    points /= 10
    radii /= 10

    start = time.time()
    pt = PowerTesselation(points, weights=radii)
    print('time =', time.time() - start)

    fig, ax = plt.subplots()
    voronoi_plot_2d(pt, ax=ax)
    ax.scatter(pt.vertices[:, 0], pt.vertices[:, 1], color='r', zorder=100)
    for (x, y), r in zip(points, radii):
        circle = plt.Circle((x, y), r, color='r', fill=False)
        ax.add_artist(circle)
    ax.set_xlim(0, 1)
    ax.set_ylim(0, 1)
    ax.set_aspect('equal')
    plt.show()


def test_voronoi_2d():
    points, radii = generateTesC(10, 10, seed=1234)
    points /= 10
    radii /= 10
    start = time.time()
    vor = Voronoi(mirror_dataBeam(points, points.shape[1], points.shape[1]*[1]))
    print('time =', time.time() - start)
    #print('#'*50)
    #print('points:', vor.points)
    #print('vertices:', vor.vertices)
    #print('ridge_points:', vor.ridge_points)
    #print('ridge_vertices:', vor.ridge_vertices)
    #print('regions:', vor.regions)
    #print('point_region:', vor.point_region)
    if True:
        fig, ax = plt.subplots()
        voronoi_plot_2d(vor, ax=ax)
        ax.set_xlim(0, 1)
        ax.set_ylim(0, 1)
        ax.set_aspect('equal')
        plt.show()


def test_power_voronoi_3d():
    points, radii = generateTesC(10, 10, 10, seed=1234)
    points /= 10
    radii /= 10

    start = time.time()
    pt = PowerTesselation(points, weights=radii)
    print('time =', time.time() - start)

    voronoi_plot_3d_vtk(pt)


def test_voronoi_3d():
    points, radii = generateTesC(10, 10, 10, seed=1234)
    points /= 10
    radii /= 10
    start = time.time()
    vor = Voronoi(mirror_dataBeam(points, points.shape[1], points.shape[1]*[1]))
    #vor = PowerTesselation(points)
    print('time =', time.time() - start)
    #print('#'*50)
    #print('points:', vor.points)
    #print('vertices:', vor.vertices)
    #print('ridge_points:', vor.ridge_points)
    #print('ridge_vertices:', vor.ridge_vertices)
    #print('regions:', vor.regions)
    #print('point_region:', vor.point_region)

    voronoi_plot_3d_vtk(vor)

if __name__ == '__main__':
    test_power_voronoi_2d()
    test_voronoi_2d()
    test_power_voronoi_3d()
    test_voronoi_3d()
