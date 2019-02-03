import sys
sys.path.append('../')

from pydmga import draw
from pydmga.draw import assets
from pydmga.draw import geometry
from math import sqrt
import random
#balls = [geometry.Ball(0, 5.0, 5.0, 5.0, 2.0)]
balls = []
# balls = [geometry.Ball(0, 5.0, 5.0, 5.0, 2.0),
#         geometry.Ball(1, 5.0-1.0, 5.0, 5.0, 1.5),
#         geometry.Ball(2, 5.0+1.0, 5.0, 5.0, 1.5),
#         geometry.Ball(3, 5.0+0.5, 5.0+0.86, 5.0, 1.5),
#         geometry.Ball(4, 5.0+0.5, 5.0-0.86, 5.0, 1.5),
#         geometry.Ball(5, 5.0-0.5, 5.0+0.86, 5.0, 1.5),
#         geometry.Ball(6, 5.0-0.5, 5.0-0.86, 5.0, 1.5),
#         geometry.Ball(7, 5.0, 5.0, 5.0+1.5, 0.5),
#         geometry.Ball(8, 5.0, 5.0, 5.0-1.5, 1.8)]

r = draw.VoroRender(geometry.ArrayInput(balls, []))
r.grow_particles = True

colors = [(1.0, 0.0, 0.0, 1.0), (0.0, 1.0, 0.0, 1.0), (0.0, 0.0, 1.0, 1.0),
          (1.0, 1.0, 0.0, 1.0), (0.0, 1.0, 1.0, 1.0), (1.0, 0.0, 1.0, 1.0)]
colors_t = [(1.0, 0.0, 0.0, 0.01), (0.0, 1.0, 0.0, 0.01), (0.0, 0.0, 1.0, 0.01),
            (1.0, 1.0, 0.0, 0.01), (0.0, 1.0, 1.0, 0.01), (1.0, 0.0, 1.0, 0.01)]
file = open('../../bin/test/sasa_test_1.txt', 'r')
exec(file)
voronoi = sasa['voronoi']
# for i, side in enumerate([voronoi[1], voronoi[2], voronoi[4], voronoi[5], voronoi[6], voronoi[7]]):
for i, side in enumerate(voronoi):
    r.add(draw.ConvexPolygons(reversed(side), colors_t[i % 6]))
    r.add(draw.ConvexPolygons(side, colors_t[i % 6]))
centers = sasa['centers']
# for i, side_center in enumerate([centers[1], centers[2], centers[4], centers[5], centers[6], centers[7]]):
# for i, side_center in enumerate(centers):
#     pc = draw.PointCloud().add_vertex(side_center, 3.0, colors[i % 6])
#     r.add(pc)

# r.add(draw.PointCloud().add_vertices(sasa['centers'], 3.0, (1.0, 0.0,0.0,1.0)))
# r.add(draw.PointCloud().add_vertices(sasa['vertices'], 3.0, (1.0, 1.0,0.0,1.0)))
border = sasa['borders']
# for i, border in enumerate([border[0], border[2], border[4]]):
for i, b in enumerate(border):
    border_arcs = draw.SphereArcs((5.0, 5.0, 5.0), 2.0, 2, colors[i % 6])
    border_arcs.add_arc(b[1], b[0], b[2])
    r.add(draw.PointCloud().add_vertices([b[1]], 3.0, colors[i % 6]))
    r.add(border_arcs)

polygons = sasa['polygons']
ga = draw.GreatArcs((5.0, 5.0, 5.0), 1.0, 1.0, (0.7, 1.0, 0.2, 1.0))
for poly in polygons:
    ga.add_arcs(poly, False)
r.add(ga)


r.add(draw.Lines([(5.3165, 4.5, 6.28446), (4.5, 4.5, 5.70711),
                  (4.5, 3.6835, 5.12976)], 2.0, (1.0, 1.0, 0.0, 0.5)))
r.add(draw.Lines([(4, 5.8165, 4.42265), (4, 5, 5),
                  (4, 4.1835, 4.42265)], 2.0, (1.0, 1.0, 0.0, 0.5)))

# for i, poly in enumerate([polygons[0]]):
#     poly_arcs = draw.SphereArcs((5.0, 5.0, 5.0), 2.0, 2, colors[i % 6])
#     poly_arcs.add_cutting_plane_arc((5.0, 5.0, 5.0), poly[0], poly[1])
#     poly_arcs.add_cutting_plane_arc((5.0, 5.0, 5.0), poly[1], poly[2])
#     #r.add(draw.PointCloud().add_vertices([b[1]], 3.0, colors[i % 6]))
#     r.add(poly_arcs)


def translate(myNodePath, x, y, z):
    for child in myNodePath.getChildren():
        child.setPos(x, y, z)
        translate(child, x, y, z)


translate(r.render, -5.0, -5.0, -5.0)
r.add(draw.Asset(assets.sphere, (0.0, 0.0, 0.0), 1.0, (1.0, 1.0, 1.0, 0.3)))

# r.add(draw.ConvexPolygons([(-1.0, -1.0, 0.0),(-1.0, 1.0, 0.0),(1.0, 1.0, 0.5)], (1.0, 0.0, 0.0, 0.5)))
# r.add(draw.ConvexPolygons([(-1.0, -1.0, 0.0),(1.0, 1.0, 0.5),(-1.0, 1.0, 0.0)], (1.0, 0.0, 0.0, 0.5)))
# r.add(draw.Outline([(-1.0, -1.0, 0.0),(1.0, 1.0, 0.5),(-1.0, 1.0, 0.0)], 2.0, (1.0, 1.0, 0.0, 0.5)))
# v = [(random.random() * 10.0 - 5.0, random.random() * 10.0 - 5.0, random.random() * 10.0 - 5.0) for i in range(100)]
# r.add(draw.PointCloud().add_vertices(v, 2.0, (0.0,1.0,1.0,1.0)))
# r.add(draw.Asset(assets.sphere, (0.0, 0.0, 0.0), 5.0, (1.0, 1.0, 1.0, 0.5)))
# ga = draw.GreatArcs((0.0, 0.0, 0.0), 5.0).add_arcs([(5.0,0.0,0.0),(0.0,5.0,0.0),(0.0,0.0,5.0)])
# r.add(ga)
# r.add(draw.PointCloud().add_vertices(ga.vertices, 12.0, (0.0, 1.0,0.0,1.0)))
# #r.add(draw.SphereArcs((0, 0, 0), 5, 1, (1,0,0,1)).add_arc((-4, 0.0, 0.0),(-4,3,0),(-4,0,3)))
# #r.add(draw.SphereArcs((0, 0, 0), 5, 2, (0,1,0,1)).add_arc((-4, 0.0, 0.0),(-4,0,3),(-4,1,sqrt(8))))
# r.add(draw.SphereArcs((0, 0, 0), 5, 2, (0,1,0,1)).add_cutting_plane_arc((-1, 0.0, 0.0), 4.0, (-4,0,3),(-4,3,0)))
# r.add(draw.SphereArcs((0, 0, 0), 5, 2, (0,0,1,1)).add_cutting_plane_arc((-1, 0.0, 0.0), 4.0, (-4,3,0),(-4,0,3)))
# r.add(draw.SphereArcs((0, 0, 0), 5, 2, (0,1,0,1)).add_cutting_plane_arc((0, 1.0, 0.0), 4.0, (0,4,3),(0,4,-3)))
# r.add(draw.SphereArcs((0, 0, 0), 5, 2, (0,0,1,1)).add_cutting_plane_arc((0, 1.0, 0.0), 4.0, (0,4,-3),(0,4,3)))
# r.add(draw.SphereArcs((0, 0, 0), 5, 2, (1,0,0,1)).add_cutting_plane_arc((0, 0.0, 1.0), 0.0, (0,-5,0),(-5,0,0)))
# r.add(draw.SphereArcs((0, 0, 0), 5, 2, (0.5,1,0,1)).add_cutting_plane_arc((0, 0.0, 1.0), 0.0, (-5,0,0), (0,-5,0)))
# r.add(draw.PointCloud().add_vertices([(-4,0,3),(-4,3,0),(0,4,-3),(0,4,3)], 8.0, (0.0, 1.0,1.0,1.0)))
# r.add(draw.PointCloud().add_vertices([(0,-5,0),(-5,0,0)], 8.0, (1.0, 0.0,0.0,1.0)))
# r.add(draw.Asset(assets.cube, (-5.0, 5.0, 5.0), 1.0, (0.0, 0.0, 1.0, 0.5)))
r.run()
