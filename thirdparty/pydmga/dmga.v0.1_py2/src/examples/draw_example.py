import sys
sys.path.append('../')

from pydmga import draw
from pydmga.draw import render
from pydmga.draw import assets
from math import sqrt
import random

r = draw.Render(render.Color.WHITE_OPAQUE)
r.add(draw.ConvexPolygons([(-1.0, -1.0, 0.0),(-1.0, 1.0, 0.0),(1.0, 1.0, 0.5)], (1.0, 0.0, 0.0, 0.5)))
r.add(draw.ConvexPolygons([(-1.0, -1.0, 0.0),(1.0, 1.0, 0.5),(-1.0, 1.0, 0.0)], (1.0, 0.0, 0.0, 0.5)))
r.add(draw.Outline([(-1.0, -1.0, 0.0),(1.0, 1.0, 0.5),(-1.0, 1.0, 0.0)], 2.0, (1.0, 1.0, 0.0, 0.5)))
v = [(random.random() * 10.0 - 5.0, random.random() * 10.0 - 5.0, random.random() * 10.0 - 5.0) for i in range(100)]
for i in range(len(v)):
    d = sqrt((v[i][0])**2 + (v[i][1])**2 + (v[i][2])**2);
    if (d > 5.0):
        v[i] = (5.0 * v[i][0] / d, 5.0 * v[i][1] / d, 5.0 * v[i][2] / d)
r.add(draw.PointCloud().add_vertices(v, 5.0, (0.5,0.5,0.5,1.0)))
r.add(draw.Asset(assets.sphere, (0.0, 0.0, 0.0), 5.0, render.Color.BLUE_QUARTER))
ga = draw.GreatArcs((0.0, 0.0, 0.0), 5.0).add_arcs([(5.0,0.0,0.0),(0.0,5.0,0.0),(0.0,0.0,5.0)])
r.add(ga)
r.add(draw.PointCloud().add_vertices(ga.vertices, 12.0, (0.0, 1.0,0.0,1.0)))
r.add(draw.SphereArcs((0, 0, 0), 5, 2, (0,1,0,1)).add_cutting_plane_arc((-1, 0.0, 0.0), 4.0, (-4,0,3),(-4,3,0)))
r.add(draw.SphereArcs((0, 0, 0), 5, 2, (0,0,1,1)).add_cutting_plane_arc((-1, 0.0, 0.0), 4.0, (-4,3,0),(-4,0,3)))
r.add(draw.SphereArcs((0, 0, 0), 5, 2, (0,1,0,1)).add_cutting_plane_arc((0, 1.0, 0.0), 4.0, (0,4,3),(0,4,-3)))
r.add(draw.SphereArcs((0, 0, 0), 5, 2, (0,0,1,1)).add_cutting_plane_arc((0, 1.0, 0.0), 4.0, (0,4,-3),(0,4,3)))
r.add(draw.SphereArcs((0, 0, 0), 5, 2, (1,0,0,1)).add_cutting_plane_arc((0, 0.0, 1.0), 0.0, (0,-5,0),(-5,0,0)))
r.add(draw.SphereArcs((0, 0, 0), 5, 2, (0.5,1,0,1)).add_cutting_plane_arc((0, 0.0, 1.0), 0.0, (-5,0,0), (0,-5,0)))
r.add(draw.PointCloud().add_vertices([(-4,0,3),(-4,3,0),(0,4,-3),(0,4,3)], 8.0, (0.0, 1.0,1.0,1.0)))
r.add(draw.PointCloud().add_vertices([(0,-5,0),(-5,0,0)], 8.0, (1.0, 0.0,0.0,1.0)))
#r.add(draw.Asset(assets.cube, (-5.0, 5.0, 5.0), 1.0, render.Color.GREEN_OPAQUE))
r.run()
