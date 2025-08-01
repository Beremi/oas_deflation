import sys
sys.path.append('../')

from pydmga.geometry import CastCylinderGeometry
from pydmga.container import Container
from pydmga.diagram import Diagram
from pydmga.draw import assets
from pydmga.draw import render
from pydmga.draw.render import Asset
from pydmga.draw.render import Color
from random import random
from random import seed
from math import pi

(R, H) = (5.0, 10.0)
count = 64
data = [(i, 2.0 * R * random() - R, 2.0 * R * random() -
         R, H * random(), 0.1) for i in range(count)]

geometry = CastCylinderGeometry(R, H, 128)
container = Container(geometry)
container.add(data)

diagram = Diagram(container)

display = render.Render(Color.WHITE_OPAQUE)


def is_divider(p1, p2):
    return ((p1[0] - p2[0])**2 + (p1[1] - p2[1])**2) < 0.0001


colors = [(1, 0, 0, 0.3), (0, 1, 0, 0.3), (0, 0, 1, 0.3),
          (1, 1, 0, 0.3), (0, 1, 1, 0.3), (1, 0, 1, 0.3)]
all_area = 0
for i, cell in enumerate(diagram):
    (id, x, y, z, r) = container.get(i)
    display.add(render.Asset(assets.sphere, (x, y, z), r, Color.BLUE_QUARTER))
    for j, side in enumerate(cell.sides):
        if (geometry.on_boundary(side.neighbour)):
            polygon = side.as_coords()
            # display.add(render.Outline(polygon))
            for k in range(len(polygon)):  # not quite that I wanted, but...
                (p, q) = (polygon[k], polygon[(k+1) % len(polygon)])
                if not(is_divider(p, q)):
                    display.add(render.Lines([p, q], 1.0))
            display.add(render.ConvexPolygons(
                polygon, colors[i % len(colors)]))
            display.add(render.ConvexPolygons(
                reversed(polygon), colors[i % len(colors)]))
            all_area += side.area()

print("All area is", geometry.adjust_area(
    all_area), "should be", (2.0 * pi * R * H))

display.run()
