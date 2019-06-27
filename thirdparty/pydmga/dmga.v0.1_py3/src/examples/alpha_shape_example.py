import sys
sys.path.append('../')

from pydmga.geometry import OrthogonalGeometry
from pydmga.container import Container
from pydmga.diagram import Diagram
from pydmga.shape import AlphaShape
from pydmga import draw
from pydmga.draw import assets
from pydmga.draw import Asset
from pydmga.draw import render
from pydmga.draw.render import Color
from random import random
from random import seed
from math import sqrt, pi, sin, cos

(box_x, box_y, box_z) = (10.0, 10.0, 10.0)
count = 10
#data = [ (i, box_x * random(), box_y * random(), box_z * random(), random()) for i in range(count)]
data = [(i, 3.0 * sin(pi * 2 * i / count), 3.0 * cos(pi * 2 * i / count),
         5.0 + random() * 2.0 - 1.0, 2.0 * random()) for i in range(count)]
data.append((count, 0.0, 0.0, 9.0, 2.1))
data.append((count, 0.0, 0.0, 1.0, 1.2))

geometry = OrthogonalGeometry(10, 10, 10, True, True, True)
container = Container(geometry)
container.add(data)

diagram = Diagram(container, True)
print("Diagram handle =", diagram._diagram_handle)
shape = AlphaShape(diagram)

display = render.Render(Color.WHITE_OPAQUE)
display.show_shape = True

max_alpha = shape.max_alpha_threshold()
alpha_for_shape = 0.0  # max_alpha / 3.0
print("Max alpha for whole thing is", max_alpha,
      "chosen alpha is", alpha_for_shape)
vertices = draw.PointCloud()
for i, cell_shape in enumerate(shape):
    print("max alpha for cell", i, "is", cell_shape.max_alpha_threshold())
    (id, x, y, z, r) = container.get(i)
    display.add(Asset(assets.sphere, (x, y, z), sqrt(
        r*r + alpha_for_shape), Color.BLUE_QUARTER))
    for side in cell_shape:
        print(side.as_coords())
    for simplex in cell_shape.shape():
        print("simplex:")
        print("\t", simplex.as_tuple())
        print("\t", simplex.as_dict())
        print("\t", simplex.as_dict_short())
        if (simplex.alpha_threshold > alpha_for_shape):
            #color = (0.3, 0.3, 1.0 - simplex.alpha_threshold / max_alpha, 0.1)
            color = Color.GRAY30_QUARTER
            if simplex.dimension == 1:
                part = simplex.diagram_part()
                print("Side is:", part.vertices)
                #display.add(draw.ConvexPolygons(part.as_coords(), color))
                #display.add(draw.ConvexPolygons(reversed(part.as_coords()), color))
            elif simplex.dimension == 2:
                part = simplex.diagram_part()
                print("Edge is:", part.as_pair())
                display.add(draw.Outline(part.as_coords(), 1.0, color))
            elif simplex.dimension == 3:
                part = simplex.diagram_part()
                print("Vertex is:", part.id, ",", part.as_coords())
                vertices.add_vertex(part.as_coords(), 2.0, color)
    for simplex in cell_shape.shape(alpha_for_shape):
        print("simplex:")
        print("\t", simplex.as_tuple())
        print("\t", simplex.as_dict())
        print("\t", simplex.as_dict_short())
        #color = (0.0, 1.0 - simplex.alpha_threshold / max_alpha, 0.0, 0.3)
        color = Color.GREEN_QUARTER
        if simplex.dimension == 1:
            part = simplex.diagram_part()
            display.add(draw.ConvexPolygons(part.as_coords(), color))
            display.add(draw.ConvexPolygons(reversed(part.as_coords()), color))
        elif simplex.dimension == 2:
            display.add(draw.Outline(
                simplex.diagram_part().as_coords(), 5.0, color))
        elif simplex.dimension == 3:
            vertices.add_vertex(
                simplex.diagram_part().as_coords(), 15.0, color)
display.add(vertices)

display.run()
