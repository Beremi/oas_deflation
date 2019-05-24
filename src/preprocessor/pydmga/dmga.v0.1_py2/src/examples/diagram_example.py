import sys
sys.path.append('../')

from pydmga.geometry import OrthogonalGeometry
from pydmga.container import Container
from pydmga.diagram import Diagram
from pydmga.draw import assets
from pydmga.draw import render
from random import random
from random import seed
from pydmga.draw.render import Color

(box_x, box_y, box_z) = (10.0, 10.0, 10.0);
count = 8;    
data = [ (i, box_x * random(), box_y * random(), box_z * random(), random()) for i in range(count)]

count = 4
data = [(0, 5.0, 5.0, 5.0, 1.0), 
        (1, 5.0-1.0, 5.0, 5.0, 1.0), 
        (2, 5.0, 5.0-1.0, 5.0, 1.0), 
        (3, 5.0, 5.0, 5.0-1.0, 1.0)]
geometry  = OrthogonalGeometry(10, 10, 10, True, True, True);

#geometry  = OrthogonalGeometry(10, 10, 10, False, False, False);
container = Container(geometry)
container.add(data)   

diagram = Diagram(container, True)

display = render.Render(Color.WHITE_OPAQUE)

volume = 0.0
colors = [(1,0,0,0.3), (0,1,0,0.3), (0,0,1,0.3), (1,1,0,0.3), (0,1,1,0.3), (1,0,1,0.3)]
for i, cell in enumerate(diagram):
    print "Entering cell", i
    (id, x, y, z, r) = container.get(i)
    display.add(render.Asset(assets.sphere, (x, y, z), r, Color.BLUE_QUARTER))
#    volume += cell.volume()
    neighbours = []
    all_area = 0
    for j, side in enumerate(cell.sides):
        polygon = side.as_coords()
        neighbours.append(side.neighbour)                
        display.add(render.Outline(polygon))
        display.add(render.ConvexPolygons(polygon, colors[i % len(colors)]))
#        print "area of side", j, "is", side.area()
        all_area += side.area()
    volume += cell.volume()
    print "Sum of area sides is", all_area, "should be", cell.area()
    print "Volume of cell is", cell.volume()
    print "Neighbours of cell are", neighbours
    print "Leaving cell", i
    
print "Volume is", volume, "should be", box_x * box_y * box_z

display.run()