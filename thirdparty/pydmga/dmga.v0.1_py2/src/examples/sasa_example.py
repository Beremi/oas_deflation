import sys
sys.path.append('../')

from pydmga.geometry import OrthogonalGeometry
from pydmga.container import Container
from pydmga.diagram import Diagram
from pydmga.shape import SASAShape
from pydmga import draw
from pydmga.draw import Asset
from pydmga.draw import assets
from pydmga.draw import render
from pydmga.draw.render import Color
from random import random
from random import seed
from math import pi

(box_x, box_y, box_z) = (10.0, 10.0, 10.0);
# count = 8;
# max_radius = 2.0;
# data = [ (i, (box_x-max_radius*2.1) * random() + max_radius*1.05, (box_y-max_radius*2.1) * random() + max_radius*1.05, (box_z-max_radius*2.1) * random() + max_radius*1.05, 1.0 + (max_radius-1.0) * random()) for i in range(count)]
# #data = [(0, 8.99062649702234, 5.293750020205118, 8.647568140450064, 2.033918072045217), (1, 2.610463170143009, 5.729168616158611, 6.947404584325251, 2.9719289851933004), (2, 8.278264394611755, 5.1658481037516655, 3.9395652765651956, 1.252755927361921), (3, 4.706306247259862, 9.494252066699207, 2.2558959426829253, 2.206051037913281), (4, 8.807451088430975, 7.729281688474261, 4.243891066173331, 2.9048347217475605), (5, 4.362624091754058, 4.1701702708018615, 1.8065528298308509, 1.4575325914570898), (6, 0.10927285408377396, 3.6039247527316496, 6.0888565111353286, 2.9772356480057707), (7, 4.532135954362886, 3.2599778709992346, 7.336914472765875, 2.920787853134759)]
# data = [(0, 5.229342494403202, 3.587544145698443, 2.432157527229299, 1.9262072656175517), (1, 7.876368141712307, 3.4228698343361477, 7.205457277647881, 1.738350400074077), (2, 5.131821158173861, 4.589239991939696, 6.193628555085105, 1.9438020282245416), (3, 5.684876367233471, 5.610585504939861, 3.7573637085390845, 1.8092105315076568), (4, 7.0633503242095035, 7.29773100735671, 4.683022261442364, 1.4145257643039395), (5, 2.8443059554330086, 6.724028128725234, 5.571359601566978, 1.2509143255447508), (6, 7.6607412363629805, 5.815090768012357, 3.065635033608226, 1.4374958110420435), (7, 6.638715988552352, 7.270022341686143, 6.887485307514176, 1.9067959290577838)]
# print "Data used:", data

# count = 5
# data = [(0, 5.0, 5.0, 5.0, 1.0), 
#         (1, 5.0-1.0, 5.0, 5.0, 1.0), 
#         (2, 5.0+1.0, 5.0, 5.0, 1.0), 
#         (1, 5.0, 5.0-1.0, 5.0, 1.0), 
#         (2, 5.0, 5.0+1.0, 5.0, 1.0)]

count = 4
data = [(0, 5.0, 5.0, 5.0, 1.0), 
        (1, 5.0-1.0, 5.0, 5.0, 1.0), 
        (2, 5.0, 5.0-1.0, 5.0, 1.0), 
        (3, 5.0, 5.0, 5.0-1.0, 1.0)]
        #(3, 5.0-0.5, 5.0-0.5, 5.0-0.5, 1.0)]

# count = 3
# data = [(0, 5.0, 5.0, 5.0, 1.0), 
#         (1, 5.0-1.0, 5.0, 5.0, 1.0), 
#         (2, 5.0, 5.0-1.0, 5.0, 1.0)]

# count = 4
# data = [(0, 5.0, 5.0, 5.0, 1.0), 
#         (1, 5.0-0.75, 5.0-0.45, 5.0, 0.3), 
#         (2, 5.0-0.45, 5.0-0.75, 5.0, 0.3),
#         (3, 5.0-0.6, 5.0-0.6, 5.0+0.6, 0.38)]

# geometry  = OrthogonalGeometry(10, 10, 10, False, False, False);
# geometry  = OrthogonalGeometry(box_x, box_y, box_z, True, True, True);
geometry  = OrthogonalGeometry(box_x, box_y, box_z, True, True, True);
container = Container(geometry)
container.add(data)   

diagram = Diagram(container)
shape = SASAShape(diagram)

display = render.Render(Color.WHITE_OPAQUE)

area = shape.sas_area()
print "SASA Area:", area
colors = [Color.RED_OPAQUE]
for ball in container:
    (id, x, y, z, r) = ball
#    display.add(Asset(assets.sphere, (x,y,z), r, Color.BLUE_QUARTER))
#     display.add(Asset(assets.sphere, (x+box_x,y,z), r, (1.0, 1.0, 1.0, 0.1)))
#     display.add(Asset(assets.sphere, (x-box_x,y,z), r, (1.0, 1.0, 1.0, 0.1)))
#     display.add(Asset(assets.sphere, (x,y+box_y,z), r, (1.0, 1.0, 1.0, 0.1)))
#     display.add(Asset(assets.sphere, (x,y-box_y,z), r, (1.0, 1.0, 1.0, 0.1)))
#     display.add(Asset(assets.sphere, (x,y,z+box_z), r, (1.0, 1.0, 1.0, 0.1)))
#     display.add(Asset(assets.sphere, (x,y,z-box_z), r, (1.0, 1.0, 1.0, 0.1)))

for i, cell_shape in enumerate(shape):
    if i == 0:
        (id, x, y, z, r) = container.get(i)
        display.add(Asset(assets.sphere, (x,y,z), r, Color.BLUE_QUARTER))
        all_area = cell_shape.sas_area() + cell_shape.excluded_area()    
        print "area is", cell_shape.sas_area(), ", excluded area is", cell_shape.excluded_area()
        print "sum is", all_area, "should be", 4.0 * pi * (r**2);
        print cell_shape.area(), ",", cell_shape.volume()
        b = 0
        for a, side in enumerate(cell_shape):
            if side.neighbour >= 0:
                b += 1
                display.add(render.Outline(side.as_coords(), 2.0, Color.BLACK_OPAQUE))
        points = draw.PointCloud()
        great_arcs = render.GreatArcs((5,5,5),1.0,2.0,Color.BLUE_OPAQUE)
        for j, arc in enumerate(cell_shape.shape()):       
            print "Arc is: ", arc
            arcs = draw.SphereArcs((x,y,z), r, 2.0, colors[j % len(colors)])
            arcs.add_arc(arc.on_plane, arc.first, arc.second)
            #display.add(render.Lines([arc.first, arc.on_plane, arc.second], 2.0, Color.RED_OPAQUE))
            great_arcs.add_arcs([arc.first, arc.on_sphere, arc.second], False)
            points.add_vertex(arc.first, 10.0, Color.GREEN_OPAQUE)
            points.add_vertex(arc.second, 10.0, Color.GREEN_OPAQUE)
            display.add(arcs)
        display.add(great_arcs)
        display.add(points)

display.run()