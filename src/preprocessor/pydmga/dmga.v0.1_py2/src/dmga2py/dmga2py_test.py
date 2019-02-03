import sys
sys.path.append('../')

import dmga2py

from pydmga.draw import geometry
from pydmga.draw import assets
from pydmga.draw import render

print "[TEST] new geometry"
g1 = dmga2py.new_orthogonal_geometry(0, 0, 0, 5.0, 5.0, 5.0, 0, 0, 0)
g2 = dmga2py.new_orthogonal_geometry(0, 0, 0, 1.0, 1.0, 1.0, 1, 0, 1)
g3 = dmga2py.new_orthogonal_geometry(0, 0, 0, 1.0, 2.0, 9.0, 1, 1, 1)

print "[TEST] free objects"
dmga2py.free_object(g2)

print "[TEST] new container"
c1 = dmga2py.new_basic_container(g1)
c3 = dmga2py.new_basic_container(g3)

print "[TEST] free objects"
dmga2py.free_object(c3)
dmga2py.free_object(g3)

print "[TEST] ADDING"
balls = []
dmga2py.basic_container_add(c1, 1, 1.0, 0.5, 1.0, 0.5)
balls.append(geometry.Ball(1, 1.0, 0.5, 1.0, 0.5))
dmga2py.basic_container_add(c1, 2, 1.1, 0.5, 1.2, 0.5)
balls.append(geometry.Ball(2, 1.1, 0.5, 1.2, 0.5))
dmga2py.basic_container_add(c1, 3, 1.0, 2.5, 1.5, 0.5)
balls.append(geometry.Ball(3, 1.0, 2.5, 1.5, 0.5))
dmga2py.basic_container_add(c1, 4, 0.9, 4.1, 3.1, 0.9)
balls.append(geometry.Ball(4, 0.9, 4.1, 3.1, 0.9))

print "[TEST] new diagram"
d1 = dmga2py.new_diagram(c1, 1)

print "[TEST] get cells"
cell1 = dmga2py.diagram_get_cell(d1, 0)
cell2 = dmga2py.diagram_get_cell(d1, 1)
cell3 = dmga2py.diagram_get_cell(d1, 2)
cell4 = dmga2py.diagram_get_cell(d1, 3)

print "[TEST] cell stats"
vol1 = dmga2py.cell_get_volume(cell1)
vol2 = dmga2py.cell_get_volume(cell2)
vol3 = dmga2py.cell_get_volume(cell3)
vol4 = dmga2py.cell_get_volume(cell4)
sides1 = dmga2py.cell_get_sides_count(cell1)
vertices1 = dmga2py.cell_get_vertex_count(cell1)
print "[TEST] cell 1 volume = ", vol1
print "[TEST] cell 2 volume = ", vol2
print "[TEST] cell 3 volume = ", vol3
print "[TEST] cell 4 volume = ", vol4
print "[TEST] summary volume = ", vol1 + vol2 + vol3 + vol4, "should be", 5*5*5

print "[TEST] cell iterator"
iter = dmga2py.cell_iterator(cell1)
print "[TEST] current vertex is", dmga2py.cell_iterator_current_vertex(iter)
print "[TEST] its coords are", dmga2py.cell_iterator_current_vertex_coords(iter)

while (dmga2py.cell_iterator_is_finished(iter) == 0):
    dmga2py.cell_iterator_mark(iter)   
    dmga2py.cell_iterator_forward(iter)
    print "[TEST] arriving at", dmga2py.cell_iterator_current_vertex(iter), dmga2py.cell_iterator_current_vertex_coords(iter)
    if (dmga2py.cell_iterator_is_marked(iter) == 1):
        dmga2py.cell_iterator_jump(iter)
        print "[TEST] new side at", dmga2py.cell_iterator_current_vertex(iter), dmga2py.cell_iterator_current_vertex_coords(iter)

print "[TEST] Render 3D"
cells = []
dmgacells = [cell1, cell2, cell3, cell4]
for dcel in dmgacells:
    vertcount = dmga2py.cell_get_vertex_count(dcel)
    vertices = []
    sides = []
    for v in range(vertcount):
        coords = dmga2py.cell_get_vertex_coords(dcel, v)
        vertices.append( (coords[0], coords[1], coords[2], {"alpha": 0.0}) )
    iter = dmga2py.cell_iterator(dcel)
    polygon = []
    while (dmga2py.cell_iterator_is_finished(iter) == 0):        
        polygon.append( (dmga2py.cell_iterator_current_vertex(iter), {"alpha": 0.0}) )
        dmga2py.cell_iterator_mark(iter) 
        dmga2py.cell_iterator_forward(iter)
        if (dmga2py.cell_iterator_is_marked(iter) == 1):
            sides.append( (polygon, {"alpha": 0.0}) );
            polygon = []
            dmga2py.cell_iterator_jump(iter) 
    cells.append( (vertices, sides) )       
data = geometry.ArrayInput(balls, cells)
rend = render.VoroRender(data)
rend.run()

print "[END OF TESTS]"
