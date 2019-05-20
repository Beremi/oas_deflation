import numpy as np
import matplotlib.pyplot as plt
import pydmga

data = [(0, 5.0, 5.0, 5.0, 1.0), 
        (1, 5.0-1.0, 5.0, 5.0, 1.0), 
        (2, 5.0, 5.0-1.0, 5.0, 1.0), 
        (3, 5.0, 5.0, 5.0-1.0, 1.0)]
geometry  = pydmga.geometry.OrthogonalGeometry(10, 10, 10, True, True, True);

#geometry  = OrthogonalGeometry(10, 10, 10, False, False, False);
container = pydmga.container.Container(geometry)
container.add(data)   


diagram = pydmga.diagram.Diagram(container, True)
print(diagram.size())
for i, cell in enumerate(diagram):
    neighbours = []
    
    for side in cell.sides:
        print(side.area())
        print(side.as_coords())

    print(cell.sides.size())

    for edge in cell.edges:
        print(edge.as_coords())
        print(edge.id)

    print(cell.edges.size())

    for vertex in cell.vertices:
        print(vertex.as_coords())
        print(vertex.id)

    print(cell.vertices.size())
        
    #for j, side in enumerate(cell.sides):
        #neighbours.append(side.neighbour)
        #print(side.as_list())
    #print(neighbours)
    
    #for vert in cell.vertices
        
