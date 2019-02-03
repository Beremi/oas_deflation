import sys
sys.path.append('../')

from pydmga.geometry import OrthogonalGeometry
from pydmga.geometry import CastSphereGeometry
from pydmga.geometry import CastCylinderGeometry
from pydmga.container import Container
from pydmga.diagram import Diagram
from pydmga.render import geometry
from pydmga.render import graphics
from pydmga.render import render
from random import random
#from dmga.io.pdb import PDBReader  

solvent_from = -1
data = []
import fileinput
i = 0
cryst_line = ""
for line in fileinput.input():
    tag = line[0:6].strip()
    #print tag
    if (tag == 'CRYST1'):
        cryst_line = line
        box_x = float(line[6:15].strip())
        box_y = float(line[15:24].strip())
        box_z = float(line[24:33].strip())
        tilt_x = float(line[33:40].strip()) 
        tilt_y = float(line[40:47].strip()) 
        tilt_z = float(line[47:54].strip()) 
        group = line[55:66] 
        z_value = int(line[66:70].strip()) 
        #print "we have box:", (box_x, box_y, box_z), (tilt_x, tilt_y, tilt_z), group, z_value
    if (tag == 'ATOM'): 
        id = int(line[6:11].strip())
        atom_name = line[12:16]
        alt_loc = line[16:17]
        res_name = line[17:20]
        chain_id = line[21:22]
        seq_no = int(line[22:27])
        x = float(line[30:38].strip())
        y = float(line[38:46].strip())
        z = float(line[46:54].strip())
        r = float(line[54:60].strip())
        # drop rest...
        if (atom_name[1] != 'H'):
            data.append( (i, x, y, z, r) )            
            if (solvent_from < 0 and res_name.strip() == 'SOL'):
                solvent_from = i        
            i += 1
        #print "We have atom:", id, atom_name, alt_loc, res_name, chain_id, seq_no, (x, y, z, r)
        
geometry = OrthogonalGeometry(box_x, box_y, box_z, True, True, True)
container = Container(geometry)
#print "Data lenght:", solvent_from, len(data)
container.add(data)
diagram = Diagram(container)

#display = render.Render()
colors = []
surface_by_part = []
surface_by_part_half = []
mark = [False for i in range(len(data)+1)]
print len(mark)
for i in range(solvent_from, len(data)):
    if (mark[i]): continue
    queue = [i]
    #print "new part at", i
    mark[i] = True        
    surface_by_part.append([])
    surface_by_part_half.append([])
    while(len(queue) > 0):
        i = queue.pop()                
        cell = diagram.get_cell(i)      
        for side in cell.sides:
            neighbour = side.neighbour
            polygon = side.as_coords()
            if (neighbour < solvent_from):                
                if (neighbour < 50 * 160):
                    surface_by_part_half[-1].append(0)
                else:
                    surface_by_part_half[-1].append(1)
                surface_by_part[-1].append(polygon) ##add new polygon to this part
            else:
                if (not(mark[neighbour])):
                    mark[neighbour] = True 
                    queue.append(neighbour)

outer = surface_by_part[0]
surface_by_part[0] = []
surface_by_part.insert(0, [])
for i in range(len(outer)):
    surface_by_part[surface_by_part_half[0][i]].append(outer[i])

import sys
filename = ".".join(sys.argv[1].split(".")[0:-1])   
print "filename is", filename 
pa_id = 100000000
#print "There are", len(surface_by_part), "parts of water in this system."
f = None
for i, part in enumerate(surface_by_part):
    file_count = 0
    for polygon in part:
        if (pa_id + len(polygon) >= 9999):
            if (f is not None): f.close()
            file_count += 1
            f = open("/home/robson/{0}_{1:03d}_{2:03d}.pdb".format(filename, i, file_count), "w")
            f.write(cryst_line)
            pa_id = 0
        #print last vertex of polygon 
        to_connect = []   
        for v in polygon:
            pa_id += 1
            to_connect.append(pa_id)
            if (i == 0):
                if v[2] < 30.0:
                    v = (v[0], v[1], v[2] + box_z)
            if (i == 1):
                if v[2] > 30.0:
                    v = (v[0], v[1], v[2] - box_z)
            f.write("HETATM{0: 5d} PS1  PSD P   1    {1: 8.3f}{2: 8.3f}{3: 8.3f}  0.00  0.00      PSDOPS  \n".format(pa_id, v[0], v[1], v[2]))
        to_connect.insert(0, pa_id)
        for j in range(len(polygon)):
            f.write("CONECT{0: 5d}{1: 5d}\n".format(to_connect[j], to_connect[j+1]))
            
                    
    pa_id = 100000000
        
if not f.closed:
    f.close()
    
#set connect_cutoff,0.0
#set connect_mode,1
   
        
#display.add_outline_polygon(polygon, color)            
#display.run()
