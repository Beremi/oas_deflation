#!/usr/bin/env python3
#-*- coding: utf-8 -*-
import numpy as np
import os
import pdb, shutil
import gmsh
import sys

def generateMesh3D():
     
    # Initialize gmsh:
    gmsh.initialize()

    # body:
    factor = 20
    point1 = gmsh.model.geo.add_point(0, 0, 0, factor)
    point2 = gmsh.model.geo.add_point(0, 100, 0, factor)
    point3 = gmsh.model.geo.add_point(0, 100, 100, factor)
    point4 = gmsh.model.geo.add_point(0, 0, 100, factor)
    point5 = gmsh.model.geo.add_point(100, 0, 0, factor)
    point6 = gmsh.model.geo.add_point(100, 100, 0, factor)
    point7 = gmsh.model.geo.add_point(100, 100, 100, factor)
    point8 = gmsh.model.geo.add_point(100, 0, 100, factor)
    line1 = gmsh.model.geo.add_line(point1,point2)
    line2 = gmsh.model.geo.add_line(point2,point3)
    line3 = gmsh.model.geo.add_line(point3,point4)
    line4 = gmsh.model.geo.add_line(point4,point1)
    line5 = gmsh.model.geo.add_line(point5,point6)
    line6 = gmsh.model.geo.add_line(point6,point7)
    line7 = gmsh.model.geo.add_line(point7,point8)
    line8 = gmsh.model.geo.add_line(point8,point5)
    line9 = gmsh.model.geo.add_line(point1,point5)
    line10 = gmsh.model.geo.add_line(point2,point6)
    line11 = gmsh.model.geo.add_line(point3,point7)
    line12 = gmsh.model.geo.add_line(point4,point8)

    face1 = gmsh.model.geo.add_curve_loop([line1, line2, line3, line4])
    face2 = gmsh.model.geo.add_curve_loop([line5, line6, line7, line8])
    face3 = gmsh.model.geo.add_curve_loop([line1, line10, -line5, -line9])
    face4 = gmsh.model.geo.add_curve_loop([line2, line11, -line6, -line10])
    face5 = gmsh.model.geo.add_curve_loop([line3, line12, -line7, -line11])
    face6 = gmsh.model.geo.add_curve_loop([line4, line9, -line8, -line12])

    pl1 = gmsh.model.geo.addPlaneSurface([face1])
    pl2 = gmsh.model.geo.addPlaneSurface([face2])
    pl3 = gmsh.model.geo.addPlaneSurface([face3])
    pl4 = gmsh.model.geo.addPlaneSurface([face4])
    pl5 = gmsh.model.geo.addPlaneSurface([face5])
    pl6 = gmsh.model.geo.addPlaneSurface([face6])

    l = gmsh.model.geo.addSurfaceLoop([pl1,pl2,pl3,pl4,pl5,pl6])
    gmsh.model.geo.addVolume([l])


    #gmsh.model.geo.recombine(face1)
    #Mesh.Smoothing = 100;
    #Mesh 1;
    #Save "TDCB.msh";

    field = gmsh.model.mesh.field
    #field.add("MathEval", 1)
    #field.setString(1, "F", "0.01*(1.0+30.*(y-x*x)*(y-x*x) + (1-x)*(1-x))")
    #field.setAsBackgroundMesh(1)


    # To generate quadrangles instead of triangles
    #gmsh.option.setNumber("Mesh.RecombineAll", 1)

     
    # Create the relevant Gmsh data structures
    # from Gmsh model.
    gmsh.model.geo.synchronize()
     
    # Generate mesh:
    gmsh.model.mesh.generate(3)
     
    # Write mesh data:
    mesh = gmsh.write("TDCB.vtk")

    # Creates  graphical user interface
    #if 'close' not in sys.argv:
    #    gmsh.fltk.run()

    gmsh.finalize()
    
    ###########################################
    file1 = open("TDCB.vtk","r")
    lines = file1.readlines()
    file1.close()

    nl = 0
    el = 0
    for i,l in enumerate(lines):
        if "POINTS" in l: nl = i
        if "CELLS" in l: el = i
    nn = int(lines[nl].split()[1])
    ne = int(lines[el].split()[1])

    nodes = np.zeros((nn,3))
    i = nl
    k = 0
    while 1:
        i+=1
        ll = lines[i].split()
        if len(ll)==3:
            nodes[k,0] = float(ll[0])
            nodes[k,1] = float(ll[1])
            nodes[k,2] = float(ll[2])
            k+=1
            if k==nn: break

    elements = np.zeros((0,4))
    i = el
    k = 0
    while 1:
        i+=1
        if i>=len(lines): break
        ll = lines[i].split()
        if len(ll)==5:
            elements = np.vstack((elements,np.flip(np.array(ll)[1:])))
            k+=1  

    elements = elements.astype(int)

    

    print(len(nodes),"nodes",len(elements),"elements")

    return nodes, elements
    

if __name__ == "__main__":
    nodes,elements = generateMesh3D()
    


