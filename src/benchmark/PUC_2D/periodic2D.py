#!/usr/bin/env python3
import numpy as np
import matplotlib
import os
import matplotlib.cm as cm
import matplotlib.pyplot as plt
from scipy.spatial import Voronoi, Delaunay

if 1:
    nsteps= 100
    saveperiod= 1
    dt = 24.*3600./nsteps
    d = 1
    w = 1
    lmin = 0.05
    nsteps = 1
    seed = 5
    Ex = 0
    Ey = 0.1
    Gxy = 0.2

    np.random.seed(seed)

    nodes = np.empty((0,2))
    ids = np.empty(0).astype(int)
    mirtype = np.empty(0).astype(int)    

    i = 0
    mastersize = np.array([w,d])

    periodicband = 3*lmin

    #NODES    
    while i<max(1000,len(nodes)):
        point = (np.random.rand(2)-0.5)*mastersize+mastersize/2.
        if len(nodes) == 0: dist = 1E6
        else: dist = min(np.sum(np.square(nodes-point),1))
        if dist>lmin*lmin:
            xplus = False
            yplus = False
            xmins = False
            ymins = False
            k = len(nodes)
            nodes = np.vstack((nodes, point))
            ids = np.hstack((ids, 0))
            mirtype = np.hstack((mirtype, 0))
            if(point[0]<periodicband): xplus = True
            if(point[0]>w-periodicband): xmins = True
            if(point[1]<periodicband): yplus = True
            if(point[1]>d-periodicband): ymins = True
            
            if xplus:
                nodes = np.vstack((nodes, point))
                nodes[-1][0] += w
                ids = np.hstack((ids, k))
                mirtype = np.hstack((mirtype, 1))
            if xmins:
                nodes = np.vstack((nodes, point))
                nodes[-1][0] -= w
                ids = np.hstack((ids, 0))   
                mirtype = np.hstack((mirtype, -1))
            if yplus:
                nodes = np.vstack((nodes, point))
                nodes[-1][1] += d
                ids = np.hstack((ids, k))  
                mirtype = np.hstack((mirtype, 2))   
            if ymins:
                nodes = np.vstack((nodes, point))
                nodes[-1][1] -= d
                ids = np.hstack((ids, 0))   
                mirtype = np.hstack((mirtype, -1))  
            if xmins and ymins:
                nodes = np.vstack((nodes, point))
                nodes[-1][0] -= w; nodes[-1][1] -= d 
                ids = np.hstack((ids, 0))
                mirtype = np.hstack((mirtype, -1))
            if xmins and yplus:
                nodes = np.vstack((nodes, point))
                nodes[-1][0] -= w; nodes[-1][1] += d 
                ids = np.hstack((ids, 0)) 
                mirtype = np.hstack((mirtype, -1))
            if xplus and ymins:
                nodes = np.vstack((nodes, point))
                nodes[-1][0] += w; nodes[-1][1] -= d 
                ids = np.hstack((ids, 0))
                mirtype = np.hstack((mirtype, -1))
            if xplus and yplus:
                nodes = np.vstack((nodes, point))
                nodes[-1][0] += w; nodes[-1][1] += d 
                ids = np.hstack((ids, k)) 
                mirtype = np.hstack((mirtype, 3))

            i = 0
        else: i += 1

    N = nodes.shape[0]

    v = Voronoi(nodes)

    valid_ridge_nodes = np.empty((0,2)).astype(int)
    valid_ridge_vertices = np.empty((0,2)).astype(int)
    for ir,r in enumerate(v.ridge_points): 
        if (mirtype[r[0]]==0 and mirtype[r[1]]>=0) or (mirtype[r[1]]==0 and mirtype[r[0]]>=0) or mirtype[r[0]]*mirtype[r[1]]==2: 
            valid_ridge_nodes = np.vstack((valid_ridge_nodes,r))
            valid_ridge_vertices = np.vstack((valid_ridge_vertices,v.ridge_vertices[ir]))
    valid_points = np.unique(valid_ridge_nodes.flatten())
    inverse = np.zeros(len(nodes)).astype(int)
    for k in range(len(valid_points)): inverse[valid_points[k]] = k  
    valid_ridge_nodes = inverse[valid_ridge_nodes]
    valid_ids = ids[valid_points]
    valid_ids = inverse[valid_ids]
    mirtype = mirtype[valid_points]
    valid_nodes = nodes[valid_points]


    valid_vertices = np.unique(valid_ridge_vertices.flatten())
    inverse = np.zeros(len(v.vertices)).astype(int)
    for k in range(len(valid_vertices)): inverse[valid_vertices[k]] = k  
    valid_ridge_vertices = inverse[valid_ridge_vertices]
    valid_vertices = v.vertices[valid_vertices]
        

    #plotting
    #"""
    fig = plt.figure(figsize=[7,7])
    ax = fig.add_axes([0.05,0.05,0.8,0.9])
    for i in valid_nodes:
        ax.plot(i[0],i[1],'ok', alpha = 0.5)
    for i in valid_ridge_nodes:
        ax.plot([valid_nodes[i[0],0],valid_nodes[i[1],0]],[valid_nodes[i[0],1],valid_nodes[i[1],1]],'k', alpha = 0.5)
    for i in valid_vertices:
        ax.plot(i[0],i[1],'ob', alpha = 0.5)
    for i in valid_ridge_vertices:
        ax.plot([valid_vertices[i[0],0],valid_vertices[i[1],0]],[valid_vertices[i[0],1],valid_vertices[i[1],1]],'k', alpha = 0.5)
    for i in range(len(mirtype)): 
        if mirtype[i]>0:
            ax.plot([valid_nodes[i,0],valid_nodes[valid_ids[i],0]],[valid_nodes[i,1],valid_nodes[valid_ids[i],1]],'g', alpha = 0.5)
    plt.show()
    #"""


    cf = open("master.inp","w") 
    cf.write("Dimension 2"+os.linesep)
    cf.write("Solver	SteadyStateNonLinearSolver	time_step	0.1	total_time	16"+os.linesep)
    cf.write("NodeFiles	1	nodes.inp"+os.linesep)
    cf.write("MatFiles	1	materials.inp"+os.linesep)
    cf.write("ElemFiles	1	elems.inp"+os.linesep)
    cf.write("PBlockFiles	1	blocks.inp"+os.linesep)
    cf.write("FunctionFiles	1	functions.inp"+os.linesep)
    cf.write("ExporterFiles	1	exporters.inp"+os.linesep)
    cf.close()

    cf = open("nodes.inp","w") 
    for i in valid_nodes:
        cf.write("Particle\t%.10e\t%.10e"%(i[0],i[1])+os.linesep)
    for i in valid_vertices:
        cf.write("AuxNode\t%.10e\t%.10e"%(i[0],i[1])+os.linesep)
    cf.close()

    cf = open("materials.inp","w") 
    cf.write("FatigueShearMaterial	E0	43.0e9	alpha	0.300000    density 2200.0 tauBar 4.0e6 Kin 0.0 gamma 10.0e6 S 0.0025e6 m 0"+os.linesep)
    #cf.write("DisMechMaterial	E0	43.0e9	alpha	0.300000    density 2200.0"+os.linesep)
    cf.close()

    cf = open("elems.inp","w") 
    nnodes = len(valid_nodes)
    for i in range(len(valid_ridge_vertices)):
        cf.write("LTCBEAM\t%d\t%d\t2\t%d\t%d\t0"%(valid_ridge_nodes[i,0],valid_ridge_nodes[i,1],valid_ridge_vertices[i,0]+nnodes,valid_ridge_vertices[i,1]+nnodes)+os.linesep)
    cf.close()

    
    cf = open("blocks.inp","w") 
    ndepend = len(np.where(mirtype>0)[0])
    #ex ey gxy sx sy sxy
    cf.write("BasicPeriodicBC\tsize\t2\t%e\t%e\tload\t%d\tex\t%d\tey\t%d\tgxy\t%d\tpairs\t%d"%(w,d,3,0,1,2,ndepend))
    for i in range(len(mirtype)): 
        if mirtype[i]>0: cf.write("\t%d\t%d"%(i,valid_ids[i]))
    cf.write(os.linesep)
    cf.close()
    
    cf = open("functions.inp","w") 
    cf.write("PWLFunction	1	0.000000	0"+os.linesep)
    cf.write("ConstSawToothFn value -5e-4 period 4"+os.linesep)
    cf.close()

    cf = open("exporters.inp","w") 
    cf.write("VTKElementExporter out  saveEvery 1e-20 0"+os.linesep)
    cf.close()

