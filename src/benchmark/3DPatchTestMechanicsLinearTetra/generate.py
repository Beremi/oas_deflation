#!/usr/bin/env python3
#-*- coding: utf-8 -*-
import numpy as np
import os
import pdb, shutil
import gmsh
import sys
from gmshpy import generateMesh3D


if __name__ == "__main__":

    folder = "CUBE"
    if not os.path.isdir(folder): os.mkdir(folder)

    nodes,elements = generateMesh3D()
    inp = open(os.path.join(folder,'nodes.inp'),'w')
    for n in nodes:
        inp.write("MechNode\t%e\t%e\t%e\n"%(n[0],n[1],n[2]))
    inp.close()
    nn = len(nodes)    

    inp = open(os.path.join(folder,'elems.inp'),'w')
    for e in elements:
        inp.write("MechanicalTet\t%d\t%d\t%d\t%d\t%d\n"%(e[0],e[1],e[2],e[3],0))  
    inp.close()

    inp = open(os.path.join(folder,'functions.inp'),'w')
    inp.write("PWLFunction\t1\t0\t0\n")
    inp.write("PWLFunction\t1\t0\t0.01\n")  
    inp.close()

    inp = open(os.path.join(folder,'master.inp'),'w')
    inp.write("Dimension\t3\n")
    inp.write("Solver\tsolver.inp\n")
    inp.write("NodeFiles\t1\tnodes.inp\n")
    inp.write("MatFiles\t1\tmats.inp\n")
    inp.write("ElemFiles\t1\telems.inp\n")
    inp.write("BCFiles\t1\tbc.inp\n")
    inp.write("FunctionFiles\t1\tfunctions.inp\n")
    inp.write("ExporterFiles\t1\texporters.inp")
    inp.close()

    inp = open(os.path.join(folder,'exporters.inp'),'w')
    inp.write("VTKElementExporter	elems	saveEvery	0.01	pointData  1	displacements	extrapolatedNodeData   2 stress strain\n")
    inp.close()

    inp = open(os.path.join(folder,'solver.inp'),'w')
    inp.write("SteadyStateLinearSolver\n")
    inp.write("time_step	1\n")
    inp.write("total_time	1\n")
    inp.close()

    inp = open(os.path.join(folder,'mats.inp'),'w')
    inp.write("TensMechMaterial	E	10e9    nu    0.3	density	1.000000\n")
    inp.close()

    indB = np.where(nodes[:,0]<1e-4)[0]
    indT = np.where(nodes[:,0]>100.-1e-4)[0]

    inp = open(os.path.join(folder,'bc.inp'),'w')
    for i in indB:  
        if i==0: 
            inp.write("NodalBC\t%d\t0\t0\t0\t-1\t-1\t-1\n"%(i))
        elif i==1: 
            inp.write("NodalBC\t%d\t0\t-1\t0\t-1\t-1\t-1\n"%(i))
        elif i==3: 
            inp.write("NodalBC\t%d\t0\t0\t-1\t-1\t-1\t-1\n"%(i))
        else:
            inp.write("NodalBC\t%d\t0\t-1\t-1\t-1\t-1\t-1\n"%(i))
    for i in indT:
        inp.write("NodalBC\t%d\t1\t-1\t-1\t-1\t-1\t-1\n"%(i))
    inp.close() 



