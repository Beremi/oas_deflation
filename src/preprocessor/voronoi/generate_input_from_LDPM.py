#!/usr/bin/env python3
#-*- coding: utf-8 -*-

# created by Jan Elias
# jan.elias@vut.cz
# Brno University of Technology,
# 2022

import numpy as np
import os

def is_between(n1, n2, x):
    v = n1-n2
    d1 = x - n1
    v = v/np.linalg.norm(v)
    d1 = d1/np.linalg.norm(d1)
    if abs(d1.dot(v))>1-1e-8:
        return True
    else: return False

def is_on_surface(n1, n2, n3, x):
    pairs = [[n1,n2, n3],[n1, n3, n2],[n2, n3, n1]]
    for i in range(1):
        v = pairs[i][0]-pairs[i][2]
        w = pairs[i][1]-pairs[i][2]
        v = v/np.linalg.norm(v)
        w = w/np.linalg.norm(w)
        norm = np.cross(v,w)
        norm = norm/np.linalg.norm(norm) 

        v = pairs[i][0]-x
        w = pairs[i][1]-x
        v = v/np.linalg.norm(v)
        w = w/np.linalg.norm(w)
        normX = np.cross(v,w)
        normX = normX/np.linalg.norm(normX)

        if abs(normX.dot(norm))>1-1e-8:
            return True
    return False

if __name__ == '__main__':

    
    LDPM_folder_name = "LDPM_generator"

    gfileending = ["-data-mesh.dat","-data-facet.dat"]
    master_name = ""
    for n in os.listdir(LDPM_folder_name):
        if n.endswith(gfileending[0]):
            master_name = n.split("-")[0]
    if master_name == "": 
        print("Error: LDPM folder ", LDPM_folder_name, " does not contains geometry file ending with ", gfileending[0])
        exit(1)
    export_folder = master_name;

    if (not os.path.isdir(export_folder)): os.mkdir(export_folder)
       
    #load nodes
    #load tetras
    with open(os.path.join(LDPM_folder_name, master_name + gfileending[0])) as f:
        while 1: 
            line = f.readline()
            if not line: break
            line = line.strip()
            if (line.startswith("Number of Vertices:")):
                nver = int(line.split()[-1])
                nodes = np.zeros((nver,3))
                i = 0
                while i < nver:
                    line = f.readline().strip()
                    if not line=="" and (line[0].isdigit() or line[0]=="-"):                        
                        nodes[i] = np.array(line.split()).astype(float) 
                        i+=1
            if (line.startswith("Number of Tets:")):
                ntets = int(line.split()[-1])
                #fi = open(os.path.join(export_folder, "elements.inp"),"w")
                tetras = np.zeros((ntets,15)).astype(int)-1
                i = 0
                while i < ntets:
                    line = f.readline()
                    if not line=="" and line[0].isdigit():
                        tetras[i,:4] = np.array(line.split()).astype(int)-1
                        i+=1
          
    centroids = np.zeros((0,3)) 
    linenodes = np.zeros((0,3))
    linecodes = np.zeros((0,2)).astype(int)
    surfnodes = np.zeros((0,3))
    surfcodes = np.zeros((0,3)).astype(int)   
  
    #load facet structures
    with open(os.path.join(LDPM_folder_name, master_name + gfileending[1])) as f:
        while 1: 
            line = f.readline()
            if not line: break
            line = line.strip()
            if (line.startswith("Tet index:")):
                teti = int(line.split()[-1])-1
                tetnodes = tetras[teti][:4]
                tetnodesX = nodes[tetnodes]
 
                for k in range(12): #12 facets
                
                    #simplex data
                    m = f.readline().split();
                    n1 = int(m[0])-1; #node numbers
                    n2 = int(m[1])-1;
                    n1num = np.where(tetnodes == n1)[0]
                    n2num = np.where(tetnodes == n2)[0]
                    if len(n1num)==1: n1num = n1num[0]
                    else:
                        print("Facet node not found in tetrahedron")
                        exit(1);
                    if len(n2num)==1: n2num = n2num[0]
                    else:
                        print("Facet node not found in tetrahedron")
                        exit(1);
                    remainnum = [x for x in range(4) if not x in [n1num, n2num]]

                    #read nodes
                    f.readline(); #local coord system
                    m = np.array(f.readline().split()).astype(float);
                    mred = [m[3:6], m[6:9]] #centroid is always first

                    #determine line node
                    isonline = np.zeros(3).astype(bool)
                    for im, mm in enumerate(mred):
                        isonline[im] = is_between(tetnodesX[n1num], tetnodesX[n2num], mm)
                        if isonline[im]: break #speedup, not safe
                    linenode = np.where(isonline)[0]
                    if len(linenode)==1: linenode = linenode[0]
                    else:
                        print("Line node not found")
                        exit(1);
                    
                    #determine surface node
                    isonsurf = np.zeros(6).astype(bool)                    
                    for im, mm in enumerate(mred):
                        if im==linenode: continue
                        isonsurf[im*2]   = is_on_surface(tetnodesX[n1num], tetnodesX[n2num], tetnodesX[remainnum[0]], mm)
                        if isonsurf[im]: break #speedup, not safe
                        isonsurf[im*2+1] = is_on_surface(tetnodesX[n1num], tetnodesX[n2num], tetnodesX[remainnum[1]], mm)
                        if isonsurf[im]: break #speedup, not safe
                    surfnode = np.where(isonsurf)[0]
                    if len(surfnode)==1: surfnode = surfnode[0]
                    else:
                        print("Surface node not found")
                        exit(1);
                    surfnum = surfnode%2
                    surfnode = int((surfnode-surfnum)/2.+0.5)

                    #centroid
                    centnode = 0
                    while surfnode==centnode or linenode==centnode: centnode +=1
    
                    #add_to_auxnodes
                    #add centroid
                    if k==0: 
                        centroids = np.vstack((centroids, m[:3]))
                        tetras[teti][14] = len(centroids)-1
                    
                    #add line node                    
                    mink = min(n1num,n2num)
                    maxk = max(n1num,n2num)
                    h = 0
                    if mink==0: h = 3+maxk
                    elif mink==1: h = 5+maxk
                    else: h = 9
                    if (tetras[teti][h]==-1): 
                        lnodes = np.sort(tetnodes[np.array([mink,maxk])])
                        ind1 = np.where(linecodes[:,0]==lnodes[0])[0]
                        ind2 = np.where(linecodes[ind1,1]==lnodes[1])[0]
                        if len(ind2)==0: 
                            linecodes = np.vstack((linecodes, lnodes))
                            linenodes = np.vstack((linenodes, mred[linenode]))                        
                            tetras[teti][h] = len(linecodes)-1
                        else:
                            #print("found line")
                            tetras[teti][h] = ind1[ind2]

                    #add surf node                    
                    snodes = np.sort(tetnodes[np.array([n1num, n2num, remainnum[surfnum]])])
                    h = n1num + n2num + remainnum[surfnum] + 7
                    if (tetras[teti][h]==-1): 
                        ind1 = np.where(surfcodes[:,0]==snodes[0])[0]
                        ind2 = np.where(surfcodes[ind1,1]==snodes[1])[0]
                        ind1 = ind1[ind2]
                        ind3 = np.where(surfcodes[ind1,2]==snodes[2])[0]
                        if len(ind3)==0: 
                            surfcodes = np.vstack((surfcodes, snodes))
                            surfnodes = np.vstack((surfnodes, mred[surfnode]))                        
                            tetras[teti][h] = len(surfcodes)-1
                        else:
                            #print("found surf")
                            tetras[teti][h] = ind1[ind3] 
                if (-1 in tetras[teti]):
                    print("Tetrahedron reader failed")
                    print(tetras[teti])
                    exit(1)
                if(teti%1000==0 and teti>0): 
                    print("processed ", teti*100./len(tetras), "%")
                    #break
                #print(tetras[teti])       
   
    nodes = nodes/1000. #mm
    auxnodes = np.vstack((centroids,linenodes,surfnodes))/1000. #mm
    lenn  = len(nodes)
    lenc = len(centroids)
    lenl = len(linenodes)
    fi = open(os.path.join(export_folder, "nodes.inp"),"w")
    for n in nodes:
        fi.write("Particle\t%.10e\t%.10e\t%.10e\n"%(n[0],n[1],n[2]))
    fi.close()
    fi = open(os.path.join(export_folder, "auxNodes.inp"),"w")
    for n in auxnodes:
        fi.write("AuxNode\t%.10e\t%.10e\t%.10e\n"%(n[0],n[1],n[2]))
    fi.close()
    fi = open(os.path.join(export_folder, "mechElems.inp"),"w")
    for t in tetras:
        fi.write("LDPMTetra\t%d\t%d\t%d\t%d"%(t[0],t[1],t[2],t[3]))
        fi.write("\t%d\t%d\t%d\t%d\t%d\t%d"%(t[4]+lenn+lenc,t[5]+lenn+lenc,t[6]+lenn+lenc,t[7]+lenn+lenc,t[8]+lenn+lenc,t[9]+lenn+lenc))
        fi.write("\t%d\t%d\t%d\t%d\t%d\t%d\n"%(t[10]+lenn+lenc+lenl,t[11]+lenn+lenc+lenl,t[12]+lenn+lenc+lenl,t[13]+lenn+lenc+lenl,t[14]+lenn,0))
    fi.close()    
                    
    print(len(tetras), "tertahedrons, ", len(nodes), "nodes, ", len(centroids), "centroids, ", len(linenodes), "line nodes, ", len(surfnodes), "surface nodes")
