#!/usr/bin/env python3
import os
import shutil
import numpy as np
import subprocess

def random_three_vector():
    phi = np.random.uniform(0,np.pi*2)
    costheta = np.random.uniform(-1,1)

    theta = np.arccos( costheta )
    x = np.sin( theta) * np.cos( phi )
    y = np.sin( theta) * np.sin( phi )
    z = np.cos( theta )
    return np.array([x,y,z])

np.random.seed(6)
normal = random_three_vector()
t1 = random_three_vector()
t1 = t1 - np.dot(normal,t1)*normal
t1 = t1/np.linalg.norm(t1)
t2 = np.cross(normal, t1)


folder = "rotated"
if os.path.isdir(folder):
    shutil.rmtree(folder)
os.mkdir(folder)

load = np.loadtxt("functions.inp",skiprows = 1, usecols = 3)
file = open(os.path.join(folder,"functions.inp"),"w")
file.write("PWLFunction	1	0   0" + os.linesep)
file.write("PWLFunction	1	0   %e"%(load[0]*normal[0]+load[1]*t1[0]+load[2]*t2[0]) + os.linesep)
file.write("PWLFunction	1	0   %e"%(load[0]*normal[1]+load[1]*t1[1]+load[2]*t2[1]) + os.linesep)
file.write("PWLFunction	1	0   %e"%(load[0]*normal[2]+load[1]*t1[2]+load[2]*t2[2]) + os.linesep)
file.write("PWLFunction	1	0   %e"%(load[3]*normal[0]+load[4]*t1[0]+load[5]*t2[0]) + os.linesep)
file.write("PWLFunction	1	0   %e"%(load[3]*normal[1]+load[4]*t1[1]+load[5]*t2[1]) + os.linesep)
file.write("PWLFunction	1	0   %e"%(load[3]*normal[2]+load[4]*t1[2]+load[5]*t2[2]) + os.linesep)
file.close()

points = np.loadtxt("nodes.inp",skiprows = 0, usecols = [1,2,3])
shift = np.random.rand(3)
file = open(os.path.join(folder,"nodes.inp"),"w")
p1 = shift + points[0,0]*normal+points[0,1]*t1+points[0,2]*t2
p2 = shift + points[1,0]*normal+points[1,1]*t1+points[1,2]*t2
file.write("Particle	%e	%e %e"%(p1[0],p1[1],p1[2]) + os.linesep)
file.write("Particle	%e	%e %e"%(p2[0],p2[1],p2[2]) + os.linesep)
file.close()

file = open(os.path.join(folder,"elems.inp"),"w")
file.write("TimoshenkoBeam3D 0 1 0 crosssection 0 zrefpoint %e %e %e"%(t2[0],t2[1],t2[2]) + os.linesep)
file.close()

for k in ["solver.inp","mechBC.inp","materials.inp","master.inp","exporters.inp","cross_sections.inp"]:
    shutil.copy(k,os.path.join(folder,k))


sampleA = subprocess.Popen(["/media/jan/Data/VYPOCTY/OAS/binaries/bin/OAS", 'master.inp'])
sampleB = subprocess.Popen(["/media/jan/Data/VYPOCTY/OAS/binaries/bin/OAS", os.path.join(folder,'master.inp')])
sampleA.communicate()
sampleB.communicate()

displA = np.loadtxt(os.path.join("results","nodes_00001.out"),skiprows = 2, usecols = [1,2,3,4,5,6])
displB = np.loadtxt(os.path.join(folder,"results","nodes_00001.out"),skiprows = 2, usecols = [1,2,3,4,5,6])
displB = np.array([displB[:3].dot(normal),displB[:3].dot(t1),displB[:3].dot(t2),displB[3:].dot(normal),displB[3:].dot(t1),displB[3:].dot(t2)])
print("displacement:")
print(np.column_stack((displA,displB,np.divide((displB-displA),displA))))
print("internal force:")
intfA = np.loadtxt(os.path.join("results","internalforces_00001.out"),skiprows = 1, usecols = [5,6,7,8,9,10])
intfB = np.loadtxt(os.path.join(folder,"results","internalforces_00001.out"),skiprows = 1, usecols = [5,6,7,8,9,10])
for k in range(len(intfA)):
    print("IP ",k)
    print(np.column_stack((intfA[k],intfB[k],np.divide((intfB[k]-intfA[k]),intfA[k]))))

#shutil.rmtree(folder)
