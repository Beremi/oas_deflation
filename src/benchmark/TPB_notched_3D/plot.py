#!/usr/bin/env python3
import numpy as np
import matplotlib
import matplotlib.cm as cm
import matplotlib.pyplot as plt
import os
from scipy import interpolate
from matplotlib.collections import LineCollection, PatchCollection

plt.rcParams.update({'font.size': 18})
plt.rcParams.update({'axes.linewidth': 2})
plt.rcParams.update({'text.usetex':True})
plt.rcParams.update({'font.family' : 'serif'})
plt.rcParams.update({'font.serif' : 'Times New Roman'})

fig = plt.figure(figsize=(5,5))
ax= fig.add_axes([0.15,0.13,0.82,0.85])

upper = np.loadtxt("results_to_check/TPB_A_exp_upper.dat")
lower = np.loadtxt("results_to_check/TPB_A_exp_lower.dat")
cus = np.loadtxt("results_to_check/TPB_A_model_Cusatis.dat")
data = np.loadtxt("results_to_check/LD.out",skiprows=1, usecols=[2,3])

f = interpolate.interp1d(lower[:,0], lower[:,1])

ax.fill_between(upper[:,0], f(upper[:,0]), upper[:,1],color="gray", label="experimental range", alpha=0.5)
ax.plot(cus[:,0],cus[:,1], ls=":", color="k",label="LDPM simulation",lw=4)
ax.plot(np.hstack((0,data[:,1]))*1000,np.hstack((0,-data[:,0]))/1000, color="k",label="our simulation",lw=4)

ax.set_xlim([0,1])
ax.set_ylim([0,3.8])
ax.set_xlabel("deflection [mm]")
ax.set_ylabel("force [kN]")

ax.legend(frameon=False)
plt.savefig("TPB_A.pdf")
plt.show()
            
            
        
