import numpy as np
import matplotlib.pyplot as plt
import scipy
import math
from IPython.display import clear_output



#
##################################################
#### MECHANICAL BEAM ELEMENT FOR LATTICE MODEL ###
class latticeBeam:
    def __init__ (self, dim, nodeAidx, nodeBidx,  materialIdx):
        self.type = 0 #latticeBeam
        self.dim = dim
        self.nodeA = nodeAidx
        self.nodeB = nodeBidx

        self.material = materialIdx

    def getString(self):
        line = '%d'%(self.nodeA)  + '\t' + '%d'%(self.nodeB) + '\t' +  '%d'%(self.material)
        return line

##################################################






##################################################
#### TRANSPORT PATH ELEMENT FOR LATTICE MODEL ####
class transportPath2d:
    def __init__ (self,  vertexAidx, vertexBidx, materialIdx):
        self.vertexA = vertexAidx
        self.vertexB = vertexBidx

        self.material = materialIdx

    def getString(self):
        line = '%d'%(self.vertexA)  + '\t' + '%d'%(self.vertexB) + '\t' + '%d'%(self.material)
        return line

##################################################



# #matType    YoungM    Poisson    TranspC    TranspS    Density
#DisMechMaterial    E0 3.000000e+10    alpha 0.300000    density 2200.000000 TrsprtMaterial capacity 11.000000    conductivity 22.000000

##################################################
#### TRANSPORT PATH ELEMENT FOR LATTICE MODEL ####
class linearElasticMaterial:
    def __init__ (self, youngModulus, poisson, transportC, transportS, density):
        self.youngModulus = youngModulus
        self.poisson = poisson
        self.transportC = transportC
        self.transportS = transportS
        self.density = density

    def getString(self):
       # line = 'LINEL\t' + '%e'%(self.youngModulus)  + '\t' + '%f'%(self.poisson) + '\t' + '%f'%(self.transportC) + '\t' + '%f'%(self.transportS) + '\t' + '%f'%(self.density)
        line = 'DisMechMaterial'       + '\t' + 'E0\t%e'%(self.youngModulus)          + '\t' + 'alpha\t%f'%(self.poisson)         + '\nTrsprtMaterial'         + '\t' + 'density\t%f'%(self.density)         + '\t' + 'capacity\t%f'%(self.transportC)         + '\t' + 'conductivity\t%f'%(self.transportS)

        return line

##################################################


##################################################
#### General function set by table ####
class generalFunc:
    def __init__(self, table):
        self.table = table

    def getString(self):
        line = 'PWLFunction\t%d'%(len(self.table))

        for i in range (len(self.table)):
            line += '\t%f'%(self.table[i][0])

        for i in range (len(self.table)):
            line += '\t%f'%(self.table[i][1])

        return line
####################################################


##################################################
#### Sine function ####
class sineFunc:
    def __init__ (self, amplitude, freq):
        self.amplitude = amplitude
        self.freq = freq
    def getString(self):
        line = 'SineFunction\t'
        line += '%f\t' (self.freq)
        line += '%f' (self.amplitude)
        return line
##################################################
