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
class transportPath:
    def __init__ (self,  vertexAidx, vertexBidx, connectedNodes, materialIdx):
        self.vertexA = vertexAidx
        self.vertexB = vertexBidx
        self.connectedNodes = connectedNodes

        self.material = materialIdx

    def getString(self):
        line = 'LTCTRSP\t%d'%(self.vertexA)  + '\t' + '%d'%(self.vertexB) +'\t%d'%(len(self.connectedNodes))
        for i in range (len(self.connectedNodes)):
            line+='\t%d'%(self.connectedNodes[i])
        line +='\t' + '%d'%(self.material)
        return line

##################################################


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
        line = 'DisMechMaterial'       + '\t' + 'E0\t%e'%(self.youngModulus)          + '\t' + 'alpha\t%f'%(self.poisson)      + '\t' + 'density\t%f'%(self.density)     + '\nTrsprtMaterial'         + '\t' + '\t' + 'capacity\t%f'%(self.transportC)         + '\t' + 'conductivity\t%f'%(self.transportS)

        return line

##################################################
