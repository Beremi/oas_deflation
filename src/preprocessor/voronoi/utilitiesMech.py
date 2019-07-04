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
    def __init__ (self,  vertexAidx, vertexBidx, connectedNds, materialIdx):
        self.vertexA = vertexAidx
        self.vertexB = vertexBidx
        self.__connectedNodes = connectedNds
        self.nds = len(self.__connectedNodes)
        self.material = materialIdx

    def addConnectedNodes (self, ridge):
        #print(ridge)
        #print(self.getString())
        """
        if (len(self.__connectedNodes) != self.nds):
            print ('ERROR ERROR ERROR ERROR ERROR')
            print ('orig: %d' %(self.nds))
            print ('curr: %d' %(len(self.__connectedNodes)))
        else:
            print('OK')
        """
        self.__connectedNodes.append(ridge[0])
        self.__connectedNodes.append(ridge[1])
        self.nds = len(self.__connectedNodes)
        #print(self.getString())

    def getString(self):
        line = 'LTCTRSP\t%d'%(self.vertexA)  + '\t' + '%d'%(self.vertexB) +'\t%d'%(len(self.__connectedNodes))
        for i in range (len(self.__connectedNodes)):
            line+='\t%d'%(self.__connectedNodes[i])
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
