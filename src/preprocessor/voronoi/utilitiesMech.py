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

        self.__connectedNodes.append(ridge[0])
        self.__connectedNodes.append(ridge[1])
        self.nds = len(self.__connectedNodes)
        #print(self.getString())

    def getString(self):

        line = 'LTCTRSP\t%d'%(self.vertexA)  + '\t' + '%d'%(self.vertexB) +'\t%d'%(self.nds)
        for i in range (self.nds):
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





class transportBC:
    def __init__(self,  nodeIdx, transportBCarray):
        self.transportBCarray = transportBCarray
        self.nodeIdx = nodeIdx
    def getTrsprtBC(self):
        return self.transportBCarray
    def getNodeIdx(self):
        return self.nodeIdx

class mechanicalBC:
    def __init__(self, dim, nodeIdx, mechBCarray):
        self.mechBCarray = mechBCarray
        self.dim = dim
        self.nodeIdx = nodeIdx

    def getDim(self):
        return self.dim
    def getMechBC(self):
        return self.mechBCarray
    def getNodeIdx(self):
        return self.nodeIdx

    def printProps(self):
        print ('Node %d' %(self.nodeIdx))

        if (self.dim ==2):
            print('TransX: %d' %self.mechBCarray[0])
            print('TransY: %d' %self.mechBCarray[1])
            print('RotZ: %d' %self.mechBCarray[2])

        if (self.dim ==3):
            print('TransX: %d' %self.mechBCarray[0])
            print('TransY: %d' %self.mechBCarray[1])
            print('TransZ: %d' %self.mechBCarray[2])
            print('RotX: %d' %self.mechBCarray[3])
            print('RotY: %d' %self.mechBCarray[4])
            print('RotZ: %d' %self.mechBCarray[5])



class mechanicalIC:
    def __init__(self, dim, nodeIdx, mechICArray):
        self.mechICArray = mechICArray
        self.dim = dim
        self.nodeIdx = nodeIdx

    def getNodeIdx(self):
        return self.nodeIdx
    def getMechIC(self):
        return self.mechICArray


class transportIC:
    def __init__(self, vrtxIdx, pressure):
        self.vrtxIdx = vrtxIdx
        self.pressure = pressure

    def getVrtxIdx(self):
        return self.vrtxIdx

    def getPressure(self):
        return self.pressure
