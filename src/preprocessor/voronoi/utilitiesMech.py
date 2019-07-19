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
        self.connectedNodes = connectedNds
        self.nds = len(self.connectedNodes)
        self.material = materialIdx

    def addConnectedNodes (self, ridge):
        #print(ridge)
        #print(self.getString())

        self.connectedNodes.append(ridge[0])
        self.connectedNodes.append(ridge[1])
        self.nds = len(self.connectedNodes)
        #print(self.getString())

    def addSingleConnectedNode (self, idx):
        self.connectedNodes.append(idx)
        self.nds = len(self.connectedNodes)

    def getString(self):
        ndNr = self.nds /2

        line = 'LTCTRSP\t%d'%(self.vertexA)  + '\t' + '%d'%(self.vertexB) +'\t%d'%(ndNr)
        for i in range (self.nds):
            line+='\t%d'%(self.connectedNodes[i])
        line +='\t' + '%d'%(self.material)
        return line

    def getReducedString(self):
        ndNr = self.nds /2
        line = 'LTCTRSP\t%d'%(self.vertexA)  + '\t' + '%d'%(self.vertexB) +'\t%d'%(ndNr)
        for i in range (0, self.nds, 2):
            line+='\t%d'%(self.connectedNodes[i])
        #line+='\t%d'%(self.connectedNodes[len(self.connectedNodes)-1])
        line +='\t' + '%d'%(self.material)
        return line
    def printConnectedNodes(self):
        print (self.connectedNodes)

##################################################


##################################################
#### Linear ealastic materials ####
class linearElasticMaterial:
    def __init__ (self, youngModulus, poisson,  density):
        self.youngModulus = youngModulus
        self.poisson = poisson
        self.density = density

    def getString(self):
        line = 'DisMechMaterial'       + '\t' + 'E0\t%e'%(self.youngModulus)          + '\t' + 'alpha\t%f'%(self.poisson)      + '\t' + 'density\t%f'%(self.density)
        return line

##################################################


##################################################
#### Mars material ####
class MarsMaterial:
    def __init__ (self, youngModulus, poisson, density, ft, Gt):
        self.youngModulus = youngModulus
        self.poisson = poisson
        self.density = density
        self.ft = ft
        self.Gt = Gt
    def getString (self):
        line = 'MarsMaterial\t' +'E0\t%e'%(self.youngModulus)   + '\t' + 'alpha\t%f'%(self.poisson)      + '\t' + 'density\t%f'%(self.density)     +'\t' + 'ft\t%f' %(self.ft) +'\t' + 'Gt\t%f' %(self.Gt)
        return line
##################################################

# FatigueShearMaterial	E0	43.0e9	alpha	0.300000    density 2200.0 tauBar 4.0e6 Kin 0.0 gamma 10.0e6 S 0.0025e6 m 0
##################################################
#### Fatigue material shear ####
class FatigueMaterial:
    def __init__ (self, youngModulus, poisson, density, tauBar, Kin, gamma, S, m):
        self.youngModulus = youngModulus
        self.poisson = poisson
        self.density = density
        self.tauBar = tauBar
        self.Kin = Kin
        self.gamma = gamma
        self.S = S
        self.m = m
    def getString (self):
        line = 'FatigueShearMaterial\t' +'E0\t%e'%(self.youngModulus)   + '\t' + 'alpha\t%f'%(self.poisson)      + '\t' + 'density\t%f'%(self.density)     +'\t' + 'tauBar\t%f' %(self.tauBar) +'\t' + 'Kin\t%f' %(self.Kin) +'\t' + 'gamma\t%f' %(self.gamma) +'\t' + 'S\t%f' %(self.S) +'\t' + 'm\t%f' %(self.m)
        return line
##################################################

##################################################
#### Transport material ####
class TransportMaterial:
    def __init__ (self, transportC, transportS):
        self.transportC = transportC
        self.transportS = transportS

    def getString (self):
        line = 'TrsprtMaterial'+ '\t' + 'capacity\t%f'%(self.transportC)         + '\t' + 'conductivity\t%f'%(self.transportS)
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
