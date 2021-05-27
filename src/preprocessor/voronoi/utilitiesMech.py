import numpy as np
import matplotlib.pyplot as plt
import scipy
import math
from IPython.display import clear_output

class MeasuringGauge:
    def __init__ (self, coordsA, coordsB, name, rotation):
        self.coordsA = coordsA
        self.coordsB = coordsB
        self.name = name
        self.rotation = rotation


class RigidPlate:
    def __init__ (self, govNodeIdx, dim, limits, radial = -1, innerRad = None, directIdcs = False):
        self.govNodeIdx = govNodeIdx
        self.dim = dim
        self.limits = limits
        self.radial = radial
        self.innerRad = innerRad
        # if radial, limits -> center, cylidnerRadius, cylinderThickness
        self.nodes = []
        self.directNodes = []
        self.directIdcs = directIdcs

    def setDirectNodes(self, ndsList):
        for n in ndsList:
            self.directNodes.append(n)

    def renumberVertices(self, dictionary, vertIdxStart):
        #print(vertIdxStart)
        renumberedVertices= []
        for i in range (len(self.directNodes)):
            key = dictionary.get(self.directNodes[i])
            if key:
                ##print('N %d  >  key %d' %(self.directNodes[i], key))
                renumberedVertices.append(dictionary[self.directNodes[i]] + vertIdxStart)



        self.directNodes = renumberedVertices

    def getString (self):
        if self.directIdcs == True:
            line = 'RigidPlate\t%d\t' %self.govNodeIdx
            line += '%d\t' %(len(self.directNodes)) #numSlaves
            for n in self.directNodes:
                line +='%d\t' %n


        if self.radial==0 and self.directIdcs == False:
            ## type      masterNode  numSlaves slave1 slave2 ...
            #RigidPlate  1656 18 18 19 20 21 22 23 24 25 26 27 28 29 30 31 32 33 34 35
            line = 'RigidPlate\t%d\t' %self.govNodeIdx
            line += '%d\t' %(len(self.nodes)) #numSlaves
            for n in self.nodes:
                line +='%d\t' %n

        elif self.directIdcs == False:
            line = 'CoordRigidPlate\t%d\t' %self.govNodeIdx
            line += '%e\t' %self.limits[0] #xmin
            line += '%e\t' %self.limits[1] #xmax
            line += '%e\t' %self.limits[2] #ymin
            line += '%e\t' %self.limits[3] #ymax
            if (self.dim == 3):
                line += '%e\t' %self.limits[4] #zmin
                line += '%e\t' %self.limits[5] #zmax

        return line

    def getNodesAffected (self, nodes):
        nodesAffected = []

        if self.directIdcs == True:
            nodesAffected = self.directNodes

        if self.radial == -1 and self.directIdcs == False:
            for n in range(len(nodes)):
                coords = nodes[n][0:self.dim]
                if self.dim==2:
                    if (coords[0]>self.limits[0] and coords[0]<self.limits[1] and
                        coords[1]>self.limits[2] and coords[1]<self.limits[3]):
                        nodesAffected.append(n)
                if self.dim==3:
                    if (coords[0]>self.limits[0] and coords[0]<self.limits[1] and
                        coords[1]>self.limits[2] and coords[1]<self.limits[3] and
                        coords[2]>self.limits[4] and coords[2]<self.limits[5]):
                        nodesAffected.append(n)

        if self.radial == 0 and self.directIdcs == False:
            for n in range(len(nodes)):
                coords = nodes[n][0:self.dim]
                cylinderThickness = self.limits[4]
                cylinderRad = self.limits[3]
                #print(self.innerRad)
                if self.dim==3:
                    xbounds = np.array([self.limits[0]-cylinderThickness/2, self.limits[0]+cylinderThickness/2])
                    radDist = np.linalg.norm(coords[1:3])
                    if (coords[0]>xbounds[0] and coords[0]<xbounds[1] and  radDist<=cylinderRad):
                        if self.innerRad != None:
                            if radDist>=self.innerRad:
                                nodesAffected.append(n)
                            #else:
                                #print('node bad')
                        else:
                            nodesAffected.append(n)

        self.nodes = nodesAffected


        return nodesAffected



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

    def getString(self, coupled=False):
        ndNr = self.nds
        if not coupled:
            line = 'LTCTRSP\t%d'%(self.vertexA)  + '\t' + '%d'%(self.vertexB) +'\t%d'%(ndNr)
        else:
            line = 'LTCTRSPCoupled\t%d'%(self.vertexA)  + '\t' + '%d'%(self.vertexB) +'\t%d'%(ndNr)

        for i in range (self.nds):
            line+='\t%d'%(self.connectedNodes[i])
        line +='\t' + '%d'%(self.material)
        return line

    def getStringyString(self,nodes, oldAux):
        ndNr = self.nds
        line = 'LTCTRSP\t%d'%(self.vertexA)  + '\t' + '%d'%(self.vertexB) +'\t%d'%(ndNr)
        for i in range (self.nds):
            if (self.connectedNodes[i]<nodes):
                line+=' Node %d' %self.connectedNodes[i]
            elif (self.connectedNodes[i]<nodes+oldAux):
                line+=' OldA %d' %self.connectedNodes[i]
            elif (self.connectedNodes[i]>=nodes+oldAux):
                line+=' NewA %d' %self.connectedNodes[i]

        #line +='%s' %self.connectedNodes

        return line

    def getReducedString(self, coupled=False):
        #oprava
        if (self.connectedNodes[0]!=self.connectedNodes[-1]):
            self.connectedNodes.append(self.connectedNodes[-1])
            self.connectedNodes.append(self.connectedNodes[0])
        self.nds = len(self.connectedNodes)

        ndNr = self.nds /2
        if not coupled:
            line = 'LTCTRSP\t%d'%(self.vertexA)  + '\t' + '%d'%(self.vertexB) +'\t%d'%(ndNr)
        else:
            line = 'LTCTRSPCoupled\t%d'%(self.vertexA)  + '\t' + '%d'%(self.vertexB) +'\t%d'%(ndNr)

        for i in range (0, self.nds, 2):
            line+='\t%d'%(self.connectedNodes[i])

        if (ndNr==2) :
            #print()
            #print (self.connectedNodes)
            #print (line)
            #element to be removed
            line = "#" + line

        #print()
        #print (self.connectedNodes)
        #print (line)

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

#E0	35e9	alpha	0.300000    density 2200.0 fc 200e6 ft 35e6 KinN 4e9 gammaN 20e9 m -0.2e-6 Ad 4000e-6 tauBar 4.0e6 Kin 0.0 gamma 10.0e6 S 0.00025e6 a 0
##################################################
#### Fatigue material shear ####
class FatigueMaterial:
    def __init__ (self, youngModulus, poisson, density, fc, ft, KinN, gammaN, m, Ad, tauBar, Kin, gamma, S, a):
        self.youngModulus = youngModulus
        self.poisson = poisson
        self.density = density
        self.fc = fc
        self.ft = ft
        self.KinN = KinN
        self.gammaN = gammaN
        self.m = m
        self.Ad = Ad
        self.tauBar = tauBar
        self.Kin = Kin
        self.gamma = gamma
        self.S = S
        self.a = a

    def getString (self):
        line = 'FatigueMaterial\t'
        line += 'E0\t%e'        %(self.youngModulus)    + '\t'
        line += 'alpha\t%f'     %(self.poisson)         + '\t'
        line += 'density\t%f'   %(self.density)         + '\t'
        line += 'fc\t%e'        %(self.fc)              + '\t'
        line += 'ft\t%e'        %(self.ft)              + '\t'
        line += 'KinN\t%e'      %(self.KinN)            + '\t'
        line += 'gammaN\t%e'    %(self.gammaN)          + '\t'
        line += 'm\t%f'         %(self.m)               + '\t'
        line += 'Ad\t%f'        %(self.Ad)              + '\t'
        line += 'tauBar\t%e'    %(self.tauBar)          + '\t'
        line += 'Kin\t%e'       %(self.KinN)            + '\t'
        line += 'gamma\t%e'     %(self.gammaN)          + '\t'
        line += 'S\t%e'         %(self.S)               + '\t'
        line += 'a\t%f'         %(self.a)               + '\t'

        return line
##################################################

##################################################
#### Transport material ####
class TransportMaterial:
    def __init__ (self, viscosity, permeability, density, capacity, crack_turtuosity, biot_coeff=0, coupled=False):
        self.viscosity = viscosity
        self.permeability = permeability
        self.density = density
        self.capacity = capacity
        self.crack_turtuosity = crack_turtuosity
        self.coupled = coupled
        self.biot_coeff = biot_coeff

    def getString (self):
        if self.coupled:
            line = 'TrsprtCoupledMaterial'+ '\t' + 'capacity\t%e'%(self.capacity)  + '\t' + 'density\t%e'%(self.density)  + '\t' + 'permeability\t%e'%(self.permeability)  + '\t' + 'viscosity\t%e'%(self.viscosity)  + '\t' +        'crack_turtuosity\t%e'%(self.crack_turtuosity) +   '\t' +        'biot_coeff\t%e'%(self.biot_coeff)
        else:
            line = 'TrsprtMaterial'+ '\t' + 'capacity\t%e'%(self.capacity)  + '\t' + 'density\t%e'%(self.density)  + '\t' + 'permeability\t%e'%(self.permeability)  + '\t' + 'viscosity\t%e'%(self.viscosity)  + '\t' +        'crack_turtuosity\t%e'%(self.crack_turtuosity)

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
