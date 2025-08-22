import numpy as np
import random
import utilitiesGeom
import utilitiesMech
import utilitiesNumeric
import pointGenerators
import voronoi
import math
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
#matplotlib.matplotlib_fname()
#import voronoi_viewer


from scipy.spatial import Voronoi
from scipy.spatial import voronoi_plot_2d
from scipy.spatial import Delaunay
#import tkinter
SHOW_PLOT = False


def assembleMeasuringGauges(type, D=-1, thickness = 0.1, maxLim = None, expansionRingsProps=[], notch=[],props=[]):
    measuringGauges = []
    if (type == 'dogbone2d'):
        #if (D==0.1):
        #total length LS
        coordsA = np.array([ D/2, 0])
        coordsB = np.array([ D/2, 6/4*D ])
        measuringGauges.append(utilitiesMech.MeasuringGauge(coordsA, coordsB, 'totalLS', False))
        #mid LS
        coordsA = np.array([ D/2, 3/4*D-D*0.6/2 ])
        coordsB = np.array([ D/2, 3/4*D+D*0.6/2 ])
        measuringGauges.append(utilitiesMech.MeasuringGauge(coordsA, coordsB, 'midLS', False))
        #left LC
        coordsA = np.array([ 0.2*D, 3/4*D-D*0.6/2 ])
        coordsB = np.array([ 0.2*D, 3/4*D+D*0.6/2 ])
        measuringGauges.append(utilitiesMech.MeasuringGauge(coordsA, coordsB, 'leftLS', False))
        #right LC
        coordsA = np.array([ D-0.2*D, 3/4*D-D*0.6/2 ])
        coordsB = np.array([ D-0.2*D, 3/4*D+D*0.6/2 ])
        measuringGauges.append(utilitiesMech.MeasuringGauge(coordsA, coordsB, 'rightLS', False))

    if (type == 'dogbone3d'):
        #total length LS front
        coordsA = np.array([ D/2, 0, thickness])
        coordsB = np.array([ D/2, 6/4*D, thickness ])
        measuringGauges.append(utilitiesMech.MeasuringGauge(coordsA, coordsB, 'totalLS_front', False))
        #mid LS front
        coordsA = np.array([ D/2, 3/4*D-D*0.6/2, thickness  ])
        coordsB = np.array([ D/2, 3/4*D+D*0.6/2, thickness  ])
        measuringGauges.append(utilitiesMech.MeasuringGauge(coordsA, coordsB, 'midLS_front', False))
        #left LC front
        coordsA = np.array([ 0.2*D, 3/4*D-D*0.6/2, thickness ])
        coordsB = np.array([ 0.2*D, 3/4*D+D*0.6/2, thickness ])
        measuringGauges.append(utilitiesMech.MeasuringGauge(coordsA, coordsB, 'leftLS_front', False))
        #right LC front
        coordsA = np.array([ D-0.2*D, 3/4*D-D*0.6/2, thickness ])
        coordsB = np.array([ D-0.2*D, 3/4*D+D*0.6/2, thickness ])
        measuringGauges.append(utilitiesMech.MeasuringGauge(coordsA, coordsB, 'rightLS_front', False))

        #total length LS back
        coordsA = np.array([ D/2, 0, 0])
        coordsB = np.array([ D/2, 6/4*D, 0])
        measuringGauges.append(utilitiesMech.MeasuringGauge(coordsA, coordsB, 'totalLS_back', False))
        #mid LS back
        coordsA = np.array([ D/2, 3/4*D-D*0.6/2, 0  ])
        coordsB = np.array([ D/2, 3/4*D+D*0.6/2, 0  ])
        measuringGauges.append(utilitiesMech.MeasuringGauge(coordsA, coordsB, 'midLS_back', False))
        #left LC back
        coordsA = np.array([ 0.2*D, 3/4*D-D*0.6/2, 0])
        coordsB = np.array([ 0.2*D, 3/4*D+D*0.6/2, 0])
        measuringGauges.append(utilitiesMech.MeasuringGauge(coordsA, coordsB, 'leftLS_back', False))
        #right LC back
        coordsA = np.array([ D-0.2*D, 3/4*D-D*0.6/2, 0 ])
        coordsB = np.array([ D-0.2*D, 3/4*D+D*0.6/2, 0 ])
        measuringGauges.append(utilitiesMech.MeasuringGauge(coordsA, coordsB, 'rightLS_back', False))

        #total length inside speciment
        coordsA = np.array([D - 0.2 * D, 3 / 4 * D - D * 0.6 / 2, thickness / 2])
        coordsB = np.array([D - 0.2 * D, 3 / 4 * D + D * 0.6 / 2, thickness / 2])
        measuringGauges.append(utilitiesMech.MeasuringGauge(coordsA, coordsB, 'total_LS_inside', False))

    if (type == 'dogbone2dStrip'):
        #if (D==0.1):
        #total length LS
        coordsA = np.array([ D/2, 0])
        coordsB = np.array([ D/2, 1 ])
        measuringGauges.append(utilitiesMech.MeasuringGauge(coordsA, coordsB, 'totalLS', False))

    if (type == 'dogbone2dBand'):
        #if (D==0.1):
        #total length LS
        coordsA = np.array([ maxLim[0]/2, 0])
        coordsB = np.array([ maxLim[0]/2, maxLim[1]])
        measuringGauges.append(utilitiesMech.MeasuringGauge(coordsA, coordsB, 'totalLS', False))

    if (type == 'coupledPress2d'):
        #if (D==0.1):
        #total length LS
        coordsA = np.array([ maxLim[0]/2, 0])
        coordsB = np.array([ maxLim[0]/2, maxLim[1] ])
        measuringGauges.append(utilitiesMech.MeasuringGauge(coordsA, coordsB, 'total', False))


    if(type=='2d_corrosionRebar'):
        rebarCount = expansionRingsProps[0]
        rebarDepth = expansionRingsProps[1]
        rebarDiameter = expansionRingsProps[2]
        maxLim = expansionRingsProps[3]

        for i in range (rebarCount):
            if (rebarCount==1):
                #puvodni poloha rebars polovina od kraje
                centre = np.array([ (maxLim[0]/rebarCount)*(i+0.5), maxLim[1]-rebarDepth  ])
            else:
                #poloha rebars presne jak je ve clanku
                centre = np.array([ (0.058 + (maxLim[0]-0.116)/(rebarCount-1)*i), maxLim[1]-rebarDepth  ])

            coordsA = np.array([ centre[0]-0.005, maxLim[1] ])
            coordsB = np.array([ centre[0]+0.005, maxLim[1] ])
            name = 'surf-Rebar#%d'%i
            measuringGauges.append(utilitiesMech.MeasuringGauge(coordsA, coordsB, name, False))

            coordsA = np.array([ centre[0]-rebarDiameter/2, centre[1] ])
            coordsB = np.array([ centre[0]+rebarDiameter/2, centre[1] ])
            name = 'horiz-Rebar#%d'%i
            measuringGauges.append(utilitiesMech.MeasuringGauge(coordsA, coordsB, name, False))

            coordsA = np.array([ centre[0], centre[1]-rebarDiameter/2 ])
            coordsB = np.array([ centre[0], centre[1]+rebarDiameter/2 ])
            name = 'vert-Rebar#%d'%i
            measuringGauges.append(utilitiesMech.MeasuringGauge(coordsA, coordsB, name, False))

    if(type=='reinhardt3d'):
        #total length
        coordsA = np.array([ 0, maxLim[1]/2, maxLim[2]/2])
        coordsB = np.array([ maxLim[0], maxLim[1]/2, maxLim[2]/2 ])
        measuringGauges.append(utilitiesMech.MeasuringGauge(coordsA, coordsB, 'total', False))


    if(type=='2d_singleSpring'):
        #total length
        coordsA = np.array([ 0, 0 ])
        coordsB = np.array([ maxLim[0], 0 ])
        measuringGauges.append(utilitiesMech.MeasuringGauge(coordsA, coordsB, 'total', True))

    if(type=='2d_RWTH'):
        #total length
        coordsA = np.array([ 0, 0 ])
        coordsB = np.array([ maxLim[0]/2,  maxLim[1] ])
        measuringGauges.append(utilitiesMech.MeasuringGauge(coordsA, coordsB, 'total', True))

    if(type=='cylinder3d'):
        #total length
        coordsA = np.array([ 0, 0, 0])
        coordsB = np.array([ maxLim[0], 0, 0 ])
        measuringGauges.append(utilitiesMech.MeasuringGauge(coordsA, coordsB, 'total', True))

    if(type=='tubesplit'):
        #total length
        coordsA = np.array([ 0, 0, 0])
        coordsB = np.array([ maxLim[0], 0, 0 ])
        measuringGauges.append(utilitiesMech.MeasuringGauge(coordsA, coordsB, 'totalX', True))
        coordsA = np.array([ maxLim[0]/2, -maxLim[1], 0])
        coordsB = np.array([ maxLim[0]/2, maxLim[1], 0 ])
        measuringGauges.append(utilitiesMech.MeasuringGauge(coordsA, coordsB, 'totalY', True))
        coordsA = np.array([ maxLim[0]/2, 0,-maxLim[2]])
        coordsB = np.array([ maxLim[0]/2, 0, maxLim[2]])
        measuringGauges.append(utilitiesMech.MeasuringGauge(coordsA, coordsB, 'totalZ', True))
        rad = props[0]
        thickness = props[1]
        width = props[2]
        coordsA = np.array([ maxLim[0]/2, width/2, -(rad-thickness) ])
        coordsB = np.array([ maxLim[0]/2, -width/2, -(rad-thickness) ])
        measuringGauges.append(utilitiesMech.MeasuringGauge(coordsA, coordsB, 'leftnotch', True))
        coordsA = np.array([ maxLim[0]/2, width/2, (rad-thickness) ])
        coordsB = np.array([ maxLim[0]/2, -width/2, (rad-thickness) ])
        measuringGauges.append(utilitiesMech.MeasuringGauge(coordsA, coordsB, 'rightnotch', True))

    if(type=='3pb2d'):
        coordsA = np.array([ maxLim[0]/2, maxLim[1]])
        coordsB = np.array([ 0, 0 ])
        measuringGauges.append(utilitiesMech.MeasuringGauge(coordsA, coordsB, 'totalLS', False))
        coordsA = np.array([ maxLim[0]/2-1e3, 0])
        coordsA = np.array([ maxLim[0]/2+1e3, 0])
        measuringGauges.append(utilitiesMech.MeasuringGauge(coordsA, coordsB, 'notchOpening', False))

    if(type=='box3d'):
        coordsA = np.array([ 0, maxLim[1]/2, maxLim[2]/2])
        coordsB = np.array([ maxLim[0], maxLim[1]/2, maxLim[2]/2])
        measuringGauges.append(utilitiesMech.MeasuringGauge(coordsA, coordsB, 'contrXmid', False))
        coordsA = np.array([ maxLim[0]/2, 0, maxLim[2]/2])
        coordsB = np.array([ maxLim[0]/2, maxLim[1], maxLim[2]/2])
        measuringGauges.append(utilitiesMech.MeasuringGauge(coordsA, coordsB, 'contrYmid', False))
        coordsA = np.array([ maxLim[0]/2, maxLim[1]/2, 0])
        coordsB = np.array([ maxLim[0]/2, maxLim[1]/2, maxLim[2]])
        measuringGauges.append(utilitiesMech.MeasuringGauge(coordsA, coordsB, 'contrZmid', False))




    if(type=='3pb3d'):
        coordsA = np.array([ maxLim[0]/2, maxLim[1], maxLim[2]/2])
        coordsB = np.array([ 0, 0, 0 ])
        measuringGauges.append(utilitiesMech.MeasuringGauge(coordsA, coordsB, 'totalTopMid', False))
        coordsA = np.array([ maxLim[0]/2, 0, maxLim[2]/2])
        coordsB = np.array([ 0, 0, 0 ])
        measuringGauges.append(utilitiesMech.MeasuringGauge(coordsA, coordsB, 'totalBotMid', False))
        coordsA = np.array([ maxLim[0]/2-maxLim[0]/3, 0, maxLim[2]/2])
        coordsB = np.array([ maxLim[0]/2+maxLim[0]/3, 0, maxLim[2]/2])
        measuringGauges.append(utilitiesMech.MeasuringGauge(coordsA, coordsB, 'notchOpening', False))
        if len(notch)>0:
            notchHeight=notch[0]
            notchWidth=notch[1]
            coordsA = np.array([ maxLim[0]/2-notchWidth, maxLim[1]*notchHeight, maxLim[2]/2])
            coordsB = np.array([ maxLim[0]/2+notchWidth, maxLim[1]*notchHeight, maxLim[2]/2])
            measuringGauges.append(utilitiesMech.MeasuringGauge(coordsA, coordsB, 'CTOD', False))

    if(type=='3d_brazilianDisc'):
        #total length
        coordsA = np.array([ maxLim[0]/2, maxLim[1]/2, 0])
        coordsB = np.array([ maxLim[0]/2, -maxLim[1]/2, 0 ])
        measuringGauges.append(utilitiesMech.MeasuringGauge(coordsA, coordsB, 'totalLoadDir', False))
        coordsA = np.array([ maxLim[0]/2, 0, maxLim[2]/2])
        coordsB = np.array([ maxLim[0]/2, 0, -maxLim[2]/2 ])
        measuringGauges.append(utilitiesMech.MeasuringGauge(coordsA, coordsB, 'totalAcrossDir', False))
        coordsA = np.array([ maxLim[0]/2, 0, 0.0508/2])
        coordsB = np.array([ maxLim[0]/2, 0, -0.0508/2 ])
        measuringGauges.append(utilitiesMech.MeasuringGauge(coordsA, coordsB, 'crackAcrossDir', False))

    return measuringGauges
    """
    if (type == 'dogboneD_2d'):
        #gauges D
        if (D==0.4):
            govNodes.append( np.array([ D/2,        3/4*D-0.240/2 ])  )#mid LS
            govNodes.append( np.array([ D/2,        3/4*D+0.240/2 ])  )
            govNodes.append( np.array([ D/2-0.2*D,  3/4*D-0.240/2 ])  )#left LC
            govNodes.append( np.array([ D/2-0.2*D,  3/4*D+0.240/2 ])  )
            govNodes.append( np.array([ D/2+0.2*D,  3/4*D-0.240/2 ])  )#right LC
            govNodes.append( np.array([ D/2+0.2*D,  3/4*D+0.240/2 ])  )
            for i in range (6):
                govNodesMechBC.append(utilitiesMech.mechanicalBC(dim, -1, blankMechBC))
    """


def assembleMaterialZones (elaX, dim, model='box', maxLim=None, D=None, thickness=None, limits=None, limits1=None, rebarDepth=None, rebarDiameter=None, rebarCount=None, minDist=None, notch=None, weakboundary = False):
    #limits = xmin, ymin, zmin, xmax, ymax, zmax
    materialZones = []
    #matZone 1
    matZ = []
    #matZone 2
    matZ1 = []

    matZ2 = []
    if (model=='tdcb'):
        if (dim==2):
            boundA = np.array(  [ -1e-8             , -1e-8          ] )
            matZ.append (boundA)
            boundB = np.array(  [ elaX    , maxLim[1] + 1e-8] )
            matZ.append (boundB)
            boundA1 = np.array(  [ maxLim[0]-elaX , - 1e-8] )
            matZ.append (boundA1)
            boundB1 = np.array(  [ maxLim[0] + 1e-8 , maxLim[1] + 1e8]  )
            matZ.append (boundB1)
            materialZones.append(matZ)


            boundA2 = np.array(  [ -1e-8    ,  -1e-8] )
            boundB2 = np.array(  [ maxLim[0]+1e-8   ,  maxLim[1]*0.37 ])
            matZ1.append(boundA2)
            matZ1.append(boundB2)

            boundA3 = np.array(  [ -1e-8    ,  maxLim[1]*0.79] )
            boundB3 = np.array(  [ maxLim[0]+1e-8   ,  maxLim[1]*1.1 ])
            matZ1.append(boundA3)
            matZ1.append(boundB3)
            materialZones.append(matZ1)

        if (dim==3):
            boundA = np.array(  [ maxLim[0]/2-elaX             , maxLim[1]*0.37*0.95   , -1e-8] )
            matZ.append (boundA)
            boundB = np.array(  [ maxLim[0]/2+elaX    , maxLim[1]*0.8, maxLim[2]+1e-8] )
            matZ.append (boundB)
            boundA1 = np.array(  [ maxLim[0]/2-elaX             , maxLim[1]*0.37*0.95   , -1e-8 ] )
            matZ.append (boundA1)
            boundB1 = np.array(  [ maxLim[0]/2+elaX    , maxLim[1]*0.8, maxLim[2]+1e-8] )
            matZ.append (boundB1)
            materialZones.append(matZ)

    if (model=='tdcbfem'):
        if (dim==3):
            boundA = np.array(  [ maxLim[0]+1e-5             , maxLim[2]-1e-5     , -1e-8] )
            matZ.append (boundA)
            boundB = np.array(  [ maxLim[1]-1e-5    , maxLim[3]-1e-5  , 1e8] )
            matZ.append (boundB)
            boundA1 = np.array(  [ maxLim[0]+1e-5             , maxLim[2]-1e-5     , -1e-8] )
            matZ.append (boundA1)
            boundB1 = np.array(  [ maxLim[1]-1e-5    , maxLim[3]-1e-5  , 1e8] )
            matZ.append (boundB1)
            materialZones.append(matZ)

    if (model=='hangingfraczone'):
        if (dim==2):
            boundA = np.array(  [ maxLim[0]+10e-5             , maxLim[2]-10e-5 ] )
            matZ.append (boundA)
            boundB = np.array(  [ maxLim[1]-10e-5    , maxLim[3]-10e-5 ] )
            matZ.append (boundB)
            boundA1 = np.array(  [ maxLim[0]+10e-5             , maxLim[2]-10e-5    ] )
            matZ.append (boundA1)
            boundB1 = np.array(  [ maxLim[1]-10-5    , maxLim[3]-10e-5 ] )
            matZ.append (boundB1)
            materialZones.append(matZ)

        if (dim==3):
            boundA = np.array(  [ maxLim[0]+1e-5             , maxLim[2]-1e-5     , -1e-8] )
            matZ.append (boundA)
            boundB = np.array(  [ maxLim[1]-1e-5    , maxLim[3]-1e-5  , 1e8] )
            matZ.append (boundB)
            boundA1 = np.array(  [ maxLim[0]+1e-5             , maxLim[2]-1e-5     , -1e-8] )
            matZ.append (boundA1)
            boundB1 = np.array(  [ maxLim[1]-1e-5    , maxLim[3]-1e-5  , 1e8] )
            matZ.append (boundB1)
            materialZones.append(matZ)



    if (model=='box'):
        if (dim==2):
            boundA = np.array(  [ -1e-8             , -1e-8          ] )
            matZ.append (boundA)
            boundB = np.array(  [ elaX    , maxLim[1] + 1e-8] )
            matZ.append (boundB)
            boundA1 = np.array(  [ maxLim[0]-elaX , - 1e-8] )
            matZ.append (boundA1)
            boundB1 = np.array(  [ maxLim[0] + 1e-8 , maxLim[1] + 1e8]  )
            matZ.append (boundB1)
            materialZones.append(matZ)
        if (dim==3):
            boundA = np.array(  [ -1e-8   -maxLim[0]          , -1e-8    -maxLim[1]         , -1e8 -maxLim[2]] )
            matZ.append (boundA)
            boundB = np.array(  [ elaX    , maxLim[1] + 1e8   , maxLim[2] + 1e8  ] )
            matZ.append (boundB)
            boundA1 = np.array(  [ maxLim[0]-elaX , - 1e-8  -maxLim[1]      , -1e8-maxLim[2]] )
            matZ.append (boundA1)
            boundB1 = np.array(  [ maxLim[0] + 1e-8 , maxLim[1] + 1e8   , maxLim[2] + 1e8 ]  )
            matZ.append (boundB1)
            materialZones.append(matZ)

    if (model=='dogbone'):
        if (dim==2):

            if weakboundary:
                radius = 0.725 * D
                b_radius = 0.725 * D + elaX

                # left boundary
                center = np.array([0.8 * D + radius, 3/4 * D])
                matZ = []
                matZ.append('wb')
                matZ.append(b_radius)
                matZ.append(center)
                materialZones.append(matZ)

                # right boundary
                center = np.array([0.2 * D - radius, 3/4 * D])
                matZ = []
                matZ.append('wb')
                matZ.append(b_radius)
                matZ.append(center)
                materialZones.append(matZ)

            else:
                boundA = np.array(  [ -1e-8    , -1e-8  ] )
                matZ.append (boundA)
                boundB = np.array(  [ D+1e-8   ,  elaX] )
                matZ.append (boundB)
                boundA1 = np.array(  [ -1e-8, 6/4*D - elaX] )
                matZ.append (boundA1)
                boundB1 = np.array(  [ D  ,  6/4*D+1e-8]  )
                matZ.append (boundB1)
                materialZones.append(matZ)

        if (dim==3):
            boundA = np.array(  [ -1e-8    , -1e-8, -1e-8  ] )
            matZ.append (boundA)
            boundB = np.array(  [ D+1e-8   ,  elaX, thickness+1e-8] )
            matZ.append (boundB)
            boundA1 = np.array(  [ -1e-8   , 6/4*D - elaX, -1e-8]  )
            matZ.append (boundA1)
            boundB1 = np.array(  [ 2*D, 2*D, thickness*2] )
            matZ.append (boundB1)
            materialZones.append(matZ)

    if (model=='dogboneStrip'):
        if (dim==2):
            boundA = np.array(  [ -1e-8    , -1e-8  ] )
            matZ.append (boundA)
            boundB = np.array(  [ D+1e-8   ,  3/4 * D + 1e-4  ] )
            matZ.append (boundB)
            boundA1 = np.array(  [ -1e-8, 6/4*D - 3/4 * D - 1e-4 ] )
            matZ.append (boundA1)
            boundB1 = np.array(  [ D  ,  6/4*D+1e-8]  )
            matZ.append (boundB1)
            materialZones.append(matZ)

    if (model=='dogboneBand'):
        if (dim==2):
            elasticHeight = maxLim[0]
            dmgHeight = maxLim[1]
            width = maxLim[2]


            boundA = np.array(  [ -1e-8    , -1e-8  ] )
            matZ.append (boundA)
            boundB = np.array(  [ width+1e-8   ,  elasticHeight+ 1e-8  ] )
            matZ.append (boundB)

            boundA1 = np.array(  [ -1e-8, elasticHeight+dmgHeight ] )
            matZ.append (boundA1)
            boundB1 = np.array(  [ width+1e-8  ,  elasticHeight+dmgHeight+elasticHeight+1e-8  ]  )
            matZ.append (boundB1)
            materialZones.append(matZ)


    # 3/4 * D - indent -minDist/2]

    if (model =='4pb2d'):

        boundA = np.array(  [ -1    ,  -1   ] )
        matZ.append (boundA)
        boundB = np.array(  [  maxLim[0]/5   ,  maxLim[1]/4 ] )
        matZ.append (boundB)


        boundA = np.array(  [ maxLim[0]*3/5    ,  -1   ] )
        matZ.append (boundA)
        boundB = np.array(  [  maxLim[0]+1   ,  maxLim[1]/4 ] )
        matZ.append (boundB)
        materialZones.append(matZ)


        boundA2 = np.array(  [ maxLim[0]/3-maxLim[0]/15     ,  maxLim[1]-maxLim[0]/15] )
        boundB2 = np.array(  [ maxLim[0]/3+maxLim[0]/15     ,  maxLim[1]+1 ])
        matZ1.append(boundA2)
        matZ1.append(boundB2)
        #materialZones.append(matZ2)

        boundA3 = np.array(  [ maxLim[0]/3*2-maxLim[0]/15     ,  maxLim[1]-maxLim[0]/15] )
        boundB3 = np.array(  [ maxLim[0]/3*2+maxLim[0]/15    ,  maxLim[1]+1] )
        matZ1.append(boundA3)
        matZ1.append(boundB3)
        materialZones.append(matZ1)


    if (model =='3pb3dX'):
        if limits is not None:
            boundA = np.array(  [ limits[0]    ,  limits[1],  limits[2]  ] )
            matZ.append (boundA)
            boundB = np.array(  [  limits[3]   ,   limits[4],  limits[5] ] )
            matZ.append (boundB)
            #materialZones.append(matZ)

            boundA1 = np.array(  [ limits1[0]    ,  limits1[1],  limits1[2]  ] )
            boundB1 = np.array(  [ limits1[3]    ,  limits1[4],  limits1[5]  ] )
            matZ1.append(boundA1)
            matZ1.append(boundB1)


        boundA2 = np.array(  [ -100    ,  maxLim[1]-1e-3,  -100  ] )
        boundB2 = np.array(  [ 100   ,  100,  100  ] )
        matZ2.append(boundA2)
        matZ2.append(boundB2)
        #matZ2.append(boundA2)
        #matZ2.append(boundB2)
        materialZones.append(matZ2)

    if (model=='3dcylinder'):
        boundA = np.array(  [ -1e-8   -maxLim[0]          , -1e-8    -maxLim[1]         , -1e8 -maxLim[2]] )
        matZ.append (boundA)
        boundB = np.array(  [ elaX*maxLim[0]     , maxLim[1] + 1e8   , maxLim[2] + 1e8  ] )
        matZ.append (boundB)

        boundA1 = np.array(  [ maxLim[0]-elaX*maxLim[0]  , - 1e-8  -maxLim[1]      , -1e8-maxLim[2]] )
        matZ.append (boundA1)
        boundB1 = np.array(  [ maxLim[0] + 1e-8 , maxLim[1] + 1e8   , maxLim[2] + 1e8 ]  )
        matZ.append (boundB1)
        materialZones.append(matZ)

    #centre = np.array([ (maxLim[0]/rebarCount)*(r+0.5), maxLim[1]-rebarDepth  ])
    if (model=='2d_corrosionRebar'):
        for r in range(rebarCount):
            radius = rebarDiameter/2 + 0.0002
            if (rebarCount==1):
                #puvodni poloha rebars polovina od kraje
                center = np.array([ (maxLim[0]/rebarCount)*(r+0.5), maxLim[1]-rebarDepth  ])
            else:
                #poloha rebars presne jak je ve clanku
                center = np.array([ (0.058 + (maxLim[0]-0.116)/(rebarCount-1)*r), maxLim[1]-rebarDepth  ])
            matZ = []
            matZ.append('circle')
            matZ.append(radius)
            matZ.append(center)
            matZ.append(rebarCount)
            materialZones.append(matZ)

    if (model=='3d_corrosionRebar'):
        for r in range(rebarCount):
            radius = rebarDiameter/2 + 0.0002
            if (rebarCount==1):
                #puvodni poloha rebars polovina od kraje
                center = np.array([ (maxLim[0]/rebarCount)*(r+0.5), maxLim[1]-rebarDepth  ])
            else:
                #poloha rebars presne jak je ve clanku
                center = np.array([ (0.058 + (maxLim[0]-0.116)/(rebarCount-1)*r), maxLim[1]-rebarDepth  ])
            matZ = []
            matZ.append('circle')
            matZ.append(radius)
            matZ.append(center)
            matZ.append(rebarCount)
            materialZones.append(matZ)


    return materialZones




def createSingleSpringTestModel(length, master_folder):
    print ('Creating single spring test model.')
    #defining functions
    functions = []
    #### Defining functions
    #0 sine func
    #fn = utilitiesNumeric.sineFunc(10,22)

    #constant zero
    fn = utilitiesNumeric.constantFunc(0)
    functions.append (fn)


    fn2 = utilitiesNumeric.sawToothConstFunc(value = 1, period = 1, sym =1)
    functions.append (fn2)

    """
    fn4 = utilitiesNumeric.varyingSawToothFunction(fn1, fn2)
    functions.append(fn4)
    """
    dim = 2
    idt = length /2
    maxLim = np.array([  length*2    ,   2  ])

    node_coords, mechBC_merged,  govNodes, govNodesMechBC, rigidPlates = assembleTwoNodeSpringTest(maxLim, idt)

    print('Conducting Voronoi tesselation...', end = '')
    vor, regions, vertices, polygons, areas, centroids, points = utilitiesNumeric.runMirroredVoronoi (node_coords, 2, maxLim)
    print('done.')

    # if SHOW_PLOT:
    #     fig = voronoi.voronoi_plot_2d(vor, show_vertices = True)
    #     plt.show()

    #print(node_coords)
    node_coords = np.asarray(node_coords)
    # if SHOW_PLOT:
    #     fig, ax = plt.subplots()
    #     ax.scatter(node_coords[:,0], node_coords[:,1])
    #     plt.show()

    return node_coords, mechBC_merged,  vor, areas, functions,govNodes, govNodesMechBC, rigidPlates


def create_2d_SingleContat(length, master_folder,maxLim):
    print ('Creating 2d singlecontact model.')
    functions = []

    #constant zero
    fn = utilitiesNumeric.constantFunc(0)
    functions.append (fn)


    fn2 = utilitiesNumeric.sawToothConstFunc(value = 1, period = 1, sym =1)
    functions.append (fn2)

    dim = len(maxLim)
    idt = length /2


    node_coords, mechBC_merged,  govNodes, govNodesMechBC, rigidPlates = assemble2D_singleContact(maxLim, idt)

    print('Conducting Voronoi tesselation...', end = '')
    vor, regions, vertices, polygons, areas, centroids, points = utilitiesNumeric.runMirroredVoronoi (node_coords, 2, maxLim)
    print('done.')

    node_coords = np.asarray(node_coords)

    return node_coords, mechBC_merged,  vor, areas, functions,[], [], rigidPlates

def create_3d_SingleContat(length, master_folder,maxLim):
    print ('Creating 3d singlecontact model.')
    functions = []

    #constant zero
    fn = utilitiesNumeric.constantFunc(0)
    functions.append (fn)


    fn2 = utilitiesNumeric.sawToothConstFunc(value = 1, period = 1, sym =1)
    functions.append (fn2)

    dim = len(maxLim)
    idt = length /2


    node_coords, mechBC_merged,  govNodes, govNodesMechBC, rigidPlates = assemble3D_singleContact(maxLim, idt)
    #print(node_coords)
    print('Conducting Voronoi tesselation...', end = '')
    vor, areas = utilitiesNumeric.runMirroredVoronoi (node_coords, dim, maxLim)
    print('done.')

    node_coords = np.asarray(node_coords)
    #print(mechBC_merged)

    return node_coords, mechBC_merged,  vor, areas, functions,[], [], rigidPlates


def createDiamondTestModel(width, height):
    print('Creating diamond test model.')
    #defining functions
    functions = []
    #### Defining functions

    #constant zero
    fn = utilitiesNumeric.constantFunc(0)
    functions.append (fn)
    #constant load
    fn1 = utilitiesNumeric.constantFunc(-1e-3)
    functions.append (fn1)

    dim = 2
    idtW = 1e-10
    idtH = 1e-10

    maxLim = np.array([  width   ,   height ])
    #shifts = -maxLim / 2
    shifts = np.zeros(2)

    node_coords,  mechBC_merged = assembleDiamondTest(maxLim, idtW, idtH)

    #print(node_coords)

    print('Conducting Voronoi tesselation...', end = '')
    vor, regions, vertices, polygons, areas, centroids, points = utilitiesNumeric.runMirroredVoronoi (node_coords, 2, maxLim, shifts=shifts)
    print('done.')

    #print(vor.points)
    #print(areas)

    if SHOW_PLOT:
        fig = voronoi.voronoi_plot_2d(vor, show_vertices = True)
        plt.show()

    #print (points)

    idt = 1e-3
    ### indirect setting of transportBCs by spatial selection of vertices
    transportBC_merged = []
    ### selecting vertices in model
    noTrsprtBC = np.array([ -1 , -1 ])
    boundA = np.array(  [ shifts[0] - idt , shifts[1] -idt ] )
    boundB = np.array(  [ maxLim[0] +shifts[0] + idt , maxLim[1] +shifts[0]  + idt  ]  )
    allVrtcs = utilitiesGeom.returnSelectedPts(boundA, boundB, vor.vertices)

    #print(allVrtcs)
    for i in range (len(allVrtcs)):
        trsBC = utilitiesMech.transportBC(allVrtcs[i], noTrsprtBC)
        transportBC_merged.append(trsBC)


    return node_coords, mechBC_merged, transportBC_merged, vor, areas, functions



def create2dSSBeamUnifLoad(maxLim, minDist, trials, notch = -1,      loadWidth = 1, fracZoneWidth = 0.15,                           orthogonalFracZone=False, notchWidth =-1,                           node_coords_init=None,                                   activeTransport=False,                           coupled = False, specifiedNodes=[], loading="3pb", gapWidth=0):
    print('Creating 2d simply supported beam, uniform load.')

    #notchWidth = 1.5e-3 /2
    if notchWidth == -1:
        notchWidth = minDist/2
    else:
        notchWidth /= 2

    node_coords, mechBC_merged, mechInitC_merged, notches, govNodes, govNodesMechBC, rigidPlates  = assemble2DSSBeamBending(maxLim, minDist, trials, notch, loadWidth, fracZoneWidth, orthogonalFracZone=orthogonalFracZone, notchWidth=notchWidth, node_coords_init=node_coords_init,  coupled=coupled, specifiedNodes=specifiedNodes, loading=loading, gapWidth=gapWidth);

    if notch == 0:
        notch = None
    else:
        notch=[notch,notchWidth]

    print('Conducting Voronoi tesselation...', end = '')
    vor, regions, vertices, polygons, areas, centroids, points = utilitiesNumeric.runMirroredVoronoi (node_coords, 2, maxLim,notch=notch)
    print('done.')


    ########################################################################
    functions = []
    #### Defining functions
    #0 constant zero
    fn = utilitiesNumeric.constantFunc(0)
    functions.append (fn)

    #1 loading function
    func1 = []
    func1.append( np.array([0,0]) )
    func1.append( np.array([0.1,0]) )
    func1.append( np.array([1, -1e-3]) )
    fn1 = utilitiesNumeric.generalFunc(func1)
    functions.append (fn1)

    if gapWidth>0:
        func1 = []
        func1.append( np.array([0,0]) )
        func1.append( np.array([0.1,100]) )
        func1.append( np.array([1, 100]) )
        fn1 = utilitiesNumeric.generalFunc(func1)
        functions.append (fn1)


    ########################################################################
    ### indirect setting of transportBCs by spatial selection of vertices
    transportBC_merged = []


    #"""
    ### selecting vertices on the bottom surface
    boundA = np.zeros(2)-1e-8
    boundB = np.array([maxLim[0], 1e-8])
    faces1 = utilitiesGeom.returnSelectedPts(boundA, boundB, vor.vertices)
    vert = vor.vertices[faces1,:]

    fn1 = utilitiesNumeric.constantFunc(100)
    functions.append (fn1)

    for i,k in enumerate(faces1):
        trsBC = utilitiesMech.transportBC(k,[2,-1])
        transportBC_merged.append(trsBC)

    #"""
    ### selecting vertices on the top surface
    boundA = np.array([-1e-8, maxLim[1]-1e-8])
    boundB = np.array([maxLim[0], maxLim[1]+1e-8])
    faces1 = utilitiesGeom.returnSelectedPts(boundA, boundB, vor.vertices)
    vert = vor.vertices[faces1,:]

    fn1 = utilitiesNumeric.constantFunc(0)
    functions.append (fn1)

    for i,k in enumerate(faces1):
        trsBC = utilitiesMech.transportBC(k,[3,-1])
        transportBC_merged.append(trsBC)

    #"""

    #print(transportBC_merged)


    #return node_coords, mechBC_merged, transportBC_merged,  notches, govNodes, govNodesMechBC, rigidPlates, vor, areas, functions

    return node_coords, mechBC_merged, [], vor, areas, functions, notches, govNodes,    govNodesMechBC, rigidPlates, transportBC_merged




def create2dCantileverBending(maxLim, minDist, trials, powerTes, nodefile = None ):

    print('Creating 2d cantilever, bending.')
    ### sampling of nodes
    ### direct setting of mechanicalBCs
    node_coords, radii, mechBC_merged, mechIC_merged  = assemble2DCantileverBending(maxLim, minDist, trials, powerTes, nodefile = nodefile );

    if not powerTes:
        print('Conducting Voronoi tesselation...', end = '')
        vor, regions, vertices, polygons, areas, centroids, points = utilitiesNumeric.runMirroredVoronoi (node_coords, 2, maxLim)
    else:
        print('Conducting Power tesselation...', end = '')
        vor, regions, vertices, polygons, areas, centroids, points = utilitiesNumeric.runMirroredPower(node_coords, radii, 2, maxLim)
    print('done.')


    ########################################################################
    functions = []
    ### creating functions
    #### Defining functions
    #0 constant zero
    fn = utilitiesNumeric.constantFunc(0)
    functions.append (fn)

    #1 loading function, single force top right, bilinear
    func1 = []
    func1.append( np.array([0,0]) )
    func1.append( np.array([50, -50e4]) )
    fn1 = utilitiesNumeric.generalFunc(func1)
    functions.append (fn1)

    #transport function, leftFace, constant
    fn2 = utilitiesNumeric.constantFunc(20)
    functions.append (fn2)

    #transport function, rightFace, bilinear
    func3 = []
    func3.append( np.array([0,0]) )
    func3.append( np.array([50, 500]) )
    fn3 = utilitiesNumeric.generalFunc(func3)
    functions.append (fn3)

    ########################################################################
    ### indirect setting of transportBCs by spatial selection of vertices
    transportBC_merged = []
    transportIC_merged = []
    ### selecting vertices on the left surface
    leftFaceBC = np.array([2,-1])
    #leftFaceIC = 25.6
    boundA = np.array(  [-1e-8 , maxLim[1]/10*8] )
    boundB = np.array(  [ 1e-8 , maxLim[1]]  )
    leftFace = utilitiesGeom.returnSelectedPts(boundA, boundB, vor.vertices)
    #print(leftFace)
    for i in range (len(leftFace)):
        trsBC = utilitiesMech.transportBC(leftFace[i], leftFaceBC)
        transportBC_merged.append(trsBC)
        #trsIC = utilitiesGeom.transportIC(leftFace[i], leftFaceIC)
        #transportIC_merged.append(trsIC)

    ### selecting vertices on the right surface
    rightFaceBC = np.array([3,-1])
    boundA = np.array(  [maxLim[0] - 1e-8 , - 1e-8] )
    boundB = np.array(  [maxLim[0] + 1e-8 , maxLim[1]/10*2 ] )
    rightFace = utilitiesGeom.returnSelectedPts(boundA, boundB, vor.vertices)
    #print(rightFace)
    for i in range (len(rightFace)):
        trsBC = utilitiesMech.transportBC(rightFace[i], rightFaceBC)
        transportBC_merged.append(trsBC)

    return node_coords, mechBC_merged, mechIC_merged, transportBC_merged, transportIC_merged, vor, areas, functions, radii


def create2dCantileverUniTens(maxLim, minDist, trials):
    print('Creating 2d cantilever, uniform tension.')
    ### sampling of nodes
    ### direct setting of mechanicalBCs
    node_coords, mechBC_merged, mechIC_merged  = assemble2DCantileverUniTens(maxLim, minDist, trials );

    print('Conducting Voronoi tesselation...', end = '')
    vor, regions, vertices, polygons, areas, centroids, points = utilitiesNumeric.runMirroredVoronoi (node_coords, 2, maxLim)
    print('done.')


    ########################################################################
    functions = []
    #### Defining functions
    #0 constant zero
    fn = utilitiesNumeric.constantFunc(0)
    functions.append (fn)

    #1 loading function, single force top right, bilinear
    fn1 = utilitiesNumeric.sawToothConstFunc(value = -0.5e-4, period = 10)
    functions.append (fn1)

    #transport function, leftFace, constant
    fn2 = utilitiesNumeric.constantFunc(0)
    functions.append (fn2)

    #transport function, rightFace, bilinear
    fn3 = utilitiesNumeric.constantFunc(100)
    functions.append (fn3)



    ########################################################################
    ### indirect setting of transportBCs by spatial selection of vertices
    transportBC_merged = []
    transportIC_merged = []
    ### selecting vertices on the left surface
    leftFaceBC = np.array([2,-1])
    #leftFaceIC = 25.6
    boundA = np.array(  [-1e-8 , maxLim[1]/10*8] )
    boundB = np.array(  [ 1e-8 , maxLim[1]]  )
    leftFace = utilitiesGeom.returnSelectedPts(boundA, boundB, vor.vertices)
    #print(leftFace)
    for i in range (len(leftFace)):
        trsBC = utilitiesMech.transportBC(leftFace[i], leftFaceBC)
        transportBC_merged.append(trsBC)
        #trsIC = utilitiesGeom.transportIC(leftFace[i], leftFaceIC)
        #transportIC_merged.append(trsIC)

    ### selecting vertices on the right surface
    rightFaceBC = np.array([3,-1])
    boundA = np.array(  [maxLim[0] - 1e-8 , - 1e-8] )
    boundB = np.array(  [maxLim[0] + 1e-8 , maxLim[1]/10*2 ] )
    rightFace = utilitiesGeom.returnSelectedPts(boundA, boundB, vor.vertices)
    #print(rightFace)
    for i in range (len(rightFace)):
        trsBC = utilitiesMech.transportBC(rightFace[i], rightFaceBC)
        transportBC_merged.append(trsBC)

    return node_coords, mechBC_merged, mechIC_merged, transportBC_merged, transportIC_merged, vor, areas, functions


def create2dbeamConfinedPress(maxLim, minDist, trials ):
    print('Creating 2d cantilever, uniform tension.')
    ### sampling of nodes
    ### direct setting of mechanicalBCs
    node_coords, mechBC_merged, mechIC_merged  = assemble2dbeamConfinedPress(maxLim, minDist, trials );

    print('Conducting Voronoi tesselation...', end = '')
    vor, regions, vertices, polygons, areas, centroids, points = utilitiesNumeric.runMirroredVoronoi (node_coords, 2, maxLim)
    print('done.')


    ########################################################################
    functions = []
    #### Defining functions
    #0 constant zero
    fn = utilitiesNumeric.constantFunc(0)
    functions.append (fn)

    #1 loading function, single force top right, bilinear
    #fn1 = utilitiesNumeric.constantFunc(1e-3)#utilitiesNumeric.sawToothConstFunc(value = -0.5e-4, period = 10)
    fn1 = utilitiesNumeric.sawToothConstFunc(value = -0.5e-4, period = 10)
    functions.append (fn1)
    #functions.append (fn1)

    #transport function, leftFace, constant
    fn2 = utilitiesNumeric.constantFunc(0)
    functions.append (fn2)

    #transport function, rightFace, bilinear
    fn3 = utilitiesNumeric.constantFunc(100)
    functions.append (fn3)



    ########################################################################
    ### indirect setting of transportBCs by spatial selection of vertices
    transportBC_merged = []
    transportIC_merged = []
    ### selecting vertices on the left surface
    leftFaceBC = np.array([2,-1])
    #leftFaceIC = 25.6
    boundA = np.array(  [-1e-8 , maxLim[1]/10*8] )
    boundB = np.array(  [ 1e-8 , maxLim[1]]  )
    leftFace = utilitiesGeom.returnSelectedPts(boundA, boundB, vor.vertices)
    #print(leftFace)
    for i in range (len(leftFace)):
        trsBC = utilitiesMech.transportBC(leftFace[i], leftFaceBC)
        transportBC_merged.append(trsBC)
        #trsIC = utilitiesGeom.transportIC(leftFace[i], leftFaceIC)
        #transportIC_merged.append(trsIC)

    ### selecting vertices on the right surface
    rightFaceBC = np.array([3,-1])
    boundA = np.array(  [maxLim[0] - 1e-8 , - 1e-8] )
    boundB = np.array(  [maxLim[0] + 1e-8 , maxLim[1]/10*2 ] )
    rightFace = utilitiesGeom.returnSelectedPts(boundA, boundB, vor.vertices)
    #print(rightFace)
    for i in range (len(rightFace)):
        trsBC = utilitiesMech.transportBC(rightFace[i], rightFaceBC)
        transportBC_merged.append(trsBC)

    return node_coords, mechBC_merged, mechIC_merged, transportBC_merged, transportIC_merged, vor, areas, functions



def createPatchTestTransport(maxLim, minDist, trials, dim, powerTes):
    print('Creating patch test')
    ### sampling of nodes
    ### direct setting of mechanicalBCs
    node_coords, radii, mechBC_merged, mechIC_merged  = assemblePatchTestTransport(maxLim, minDist, trials, dim);

    print('Conducting Voronoi tesselation...', end = '')
    if not powerTes:
        if (dim==2):
            vor, regions, vertices, polygons, areas, centroids, points = utilitiesNumeric.runMirroredVoronoi (node_coords, dim, maxLim)
        else:
            vor, areas = utilitiesNumeric.runMirroredVoronoi (node_coords, dim, maxLim)
    else:
        if (dim==2):
            vor, regions, vertices, polygons, areas, centroids, points = utilitiesNumeric.runMirroredPower(node_coords, radii, 2, maxLim)
        else:
            vor, areas = utilitiesNumeric.runMirroredPower(node_coords, radii, 3, maxLim)
    print('done.')

    # if SHOW_PLOT:
    #     fig = voronoi_plot_2d(vor, show_vertices=True, line_colors='orange',  line_width=2, line_alpha=0.6, point_size=2)
    #     plt.show()

    ########################################################################
    ### indirect setting of transportBCs by spatial selection of vertices
    transportBC_merged = []
    functions = []

    #"""
    ### selecting vertices on the left surface
    boundA = np.zeros(dim)-1e-8
    boundB = maxLim + 1e-8
    faces1 = utilitiesGeom.returnSelectedPts(boundA, boundB, vor.vertices)
    vert = vor.vertices[faces1,:]
    boundA = np.zeros(dim)+1e-8
    boundB = maxLim - 1e-8
    faces0 = utilitiesGeom.excludeSelectedPts(boundA, boundB, vert)
    faces = faces1[faces0]

    for i,k in enumerate(faces):
        fn1 = utilitiesNumeric.constantFunc(np.sin(vor.vertices[k,0])*np.exp(vor.vertices[k,1]))
        functions.append (fn1)
        trsBC = utilitiesMech.transportBC(k,[i,-1])
        transportBC_merged.append(trsBC)
    #"""
    """
    #transport function, leftFace, constant
    fn2 = utilitiesNumeric.constantFunc(0)
    functions.append (fn2)
    fn3 = utilitiesNumeric.constantFunc(1)
    functions.append (fn3)


    ########################################################################
    ### indirect setting of transportBCs by spatial selection of vertices
    transportBC_merged = []
    transportIC_merged = []
    ### selecting vertices on the left surface
    leftFaceBC = np.array([0,-1])
    boundA = np.array(  [-1e-8 , 0] )
    boundB = np.array(  [ 1e-8 , maxLim[1]]  )
    leftFace = utilitiesGeom.returnSelectedPts(boundA, boundB, vor.vertices)
    #print(leftFace)
    for i in range (len(leftFace)):
        trsBC = utilitiesMech.transportBC(leftFace[i], leftFaceBC)
        transportBC_merged.append(trsBC)

    ### selecting vertices on the right surface
    rightFaceBC = np.array([1,-1])
    boundA = np.array(  [maxLim[0] - 1e-8, 0] )
    boundB = np.array(  [maxLim[0] + 1e-8 , maxLim[1]]  )
    rightFace = utilitiesGeom.returnSelectedPts(boundA, boundB, vor.vertices)
    #print(rightFace)
    for i in range (len(rightFace)):
        trsBC = utilitiesMech.transportBC(rightFace[i], rightFaceBC)
        transportBC_merged.append(trsBC)
    """

    return node_coords, [], transportBC_merged, vor, areas, functions, radii


def createCoupledCompression(maxLim, minDist, trials, dim, powerTes):
    print('Creating coupled compression test')
    node_coords, radii, mechBC_merged, mechIC_merged  = assemblePatchTestTransport(maxLim, minDist, trials, dim);
    node_coords_list = node_coords.tolist()

    indent = 1e-8
    nodeA = np.array([indent, indent])
    nodeB = np.array([indent, maxLim[1]-indent])
    pointGenerators.generateNodesLine2dRand(nodeA, nodeB, minDist, dim, node_coords_list,  trials, True, True)
    nodeA = np.array([maxLim[0]-indent, indent])
    nodeB = np.array([maxLim[0]-indent, maxLim[1]-indent])
    pointGenerators.generateNodesLine2dRand(nodeA, nodeB, minDist, dim, node_coords_list,  trials, True, True)
    node_coords = np.array(node_coords_list)
    radii = np.hstack((radii, np.zeros(len(node_coords)-len(radii))))
    #pointGenerators.generateNodesOrtoSurface3dRand(nodeA, nodeB, minDist, dim, node_coords, trials)


    print('Conducting Voronoi tesselation...', end = '')
    if not powerTes:
        if (dim==2):
            vor, regions, vertices, polygons, areas, centroids, points = utilitiesNumeric.runMirroredVoronoi (node_coords, dim, maxLim)
        else:
            vor, areas = utilitiesNumeric.runMirroredVoronoi (node_coords, dim, maxLim)
    else:
        if (dim==2):
            vor, regions, vertices, polygons, areas, centroids, points = utilitiesNumeric.runMirroredPower(node_coords, radii, 2, maxLim)
        else:
            vor, areas = utilitiesNumeric.runMirroredPower(node_coords, radii, 3, maxLim)
    print('done.')

    ########################################################################
    ### indirect setting of transportBCs by spatial selection of vertices
    transportBC_merged = []
    mechanicalBC_merged = []
    functions = []

    func = []
    func.append( np.array([0,0]) )
    func.append( np.array([1,-1e-3]) )

    functions.append ( utilitiesNumeric.constantFunc( 0. ) ) #constant BC
    functions.append ( utilitiesNumeric.generalFunc(func) ) #top movement
    functions.append ( utilitiesNumeric.constantFunc( 0. ) ) #bottom pressure
    functions.append ( utilitiesNumeric.constantFunc( 1e6 ) ) #top pressure

    ### selecting vertices and nodes on the top and bottom surface
    boundA = np.zeros(dim)-1e-7
    boundB = np.copy(maxLim)+1e-7
    boundB[0] = 1e-7
    botVert = utilitiesGeom.returnSelectedPts(boundA, boundB, vor.vertices)
    botNode = utilitiesGeom.returnSelectedPts(boundA, boundB, node_coords)
    boundA[0] += maxLim[0]
    boundB[0] += maxLim[0]
    topVert = utilitiesGeom.returnSelectedPts(boundA, boundB, vor.vertices)
    topNode = utilitiesGeom.returnSelectedPts(boundA, boundB, node_coords)

    nodeBC1 = -np.ones(6*dim-6).astype(int)
    nodeBC1[0] = 0
    nodeBC1X = np.copy(nodeBC1)
    nodeBC1X[1] = 0
    for i,k in enumerate(botNode):
        if (i==0): mBC = utilitiesMech.mechanicalBC(dim, k, nodeBC1X)
        else: mBC = utilitiesMech.mechanicalBC(dim, k, nodeBC1)
        mechanicalBC_merged.append(mBC)
    nodeBC2 = np.copy(nodeBC1)
    nodeBC2[0] = 1
    for i,k in enumerate(topNode):
        mBC = utilitiesMech.mechanicalBC(dim, k, nodeBC2)
        mechanicalBC_merged.append(mBC)

    for i,k in enumerate(botVert):
        trsBC = utilitiesMech.transportBC(k,[2,-1])
        transportBC_merged.append(trsBC)
    for i,k in enumerate(topVert):
        trsBC = utilitiesMech.transportBC(k,[3,-1])
        transportBC_merged.append(trsBC)

    return node_coords, mechanicalBC_merged, transportBC_merged, vor, areas, functions, radii

def create2DUniaxialTension(maxLim, minDist, trials, dim, powerTes):
    print('Creating uniaxial Tension')
    ### sampling of nodes
    ### direct setting of mechanicalBCs
    node_coords, radii, mechBC_merged, mechIC_merged  = assemblePatchTestTransport(maxLim, minDist, trials, dim);

    print('Conducting Voronoi tesselation...', end = '')
    if not powerTes:
        if (dim==2):
            vor, regions, vertices, polygons, areas, centroids, points = utilitiesNumeric.runMirroredVoronoi (node_coords, dim, maxLim)
        else:
            vor, areas = utilitiesNumeric.runMirroredVoronoi (node_coords, dim, maxLim)
    else:
        if (dim==2):
            vor, regions, vertices, polygons, areas, centroids, points = utilitiesNumeric.runMirroredPower(node_coords, radii, 2, maxLim)
        else:
            vor, areas = utilitiesNumeric.runMirroredPower(node_coords, radii, 3, maxLim)
    print('done.')

    # if SHOW_PLOT:
    #     fig = voronoi_plot_2d(vor, show_vertices=True, line_colors='orange',  line_width=2, line_alpha=0.6, point_size=2)
    #     plt.show()

    ########################################################################
    ### indirect setting of transportBCs by spatial selection of vertices
    transportBC_merged = []
    functions = []

    #"""
    ### selecting vertices on the left surface
    boundA = np.zeros(dim)-1e-8
    boundB = maxLim + 1e-8
    faces1 = utilitiesGeom.returnSelectedPts(boundA, boundB, vor.vertices)
    vert = vor.vertices[faces1,:]
    boundA = np.zeros(dim)+1e-8
    boundB = maxLim - 1e-8
    faces0 = utilitiesGeom.excludeSelectedPts(boundA, boundB, vert)
    faces = faces1[faces0]

    for i,k in enumerate(faces):
        fn1 = utilitiesNumeric.constantFunc(np.sin(vor.vertices[k,0])*np.exp(vor.vertices[k,1]))
        functions.append (fn1)
        trsBC = utilitiesMech.transportBC(k,[i,-1])
        transportBC_merged.append(trsBC)
    #"""
    """
    #transport function, leftFace, constant
    fn2 = utilitiesNumeric.constantFunc(0)
    functions.append (fn2)
    fn3 = utilitiesNumeric.constantFunc(1)
    functions.append (fn3)


    ########################################################################
    ### indirect setting of transportBCs by spatial selection of vertices
    transportBC_merged = []
    transportIC_merged = []
    ### selecting vertices on the left surface
    leftFaceBC = np.array([0,-1])
    boundA = np.array(  [-1e-8 , 0] )
    boundB = np.array(  [ 1e-8 , maxLim[1]]  )
    leftFace = utilitiesGeom.returnSelectedPts(boundA, boundB, vor.vertices)
    #print(leftFace)
    for i in range (len(leftFace)):
        trsBC = utilitiesMech.transportBC(leftFace[i], leftFaceBC)
        transportBC_merged.append(trsBC)

    ### selecting vertices on the right surface
    rightFaceBC = np.array([1,-1])
    boundA = np.array(  [maxLim[0] - 1e-8, 0] )
    boundB = np.array(  [maxLim[0] + 1e-8 , maxLim[1]]  )
    rightFace = utilitiesGeom.returnSelectedPts(boundA, boundB, vor.vertices)
    #print(rightFace)
    for i in range (len(rightFace)):
        trsBC = utilitiesMech.transportBC(rightFace[i], rightFaceBC)
        transportBC_merged.append(trsBC)
    """

    return node_coords, [], transportBC_merged, vor, areas, functions, radii


def create2dCoupledPress(maxLim, minDist, trials):
    dim = 2
    print('Creating coupled artificial crack')
    ### sampling of nodes
    ### direct setting of mechanicalBCs
    node_coords, mechBC_merged, mechInitC_merged, govNodes, govNodesMechBC, rigidPlates = assemble2dCoupledPress(maxLim, minDist, trials);

    print('Conducting Voronoi tesselation...', end = '')
    vor, regions, vertices, polygons, areas, centroids, points = utilitiesNumeric.runMirroredVoronoi (node_coords, dim, maxLim)
    print('done.')

    # if SHOW_PLOT:
    #     fig = voronoi_plot_2d(vor, show_vertices=True, line_colors='orange',  line_width=2, line_alpha=0.6, point_size=2)
    #     plt.show()

    functions = []

    #### Defining functions
    #0 constant zero
    fn = utilitiesNumeric.constantFunc(0)
    functions.append (fn)

    #1 loading function
    func2 = []
    func2.append( np.array([0,0]) )
    func2.append( np.array([1, 10e-3]) )
    fn2 = utilitiesNumeric.generalFunc(func2)
    functions.append (fn2)

    func3 = utilitiesNumeric.constantFunc(1)
    functions.append (func3)

    ########################################################################
    ### indirect setting of transportBCs by spatial selection of vertices
    transportBC_merged = []

    #"""
    ### indirect setting of transportBCs by spatial selection of vertices
    transportBC_merged = []
    transportIC_merged = []
    ### selecting vertices on the left surface
    boundA = np.array(  [-1e-4 , 1e-4] )
    boundB = np.array(  [ 1e-4 , maxLim[1]-1e-4]  )
    leftFace = utilitiesGeom.returnSelectedPts(boundA, boundB, vor.vertices)




    ### selecting vertices on the right surface
    boundA = np.array(  [ maxLim[0]-1e-4 , 1e-4] )
    boundB = np.array(  [ maxLim[0]+1e-4 , maxLim[1]-1e-4]  )
    rightFace = utilitiesGeom.returnSelectedPts(boundA, boundB, vor.vertices)

    ### selecting vertices on the top surface
    boundA = np.array(  [ -1e-5 , -1e-5] )
    boundB = np.array(  [ maxLim[0]+1e-5 , 1e-5]  )
    topFace = utilitiesGeom.returnSelectedPts(boundA, boundB, vor.vertices)

    ### selecting vertices on the bot surface
    boundA = np.array(  [ -1e-5 , maxLim[1]-1e-5] )
    boundB = np.array(  [ maxLim[0]+1e-5 , maxLim[1]+1e-5]  )
    botFace = utilitiesGeom.returnSelectedPts(boundA, boundB, vor.vertices)

    print ('vrtx %s' %leftFace)
    print ('vrtx %s' %rightFace)
    print ('vrtx %s' %topFace)
    print ('vrtx %s' %botFace)
    #setting of transport bc by using rigid plate
    rigidPlatesTrspt = []
    govNodesTrspt = []
    govNodesTrsptBC = []

    trsptLeftRigidPlate = utilitiesMech.RigidPlate(-len(govNodesTrspt)-1, 3, None, directIdcs = True)
    trsptLeftRigidPlate.setDirectNodes(leftFace)
    trsptLeftRigidPlateMechBC = np.array([0,-1])
    rigidPlatesTrspt.append(trsptLeftRigidPlate)
    govNodesTrspt.append(np.array([ 0, 0]))
    govNodesTrsptBC.append(utilitiesMech.transportBC(govNodesTrspt[-1], trsptLeftRigidPlateMechBC))

    trsptRightRigidPlate = utilitiesMech.RigidPlate(-len(govNodesTrspt)-1, 3, None, directIdcs = True)
    trsptRightRigidPlate.setDirectNodes(rightFace)
    trsptRightRigidPlateMechBC = np.array([2,-1])
    rigidPlatesTrspt.append(trsptRightRigidPlate)
    govNodesTrspt.append(np.array([ -1, -1]))
    govNodesTrsptBC.append(utilitiesMech.transportBC(govNodesTrspt[-1], trsptRightRigidPlateMechBC))


    trsptTopRigidPlate = utilitiesMech.RigidPlate(-len(govNodesTrspt)-1, 3, None, directIdcs = True)
    trsptTopRigidPlate.setDirectNodes(topFace)
    trsptTopRigidPlateMechBC = np.array([-1,0])
    rigidPlatesTrspt.append(trsptTopRigidPlate)
    govNodesTrspt.append(np.array([ -1, 2]))
    govNodesTrsptBC.append(utilitiesMech.transportBC(govNodesTrspt[-1], trsptTopRigidPlateMechBC))


    trsptBotRigidPlate = utilitiesMech.RigidPlate(-len(govNodesTrspt)-1, 3, None, directIdcs = True)
    trsptBotRigidPlate.setDirectNodes(botFace)
    trsptBotRigidPlateMechBC = np.array([-1,0])
    rigidPlatesTrspt.append(trsptBotRigidPlate)
    govNodesTrspt.append(np.array([ -1, 3]))
    govNodesTrsptBC.append(utilitiesMech.transportBC(govNodesTrspt[-1], trsptBotRigidPlateMechBC))


    return node_coords, mechBC_merged, transportBC_merged,   govNodes, govNodesMechBC, rigidPlates, vor, areas, functions,  rigidPlatesTrspt, govNodesTrspt, govNodesTrsptBC




def create2dCoupledArtificialCrack(maxLim, minDist, trials, notchH):
    dim = 2
    print('Creating coupled artificial crack')
    ### sampling of nodes
    ### direct setting of mechanicalBCs
    slitWidth = minDist*0.8
    node_coords, mechBC_merged, mechInitC_merged, notches, govNodes, govNodesMechBC, rigidPlates = assemble2dCoupledArtificialCrack(maxLim, minDist, trials, slitWidth, notchH);

    print('Conducting Voronoi tesselation...', end = '')
    vor, regions, vertices, polygons, areas, centroids, points = utilitiesNumeric.runMirroredVoronoi (node_coords, dim, maxLim)
    print('done.')

    # if SHOW_PLOT:
    #     fig = voronoi_plot_2d(vor, show_vertices=True, line_colors='orange',  line_width=2, line_alpha=0.6, point_size=2)
    #     plt.show()

    functions = []

    #### Defining functions
    #0 constant zero
    fn = utilitiesNumeric.constantFunc(0)
    functions.append (fn)

    #1 loading function
    func2 = []
    func2.append( np.array([0,0]) )
    func2.append( np.array([1, -1e-3]) )
    fn2 = utilitiesNumeric.generalFunc(func2)
    functions.append (fn2)

    func3 = []
    func3.append( np.array([0,0]) )
    func3.append( np.array([1, 1e-3]) )
    fn3 = utilitiesNumeric.generalFunc(func3)
    functions.append (fn3)

    ########################################################################
    ### indirect setting of transportBCs by spatial selection of vertices
    transportBC_merged = []


    #"""
    ### selecting vertices on the bottom surface
    boundA = np.zeros(2)-1e-8
    boundB = np.array([maxLim[0], 1e-8])
    faces1 = utilitiesGeom.returnSelectedPts(boundA, boundB, vor.vertices)
    vert = vor.vertices[faces1,:]

    fn1 = utilitiesNumeric.constantFunc(100)
    functions.append (fn1)

    for i,k in enumerate(faces1):
        trsBC = utilitiesMech.transportBC(k,[3,-1])
        transportBC_merged.append(trsBC)

    #"""
    ### selecting vertices on the top surface
    boundA = np.array([-1e-8, maxLim[1]-1e-8])
    boundB = np.array([maxLim[0], maxLim[1]+1e-8])
    faces1 = utilitiesGeom.returnSelectedPts(boundA, boundB, vor.vertices)
    vert = vor.vertices[faces1,:]

    fn1 = utilitiesNumeric.constantFunc(0)
    functions.append (fn1)

    for i,k in enumerate(faces1):
        trsBC = utilitiesMech.transportBC(k,[4,-1])
        transportBC_merged.append(trsBC)

    #"""


    return node_coords, mechBC_merged, transportBC_merged,  notches, govNodes, govNodesMechBC, rigidPlates, vor, areas, functions






def create3dCoupledArtificialCrack(maxLim, minDist, trials, notchH):
    dim = 3
    print('Creating coupled artificial crack')
    ### sampling of nodes
    ### direct setting of mechanicalBCs
    slitWidth = minDist*1
    node_coords, mechBC_merged, mechInitC_merged, notches, govNodes, govNodesMechBC, rigidPlates = assemble3dCoupledArtificialCrack(maxLim, minDist, trials, slitWidth, notchH);

    print('Conducting Voronoi tesselation...', end = '')
    vor, volumes = utilitiesNumeric.runMirroredVoronoi (node_coords, dim, maxLim)
    print('done.')

    # if SHOW_PLOT:
    #     fig = voronoi_plot_2d(vor, show_vertices=True, line_colors='orange',  line_width=2, line_alpha=0.6, point_size=2)
    #     plt.show()

    functions = []

    #### Defining functions
    #0 constant zero
    fn = utilitiesNumeric.constantFunc(0)
    functions.append (fn)

    #1 loading function
    func2 = []
    func2.append( np.array([0,0]) )
    func2.append( np.array([1, -1e-3]) )
    fn2 = utilitiesNumeric.generalFunc(func2)
    functions.append (fn2)

    func3 = []
    func3.append( np.array([0,0]) )
    func3.append( np.array([1, 1e-3]) )
    fn3 = utilitiesNumeric.generalFunc(func3)
    functions.append (fn3)

    ########################################################################
    ### indirect setting of transportBCs by spatial selection of vertices
    transportBC_merged = []


    #"""
    ### selecting vertices on the bottom surface
    boundA = np.zeros(3)-1e-8
    boundB = np.array([maxLim[0], 1e-8, maxLim[2]+1e-8])
    faces1 = utilitiesGeom.returnSelectedPts(boundA, boundB, vor.vertices)
    vert = vor.vertices[faces1,:]

    fn1 = utilitiesNumeric.constantFunc(100)
    functions.append (fn1)
    """
    for i,k in enumerate(faces1):
        trsBC = utilitiesMech.transportBC(k,[3,-1])
        transportBC_merged.append(trsBC)

    #"""
    ### selecting vertices on the top surface
    boundA = np.array([-1e-8, maxLim[1]-1e-8, -1e-8])
    boundB = np.array([maxLim[0], maxLim[1]+1e-8, maxLim[2]+1e-8])
    faces2 = utilitiesGeom.returnSelectedPts(boundA, boundB, vor.vertices)
    vert = vor.vertices[faces2,:]

    fn1 = utilitiesNumeric.constantFunc(0)
    functions.append (fn1)
    """
    for i,k in enumerate(faces2):
        trsBC = utilitiesMech.transportBC(k,[4,-1])
        transportBC_merged.append(trsBC)

    #"""

    govNodesTrspt = []
    govNodesTrsptBC = []
    rigidPlatesTrspt = []

    trsptLeftRigidPlate = utilitiesMech.RigidPlate(-len(govNodesTrspt)-1, 3, None, directIdcs = True)
    trsptLeftRigidPlate.setDirectNodes(faces1)
    trsptLeftRigidPlateMechBC = np.array([3,-1])
    rigidPlatesTrspt.append(trsptLeftRigidPlate)
    govNodesTrspt.append(np.array([ 0, 0, 0]))
    govNodesTrsptBC.append(utilitiesMech.transportBC(govNodesTrspt[-1], trsptLeftRigidPlateMechBC))

    trsptRightRigidPlate = utilitiesMech.RigidPlate(-len(govNodesTrspt)-1, 3, None, directIdcs = True)
    trsptRightRigidPlate.setDirectNodes(faces2)
    trsptRightRigidPlateMechBC = np.array([4,-1])
    rigidPlatesTrspt.append(trsptRightRigidPlate)
    govNodesTrspt.append(np.array([ -1, -1, -1]))
    govNodesTrsptBC.append(utilitiesMech.transportBC(govNodesTrspt[-1], trsptRightRigidPlateMechBC))

    #return node_coords, mechBC_merged, transportBC_merged,  govNodes, govNodesMechBC, rigidPlates, vor, functions, rigidPlatesTrspt, govNodesTrspt, govNodesTrsptBC

    return node_coords, mechBC_merged, transportBC_merged, govNodes, govNodesMechBC, rigidPlates, vor, volumes, functions, rigidPlatesTrspt, govNodesTrspt, govNodesTrsptBC




def create2d_CFRAC_Clover(maxLim, minDist, trials, holeMinDist, holeDiameter,baseMinDist,fineWidth,fineHeight, interfaceMinDist=-1, roughMinDistCoef=1):
    print('Creating CFRAC clover model...')
    dim=2

    ### sampling of nodes
    ### direct setting of mechanicalBCs
    sampleBorders = True
    node_coords, mechBC_merged, mechInitC_merged,  govNodes, govNodesMechBC, rigidPlates, functions  = assemble2d_CFRAC_Clover(maxLim, minDist, trials, holeMinDist, holeDiameter,baseMinDist,fineWidth,fineHeight)


    print('Conducting Voronoi tesselation...', end = '')
    vor, regions, vertices, polygons, areas, centroids, points = utilitiesNumeric.runMirroredVoronoiClover (node_coords, dim, maxLim, holeDiameter)
    print('done.')


    #### Defining functions
    #0 constant zero
    fn = utilitiesNumeric.constantFunc(0)
    functions.append (fn)

    #2 mech loading function
    func2 = []
    func2.append( np.array([0,0]) )
    func2.append( np.array([1, 1e-4]) )
    fn2 = utilitiesNumeric.generalFunc(func2)
    functions.append (fn2)


    return node_coords, mechBC_merged, [], govNodes, govNodesMechBC, rigidPlates, vor, areas, functions



def create3d_CFRAC_Clover(maxLim, minDist, trials, holeMinDist, holeDiameter, interfaceMinDist=-1, roughMinDistCoef=1):
    print('Creating CFRAC 3d clover model...')
    dim=3

    ### sampling of nodes
    ### direct setting of mechanicalBCs
    sampleBorders = True
    node_coords, mechBC_merged, mechInitC_merged,  govNodes, govNodesMechBC, rigidPlates, functions  = assemble3d_CFRAC_Clover(maxLim, minDist, trials, holeMinDist, holeDiameter)

    print('Conducting Voronoi tesselation...', end = '')
    vor, volumes= utilitiesNumeric.runMirroredVoronoiClover (node_coords, dim, maxLim, holeDiameter)
    print('done.')


    #### Defining functions
    #0 constant zero
    fn = utilitiesNumeric.constantFunc(0)
    functions.append (fn)

    #2 mech loading function
    func2 = []
    func2.append( np.array([0,0]) )
    func2.append( np.array([1, 1e-4]) )
    fn2 = utilitiesNumeric.generalFunc(func2)
    functions.append (fn2)


    return node_coords, mechBC_merged, [], govNodes, govNodesMechBC, rigidPlates, vor, [], functions





def create2d_CFRAC_TDCB(maxLim, minDist, trials, holeMinDist, holeDiameter, interfaceMinDist=-1, roughMinDistCoef=1, elazonewidth=10*0.001):
    print('Creating CFRAC TDCB model...')
    dim=2

    ### sampling of nodes
    ### direct setting of mechanicalBCs
    sampleBorders = True
    node_coords, mechBC_merged, mechInitC_merged,  govNodes, govNodesMechBC, rigidPlates, functions,notches,node_indices  = assemble2d_CFRAC_TDCB(maxLim, minDist, trials, holeMinDist, holeDiameter,roughMinDistCoef,elazonewidth)


    print('Conducting Voronoi tesselation...', end = '')
    vor, regions, vertices, polygons, areas, centroids, points = utilitiesNumeric.runMirroredVoronoiTDCB (node_coords, dim, maxLim, holeDiameter)
    print('done.')


    #### Defining functions
    #0 constant zero
    fn = utilitiesNumeric.constantFunc(0)
    functions.append (fn)

    #2 mech loading function
    func2 = []
    func2.append( np.array([0,0]) )
    func2.append( np.array([1, 5e-5]) )
    fn2 = utilitiesNumeric.generalFunc(func2)
    functions.append (fn2)


    return node_coords, mechBC_merged, [], govNodes, govNodesMechBC, rigidPlates, vor, areas, functions,notches,node_indices



def create2d_Hanging_FracZone(maxLim, minDist, trials):
    print('Creating 2d frac zone...')
    dim=2

    ### sampling of nodes
    ### direct setting of mechanicalBCs
    node_coords, mechBC_merged, mechInitC_merged,  govNodes, govNodesMechBC, rigidPlates, functions  = assemble2d_Hanging_FracZone(maxLim, minDist, trials)

    print('Conducting Voronoi tesselation...', end = '')
    vor, regions, vertices, polygons, areas, centroids, points = utilitiesNumeric.runMirroredVoronoi (node_coords, dim, maxLim)
    print('done.')


    #### Defining functions
    #0 constant zero
    fn = utilitiesNumeric.constantFunc(0)
    functions.append (fn)

    #2 mech loading function
    func2 = []
    func2.append( np.array([0,0]) )
    func2.append( np.array([1, 3e-5]) )
    fn2 = utilitiesNumeric.generalFunc(func2)
    functions.append (fn2)


    return node_coords, mechBC_merged, [], govNodes, govNodesMechBC, rigidPlates, vor, [], functions

def create2d_box_with_periodic_nodes(maxLim, minDist, trials):
    print('Creating 2d frac zone...')
    dim=2

    ### sampling of nodes
    ### direct setting of mechanicalBCs
    node_coords, mechBC_merged, mechInitC_merged,  govNodes, govNodesMechBC, rigidPlates, functions, radii  = assemble2d_box_with_periodic_nodes(maxLim, minDist, trials)

    print('Conducting Power tesselationSSSSSSS...', end = '')
    vor, regions, vertices, polygons, areas, centroids, points = utilitiesNumeric.runMirroredVoronoi (node_coords, dim, maxLim)
    #vor, regions, vertices, polygons, areas, centroids, points = utilitiesNumeric.runMirroredPower(node_coords, radii, dim, maxLim)
    print('done.')


    #### Defining functions
    #0 constant zero
    fn = utilitiesNumeric.constantFunc(0)
    functions.append (fn)

    #2 mech loading function
    func2 = []
    func2.append( np.array([0,0]) )
    func2.append( np.array([1, 3e-5]) )
    fn2 = utilitiesNumeric.generalFunc(func2)
    functions.append (fn2)


    return node_coords, mechBC_merged, [], govNodes, govNodesMechBC, rigidPlates, vor, [], functions, radii



def create3d_Hanging_FracZone(maxLim, minDist, topsupwidth, trials):
    print('Creating 3d frac zone...')
    dim=3

    ### sampling of nodes
    ### direct setting of mechanicalBCs
    node_coords, mechBC_merged, mechInitC_merged,  govNodes, govNodesMechBC, rigidPlates, functions, radii  = assemble3d_Hanging_FracZone(maxLim, minDist, topsupwidth, trials)

    #print('Conducting Voronoi tesselation...', end = '')
    #vor, volumes= utilitiesNumeric.runMirroredVoronoi (node_coords, dim, maxLim)
    #print('done.')

    print('Conducting Power tesselationSSSSSSS...', end = '')
    #vor, regions, vertices, polygons, areas, centroids, points = utilitiesNumeric.runMirroredVoronoi (node_coords, dim, maxLim)
    vor, areas = utilitiesNumeric.runMirroredPower(node_coords, radii, dim, maxLim)
    print('done.')


    #### Defining functions
    #0 constant zero
    fn = utilitiesNumeric.constantFunc(0)
    functions.append (fn)

    #2 mech loading function
    func2 = []
    func2.append( np.array([0,0]) )
    func2.append( np.array([1, 3e-5]) )
    fn2 = utilitiesNumeric.generalFunc(func2)
    functions.append (fn2)

    return node_coords, mechBC_merged, [], govNodes, govNodesMechBC, rigidPlates, vor, [], functions, radii





def create3d_CFRAC_TDCB(maxLim, minDist, trials, holeMinDist, holeDiameter, interfaceMinDist=-1, roughMinDistCoef=1, elazonewidth=10*0.001, notchWidth=0,fracZoneHeight=0.8,fracZoneOverhang=0.2, fem=False):
    print('Creating CFRAC 3d TDCB model...')
    dim=3

    if notchWidth == 0:
        notchWidth = minDist

    ### sampling of nodes
    ### direct setting of mechanicalBCs
    sampleBorders = True
    node_coords, mechBC_merged, mechInitC_merged,  govNodes, govNodesMechBC, rigidPlates, functions,notches,node_indices  = assemble3d_CFRAC_TDCB(maxLim, minDist, trials, holeMinDist, holeDiameter, roughMinDistCoef,elazonewidth, notchWidth, fracZoneOverhang, fracZoneHeight, fem=fem)

    print('Conducting Voronoi tesselation...', end = '')
    vor, volumes= utilitiesNumeric.runMirroredVoronoiTDCB (node_coords, dim, maxLim, holeDiameter)
    print('done.')


    #### Defining functions
    #0 constant zero
    fn = utilitiesNumeric.constantFunc(0)
    functions.append (fn)

    #2 mech loading function
    func2 = []
    func2.append( np.array([0,0]) )
    func2.append( np.array([1, 3e-5]) )
    fn2 = utilitiesNumeric.generalFunc(func2)
    functions.append (fn2)


    return node_coords, mechBC_merged, [], govNodes, govNodesMechBC, rigidPlates, vor, [], functions,notches,node_indices




def create2dCorrosionRebar(maxLim, minDist, trials, rebarMinDist, rebarDiameter, rebarCount, rebarDepth, node_coords_init=None, interfaceMinDist=-1, roughMinDistCoef=1, adaptivityReady=False, fineRingThickness=1, fineRegDepth=0.5, gradientRegDepth=0.2):
    print('Creating corrosion rebar model...')
    dim=2
    print(roughMinDistCoef)

    ### sampling of nodes
    ### direct setting of mechanicalBCs
    sampleBorders = True
    node_coords, mechBC_merged, mechInitC_merged,  govNodes, govNodesMechBC, rigidPlates, functions, interfaceNodeIndices  = assemble2dCorrosionRebar(maxLim, minDist, trials, rebarMinDist, interfaceMinDist, rebarDiameter, rebarCount, rebarDepth, sampleBorders, node_coords_init=node_coords_init, roughMinDistCoef=roughMinDistCoef, adaptivityReady=adaptivityReady, fineRingThickness=fineRingThickness, fineRegDepth=fineRegDepth, gradientRegDepth=gradientRegDepth)


    print('Conducting Voronoi tesselation...', end = '')
    vor, regions, vertices, polygons, areas, centroids, points = utilitiesNumeric.runMirroredVoronoiRebars (node_coords, dim, maxLim, rebarDiameter, rebarDepth, rebarCount)
    print('done.')


    #### Defining functions
    #0 constant zero
    fn = utilitiesNumeric.constantFunc(0)
    functions.append (fn)

    #1 mech loading function
    fn = utilitiesNumeric.constantFunc(0)
    functions.append (fn)

    #2 mech loading function
    func2 = []
    func2.append( np.array([0,0]) )
    func2.append( np.array([1, 1e-2]) )
    fn2 = utilitiesNumeric.generalFunc(func2)
    functions.append (fn2)

    ########################################################################
    ### indirect setting of transportBCs by spatial selection of vertices
    transportBC_merged = []
    govNodesTrspt = []
    govNodesTrsptBC = []
    rigidPlatesTrspt = []



    ### selecting vertices on the bottom surface
    boundA = np.zeros(2)-1e-8
    boundB = np.array([maxLim[0], 1e-8])
    faces1 = utilitiesGeom.returnSelectedPts(boundA, boundB, vor.vertices)
    vert = vor.vertices[faces1,:]

    #fn1 = utilitiesNumeric.constantFunc(0)
    #functions.append (fn1)

    trsptBottomRigidPlate = utilitiesMech.RigidPlate(-len(govNodesTrspt)-1, 3, None, directIdcs = True)
    trsptBottomRigidPlate.setDirectNodes(faces1)
    trsptBottomRigidPlateMechBC = np.array([0,-1])
    rigidPlatesTrspt.append(trsptBottomRigidPlate)
    govNodesTrspt.append(np.array([ 0, 0]))
    govNodesTrsptBC.append(utilitiesMech.transportBC(govNodesTrspt[-1], trsptBottomRigidPlateMechBC))



    ### selecting vertices on the top surface
    boundA = np.array([-1e-8, maxLim[1]-1e-8])
    boundB = np.array([maxLim[0], maxLim[1]+1e-8])
    faces2 = utilitiesGeom.returnSelectedPts(boundA, boundB, vor.vertices)
    vert = vor.vertices[faces1,:]

    trsptTopRigidPlate = utilitiesMech.RigidPlate(-len(govNodesTrspt)-1, 2, None, directIdcs = True)
    trsptTopRigidPlate.setDirectNodes(faces2)
    trsptTopRigidPlateMechBC = np.array([0,-1])
    rigidPlatesTrspt.append(trsptTopRigidPlate)
    govNodesTrspt.append(np.array([ -1, 1]))
    govNodesTrsptBC.append(utilitiesMech.transportBC(govNodesTrspt[-1], trsptTopRigidPlateMechBC))



    ### selecting vertices on the left surface
    boundA = np.array([-1e-8, 1e-8])
    boundB = np.array([1e-8, maxLim[1]-1e-8])
    faces3 = utilitiesGeom.returnSelectedPts(boundA, boundB, vor.vertices)
    vert = vor.vertices[faces1,:]

    trsptRightRigidPlate = utilitiesMech.RigidPlate(-len(govNodesTrspt)-1, 2, None, directIdcs = True)
    trsptRightRigidPlate.setDirectNodes(faces3)
    trsptRightRigidPlateMechBC = np.array([0,-1])
    rigidPlatesTrspt.append(trsptRightRigidPlate)
    govNodesTrspt.append(np.array([ -1, 2]))
    govNodesTrsptBC.append(utilitiesMech.transportBC(govNodesTrspt[-1], trsptRightRigidPlateMechBC))



    ### selecting vertices on the left surface
    boundA = np.array([maxLim[0]-1e-8, 1e-8])
    boundB = np.array([maxLim[0]+1e-8, maxLim[1]-1e-8])
    faces4 = utilitiesGeom.returnSelectedPts(boundA, boundB, vor.vertices)
    vert = vor.vertices[faces1,:]

    trsptLeftRigidPlate = utilitiesMech.RigidPlate(-len(govNodesTrspt)-1, 2, None, directIdcs = True)
    trsptLeftRigidPlate.setDirectNodes(faces4)
    trsptLeftRigidPlateMechBC = np.array([0,-1])
    rigidPlatesTrspt.append(trsptLeftRigidPlate)
    govNodesTrspt.append(np.array([ -1, -1]))
    govNodesTrsptBC.append(utilitiesMech.transportBC(govNodesTrspt[-1], trsptLeftRigidPlateMechBC))




    funcTrsprt = []
    funcTrsprt.append( np.array([0,0]) )
    funcTrsprt.append( np.array([1, 100]) )
    fnTrsprt = utilitiesNumeric.generalFunc(funcTrsprt)
    functions.append (fnTrsprt)

    funcTrsprt1 = []
    funcTrsprt1.append( np.array([0,1]) )
    #funcTrsprt1.append( np.array([1, 1]) )
    fnTrsprt1 = utilitiesNumeric.generalFunc(funcTrsprt1)
    functions.append (fnTrsprt1)

    ### selecting vertices on rebar interfaces
    clipRadius = rebarDiameter/2*1.01
    for r in range (rebarCount):
        interfaceVertices = []
        for i in range (len(vor.vertices)):
            if (rebarCount==1):
                #puvodni poloha rebars polovina od kraje
                center = np.array([ (maxLim[0]/rebarCount)*(r+0.5), maxLim[1]-rebarDepth  ])
            else:
                #poloha rebars presne jak je ve clanku
                center = np.array([ (0.058 + (maxLim[0]-0.116)/(rebarCount-1)*r), maxLim[1]-rebarDepth  ])

            vertex = vor.vertices[i]
            if (np.linalg.norm(vertex-center) < clipRadius):
                interfaceVertices.append(i)

        #print(interfaceVertices)
        trsptRebarRigidPlate = utilitiesMech.RigidPlate(-len(govNodesTrspt)-1, 4, None, directIdcs = True)
        trsptRebarRigidPlate.setDirectNodes(interfaceVertices.copy())
        trsptRebarRigidPlateMechBC = np.array([2,-1])
        rigidPlatesTrspt.append(trsptRebarRigidPlate)
        govNodesTrspt.append(center-1e-4)
        govNodesTrsptBC.append(utilitiesMech.transportBC(govNodesTrspt[-1], trsptRebarRigidPlateMechBC))



    return node_coords, mechBC_merged, transportBC_merged, govNodes, govNodesMechBC, rigidPlates, vor, areas, functions, rigidPlatesTrspt, govNodesTrspt, govNodesTrsptBC, interfaceNodeIndices





def create3dCorrosionRebar(maxLim, minDist, trials, rebarMinDist, rebarDiameter, rebarCount, rebarDepth, node_coords_init=None, interfaceMinDist=-1, roughMinDistCoef=1, adaptivityReady=False, fineRingThickness=1, fineRegDepth=0.5, gradientRegDepth=0.2):
    print('Creating 3d corrosion rebar model...')
    dim=3

    ### sampling of nodes
    ### direct setting of mechanicalBCs
    sampleBorders = True
    node_coords, mechBC_merged, mechInitC_merged,  govNodes, govNodesMechBC, rigidPlates, functions, interfaceNodeIndices  = assemble3dCorrosionRebar(maxLim, minDist, trials, rebarMinDist, interfaceMinDist, rebarDiameter, rebarCount, rebarDepth, sampleBorders, node_coords_init=node_coords_init, roughMinDistCoef=roughMinDistCoef, adaptivityReady=adaptivityReady, fineRingThickness=fineRingThickness, fineRegDepth=fineRegDepth, gradientRegDepth=gradientRegDepth)

    print('Conducting Voronoi tesselation...', end = '')
    #vor, regions, vertices, polygons, areas, centroids, points = utilitiesNumeric.runMirroredVoronoi (node_coords, dim, maxLim)
    vor, areas = utilitiesNumeric.runMirroredVoronoiRebars (node_coords, dim, maxLim, rebarDiameter, rebarDepth, rebarCount)
    print('done.')

    # if SHOW_PLOT:
    #     fig = voronoi_plot_2d(vor, show_vertices=True, line_colors='orange',  line_width=2, line_alpha=0.6, point_size=2)
    #     plt.show()

    #functions = []

    #### Defining functions
    #0 constant zero
    fn = utilitiesNumeric.constantFunc(0)
    functions.append (fn)

    #1 mech loading function
    fn = utilitiesNumeric.constantFunc(0)
    functions.append (fn)

    #2 mech loading function
    func2 = []
    func2.append( np.array([0,0]) )
    func2.append( np.array([1, 1e-2]) )
    fn2 = utilitiesNumeric.generalFunc(func2)
    functions.append (fn2)

    ########################################################################
    ### indirect setting of transportBCs by spatial selection of vertices
    transportBC_merged = []
    govNodesTrspt = []
    govNodesTrsptBC = []
    rigidPlatesTrspt = []



    ### selecting vertices on the bottom surface

    boundA = np.zeros(3)-1e-8
    boundB = np.array([maxLim[0]+1e-8, 1e-8, maxLim[2]+1e-8])
    faces1 = utilitiesGeom.returnSelectedPts(boundA, boundB, vor.vertices)
    vert = vor.vertices[faces1,:]

    #fn1 = utilitiesNumeric.constantFunc(0)
    #functions.append (fn1)

    trsptBottomRigidPlate = utilitiesMech.RigidPlate(-len(govNodesTrspt)-1, 3, None, directIdcs = True)
    trsptBottomRigidPlate.setDirectNodes(faces1)
    trsptBottomRigidPlateMechBC = np.array([0,-1])
    rigidPlatesTrspt.append(trsptBottomRigidPlate)
    govNodesTrspt.append(np.array([ 0, 0,0]))
    govNodesTrsptBC.append(utilitiesMech.transportBC(govNodesTrspt[-1], trsptBottomRigidPlateMechBC))



    ### selecting vertices on the top surface
    boundA = np.array([-1e-8, maxLim[1]-1e-8, -1e-8])
    boundB = np.array([maxLim[0]+1e-8, maxLim[1]+1e-8, maxLim[2]+1e-8])
    faces2 = utilitiesGeom.returnSelectedPts(boundA, boundB, vor.vertices)
    vert = vor.vertices[faces1,:]

    trsptTopRigidPlate = utilitiesMech.RigidPlate(-len(govNodesTrspt)-1, 3, None, directIdcs = True)
    trsptTopRigidPlate.setDirectNodes(faces2)
    trsptTopRigidPlateMechBC = np.array([-1,-1])
    rigidPlatesTrspt.append(trsptTopRigidPlate)
    govNodesTrspt.append(np.array([ 0, 0,1]))
    govNodesTrsptBC.append(utilitiesMech.transportBC(govNodesTrspt[-1], trsptTopRigidPlateMechBC))



    ### selecting vertices on the left surface
    boundA = np.array([-1e-8, 1e-7, -1e-8])
    boundB = np.array([1e-8, maxLim[1]-1e-7, maxLim[2]+1e-8])
    faces3 = utilitiesGeom.returnSelectedPts(boundA, boundB, vor.vertices)
    vert = vor.vertices[faces1,:]

    trsptRightRigidPlate = utilitiesMech.RigidPlate(-len(govNodesTrspt)-1, 3, None, directIdcs = True)
    trsptRightRigidPlate.setDirectNodes(faces3)
    trsptRightRigidPlateMechBC = np.array([-1,-1])
    rigidPlatesTrspt.append(trsptRightRigidPlate)
    govNodesTrspt.append(np.array([ 0, 0,2]))
    govNodesTrsptBC.append(utilitiesMech.transportBC(govNodesTrspt[-1], trsptRightRigidPlateMechBC))



    ### selecting vertices on the left surface
    boundA = np.array([maxLim[0]-1e-8, 1e-7, -1e-8])
    boundB = np.array([maxLim[0]+1e-8, maxLim[1]-1e-7,  maxLim[2]+1e-8])
    faces4 = utilitiesGeom.returnSelectedPts(boundA, boundB, vor.vertices)
    vert = vor.vertices[faces1,:]

    trsptLeftRigidPlate = utilitiesMech.RigidPlate(-len(govNodesTrspt)-1, 3, None, directIdcs = True)
    trsptLeftRigidPlate.setDirectNodes(faces4)
    trsptLeftRigidPlateMechBC = np.array([-1,-1])
    rigidPlatesTrspt.append(trsptLeftRigidPlate)
    govNodesTrspt.append(np.array([ 0, 3, 0]))
    govNodesTrsptBC.append(utilitiesMech.transportBC(govNodesTrspt[-1], trsptLeftRigidPlateMechBC))




    funcTrsprt = []
    funcTrsprt.append( np.array([0,0]) )
    funcTrsprt.append( np.array([1, 100]) )
    fnTrsprt = utilitiesNumeric.generalFunc(funcTrsprt)
    functions.append (fnTrsprt)

    funcTrsprt1 = []
    funcTrsprt1.append( np.array([0,1]) )
    #funcTrsprt1.append( np.array([1, 1]) )
    fnTrsprt1 = utilitiesNumeric.generalFunc(funcTrsprt1)
    functions.append (fnTrsprt1)

    ### selecting vertices on rebar interfaces
    clipRadius = rebarDiameter/2 + minDist/3

    interfaceVertexIndices = []
    for r in range (rebarCount):
        interfaceVertices = []
        for i in range (len(vor.vertices)):
            if (rebarCount==1):
                #puvodni poloha rebars polovina od kraje
                center = np.array([ (maxLim[0]/rebarCount)*(r+0.5), maxLim[1]-rebarDepth ,0 ])
            else:
                #poloha rebars presne jak je ve clanku
                center = np.array([ (0.058 + (maxLim[0]-0.116)/(rebarCount-1)*r), maxLim[1]-rebarDepth,0  ])

            vertex = vor.vertices[i]
            if (np.linalg.norm(vertex[0:2]-center[0:2]) < clipRadius):
                interfaceVertices.append(i)

        interfaceVertexIndices.append(interfaceVertices.copy())

        #print(interfaceVertices)
        trsptRebarRigidPlate = utilitiesMech.RigidPlate(-len(govNodesTrspt)-1, 3, None, directIdcs = True)
        trsptRebarRigidPlate.setDirectNodes(interfaceVertices.copy())
        trsptRebarRigidPlateMechBC = np.array([4,-1])
        rigidPlatesTrspt.append(trsptRebarRigidPlate)
        govNodesTrspt.append(center-1e-5)
        govNodesTrsptBC.append(utilitiesMech.transportBC(govNodesTrspt[-1], trsptRebarRigidPlateMechBC))

    interfaceVertexIndices = np.asarray(interfaceVertexIndices)
    return node_coords, mechBC_merged, transportBC_merged, govNodes, govNodesMechBC, rigidPlates, vor, areas, functions, rigidPlatesTrspt, govNodesTrspt, govNodesTrsptBC, interfaceVertexIndices





def createCoupledBrazilianDisc(center, cylinderRad, cylinderHeight,  minDist, trials, powerTes=False):
    dim = 3
    print('Creating coupled brazilian disc...')
    ### sampling of nodes
    ### direct setting of mechanicalBCs
    node_coords, mechBC_merged, mechInitC_merged, govNodes, govNodesMechBC, rigidPlates, radii  = assembleCoupledBrazilianDisc(center, cylinderRad, cylinderHeight, minDist, trials, 0, powerTes = powerTes )


    fig = plt.figure()
    ax = Axes3D(fig)
    #ax.auto_scale_xyz([0, maxLim[0]], [0, maxLim[1]], [0, maxLim[2]])
    ax.scatter(node_coords[:,0], node_coords[:,1], node_coords[:,2])
    plt.show()


    utilitiesGeom.saveNodes('', node_coords, "node",dim, "nodes.inp")




    print('Conducting Voronoi tesselation...', end='')
    vor, volumes = utilitiesNumeric.runCylinderMirroredVoronoi (node_coords, center, cylinderRad, cylinderHeight, 0, quarter=False, weights=radii)
    print('done.')


    functions = []

    #### Defining functions
    #0 constant zero
    fn = utilitiesNumeric.constantFunc(0)
    functions.append (fn)

    #1 loading function
    func2 = []
    func2.append( np.array([0,0]) )
    func2.append( np.array([1, -1e-3]) )
    fn2 = utilitiesNumeric.generalFunc(func2)
    functions.append (fn2)

    #2 unitary pressure
    fn3 = utilitiesNumeric.constantFunc(1)
    functions.append (fn3)


    #"""
    ### indirect setting of transportBCs by spatial selection of vertices
    transportBC_merged = []
    transportIC_merged = []
    ### selecting vertices on the left surface
    leftFaceBC = np.array([0,-1])
    boundA = np.array(  [-1e-4 , -cylinderRad*2, -cylinderRad*2] )
    boundB = np.array(  [ 1e-4 , cylinderRad*2, cylinderRad*2]  )
    leftFace = utilitiesGeom.returnSelectedPts(boundA, boundB, vor.vertices)


    ### selecting vertices on the right surface
    rightFaceBC = np.array([2,-1])
    boundA = np.array(  [cylinderHeight-1e-4 , -cylinderRad*2, -cylinderRad*2] )
    boundB = np.array(  [cylinderHeight+1e-4 , cylinderRad*2, cylinderRad*2]  )
    rightFace = utilitiesGeom.returnSelectedPts(boundA, boundB, vor.vertices)


    #setting of transport bc by using rigid plates
    rigidPlatesTrspt = []
    govNodesTrspt = []
    govNodesTrsptBC = []

    trsptLeftRigidPlate = utilitiesMech.RigidPlate(-len(govNodesTrspt)-1, 3, None, directIdcs = True)
    trsptLeftRigidPlate.setDirectNodes(leftFace)
    trsptLeftRigidPlateMechBC = np.array([0,-1])
    rigidPlatesTrspt.append(trsptLeftRigidPlate)
    govNodesTrspt.append(np.array([ 0, 0, 0]))
    govNodesTrsptBC.append(utilitiesMech.transportBC(govNodesTrspt[-1], trsptLeftRigidPlateMechBC))

    trsptRightRigidPlate = utilitiesMech.RigidPlate(-len(govNodesTrspt)-1, 3, None, directIdcs = True)
    trsptRightRigidPlate.setDirectNodes(rightFace)
    trsptRightRigidPlateMechBC = np.array([2,-1])
    rigidPlatesTrspt.append(trsptRightRigidPlate)
    govNodesTrspt.append(np.array([ -1, -1, -1]))
    govNodesTrsptBC.append(utilitiesMech.transportBC(govNodesTrspt[-1], trsptRightRigidPlateMechBC))

    return node_coords, mechBC_merged, transportBC_merged,  govNodes, govNodesMechBC, rigidPlates, vor, [], functions, rigidPlatesTrspt, govNodesTrspt, govNodesTrsptBC, radii

def create2dPeriodicShear(maxLim, minDist, trials, powerTes ):
    print('Creating 2d periodic rectangle, shear loaded.')
    ### sampling of nodes
    ### direct setting of mechanicalBCs
    node_coords, mechBC_merged, mechInitC_merged, radii = asssemble2dPeriodicShear(maxLim, minDist, trials, powerTes );

    #print ('Conducting Voronoi tesselation...', end ='')
    #vor = Voronoi(node_coords)
    #regions, vertices, polygons, areas, centroids, points = voronoi.voronoi_2d(vor, maxLim)
    #print('done.')

    if powerTes == False:
        vor = Voronoi(node_coords)
        regions, vertices, polygons, areas, centroids, points = voronoi.voronoi_2d(vor, maxLim)
    else:
        vor, regions, vertices, polygons, areas, centroids, points = utilitiesNumeric.runPowerPlain(node_coords, radii, 2, maxLim)
    print('done.')


    # if SHOW_PLOT:
    #     fig = voronoi_plot_2d(vor, show_vertices=False, line_colors='orange',  line_width=2, line_alpha=0.6, point_size=2)
    #     plt.show()

    ########################################################################
    functions = []
    #### Defining functions
    #0 constant zero
    fn = utilitiesNumeric.constantFunc(0)
    functions.append (fn)

    #1 loading function, single force top right, bilinear
    fn1 = utilitiesNumeric.sawToothConstFunc(value = -5e-4, period = 4)
    functions.append (fn1)

    mechIC_merged = []

    ########################################################################
    ### indirect setting of transportBCs by spatial selection of vertices
    transportBC_merged = []
    transportIC_merged = []
    ### selecting vertices on the left surface
    leftFaceBC = np.array([2,-1])
    #leftFaceIC = 25.6
    boundA = np.array(  [-1e-8 , maxLim[1]/10*8] )
    boundB = np.array(  [ 1e-8 , maxLim[1]]  )
    leftFace = utilitiesGeom.returnSelectedPts(boundA, boundB, vor.vertices)
    #print(leftFace)
    for i in range (len(leftFace)):
        trsBC = utilitiesMech.transportBC(leftFace[i], leftFaceBC)
        transportBC_merged.append(trsBC)
        #trsIC = utilitiesGeom.transportIC(leftFace[i], leftFaceIC)
        #transportIC_merged.append(trsIC)

    ### selecting vertices on the right surface
    rightFaceBC = np.array([3,-1])
    boundA = np.array(  [maxLim[0] - 1e-8 , - 1e-8] )
    boundB = np.array(  [maxLim[0] + 1e-8 , maxLim[1]/10*2 ] )
    rightFace = utilitiesGeom.returnSelectedPts(boundA, boundB, vor.vertices)
    #print(rightFace)
    for i in range (len(rightFace)):
        trsBC = utilitiesMech.transportBC(rightFace[i], rightFaceBC)
        transportBC_merged.append(trsBC)

    return node_coords, mechBC_merged, mechIC_merged, transportBC_merged, transportIC_merged, vor, areas, functions, radii



def create3dPeriodicShear(maxLim, minDist, trials, powerTes ):
    print('Creating 3d periodic rectangle, shear loaded.')
    ### sampling of nodes
    ### direct setting of mechanicalBCs
    node_coords, mechBC_merged, mechInitC_merged, radii = asssemble3dPeriodicRectangle(maxLim, minDist, trials, powerTes );

    print ('Conducting Voronoi tesselation...', end ='')
    if powerTes == False:
        vor = Voronoi(node_coords)
        volumes = voronoi.voronoi_3d(vor, maxLim)
    else:
        vor, volumes = utilitiesNumeric.runPowerPlain(node_coords, radii, 3, maxLim)
    print('done.')

    ########################################################################
    functions = []
    #### Defining functions
    #0 constant zero
    fn = utilitiesNumeric.constantFunc(0)
    functions.append (fn)

    #1 loading function, single force top right, bilinear
    fn1 = utilitiesNumeric.sawToothConstFunc(value = -5e-4, period = 4)
    #functions.append (fn1)

    func2 = []
    func2.append( np.array([0,0]) )
    func2.append( np.array([1, -1e-3]) )
    fn2 = utilitiesNumeric.generalFunc(func2)
    functions.append (fn2)

    mechIC_merged = []

    ########################################################################
    ### indirect setting of transportBCs by spatial selection of vertices
    transportBC_merged = []
    transportIC_merged = []
    ### selecting vertices on the left surface
    """
    leftFaceBC = np.array([2,-1])
    #leftFaceIC = 25.6
    boundA = np.array(  [-1e-8 , maxLim[1]/10*8] )
    boundB = np.array(  [ 1e-8 , maxLim[1]]  )
    leftFace = utilitiesGeom.returnSelectedPts(boundA, boundB, vor.vertices)
    #print(leftFace)
    for i in range (len(leftFace)):
        trsBC = utilitiesMech.transportBC(leftFace[i], leftFaceBC)
        transportBC_merged.append(trsBC)
        #trsIC = utilitiesGeom.transportIC(leftFace[i], leftFaceIC)
        #transportIC_merged.append(trsIC)

    ### selecting vertices on the right surface
    rightFaceBC = np.array([3,-1])
    boundA = np.array(  [maxLim[0] - 1e-8 , - 1e-8] )
    boundB = np.array(  [maxLim[0] + 1e-8 , maxLim[1]/10*2 ] )
    rightFace = utilitiesGeom.returnSelectedPts(boundA, boundB, vor.vertices)
    #print(rightFace)
    for i in range (len(rightFace)):
        trsBC = utilitiesMech.transportBC(rightFace[i], rightFaceBC)
        transportBC_merged.append(trsBC)

    #return node_coords, mechBC_merged, mechIC_merged, transportBC_merged, transportIC_merged, vor, areas, functions, nodePositions, coupledNodes, mirtype
    """

    return node_coords, mechBC_merged, mechIC_merged, transportBC_merged, transportIC_merged, vor, volumes, functions, radii


def create2dCoupledRVE(maxLim, minDist, trials, powerTes ):
    print('Creating 2d periodic RVE.')
    ### sampling of nodes
    ### direct setting of mechanicalBCs
    node_coords, mechBC_merged, mechInitC_merged, radii = asssemble2dPeriodicShear(maxLim, minDist, trials, powerTes );

    print ('Conducting Voronoi tesselation...', end ='')
    vor = Voronoi(node_coords)
    regions, vertices, polygons, areas, centroids, points = voronoi.voronoi_2d(vor, maxLim)
    print('done.')


    # if SHOW_PLOT:
    #     fig = voronoi_plot_2d(vor, show_vertices=False, line_colors='orange',  line_width=2, line_alpha=0.6, point_size=2)
    #     plt.show()

    ########################################################################
    functions = []
    #### Defining functions
    #0 constant zero
    fn = utilitiesNumeric.constantFunc(0)
    functions.append (fn)

    #1 loading function, single force top right, bilinear
    fn1 = utilitiesNumeric.sawToothConstFunc(value = -5e-4, period = 4)
    functions.append (fn1)

    mechIC_merged = []

    ########################################################################
    ### indirect setting of transportBCs by spatial selection of vertices
    transportBC_merged = []
    transportIC_merged = []
    ### selecting vertices on the left surface
    leftFaceBC = np.array([2,-1])
    #leftFaceIC = 25.6
    boundA = np.array(  [-1e-8 , maxLim[1]/10*8] )
    boundB = np.array(  [ 1e-8 , maxLim[1]]  )
    leftFace = utilitiesGeom.returnSelectedPts(boundA, boundB, vor.vertices)
    #print(leftFace)
    for i in range (len(leftFace)):
        trsBC = utilitiesMech.transportBC(leftFace[i], leftFaceBC)
        transportBC_merged.append(trsBC)
        #trsIC = utilitiesGeom.transportIC(leftFace[i], leftFaceIC)
        #transportIC_merged.append(trsIC)

    ### selecting vertices on the right surface
    rightFaceBC = np.array([3,-1])
    boundA = np.array(  [maxLim[0] - 1e-8 , - 1e-8] )
    boundB = np.array(  [maxLim[0] + 1e-8 , maxLim[1]/10*2 ] )
    rightFace = utilitiesGeom.returnSelectedPts(boundA, boundB, vor.vertices)
    #print(rightFace)
    for i in range (len(rightFace)):
        trsBC = utilitiesMech.transportBC(rightFace[i], rightFaceBC)
        transportBC_merged.append(trsBC)

    print("BOUNDARY CONDITIONS",transportBC_merged)

    return node_coords, mechBC_merged, mechIC_merged, transportBC_merged, transportIC_merged, vor, areas, functions, radii


def create2dCircRVE(maxLim, minDist, trials, powerTes):
    print('Creating 2d circular periodic RVE.')

    ### sampling of nodes
    all_node_coords, all_node_coords_polar, mechBC_merged, mechInitC_merged, radii  = asssemble2dCircRVE(maxLim, minDist, trials, powerTes )

    if powerTes==False:
        print('Conducting Voronoi tesselation...', end = '')
        vor = utilitiesNumeric.Voronoi (all_node_coords)
        regions, vertices, polygons, areas, centroids, points = voronoi.voronoi_2d(vor, maxLim)
    else:
        print('Conducting Power tesselation...', end = '')
        vor, regions, vertices, polygons, areas, centroids, points = utilitiesNumeric.runPowerPlain(all_node_coords, radii, 2, [-maxLim/2, maxLim/2])
    print('done.')

    
    # fig = voronoi_plot_2d(vor, show_vertices=False, line_colors='orange',  line_width=2, line_alpha=0.6, point_size=2)
    # plt.show()

    ########################################################################
    functions = []
    #### Defining functions
    #0 constant zero
    fn = utilitiesNumeric.constantFunc(0)
    functions.append (fn)

    #1 loading function, single force top right, bilinear
    fn1 = utilitiesNumeric.sawToothConstFunc(value = -5e-4, period = 4)
    functions.append (fn1)

    mechIC_merged = []

    ########################################################################
    ### transport not needed so far
    transportBC_merged = []
    transportIC_merged = []

    return all_node_coords, all_node_coords_polar, mechBC_merged, mechIC_merged, transportBC_merged, transportIC_merged, vor, areas, functions, radii


def assembleTwoNodeSpringTest (maxLim, idt):
    node_coords = []
    mechBC_merged = []

    dim = 2

    govNodes = None
    govNodesMechBC = None
    rigidPlates = None
    """
    indentRP = 1e-1
    rightRigidPlateMechBC = np.array([1,2,0, -1,-1,-1])
    rightRigidPlate = utilitiesMech.RigidPlate(-2, 2,
    np.array([ maxLim[0] - idt -indentRP,
     maxLim[0] - idt+indentRP,
     -1000,
      1000 ]))
    rigidPlates.append(rightRigidPlate)
    govNodes.append(np.array([ maxLim[0] - idt+indentRP, 0]))
    govNodesMechBC.append(utilitiesMech.mechanicalBC(dim, -2, rightRigidPlateMechBC))
    #"""

    nodeA = np.array ( [ 0 + idt , maxLim[1]/2 ] )
    pointGenerators.generateSingleNode(nodeA, dim, node_coords)
    mBC = utilitiesMech.mechanicalBC(dim, 0, np.array([0,0,0, -1,-1,-1]))
    mechBC_merged.append(mBC)

    nodeB = np.array ( [ maxLim[0] - idt, maxLim[1]/2 ] )
    pointGenerators.generateSingleNode(nodeB, dim, node_coords)
    mBC = utilitiesMech.mechanicalBC(dim, 1, np.array([-1,-1,0, 1,1,-1]))
    mechBC_merged.append(mBC)

    return node_coords, mechBC_merged,  govNodes, govNodesMechBC, rigidPlates


def assemble2D_singleContact (maxLim, idt):
    node_coords = []
    mechBC_merged = []

    dim = 2
    govNodes = []
    govNodesMechBC = []
    rigidPlates = []

    nodeA = np.array ( [ 0 + idt , maxLim[1]/2 ] )
    pointGenerators.generateSingleNode(nodeA, dim, node_coords)
    mBC = utilitiesMech.mechanicalBC(dim, 0, np.array([0,0,0, -1,-1,-1]))
    mechBC_merged.append(mBC)

    nodeB = np.array ( [ maxLim[0] - idt, maxLim[1]/2 ] )
    pointGenerators.generateSingleNode(nodeB, dim, node_coords)
    mBC = utilitiesMech.mechanicalBC(dim, 1, np.array([-1,-1,0, 1,1,-1]))
    mechBC_merged.append(mBC)

    return node_coords, mechBC_merged,  govNodes, govNodesMechBC, rigidPlates

def assemble3D_singleContact (maxLim, idt):
    node_coords = []
    mechBC_merged = []

    dim = 3
    govNodes = []
    govNodesMechBC = []
    rigidPlates = []

    nodeA = np.array ( [ 0 + idt , maxLim[1]/2,maxLim[2]/2 ] )
    pointGenerators.generateSingleNode(nodeA, dim, node_coords)
    mBC = utilitiesMech.mechanicalBC(dim, 0, np.array([0,0,0, 0,0,0, -1,-1,-1, -1,-1,-1]))
    mechBC_merged.append(mBC)

    nodeB = np.array ( [ maxLim[0] - idt, maxLim[1]/2,maxLim[2]/2 ] )
    pointGenerators.generateSingleNode(nodeB, dim, node_coords)
    mBC = utilitiesMech.mechanicalBC(dim, 1, np.array([1,0,0, 0,0,0, -1,-1,-1, -1,-1,-1]))
    mechBC_merged.append(mBC)

    return node_coords, mechBC_merged,  govNodes, govNodesMechBC, rigidPlates


def create2dDogBone(minDist, trials, D=1.0, excentricity = 50, symmetric=False, edgeMinDistCoef=1.0, roughDogBone=0, roughEdgeDogbone = 0, roughMinDistCoef=1, interLayerThickness=2, powerTes = False, weakboundary = 0):
    print('Creating 2d dog bone....')
    #


    node_coords_all, node_indices_dogbone, mechBC_merged, mechInitC_merged, node_count, govNodes, govNodesMechBC, rigidPlates, radii  = assemble2dDogBone(D, minDist, trials, excentricity = excentricity, symmetric = symmetric, edgeMinDistCoef=edgeMinDistCoef, roughDogBone=roughDogBone, roughEdgeDogbone=roughEdgeDogbone, roughMinDistCoef=roughMinDistCoef, interLayerThickness=interLayerThickness, powerTes = powerTes, weakboundary = weakboundary);

    node_coords_all = np.asarray(node_coords_all)
    radii = np.asarray(radii)

    """
    fig, ax = plt.subplots()
    ax.scatter(node_coords_all[:,0], node_coords_all[:,1])
    ax.scatter(node_coords_all[node_indices_dogbone,0], node_coords_all[node_indices_dogbone,1])
    plt.show()
    #"""


    if powerTes:
        print('Conducting Power tesselation...', end = '')
        # vor = utilitiesNumeric.runMirroredVoronoiDogBone(node_coords_all, 2, D)
        vor = utilitiesNumeric.runMirroredPowerDogBone(node_coords_all, 2, D, radii=radii)
        print('done.')

    else:
        print('Conducting Voronoi tesselation...', end = '')
        vor = utilitiesNumeric.runMirroredVoronoiDogBone(node_coords_all, 2, D)
        print('done.')


    node_coords = np.copy(node_coords_all)
    # areas = []
    # for i in range (node_count): areas.append(0)
    # areas = np.asarray(areas)
    areas = np.zeros(node_count)

    # if SHOW_PLOT:
    #     fig = voronoi_plot_2d(vor, show_vertices=True, line_colors='orange',line_width=2, line_alpha=0.6, point_size=2)
    #     plt.show()

    ########################################################################
    functions = []
    #### Defining functions
    #0 constant zero
    fn = utilitiesNumeric.constantFunc(0)
    functions.append (fn)

    #1 loading function
    func1 = []
    func1.append( np.array([0,0]) )
    func1.append( np.array([1, -1e-3]) )
    fn1 = utilitiesNumeric.generalFunc(func1)
    functions.append (fn1)

    ########################################################################
    ### indirect setting of transportBCs by spatial selection of vertices
    transportBC_merged = []
    transportIC_merged = []


    return node_coords, mechBC_merged, mechInitC_merged, transportBC_merged, transportIC_merged, vor, areas, functions,  govNodes, govNodesMechBC, rigidPlates, node_indices_dogbone



def create2dDogBoneStrip(minDist, D=1.0, excentricity = 50, randomStrip = 0):
    print('Creating 2d dog bone strip....')
    #
    dim=2

    node_coords, mechBC_merged, mechInitC_merged, node_count, govNodes, govNodesMechBC, rigidPlates  = assemble2dDogBoneStrip(D, minDist, excentricity = excentricity, randomStrip = randomStrip);

    node_coords = np.asarray(node_coords)


    fig, ax = plt.subplots()
    ax.scatter(node_coords[:,0], node_coords[:,1])
    plt.show()


    print('Conducting Voronoi tesselation...', end = '')
    maxLim = np.array((0.8*D, minDist))
    vor, regions, vertices, polygons, areas, centroids, points = utilitiesNumeric.runMirroredVoronoi (node_coords, dim, maxLim)
    print('done.')

    node_coords = np.copy(node_coords)
    areas = []
    for i in range (node_count): areas.append(0)
    areas = np.asarray(areas)


    # if SHOW_PLOT:
    #     fig = voronoi_plot_2d(vor, show_vertices=True, line_colors='orange',line_width=2, line_alpha=0.6, point_size=2)
    #     plt.show()

    ########################################################################
    functions = []
    #### Defining functions
    #0 constant zero
    fn = utilitiesNumeric.constantFunc(0)
    functions.append (fn)

    #1 loading function
    func1 = []
    func1.append( np.array([0,0]) )
    func1.append( np.array([1, -1e-3]) )
    fn1 = utilitiesNumeric.generalFunc(func1)
    functions.append (fn1)

    ########################################################################
    ### indirect setting of transportBCs by spatial selection of vertices
    transportBC_merged = []
    transportIC_merged = []


    return node_coords, mechBC_merged, mechInitC_merged, transportBC_merged, transportIC_merged, vor, areas, functions,  govNodes, govNodesMechBC, rigidPlates



def create2dDogBoneBand(maxLim, minDist,  roughMinDistCoef=1, elasticHeightCoef=1):
    print('Creating 2d dog bone strip....')
    #
    dim=2


    node_coords, mechBC_merged, mechInitC_merged, node_count, govNodes, govNodesMechBC, rigidPlates  = assemble2dDogBoneBand(maxLim, minDist, roughMinDistCoef=roughMinDistCoef, elasticHeightCoef=elasticHeightCoef);



    """
    node_coords = np.asarray(node_coords)
    fig, ax = plt.subplots()
    ax.scatter(node_coords[:,0], node_coords[:,1])
    plt.show()
    """

    print('Conducting Voronoi tesselation...', end = '')
    #maxLim = np.array((0.8*D, minDist))
    vor, regions, vertices, polygons, areas, centroids, points = utilitiesNumeric.runMirroredVoronoi (node_coords, dim, maxLim)
    print('done.')

    node_coords = np.copy(node_coords)
    areas = []
    for i in range (node_count): areas.append(0)
    areas = np.asarray(areas)


    # if SHOW_PLOT:
    #     fig = voronoi_plot_2d(vor, show_vertices=True, line_colors='orange',line_width=2, line_alpha=0.6, point_size=2)
    #     plt.show()

    ########################################################################
    functions = []
    #### Defining functions
    #0 constant zero
    fn = utilitiesNumeric.constantFunc(0)
    functions.append (fn)

    #1 loading function
    func1 = []
    func1.append( np.array([0,0]) )
    func1.append( np.array([1, 1e-3]) )
    fn1 = utilitiesNumeric.generalFunc(func1)
    functions.append (fn1)

    ########################################################################
    ### indirect setting of transportBCs by spatial selection of vertices
    transportBC_merged = []
    transportIC_merged = []


    return node_coords, mechBC_merged, mechInitC_merged, transportBC_merged, transportIC_merged, vor, areas, functions,  govNodes, govNodesMechBC, rigidPlates


def create3dDogBone(minDist, trials,thickness, D=1.0,D0=None,H=None,H0=None, excentricity_X = 20, excentricity_Z = 0, symmetric=False):
    print('Creating 3sd dog bone....')
    #
    node_coords_all, node_indices_dogbone, mechBC_merged, mechInitC_merged, node_count,govNodes, govNodesMechBC, rigidPlates,radii  = assemble3dDogBone(D, minDist, trials,thickness,D0=D0,H=H,H0=H0, excentricity_X= excentricity_X, excentricity_Z= excentricity_Z, symmetric = symmetric)

    node_coords_all = np.asarray(node_coords_all)
    """
    if SHOW_PLOT:
        fig = plt.figure()
        ax = fig.add_subplot(111, projection='3d')
        ax.scatter(node_coords_all[:,0], node_coords_all[:,1], node_coords_all[:,2], c = 'b', marker='o')
        plt.show()
    """
    print('Conducting Voronoi tesselation...', end = '')
    vor = utilitiesNumeric.runMirroredPowerDogBone_old(node_coords_all, 3, D,H=H or 6/4*D, thickness = thickness,radii=radii)
    print('done.')

    node_coords_all = node_coords_all[0:node_count]
    areas = []
    for i in range (node_count): areas.append(0)
    areas = np.asarray(areas)


    ########################################################################
    functions = []
    #0 constant zero
    fn = utilitiesNumeric.constantFunc(0)
    functions.append (fn)

    #1 loading function
    func1 = []
    func1.append( np.array([0,0]) )
    func1.append( np.array([1, -1e-4]) )
    fn1 = utilitiesNumeric.generalFunc(func1)
    functions.append (fn1)

    ########################################################################
    ### indirect setting of transportBCs by spatial selection of vertices
    transportBC_merged = []
    transportIC_merged = []


    return node_coords_all, mechBC_merged, mechInitC_merged, transportBC_merged, transportIC_merged, vor, areas, functions, govNodes, govNodesMechBC, rigidPlates, node_indices_dogbone








def assembleDiamondTest (maxLim, idtW, idtH):
    node_coords = []
    mechBC_merged = []
    dim = 2

    #left mid
    nodeA = np.array ( [ 0 + idtW , maxLim[1]/2 ] )
    nodeAmechBC = np.array([-1,-1,-1 , -1,-1 ,-1])
    pointGenerators.generateSingleNode(nodeA, 2, node_coords)
    mBC = utilitiesMech.mechanicalBC(dim, 0, nodeAmechBC)
    mechBC_merged.append(mBC)

    #right mid
    nodeB = np.array ( [ maxLim[0] - idtW, maxLim[1]/2 ] )
    nodeBmechBC = np.array([-1,-1,-1 , -1,-1,-1])
    pointGenerators.generateSingleNode(nodeB, 2, node_coords)
    mBC = utilitiesMech.mechanicalBC(dim, 1, nodeBmechBC)
    mechBC_merged.append(mBC)

    #top
    nodeC = np.array ( [ maxLim[0]/2, maxLim[1] -idtH] )
    nodeCmechBC = np.array([0, -1 , 0 , -1, 1, -1])
    pointGenerators.generateSingleNode(nodeC, 2, node_coords)
    mBC = utilitiesMech.mechanicalBC(dim, 2, nodeCmechBC)
    mechBC_merged.append(mBC)

    #bottom
    nodeD = np.array ( [ maxLim[0]/2, idtH] )
    nodeDmechBC = np.array([0,0,0,  -1, -1, -1])
    pointGenerators.generateSingleNode(nodeD, 2, node_coords)
    mBC = utilitiesMech.mechanicalBC(dim, 3, nodeDmechBC)
    mechBC_merged.append(mBC)


    return node_coords,  mechBC_merged






def create3dSSBeamUnifLoad(maxLim, minDist, trials, notch = -1, loadWidth = 1, fracZoneWidth = 0.15, orthogonalFracZone = False, notchWidth = -1, activeMechanics = True, activeTransport=False, coupled=False, node_coords_init=None, specifiedNodes=[], roughMinDistCoef=1, supportDivision=10,gapWidth=0,blank=0,powerTes=False,supportWidth=None,span=0,gradientZoneWidth=0):

    radii = []
    print('Creating 3d simply supported beam, uniform load.')
    #govNodes, rigidPlates
    if not blank:
        node_coords, mechBC_merged, mechInitC_merged, notches, govNodes, govNodesMechBC, rigidPlates,radii  = assemble3DSSBeamBending(maxLim, minDist, trials, notch, loadWidth, fracZoneWidth=fracZoneWidth, orthogonalFracZone=orthogonalFracZone, notchWidth = notchWidth, coupled=coupled, node_coords_init=node_coords_init, specifiedNodes=specifiedNodes, roughMinDistCoef=roughMinDistCoef, supportDivision=supportDivision,gapWidth=gapWidth,blank=blank,powerTes=powerTes,supportWidth=supportWidth,span=span,gradientZoneWidth=gradientZoneWidth);
        node_coords = np.asarray(node_coords)
    else:
        node_coords, mechBC_merged, mechInitC_merged, notches, govNodes, govNodesMechBC, rigidPlates,radii  = assemble3DSSBeamBendingBlank(maxLim, minDist, trials, notch, loadWidth, fracZoneWidth=fracZoneWidth, orthogonalFracZone=orthogonalFracZone, notchWidth = notchWidth, coupled=coupled, node_coords_init=node_coords_init, specifiedNodes=specifiedNodes, roughMinDistCoef=roughMinDistCoef, supportDivision=supportDivision,gapWidth=gapWidth,blank=blank,powerTes=powerTes,supportWidth=supportWidth,span=span,gradientZoneWidth=gradientZoneWidth);
        node_coords = np.asarray(node_coords)
        
    for i in node_coords:
        if i[1]<0:
            print(i)
            return

    """
    if SHOW_PLOT:
        fig = plt.figure()
        ax = Axes3D(fig)
        ax.auto_scale_xyz([0, maxLim[0]], [0, maxLim[1]], [0, maxLim[2]])
        ax.scatter(node_coords[:,0], node_coords[:,1], node_coords[:,2])
        plt.show()
    """
    print('Conducting Voronoi tesselation...', end = '')


    if notch == 0:
        notch = None
    else:
        notch=[notch,notchWidth]

    if powerTes==False:
        vor, volumes = utilitiesNumeric.runMirroredVoronoi (node_coords, 3, maxLim,notch=notch)
        radii=[]
    else:
        vor, volumes = utilitiesNumeric.runMirroredPower(node_coords, radii, 3, maxLim,notch=notch)


    ########################################################################
    functions = []
    #### Defining functions
    #0 constant zero
    fn = utilitiesNumeric.constantFunc(0)
    functions.append (fn)

    if gapWidth >0:
        func1 = []
        func1.append( np.array([0,0]) )
        func1.append( np.array([0.1,0]) )
        func1.append( np.array([1, -0.35e-3]) )
        fn1 = utilitiesNumeric.generalFunc(func1)
        functions.append (fn1)

        func1 = []
        func1.append( np.array([0,0]) )
        func1.append( np.array([0.1,100]) )
        func1.append( np.array([1,  100]) )
        fn1 = utilitiesNumeric.generalFunc(func1)
        functions.append (fn1)
    else:
        func1 = []
        func1.append( np.array([0,0]) )
        func1.append( np.array([1, -1e-4]) )
        fn1 = utilitiesNumeric.generalFunc(func1)
        functions.append (fn1)

        fn2 = utilitiesNumeric.constantFunc(100)
        functions.append (fn2)

    transportBC_merged = []
    govNodesTrspt = []
    govNodesTrsptBC = []
    rigidPlatesTrspt = []

    if coupled == True:
        ### selecting vertices on the bottom surface
        botFaceBC = np.array([2,-1])
        boundA = np.array(  [-1e-4 , -1e-5, -1e-5] )
        boundB = np.array(  [ maxLim[0]+1e-4 , 1e-5, maxLim[2]+1e-5]  )
        botFace = utilitiesGeom.returnSelectedPts(boundA, boundB, vor.vertices)

        ### selecting vertices inside notch
        try:
            notchFaceBC = np.array([2,-1])
            boundA = np.array(  [ maxLim[0]/2-notchWidth/3 , -1e-5, -1e-5] )
            boundB = np.array(  [ maxLim[0]/2+notchWidth/3 , maxLim[1]*notch+minDist*0.9, maxLim[2]+1e-5]  )
            notchFace = utilitiesGeom.returnSelectedPts(boundA, boundB, vor.vertices)
        except:
            notchFace = []

        for i in range (len(notchFace)):
            if not notchFace[i] in botFace:
                botFace = np.hstack((botFace, notchFace[i]))
            #else:
                #print('Point is already a part of the bottom surface.')

        trsptBottomRigidPlate = utilitiesMech.RigidPlate(-len(govNodesTrspt)-1, 3, None, directIdcs = True)
        trsptBottomRigidPlate.setDirectNodes(botFace)
        trsptBottomRigidPlateMechBC = botFaceBC
        rigidPlatesTrspt.append(trsptBottomRigidPlate)
        govNodesTrspt.append(np.array([ 0, 0,0]))
        govNodesTrsptBC.append(utilitiesMech.transportBC(govNodesTrspt[-1], trsptBottomRigidPlateMechBC))

        #for i in range (len(botFace)):
        #    trsBC = utilitiesMech.transportBC(botFace[i], botFaceBC)
        #    transportBC_merged.append(trsBC)

        ### selecting vertices on the top surface
        topFaceBC = np.array([0,-1])
        boundA = np.array(  [-1e-4 , maxLim[1]-1e-5, -1e-5] )
        boundB = np.array(  [ maxLim[0]+1e-4 , maxLim[1]+1e-5, maxLim[2]+1e-5]  )
        topFace = utilitiesGeom.returnSelectedPts(boundA, boundB, vor.vertices)

        trsptTopRigidPlate = utilitiesMech.RigidPlate(-len(govNodesTrspt)-1, 2, None, directIdcs = True)
        trsptTopRigidPlate.setDirectNodes(topFace)
        trsptTopRigidPlateMechBC = topFaceBC
        rigidPlatesTrspt.append(trsptTopRigidPlate)
        govNodesTrspt.append(np.array([ -1, 1,0]))
        govNodesTrsptBC.append(utilitiesMech.transportBC(govNodesTrspt[-1], trsptTopRigidPlateMechBC))

        #for i in range (len(topFace)):
        #    trsBC = utilitiesMech.transportBC(topFace[i], topFaceBC)
        #    transportBC_merged.append(trsBC)



    #return node_coords, mechBC_merged, mechInitC_merged,  vor, volumes, functions, notches, govNodes, govNodesMechBC, rigidPlates, transportBC_merged, transportIC_merged
    #return node_coords, mechBC_merged, transportBC_merged, govNodes, govNodesMechBC, rigidPlates, vor, areas, functions, rigidPlatesTrspt, govNodesTrspt, govNodesTrsptBC
    #return node_coords, mechBC_merged, mechInitC_merged,  vor, volumes, functions, notches, govNodes, govNodesMechBC, rigidPlates, rigidPlatesTrspt, govNodesTrspt, govNodesTrsptBC, transportBC_merged
    return node_coords, mechBC_merged, transportBC_merged,  govNodes, govNodesMechBC, rigidPlates, vor, [], functions, rigidPlatesTrspt, govNodesTrspt, govNodesTrsptBC, radii,notches

def create3dCube(maxLim, minDist, trials, powerTes, coupled=False, node_coords_init=None, periodic = False ):
    print('Creating 3d cube. Power tesselation: %s' %powerTes)
    #govNodes, rigidPlates
    dim = 3
    node_coords, mechBC_merged, mechInitC_merged,  govNodes, govNodesMechBC, rigidPlates,radii  = assemble3Dcube(maxLim, minDist, trials, powerTes, coupled=coupled, node_coords_init=node_coords_init, periodic = periodic);
    node_coords = np.asarray(node_coords)
    """
    if SHOW_PLOT:
        fig = plt.figure()
        ax = Axes3D(fig)
        ax.auto_scale_xyz([0, maxLim[0]], [0, maxLim[1]], [0, maxLim[2]])
        ax.scatter(node_coords[:,0], node_coords[:,1], node_coords[:,2])
        plt.show()
    """

    if powerTes==False:
        print('Conducting Voronoi tesselation...', end = '')
        vor, volumes = utilitiesNumeric.runMirroredVoronoi (node_coords, 3, maxLim)
    else:
        print('Conducting Power tesselation...', end = '')
        vor, volumes = utilitiesNumeric.runMirroredPower(node_coords, radii, 3, maxLim)
    print('done.')


    ########################################################################
    functions = []
    #### Defining functions
    #0 constant zero
    fn = utilitiesNumeric.constantFunc(0)
    functions.append (fn)

    #1 loading function, top surface,
    func1 = []
    func1.append( np.array([0,0]) )
    func1.append( np.array([1, 1e-4]) )
    fn1 = utilitiesNumeric.generalFunc(func1)
    functions.append (fn1)

    fn2 = utilitiesNumeric.constantFunc(100)
    functions.append (fn2)

    transportBC_merged = []
    transportIC_merged = []
    if coupled == True:

        ### selecting vertices on the bottom surface
        botFaceBC = np.array([2,-1])
        boundA = np.array(  [-1e-4 , -1e-5, -1e-5] )
        boundB = np.array(  [ maxLim[0]+1e-4 , 1e-5, maxLim[2]+1e-5]  )
        botFace = utilitiesGeom.returnSelectedPts(boundA, boundB, vor.vertices)

        for i in range (len(botFace)):
            trsBC = utilitiesMech.transportBC(botFace[i], botFaceBC)
            transportBC_merged.append(trsBC)

        ### selecting vertices on the top surface
        topFaceBC = np.array([0,-1])
        boundA = np.array(  [-1e-4 , maxLim[1]-1e-5, -1e-5] )
        boundB = np.array(  [ maxLim[0]+1e-4 , maxLim[1]+1e-5, maxLim[2]+1e-5]  )
        topFace = utilitiesGeom.returnSelectedPts(boundA, boundB, vor.vertices)

        for i in range (len(topFace)):
            trsBC = utilitiesMech.transportBC(topFace[i], topFaceBC)
            transportBC_merged.append(trsBC)

    return node_coords, mechBC_merged, mechInitC_merged,  vor, volumes, functions, govNodes, govNodesMechBC, rigidPlates, transportBC_merged, transportIC_merged


def create3dBalbet(maxLim, minDist, trials, powerTes, shotRadius=0.02, shotGradientRadius=0.01, roughMinDistCoef=2,coupled=False, node_coords_init=None ):
    print('Creating 3d balbet.')
    #govNodes, rigidPlates
    dim = 3
    node_coords, mechBC_merged, mechInitC_merged,  govNodes, govNodesMechBC, rigidPlates,radii  = assemble3DBalbet(maxLim, minDist, trials, powerTes, shotRadius=shotRadius, shotGradientRadius=shotGradientRadius, roughMinDistCoef=roughMinDistCoef, coupled=coupled, node_coords_init=node_coords_init);
    node_coords = np.asarray(node_coords)

    print('Conducting Voronoi tesselation...', end = '')
    vor, volumes = utilitiesNumeric.runMirroredVoronoi (node_coords, 3, maxLim)
    print('done.')

    ########################################################################
    functions = []
    #### Defining functions
    #0 constant zero
    fn = utilitiesNumeric.constantFunc(0)
    functions.append (fn)

    #shot
    func1 = []
    func1.append( np.array([0,0]) )
    func1.append( np.array([1, 1e-3]) )
    fn1 = utilitiesNumeric.generalFunc(func1)
    functions.append (fn1)


    transportBC_merged = []
    transportIC_merged = []


    return node_coords, mechBC_merged, mechInitC_merged,  vor, volumes, functions, govNodes, govNodesMechBC, rigidPlates, transportBC_merged, transportIC_merged


def create3dConsolidation(maxLim, minDist, trials, powerTes, coupled=False, node_coords_init=None ):
    print('Creating 3d cube. Power tesselation: %s' %powerTes)
    #govNodes, rigidPlates
    dim = 3
    node_coords, mechBC_merged, mechInitC_merged,  govNodes, govNodesMechBC, rigidPlates, radii  = assemble3Dcube(maxLim, minDist, trials, powerTes, coupled=coupled, node_coords_init=node_coords_init);
    node_coords = np.asarray(node_coords)
    """
    if SHOW_PLOT:
        fig = plt.figure()
        ax = Axes3D(fig)
        ax.auto_scale_xyz([0, maxLim[0]], [0, maxLim[1]], [0, maxLim[2]])
        ax.scatter(node_coords[:,0], node_coords[:,1], node_coords[:,2])
        plt.show()
    """
    print('Conducting Voronoi tesselation...', end = '')
    vor, volumes = utilitiesNumeric.runMirroredVoronoi (node_coords, 3, maxLim)
    print('done.')

    ########################################################################
    functions = []
    #### Defining functions
    #0 constant zero
    fn = utilitiesNumeric.constantFunc(0)
    functions.append (fn)

    #1 loading function, top surface,
    func1 = []
    func1.append( np.array([0,0]) )
    func1.append( np.array([1, -1e-2]) )
    fn1 = utilitiesNumeric.generalFunc(func1)
    functions.append (fn1)

    fn2 = utilitiesNumeric.constantFunc(100)
    functions.append (fn2)

    fn3 = utilitiesNumeric.constantFunc(0)
    functions.append (fn3)

    transportBC_merged = []
    transportIC_merged = []

    """
    ### selecting vertices on the bottom surface
    botFaceBC = np.array([2,-1])
    boundA = np.array(  [-1e-4 , -1e-5, -1e-5] )
    boundB = np.array(  [ maxLim[0]+1e-4 , 1e-5, maxLim[2]+1e-5]  )
    botFace = utilitiesGeom.returnSelectedPts(boundA, boundB, vor.vertices)

    for i in range (len(botFace)):
        trsBC = utilitiesMech.transportBC(botFace[i], botFaceBC)
        transportBC_merged.append(trsBC)
    """

    ### selecting vertices on the right surface
    rightFaceBC = np.array([2,-1])
    boundA = np.array(  [maxLim[0]-1e-4 , -1e-5, -1e-5] )
    boundB = np.array(  [ maxLim[0]+1e-4 , maxLim[1]+1e-5, maxLim[2]+1e-5]  )
    rightFace = utilitiesGeom.returnSelectedPts(boundA, boundB, vor.vertices)
    for i in range (len(rightFace)):
        trsBC = utilitiesMech.transportBC(rightFace[i], rightFaceBC)
        transportBC_merged.append(trsBC)


    #block transversal movements
    FaceBC = np.array([-1,3,3,-1,-1,-1,-1,-1,-1,-1,-1,-1])
    boundA = np.array(  [-1e-4 , -1e-5, -1e-5] )
    boundB = np.array(  [ maxLim[0]+1e-4 , 1e-5, maxLim[2]+1e-5]  )
    faceBot = utilitiesGeom.returnSelectedPts(boundA, boundB, node_coords)
    boundA = np.array(  [-1e-4 , maxLim[1]-1e-5, -1e-5] )
    boundB = np.array(  [ maxLim[0]+1e-4 , maxLim[1]+1e-5, maxLim[2]+1e-5]  )
    faceTop = utilitiesGeom.returnSelectedPts(boundA, boundB, node_coords)
    boundA = np.array(  [-1e-4 , -1e-5, -1e-5] )
    boundB = np.array(  [ maxLim[0]+1e-4 , maxLim[1]+1e-5, 1e-5]  )
    faceFro = utilitiesGeom.returnSelectedPts(boundA, boundB, node_coords)
    boundA = np.array(  [-1e-4 , -1e-5, maxLim[2]-1e-5] )
    boundB = np.array(  [ maxLim[0]+1e-4 , maxLim[1]+1e-5, maxLim[2]+1e-5]  )
    faceBac = utilitiesGeom.returnSelectedPts(boundA, boundB, node_coords)
    face = np.hstack((faceBot,faceTop,faceFro,faceBac))

    for i in range (len(face)):
        mechBC = utilitiesMech.mechanicalBC(3, face[i], FaceBC)
        mechBC_merged.append(mechBC)


    return node_coords, mechBC_merged, mechInitC_merged,  vor, volumes, functions, govNodes, govNodesMechBC, rigidPlates, transportBC_merged, transportIC_merged, radii




def create3dDam(maxLim, minDist, trials, Xtop):
    print('Creating 3d Dam....')
    #
    node_coords, radii, mechBC_merged, mechInitC_merged  = assemble3dDam(maxLim, minDist, trials, Xtop);

    node_coords = np.asarray(node_coords)
    node_count = len(node_coords)
    """
    if SHOW_PLOT:
        fig = plt.figure()
        ax = fig.add_subplot(111, projection='3d')
        ax.scatter(node_coords[:,0], node_coords[:,1], node_coords[:,2], c = 'b', marker='o')
        plt.show()
    """
    print('Conducting Voronoi tesselation...', end = '')
    vor = utilitiesNumeric.runMirroredPowerDam(node_coords, radii, 3, maxLim, Xtop)
    print('done.')

    node_coords = node_coords[0:node_count]
    areas = []
    for i in range (node_count): areas.append(0)
    areas = np.asarray(areas)

    # if SHOW_PLOT:
    #     fig = voronoi_plot_2d(vor, show_vertices=False, line_colors='orange',line_width=2, line_alpha=0.6, point_size=2)
    #     plt.show()

    ########################################################################
    functions = []
    #### Defining functions
    #0 constant zero
    fn = utilitiesNumeric.constantFunc(0)
    functions.append (fn)

    ########################################################################
    ### indirect setting of transportBCs by spatial selection of vertices
    transportBC_merged = []
    transportIC_merged = []
    ### selecting vertices on the left surface
    #right face
    alpha = np.arctan( (maxLim[0] - Xtop)/maxLim[2] )
    planenorm = np.array([np.cos(alpha), 0., np.sin(alpha)])
    planeconst = -planenorm[0]*maxLim[0] - planenorm[1]*maxLim[1]
    rightFace = np.where( abs( np.dot(vor.vertices,planenorm) + planeconst)<1e-8)[0]
    for i in range (len(rightFace)):
        trsBC = utilitiesMech.transportBC(rightFace[i], [0,-1])
        transportBC_merged.append(trsBC)

    leftFaceBC = np.array([2,-1])
    boundA = np.array(  [-1e-8 , 0, 0] )
    boundB = np.array(  [ 1e-8 , maxLim[1], maxLim[2]]  )
    leftFace = utilitiesGeom.returnSelectedPts(boundA, boundB, vor.vertices)
    #print(leftFace)
    for i in range (len(leftFace)):
        trsBC = utilitiesMech.transportBC(leftFace[i], [i+1,-1])
        transportBC_merged.append(trsBC)
        func1 = []
        func1.append( np.array([0,0]) )
        func1.append( np.array([1, (maxLim[2] - vor.vertices[leftFace[i],2])*1000.*9.81 ]) )
        fn1 = utilitiesNumeric.generalFunc(func1)
        functions.append (fn1)


    return node_coords, mechBC_merged, mechInitC_merged, transportBC_merged, transportIC_merged, vor, areas, functions



def create3dReinhardtTension(maxLim, minDist, trials, fracZoneWidth = 0.15 ,roughCoef=1,gapWidth=0):
    print('Creating 3d simply supported beam, uniform load.')
    #govNodes, rigidPlates
    node_coords, mechBC_merged, mechInitC_merged, notchNodes, govNodes, govNodesMechBC, rigidPlates  = assemble3DReinhardtTension(maxLim, minDist, trials, fracZoneWidth=fracZoneWidth,roughCoef=roughCoef,gapWidth=gapWidth);
    node_coords = np.asarray(node_coords)
    """
    if SHOW_PLOT:
        fig = plt.figure()
        ax = Axes3D(fig)
        ax.auto_scale_xyz([0, maxLim[0]], [0, maxLim[1]], [0, maxLim[2]])
        ax.scatter(node_coords[:,0], node_coords[:,1], node_coords[:,2])
        plt.show()
    """
    print('Conducting Voronoi tesselation...', end = '')
    vor, volumes = utilitiesNumeric.runMirroredVoronoi (node_coords, 3, maxLim)
    print('done.')

    ########################################################################
    functions = []
    #### Defining functions
    #0 constant zero
    fn = utilitiesNumeric.constantFunc(0)
    functions.append (fn)

    if gapWidth >0:
        func1 = []
        func1.append( np.array([0,0]) )
        func1.append( np.array([0.1,0]) )
        func1.append( np.array([1, 0.005]) )
        fn1 = utilitiesNumeric.generalFunc(func1)
        functions.append (fn1)

        func1 = []
        func1.append( np.array([0,0]) )
        func1.append( np.array([0.1,-100]) )
        func1.append( np.array([1, -100]) )
        fn1 = utilitiesNumeric.generalFunc(func1)
        functions.append (fn1)

        func1 = []
        func1.append( np.array([0,0]) )
        func1.append( np.array([0.1,100]) )
        func1.append( np.array([1,  100]) )
        fn1 = utilitiesNumeric.generalFunc(func1)
        functions.append (fn1)
    else:
        func1 = []
        func1.append( np.array([0,0]) )
        func1.append( np.array([1, 10000]) )
        fn1 = utilitiesNumeric.generalFunc(func1)
        functions.append (fn1)




    return node_coords, mechBC_merged, mechInitC_merged, vor, volumes, functions, notchNodes, govNodes, govNodesMechBC, rigidPlates




def create3dCantileverBending(maxLim, minDist, trials ):
    ### sampling of nodes
    ### direct setting of mechanicalBCs
    node_coords, mechBC_merged, mechIC_merged  = assemble3dCantileverBending(maxLim, minDist, trials )


    print('Conducting Voronoi tesselation...', end='')
    ### conducting Voronoi tesselation
    vor, volumes = utilitiesNumeric.runMirroredVoronoi (node_coords, 3, maxLim)
    ### extracting characteristics of the Vor diagram

    print('done.')

    """
    if SHOW_PLOT:
        fig = plt.figure()
        ax = plt.axes(projection='3d')

        d = np.asarray(node_coords)
        xcoords = d[:,0]
        ycoords = d[:,1]
        zcoords = d[:,2]

        ax.scatter3D(xcoords,ycoords,zcoords)
        ax.set_xlim3d(-maxLim[0]*1,2*maxLim[0])
        ax.set_ylim3d(-maxLim[1]*1,2*maxLim[1])
        ax.set_zlim3d(-maxLim[2]*1,2*maxLim[2])
        plt.show()
    """

    ########################################################################
    functions = []
    ### creating functions
    #### Defining functions
    #0 constant zero
    fn = utilitiesNumeric.constantFunc(0)
    functions.append (fn)

    #1 loading function. Const

    func = []
    func.append( np.array([0,0]) )
    func.append( np.array([10, -1e-2]) )
    fn1= utilitiesNumeric.generalFunc(func)
    #fn1 = utilitiesNumeric.sawToothConstFunc(value = -1e-1, period = 10, sym = 1)
    functions.append (fn1)

    #transport function, leftFace, constant
    fn2 = utilitiesNumeric.constantFunc(20)
    functions.append (fn2)

    #transport function, rightFace, bilinear
    func3 = []
    func3.append( np.array([0,0]) )
    func3.append( np.array([50, 500]) )
    fn3 = utilitiesNumeric.generalFunc(func3)
    functions.append (fn3)

    ########################################################################
    ### indirect setting of transportBCs by spatial selection of vertices
    transportBC_merged = []
    transportIC_merged = []
    ### selecting vertices on the left surface
    leftFaceBC = np.array([2,-1])
    boundA = np.array(  [-1e-8 , 0] )
    boundB = np.array(  [ 1e-8 , maxLim[1]]  )
    leftFace = utilitiesGeom.returnSelectedPts(boundA, boundB, vor.vertices)
    #print(leftFace)
    for i in range (len(leftFace)):
        trsBC = utilitiesMech.transportBC(leftFace[i], leftFaceBC)
        transportBC_merged.append(trsBC)

    ### selecting vertices on the right surface
    rightFaceBC = np.array([3,-1])
    boundA = np.array(  [maxLim[0] - 1e-8, 0] )
    boundB = np.array(  [maxLim[0] + 1e-8 , maxLim[1]]  )
    rightFace = utilitiesGeom.returnSelectedPts(boundA, boundB, vor.vertices)
    #print(rightFace)
    for i in range (len(rightFace)):
        trsBC = utilitiesMech.transportBC(rightFace[i], rightFaceBC)
        transportBC_merged.append(trsBC)

    return node_coords, mechBC_merged, mechIC_merged, transportBC_merged, transportIC_merged, vor, volumes, functions








def create3dCantileverUniPressFree(maxLim, minDist, trials ):
    ### sampling of nodes
    ### direct setting of mechanicalBCs
    node_coords, mechBC_merged, mechIC_merged  = assemble3dCantileverUniPressFree(maxLim, minDist, trials )


    print('Conducting Voronoi tesselation...', end='')
    ### conducting Voronoi tesselation
    vor, volumes = utilitiesNumeric.runMirroredVoronoi (node_coords, 3, maxLim)
    ### extracting characteristics of the Vor diagram

    print('done.')


    ########################################################################
    functions = []
    ### creating functions
    #### Defining functions
    #0 constant zero
    fn = utilitiesNumeric.constantFunc(0)
    functions.append (fn)

    #1 loading function. Const

    func = []
    func.append( np.array([0,0]) )
    func.append( np.array([10, -1e-2]) )
    fn1= utilitiesNumeric.generalFunc(func)
    #fn1 = utilitiesNumeric.sawToothConstFunc(value = -1e-1, period = 10, sym = 1)
    functions.append (fn1)

    #transport function, leftFace, constant
    fn2 = utilitiesNumeric.constantFunc(20)
    functions.append (fn2)

    #transport function, rightFace, bilinear
    func3 = []
    func3.append( np.array([0,0]) )
    func3.append( np.array([10, 500]) )
    fn3 = utilitiesNumeric.generalFunc(func3)
    functions.append (fn3)

    ########################################################################
    ### indirect setting of transportBCs by spatial selection of vertices
    transportBC_merged = []
    transportIC_merged = []
    ### selecting vertices on the left surface
    leftFaceBC = np.array([2,-1])
    boundA = np.array(  [-1e-8 , -1e-8, -1e-8 ] )
    boundB = np.array(  [ 1e-8 , maxLim[1]+1e8,  maxLim[2]+1e8]  )
    leftFace = utilitiesGeom.returnSelectedPts(boundA, boundB, vor.vertices)
    #print(leftFace)
    for i in range (len(leftFace)):
        trsBC = utilitiesMech.transportBC(leftFace[i], leftFaceBC)
        transportBC_merged.append(trsBC)

    ### selecting vertices on the right surface
    rightFaceBC = np.array([3,-1])
    boundA = np.array(  [maxLim[0] - 1e-8, -1e8, -1e8] )
    boundB = np.array(  [maxLim[0] + 1e-8 , maxLim[1]+1e8, maxLim[2]+1e8 ]  )
    rightFace = utilitiesGeom.returnSelectedPts(boundA, boundB, vor.vertices)
    #print(rightFace)
    for i in range (len(rightFace)):
        trsBC = utilitiesMech.transportBC(rightFace[i], rightFaceBC)
        transportBC_merged.append(trsBC)


    return node_coords, mechBC_merged, mechIC_merged, transportBC_merged, transportIC_merged, vor, volumes, functions









def create3dCantileverUniPressConfined(maxLim, minDist, trials ):
    ### sampling of nodes
    ### direct setting of mechanicalBCs
    node_coords, mechBC_merged, mechIC_merged  = assemble3dCantileverUniPressConfined(maxLim, minDist, trials )


    print('Conducting Voronoi tesselation...', end='')
    ### conducting Voronoi tesselation
    vor, volumes = utilitiesNumeric.runMirroredVoronoi (node_coords, 3, maxLim)
    ### extracting characteristics of the Vor diagram

    print('done.')


    ########################################################################
    functions = []
    ### creating functions
    #### Defining functions
    #0 constant zero
    fn = utilitiesNumeric.constantFunc(0)
    functions.append (fn)

    #1 loading function. Const
    func = []
    func.append( np.array([0,0]) )
    func.append( np.array([10, -1e-2]) )
    fn1= utilitiesNumeric.generalFunc(func)
    #fn1 = utilitiesNumeric.sawToothConstFunc(value = -1e-1, period = 10, sym = 1)
    functions.append (fn1)

    #transport function, leftFace, constant
    fn2 = utilitiesNumeric.constantFunc(20)
    functions.append (fn2)

    #transport function, rightFace, bilinear
    func3 = []
    func3.append( np.array([0,0]) )
    func3.append( np.array([10, 500]) )
    fn3 = utilitiesNumeric.generalFunc(func3)
    functions.append (fn3)

    ########################################################################
    ### indirect setting of transportBCs by spatial selection of vertices
    transportBC_merged = []
    transportIC_merged = []
    ### selecting vertices on the left surface
    leftFaceBC = np.array([2,-1])
    boundA = np.array(  [-1e-8 , -1e-8, -1e-8 ] )
    boundB = np.array(  [ 1e-8 , maxLim[1]+1e8,  maxLim[2]+1e8]  )
    leftFace = utilitiesGeom.returnSelectedPts(boundA, boundB, vor.vertices)
    #print(leftFace)
    for i in range (len(leftFace)):
        trsBC = utilitiesMech.transportBC(leftFace[i], leftFaceBC)
        transportBC_merged.append(trsBC)

    ### selecting vertices on the right surface
    rightFaceBC = np.array([3,-1])
    boundA = np.array(  [maxLim[0] - 1e-8, -1e8, -1e8] )
    boundB = np.array(  [maxLim[0] + 1e-8 , maxLim[1]+1e8, maxLim[2]+1e8 ]  )
    rightFace = utilitiesGeom.returnSelectedPts(boundA, boundB, vor.vertices)
    #print(rightFace)
    for i in range (len(rightFace)):
        trsBC = utilitiesMech.transportBC(rightFace[i], rightFaceBC)
        transportBC_merged.append(trsBC)

    return node_coords, mechBC_merged, mechIC_merged, transportBC_merged, transportIC_merged, vor, volumes, functions





def create3dcylinderUniPressConfined(center, radius, height, minDist, trials, directionDim ):
    ### sampling of nodes
    ### direct setting of mechanicalBCs
    node_coords, mechBC_merged, mechIC_merged  = assemble3dcylinderUniPressConfined(center, radius, height, minDist, trials, directionDim )

    #print(*node_coords, sep='\n')

    #node_coords = np.asarray(node_coords)
    """
    if SHOW_PLOT:
        fig = plt.figure()
        ax = Axes3D(fig)
        ax.scatter(node_coords[:,0], node_coords[:,1], node_coords[:,2])
        plt.show()
    """
    print('Conducting Voronoi tesselation...', end='')
    ### conducting Voronoi tesselation
    vor, volumes = utilitiesNumeric.runCylinderMirroredVoronoi (node_coords, center, radius, height, directionDim)
    ### extracting characteristics of the Vor diagram
    print('done.')

    ########################################################################
    functions = []
    ### creating functions
    #### Defining functions
    #0 constant zero
    fn = utilitiesNumeric.constantFunc(0)
    functions.append (fn)

    #1 loading function. Const
    func = []
    func.append( np.array([0,0]) )
    func.append( np.array([10, -1e-2]) )
    fn1= utilitiesNumeric.generalFunc(func)
    #fn1 = utilitiesNumeric.sawToothConstFunc(value = -1e-1, period = 10, sym = 1)
    functions.append (fn1)

    #transport function, leftFace, constant
    fn2 = utilitiesNumeric.constantFunc(20)
    functions.append (fn2)

    #transport function, rightFace, bilinear
    func3 = []
    func3.append( np.array([0,0]) )
    func3.append( np.array([10, 500]) )
    fn3 = utilitiesNumeric.generalFunc(func3)
    functions.append (fn3)

    ########################################################################
    ### indirect setting of transportBCs by spatial selection of vertices
    transportBC_merged = []
    transportIC_merged = []
    ### selecting vertices on the left surface
    leftFaceBC = np.array([2,-1])
    boundA = np.array(  [-1e-8 , center[1]-radius, center[2]-radius] )
    boundB = np.array(  [ 1e-8 , center[1]+radius,  center[2]+radius]  )
    leftFace = utilitiesGeom.returnSelectedPts(boundA, boundB, vor.vertices)
    #print(leftFace)
    for i in range (len(leftFace)):
        trsBC = utilitiesMech.transportBC(leftFace[i], leftFaceBC)
        transportBC_merged.append(trsBC)

    ### selecting vertices on the right surface
    rightFaceBC = np.array([3,-1])
    boundA = np.array(  [height-1e-8 , center[1]-radius, center[2]-radius] )
    boundB = np.array(  [height+1e-8 , center[1]+radius,  center[2]+radius]  )
    rightFace = utilitiesGeom.returnSelectedPts(boundA, boundB, vor.vertices)
    #print(rightFace)
    for i in range (len(rightFace)):
        trsBC = utilitiesMech.transportBC(rightFace[i], rightFaceBC)
        transportBC_merged.append(trsBC)
        print(i)


    return node_coords, mechBC_merged, mechIC_merged, transportBC_merged, transportIC_merged, vor, volumes, functions




def create3dcylinderUniPressFree(center, radius, height, minDist, trials, directionDim ):
    ### sampling of nodes
    ### direct setting of mechanicalBCs
    node_coords, mechBC_merged, mechIC_merged  = assemble3dcylinderUniPressFree(center, radius, height, minDist, trials, directionDim )

    #print(*node_coords, sep='\n')

    #node_coords = np.asarray(node_coords)
    """
    if SHOW_PLOT:
        fig = plt.figure()
        ax = Axes3D(fig)
        ax.scatter(node_coords[:,0], node_coords[:,1], node_coords[:,2])
        plt.show()
    """
    print('Conducting Voronoi tesselation...', end='')
    ### conducting Voronoi tesselation
    vor, volumes = utilitiesNumeric.runCylinderMirroredVoronoi (node_coords, center, radius, height, directionDim)
    ### extracting characteristics of the Vor diagram
    print('done.')

    ########################################################################
    functions = []
    ### creating functions
    #### Defining functions
    #0 constant zero
    fn = utilitiesNumeric.constantFunc(0)
    functions.append (fn)

    #1 loading function. Const
    func = []
    func.append( np.array([0,0]) )
    func.append( np.array([10, -1e-2]) )
    fn1= utilitiesNumeric.generalFunc(func)
    #fn1 = utilitiesNumeric.sawToothConstFunc(value = -1e-1, period = 10, sym = 1)
    functions.append (fn1)

    #transport function, leftFace, constant
    fn2 = utilitiesNumeric.constantFunc(20)
    functions.append (fn2)

    #transport function, rightFace, bilinear
    func3 = []
    func3.append( np.array([0,0]) )
    func3.append( np.array([10, 500]) )
    fn3 = utilitiesNumeric.generalFunc(func3)
    functions.append (fn3)


    ########################################################################
    ### indirect setting of transportBCs by spatial selection of vertices
    transportBC_merged = []
    transportIC_merged = []
    ### selecting vertices on the left surface
    leftFaceBC = np.array([2,-1])
    boundA = np.array(  [-1e-8 , center[1]-radius, center[2]-radius] )
    boundB = np.array(  [ 1e-8 , center[1]+radius,  center[2]+radius]  )
    leftFace = utilitiesGeom.returnSelectedPts(boundA, boundB, vor.vertices)
    #print(leftFace)
    for i in range (len(leftFace)):
        trsBC = utilitiesMech.transportBC(leftFace[i], leftFaceBC)
        transportBC_merged.append(trsBC)

    ### selecting vertices on the right surface
    rightFaceBC = np.array([3,-1])
    boundA = np.array(  [height-1e-8 , center[1]-radius, center[2]-radius] )
    boundB = np.array(  [height+1e-8 , center[1]+radius,  center[2]+radius]  )
    rightFace = utilitiesGeom.returnSelectedPts(boundA, boundB, vor.vertices)
    #print(rightFace)
    for i in range (len(rightFace)):
        trsBC = utilitiesMech.transportBC(rightFace[i], rightFaceBC)
        transportBC_merged.append(trsBC)
        print(i)


    return node_coords, mechBC_merged, mechIC_merged, transportBC_merged, transportIC_merged, vor, volumes, functions





def create3dcylinderTorsionFree(center, radius, height, minDist, trials, directionDim, powerTes=False ):

    ########################################################################
    functions = []
    ### creating functions
    #### Defining functions
    #0 constant zero
    fn = utilitiesNumeric.constantFunc(0)
    functions.append (fn)



    ### sampling of nodes
    ### direct setting of mechanicalBCs
    node_coords, mechBC_merged, mechIC_merged  = assemble3dcylinderTorsionFree(center, radius, height, minDist, trials, directionDim, functions )

    #node_coords, mechBC_merged, mechInitC_merged, govNodes, govNodesMechBC, rigidPlates, radii  = assembleCoupledBrazilianDisc(center, cylinderRad, cylinderHeight, minDist, trials, 0, powerTes = powerTes )


    #print(*node_coords, sep='\n')

    #node_coords = np.asarray(node_coords)

    """
    if SHOW_PLOT:
        fig = plt.figure()
        ax = Axes3D(fig)
        ax.scatter(node_coords[:,0], node_coords[:,1], node_coords[:,2])
        plt.show()
    """

    print('Conducting Voronoi tesselation...', end='')
    ### conducting Voronoi tesselation
    vor, volumes = utilitiesNumeric.runCylinderMirroredVoronoi (node_coords, center, radius, height, directionDim)
    ### extracting characteristics of the Vor diagram
    print('done.')

    ########################################################################
    ### indirect setting of transportBCs by spatial selection of vertices
    transportBC_merged = []
    transportIC_merged = []
    ### selecting vertices on the left surface
    leftFaceBC = np.array([2,-1])
    boundA = np.array(  [-1e-8 , center[1]-radius,  center[2]-radius] )
    boundB = np.array(  [ 1e-8 , center[1]+radius,  center[2]+radius]  )
    leftFace = utilitiesGeom.returnSelectedPts(boundA, boundB, vor.vertices)
    #print(leftFace)
    for i in range (len(leftFace)):
        trsBC = utilitiesMech.transportBC(leftFace[i], leftFaceBC)
        transportBC_merged.append(trsBC)

    ### selecting vertices on the right surface
    rightFaceBC = np.array([3,-1])
    boundA = np.array(  [height-1e-8 , center[1]-radius, center[2]-radius] )
    boundB = np.array(  [height+1e-8 , center[1]+radius,  center[2]+radius]  )
    rightFace = utilitiesGeom.returnSelectedPts(boundA, boundB, vor.vertices)
    #print(rightFace)
    for i in range (len(rightFace)):
        trsBC = utilitiesMech.transportBC(rightFace[i], rightFaceBC)
        transportBC_merged.append(trsBC)
        print(i)


    return node_coords, mechBC_merged, mechIC_merged, transportBC_merged, transportIC_merged, vor, volumes, functions


def create3dcylinderTorsionPressFree(center, radius, height, minDist, trials, directionDim, activeTransport, powerTes = False):

    ########################################################################
    functions = []
    ### creating functions
    #### Defining functions
    #0 constant zero
    fn = utilitiesNumeric.constantFunc(0)
    functions.append (fn)

    #1 loading function, pressure X
    func1 = []
    func1.append( np.array([0,0]) )
    func1.append( np.array([0.5, -1e-3]) )
    func1.append( np.array([1, -1e-3]) )
    fn1 = utilitiesNumeric.generalFunc(func1)
    functions.append (fn1)

    #1 loading function, rotation X
    func2 = []
    func2.append( np.array([0,0]) )
    func2.append( np.array([0.5, 0]) )
    func2.append( np.array([1, 1e-3]) )
    fn2 = utilitiesNumeric.generalFunc(func2)
    functions.append (fn2)


    ### sampling of nodes
    ### direct setting of mechanicalBCs
    node_coords, mechBC_merged, govNodes, govNodesMechBC, rigidPlates, radii  = assemble3dcylinderTorsionPressFree(center, radius, height, minDist, trials, directionDim, functions, powerTes=powerTes )

    """
    node_coords = np.asarray(node_coords)
    if SHOW_PLOT:
        fig = plt.figure()
        ax = Axes3D(fig)
        ax.scatter(node_coords[:,0], node_coords[:,1], node_coords[:,2])
        plt.show()
    """

    print('Conducting Voronoi tesselation...', end='')
    ### conducting Voronoi tesselation
    vor, volumes = utilitiesNumeric.runCylinderMirroredVoronoi (node_coords, center, radius, height, directionDim)
    ### extracting characteristics of the Vor diagram
    print('done.')



    ########################################################################
    ### indirect setting of transportBCs by spatial selection of vertices
    transportBC_merged = []

    transportIC_merged = []

    govNodesTrspt = []
    govNodesTrsptBC = []
    rigidPlatesTrspt = []

    maxLim = np.array([height, 2*radius, 2*radius])

    if activeTransport:
        #"""
        ### selecting vertices on the bottom surface
        boundA = np.array([-minDist/10, -100, -100])
        boundB = np.array([minDist/10, 100, 100])
        faces1 = utilitiesGeom.returnSelectedPts(boundA, boundB, vor.vertices)
        vert = vor.vertices[faces1,:]



        fn1 = utilitiesNumeric.constantFunc(100)
        functions.append (fn1)

        ### selecting vertices on the top surface
        boundA = np.array([maxLim[0]-minDist/10, -100, -100])
        boundB = np.array([maxLim[0]+minDist/10, 100, 100])
        faces2 = utilitiesGeom.returnSelectedPts(boundA, boundB, vor.vertices)
        vert = vor.vertices[faces2,:]

        fn1 = utilitiesNumeric.constantFunc(0)
        functions.append (fn1)



        print(faces1)
        print(faces2)

        trsptLeftRigidPlate = utilitiesMech.RigidPlate(-len(govNodesTrspt)-1, 3, None, directIdcs = True)
        trsptLeftRigidPlate.setDirectNodes(faces1)
        trsptLeftRigidPlateMechBC = np.array([3,-1])
        rigidPlatesTrspt.append(trsptLeftRigidPlate)
        govNodesTrspt.append(np.array([ 0, 0, 0]))
        govNodesTrsptBC.append(utilitiesMech.transportBC(govNodesTrspt[-1], trsptLeftRigidPlateMechBC))

        trsptRightRigidPlate = utilitiesMech.RigidPlate(-len(govNodesTrspt)-1, 3, None, directIdcs = True)
        trsptRightRigidPlate.setDirectNodes(faces2)
        trsptRightRigidPlateMechBC = np.array([4,-1])
        rigidPlatesTrspt.append(trsptRightRigidPlate)
        govNodesTrspt.append(np.array([ height, 0, 0]))
        govNodesTrsptBC.append(utilitiesMech.transportBC(govNodesTrspt[-1], trsptRightRigidPlateMechBC))


    return node_coords, mechBC_merged,  vor, volumes, functions,  govNodes, govNodesMechBC, rigidPlates, transportBC_merged, govNodesTrspt, rigidPlatesTrspt, govNodesTrsptBC, radii

    #return node_coords, mechBC_merged, transportBC_merged, govNodes, govNodesMechBC, rigidPlates, vor, volumes, functions, rigidPlatesTrspt, govNodesTrspt, govNodesTrsptBC



def create3dRWTHShearCylinder(center, radius, height, minDist, trials, notchRadLeft, notchRadRight, notchWidthLeft, notchWidthRight, notchDepth, quarter=False):
    directionDim=0
    ########################################################################
    functions = []
    ### creating functions
    #### Defining functions
    #0 constant zero
    fn = utilitiesNumeric.constantFunc(0)
    functions.append (fn)

    #1 loading function, pressure X
    func1 = []
    func1.append( np.array([0,0]) )
    func1.append( np.array([1, -1e-3]) )
    fn1 = utilitiesNumeric.generalFunc(func1)
    functions.append (fn1)




    ### sampling of nodes
    ### direct setting of mechanicalBCs
    node_coords, mechBC_merged, govNodes, govNodesMechBC, rigidPlates, notches  = assemble3dRWTHShearCylinder(center, radius, height, minDist, trials, directionDim, functions, notchRadLeft, notchRadRight, notchWidthLeft, notchWidthRight,notchDepth, quarter = quarter )

    """
    node_coords = np.asarray(node_coords)
    if SHOW_PLOT:
        fig = plt.figure()
        ax = Axes3D(fig)
        ax.scatter(node_coords[:,0], node_coords[:,1], node_coords[:,2])
        plt.show()
    """

    print('Conducting Voronoi tesselation...', end='')
    ### conducting Voronoi tesselation
    vor, volumes = utilitiesNumeric.runCylinderMirroredVoronoi (node_coords, center, radius, height, directionDim, quarter = quarter)
    ### extracting characteristics of the Vor diagram
    print('done.')



    return node_coords, mechBC_merged,  vor, volumes, functions,govNodes, govNodesMechBC, rigidPlates, notches



def create2dRWTHShearCylinder(radius, height, minDist, trials, innerRadTop, innerRadBottom, notchWidth, notchDepth):
    directionDim=0
    ########################################################################
    functions = []
    ### creating functions
    #### Defining functions
    #0 constant zero
    fn = utilitiesNumeric.constantFunc(0)
    functions.append (fn)

    #1 loading function, pressure X
    func1 = []
    func1.append( np.array([0,0]) )
    func1.append( np.array([1, -1e-3]) )
    fn1 = utilitiesNumeric.generalFunc(func1)
    functions.append (fn1)


    maxLim = np.array([radius*2, height])

    ### sampling of nodes
    ### direct setting of mechanicalBCs
    node_coords, mechBC_merged, mechInitC_merged, notches, govNodes, govNodesMechBC, rigidPlates = assemble2dRWTHShearCylinder(maxLim, minDist, trials, innerRadTop, innerRadBottom, notchWidth, notchDepth)



    """
    node_coords = np.asarray(node_coords)
    fig, ax = plt.subplots()
    ax.scatter(node_coords[:,0], node_coords[:,1])
    if SHOW_PLOT:
        plt.show()
    #"""

    print('Conducting Voronoi tesselation...', end='')
    ### conducting Voronoi tesselation
    vor, regions, vertices, polygons, areas, centroids, points = utilitiesNumeric.runMirroredVoronoi (node_coords, 2, maxLim)
    ### extracting characteristics of the Vor diagram
    print('done.')



    return node_coords, mechBC_merged,  vor, areas, functions,govNodes, govNodesMechBC, rigidPlates, notches





def create3dtubeTorsionFree(center, radius, height, thickness, minDist, trials, directionDim, rotationAngle = 0.001 ):

    ########################################################################
    functions = []
    ### creating functions
    #### Defining functions
    #0 constant zero
    fn = utilitiesNumeric.constantFunc(0)
    functions.append (fn)

    directionDim = int(directionDim)
    ### sampling of nodes
    ### direct setting of mechanicalBCs
    node_coords, mechBC_merged, mechIC_merged = assemble3dtubeTorsionFree(center, radius, height, thickness, minDist, trials, directionDim, functions, rotationAngle)
    #node_coords, mechBC_merged, mechIC_merged = assemble3dslimTubeTorsionFree(center, radius, height, thickness, minDist, trials, directionDim, functions )

    #print(*node_coords, sep='\n')


    node_coords = np.asarray(node_coords)
    """
    if SHOW_PLOT:
        fig = plt.figure()
        ax = Axes3D(fig)
        ax.scatter(node_coords[:,0], node_coords[:,1], node_coords[:,2])
        plt.show()
    """
    print('Conducting Voronoi tesselation...', end='')
    ### conducting Voronoi tesselation
    vor, volumes = utilitiesNumeric.runTubeMirroredVoronoi (node_coords, center, radius, height, thickness, directionDim)
    ### extracting characteristics of the Vor diagram
    print('done.')

    ########################################################################
    ### indirect setting of transportBCs by spatial selection of vertices
    transportBC_merged = []
    transportIC_merged = []
    ### selecting vertices on the left surface
    leftFaceBC = np.array([2,-1])
    boundA = np.array(  [-1e-8 , center[1]-radius,  center[2]-radius] )
    boundB = np.array(  [ 1e-8 , center[1]+radius,  center[2]+radius]  )
    leftFace = utilitiesGeom.returnSelectedPts(boundA, boundB, vor.vertices)
    #print(leftFace)
    for i in range (len(leftFace)):
        trsBC = utilitiesMech.transportBC(leftFace[i], leftFaceBC)
        transportBC_merged.append(trsBC)

    ### selecting vertices on the right surface
    rightFaceBC = np.array([3,-1])
    boundA = np.array(  [height-1e-8 , center[1]-radius, center[2]-radius] )
    boundB = np.array(  [height+1e-8 , center[1]+radius,  center[2]+radius]  )
    rightFace = utilitiesGeom.returnSelectedPts(boundA, boundB, vor.vertices)
    #print(rightFace)
    for i in range (len(rightFace)):
        trsBC = utilitiesMech.transportBC(rightFace[i], rightFaceBC)
        transportBC_merged.append(trsBC)
        print(i)


    return node_coords, mechBC_merged, mechIC_merged, transportBC_merged, transportIC_merged, vor, volumes, functions




def create3dBiparvaTubeTransport( radius, height, thickness, minDist, trials, maxLim):
    center = np.zeros((3))
    ########################################################################
    functions = []
    fn = utilitiesNumeric.constantFunc(0)
    functions.append (fn)


    func1 = []
    func1.append( np.array([0,0]) )
    func1.append( np.array([1, -1e-3]) )
    fn1 = utilitiesNumeric.generalFunc(func1)
    functions.append (fn1)

    fn1 = utilitiesNumeric.constantFunc(100)
    functions.append (fn1)

    ### sampling of nodes
    node_coords, mechBC_merged,  govNodes, govNodesMechBC, rigidPlates = assemble3dBiparvaTubeTransport(center, radius, height, thickness, minDist, trials)
    node_coords = np.asarray(node_coords)

    print (len(node_coords))

    fig = plt.figure()
    ax = Axes3D(fig)
    ax.scatter(node_coords[:,0], node_coords[:,1], node_coords[:,2])
    plt.show()

    print('Conducting Voronoi tesselation...', end='')
    directionDim = 0
    vor, volumes = utilitiesNumeric.runTubeMirroredVoronoi (node_coords, center, radius, height, thickness, directionDim)
    print('done.')

    ########################################################################
    ### indirect setting of transportBCs by spatial selection of vertices
    transportBC_merged = []
    transportIC_merged = []

    govNodesTrspt = []
    govNodesTrspt.append(np.array([ 0, 0, 0]))


    govNodesTrsptBC = []
    rigidPlatesTrspt = []

    modelVertices = utilitiesGeom.returnSelectedPtsRadial (radius-thickness-1e-3 , radius+1e-3 , vor.vertices)
    ### selecting vertices on the outer surface

    outerFaceBC = np.array([2,-1])
    outerFace = utilitiesGeom.returnSelectedPtsRadial (radius-minDist/2 , radius+minDist/2 , vor.vertices, xmin = 1e-5, xmax = height - 1e-5)
    trsptOuterRigidPlate = utilitiesMech.RigidPlate(-len(govNodesTrspt)-1, 3, None, directIdcs = True)
    trsptOuterRigidPlate.setDirectNodes(outerFace)
    rigidPlatesTrspt.append(trsptOuterRigidPlate)
    govNodesTrsptBC.append(utilitiesMech.transportBC(govNodesTrspt[-1], outerFaceBC))

    #for i in range (len(outerFace)):
    #   trsBC = utilitiesMech.transportBC(outerFace[i], outerFaceBC)
    #   transportBC_merged.append(trsBC)

    if SHOW_PLOT:
        fig = plt.figure()
        ax = Axes3D(fig)
        ax.scatter(vor.vertices[modelVertices,0], vor.vertices[modelVertices,1], vor.vertices[modelVertices,2])
        ax.scatter(vor.vertices[outerFace,0], vor.vertices[outerFace,1], vor.vertices[outerFace,2])
        plt.show()

    govNodesTrspt.append(np.array([ -1, -1, -1]))
    innerFaceBC = np.array([0,-1])
    innerFace = utilitiesGeom.returnSelectedPtsRadial ((radius-thickness)-minDist/2 , (radius-thickness)+minDist/2, vor.vertices,xmin = 1e-5, xmax = height - 1e-5)
    trsptInnerRigidPlate = utilitiesMech.RigidPlate(-len(govNodesTrspt)-1, 3, None, directIdcs = True)
    trsptInnerRigidPlate.setDirectNodes(innerFace)
    rigidPlatesTrspt.append(trsptInnerRigidPlate)
    govNodesTrsptBC.append(utilitiesMech.transportBC(govNodesTrspt[-1], innerFaceBC))


    govNodesTrspt.append(np.array([ 1, -1, -1]))
    topFaceBC = np.array([-1,0])
    boundA = np.array(  [-1e-5 , -radius*10, -radius*10] )
    boundB = np.array(  [ 1e-5 , radius*10, radius*10] )
    topFace = utilitiesGeom.returnSelectedPts(boundA, boundB, vor.vertices)
    print(topFace)
    trsptTopRigidPlate = utilitiesMech.RigidPlate(-len(govNodesTrspt)-1, 3, None, directIdcs = True)
    trsptTopRigidPlate.setDirectNodes(topFace)
    rigidPlatesTrspt.append(trsptTopRigidPlate)
    govNodesTrsptBC.append(utilitiesMech.transportBC(govNodesTrspt[-1], topFaceBC))

    govNodesTrspt.append(np.array([ 2, -1, -1]))
    boundA = np.array(  [height-1e-5 , -radius*10, -radius*10] )
    boundB = np.array(  [height+1e-5  , radius*10, radius*10] )
    botFace = utilitiesGeom.returnSelectedPts(boundA, boundB, vor.vertices)
    print(botFace)
    trsptBotRigidPlate = utilitiesMech.RigidPlate(-len(govNodesTrspt)-1, 3, None, directIdcs = True)
    trsptBotRigidPlate.setDirectNodes(botFace)
    rigidPlatesTrspt.append(trsptBotRigidPlate)
    govNodesTrsptBC.append(utilitiesMech.transportBC(govNodesTrspt[-1], topFaceBC))



    #for i in range (len(innerFace)):
    #    trsBC = utilitiesMech.transportBC(innerFace[i], innerFaceBC)
    #    transportBC_merged.append(trsBC)

    if SHOW_PLOT:
        fig = plt.figure()
        ax = Axes3D(fig)
        ax.scatter(vor.vertices[modelVertices,0], vor.vertices[modelVertices,1], vor.vertices[modelVertices,2])
        ax.scatter(vor.vertices[innerFace,0], vor.vertices[innerFace,1], vor.vertices[innerFace,2])
        plt.show()



    radii = np.zeros((len(node_coords))) + minDist
    return node_coords, mechBC_merged, govNodes, govNodesMechBC, rigidPlates, transportBC_merged, vor, volumes, functions, radii, rigidPlatesTrspt, govNodesTrspt, govNodesTrsptBC





def create3dTubeSplit( radius, height, thickness, minDist, trials, maxLim,notchH=0 ,notchWidth=0):
    center = np.zeros((3))
    ########################################################################
    functions = []
    fn = utilitiesNumeric.constantFunc(0)
    functions.append (fn)


    func1 = []
    func1.append( np.array([0,0]) )
    func1.append( np.array([1, 1e-3]) )
    fn1 = utilitiesNumeric.generalFunc(func1)
    functions.append (fn1)

    fn1 = utilitiesNumeric.constantFunc(100)
    functions.append (fn1)

    ### sampling of nodes
    node_coords, mechBC_merged,  govNodes, govNodesMechBC, rigidPlates,notches = assemble3dTubeSplit(center, radius, height, thickness, minDist, trials,notchH=notchH, notchWidth=notchWidth)
    node_coords = np.asarray(node_coords)

    print (len(node_coords))

    fig = plt.figure()
    ax = Axes3D(fig)
    ax.scatter(node_coords[:,0], node_coords[:,1], node_coords[:,2])
    plt.show()

    print('Conducting Voronoi tesselation...', end='')
    directionDim = 0
    vor, volumes = utilitiesNumeric.runTubeMirroredVoronoi (node_coords, center, radius, height, thickness, directionDim)
    print('done.')

    ########################################################################
    ### indirect setting of transportBCs by spatial selection of vertices
    transportBC_merged = []
    transportIC_merged = []

    govNodesTrspt = []


    govNodesTrsptBC = []
    rigidPlatesTrspt = []




    radii = np.zeros((len(node_coords))) + minDist
    return node_coords, mechBC_merged, govNodes, govNodesMechBC, rigidPlates, transportBC_merged, vor, volumes, functions, radii, rigidPlatesTrspt, govNodesTrspt, govNodesTrsptBC,notches







def create3dTubeInnerPressure( radius, height, thickness, minDist, trials, maxLim, activeTransport, powerTes):

    center = np.zeros((3))
    ########################################################################
    functions = []
    fn = utilitiesNumeric.constantFunc(0)
    functions.append (fn)


    func1 = []
    func1.append( np.array([0,0]) )
    func1.append( np.array([1, -1e-3]) )
    fn1 = utilitiesNumeric.generalFunc(func1)
    functions.append (fn1)

    #2 mech loading function
    func2 = []
    func2.append( np.array([0,0]) )
    func2.append( np.array([1, 1e-2]) )
    fn2 = utilitiesNumeric.generalFunc(func2)
    functions.append (fn2)

    ### sampling of nodes
    node_coords, mechBC_merged,  govNodes, govNodesMechBC, rigidPlates, radii = assemble3dTubeInnerPressure(center, radius, height, thickness, minDist, trials, powerTes)
    node_coords = np.asarray(node_coords)

    print (len(node_coords))

    #fig = plt.figure()
    #ax = Axes3D(fig)
    #ax.scatter(node_coords[:,0], node_coords[:,1], node_coords[:,2])
    #plt.show()

    print('Conducting Voronoi tesselation...', end='')
    directionDim = 0
    if powerTes:
        vor, volumes = utilitiesNumeric.runTubeMirroredPower (node_coords, center, radius, height, thickness, directionDim, radii)
    else:
        vor, volumes = utilitiesNumeric.runTubeMirroredVoronoi (node_coords, center, radius, height, thickness, directionDim)
        radii = []
    print('done.')

    ########################################################################
    ### indirect setting of transportBCs by spatial selection of vertices
    transportBC_merged = []
    transportIC_merged = []

    govNodesTrspt = []
    interfaceVertexIndices = []


    govNodesTrsptBC = []
    rigidPlatesTrspt = []

    if activeTransport:
        govNodesTrspt.append(np.array([ 0, 0, 0]))

        modelVertices = utilitiesGeom.returnSelectedPtsRadial (radius-thickness-1e-3 , radius+1e-3 , vor.vertices)
        ### selecting vertices on the outer surface

        outerFaceBC = np.array([0,-1])
        outerFace = utilitiesGeom.returnSelectedPtsRadial (radius-minDist/2 , radius+minDist/2 , vor.vertices, xmin = -1e-5, xmax = height + 1e-5) # JE: I need also edges
        trsptOuterRigidPlate = utilitiesMech.RigidPlate(-len(govNodesTrspt)-1, 3, None, directIdcs = True)
        trsptOuterRigidPlate.setDirectNodes(outerFace)
        rigidPlatesTrspt.append(trsptOuterRigidPlate)
        govNodesTrsptBC.append(utilitiesMech.transportBC(govNodesTrspt[-1], outerFaceBC))

        #for i in range (len(outerFace)):
        #   trsBC = utilitiesMech.transportBC(outerFace[i], outerFaceBC)
        #   transportBC_merged.append(trsBC)

        if SHOW_PLOT:
            fig = plt.figure()
            ax = Axes3D(fig)
            ax.scatter(vor.vertices[modelVertices,0], vor.vertices[modelVertices,1], vor.vertices[modelVertices,2])
            ax.scatter(vor.vertices[outerFace,0], vor.vertices[outerFace,1], vor.vertices[outerFace,2])
            plt.show()

        govNodesTrspt.append(np.array([ -1, -1, -1]))
        innerFaceBC = np.array([2,-1])
        innerFace = utilitiesGeom.returnSelectedPtsRadial ((radius-thickness)-minDist/2 , (radius-thickness)+minDist/2, vor.vertices, xmin = -1e-5, xmax = height + 1e-5) # JE: I need also edges
        trsptInnerRigidPlate = utilitiesMech.RigidPlate(-len(govNodesTrspt)-1, 3, None, directIdcs = True)
        trsptInnerRigidPlate.setDirectNodes(innerFace)
        rigidPlatesTrspt.append(trsptInnerRigidPlate)
        govNodesTrsptBC.append(utilitiesMech.transportBC(govNodesTrspt[-1], innerFaceBC))



        interfaceVertexIndices.append(innerFace)
        interfaceVertexIndices = np.asarray(interfaceVertexIndices)

        #TOP and BOTTOM face not needed
        """
        govNodesTrspt.append(np.array([ 1, -1, -1]))
        topFaceBC = np.array([-1,0])
        boundA = np.array(  [-1e-5 , -radius*10, -radius*10] )
        boundB = np.array(  [ 1e-5 , radius*10, radius*10] )
        topFace = utilitiesGeom.returnSelectedPts(boundA, boundB, vor.vertices)
        print(topFace)
        trsptTopRigidPlate = utilitiesMech.RigidPlate(-len(govNodesTrspt)-1, 3, None, directIdcs = True)
        trsptTopRigidPlate.setDirectNodes(topFace)
        rigidPlatesTrspt.append(trsptTopRigidPlate)
        govNodesTrsptBC.append(utilitiesMech.transportBC(govNodesTrspt[-1], topFaceBC))

        govNodesTrspt.append(np.array([ 2, -1, -1]))
        boundA = np.array(  [height-1e-5 , -radius*10, -radius*10] )
        boundB = np.array(  [height+1e-5  , radius*10, radius*10] )
        botFace = utilitiesGeom.returnSelectedPts(boundA, boundB, vor.vertices)
        print(botFace)
        trsptBotRigidPlate = utilitiesMech.RigidPlate(-len(govNodesTrspt)-1, 3, None, directIdcs = True)
        trsptBotRigidPlate.setDirectNodes(botFace)
        rigidPlatesTrspt.append(trsptBotRigidPlate)
        govNodesTrsptBC.append(utilitiesMech.transportBC(govNodesTrspt[-1], topFaceBC))
        """


        #for i in range (len(innerFace)):
        #    trsBC = utilitiesMech.transportBC(innerFace[i], innerFaceBC)
        #    transportBC_merged.append(trsBC)

    """    if SHOW_PLOT:
            fig = plt.figure()
            ax = Axes3D(fig)
            ax.scatter(vor.vertices[modelVertices,0], vor.vertices[modelVertices,1], vor.vertices[modelVertices,2])
            ax.scatter(vor.vertices[innerFace,0], vor.vertices[innerFace,1], vor.vertices[innerFace,2])
            plt.show()
    """



    return node_coords, mechBC_merged, govNodes, govNodesMechBC, rigidPlates, transportBC_merged, vor, volumes, functions, radii, rigidPlatesTrspt, govNodesTrspt, govNodesTrsptBC, interfaceVertexIndices










#
######## METHOD FOR CREATING OF A 2D SUPPORTED CANTILEVER MODEL
def assemble2DCantileverBending (maxLim, minDist, trials, powerTes, nodefile=None):
    dim = 2
    #lists for the model
    node_coords = []
    mechBC_merged = []
    mechIC_merged = []

    #an indent due to mirroring of the data for voronoi tess.
    indent = 1e-8

    ###############generating of nodes, supported line left vertical ###############
    #mech bc
    lineBC = np.array([0,0,0,-1,-1,-1])

    #defining points of the line
    nodeA = np.array([indent, indent])
    nodeB = np.array([indent, maxLim[1]-indent])

    oldLen = len(node_coords)    
    pointGenerators.generateNodesLine2dRand(nodeA, nodeB, minDist, dim, node_coords,  trials, True, False)
    nrOfPoints =  (len(node_coords)) - oldLen
    #print (nrOfPoints)

    #adding mech boundary conditions
    for n in range ( nrOfPoints ):
        mBC = utilitiesMech.mechanicalBC(dim, oldLen + n, lineBC)
        mechBC_merged.append(mBC)
        #print('adding')

    ###############generating a single point top right (a line of zero length) ###############
    lineBC = np.array([-1,-1,-1,-1, 1 ,-1])

    #defining points of the line
    nodeA = np.array([maxLim[0] - indent, maxLim[1] - indent])
    nodeB = np.array([maxLim[0] - indent, indent])

    oldLen = len(node_coords)
    pointGenerators.generateNodesLine2dRand(nodeA, nodeB, minDist, dim, node_coords, trials, False, False)
    #pointGenerators.generateSingleNode (nodeA, dim, node_coords)
    nrOfPoints =  (len(node_coords)) - oldLen
    #print (nrOfPoints)

    
    radii = np.zeros(len(node_coords))
    #adding mech boundary conditions
    for n in range ( nrOfPoints ):
        mBC = utilitiesMech.mechanicalBC(dim, oldLen + n, lineBC)
        mechBC_merged.append(mBC)
        #print('adding')

    if nodefile:
        node_coords_file = np.loadtxt(nodefile)         
        node_coords = np.array(node_coords)
        node_coords = np.vstack((node_coords,node_coords_file[:,:2]))
        if powerTes:           
            radii = np.concatenate((radii,node_coords_file[:,2]), axis=0)
        else: radii = np.zeros(len(node_coords))
    else:

        ##########################################generating of points, homogeneous volume
        rectBC = np.array([-1,-1,-1,-1,-1,-1])
        #rect
        oldLen = len(node_coords)
        if powerTes:
            node_coords, radii = pointGenerators.generateParticlesRect(maxLim, minDist*0.4, minDist, 0.8, 2, trials, np.array(node_coords), np.array(radii), allow_domain_overlap = False, periodic_distance=False)
        else:
            pointGenerators.generateNodesRect(maxLim, minDist, dim, trials, node_coords)
            radii = np.zeros(len(node_coords))



        #
        newLen = len(node_coords)-1

   # print (nrOfPoints)
  #  mBC = utilitiesGeom.mechanicalBC(dim, kvadrBC, oldLen, newLen)
   # mechBC_merged.append(mBC)
    ####################################################################################################

    return node_coords, radii, mechBC_merged, mechIC_merged






#
######## METHOD FOR CREATING OF A 2D SUPPORTED CANTILEVER MODEL
def assemble2DCantileverUniTens (maxLim, minDist, trials):
    dim = 2
    #lists for the model
    node_coords = []
    mechBC_merged = []
    mechIC_merged = []

    #an indent due to mirroring of the data for voronoi tess.
    indent = 1e-8

    nodeA = np.array ( [  indent, indent ] )
    nodeAmechBC = np.array([0, 0 , -1 ,      -1 , -1 , -1])
    pointGenerators.generateSingleNode(nodeA, dim, node_coords)
    mechBC_merged.append( utilitiesMech.mechanicalBC(dim, 0, nodeAmechBC))


    nodeB = np.array ( [  indent, maxLim[1]-indent ] )
    nodeBmechBC = np.array([0, -1 , -1 ,      -1 , -1 , -1])
    pointGenerators.generateSingleNode(nodeB, dim, node_coords)
    mechBC_merged.append(utilitiesMech.mechanicalBC(dim, 1, nodeBmechBC))

    ###############generating of nodes, supported line left vertical ###############
    #mech bc
    lineBC = np.array([0,-1,-1,-1,-1,-1])

    #defining points of the line
    nodeA = np.array([indent, indent])
    nodeB = np.array([indent, maxLim[1]-indent])

    oldLen = len(node_coords)
    pointGenerators.generateNodesLine2dRand(nodeA, nodeB, minDist, dim, node_coords,  trials, False, True)
    nrOfPoints =  (len(node_coords)) - oldLen

    #adding mech boundary conditions
    for n in range ( nrOfPoints ):
        mBC = utilitiesMech.mechanicalBC(dim, oldLen + n, lineBC)
        mechBC_merged.append(mBC)


    ###############generating a nodes on right face, loaded by uni tens in X) ###############
    lineBC = np.array([1,-1,-1,   -1, -1 ,-1])

    #defining points of the line
    nodeA = np.array([maxLim[0] - indent, indent])
    nodeB = np.array([maxLim[0] - indent, maxLim[1]-indent])

    oldLen = len(node_coords)
    pointGenerators.generateNodesLine2dRand(nodeA, nodeB, minDist, dim, node_coords,  trials, True, True)
    nrOfPoints =  (len(node_coords)) - oldLen
    #print (nrOfPoints)

    #adding mech boundary conditions
    for n in range ( nrOfPoints ):
        mBC = utilitiesMech.mechanicalBC(dim, oldLen + n, lineBC)
        mechBC_merged.append(mBC)
        #print('adding')



    ###############generating of nodes, supported top  vertical ###############
    pointGenerators.generateNodesLine2dRand(nodeA, nodeB, minDist, dim, node_coords,  trials, False, False)

    ###############generating of nodes, supported bottom  vertical ###############
    pointGenerators.generateNodesLine2dRand(nodeA, nodeB, minDist, dim, node_coords,  trials, False, False)





    ##########################################generating of points, homogeneous volume
    rectBC = np.array([-1,-1,-1,-1,-1,-1])
    #rect
    oldLen = len(node_coords)
    pointGenerators.generateNodesRect(maxLim, minDist, dim, trials, node_coords)
    #
    newLen = len(node_coords)-1

   # print (nrOfPoints)
  #  mBC = utilitiesGeom.mechanicalBC(dim, kvadrBC, oldLen, newLen)
   # mechBC_merged.append(mBC)
    ####################################################################################################

    return node_coords,  mechBC_merged, mechIC_merged


def assemble2dbeamConfinedPress (maxLim, minDist, trials):
    dim = 2
    #lists for the model
    node_coords = []
    mechBC_merged = []
    mechIC_merged = []

    #an indent due to mirroring of the data for voronoi tess.
    indent = 1e-8

    nodeA = np.array ( [ 0 + indent, indent ] )
    nodeAmechBC = np.array([0, 0 , -1 ,      -1 , -1 , -1])
    pointGenerators.generateSingleNode(nodeA, dim, node_coords)
    mechBC_merged.append( utilitiesMech.mechanicalBC(dim, 0, nodeAmechBC))


    nodeB = np.array ( [ 0 + indent, maxLim[1]-indent ] )
    nodeBmechBC = np.array([0, 0 , -1 ,      -1 , -1 , -1])
    pointGenerators.generateSingleNode(nodeB, dim, node_coords)
    mechBC_merged.append(utilitiesMech.mechanicalBC(dim, 1, nodeBmechBC))

    ###############generating of nodes, supported line left vertical ###############
    #mech bc
    lineBC = np.array([0,-1,-1,-1,-1,-1])

    #defining points of the line
    nodeA = np.array([indent, indent])
    nodeB = np.array([indent, maxLim[1]-indent])

    oldLen = len(node_coords)
    pointGenerators.generateNodesLine2dRand(nodeA, nodeB, minDist, dim, node_coords,  trials, False, True)
    nrOfPoints =  (len(node_coords)) - oldLen

    #adding mech boundary conditions
    for n in range ( nrOfPoints ):
        mBC = utilitiesMech.mechanicalBC(dim, oldLen + n, lineBC)
        mechBC_merged.append(mBC)

    ###############generating a nodes on right face, loaded by uni tens in X) ###############
    lineBC = np.array([1, 0 ,-1,   -1, -1 ,-1])

    #defining points of the line
    nodeA = np.array([maxLim[0] - indent, indent])
    nodeB = np.array([maxLim[0] - indent, maxLim[1]-indent])

    oldLen = len(node_coords)
    pointGenerators.generateNodesLine2dRand(nodeA, nodeB, minDist, dim, node_coords,  trials, True, True)
    nrOfPoints =  (len(node_coords)) - oldLen
    #print (nrOfPoints)

    #adding mech boundary conditions
    for n in range ( nrOfPoints ):
        mBC = utilitiesMech.mechanicalBC(dim, oldLen + n, lineBC)
        mechBC_merged.append(mBC)
        #print('adding')

    ###############generating of nodes, supported top  vertical ###############
    #mech bc
    lineBC = np.array([-1,0,-1,-1,-1,-1])

    #defining points of the line
    nodeA = np.array([indent, maxLim[1]-indent])
    nodeB = np.array([maxLim[0]+indent, maxLim[1]-indent])

    oldLen = len(node_coords)
    pointGenerators.generateNodesLine2dRand(nodeA, nodeB, minDist, dim, node_coords,  trials, False, False)
    nrOfPoints =  (len(node_coords)) - oldLen

    #adding mech boundary conditions
    for n in range ( nrOfPoints ):
        mBC = utilitiesMech.mechanicalBC(dim, oldLen + n, lineBC)
        mechBC_merged.append(mBC)

    ###############generating of nodes, supported bottom  vertical ###############
    #mech bc
    lineBC = np.array([-1,0,-1,-1,-1,-1])

    #defining points of the line
    nodeA = np.array([indent, indent])
    nodeB = np.array([maxLim[0]+indent, indent])

    oldLen = len(node_coords)
    pointGenerators.generateNodesLine2dRand(nodeA, nodeB, minDist, dim, node_coords,  trials, False, False)
    nrOfPoints =  (len(node_coords)) - oldLen

    #adding mech boundary conditions
    for n in range ( nrOfPoints ):
        mBC = utilitiesMech.mechanicalBC(dim, oldLen + n, lineBC)
        mechBC_merged.append(mBC)




    ##########################################generating of points, homogeneous volume
    rectBC = np.array([-1,-1,-1,-1,-1,-1])
    #rect
    oldLen = len(node_coords)
    pointGenerators.generateNodesRect(maxLim, minDist, dim, trials, node_coords)
    #
    newLen = len(node_coords)-1

   # print (nrOfPoints)
  #  mBC = utilitiesGeom.mechanicalBC(dim, kvadrBC, oldLen, newLen)
   # mechBC_merged.append(mBC)
    ####################################################################################################

    return node_coords,  mechBC_merged, mechIC_merged


def assemblePatchTestTransport (maxLim, minDist, trials, dim):
    #lists for the model
    node_coords = np.zeros((0,dim))
    radii = np.zeros(0)
    mechBC_merged = []
    mechIC_merged = []

    ##########################################generating of points, homogeneous volume
    rectBC = np.array([-1,-1,-1,-1,-1,-1])
    #rect
    oldLen = len(node_coords)
    #pointGenerators.generateNodesRect(maxLim, minDist, dim, trials, node_coords)
    node_coords, radii = pointGenerators.generateParticlesRect(maxLim, minDist/4., minDist, 0.8, dim, trials, node_coords, radii)
    #
    newLen = len(node_coords)-1

   # print (nrOfPoints)
  #  mBC = utilitiesGeom.mechanicalBC(dim, kvadrBC, oldLen, newLen)
   # mechBC_merged.append(mBC)
    ####################################################################################################









    return node_coords, radii, mechBC_merged, mechIC_merged



def assemble2DSSBeamBending (maxLim, minDist, trials, notch, loadWidth,
                             fracZoneWidth,  orthogonalFracZone=False,
                             notchWidth = -1, node_coords_init=None, coupled=False, specifiedNodes=[], loading='3pb',gapWidth=0):
    dim = 2
    #lists for the model
    if node_coords_init is None:
        node_coords = []
    else:
        node_coords = node_coords_init

    #lists for the model

    mechBC_merged = []
    mechInitC_merged = []
    govNodes = []
    govNodesMechBC = []
    rigidPlates = []

    indent = 1e-8

    node_coords.append(np.array([maxLim[0]/2-0.01, maxLim[1] - indent]))
    node_coords.append(np.array([maxLim[0]/2+0.01, maxLim[1] - indent]))

    #an indent due to mirroring of the data for voronoi tess.
    notches=[]

    # the following is for remesher that works so far just for the bema WITHOUT notch
    if node_coords_init is None:
        if (len(specifiedNodes)>0):
            print ('appending specified nodes...')
            for node in specifiedNodes:
                node_coords.append((node))


        #generating notch points
        if (notch > 0):
            notchSide0 = []
            nodeA = np.array([maxLim[0]/2-notchWidth, indent])
            nodeB = np.array([maxLim[0]/2-notchWidth, maxLim[1]*notch-minDist])
            oldLen = len(node_coords)
            pointGenerators.generateNodesLine2dRand(nodeA, nodeB, minDist, dim, node_coords, trials, catchCorners=True, equidist=True)
            for i in range (oldLen, len(node_coords), 1):
                notchSide0.append(i)

            notchSide1 = []
            nodeA = np.array([maxLim[0]/2+notchWidth, indent])
            nodeB = np.array([maxLim[0]/2+notchWidth, maxLim[1]*notch-minDist])
            oldLen = len(node_coords)
            pointGenerators.generateNodesLine2dRand(nodeA, nodeB, minDist, dim, node_coords, trials, catchCorners=True, equidist=True)
            for i in range (oldLen, len(node_coords), 1):
                notchSide1.append(i)

            notchA = []
            notchA.append(notchSide0)
            notchA.append(notchSide1)
            notches.append(notchA)

            if not orthogonalFracZone:
                node_coords.append(np.array([maxLim[0]/2+notchWidth, maxLim[1]*notch]))
                node_coords.append(np.array([maxLim[0]/2-notchWidth, maxLim[1]*notch]))

    #width of the supports
    supportWidth = maxLim[0] / 21

    print('maxlim ',maxLim)
    ##################### CONSTRAINTS AND RIGID PLATES
    #rigid plate left support
    indentRP = 1e-3
    leftRigidPlateMechBC = np.array([0,0,-1, -1,-1,-1])
    leftRigidPlate = utilitiesMech.RigidPlate(-1, 2, np.array([ -indentRP, supportWidth/2+indentRP, -indentRP, 2*indentRP ]))
    rigidPlates.append(leftRigidPlate)
    govNodes.append(np.array([ indent+supportWidth/2, indent ]))
    govNodesMechBC.append(utilitiesMech.mechanicalBC(dim, -1, leftRigidPlateMechBC))
    #rigid plate left support
    rightRigidPlateMechBC = np.array([-1,0,-1, -1,-1,-1])
    rightRigidPlate = utilitiesMech.RigidPlate(-2, 2, np.array([ maxLim[0] - 2*indentRP-supportWidth/2, maxLim[0] + indentRP, -indentRP, 2*indentRP ]))
    rigidPlates.append(rightRigidPlate)
    govNodes.append(np.array([ maxLim[0]-indent-supportWidth/2, indent ]))
    govNodesMechBC.append(utilitiesMech.mechanicalBC(dim, -2, rightRigidPlateMechBC))

    print('support plate ',np.array([ -indentRP, supportWidth/2+indentRP, -indentRP, 2*indentRP ]))
    print('support plate ',np.array([ maxLim[0] - 2*indentRP-supportWidth/2, maxLim[0] + indentRP, -indentRP, 2*indentRP ]))



    if loading =='3pb':
        if gapWidth>0:
            topRigidPlateMechBC = np.array([-1,1,-1, -1,-1,-1])
            topRigidPlate = utilitiesMech.RigidPlate(-3, 2,
            np.array([
            0.5*maxLim[0]-notchWidth-gapWidth,
            0.5*maxLim[0]-notchWidth,
            maxLim[1] - 2*indentRP,
            maxLim[1] + 2*indentRP  ]))
            rigidPlates.append(topRigidPlate)
            govNodes.append(np.array([ maxLim[0]/2-notchWidth-gapWidth/2, maxLim[1]-indent ]))
            govNodesMechBC.append(utilitiesMech.mechanicalBC(dim, -3, topRigidPlateMechBC))

            rtopRigidPlateMechBC = np.array([-1,1,-1, -1,-1,-1])
            rtopRigidPlate = utilitiesMech.RigidPlate(-4, 2,
            np.array([
            0.5*maxLim[0]+notchWidth,
            0.5*maxLim[0]+notchWidth+gapWidth,
            maxLim[1] - 2*indentRP,
            maxLim[1] + 2*indentRP  ]))
            rigidPlates.append(rtopRigidPlate)
            govNodes.append(np.array([ maxLim[0]/2+notchWidth+gapWidth/2, maxLim[1]-indent ]))
            govNodesMechBC.append(utilitiesMech.mechanicalBC(dim, -4, rtopRigidPlateMechBC))

            topRigidPlateMechBC = np.array([-1,-1,-1, -1,2,-1])
            topRigidPlate = utilitiesMech.RigidPlate(-5, 2,
            np.array([
            0.5*maxLim[0]-notchWidth-gapWidth,
            0.5*maxLim[0]-notchWidth,
            - 2*indentRP,
             2*indentRP  ]))
            rigidPlates.append(topRigidPlate)
            govNodes.append(np.array([ maxLim[0]/2-notchWidth-gapWidth/2, -indent ]))
            govNodesMechBC.append(utilitiesMech.mechanicalBC(dim, -5, topRigidPlateMechBC))

            rtopRigidPlateMechBC = np.array([-1,-1,-1, -1,2,-1])
            rtopRigidPlate = utilitiesMech.RigidPlate(-6, 2,
            np.array([
            0.5*maxLim[0]+notchWidth,
            0.5*maxLim[0]+notchWidth+gapWidth,
             - 2*indentRP,
             + 2*indentRP  ]))
            rigidPlates.append(rtopRigidPlate)
            govNodes.append(np.array([ maxLim[0]/2+notchWidth+gapWidth/2, -indent ]))
            govNodesMechBC.append(utilitiesMech.mechanicalBC(dim, -6, rtopRigidPlateMechBC))


        else:
            topRigidPlateMechBC = np.array([-1,1,-1, -1,-1,-1])
            topRigidPlate = utilitiesMech.RigidPlate(-3, 2,
            np.array([
            0.5*maxLim[0]*(1-loadWidth)-indentRP,
            maxLim[0]  - 0.5*maxLim[0]*(1-loadWidth)+indentRP,
            maxLim[1] - 2*indentRP,
            maxLim[1] + 2*indentRP  ]))
            rigidPlates.append(topRigidPlate)
            govNodes.append(np.array([ maxLim[0]/2, maxLim[1]-indent ]))
            govNodesMechBC.append(utilitiesMech.mechanicalBC(dim, -3, topRigidPlateMechBC))


    if loading =='4pb':
        #rigid plate top load
        topRigidPlateMechBC = np.array([-1,-1,-1, -1,5,-1])
        topRigidPlate = utilitiesMech.RigidPlate(-3, 2,
        np.array([
         maxLim[0]/3*(1-loadWidth)-indentRP,
         maxLim[0]/3*(1+loadWidth)+indentRP,
        maxLim[1] - 2*indentRP,
        maxLim[1] + 2*indentRP  ]))
        rigidPlates.append(topRigidPlate)
        govNodes.append(np.array([ maxLim[0]/3, maxLim[1]-indent ]))
        govNodesMechBC.append(utilitiesMech.mechanicalBC(dim, -3, topRigidPlateMechBC))

        ltopRigidPlateMechBC = np.array([-1,-1,-1, -1,5,-1])
        ltopRigidPlate = utilitiesMech.RigidPlate(-4, 2,
        np.array([
        maxLim[0]/3*2*(1-loadWidth)-indentRP,
        maxLim[0]/3*2*(1+loadWidth)+indentRP,
        maxLim[1] - 2*indentRP,
        maxLim[1] + 2*indentRP  ]))
        rigidPlates.append(ltopRigidPlate)
        govNodes.append(np.array([ maxLim[0]/3*2, maxLim[1]-indent ]))
        govNodesMechBC.append(utilitiesMech.mechanicalBC(dim, -4, ltopRigidPlateMechBC))

        node_coords.append(np.array([maxLim[0]/3 ,maxLim[1]-indent/2]))
        node_coords.append(np.array([maxLim[0]/3+loadWidth/3 ,maxLim[1]-indent/2]))
        node_coords.append(np.array([maxLim[0]/3-loadWidth/3 ,maxLim[1]-indent/2]))

        node_coords.append(np.array([maxLim[0]/3*2, maxLim[1]-indent/2]))
        node_coords.append(np.array([maxLim[0]/3*2+loadWidth/3 ,maxLim[1]-indent/2]))
        node_coords.append(np.array([maxLim[0]/3*2-loadWidth/3 ,maxLim[1]-indent/2]))


    if node_coords_init is None:
        ###############generating of nodes, left horizontal support ###############
        #defining points of the line
        nodeA = np.array([indent, indent])
        nodeB = np.array([indent + supportWidth, indent])
        pointGenerators.generateNodesLine2dRand(nodeA, nodeB, minDist, dim, node_coords, trials, True, False)
        ###############generating of nodes, right horizontal support ###############
        #defining points of the line
        nodeA = np.array([maxLim[0] - supportWidth -indent, indent])
        nodeB = np.array([maxLim[0] - indent, indent])
        pointGenerators.generateNodesLine2dRand(nodeA, nodeB, minDist, dim, node_coords, trials, True, False)
        ############### loaded top face ###############
        if gapWidth>0:
            nodeA =  np.array([0.5*maxLim[0]-notchWidth-gapWidth, maxLim[1] - indent])
            nodeB =  np.array([0.5*maxLim[0]-notchWidth-gapWidth , maxLim[1] - indent])
            pointGenerators.generateNodesLine2dRand(nodeA, nodeB, minDist, dim, node_coords,  trials, False, False)
        else:
            nodeA =  np.array([indent + 0.5*maxLim[0]*(1-loadWidth) , maxLim[1] - indent])
            nodeB =  np.array([maxLim[0] - indent - 0.5*maxLim[0]*(1-loadWidth) , maxLim[1] - indent])
            #pointGenerators.generateNodesLine2dRand(nodeA, nodeB, minDist, dim, node_coords,  trials, False, False)
            #######top of frac zone
            nodeA =  np.array([0.5*maxLim[0] - 0.5*maxLim[1]*(1-notch) , maxLim[1] - indent])
            nodeB =  np.array([0.5*maxLim[0] + 0.5*maxLim[1]*(1-notch), maxLim[1] - indent])
            #pointGenerators.generateNodesLine2dRand(nodeA, nodeB, minDist, dim, node_coords,  trials, False, False)


        ########################################## rest of  faces
        nodeA =  np.array([indent + 0.5*maxLim[0]*(1-loadWidth) , maxLim[1] - indent])
        nodeB =  np.array([indent, maxLim[1] - indent])
        pointGenerators.generateNodesLine2dRand(nodeA, nodeB, minDist*2, dim, node_coords,  trials, False, False)
        nodeA =  np.array([maxLim[0] - indent  , maxLim[1] - indent])
        nodeB =  np.array([maxLim[0] - indent - 0.5*maxLim[0]*(1-loadWidth) , maxLim[1] - indent])
        pointGenerators.generateNodesLine2dRand(nodeA, nodeB, minDist*2, dim, node_coords,  trials, False, False)

        #bottom line
        if gapWidth>0:
            oldl = len(node_coords)
            nodeA =  np.array([maxLim[0]/2-notchWidth-gapWidth  ,  indent])
            nodeB =  np.array([maxLim[0]/2-notchWidth ,  indent])
            pointGenerators.generateNodesLine2dRand(nodeA, nodeB, minDist, dim, node_coords,  trials, True,True)
            nodeA =  np.array([maxLim[0]/2+notchWidth ,  indent])
            nodeB =  np.array([maxLim[0]/2+notchWidth+gapWidth ,  indent])
            pointGenerators.generateNodesLine2dRand(nodeA, nodeB, minDist, dim, node_coords,  trials, True,True)
            #newl = len(node_coords)
            #for i in range(newl-oldl):
            #        mechBC = np.array([-1,-1,-1,-1,2,-1])
            #        mBC = utilitiesMech.mechanicalBC(dim, oldl+i, mechBC)
            #            mechBC_merged.append(mBC)

        nodeA =  np.array([indent  ,  indent])
        nodeB =  np.array([maxLim[0]/2-notchWidth ,  indent])
        pointGenerators.generateNodesLine2dRand(nodeA, nodeB, minDist*2, dim, node_coords,  trials, False, False)
        nodeA =  np.array([maxLim[0]-indent  ,  indent])
        nodeB =  np.array([maxLim[0]/2+notchWidth ,  indent])
        pointGenerators.generateNodesLine2dRand(nodeA, nodeB, minDist*2, dim, node_coords,  trials, False, False)

        nodeA =  np.array([indent, indent])
        nodeB =  np.array([indent, maxLim[1] - indent])
        pointGenerators.generateNodesLine2dRand(nodeA, nodeB, minDist*2, dim, node_coords,  trials, False, False)
        nodeA =  np.array([maxLim[0]-indent, indent])
        nodeB =  np.array([maxLim[0]-indent, maxLim[1] - indent])
        pointGenerators.generateNodesLine2dRand(nodeA, nodeB, minDist*2, dim, node_coords,  trials, False, False)

        """
        nodeA =  np.array([indent  ,  indent])
        nodeB =  np.array([indent ,  maxLim[1]-indent])
        pointGenerators.generateNodesLine2dRand(nodeA, nodeB, minDist/1, dim, node_coords,  trials, False, False)
        nodeA =  np.array([maxLim[0] - indent   ,  indent])
        nodeB =  np.array([maxLim[0] - indent  ,  maxLim[1]-indent])
        pointGenerators.generateNodesLine2dRand(nodeA, nodeB, minDist/1, dim, node_coords,  trials, False, False)
        #"""
        ##########################################generating of points, fracture zone

        maxLimF = np.array([
        maxLim[0] - indent - 0.5*maxLim[0]*(1-fracZoneWidth),
        maxLim[1] - indent,
        indent + 0.5*maxLim[0]*(1-fracZoneWidth),
        maxLim[1]*notch*0.8])

        if gapWidth>0:
            maxLimF = np.array([
            maxLim[0]/2 +notchWidth+gapWidth*2,
            maxLim[1] - indent,
            maxLim[0]/2 -notchWidth-gapWidth*2,
            indent])


        if not orthogonalFracZone:
            if not fracZoneWidth ==1:
                pointGenerators.generateNodesRect(maxLimF, minDist, dim, trials, node_coords, useLowBound=True)
            else:
                pointGenerators.generateNodesRect(maxLim, minDist, dim, trials, node_coords)
        else:
            maxLimF = np.array([
            maxLim[0],
            maxLim[1] - indent,
            maxLim[0]*0.5 + maxLim[1]*(1-notch)/2,
            maxLim[1]*notch])

            pointGenerators.generateOrtogrid(maxLimF, minDist, dim, node_coords, maxLim[1]*(1-notch))
            #fracZoneWidth*maxLim[0])


        ## notch faces

        """
        maxLimF = np.array([
        maxLim[0]/2 - 0.5*maxLim[1]*(1-notch)*1.5,
        maxLim[1],
        maxLim[0]/2 + 0.5*maxLim[1]*(1-notch)*1.5,
        indent+maxLim[1]*notch/2])
        pointGenerators.generateNodesRect(maxLimF, minDist, dim, trials, node_coords, useLowBound=True)
        """

        """
        maxLimF = np.array([
        maxLim[0]/2 - 0.5*maxLim[1]*(1-notch)*2.5,
        maxLim[1],
        maxLim[0]/2 + 0.5*maxLim[1]*(1-notch)*2.5,
        indent+maxLim[1]*notch/2])
        pointGenerators.generateNodesRect(maxLimF, minDist*1.5, dim, trials, node_coords, useLowBound=True)

        maxLimF = np.array([
        maxLim[0] - indent - 0.5*maxLim[0]*(1-fracZoneWidth) - maxLim[0]*fracZoneWidth,
        maxLim[1]-indent,
        indent + 0.5*maxLim[0]*(1-fracZoneWidth) - maxLim[0]*fracZoneWidth,
        indent])
        pointGenerators.generateNodesRect(maxLimF, minDist/2, dim, trials, node_coords, useLowBound=True)
        #"""

        ##########################################generating of points, left support
        maxLimF = np.array([
        supportWidth*2,
        supportWidth*2,
        indent,
        indent])
        pointGenerators.generateNodesRect(maxLimF, minDist*1.2, dim, trials, node_coords, useLowBound=True)
        ##########################################generating of points, right support
        maxLimF = np.array([
        maxLim[0],
        supportWidth*2,
        maxLim[0]-supportWidth*2,
        indent])
        pointGenerators.generateNodesRect(maxLimF, minDist*1.2, dim, trials, node_coords, useLowBound=True)


        #rect
        pointGenerators.generateNodesRect(maxLim, minDist*2.5, dim, trials, node_coords)


    """
    node_coords = np.asarray(node_coords)
    plt.plot(node_coords[:,0], node_coords[:,1], 'o', color='black');
    plt.show()
    """

    return node_coords, mechBC_merged, mechInitC_merged, notches, govNodes, govNodesMechBC, rigidPlates


def assemble2d_CFRAC_Clover(maxLim, minDist, trials, holeMinDist, holeDiameter, baseMinDist,fineWidth,fineHeight):
    print('2d clover')
    dim = 2
    node_coords = []
    mechBC_merged = []
    mechInitC_merged = []
    govNodes = []
    govNodesMechBC = []
    rigidPlates = []
    functions = []

    #an indent due to mirroring of the data for voronoi tess.
    indent = 1e-8
    indentRP = 1e-7

    #sampling holes
    for r in range (4):
        oldLen = len(node_coords)
        interface = []

        print ('Hole #%d' %r)
        if r == 0:
            centre = np.array([ maxLim[0]/2-0.011, maxLim[1]/2 ])
        if r == 1:
            centre = np.array([ maxLim[0]/2+0.011, maxLim[1]/2 ])
        if r == 2:
            centre = np.array([ maxLim[0]/2, maxLim[1]/2-0.011 ])
        if r == 3:
            centre = np.array([ maxLim[0]/2, maxLim[1]/2+0.011 ])

        circleLength = 2*np.pi*holeDiameter/2
        nrNodes = int ( circleLength / holeMinDist )
        #print ('nrnodes: %d' %nrNodes)
        #deleno ctyrma
        nrNodes = (int (nrNodes / 4) +1 ) * 4

        nodesOld = len(node_coords)
        pointGenerators.generateNodesCircle2dRand(centre, holeDiameter/2, holeMinDist, node_coords, trials, equiAngNodes=nrNodes )
        newNodes = len(node_coords) - nodesOld

    sampleBorders=False
    if sampleBorders:
        #top
        nodeA = np.array([indent, maxLim[1]-indent])
        nodeB = np.array([maxLim[0]-indent, maxLim[1]-indent])
        pointGenerators.generateNodesLine2dRand(nodeA, nodeB, minDist, dim, node_coords, trials, catchCorners=True, equidist=False)

        #bottom
        nodeA = np.array([indent, indent])
        nodeB = np.array([maxLim[0]-indent, indent])
        pointGenerators.generateNodesLine2dRand(nodeA, nodeB, minDist, dim, node_coords, trials, catchCorners=False, equidist=False)

        #left
        nodeA = np.array([indent, indent])
        nodeB = np.array([indent, maxLim[1]-indent])
        pointGenerators.generateNodesLine2dRand(nodeA, nodeB, minDist, dim, node_coords, trials, catchCorners=False, equidist=False)


        #right
        nodeA = np.array([maxLim[0]-indent, indent])
        nodeB = np.array([maxLim[0]-indent, maxLim[1]-indent])
        pointGenerators.generateNodesLine2dRand(nodeA, nodeB, minDist, dim, node_coords, trials, catchCorners=False, equidist=False)


    leftRigidPlateMechBC = np.array([0, 0,-1,   -1,-1,-1])
    leftRigidPlate = utilitiesMech.RigidPlate(-1, 2, np.array([-indentRP, 0.015,0.007, 0.025 ]))
    rigidPlates.append(leftRigidPlate)
    govNodes.append(np.array([ 0.0075, (0.025+0.007)/2 ]))
    govNodesMechBC.append(utilitiesMech.mechanicalBC(dim, -1, leftRigidPlateMechBC))

    rightRigidPlateMechBC = np.array([1, 0,-1,   -1,-1,-1])
    rightRigidPlate = utilitiesMech.RigidPlate(-1, 2, np.array([maxLim[0]-0.015,maxLim[0]+indentRP,    0.007, 0.025 ]))
    rigidPlates.append(rightRigidPlate)
    govNodes.append(np.array([ maxLim[0]-0.0075, (0.025+0.007)/2 ]))
    govNodesMechBC.append(utilitiesMech.mechanicalBC(dim, -2, rightRigidPlateMechBC))







    #pointGenerators.generateNodesRect(maxLim, minDist, 2, trials, node_coords)

    #original model
    #left gradient
    interBounds = np.array([     indent,      indent , maxLim[0]/3,              maxLim[1]])
    pointGenerators.generateNodesRect(interBounds, minDist, dim, trials, node_coords, useLowBound=True, topMinDist = baseMinDist, bottomMinDist = baseMinDist/2, gradienDirection=0)
    interBounds = np.array([     maxLim[0]/3,      indent , maxLim[0]/3*2,              maxLim[1]])
    pointGenerators.generateNodesRect(interBounds, baseMinDist/2, dim, trials, node_coords, useLowBound=True)
    #right gradient
    interBounds = np.array([     maxLim[0]/3*2,      indent , maxLim[0],              maxLim[1]])
    pointGenerators.generateNodesRect(interBounds, minDist, dim, trials, node_coords, useLowBound=True, topMinDist = baseMinDist/2, bottomMinDist = baseMinDist, gradienDirection=0)

    # sparse model
    #bounds = np.array([     indent,      indent , maxLim[0],              maxLim[1]])
    #pointGenerators.generateNodesRect(bounds, baseMinDist, dim, trials, node_coords, useLowBound=True)

    #left gradient
    interBounds = np.array([     maxLim[0]/2-fineWidth-5*baseMinDist,      indent , maxLim[0]/2-fineWidth,              indent + fineHeight ])
    pointGenerators.generateNodesRect(interBounds, minDist, dim, trials, node_coords, useLowBound=True, topMinDist = baseMinDist, bottomMinDist = minDist, gradienDirection=0)

    #right gradient
    interBounds = np.array([     maxLim[0]/2+fineWidth,      indent , maxLim[0]/2+fineWidth+5*baseMinDist,              indent + fineHeight ])
    pointGenerators.generateNodesRect(interBounds, minDist, dim, trials, node_coords, useLowBound=True, topMinDist = minDist, bottomMinDist = baseMinDist, gradienDirection=0)

    #fine center
    bounds = np.array([      maxLim[0]/2-fineWidth,      indent ,  maxLim[0]/2+fineWidth,              indent + fineHeight])
    pointGenerators.generateNodesRect(bounds, minDist, dim, trials, node_coords, useLowBound=True)

    #top gradient
    interBounds = np.array([      maxLim[0]/2-fineWidth,      indent + fineHeight ,  maxLim[0]/2+fineWidth,              fineHeight+5*baseMinDist])
    pointGenerators.generateNodesRect(interBounds, minDist, dim, trials, node_coords, useLowBound=True, topMinDist = minDist , bottomMinDist = baseMinDist, gradienDirection=1)

    #remove points from rebars
    newNodes = []
    for i in range (len(node_coords)):
        app = True
        for r in range(4):

            if r == 0:
                centre = np.array([ maxLim[0]/2-0.011, maxLim[1]/2 ])
                if (np.linalg.norm(centre - node_coords[i]) < holeDiameter/2*0.999):
                    app = False
            if r == 1 and app==True:
                centre = np.array([ maxLim[0]/2+0.011, maxLim[1]/2])
                if not (np.linalg.norm(centre - node_coords[i]) > holeDiameter/2*0.999):
                    app = False
            if r == 2and app==True:
                centre = np.array([ maxLim[0]/2, maxLim[1]/2-0.011 ])
                if not (np.linalg.norm(centre - node_coords[i]) > holeDiameter/2*0.999):
                    app = False
            if r == 3and app==True:
                centre = np.array([ maxLim[0]/2, maxLim[1]/2+0.011 ])
                if not (np.linalg.norm(centre - node_coords[i]) > holeDiameter/2*0.999):
                    app = False

        if app==True:
            newNodes.append(node_coords[i])

    node_coords = newNodes.copy()



    node_coords = np.asarray(node_coords)
    if False:
        plt.plot(node_coords[:,0], node_coords[:,1], 'o', color='black');
        plt.show()


    return node_coords, mechBC_merged, mechInitC_merged,  govNodes, govNodesMechBC, rigidPlates, functions



def assemble3d_CFRAC_Clover(maxLim, minDist, trials, holeMinDist, holeDiameter):

    dim = 3
    node_coords = []
    mechBC_merged = []
    mechInitC_merged = []
    govNodes = []
    govNodesMechBC = []
    rigidPlates = []
    functions = []

    #an indent due to mirroring of the data for voronoi tess.
    indent = 1e-8
    indentRP = 1e-7

    #sampling holes
    for r in range (4):
        oldLen = len(node_coords)
        interface = []

        print ('Hole 3d #%d' %r)
        if r == 0:
            centre = np.array([ maxLim[0]/2-0.011, maxLim[1]/2, indent ])
        if r == 1:
            centre = np.array([ maxLim[0]/2+0.011, maxLim[1]/2,indent ])
        if r == 2:
            centre = np.array([ maxLim[0]/2, maxLim[1]/2-0.011,indent ])
        if r == 3:
            centre = np.array([ maxLim[0]/2, maxLim[1]/2+0.011,indent ])

        circleLength = 2*np.pi*holeDiameter/2
        nrNodes = int ( circleLength / holeMinDist )
        #print ('nrnodes: %d' %nrNodes)
        #deleno ctyrma
        nrNodes = (int (nrNodes / 4) +1 ) * 4

        nodesOld = len(node_coords)
        pointGenerators.generateNodesOrtoCilinderSurf3dRand(centre, holeDiameter/2,  maxLim[2], 2 ,       minDist, node_coords, trials)

        newNodes = len(node_coords) - nodesOld

    sampleBorders=False
    if sampleBorders:
        #front
        nodeA = np.array([indent, indent, indent])
        nodeB = np.array([maxLim[0]-indent, maxLim[1]-indent, maxLim[2]-indent])
        pointGenerators.generateNodesOrtoSurface3dRand(nodeA, nodeB, minDist, dim, node_coords, trials)

        #back
        nodeA = np.array([indent, indent, maxLim[2]-indent])
        nodeB = np.array([maxLim[0]-indent, maxLim[1]-indent, maxLim[2]-indent])
        pointGenerators.generateNodesOrtoSurface3dRand(nodeA, nodeB, minDist, dim, node_coords, trials)

        #top
        nodeA = np.array([indent, maxLim[1]-indent, indent])
        nodeB = np.array([maxLim[0]-indent, maxLim[1]-indent, maxLim[2]-indent])
        pointGenerators.generateNodesOrtoSurface3dRand(nodeA, nodeB, minDist, dim, node_coords, trials)

        #bottom
        nodeA = np.array([indent, indent,indent])
        nodeB = np.array([maxLim[0]-indent, indent,maxLim[2]-indent])
        pointGenerators.generateNodesOrtoSurface3dRand(nodeA, nodeB, minDist, dim, node_coords, trials)

        #left
        nodeA = np.array([indent, indent,indent])
        nodeB = np.array([indent, maxLim[1]-indent, maxLim[2]-indent])
        pointGenerators.generateNodesOrtoSurface3dRand(nodeA, nodeB, minDist, dim, node_coords, trials)

        #right
        nodeA = np.array([maxLim[0]-indent, indent,indent])
        nodeB = np.array([maxLim[0]-indent, maxLim[1]-indent,maxLim[2]-indent])
        pointGenerators.generateNodesOrtoSurface3dRand(nodeA, nodeB, minDist, dim, node_coords, trials)

    leftRigidPlateMechBC = np.array([0, 0,0,   -1,-1,-1, -1,-1,-1, -1,-1,-1])
    leftRigidPlate = utilitiesMech.RigidPlate(-1, 3, np.array([-indentRP, 0.015,0.007, 0.025, -indentRP, maxLim[2]+indentRP ]))
    rigidPlates.append(leftRigidPlate)
    govNodes.append(np.array([ 0.0075, (0.025+0.007)/2,maxLim[2]/2 ]))
    govNodesMechBC.append(utilitiesMech.mechanicalBC(dim, -1, leftRigidPlateMechBC))

    rightRigidPlateMechBC = np.array([1, 0,0,   -1,-1,-1, -1,-1,-1, -1,-1,-1])
    rightRigidPlate = utilitiesMech.RigidPlate(-1, 3, np.array([maxLim[0]-0.015,maxLim[0]+indentRP,    0.007, 0.025, -indentRP, maxLim[2]+indentRP ]))
    rigidPlates.append(rightRigidPlate)
    govNodes.append(np.array([ maxLim[0]-0.0075, (0.025+0.007)/2,maxLim[2]/2 ]))
    govNodesMechBC.append(utilitiesMech.mechanicalBC(dim, -2, rightRigidPlateMechBC))

    #left
    interBounds = np.array([     indent,      indent , indent, maxLim[0]/3,   maxLim[1], maxLim[2]])
    pointGenerators.generateNodesRect(interBounds, minDist, dim, trials, node_coords, useLowBound=True, topMinDist = minDist*2, bottomMinDist = minDist, gradienDirection=0)

    interBounds = np.array([     maxLim[0]/3,      indent , indent, maxLim[0]/3*2,   maxLim[1],  maxLim[2]])
    pointGenerators.generateNodesRect(interBounds, minDist, dim, trials, node_coords, useLowBound=True)

    interBounds = np.array([     maxLim[0]/3*2,      indent , indent, maxLim[0],  maxLim[1],  maxLim[2]])
    pointGenerators.generateNodesRect(interBounds, minDist, dim, trials, node_coords, useLowBound=True, topMinDist = minDist, bottomMinDist = minDist*2, gradienDirection=0)


    #remove points from rebars
    newNodes = []
    for i in range (len(node_coords)):
        app = True
        for r in range(4):

            if r == 0:
                centre = np.array([ maxLim[0]/2-0.011, maxLim[1]/2 ])
                if (np.linalg.norm(centre [0:2]- node_coords[i][0:2]) < holeDiameter/2*0.999):
                    app = False
            if r == 1 and app==True:
                centre = np.array([ maxLim[0]/2+0.011, maxLim[1]/2])
                if not (np.linalg.norm(centre[0:2] - node_coords[i][0:2]) > holeDiameter/2*0.999):
                    app = False
            if r == 2and app==True:
                centre = np.array([ maxLim[0]/2, maxLim[1]/2-0.011 ])
                if not (np.linalg.norm(centre[0:2] - node_coords[i][0:2]) > holeDiameter/2*0.999):
                    app = False
            if r == 3and app==True:
                centre = np.array([ maxLim[0]/2, maxLim[1]/2+0.011 ])
                if not (np.linalg.norm(centre[0:2] - node_coords[i][0:2]) > holeDiameter/2*0.999):
                    app = False

        if app==True:
            newNodes.append(node_coords[i])

    node_coords = newNodes.copy()



    node_coords = np.asarray(node_coords)
    if False:
        plt.plot(node_coords[:,0], node_coords[:,1], 'o', color='black');
        plt.show()


    return node_coords, mechBC_merged, mechInitC_merged,  govNodes, govNodesMechBC, rigidPlates, functions




def assemble2d_CFRAC_TDCB(maxLim, minDist, trials, holeMinDist, holeDiameter,roughMinDistCoef,elazonewidth):
    roughCoef = roughMinDistCoef
    dim = 2
    node_coords = []
    mechBC_merged = []
    mechInitC_merged = []
    govNodes = []
    govNodesMechBC = []
    rigidPlates = []
    functions = []

    #an indent due to mirroring of the data for voronoi tess.
    indent = 1e-8
    indentRP = 1e-7

    #sampling holes
    idcs_leftHole = []
    idcs_rightHole = []
    for r in range (2):
        print ('Hole #%d' %r)
        if r == 0:
            centre = np.array([ maxLim[0]/2-0.01, 0.011 ])
        if r == 1:
            centre = np.array([ maxLim[0]/2+0.01, 0.011 ])

        circleLength = 2*np.pi*holeDiameter/2
        nrNodes = int ( circleLength / holeMinDist )
        #print ('nrnodes: %d' %nrNodes)
        #deleno ctyrma
        nrNodes = (int (nrNodes / 4) +1 ) * 4

        nodesOld = len(node_coords)
        pointGenerators.generateNodesCircle2dRand(centre, holeDiameter/2, holeMinDist, node_coords, trials, equiAngNodes=nrNodes )
        newNodes = len(node_coords) - nodesOld
        if r == 0:
            idcs_leftHole = np.arange(0,newNodes)+nodesOld
            idcs_leftHole = idcs_leftHole.tolist()
            for h in idcs_leftHole:
                if (node_coords[h][0] > (centre[0]-0.7071*holeDiameter/2)) :
                    idcs_leftHole.remove(h)

        if r == 1:
            idcs_rightHole = np.arange(0,newNodes)+nodesOld
            idcs_rightHole = idcs_rightHole.tolist()
            for h in idcs_rightHole:
                if (node_coords[h][0] < (centre[0]+0.7071*holeDiameter/2)) :
                    idcs_rightHole.remove(h)

    leftRigidPlateMechBC = np.array([0, 0,-1,   -1,-1,-1])
    leftRigidPlate = utilitiesMech.RigidPlate(-1, 2, None, directIdcs=True )
    leftRigidPlate.setDirectNodes(idcs_leftHole)
    rigidPlates.append(leftRigidPlate)
    govNodes.append(np.array([ maxLim[0]/2-0.01, 0.011 ]))
    govNodesMechBC.append(utilitiesMech.mechanicalBC(dim, -1, leftRigidPlateMechBC))

    rightRigidPlateMechBC = np.array([1, 0,-1,   -1,-1,-1])
    rightRigidPlate = utilitiesMech.RigidPlate(-1, 2, None, directIdcs=True )
    rightRigidPlate.setDirectNodes(idcs_rightHole)
    rigidPlates.append(rightRigidPlate)
    govNodes.append(np.array([ maxLim[0]/2+0.01, 0.011 ]))
    govNodesMechBC.append(utilitiesMech.mechanicalBC(dim, -2, rightRigidPlateMechBC))

    notches = []
    #generating notch points
    notch = 0.37
    if (notch > 0):
        nTop = maxLim[1]*notch+indent-minDist/2
        if nTop < minDist:
            nTop = maxLim[1]*notch+indent
        notchWidth = 0.00125
        notchWidth = minDist/2
        notchSide0 = []
        nodeA = np.array([maxLim[0]/2-notchWidth, indent])
        nodeB = np.array([maxLim[0]/2-notchWidth, nTop])
        oldLen = len(node_coords)
        pointGenerators.generateNodesLine2dRand(nodeA, nodeB, minDist, dim, node_coords, trials, catchCorners=True, equidist=True)
        for i in range (oldLen, len(node_coords), 1):
            notchSide0.append(i)

        notchSide1 = []
        nodeA = np.array([maxLim[0]/2+notchWidth, indent])
        nodeB = np.array([maxLim[0]/2+notchWidth, nTop])
        oldLen = len(node_coords)
        pointGenerators.generateNodesLine2dRand(nodeA, nodeB, minDist, dim, node_coords, trials, catchCorners=True, equidist=True)
        for i in range (oldLen, len(node_coords), 1):
            notchSide1.append(i)

        notchA = []
        notchA.append(notchSide0)
        notchA.append(notchSide1)
        notches.append(notchA)

        node_coords.append(np.array([maxLim[0]/2, nTop+minDist/2]))
        #node_coords.append(np.array([maxLim[0]/2-notchWidth, maxLim[1]*notch]))

    rougherCoef = 1.2

    sampleBorders=True
    if sampleBorders:
        #top
        nodeA = np.array([indent, maxLim[1]-indent])
        nodeB = np.array([maxLim[0]-indent, maxLim[1]-indent])
        pointGenerators.generateNodesLine2dRand(nodeA, nodeB, minDist*roughCoef, dim, node_coords, trials, catchCorners=False, equidist=False)
        #bottom
        nodeA = np.array([indent, indent])
        nodeB = np.array([maxLim[0]-indent, indent])
        pointGenerators.generateNodesLine2dRand(nodeA, nodeB, minDist*roughCoef*rougherCoef, dim, node_coords, trials, catchCorners=False, equidist=False)

        #left
        nodeA = np.array([indent, maxLim[1]-indent])
        nodeB = np.array([0.015, indent])
        pointGenerators.generateNodesLine2dRand(nodeA, nodeB, minDist*roughCoef*rougherCoef, dim, node_coords, trials, catchCorners=False, equidist=True)
        nodeA = np.array([0, maxLim[1]-indent])
        nodeB = np.array([0.015-indent, indent])
        pointGenerators.generateNodesLine2dRand(nodeA, nodeB, minDist*roughCoef*rougherCoef, dim, node_coords, trials, catchCorners=False, equidist=True)

        #right
        nodeA = np.array([0.075, indent])
        nodeB = np.array([maxLim[0]-indent, maxLim[1]-indent])
        pointGenerators.generateNodesLine2dRand(nodeA, nodeB, minDist*roughCoef*rougherCoef, dim, node_coords, trials, catchCorners=False, equidist=True)
        #right
        nodeA = np.array([0.075+indent, indent])
        nodeB = np.array([maxLim[0], maxLim[1]-indent])
        pointGenerators.generateNodesLine2dRand(nodeA, nodeB, minDist*roughCoef*rougherCoef, dim, node_coords, trials, catchCorners=False, equidist=True)



    #pointGenerators.generateNodesRect(maxLim, minDist, 2, trials, node_coords)

    fineCenterWidth = elazonewidth
    #bottom up to notch
    interBounds = np.array([     indent,      indent , maxLim[0]-indent,              maxLim[1]*(notch-2.*minDist)])
    pointGenerators.generateNodesRect(interBounds, minDist, dim, trials, node_coords, useLowBound=True, topMinDist = minDist*roughCoef, bottomMinDist = minDist*roughCoef, gradienDirection=1)

    roughtop = 0.8
    #left
    interBounds = np.array([     indent,       indent , maxLim[0]/2,              maxLim[1]*roughtop])
    pointGenerators.generateNodesRect(interBounds, minDist, dim, trials, node_coords, useLowBound=True, topMinDist = minDist*roughCoef, bottomMinDist = minDist*roughCoef/5, gradienDirection=0)
    interBounds = np.array([     maxLim[0]/3,       maxLim[1]*notch*0.9 , maxLim[0]/2-fineCenterWidth,              maxLim[1]*roughtop])
    pointGenerators.generateNodesRect(interBounds, minDist, dim, trials, node_coords, useLowBound=True, topMinDist = minDist*roughCoef, bottomMinDist = minDist, gradienDirection=0)
    interBounds = np.array([     indent,      indent , maxLim[0]/2-fineCenterWidth,              maxLim[1]*roughtop])
    pointGenerators.generateNodesRect(interBounds, minDist, dim, trials, node_coords, useLowBound=True, topMinDist = minDist*roughCoef*rougherCoef, bottomMinDist = minDist*roughCoef, gradienDirection=0)

    #right
    interBounds = np.array([     maxLim[0]/2,      indent, maxLim[0],              maxLim[1]*roughtop])
    pointGenerators.generateNodesRect(interBounds, minDist, dim, trials, node_coords, useLowBound=True, topMinDist = minDist*roughCoef/5, bottomMinDist = minDist*roughCoef, gradienDirection=0)
    interBounds = np.array([     maxLim[0]/2+fineCenterWidth,      maxLim[1]*notch*0.9 , maxLim[0]/3*2,              maxLim[1]*roughtop])
    pointGenerators.generateNodesRect(interBounds, minDist, dim, trials, node_coords, useLowBound=True, topMinDist = minDist, bottomMinDist = minDist*roughCoef, gradienDirection=0)
    interBounds = np.array([     maxLim[0]/2+fineCenterWidth,      indent , maxLim[0],              maxLim[1]*roughtop])
    pointGenerators.generateNodesRect(interBounds, minDist, dim, trials, node_coords, useLowBound=True, topMinDist = minDist*roughCoef, bottomMinDist = minDist*roughCoef*rougherCoef, gradienDirection=0)

    #center
    interBounds = np.array([     maxLim[0]/2 - fineCenterWidth,      maxLim[1]*(notch-0.02) , maxLim[0]/2+ fineCenterWidth,              maxLim[1]*roughtop])
    pointGenerators.generateNodesRect(interBounds, minDist, dim, trials, node_coords, useLowBound=True, topMinDist = minDist, bottomMinDist = minDist, gradienDirection=0)
    #center
    interBounds = np.array([     maxLim[0]/2 - fineCenterWidth,      maxLim[1]*roughtop , maxLim[0]/2 + fineCenterWidth    ,         maxLim[1]])
    pointGenerators.generateNodesRect(interBounds, minDist, dim, trials, node_coords, useLowBound=True, topMinDist = minDist*roughCoef/2, bottomMinDist = minDist*roughCoef/2, gradienDirection=0)
    interBounds = np.array([     indent,      maxLim[1]*roughtop , maxLim[0]/2 - fineCenterWidth     ,         maxLim[1]])
    pointGenerators.generateNodesRect(interBounds, minDist, dim, trials, node_coords, useLowBound=True, topMinDist = minDist*roughCoef, bottomMinDist = minDist*roughCoef/3, gradienDirection=0)
    interBounds = np.array([        maxLim[0]/2 + fineCenterWidth,      maxLim[1]*roughtop , maxLim[0]     ,         maxLim[1]])
    pointGenerators.generateNodesRect(interBounds, minDist, dim, trials, node_coords, useLowBound=True, topMinDist = minDist*roughCoef/3, bottomMinDist = minDist*roughCoef, gradienDirection=0)

    #remove points from rebars
    newNodes = []
    for i in range (len(node_coords)):
        app = True
        for r in range(2):
            if r == 0:
                centre = np.array([ maxLim[0]/2-0.01, 0.011])
                if (np.linalg.norm(centre - node_coords[i]) < holeDiameter/2*0.999):
                    app = False
            if r == 1 and app==True:
                centre = np.array([ maxLim[0]/2+0.01, 0.011])
                if not (np.linalg.norm(centre - node_coords[i]) > holeDiameter/2*0.999):
                    app = False

        if app==True:
            newNodes.append(node_coords[i])

    node_coords = newNodes.copy()



    #identifying nodes within modelu
    model_indices = []
    for i,p in enumerate(node_coords):
        app = True
        if p[0]<0.015:
            if not (p[1]>=(0.1-(6+2/3)*p[0]) ):
                app = False
        if p[0]>0.075:
            if not (p[1]>=(6+2/3)*(p[0]-0.075) ):
                app = False

        if app == True:
            model_indices.append(i)


    node_coords = np.asarray(node_coords)
    if False:
        plt.plot(node_coords[model_indices,0], node_coords[model_indices,1], 'o', color='black');
        plt.show()


    return node_coords, mechBC_merged, mechInitC_merged,  govNodes, govNodesMechBC, rigidPlates, functions, notches,model_indices



def assemble2d_Hanging_FracZone(maxLim, minDist, trials):
    dim = 2
    node_coords = []
    mechBC_merged = []
    mechInitC_merged = []
    govNodes = []
    govNodesMechBC = []
    rigidPlates = []
    functions = []

    indent = 1e-8

    print('bottom bound surface for fem frac zone')
    nodeA = np.array([indent, indent])
    nodeB = np.array([maxLim[0]-indent, indent])
    pointGenerators.generateNodesOrtoSurface3dRand(nodeA, nodeB, minDist*0.9, dim, node_coords, trials*2)
    print('top bound surface for fem frac zone')
    nodeA = np.array([indent, maxLim[1]-indent])
    nodeB = np.array([maxLim[0]-indent, maxLim[1]-indent])
    pointGenerators.generateNodesOrtoSurface3dRand(nodeA, nodeB, minDist*0.9, dim, node_coords, trials*2)
    print('left bound surface for fem frac zone')
    nodeA = np.array([indent, indent])
    nodeB = np.array([indent, maxLim[1]-indent])
    pointGenerators.generateNodesOrtoSurface3dRand(nodeA, nodeB, minDist*0.9, dim, node_coords, trials*2)
    print('right bound surface for fem frac zone')
    nodeA = np.array([maxLim[0]-indent, indent])
    nodeB = np.array([maxLim[0]-indent, maxLim[1]-indent])
    pointGenerators.generateNodesOrtoSurface3dRand(nodeA, nodeB, minDist*0.9, dim, node_coords, trials*2)

    print('fine fracture zone volume')
    pointGenerators.generateNodesRect(maxLim, minDist, dim, trials, node_coords)

    leftRigidPlateMechBC = np.array([0, 0,0,   -1,-1,-1])
    leftRigidPlate = utilitiesMech.RigidPlate(-1, 2, None, directIdcs=True )
    rigidPlates.append(leftRigidPlate)
    govNodes.append(np.array([ maxLim[0]/2-0.01, 0.011]))
    govNodesMechBC.append(utilitiesMech.mechanicalBC(dim, -1, leftRigidPlateMechBC))

    rightRigidPlateMechBC = np.array([1, 0,0,   -1,-1,-1])
    rightRigidPlate = utilitiesMech.RigidPlate(-1, 2, None, directIdcs=True )
    rigidPlates.append(rightRigidPlate)
    govNodes.append(np.array([ maxLim[0]/2+0.01, 0.011]))
    govNodesMechBC.append(utilitiesMech.mechanicalBC(dim, -2, rightRigidPlateMechBC))

    return node_coords, mechBC_merged, mechInitC_merged,  govNodes, govNodesMechBC, rigidPlates, functions

def assemble2d_box_with_periodic_nodes(maxLim, minDist, trials):
    dim = 2
    node_coords = []
    mechBC_merged = []
    mechInitC_merged = []
    govNodes = []
    govNodesMechBC = []
    rigidPlates = []
    functions = []

    indent = 1e-8

    oldlen = len(node_coords)
    nodeA = np.array([minDist*0.4, indent])
    nodeB = np.array([maxLim[0]-minDist*0.4, indent])
    pointGenerators.generateNodesLine2dRand(nodeA, nodeB, minDist*0.5, dim, node_coords, trials*20, catchCorners=True, equidist=False)
    newnodes = np.copy(node_coords[oldlen:])
    for k in newnodes:
        k[1] = maxLim[1]-indent
        node_coords.append(k)
    nodeA = np.array([indent, minDist*0.4])
    nodeB = np.array([indent, maxLim[1]-minDist*0.4])
    oldlen = len(node_coords)
    pointGenerators.generateNodesLine2dRand(nodeA, nodeB, minDist*0.5, dim, node_coords, trials*20, catchCorners=True, equidist=False)
    newnodes = np.copy(node_coords[oldlen:])
    for k in newnodes:
        k[0] = maxLim[0]-indent
        node_coords.append(k)

    node_coords = np.array(node_coords)
    radii = np.zeros(len(node_coords))+minDist*0.4

    print("MAX LIM", maxLim)
    #pointGenerators.generateNodesRect(maxLim, minDist, dim, trials, node_coords)
    node_coords, radii = pointGenerators.generateParticlesRect(maxLim, minDist*0.4, minDist, 0.8, 2, trials, node_coords, np.array(radii), allow_domain_overlap = False, periodic_distance=False)

    leftRigidPlateMechBC = np.array([0, 0,0,   -1,-1,-1])
    leftRigidPlate = utilitiesMech.RigidPlate(-1, 2, None, directIdcs=True )
    rigidPlates.append(leftRigidPlate)
    govNodes.append(np.array([ maxLim[0]/2-0.01, 0.011]))
    govNodesMechBC.append(utilitiesMech.mechanicalBC(dim, -1, leftRigidPlateMechBC))

    rightRigidPlateMechBC = np.array([1, 0,0,   -1,-1,-1])
    rightRigidPlate = utilitiesMech.RigidPlate(-1, 2, None, directIdcs=True )
    rigidPlates.append(rightRigidPlate)
    govNodes.append(np.array([ maxLim[0]/2+0.01, 0.011]))
    govNodesMechBC.append(utilitiesMech.mechanicalBC(dim, -2, rightRigidPlateMechBC))

    return node_coords, mechBC_merged, mechInitC_merged,  govNodes, govNodesMechBC, rigidPlates, functions, radii




def assemble3d_Hanging_FracZone(maxLim, minDist, topsupwidth, trials):
    dim = 3
    node_coords = []
    mechBC_merged = []
    mechInitC_merged = []
    govNodes = []
    govNodesMechBC = []
    rigidPlates = []
    functions = []

    indent = 1e-8

    print('bottom bound surface for fem frac zone')
    nodeA = np.array([indent, indent, indent])
    nodeB = np.array([maxLim[0]-indent, indent, maxLim[2]-indent])
    pointGenerators.generateNodesOrtoSurface3dRand(nodeA, nodeB, minDist*0.9, dim, node_coords, trials*2)
    print('top bound surface for fem frac zone')
    nodeA = np.array([maxLim[0]/2. - topsupwidth/2. + indent, maxLim[1]-indent, indent])
    nodeB = np.array([maxLim[0]/2. + topsupwidth/2. - indent, maxLim[1]-indent, maxLim[2]-indent])
    pointGenerators.generateNodesOrtoSurface3dRand(nodeA, nodeB, minDist*0.9, dim, node_coords, trials*2)
    print('left bound surface for fem frac zone')
    nodeA = np.array([indent, indent, indent])
    nodeB = np.array([indent, maxLim[1]-indent, maxLim[2]-indent])
    pointGenerators.generateNodesOrtoSurface3dRand(nodeA, nodeB, minDist*0.9, dim, node_coords, trials*2)
    print('right bound surface for fem frac zone')
    nodeA = np.array([maxLim[0]-indent, indent, indent])
    nodeB = np.array([maxLim[0]-indent, maxLim[1]-indent, maxLim[2]-indent])
    pointGenerators.generateNodesOrtoSurface3dRand(nodeA, nodeB, minDist*0.9, dim, node_coords, trials*2)


    print('fine fracture zone volume')
    node_coords = np.array(node_coords)
    radii = np.zeros(len(node_coords))+minDist*0.4
    #pointGenerators.generateNodesRect(maxLim, minDist, dim, trials, node_coords)
    node_coords, radii = pointGenerators.generateParticlesRect(maxLim, minDist*0.4, minDist, 0.8, 3, trials, node_coords, np.array(radii), allow_domain_overlap = False, periodic_distance=False)

    leftRigidPlateMechBC = np.array([0, 0,0,   -1,-1,-1, -1,-1,-1, -1,-1,-1])
    leftRigidPlate = utilitiesMech.RigidPlate(-1, 3, None, directIdcs=True )
    rigidPlates.append(leftRigidPlate)
    govNodes.append(np.array([ maxLim[0]/2-0.01, 0.011 ,maxLim[2]/2]))
    govNodesMechBC.append(utilitiesMech.mechanicalBC(dim, -1, leftRigidPlateMechBC))

    rightRigidPlateMechBC = np.array([1, 0,0,   -1,-1,-1, -1,-1,-1, -1,-1,-1])
    rightRigidPlate = utilitiesMech.RigidPlate(-1, 3, None, directIdcs=True )
    rigidPlates.append(rightRigidPlate)
    govNodes.append(np.array([ maxLim[0]/2+0.01, 0.011,maxLim[2]/2 ]))
    govNodesMechBC.append(utilitiesMech.mechanicalBC(dim, -2, rightRigidPlateMechBC))

    return node_coords, mechBC_merged, mechInitC_merged,  govNodes, govNodesMechBC, rigidPlates, functions, radii



def assemble3d_CFRAC_TDCB(maxLim, minDist, trials, holeMinDist, holeDiameter, roughMinDistCoef,elazonewidth, notchWidth, fracZoneOverhang=0.1, fracZoneHeight=0.8, fem=False):
    dim = 3
    node_coords = []
    mechBC_merged = []
    mechInitC_merged = []
    govNodes = []
    govNodesMechBC = []
    rigidPlates = []
    functions = []


    rougherCoef = 1.2
    roughCoef = 0.004/minDist#roughMinDistCoef

    notch = 0.37
    roughtop = (fracZoneHeight)
    fracW = elazonewidth

    #an indent due to mirroring of the data for voronoi tess.
    indent = 1e-8
    indentRP = 1e-7

    #sampling holes
    idcs_leftHole = []
    idcs_rightHole = []
    if fem == False:
        for r in range (2):
            print ('Hole #%d' %r)
            if r == 0:
                centre = np.array([ maxLim[0]/2-0.01, 0.011, indent ])
            if r == 1:
                centre = np.array([ maxLim[0]/2+0.01, 0.011, indent ])

            circleLength = 2*np.pi*holeDiameter/2
            nrNodes = int ( circleLength / holeMinDist )
            #print ('nrnodes: %d' %nrNodes)
            #deleno ctyrma
            nrNodes = (int (nrNodes / 4) +1 ) * 4

            nodesOld = len(node_coords)
            pointGenerators.generateNodesOrtoCilinderSurf3dRand(centre, holeDiameter/2,  maxLim[2]-indent, 2 ,       holeMinDist, node_coords, trials)
            newNodes = len(node_coords) - nodesOld
            if r == 0:
                idcs_leftHole = np.arange(0,newNodes)+nodesOld
                idcs_leftHole = idcs_leftHole.tolist()
                for h in idcs_leftHole:
                    if (node_coords[h][0] > (centre[0]-0.7071*holeDiameter/2)) :
                        idcs_leftHole.remove(h)

            if r == 1:
                idcs_rightHole = np.arange(0,newNodes)+nodesOld
                idcs_rightHole = idcs_rightHole.tolist()
                for h in idcs_rightHole:
                    if (node_coords[h][0] < (centre[0]+0.7071*holeDiameter/2)) :
                        idcs_rightHole.remove(h)



    leftRigidPlateMechBC = np.array([0, 0,0,   -1,-1,-1, -1,-1,-1, -1,-1,-1])
    leftRigidPlate = utilitiesMech.RigidPlate(-1, 3, None, directIdcs=True )
    leftRigidPlate.setDirectNodes(idcs_leftHole)
    rigidPlates.append(leftRigidPlate)
    govNodes.append(np.array([ maxLim[0]/2-0.01, 0.011 ,maxLim[2]/2]))
    govNodesMechBC.append(utilitiesMech.mechanicalBC(dim, -1, leftRigidPlateMechBC))

    rightRigidPlateMechBC = np.array([1, 0,0,   -1,-1,-1, -1,-1,-1, -1,-1,-1])
    rightRigidPlate = utilitiesMech.RigidPlate(-1, 3, None, directIdcs=True )
    rightRigidPlate.setDirectNodes(idcs_rightHole)
    rigidPlates.append(rightRigidPlate)
    govNodes.append(np.array([ maxLim[0]/2+0.01, 0.011,maxLim[2]/2 ]))
    govNodesMechBC.append(utilitiesMech.mechanicalBC(dim, -2, rightRigidPlateMechBC))


    if fem == True:
        print('top bound surface for fem frac zone')
        nodeA = np.array([maxLim[0]/2-fracW, maxLim[1]*roughtop, indent])
        nodeB = np.array([maxLim[0]/2+fracW, maxLim[1]*roughtop, maxLim[2]-indent])
        pointGenerators.generateNodesOrtoSurface3dRand(nodeA, nodeB, minDist*0.9, dim, node_coords, trials*2)
        print('left bound surface for fem frac zone')
        nodeA = np.array([maxLim[0]/2-fracW, maxLim[1]*(notch-notch*fracZoneOverhang), indent])
        nodeB = np.array([maxLim[0]/2-fracW, maxLim[1]*roughtop, maxLim[2]-indent])
        pointGenerators.generateNodesOrtoSurface3dRand(nodeA, nodeB, minDist*0.9, dim, node_coords, trials*2)
        print('right bound surface for fem frac zone')
        nodeA = np.array([maxLim[0]/2+fracW, maxLim[1]*(notch-notch*fracZoneOverhang), indent])
        nodeB = np.array([maxLim[0]/2+fracW, maxLim[1]*roughtop, maxLim[2]-indent])
        pointGenerators.generateNodesOrtoSurface3dRand(nodeA, nodeB, minDist*0.9, dim, node_coords, trials*2)



    notches = []
    #generating notch points
    if (notch > 0 and fracZoneOverhang>0):

        nTop = maxLim[1]*notch+indent-minDist/2
        if nTop < minDist:
            nTop = maxLim[1]*notch+indent

        notchMinDist= np.amin([notchWidth*0.3,maxLim[2]*0.2])

        nodeA = np.array([maxLim[0]/2, nTop+minDist/2, indent])
        nodeB = np.array([maxLim[0]/2, nTop+minDist/2, maxLim[2]-indent])
        pointGenerators.generateNodesLine3dRand(nodeA, nodeB, minDist/2, dim, node_coords, trials, catchCorners=True, equidist=True)

        #if notchMinDist < 0.001:
        #    notchMinDist = 0.001
        notchSide0 = []
        oldLen = len(node_coords)

        if fem == False:
            nodeA = np.array([maxLim[0]/2-notchWidth/2, indent, indent])
            nodeB = np.array([maxLim[0]/2-notchWidth/2, indent, maxLim[2]-indent])
            #pointGenerators.generateNodesLine3dRand(nodeA, nodeB, notchMinDist, dim, node_coords, trials*2, catchCorners=True, equidist=True)

        nodeA = np.array([maxLim[0]/2-notchWidth/2, nTop, indent])
        nodeB = np.array([maxLim[0]/2-notchWidth/2, nTop, maxLim[2]-indent])
        pointGenerators.generateNodesLine3dRand(nodeA, nodeB, notchMinDist, dim, node_coords, trials*2, catchCorners=True, equidist=True)

        nodeA = np.array([maxLim[0]/2-notchWidth/2, indent, indent])
        nodeB = np.array([maxLim[0]/2-notchWidth/2, maxLim[1]*notch-minDist, maxLim[2]-indent])
        if fem == False:
            pointGenerators.generateNodesOrtoSurface3dRand(nodeA, nodeB, notchMinDist, dim, node_coords, trials*2,minDistAmongNewPoints=True)
        nodeA = np.array([maxLim[0]/2-notchWidth/2, maxLim[1]*(notch-notch*fracZoneOverhang), indent])
        nodeB = np.array([maxLim[0]/2-notchWidth/2, maxLim[1]*notch-minDist, maxLim[2]-indent])
        pointGenerators.generateNodesOrtoSurface3dRand(nodeA, nodeB, minDist, dim, node_coords, trials*2)

        if fem==True:
            nodeA = np.array([maxLim[0]/2-fracW, maxLim[1]*(notch-notch*fracZoneOverhang), indent])
            nodeB = np.array([maxLim[0]/2-notchWidth/2, maxLim[1]*(notch-notch*fracZoneOverhang), maxLim[2]-indent])
            pointGenerators.generateNodesOrtoSurface3dRand(nodeA, nodeB, notchMinDist, dim, node_coords, trials*2)

        print('fracture zone overhang left')
        interBounds = np.array([     maxLim[0]/2-fracW,  maxLim[1]*notch , indent, maxLim[0]/2-notchWidth/2,   maxLim[1]*(notch-notch*fracZoneOverhang), maxLim[2]-indent])
        pointGenerators.generateNodesRect(interBounds, minDist, dim, trials, node_coords, useLowBound=True, topMinDist = minDist, bottomMinDist = minDist, gradienDirection=0)

        if fem == False:
            print('from bottom to notch left')
            interBounds = np.array([     indent,      indent , indent, maxLim[0]/2-notchWidth/2, maxLim[1]*notch*0.9, maxLim[2]-indent])
            pointGenerators.generateNodesRect(interBounds, minDist, dim, trials, node_coords, useLowBound=True, topMinDist = minDist*roughCoef/2, bottomMinDist = minDist*roughCoef/2, gradienDirection=1)

        for i in range (oldLen, len(node_coords), 1):
            notchSide0.append(i)


        notchSide1 = []
        oldLen = len(node_coords)

        if fem == False:
            nodeA = np.array([maxLim[0]/2+notchWidth/2, indent, indent])
            nodeB = np.array([maxLim[0]/2+notchWidth/2, indent, maxLim[2]-indent])
            #pointGenerators.generateNodesLine3dRand(nodeA, nodeB, notchMinDist, dim, node_coords, trials*2, catchCorners=True, equidist=True)

        nodeA = np.array([maxLim[0]/2+notchWidth/2, nTop, indent])
        nodeB = np.array([maxLim[0]/2+notchWidth/2, nTop, maxLim[2]-indent])
        pointGenerators.generateNodesLine3dRand(nodeA, nodeB, notchMinDist, dim, node_coords, trials*2, catchCorners=True, equidist=True)

        if fem == False:
            nodeA = np.array([maxLim[0]/2+notchWidth/2, indent, indent])
            nodeB = np.array([maxLim[0]/2+notchWidth/2, maxLim[1]*notch-minDist, maxLim[2]-indent])
            pointGenerators.generateNodesOrtoSurface3dRand(nodeA, nodeB, notchMinDist, dim, node_coords, trials*2,minDistAmongNewPoints=True)

        nodeA = np.array([maxLim[0]/2+notchWidth/2, maxLim[1]*(notch-notch*fracZoneOverhang), indent])
        nodeB = np.array([maxLim[0]/2+notchWidth/2, maxLim[1]*notch-minDist, maxLim[2]-indent])
        pointGenerators.generateNodesOrtoSurface3dRand(nodeA, nodeB, minDist, dim, node_coords, trials*2)

        if fem == True:
            nodeA = np.array([maxLim[0]/2+notchWidth/2, maxLim[1]*(notch-notch*fracZoneOverhang), indent])
            nodeB = np.array([maxLim[0]/2+fracW, maxLim[1]*(notch-notch*fracZoneOverhang), maxLim[2]-indent])
            pointGenerators.generateNodesOrtoSurface3dRand(nodeA, nodeB, notchMinDist, dim, node_coords, trials*2)

        print('fracture zone overhang right')
        interBounds = np.array([     maxLim[0]/2+notchWidth/2,  maxLim[1]*notch , indent, maxLim[0]/2+fracW ,   maxLim[1]*(notch-notch*fracZoneOverhang), maxLim[2]-indent])
        pointGenerators.generateNodesRect(interBounds, minDist, dim, trials, node_coords, useLowBound=True, topMinDist = minDist, bottomMinDist = minDist, gradienDirection=0)

        if fem == False:
            print('from bottom to notch right')
            interBounds = np.array([     maxLim[0]/2+notchWidth/2,      indent , indent, maxLim[0], maxLim[1]*notch*0.9, maxLim[2]-indent])
            pointGenerators.generateNodesRect(interBounds, minDist, dim, trials, node_coords, useLowBound=True, topMinDist = minDist*roughCoef/2, bottomMinDist = minDist*roughCoef/2, gradienDirection=1)



        for i in range (oldLen, len(node_coords), 1):
            notchSide1.append(i)

        notchA = []
        notchA.append(notchSide0)
        notchA.append(notchSide1)
        notches.append(notchA)

    else:
       print('bottom bound surface for fem frac zone')
       nodeA = np.array([maxLim[0]/2-fracW, maxLim[1]*(notch-notch*fracZoneOverhang), indent])
       nodeB = np.array([maxLim[0]/2+fracW, maxLim[1]*(notch-notch*fracZoneOverhang), maxLim[2]-indent])
       pointGenerators.generateNodesOrtoSurface3dRand(nodeA, nodeB, minDist*0.9, dim, node_coords, trials*2)

    print('Generating volumes...')


    #left
    if fem == False:
        print('left gradient')
        oldnodes=len(node_coords)
        interBounds = np.array([     indent,      maxLim[1]*notch*0.9 , indent, maxLim[0]/2,   maxLim[1], maxLim[2]-indent])
        pointGenerators.generateNodesRect(interBounds, minDist, dim, trials, node_coords, useLowBound=True, topMinDist = minDist*roughCoef, bottomMinDist = minDist*roughCoef/3, gradienDirection=0)

    #right
    if fem == False:
        print('right gradient')
        interBounds = np.array([     maxLim[0]/2,      maxLim[1]*notch*0.9 , indent, maxLim[0],   maxLim[1], maxLim[2]-indent])
        pointGenerators.generateNodesRect(interBounds, minDist, dim, trials, node_coords, useLowBound=True, topMinDist = minDist*roughCoef/3, bottomMinDist = minDist*roughCoef, gradienDirection=0)

    if fem == False:
        print('from top to crack')
        interBounds = np.array([     maxLim[0]/2-fracW,      maxLim[1]*roughtop , indent, maxLim[0]/2+fracW,   maxLim[1], maxLim[2]-indent])
        pointGenerators.generateNodesRect(interBounds, minDist, dim, trials, node_coords, useLowBound=True, topMinDist = minDist, bottomMinDist = minDist*roughCoef/2, gradienDirection=1)

    #center
    print('fine fracture zone')
    interBounds = np.array([     maxLim[0]/2-fracW,  maxLim[1]*notch , indent, maxLim[0]/2+fracW,   maxLim[1]*roughtop, maxLim[2]-indent])
    pointGenerators.generateNodesRect(interBounds, minDist, dim, trials, node_coords, useLowBound=True, topMinDist = minDist, bottomMinDist = minDist, gradienDirection=0)


    #remove points from rebars
    newNodes = []
    for i in range (len(node_coords)):
        app = True
        for r in range(2):
            if r == 0:
                centre = np.array([ maxLim[0]/2-0.01, 0.011])
                if (np.linalg.norm(centre - node_coords[i][0:2]) < holeDiameter/2*0.999):
                    app = False
            if r == 1 and app==True:
                centre = np.array([ maxLim[0]/2+0.01, 0.011])
                if not (np.linalg.norm(centre - node_coords[i][0:2]) > holeDiameter/2*0.999):
                    app = False

        if app==True:
            newNodes.append(node_coords[i])

    node_coords = newNodes.copy()



    #identifying nodes within modelu
    model_indices = []
    for i,p in enumerate(node_coords):
        app = True
        if p[0]<0.015:
            if not (p[1]>=(0.1-(6+2/3)*p[0]) ):
                app = False
        if p[0]>0.075:
            if not (p[1]>=(6+2/3)*(p[0]-0.075) ):
                app = False

        if app == True:
            model_indices.append(i)


    print('model points %s' %len(model_indices))
    node_coords = np.asarray(node_coords)
    if False:
        plt.plot(node_coords[model_indices,0], node_coords[model_indices,1], 'o', color='black');
        plt.show()

        plt.plot(node_coords[model_indices,2], node_coords[model_indices,1], 'o', color='black');
        plt.show()

    print('Voronoi...')
    return node_coords, mechBC_merged, mechInitC_merged,  govNodes, govNodesMechBC, rigidPlates, functions, notches,model_indices








##############################x
def assemble2dCorrosionRebar(maxLim, minDist, trials, rebarMinDist, interfaceMinDist, rebarDiameter, rebarCount, rebarDepth, sampleBorders, node_coords_init=None, roughMinDistCoef=1, adaptivityReady=False, fineRingThickness=1, fineRegDepth=0.5, gradientRegDepth=0.2):
    dim = 2

    print (roughMinDistCoef)

    if node_coords_init is None:
        node_coords = []
    else:
        node_coords = node_coords_init


    mechBC_merged = []
    mechInitC_merged = []
    govNodes = []
    govNodesMechBC = []
    rigidPlates = []
    expansionRings = []
    functions = []

    #an indent due to mirroring of the data for voronoi tess.
    indent = 1e-8
    indentRP = 1e-7


    if node_coords_init is None:
        node_coords.append(np.array([indent, 2*indentRP]))
    mechBC = np.array([0,0,-1,-1,-1,-1])
    mBC = utilitiesMech.mechanicalBC(dim, 0, mechBC)
    mechBC_merged.append(mBC)

    if node_coords_init is None:
        #protože když je uplně v koutě, tak casto neni propojeny se dvema elementy a nema tuhost v rotaci
        node_coords.append(np.array([maxLim[0]*0.9, 2*indentRP]))
    mechBC = np.array([-1,0,-1,-1,-1,-1])
    mBC = utilitiesMech.mechanicalBC(dim, 1, mechBC)
    mechBC_merged.append(mBC)

    interfaceNodeIndices = []

    #sampling interfaces
    for r in range (rebarCount):
        if node_coords_init is None:
            oldLen = len(node_coords)
            interface = []

            print ('Interface #%d' %r)
            #rebar edge
            if (rebarCount==1):
                #puvodni poloha rebars polovina od kraje
                centre = np.array([ (maxLim[0]/rebarCount)*(r+0.5), maxLim[1]-rebarDepth  ])
            else:
                #poloha rebars presne jak je ve clanku
                centre = np.array([ (0.058 + (maxLim[0]-0.116)/(rebarCount-1)*r), maxLim[1]-rebarDepth  ])
            #print(centre)

            circleLength = 2*np.pi*rebarDiameter/2
            nrNodes = int ( circleLength / interfaceMinDist )
            #print ('nrnodes: %d' %nrNodes)
            #deleno ctyrma
            nrNodes = (int (nrNodes / 4) +1 ) * 4

            nodesOld = len(node_coords)
            pointGenerators.generateNodesCircle2dRand(centre, rebarDiameter/2, interfaceMinDist, node_coords, trials, equiAngNodes=nrNodes )
            newNodes = len(node_coords) - nodesOld

            for n in range ( newNodes ):
                #print(nodesOld + n)
                lnfc = len(functions)
                mechBC = np.array([lnfc, lnfc+1, -1, -1 , -1, -1,])
                point = node_coords[nodesOld + n]

                directionVec = (point - centre) / np.linalg.norm(point - centre)
                uCoef = 0.01
                directionVec *= uCoef

                funcX = []
                funcX.append( np.array([0,0]) )
                funcX.append( np.array([1, directionVec[0] ]) )
                fnX = utilitiesNumeric.generalFunc(funcX)
                #functions.append(fnX)
                #print(funcX)

                funcY = []
                funcY.append( np.array([0,0]) )
                funcY.append( np.array([1, directionVec[1] ]) )
                fnY = utilitiesNumeric.generalFunc(funcY)
                #functions.append(fnY)


                mBC = utilitiesMech.mechanicalBC(dim, nodesOld + n, mechBC)
                #mechBC_merged.append(mBC)

            newLen = len(node_coords)

            for i in range (newLen-oldLen):
                interface.append(oldLen+i)
            interfaceNodeIndices.append(interface)

        rebarBC = np.array([-1, 2,  -1,-1,-1,-1,  -1,-1,-1])
        govNodesMechBC.append(utilitiesMech.mechanicalBC(dim, (-3), rebarBC))



    if node_coords_init is None:
        if sampleBorders:

            #top
            nodeA = np.array([indent, maxLim[1]-indent])
            nodeB = np.array([maxLim[0]-indent, maxLim[1]-indent])
            if adaptivityReady:
                pointGenerators.generateNodesLine2dRand(nodeA, nodeB, minDist*0.8, dim, node_coords, trials*2, catchCorners=True, equidist=False)
            else:
                pointGenerators.generateNodesLine2dRand(nodeA, nodeB, minDist*roughMinDistCoef*0.7, dim, node_coords, trials*2, catchCorners=True, equidist=False)
            """
            topRigidPlateMechBC = np.array([-1, -1,-1,   -1,-1,-1])
            topRigidPlate = utilitiesMech.RigidPlate(-1, 2, np.array([-indentRP, maxLim[0]+indentRP, maxLim[1]-indentRP, maxLim[1]+indentRP ]))
            rigidPlates.append(topRigidPlate)
            govNodes.append(np.array([ maxLim[0]/2, maxLim[1]-indent ]))
            govNodesMechBC.append(utilitiesMech.mechanicalBC(dim, -1, topRigidPlateMechBC))
            """

            #bottom
            nodeA = np.array([indent, indent])
            nodeB = np.array([maxLim[0]-indent, indent])
            if adaptivityReady:
                pointGenerators.generateNodesLine2dRand(nodeA, nodeB, minDist*roughMinDistCoef*0.8, dim, node_coords, trials, catchCorners=False, equidist=False)
            else:
                pointGenerators.generateNodesLine2dRand(nodeA, nodeB, minDist*roughMinDistCoef*0.8, dim, node_coords, trials, catchCorners=False, equidist=False)

            """
            bottomRigidPlateMechBC = np.array([-1, 0,0,   -1,-1,-1])
            bottomRigidPlate = utilitiesMech.RigidPlate(-1, 2, np.array([-indentRP, maxLim[0]+indentRP, -indentRP, indentRP ]))
            rigidPlates.append(bottomRigidPlate)
            govNodes.append(np.array([ maxLim[0]/2, indentRP ]))
            govNodesMechBC.append(utilitiesMech.mechanicalBC(dim, -2, bottomRigidPlateMechBC))
            """



            #left
            nodeA = np.array([indent, indent])
            nodeB = np.array([indent, maxLim[1]-indent])
            #pointGenerators.generateNodesLine2dRand(nodeA, nodeB, minDist*0.8, dim, node_coords, trials, catchCorners=False, equidist=False)

            """
            leftRigidPlateMechBC = np.array([0, -1,-1,   -1,-1,-1])
            leftRigidPlate = utilitiesMech.RigidPlate(-1, 2, np.array([-indentRP, indentRP, minDist/2, maxLim[1]-minDist/2 ]))
            rigidPlates.append(leftRigidPlate)
            govNodes.append(np.array([ indentRP, maxLim[1]/2 ]))
            govNodesMechBC.append(utilitiesMech.mechanicalBC(dim, -3, leftRigidPlateMechBC))
            """

            #right
            nodeA = np.array([maxLim[0]-indent, indent])
            nodeB = np.array([maxLim[0]-indent, maxLim[1]-indent])
            #pointGenerators.generateNodesLine2dRand(nodeA, nodeB, minDist*0.8, dim, node_coords, trials, catchCorners=False, equidist=False)

            """
            rightRigidPlateMechBC = np.array([0, -1,-1,   -1,-1,-1])
            rightRigidPlate = utilitiesMech.RigidPlate(-1, 2, np.array([maxLim[0]-indentRP, maxLim[0]+indentRP, minDist/2, maxLim[1]-minDist/2 ]))
            rigidPlates.append(rightRigidPlate)
            govNodes.append(np.array([ maxLim[0], maxLim[1]/2 ]))
            govNodesMechBC.append(utilitiesMech.mechanicalBC(dim, -4, rightRigidPlateMechBC))
            """






    if node_coords_init is None:
        fineRegDepth *= maxLim[1]
        #top half rect
        fineTopBounds = np.array([     indent,          maxLim[1] -fineRegDepth,      maxLim[0],         maxLim[1]   ])

        if adaptivityReady:
            pointGenerators.generateNodesRect(fineTopBounds, minDist, dim, trials, node_coords, useLowBound=True)
        else:
            pointGenerators.generateNodesRect(fineTopBounds, minDist, dim, trials, node_coords, useLowBound=True)
        #pointGenerators.generateNodesRect(maxLim, minDist*roughMinDistCoef, dim, trials, node_coords, useLowBound=True)

        #sampling rebars
        for r in range (rebarCount):
            #rebar edge
            if (rebarCount==1):
                #puvodni poloha rebars polovina od kraje
                centre = np.array([ (maxLim[0]/rebarCount)*(r+0.5), maxLim[1]-rebarDepth  ])
            else:
                #poloha rebars presne jak je ve clanku
                centre = np.array([ (0.058 + (maxLim[0]-0.116)/(rebarCount-1)*r), maxLim[1]-rebarDepth  ])

            #fine surf above rebar
            nodeA = np.array([centre[0]-rebarDiameter, maxLim[1]-indent ])
            nodeB = np.array([centre[0]+rebarDiameter, maxLim[1]-indent ])
            #pointGenerators.generateNodesLine2dRand(nodeA, nodeB, minDist, dim, node_coords, trials, catchCorners=False, equidist=False)

            if adaptivityReady:


                #sampling fine annulus around rebar
                #pointGenerators.generateNodesOrtoAnnulus2dRand(centre, (rebarDiameter/2)*(fineRingThickness+1), (rebarDiameter/2)*(fineRingThickness),  minDist, node_coords, trials)
                #sampling gradient annulus around rebar
                pointGenerators.generateNodesOrtoAnnulus2dRand(centre,  (rebarDiameter/2)*(fineRingThickness+1), (rebarDiameter/2)*(fineRingThickness),  minDist/roughMinDistCoef, node_coords, trials, minD=minDist/roughMinDistCoef, maxD=minDist)

            #bounds = np.array([     centre[0]-rebarDiameter*2,          centre[1],      centre[0]+rebarDiameter*2,       maxLim[1]-indent     ])
            #pointGenerators.generateNodesRect(bounds, minDist*roughMinDistCoef/2, dim, trials, node_coords, useLowBound=True)

            """
            if adaptivityReady:
                #left rect
                interHeight = (maxLim[1]-fineRegDepth)  / 2
                bounds = np.array([     centre[0]-rebarDiameter/2*4,          centre[1],      centre[0]-rebarDiameter,       maxLim[1]-indent     ])
                pointGenerators.generateNodesRect(bounds, minDist, dim, trials, node_coords, useLowBound=True, topMinDist = minDist*roughMinDistCoef/2, bottomMinDist = minDist, gradienDirection=0)

                #right rect
                interHeight = (maxLim[1]-fineRegDepth)  / 2
                bounds = np.array([     centre[0]+rebarDiameter,          centre[1],      centre[0]+rebarDiameter/2*4,       maxLim[1]-indent     ])
                pointGenerators.generateNodesRect(bounds, minDist, dim, trials, node_coords, useLowBound=True, topMinDist = minDist, bottomMinDist = minDist*roughMinDistCoef/2, gradienDirection=0)


                bounds = np.array([     centre[0]-rebarDiameter,          centre[1],      centre[0]+rebarDiameter,       maxLim[1]-indent     ])
                pointGenerators.generateNodesRect(bounds, minDist, dim, trials, node_coords, useLowBound=True)
            """



            #remove points from rebars
            newNodes = []
            for i in range (len(node_coords)):
                if (np.linalg.norm(centre - node_coords[i]) > rebarDiameter/2*0.999):
                    newNodes.append(node_coords[i])

            node_coords = newNodes.copy()

            if (rebarMinDist<0):
                print ('rebars centre')
            #    node_coords.append(centre*1e-5)
            else:
                #rebar crossection
                pointGenerators.generateNodesOrtoCircle2dRand(centre, rebarDiameter/2, rebarMinDist, node_coords, trials)







        #intermediate rect
        interHeight = maxLim[1] * gradientRegDepth
        interBounds = np.array([     indent,      maxLim[1] -fineRegDepth - interHeight , maxLim[0],              maxLim[1] -fineRegDepth])
        #interBounds = np.array([     indent,     indent, maxLim[0],              maxLim[1] -fineRegDepth])
        if node_coords_init is None:
            if adaptivityReady:
                pointGenerators.generateNodesRect(interBounds, minDist, dim, trials, node_coords, useLowBound=True, topMinDist = minDist, bottomMinDist = minDist, gradienDirection=1)
            else:
                pointGenerators.generateNodesRect(interBounds, minDist, dim, trials, node_coords, useLowBound=True, topMinDist = minDist*roughMinDistCoef, bottomMinDist = minDist, gradienDirection=1)


            #bottom rough rect
            roughBottomBounds =  np.array([     indent,      indent,  maxLim[0],              maxLim[1] -fineRegDepth - interHeight ])
            if adaptivityReady:
                pointGenerators.generateNodesRect(roughBottomBounds, minDist, dim, trials, node_coords, useLowBound=True)
            else:
                pointGenerators.generateNodesRect(roughBottomBounds, minDist*roughMinDistCoef, dim, trials, node_coords, useLowBound=True)




    node_coords = np.asarray(node_coords)
    if SHOW_PLOT:
        plt.plot(node_coords[:,0], node_coords[:,1], 'o', color='black');
        plt.show()


    return node_coords, mechBC_merged, mechInitC_merged,  govNodes, govNodesMechBC, rigidPlates, functions, interfaceNodeIndices

###

def assemble3dCorrosionRebar(maxLim, minDist, trials, rebarMinDist, interfaceMinDist, rebarDiameter, rebarCount, rebarDepth, sampleBorders, node_coords_init=None, roughMinDistCoef=1, adaptivityReady=False, fineRingThickness=1, fineRegDepth=0.5, gradientRegDepth=0.2):
    dim = 3

    rebarRadius = rebarDiameter/2

    if node_coords_init is None:
        node_coords = []
    else:
        node_coords = node_coords_init


    mechBC_merged = []
    mechInitC_merged = []
    govNodes = []
    govNodesMechBC = []
    rigidPlates = []
    expansionRings = []
    functions = []

    #an indent due to mirroring of the data for voronoi tess.
    indent = 1e-8
    indentRP = 1e-6


    if node_coords_init is None:
        node_coords.append(np.array([indentRP*3, indentRP*3, indentRP*3]))
        mechBC = np.array([0,0,0,  0,0,0,  -1,-1,-1,   -1,-1,-1])
        mBC = utilitiesMech.mechanicalBC(dim, 0, mechBC)
        mechBC_merged.append(mBC)


    nodeA = np.array([indent, indent, indent])
    nodeB = np.array([maxLim[0]-indent, indent, maxLim[2]-indent])
    oldLen = len(node_coords)
    pointGenerators.generateNodesOrtoSurface3dRand(nodeA, nodeB, minDist, dim, node_coords, trials)

    ##################### CONSTRAINTS AND RIGID PLATES
    indentRP = 1e-6
    topRigidPlateMechBC = np.array([-1, -1,-1,  0, 0, 0, -1,-1,-1, -1,-1,-1])#np.array([1,0,-1])
    topRigidPlate = utilitiesMech.RigidPlate(-1, 3,    np.array([    indentRP, maxLim[0]-indentRP, -indentRP, indentRP, -indentRP,      maxLim[2]+indentRP])    )
    rigidPlates.append(topRigidPlate)
    govNodes.append(np.array([ maxLim[0]/2, indent, maxLim[2]/2 ]))
    govNodesMechBC.append(utilitiesMech.mechanicalBC(dim, -1, topRigidPlateMechBC))

    #front surf
    frontRigidPlateMechBC = np.array([0, -1,-1,  0, 0, 0, -1,-1,-1, -1,-1,-1])#np.array([1,0,-1])
    frontRigidPlate = utilitiesMech.RigidPlate(-1, 3,    np.array([    -indentRP, maxLim[0]+indentRP, -indentRP, maxLim[1]+indentRP, -indentRP,      indentRP])    )
    rigidPlates.append(frontRigidPlate)
    govNodes.append(np.array([ maxLim[0]/2, maxLim[1]/2, 0 ]))
    govNodesMechBC.append(utilitiesMech.mechanicalBC(dim, -2, frontRigidPlateMechBC))

    nodeA = np.array([indent, indent, indent])
    nodeB = np.array([maxLim[0]-indent, maxLim[1]-indent, indent])
    oldLen = len(node_coords)
    pointGenerators.generateNodesOrtoSurface3dRand(nodeA, nodeB, minDist, dim, node_coords, trials)


    #rear surf
    frontRigidPlateMechBC = np.array([0, -1,-1,  0, 0, 0, -1,-1,-1, -1,-1,-1])#np.array([1,0,-1])
    frontRigidPlate = utilitiesMech.RigidPlate(-1, 3,    np.array([    -indentRP, maxLim[0]+indentRP, -indentRP, maxLim[1]+indentRP, maxLim[2]-indentRP,     maxLim[2]+indentRP])    )
    rigidPlates.append(frontRigidPlate)
    govNodes.append(np.array([ maxLim[0]/2, maxLim[1]/2, maxLim[2] ]))
    govNodesMechBC.append(utilitiesMech.mechanicalBC(dim, -3, frontRigidPlateMechBC))

    nodeA = np.array([indent, indent, maxLim[2]-indent])
    nodeB = np.array([maxLim[0]-indent, maxLim[1]-indent, maxLim[2]-indent])
    oldLen = len(node_coords)
    pointGenerators.generateNodesOrtoSurface3dRand(nodeA, nodeB, minDist, dim, node_coords, trials)


    """
    node_coords=np.asarray(node_coords)
    fig = plt.figure()
    ax = Axes3D(fig)
    ax.scatter(node_coords[:,0], node_coords[:,1], node_coords[:,2])
    plt.show()
    """

    nodeA = np.array([indent, maxLim[1]-indent, indent])
    nodeB = np.array([maxLim[0]-indent, maxLim[1]-indent, maxLim[2]-indent])
    oldLen = len(node_coords)
    pointGenerators.generateNodesOrtoSurface3dRand(nodeA, nodeB, minDist, dim, node_coords, trials)




    if node_coords_init is None:
        if sampleBorders:
            #top
            nodeA = np.array([indent, maxLim[1]-indent,  indent])
            nodeB = np.array([maxLim[0]-indent, maxLim[1]-indent, maxLim[2]-indent])
            if adaptivityReady:
                #pointGenerators.generateNodesLine2dRand(nodeA, nodeB, minDist*0.7, dim, node_coords, trials*2, catchCorners=True, equidist=False)
                pointGenerators.generateNodesOrtoSurface3dRand(nodeA, nodeB, minDist*0.7, dim, node_coords, trials*2)
            else:
                #pointGenerators.generateNodesLine2dRand(nodeA, nodeB, minDist*roughMinDistCoef*0.7, dim, node_coords, trials*2, catchCorners=True, equidist=False)
                pointGenerators.generateNodesOrtoSurface3dRand(nodeA, nodeB, minDist*roughMinDistCoef*0.7, dim, node_coords, trials*2)

            #bottom
            nodeA = np.array([indent, indent, indent])
            nodeB = np.array([maxLim[0]-indent, indent, maxLim[2]-indent])
            if adaptivityReady:
                pointGenerators.generateNodesOrtoSurface3dRand(nodeA, nodeB, minDist*0.7, dim, node_coords, trials*2)
            else:
                pointGenerators.generateNodesOrtoSurface3dRand(nodeA, nodeB, minDist*roughMinDistCoef*0.7, dim, node_coords, trials*2)


    #rebar

    interfaceNodeIndices = []

    #sampling interfaces
    for r in range (rebarCount):
        if node_coords_init is None:
            oldLen = len(node_coords)
            interface = []

            print ('Interface #%d' %r)
            #rebar edge
            if (rebarCount==1):
                #puvodni poloha rebars polovina od kraje
                centre = np.array([ (maxLim[0]/rebarCount)*(r+0.5), maxLim[1]-rebarDepth, indent  ])
            else:
                #poloha rebars presne jak je ve clanku
                centre = np.array([ (0.058 + (maxLim[0]-0.116)/(rebarCount-1)*r), maxLim[1]-rebarDepth , indent ])
            #print(centre)

            circleLength = 2*np.pi*rebarDiameter/2
            nrNodes = int ( circleLength / interfaceMinDist )
            #print ('nrnodes: %d' %nrNodes)
            #deleno ctyrma
            nrNodes = (int (nrNodes / 4) +1 ) * 4

            nodesOld = len(node_coords)
            #pointGenerators.generateNodesCircle2dRand(centre, rebarDiameter/2, interfaceMinDist, node_coords, trials, equiAngNodes=nrNodes )
            pointGenerators.generateNodesOrtoCilinderSurf3dRand(centre, rebarRadius-1e-5, maxLim[2]-indent-indent, 2, interfaceMinDist,  node_coords, trials, equiAngNodes=nrNodes )
            newNodes = len(node_coords) - nodesOld

            for n in range ( newNodes ):
                #print(nodesOld + n)
                lnfc = len(functions)
                mechBC = np.array([lnfc, lnfc+1, -1, -1 , -1, -1,])
                point = node_coords[nodesOld + n]

                directionVec = (point - centre) / np.linalg.norm(point - centre)
                uCoef = 0.01
                directionVec *= uCoef

                funcX = []
                funcX.append( np.array([0,0]) )
                funcX.append( np.array([1, directionVec[0] ]) )
                fnX = utilitiesNumeric.generalFunc(funcX)
                #functions.append(fnX)
                #print(funcX)

                funcY = []
                funcY.append( np.array([0,0]) )
                funcY.append( np.array([1, directionVec[1] ]) )
                fnY = utilitiesNumeric.generalFunc(funcY)
                #functions.append(fnY)


                mBC = utilitiesMech.mechanicalBC(dim, nodesOld + n, mechBC)
                #mechBC_merged.append(mBC)

        newLen = len(node_coords)

        for i in range (newLen-oldLen):
            interface.append(oldLen+i)
        interfaceNodeIndices.append(interface)

        rebarBC = np.array([2, -1,              -1,-1, -1, -1,  -1,-1,-1,   -1,-1,-1])
        govNodesMechBC.append(utilitiesMech.mechanicalBC(dim, (-3), rebarBC))



    if node_coords_init is None:
        fineRegDepth *= maxLim[1]
        #top half rect
        fineTopBounds = np.array([     indent,          maxLim[1] -fineRegDepth,    indent,   maxLim[0],         maxLim[1], maxLim[2]  ])

        if adaptivityReady:
            pointGenerators.generateNodesRect(fineTopBounds, minDist, dim, trials, node_coords, useLowBound=True)
        else:
            pointGenerators.generateNodesRect(fineTopBounds, minDist, dim, trials, node_coords, useLowBound=True)
        #pointGenerators.generateNodesRect(maxLim, minDist*roughMinDistCoef, dim, trials, node_coords, useLowBound=True)


        #sampling rebars
        for r in range (rebarCount):
            #rebar edge
            if (rebarCount==1):
                #puvodni poloha rebars polovina od kraje
                centre = np.array([ (maxLim[0]/rebarCount)*(r+0.5), maxLim[1]-rebarDepth, indent  ])
            else:
                #poloha rebars presne jak je ve clanku
                centre = np.array([ (0.058 + (maxLim[0]-0.116)/(rebarCount-1)*r), maxLim[1]-rebarDepth, indent  ])

            #fine surf above rebar
            nodeA = np.array([centre[0]-rebarDiameter, maxLim[1]-indent, indent ])
            nodeB = np.array([centre[0]+rebarDiameter, maxLim[1]-indent, maxLim[2]-indent ])
            #pointGenerators.generateNodesLine2dRand(nodeA, nodeB, minDist, dim, node_coords, trials, catchCorners=False, equidist=False)

            if adaptivityReady:
                pointGenerators.generateNodesOrtoAnnulus2dRand(centre,  (rebarDiameter/2)*(fineRingThickness+1), (rebarDiameter/2)*(fineRingThickness),  minDist/roughMinDistCoef, node_coords, trials, minD=minDist/roughMinDistCoef, maxD=minDist)


        #intermediate rect
        interHeight = maxLim[1] * gradientRegDepth
        interBounds = np.array([     indent,      maxLim[1] -fineRegDepth - interHeight , indent, maxLim[0],              maxLim[1] -fineRegDepth   ,  maxLim[2]-indent])
        #interBounds = np.array([     indent,     indent, maxLim[0],              maxLim[1] -fineRegDepth])
        if adaptivityReady:
            pointGenerators.generateNodesRect(interBounds, minDist, dim, trials, node_coords, useLowBound=True, topMinDist = minDist, bottomMinDist = minDist, gradienDirection=1)
        else:
            pointGenerators.generateNodesRect(interBounds, minDist, dim, trials, node_coords, useLowBound=True, topMinDist = minDist*roughMinDistCoef, bottomMinDist = minDist, gradienDirection=1)


        #bottom rough rect
        roughBottomBounds =  np.array([     indent,      indent, indent,  maxLim[0],              maxLim[1] -fineRegDepth - interHeight,  maxLim[2]-indent ])
        if adaptivityReady:
            pointGenerators.generateNodesRect(roughBottomBounds, minDist, dim, trials, node_coords, useLowBound=True)
        else:
            pointGenerators.generateNodesRect(roughBottomBounds, minDist*roughMinDistCoef, dim, trials, node_coords, useLowBound=True)


        #remove points from rebars
        newNodes = []
        for i in range (len(node_coords)):
            if (np.linalg.norm(centre[0:2] - node_coords[i][0:2]) > rebarDiameter/2*0.99):
                newNodes.append(node_coords[i])

        #.copy()
    node_coords = newNodes
    node_coords = np.asarray(node_coords)



    node_coords=np.asarray(node_coords)
    fig = plt.figure()
    ax = Axes3D(fig)
    ax.scatter(node_coords[:,0], node_coords[:,1], node_coords[:,2])
    plt.show()



    """
    fig, ax = plt.subplots()
    ax.scatter(node_coords[:,0], node_coords[:,1])
    plt.show()


    fig = plt.figure()
    ax = Axes3D(fig)
    ax.scatter(node_coords[:,0], node_coords[:,1], node_coords[:,2])
    plt.show()
    #"""

    return node_coords, mechBC_merged, mechInitC_merged,  govNodes, govNodesMechBC, rigidPlates, functions, interfaceNodeIndices




def assemble2dCoupledPress(maxLim, minDist, trials):
    dim = 2
    #lists for the model
    node_coords = []
    mechBC_merged = []
    mechInitC_merged = []
    govNodes = []
    govNodesMechBC = []
    rigidPlates = []

    #an indent due to mirroring of the data for voronoi tess.
    indent = 1e-8



    node_coords.append(np.array([maxLim[0]/2, indent]))
    node_coords.append(np.array([maxLim[0]/2, maxLim[1]-indent]))
    #node_coords.append(np.array([maxLim[0]/2-notchWidth, maxLim[1]*notch/2 , maxLim[2]/2]))

    ##################### CONSTRAINTS AND RIGID PLATES
    #top plate
    indentRP = 1e-6
    topRigidPlateMechBC = np.array([0, 1,0,   -1,-1,-1])#np.array([1,0,-1])
    topRigidPlate = utilitiesMech.RigidPlate(-1, 2, np.array([ -indentRP, maxLim[0]+indentRP, -indentRP, indentRP ]))
    rigidPlates.append(topRigidPlate)
    govNodes.append(np.array([ maxLim[0]/2, 0 ]))
    govNodesMechBC.append(utilitiesMech.mechanicalBC(dim, -1, topRigidPlateMechBC))

    #bottom plate
    bottomRigidPlateMechBC = np.array([0, 0,0,   -1,-1,-1]) #np.array([2,0,-1])
    bottomRigidPlate = utilitiesMech.RigidPlate(-1, 2, np.array([ -indentRP, maxLim[0]+indentRP,-indentRP+maxLim[1], maxLim[1]+indentRP ]))
    rigidPlates.append(bottomRigidPlate)
    govNodes.append(np.array([ maxLim[0]/2, maxLim[1] ]))
    govNodesMechBC.append(utilitiesMech.mechanicalBC(dim, -2, bottomRigidPlateMechBC))


    facesMult = 1
    ########################################## rest of  faces
    nodeA =  np.array([indent , indent])
    nodeB =  np.array([maxLim[0] - indent,  indent])
    pointGenerators.generateNodesLine2dRand(nodeA, nodeB, minDist*facesMult, dim, node_coords,  trials, False, False)


    nodeA =  np.array([indent , maxLim[1] - indent])
    nodeB =  np.array([maxLim[0] - indent,  maxLim[1] -indent])
    pointGenerators.generateNodesLine2dRand(nodeA, nodeB, minDist*facesMult, dim, node_coords,  trials, False, False)


    #rect
    pointGenerators.generateNodesRect(maxLim, minDist, dim, trials, node_coords)


    return node_coords, mechBC_merged, mechInitC_merged,  govNodes, govNodesMechBC, rigidPlates








def assemble2dCoupledArtificialCrack (maxLim, minDist, trials, slitWidth, notch):
    dim = 2
    #lists for the model
    node_coords = []
    mechBC_merged = []
    mechInitC_merged = []
    govNodes = []
    govNodesMechBC = []
    rigidPlates = []

    #an indent due to mirroring of the data for voronoi tess.
    indent = 1e-8
    notchWidth = slitWidth/2
    #generating notch points

    nodeA = np.array([maxLim[0]/2-notchWidth, indent])
    nodeB = np.array([maxLim[0]/2-notchWidth, maxLim[1]*notch])
    oldLen = len(node_coords)
    pointGenerators.generateNodesLine2dRand(nodeA, nodeB, minDist, dim, node_coords, trials, catchCorners=True, equidist=True)

    nodeA = np.array([maxLim[0]/2+notchWidth, indent])
    nodeB = np.array([maxLim[0]/2+notchWidth, maxLim[1]*notch])
    oldLen = len(node_coords)
    pointGenerators.generateNodesLine2dRand(nodeA, nodeB, minDist, dim, node_coords, trials, catchCorners=True, equidist=True)

    #node_coords.append(np.array([maxLim[0]/2+notchWidth, maxLim[1]*notch]))
    #node_coords.append(np.array([maxLim[0]/2-notchWidth, maxLim[1]*notch]))


    ##################### CONSTRAINTS AND RIGID PLATES
    #rigid plate left slit face
    indentRP = 1e-9
    leftRigidPlateMechBC = np.array([1, 0,-1,   -1,-1,-1])#np.array([1,0,-1])
    leftRigidPlate = utilitiesMech.RigidPlate(-1, 2, np.array([ maxLim[0]/2-notchWidth-indentRP, maxLim[0]/2-notchWidth+indentRP, -indentRP, maxLim[1]*notch+indentRP ]))
    rigidPlates.append(leftRigidPlate)
    govNodes.append(np.array([ maxLim[0]/2-notchWidth, indent ]))
    govNodesMechBC.append(utilitiesMech.mechanicalBC(dim, -1, leftRigidPlateMechBC))

    #rigid plate left support
    rightRigidPlateMechBC = np.array([2, 0,-1,   -1,-1,-1]) #np.array([2,0,-1])
    rightRigidPlate = utilitiesMech.RigidPlate(-1, 2, np.array([ maxLim[0]/2+notchWidth-indentRP, maxLim[0]/2+notchWidth+indentRP, -indentRP, maxLim[1]*notch+indentRP ]))
    rigidPlates.append(rightRigidPlate)
    govNodes.append(np.array([ maxLim[0]/2+notchWidth, indent ]))
    govNodesMechBC.append(utilitiesMech.mechanicalBC(dim, -2, rightRigidPlateMechBC))


    facesMult = 1
    ########################################## rest of  faces
    nodeA =  np.array([indent , maxLim[1] - indent])
    nodeB =  np.array([maxLim[0] - indent, maxLim[1] - indent])
    pointGenerators.generateNodesLine2dRand(nodeA, nodeB, minDist*facesMult, dim, node_coords,  trials, False, False)

    nodeA =  np.array([indent  ,  indent])
    nodeB =  np.array([maxLim[0] - indent ,  indent])
    pointGenerators.generateNodesLine2dRand(nodeA, nodeB, minDist*facesMult, dim, node_coords,  trials, False, False)

    nodeA =  np.array([indent, indent])
    nodeB =  np.array([indent, maxLim[1] - indent])
    pointGenerators.generateNodesLine2dRand(nodeA, nodeB, minDist*facesMult, dim, node_coords,  trials, False, False)

    nodeA =  np.array([maxLim[0]-indent, indent])
    nodeB =  np.array([maxLim[0]-indent, maxLim[1] - indent])
    pointGenerators.generateNodesLine2dRand(nodeA, nodeB, minDist*facesMult, dim, node_coords,  trials, False, False)


    #rect
    pointGenerators.generateNodesRect(maxLim, minDist, dim, trials, node_coords)


    return node_coords, mechBC_merged, mechInitC_merged, [], govNodes, govNodesMechBC, rigidPlates







def assemble3dCoupledArtificialCrack (maxLim, minDist, trials, slitWidth, notch):
    dim = 3
    #lists for the model
    node_coords = []
    mechBC_merged = []
    mechInitC_merged = []
    govNodes = []
    govNodesMechBC = []
    rigidPlates = []

    #an indent due to mirroring of the data for voronoi tess.
    indent = 1e-8
    notchWidth = slitWidth/2
    #generating notch points

    node_coords.append(np.array([maxLim[0]/2-notchWidth, maxLim[1]*notch/2 , maxLim[2]/2]))
    nodeA = np.array([maxLim[0]/2-notchWidth, indent, indent])
    nodeB = np.array([maxLim[0]/2-notchWidth, maxLim[1]*notch-indent, maxLim[2]-indent])
    oldLen = len(node_coords)
    pointGenerators.generateNodesOrtoSurface3dRand(nodeA, nodeB, minDist, dim, node_coords, trials*10)
    #pointGenerators.generateNodesLine2dRand(nodeA, nodeB, minDist, dim, node_coords, trials, catchCorners=True, equidist=True)

    nodeA = np.array([maxLim[0]/2+notchWidth, indent, indent])
    nodeB = np.array([maxLim[0]/2+notchWidth, maxLim[1]*notch-indent, maxLim[2]-indent])
    oldLen = len(node_coords)
    pointGenerators.generateNodesOrtoSurface3dRand(nodeA, nodeB, minDist, dim, node_coords, trials*10)
    #pointGenerators.generateNodesLine2dRand(nodeA, nodeB, minDist, dim, node_coords, trials, catchCorners=True, equidist=True)

    #node_coords.append(np.array([maxLim[0]/2+notchWidth, maxLim[1]*notch]))
    #node_coords.append(np.array([maxLim[0]/2-notchWidth, maxLim[1]*notch]))


    ##################### CONSTRAINTS AND RIGID PLATES
    #rigid plate left slit face
    indentRP = 1e-6
    leftRigidPlateMechBC = np.array([1, 0,0,  0,0,0, -1,-1,-1, -1,-1,-1])#np.array([1,0,-1])
    leftRigidPlate = utilitiesMech.RigidPlate(-1, 3, np.array([ maxLim[0]/2-notchWidth-indentRP, maxLim[0]/2-notchWidth+indentRP, -indentRP, maxLim[1]*notch+indentRP, -indentRP, maxLim[2]+indentRP ]))
    rigidPlates.append(leftRigidPlate)
    govNodes.append(np.array([ maxLim[0]/2-notchWidth, indent, maxLim[2]/2 ]))
    govNodesMechBC.append(utilitiesMech.mechanicalBC(dim, -1, leftRigidPlateMechBC))

    #rigid plate right slit face
    rightRigidPlateMechBC = np.array([2, 0,0,    0,0,0, -1,-1,-1, -1,-1,-1]) #np.array([2,0,-1])
    rightRigidPlate = utilitiesMech.RigidPlate(-1, 3, np.array([ maxLim[0]/2+notchWidth-indentRP, maxLim[0]/2+notchWidth+indentRP, -indentRP, maxLim[1]*notch+indentRP,  -indentRP, maxLim[2]+indentRP ]))
    rigidPlates.append(rightRigidPlate)
    govNodes.append(np.array([ maxLim[0]/2+notchWidth, indent, maxLim[2]/2 ]))
    govNodesMechBC.append(utilitiesMech.mechanicalBC(dim, -2, rightRigidPlateMechBC))


    facesMult = 1
    ########################################## rest of  faces
    nodeA =  np.array([indent , maxLim[1] - indent, indent])
    nodeB =  np.array([maxLim[0] - indent, maxLim[1] - indent, maxLim[2]-indent])
    pointGenerators.generateNodesOrtoSurface3dRand(nodeA, nodeB, minDist, dim, node_coords, trials)


    nodeA =  np.array([indent ,  indent, indent])
    nodeB =  np.array([maxLim[0] - indent,  indent, maxLim[2]-indent])
    pointGenerators.generateNodesOrtoSurface3dRand(nodeA, nodeB, minDist, dim, node_coords, trials)




    #rect
    pointGenerators.generateNodesRect(maxLim, minDist, dim, trials, node_coords)

    print('generated nodes %d' %len(node_coords))
    return node_coords, mechBC_merged, mechInitC_merged, [], govNodes, govNodesMechBC, rigidPlates

def assemble2dDogBone(D, minDist, trials, excentricity = 50, symmetric=0, edgeMinDistCoef = 1.0, roughDogBone=0, roughEdgeDogbone=0, roughMinDistCoef=1, interLayerThickness=2, powerTes = False, weakboundary = 0 ):
    dim = 2
    #lists for the model
    node_coords = []
    radii = []
    mechBC_merged = []
    mechInitC_merged = []
    govNodes = []
    govNodesMechBC = []
    rigidPlates = []

    indent = 1e-6

    oldLen = len(node_coords)

    if powerTes:
        maxDiam = minDist
        minDiam = minDist/4
        boundary_dist = minDiam * 3
        boundary_radii = minDiam / 2
    else:
        boundary_dist = minDist

    #####################nodes of interest
    node_coords.append(np.array([  D/2,  indent   ])) #top mid
    node_coords.append(np.array([  D/2,  6/4*D - indent  ]))  #bottom mid
    #gauges B
    #if (D==0.1):
    node_coords.append( np.array([ D/2,         3/4*D-D*0.6/2 ])  )#mid LS
    node_coords.append( np.array([ D/2,         3/4*D+D*0.6/2 ])  )
    node_coords.append( np.array([ 0.2*D,       3/4*D-D*0.6/2 ])  ) #left LC (not LC)
    node_coords.append( np.array([ 0.2*D,       3/4*D+D*0.6/2 ])  )
    node_coords.append( np.array([ D-0.2*D,     3/4*D-D*0.6/2 ])  )#right LC (not LC)
    node_coords.append( np.array([ D-0.2*D,     3/4*D+D*0.6/2 ])  )

    #top line of dogbone
    nodeA = np.array([indent, indent])
    nodeB = np.array([D-indent, indent])
    pointGenerators.generateNodesLine2dRand(nodeA, nodeB, boundary_dist, dim,
                                            node_coords, trials, False, False)

    #bottom line of dogbone
    nodeA = np.array([indent,  6/4 * D - indent])
    nodeB = np.array([D-indent, 6/4 * D - indent])
    pointGenerators.generateNodesLine2dRand(nodeA, nodeB, boundary_dist, dim,
                                            node_coords, trials, False, False)

    ##################### CONSTRAINTS AND RIGID PLATES
    #top rigid plate
    indentRP = 1e-3
    topRigidPlateMechBC = np.array([0, 1,-1,   -1,-1,-1])
    topRigidPlate = utilitiesMech.RigidPlate(-1, 2, np.array([
    -indentRP,
    indentRP+D,
    -indentRP,
     +indentRP  ]))
    rigidPlates.append(topRigidPlate)
    govNodes.append(np.array([ D/2+D/excentricity, indent ]))
    govNodesMechBC.append(utilitiesMech.mechanicalBC(dim, -1, topRigidPlateMechBC))
    #bottom rigid plate
    bottomRigidPlateMechBC = np.array([0,0,-1,   -1,-1,-1])
    bottomRigidPlate = utilitiesMech.RigidPlate(-1, 2, np.array([
    -indentRP,
    indentRP+D,
    -indentRP+6/4 * D,
     +indentRP+6/4 * D  ]))
    rigidPlates.append(bottomRigidPlate)
    govNodes.append(np.array([ D/2+D/excentricity, 6/4 * D-indent ]))
    govNodesMechBC.append(utilitiesMech.mechanicalBC(dim, -1, bottomRigidPlateMechBC))
    #####################

    if powerTes:
        radii = [boundary_radii] * len(node_coords) # TODO:

    #node_coords, radii =
        pointGenerators.generateParticlesDogbonePow(D, minDiam, maxDiam, 1, dim, trials,
                                    node_coords, radii, allow_domain_overlap = False,
                                    periodic_distance=False,useLowBound=False)
    else:
        pointGenerators.generateParticlesDogbone(D, minDist, minDist, 1, dim, trials,
                                    node_coords, radii, allow_domain_overlap = False,
                                    periodic_distance=False,useLowBound=False)

    fig, ax = plt.subplots()
    #draw_dogbone(ax, D)

    #ax.scatter(points[:, 0], points[:, 1])
    ax.scatter(np.array(node_coords)[:, 0], np.array(node_coords)[:, 1])

    if powerTes:
        for (x, y), r in zip(node_coords, radii):
            circle = plt.Circle((x, y), r, fill = False )
            ax.add_artist(circle)
    else:
        for (x, y), r in zip(node_coords, np.ones(len(node_coords)) * minDist/2):
            circle = plt.Circle((x, y), r, fill = False )
            ax.add_artist(circle)

    print("asdasd")
    ax.set_aspect('equal')
    #plt.show()
    node_coords=np.asarray(node_coords)
    plt.scatter(node_coords[:,0],node_coords[:,1])
    plt.show()

    node_count = len (node_coords)
    return (node_coords, list(range(len(node_coords))), mechBC_merged, mechInitC_merged,
            node_count, govNodes, govNodesMechBC, rigidPlates, radii)


def assemble2dDogBone_old(D, minDist, trials, excentricity = 50, symmetric=0, edgeMinDistCoef = 1.0, roughDogBone=0, roughEdgeDogbone=0, roughMinDistCoef=1, interLayerThickness=2, powerTes = False, weakboundary = 0 ):

    if roughDogBone >0 :
        sampleCircularBorders = False
    else:
        sampleCircularBorders = True

    dim = 2
    #lists for the model
    node_coords = []
    radii = []
    mechBC_merged = []
    mechInitC_merged = []
    govNodes = []
    govNodesMechBC = []
    rigidPlates = []

    indent = 1e-6

    oldLen = len(node_coords)

    if powerTes:
        maxDiam = minDist
        minDiam = minDist/4
        boundary_dist = minDist * 2
        boundary_radii = minDist/2
    else:
        boundary_dist = minDist

    if roughDogBone > 0:
        altMinDist = roughMinDistCoef * boundary_dist
    else:
        altMinDist = boundary_dist

    if symmetric >0:
        if symmetric == 1:
            nodeA = np.array([0.2*D+2*indent,  3/4 * D - indent -minDist/2])
            nodeB = np.array([D*0.8-2*indent, 3/4 * D - indent -minDist/2])
        if symmetric == 2:
            nodeA = np.array([0.2*D+minDist/2,  3/4 * D - indent -minDist/2])
            nodeB = np.array([D*0.8-minDist/2, 3/4 * D - indent -minDist/2])
        pointGenerators.generateNodesLine2dRand(nodeA, nodeB, minDist, dim, node_coords, trials, True, True)

        if symmetric == 1:
            nodeA = np.array([0.2*D+2*indent,  3/4 * D - indent +minDist/2])
            nodeB = np.array([D*0.8-2*indent, 3/4 * D - indent +minDist/2])
        if symmetric == 2:
            nodeA = np.array([0.2*D+minDist/2,  3/4 * D - indent +minDist/2])
            nodeB = np.array([D*0.8-minDist/2, 3/4 * D - indent +minDist/2])
        pointGenerators.generateNodesLine2dRand(nodeA, nodeB, minDist, dim, node_coords, trials, True, True)



    #####################nodes of interest
    node_coords.append(np.array([  D/2,  indent   ])) #top mid
    node_coords.append(np.array([  D/2,  6/4*D - indent  ]))  #bottom mid
    #gauges B
    #if (D==0.1):
    node_coords.append( np.array([ D/2,         3/4*D-D*0.6/2 ])  )#mid LS
    node_coords.append( np.array([ D/2,         3/4*D+D*0.6/2 ])  )
    node_coords.append( np.array([ 0.2*D,       3/4*D-D*0.6/2 ])  ) #left LC
    node_coords.append( np.array([ 0.2*D,       3/4*D+D*0.6/2 ])  )
    node_coords.append( np.array([ D-0.2*D,     3/4*D-D*0.6/2 ])  )#right LC
    node_coords.append( np.array([ D-0.2*D,     3/4*D+D*0.6/2 ])  )



    #top line of dogbone
    nodeA = np.array([indent, indent])
    nodeB = np.array([D-indent, indent])
    pointGenerators.generateNodesLine2dRand(nodeA, nodeB, altMinDist*0.9, dim, node_coords, trials, False, True)

    #bottom line of dogbone
    nodeA = np.array([indent,  6/4 * D - indent])
    nodeB = np.array([D-indent, 6/4 * D - indent])
    pointGenerators.generateNodesLine2dRand(nodeA, nodeB, altMinDist*0.9, dim, node_coords, trials, False, True)



    #top left edge of dogbone
    nodeA = np.array([indent, indent])
    nodeB = np.array([indent, 1/4*D-indent])
    pointGenerators.generateNodesLine2dRand(nodeA, nodeB, altMinDist, dim, node_coords, trials, True, False)

    #top right edge of dogbone
    nodeA = np.array([D-indent, indent])
    nodeB = np.array([D-indent, 1/4*D-indent])
    pointGenerators.generateNodesLine2dRand(nodeA, nodeB, altMinDist, dim, node_coords, trials, True, False)

    #bottom left edge of dogbone
    nodeA = np.array([indent, 5/4*D-indent])
    nodeB = np.array([indent, 6/4*D-indent])
    pointGenerators.generateNodesLine2dRand(nodeA, nodeB, altMinDist, dim, node_coords, trials, True, False)

    #bottom right edge of dogbone
    nodeA = np.array([D-indent, 5/4*D-indent])
    nodeB = np.array([D-indent, 6/4*D-indent])
    pointGenerators.generateNodesLine2dRand(nodeA, nodeB, altMinDist, dim, node_coords, trials, True, False)

    #top line of dogbone
    nodeA = np.array([indent, indent])
    nodeB = np.array([D-indent, indent])
    pointGenerators.generateNodesLine2dRand(nodeA, nodeB, altMinDist, dim, node_coords, trials, False, True)

    #bottom line of dogbone
    nodeA = np.array([indent,  6/4 * D - indent])
    nodeB = np.array([D-indent, 6/4 * D - indent])
    pointGenerators.generateNodesLine2dRand(nodeA, nodeB, altMinDist, dim, node_coords, trials, False, True)

    uniquePoints =  (len(node_coords)) - oldLen

    ##################### CONSTRAINTS AND RIGID PLATES
    #top rigid plate
    indentRP = 1e-3
    topRigidPlateMechBC = np.array([0, 1,-1,   -1,-1,-1])
    topRigidPlate = utilitiesMech.RigidPlate(-1, 2, np.array([
    -indentRP,
    indentRP+D,
    -indentRP,
     +indentRP  ]))
    rigidPlates.append(topRigidPlate)
    govNodes.append(np.array([ D/2+D/excentricity, indent ]))
    govNodesMechBC.append(utilitiesMech.mechanicalBC(dim, -1, topRigidPlateMechBC))
    #bottom rigid plate
    bottomRigidPlateMechBC = np.array([0,0,-1,   -1,-1,-1])
    bottomRigidPlate = utilitiesMech.RigidPlate(-1, 2, np.array([
    -indentRP,
    indentRP+D,
    -indentRP+6/4 * D,
     +indentRP+6/4 * D  ]))
    rigidPlates.append(bottomRigidPlate)
    govNodes.append(np.array([ D/2+D/excentricity, 6/4 * D-indent ]))
    govNodesMechBC.append(utilitiesMech.mechanicalBC(dim, -1, bottomRigidPlateMechBC))
    #####################



    centreA = np.array( [-0.525 * D, 3/4 * D] )
    centreB = np.array( [ 1.525 * D, 3/4 * D] )
    if roughEdgeDogbone==0:
        #sampling on circular borders
        radius = 0.725*D
        angleLimitA =   -np.arcsin( 0.5*D / radius)
        angleLimitB =   np.arcsin( 0.5*D / radius)
        if roughDogBone >1:
            angleLimitA =   -np.arcsin(  roughDogBone * minDist / radius)
            angleLimitB =   np.arcsin(  roughDogBone * minDist / radius)
        #if symmetric == True:
        #    angleLimitB =   0
        mirroredPointsA =  pointGenerators.generateNodesCircle2dRand(centreA, radius+indent, minDist*edgeMinDistCoef, node_coords, trials, angleLimitA=angleLimitA, angleLimitB=angleLimitB, mirrorIndent = indent*10 )
        angleLimitA +=   np.pi
        angleLimitB +=   np.pi
        #if symmetric == True:
        #    angleLimitA =    np.pi
        mirroredPointsB =  pointGenerators.generateNodesCircle2dRand(centreB, radius+indent, minDist*edgeMinDistCoef, node_coords, trials, angleLimitA=angleLimitA, angleLimitB=angleLimitB, mirrorIndent = indent*10)
    else:
        sampleCircularBorders = False
        mirroredPointsA = []
        mirroredPointsB = []

    #"random sampling along the border"
    if roughEdgeDogbone==3:
        radius = 0.725*D
        radiusSpread = minDist
        angleLimitA =   -np.arcsin( 0.5*D / radius)
        angleLimitB =   np.arcsin( 0.5*D / radius)
        if roughDogBone >1:
            angleLimitA =   -np.arcsin( 2* roughDogBone * minDist / radius)
            angleLimitB =   np.arcsin( 2* roughDogBone * minDist / radius)
        pointGenerators.generateNodesCircle2dRand(centreA, radius+radiusSpread, minDist, node_coords, trials, angleLimitA=angleLimitA, angleLimitB=angleLimitB, mirrorIndent = 0, radiusSpread = radiusSpread )
        angleLimitA +=   np.pi
        angleLimitB +=   np.pi
        pointGenerators.generateNodesCircle2dRand(centreB, radius+radiusSpread, minDist, node_coords, trials, angleLimitA=angleLimitA, angleLimitB=angleLimitB, mirrorIndent =0,
        radiusSpread = radiusSpread)
    """
        node_coords = np.asarray(node_coords)
        plt.plot(node_coords[:,0], node_coords[:,1], 'o', color='black');
        plt.show()
    #"""

    # boundary with diferent discretization
    if weakboundary > 0:
        print('Assembling weak boundary...\n', end='')

        boundary_minDist = 0.5 * minDist
        radius = 0.725*D
        maxLim_left = [0.2 * D + weakboundary, 1.5 * D]
        node_coords_BoundaryRect = []
        node_coords_inBoundary = []
        pointGenerators.generateNodesRect(maxLim_left, boundary_minDist, 2, trials, node_coords_BoundaryRect)

        for i in range(len(node_coords_BoundaryRect)):
            node = node_coords_BoundaryRect[i]
            distA = np.linalg.norm( node - centreA)
            distB = np.linalg.norm( node - centreB)
            if (distA > radius and distA < (radius + weakboundary) or  distB > radius and distB < (radius + weakboundary)):
                # LEFT BOUNDARY
                node_coords.append(np.array(node))
                # RIGHT BOUNDARY
                node_coords.append(np.array(np.array([D,0]) + node * np.array([-1,1])))

    if roughDogBone == 1: #hrubsi jen obdelniky prilozek
        #top rough rectangle
        oldLen = len(node_coords)
        maxLim = np.array([  D    ,  1/4*D ])
        pointGenerators.generateNodesRect(maxLim, altMinDist, 2, trials, node_coords)
        #bottom rough rectangle
        oldLen = len(node_coords)
        maxLimF = np.array([     indent,       5/4 * D,        D,        6/4 * D   ])
        pointGenerators.generateNodesRect(maxLimF, altMinDist, 2, trials, node_coords, useLowBound=True)
        # middle fine asssemble3dPeriodicRectanglemaxLimF = np.array([     indent,       5/4 * D,        D,        6/4 * D   ])
        maxLimF = np.array([     indent,       1/4 * D,        D,        5/4 * D   ])
        pointGenerators.generateNodesRect(maxLimF, minDist, 2, trials, node_coords, useLowBound=True)
        nrOfPoints =  (len(node_coords)) - oldLen

    elif roughDogBone > 1: #hrubsi krome pruhu +-10xmindist od prostredka
        #top rough rectangle
        oldLen = len(node_coords)
        maxLim = np.array([  D    ,  3/4*D  - (roughDogBone+interLayerThickness)*minDist ])
        pointGenerators.generateNodesRect(maxLim, altMinDist, 2, trials, node_coords)
        #bottom rough rectangle
        oldLen = len(node_coords)
        maxLimF = np.array([     indent,       3/4 * D+   (roughDogBone+interLayerThickness)*minDist,        D,        6/4 * D   ])
        pointGenerators.generateNodesRect(maxLimF, altMinDist, 2, trials, node_coords, useLowBound=True)


        maxLimF = np.array([     indent,       3/4*D  - (roughDogBone+interLayerThickness)*minDist  ,        D,        3/4*D  - (roughDogBone)*minDist ])
        pointGenerators.generateNodesRect(maxLimF, minDist*4, 2, trials, node_coords, useLowBound=True, topMinDist = altMinDist, bottomMinDist = minDist, gradienDirection=1)

        maxLimF = np.array([     indent,       3/4*D  + roughDogBone*minDist  ,        D,        3/4*D  + (roughDogBone+interLayerThickness)*minDist ])
        pointGenerators.generateNodesRect(maxLimF, minDist*4, 2, trials, node_coords, useLowBound=True, topMinDist = minDist, bottomMinDist = altMinDist, gradienDirection=1)

        # middle fine
        maxLimF = np.array([     indent,       3/4*D  - (roughDogBone)*minDist,        D,        3/4*D  + (roughDogBone)*minDist  ])
        pointGenerators.generateNodesRect(maxLimF, minDist, 2, trials, node_coords, useLowBound=True)

        nrOfPoints =  (len(node_coords)) - oldLen

    else:
        #rectangle of dogbone
        oldLen = len(node_coords)
        maxLim = np.array([  D    ,  6/4*D ])
        #if symmetric == True:
        #    maxLim = np.array([  D    ,  3/4*D ])

        if powerTes:
            radii = [boundary_radii for r in range(len(node_coords))]
            node_coords, radii = pointGenerators.generateParticlesRect(maxLim, minDiam, maxDiam, 0.8, 2, trials, np.array(node_coords), np.array(radii), allow_domain_overlap = False, periodic_distance=False)
        else:
            pointGenerators.generateNodesRect(maxLim, minDist, 2, trials, node_coords)
            # print(list(node_coords))
            # pointGenerators.generateNodesRect_old(maxLim, minDist, 2, trials, list(node_coords))
        nrOfPoints =  (len(node_coords)) - oldLen



    node_coords_all = np.copy ( node_coords )
    node_coords_dogbone = []
    node_indices_dogbone = []
    radii_dogbone = []

    #dumping points outside bone

    radius = np.linalg.norm( centreB - np.array([D, 1/4*D]))
    print('Dumping points within bordering circles...', end='')
    node_coords_out = []
    for i in range(len(node_coords_all)):
        node = node_coords_all[i]
        distA = np.linalg.norm( node - centreA)
        distB = np.linalg.norm( node - centreB)
        if (distA > radius and distB > radius):
            node_indices_dogbone.append(i)
            node_coords_dogbone.append(node)
            if len(radii) > 0:
                radii_dogbone.append(radii[i])
    print('done.')
    node_coords_dogbone = np.asarray(node_coords_dogbone)
    radii_dogbone = np.asarray(radii_dogbone)


    if roughEdgeDogbone == 2 or roughEdgeDogbone == 3:
        mirrored_coords = []
        mirrored_radii = []
        #mirroring rough edge dogbone circular borders
        dogboneRadius = 0.725*D
        leftCenter = np.array( [-0.525 * D, 3/4 * D] )
        rightCenter = np.array( [ 1.525 * D, 3/4 * D] )
        for i, node in enumerate(node_coords_dogbone):
            if node[0]<D/2:
                #left half, mirroring to left center
                nodeRad = np.linalg.norm(leftCenter-node)
                distFromEdge = nodeRad - dogboneRadius
                mirroredNodeRad = nodeRad - 2*distFromEdge
                #
                relativeNodeCoords = node - leftCenter
                mirroredNodeRelativeCoords = relativeNodeCoords * (mirroredNodeRad/nodeRad)
                mirroredNodeAbsoluteCoords = mirroredNodeRelativeCoords + leftCenter
                #
                if (mirroredNodeAbsoluteCoords[0]>=0):
                    mirrored_coords.append(mirroredNodeAbsoluteCoords)
                    if len(radii) > 0:
                        mirrored_radii.append(radii_dogbone[i])

                    """
                    print()
                    print ('dist from edge %s' %distFromEdge)
                    print ('mirroredNodeRad %s' %mirroredNodeRad)
                    print ('nodeRad %s' %nodeRad)
                    print ('dogboneRadius %s' %dogboneRadius)
                    plt.plot(node_coords_dogbone[:,0], node_coords_dogbone[:,1], '.', color='black');
                    plt.plot(leftCenter[0], leftCenter[1], 'o', color='black');
                    plt.plot(mirroredNodeAbsoluteCoords[0], mirroredNodeAbsoluteCoords[1], 'o', color='red');
                    plt.plot(node[0], node[1], 'o', color='blue');
                    plt.show()
                    """
            else:
                #right half, mirroring to left center
                nodeRad = np.linalg.norm(rightCenter-node)
                distFromEdge = nodeRad - dogboneRadius
                mirroredNodeRad = nodeRad - 2*distFromEdge
                #
                relativeNodeCoords = rightCenter - node
                mirroredNodeRelativeCoords = relativeNodeCoords * (mirroredNodeRad/nodeRad)
                mirroredNodeAbsoluteCoords = rightCenter - mirroredNodeRelativeCoords
                #
                if (mirroredNodeAbsoluteCoords[0]<=D):
                    mirrored_coords.append(mirroredNodeAbsoluteCoords)
                    if len(radii) > 0:
                        mirrored_radii.append(radii_dogbone[i])

        mirrored_coords = np.asarray(mirrored_coords)
        mirrored_radii = np.asarray(mirrored_radii)

        """
        plt.plot(mirrored_coords[:,0], mirrored_coords[:,1], 'o', color='black');
        plt.plot(node_coords_dogbone[:,0], node_coords_dogbone[:,1], 'x', color='red');
        plt.show()
        #"""
        node_coords_dogbone = np.vstack((node_coords_dogbone, mirrored_coords))
        node_coords_all = np.copy(node_coords_dogbone)
        radii = np.hstack((radii_dogbone, mirrored_radii))

        node_indices_dogbone = []
        for i in range(len(node_coords_all)):
            node = node_coords_all[i]
            distA = np.linalg.norm( node - centreA)
            distB = np.linalg.norm( node - centreB)
            if (distA > radius and distB > radius):
                node_indices_dogbone.append(i)




    if sampleCircularBorders == True:
        mirroredMiddle = []
        mirroredMiddle.append(np.array([  0.2*D-1e-6,  3/4 * D - indent  +minDist/2 ]))
        mirroredMiddle.append(np.array([  D*0.8+1e-6,  3/4 * D - indent  +minDist/2 ]))
        mirroredMiddle.append(np.array([  0.2*D-1e-6,  3/4 * D - indent  -minDist/2 ]))
        mirroredMiddle.append(np.array([  D*0.8+1e-6,  3/4 * D - indent  -minDist/2 ]))

        node_coords_all = np.vstack( (node_coords_all, mirroredPointsA, mirroredPointsB) )


    # plt.plot(mirrored_coords[:, 0], mirrored_coords[:, 1], 'o', color='blue')
    # plt.plot(node_coords_all[:,0], node_coords_all[:,1], 'o', color='green');
    # plt.plot(node_coords_all[node_indices_dogbone,0], node_coords_all[node_indices_dogbone,1], 'x', color='red');
    # plt.show()
    #"""

    node_count = len (node_coords_all)
    return node_coords_all, node_indices_dogbone, mechBC_merged, mechInitC_merged, node_count, govNodes, govNodesMechBC, rigidPlates, radii




def assemble2dDogBoneStrip(D, minDist, excentricity = 50 , randomStrip = 0):

    dim = 2
    #lists for the model
    node_coords = []
    mechBC_merged = []
    mechInitC_merged = []
    govNodes = []
    govNodesMechBC = []
    rigidPlates = []

    indent = 1e-6
    trials = 500

    if randomStrip == 0:
        print ('Ekvidist')
        nodeA = np.array([0.2*D+2*indent,  indent])
        nodeB = np.array([D*0.8-2*indent, indent])
        pointGenerators.generateNodesLine2dRand(nodeA, nodeB, minDist, dim, node_coords, trials, True, True)

        nodeA = np.array([0.2*D+2*indent,  minDist-indent])
        nodeB = np.array([D*0.8-2*indent, minDist-indent])
        pointGenerators.generateNodesLine2dRand(nodeA, nodeB, minDist, dim, node_coords, trials, True, True)

    if randomStrip == 1 :
        print ('Random w corners')
        nodeA = np.array([0.2*D+2*indent,  indent])
        nodeB = np.array([D*0.8-2*indent, indent])
        pointGenerators.generateNodesLine2dRand(nodeA, nodeB, minDist, dim, node_coords, trials, True, False)

        secondLine = []
        nodeA = np.array([0.2*D+2*indent,  minDist-indent])
        nodeB = np.array([D*0.8-2*indent, minDist-indent])
        pointGenerators.generateNodesLine2dRand(nodeA, nodeB, minDist, dim, secondLine, trials, True, False)

        for p in secondLine:
            node_coords.append(p)

    if randomStrip == 2 :
        print ('Random')
        nodeA = np.array([0.2*D+2*indent,  indent])
        nodeB = np.array([D*0.8-2*indent, indent])
        pointGenerators.generateNodesLine2dRand(nodeA, nodeB, minDist, dim, node_coords, trials, False, False)

        secondLine = []
        nodeA = np.array([0.2*D+2*indent,  minDist-indent])
        nodeB = np.array([D*0.8-2*indent, minDist-indent])
        pointGenerators.generateNodesLine2dRand(nodeA, nodeB, minDist, dim, secondLine, trials, False, False)

        for p in secondLine:
            node_coords.append(p)


    ##################### CONSTRAINTS AND RIGID PLATES
    #top rigid plate
    indentRP = 1e-5
    topRigidPlateMechBC = np.array([0, 1,-1,   -1,-1,-1])
    topRigidPlate = utilitiesMech.RigidPlate(-1, 2, np.array([
    -indentRP,
    indentRP+D,
    -indentRP,
     +indentRP  ]))
    rigidPlates.append(topRigidPlate)
    govNodes.append(np.array([ D/2+D/excentricity, 0 ]))
    govNodesMechBC.append(utilitiesMech.mechanicalBC(dim, -1, topRigidPlateMechBC))
    #bottom rigid plate
    bottomRigidPlateMechBC = np.array([0,0,-1,   -1,-1,-1])
    bottomRigidPlate = utilitiesMech.RigidPlate(-1, 2, np.array([
    -indentRP,
    indentRP+D,
    -indentRP+minDist,
     +indentRP+minDist ]))
    rigidPlates.append(bottomRigidPlate)
    govNodes.append(np.array([ D/2+D/excentricity, minDist ]))
    govNodesMechBC.append(utilitiesMech.mechanicalBC(dim, -1, bottomRigidPlateMechBC))
    #####################

    """
    plt.plot(node_coords[:,0], node_coords[:,1], 'o', color='green');
    plt.show()
    """
    node_count = len (node_coords)
    return node_coords, mechBC_merged, mechInitC_merged, node_count, govNodes, govNodesMechBC, rigidPlates




def assemble2dDogBoneBand(maxLim, minDist, roughMinDistCoef=1, elasticHeightCoef=1 ):

    dim = 2
    #lists for the model
    node_coords = []
    mechBC_merged = []
    mechInitC_merged = []
    govNodes = []
    govNodesMechBC = []
    rigidPlates = []

    indent = 1e-7
    trials = 1000


    elasticHeight = maxLim[1] * elasticHeightCoef / (2*elasticHeightCoef + 1)
    dmgHeight = maxLim[1] * 1 / (2*elasticHeightCoef + 1)
    width = maxLim[0]

    #top elastic rect
    nodeA = np.array([indent,  (elasticHeight+dmgHeight+elasticHeight)-indent])
    nodeB = np.array([width-indent, (elasticHeight+dmgHeight+elasticHeight)-indent])
    pointGenerators.generateNodesLine2dRand(nodeA, nodeB, minDist, dim, node_coords, trials, True, True)
    bounds = np.array([     indent,          (elasticHeight+dmgHeight),      width-indent,      (elasticHeight+dmgHeight+elasticHeight)-indent    ])
    pointGenerators.generateNodesRect(bounds, minDist, dim, trials, node_coords, useLowBound=True, topMinDist = minDist, bottomMinDist = minDist*roughMinDistCoef, gradienDirection=1)

    #mid dmg rect
    bounds = np.array([     indent,          (elasticHeight),      width-indent,       (elasticHeight+dmgHeight)    ])
    pointGenerators.generateNodesRect(bounds, minDist, dim, trials, node_coords, useLowBound=True)

    #bottom elastic rect
    nodeA = np.array([indent, indent])
    nodeB = np.array([width-indent, indent])
    pointGenerators.generateNodesLine2dRand(nodeA, nodeB, minDist, dim, node_coords, trials, True, True)
    bounds = np.array([     indent,          indent,      width-indent,      elasticHeight  ])
    pointGenerators.generateNodesRect(bounds, minDist, dim, trials, node_coords, useLowBound=True, bottomMinDist = minDist, topMinDist = minDist*roughMinDistCoef, gradienDirection=1)



    ##################### CONSTRAINTS AND RIGID PLATES
    #top rigid plate
    indentRP = 1e-6
    topRigidPlateMechBC = np.array([0, 1,0,   -1,-1,-1])
    topRigidPlate = utilitiesMech.RigidPlate(-1, 2, np.array([
    -indentRP,
    indentRP+width,
    -indentRP+(elasticHeight+dmgHeight+elasticHeight),
     +indentRP+(elasticHeight+dmgHeight+elasticHeight)  ]))
    rigidPlates.append(topRigidPlate)
    govNodes.append(np.array([ width/2, (elasticHeight+dmgHeight+elasticHeight) ]))
    govNodesMechBC.append(utilitiesMech.mechanicalBC(dim, -1, topRigidPlateMechBC))

    #bottom rigid plate
    bottomRigidPlateMechBC = np.array([0,0,0,   -1,-1,-1])
    bottomRigidPlate = utilitiesMech.RigidPlate(-1, 2, np.array([
    -indentRP,
    indentRP+width,
    -indentRP,
     +indentRP ]))
    rigidPlates.append(bottomRigidPlate)
    govNodes.append(np.array([ width/2, 0 ]))
    govNodesMechBC.append(utilitiesMech.mechanicalBC(dim, -2, bottomRigidPlateMechBC))
    #####################

    """
    node_coords = np.asarray(node_coords)
    plt.plot(node_coords[:,0], node_coords[:,1], 'o', color='green');
    plt.show()
    """

    node_count = len (node_coords)
    return node_coords, mechBC_merged, mechInitC_merged, node_count, govNodes, govNodesMechBC, rigidPlates



def assemble3dDogBone(D, minDist, trials, thickness, D0=None,H=None,H0=None,  excentricity_X = 20, excentricity_Z = 0, symmetric=0):
    dim = 3
    #lists for the model
    node_coords = []
    mechBC_merged = []
    mechInitC_merged = []
    node_coords_out = []
    govNodes = []
    govNodesMechBC = []
    rigidPlates = []

    indent = 1e-5

    if H is None:
        H = 6/4*D

    if symmetric == 1:
        remainerX = (int(0.6 * D * 1000) % int(minDist * 1000)) * 0.001 # ! float to integer conversion due to the rounding error !
        remainerZ = ((int(thickness * 1000) - int(minDist * 1000)) % int(minDist * 1000)) * 0.001 # ! float to integer conversion due to the rounding error !
        grainsN_X = int((0.6*D - minDist - remainerX) / minDist)+1 # number of grains in x direction
        grainsN_Z = int((thickness - minDist - remainerZ) / minDist)+1 # number of grains in z direction
        distX = minDist + (remainerX / (grainsN_X - 1)) # distance between grains in x direction
        distZ = minDist + (remainerZ / (grainsN_Z - 1)) # distance between grains in z direction

        for IndexOfRow in range(grainsN_Z):
            nodeA = np.array([0.2 * D + minDist / 2 + 2 * indent, 3 / 4 * D - minDist/2 + indent, minDist / 2 + IndexOfRow * distZ])
            nodeB = np.array([0.8 * D - minDist / 2 - 2 * indent, 3 / 4 * D - minDist/2 + indent, minDist / 2 + IndexOfRow * distZ])
            pointGenerators.generateNodesLine3dRand(nodeA, nodeB, distX, dim, node_coords, trials, True, True)

            nodeC = np.array([0.2 * D + minDist / 2 + 2 * indent, 3 / 4 * D + minDist/2 - indent, minDist / 2 + IndexOfRow * distZ])
            nodeD = np.array([0.8 * D - minDist / 2 - 2 * indent, 3 / 4 * D + minDist/2 - indent, minDist / 2 + IndexOfRow * distZ])
            pointGenerators.generateNodesLine3dRand(nodeC, nodeD, distX, dim, node_coords, trials, True, True)

            #print(f"nodes coords = {IndexOfRow}{[nodeA, nodeB, nodeC, nodeD]}")

    node_coords.append(np.array([(D-indent)/2, indent, thickness/2]))
    node_coords.append(np.array([(D-indent)/2, H - indent, thickness/2]))

    ##################### CONSTRAINTS AND RIGID PLATES
    #top rigid plate
    indentRP = 1e-3
    topRigidPlateMechBC = np.array([0,1,0,  -1,-1,-1, -1,-1,-1, -1,-1,-1,])
    topRigidPlate = utilitiesMech.RigidPlate(-1, 3, np.array([
    -indentRP,
    indentRP+D,
    -indentRP,
     +indentRP,
     -indentRP,
     thickness+indentRP  ]))
    rigidPlates.append(topRigidPlate)
    govNodes.append(np.array([D / 2 + D / excentricity_X, indent, thickness / 2 + excentricity_Z]))
    govNodesMechBC.append(utilitiesMech.mechanicalBC(dim, -1, topRigidPlateMechBC))
    #bottom rigid plate
    bottomRigidPlateMechBC =  np.array([0,0,0,  -1,-1,-1, -1,-1,-1, -1,-1,-1,])
    bottomRigidPlate = utilitiesMech.RigidPlate(-1, 3, np.array([
    -indentRP,
    indentRP+D,
    -indentRP+H,
     +indentRP+H,
     -indentRP,
     thickness+indentRP  ]))
    rigidPlates.append(bottomRigidPlate)
    govNodes.append(np.array([D / 2 + D / excentricity_X, 6 / 4 * D - indent, thickness / 2 + excentricity_Z]))
    govNodesMechBC.append(utilitiesMech.mechanicalBC(dim, -1, bottomRigidPlateMechBC))
    #####################

    #top surface of dogbone
    nodeA = np.array([indent, indent, indent])
    nodeB = np.array([D-indent, indent, thickness-indent ])
    pointGenerators.generateNodesOrtoSurface3dRand(nodeA, nodeB, minDist, dim, node_coords, trials)

    #bottom surface of dogbone
    nodeA = np.array([indent,  H - indent, indent])
    nodeB = np.array([D-indent, H - indent, thickness-indent])
    pointGenerators.generateNodesOrtoSurface3dRand(nodeA, nodeB, minDist, dim, node_coords, trials)

    # block
    oldLen = len(node_coords)
    maxLim = np.array([D,  H, thickness])
    radii = np.zeros(len(node_coords))+minDist/2
    #pointGenerators.generateNodesRect(maxLim, minDist, 3, trials, node_coords)
    node_coords, radii = pointGenerators.generateParticlesRect(maxLim, minDist*0.5, minDist, 0.9, dim, trials, np.asarray(node_coords), radii, periodic_distance=False)
    
    
    if H0 is None:
        #Van Mier dogbone geometry
        centreA = np.array( [-0.525 * D, 3/4 * D, indent] )
        centreB = np.array( [ 1.525 * D, 3/4 * D, indent] )
        radius = np.linalg.norm( centreB[0:2] - np.array([D, 1/4*D]))
        dogboneRadius = radius

        leftCenter = np.array([-0.525 * D, 3 / 4 * D])
        rightCenter = np.array([1.525 * D, 3 / 4 * D])
    else:
        #Custom geometry
        def circumcenter(A, B, C):
            x1, y1 = A
            x2, y2 = B
            x3, y3 = C
            D = 2 * (x1*(y2 - y3) - x2*(y1 - y3) + x3*(y1 - y2))
            Ux = ((x1**2 + y1**2)*(y2 - y3) + (x2**2 + y2**2)*(y3 - y1) + (x3**2 + y3**2)*(y1 - y2)) / D
            Uy = ((x1**2 + y1**2)*(x3 - x2) + (x2**2 + y2**2)*(x1 - x3) + (x3**2 + y3**2)*(x2 - x1)) / D
            radius = math.sqrt((Ux - x1)**2 + (Uy - y1)**2)
            return np.asarray([Ux, Uy, thickness/2]),radius
        centreA,radius = circumcenter([0, (H-H0)/2], [0, H-(H-H0)/2], [(D-D0)/2 , H/2] )
        dogboneRadius = radius
        centreB,_ = circumcenter([D, (H-H0)/2], [D, H-(H-H0)/2], [D-(D-D0)/2 , H/2] )
        leftCenter = centreA[:2]
        rightCenter = centreB[:2]
    
    node_coords_out = []
    radii_out = []
    node_indices_dogbone = []
    i = 0

    for ni, node in enumerate(node_coords):
        distA = np.linalg.norm( node[0:2] - centreA[0:2])
        distB = np.linalg.norm( node[0:2] - centreB[0:2])
        in_dogboneBLock =  0 < node[0] < D and 0 < node[2] < thickness
        if in_dogboneBLock and (distA > radius and distB > radius):
            node_coords_out.append(node)
            radii_out.append(radii[ni])
            node_indices_dogbone.append(i)
            i += 1

    node_count = len(node_coords)

    # mirroring rough edge dogbone circular borders
    mirrored_coords = []
    mirrored_radii = []
     
    for ni,node in enumerate(node_coords_out):
        if node[0] < D / 2:
            # left half, mirroring to left center
            nodeRad = np.linalg.norm(leftCenter - node[0:2])
            distFromEdge = nodeRad - dogboneRadius
            mirroredNodeRad = nodeRad - 2 * distFromEdge
            #
            relativeNodeCoords = node[0:2] - leftCenter
            mirroredNodeRelativeCoords = relativeNodeCoords * (mirroredNodeRad / nodeRad)
            mirroredNodeAbsoluteCoords = mirroredNodeRelativeCoords + leftCenter
            mirroredNodeAbsoluteCoords = np.append(mirroredNodeAbsoluteCoords, node[2])
            #
            if (mirroredNodeAbsoluteCoords[0] >= 0):
                mirrored_coords.append(mirroredNodeAbsoluteCoords)
                mirrored_radii.append(radii[ni])
        else:
            # right half, mirroring to left center
            nodeRad = np.linalg.norm(rightCenter - node[0:2])
            distFromEdge = nodeRad - dogboneRadius
            mirroredNodeRad = nodeRad - 2 * distFromEdge
            #
            relativeNodeCoords = rightCenter - node[0:2]
            mirroredNodeRelativeCoords = relativeNodeCoords * (mirroredNodeRad / nodeRad)
            mirroredNodeAbsoluteCoords = rightCenter - mirroredNodeRelativeCoords
            mirroredNodeAbsoluteCoords = np.append(mirroredNodeAbsoluteCoords, node[2])
            #
            if (mirroredNodeAbsoluteCoords[0] <= D):
                mirrored_coords.append(mirroredNodeAbsoluteCoords)
                mirrored_radii.append(radii[ni])

    mirrored_coords = np.asarray(mirrored_coords)
    mirrored_radii = np.asarray(mirrored_radii)
    node_coords_out = np.vstack((node_coords_out, mirrored_coords))
    radii_out = np.hstack((radii_out, mirrored_radii))
    
    #plt.scatter(np.asarray(node_coords_out)[node_indices_dogbone,0],np.asarray(node_coords_out)[node_indices_dogbone,1])
    #plt.show()
    return node_coords_out, node_indices_dogbone, mechBC_merged, mechInitC_merged, node_count,govNodes, govNodesMechBC, rigidPlates, radii_out













def asssemble2dPeriodicShear (maxLim, minDist, trials, powerTes):
    dim = 2
    #lists for the model
    node_coords = []
    mechBC_merged = []
    mechInitC_merged = []

    print('assembling 2d periodic ')
    ###########generating of points in rectangle


    if powerTes == False:
        pointGenerators.generateNodesRectPeriodic(maxLim, minDist, dim, trials, node_coords)
        radii = np.zeros(len(node_coords))
    else:
        #TODO: power Tesselation
        node_coords = np.zeros((0,dim))
        radii = np.zeros(0)
        node_coords, radii = pointGenerators.generateParticlesRect(maxLim, minDist*0.4, minDist, 0.8, dim, trials, node_coords, radii, allow_domain_overlap = True, periodic_distance=True)

    node_coords = np.asarray(node_coords)
    #masters = np.ones(len(node_coords)).astype(int)*(-1)

    limit = 8*minDist
    XA = np.where(node_coords[:,0]<limit)[0]
    XB = np.where(node_coords[:,0]>maxLim[0]-limit)[0]
    YA = np.where(node_coords[:,1]<limit)[0]
    YB = np.where(node_coords[:,1]>maxLim[1]-limit)[0]

    XAYA = XA[np.where(node_coords[XA,1]<limit)[0]]
    XAYB = XA[np.where(node_coords[XA,1]>maxLim[1]-limit)[0]]
    XBYA = XB[np.where(node_coords[XB,1]<limit)[0]]
    XBYB = XB[np.where(node_coords[XB,1]>maxLim[1]-limit)[0]]

    nNds = np.vstack((
    node_coords,
    node_coords[XA] + np.array([maxLim[0], 0]),
    node_coords[YA] + np.array([0, maxLim[1]]),
    #
    node_coords[XB] + np.array([-maxLim[0], 0]),
    node_coords[YB] + np.array([0, -maxLim[1]]),
    #
    node_coords[XAYA] + np.array([maxLim[0], maxLim[1]]),
    node_coords[XBYA] + np.array([-maxLim[0], maxLim[1]]),
    node_coords[XAYB] + np.array([maxLim[0], -maxLim[1]]),
    node_coords[XBYB] + np.array([-maxLim[0], -maxLim[1]])
    ))

    #masters = np.hstack(( masters,XA,YA,XB,YB,XAYA,XBYA,XAYB,XBYB ))
    radii = np.hstack(( radii, radii[np.hstack((XA,YA,XB,YB,XAYA,XBYA,XAYB,XBYB ))]))


    nNds = np.asarray(nNds)
    # if True:
        #ax.auto_scale_xyz([-maxLim[0], 2*maxLim[0]], [-maxLim[1], 2*maxLim[1]], [-maxLim[2], 2*maxLim[2]])
        #ax.scatter(node_coords[:,0], node_coords[:,1], node_coords[:,2], color='r')
        # plt.scatter(nNds[:,0], nNds[:,1],  color='r')
        # plt.show()




    return nNds, mechBC_merged, mechInitC_merged, radii#, masters




def asssemble3dPeriodicRectangle (maxLim, minDist, trials, powerTes):
    dim = 3
    #lists for the model
    node_coords = []
    mechBC_merged = []
    mechInitC_merged = []


    print('assembling 3d periodic ')
    ###########generating of points in rectangle
    if powerTes == False:
        pointGenerators.generateNodesRectPeriodic(maxLim, minDist, dim, trials, node_coords)

    else:
        node_coords = np.zeros((0,dim))
        radii = np.zeros(len(node_coords))
        node_coords, radii = pointGenerators.generateParticlesRect(maxLim, minDist*0.4, minDist, 0.8, dim, trials, node_coords, radii, allow_domain_overlap = True, periodic_distance=True)


    node_coords = np.asarray(node_coords)
    #masters = np.ones(len(node_coords)).astype(int)*(-1)

    limit = 8*minDist
    XA = np.where(node_coords[:,0]<limit)[0]
    XB = np.where(node_coords[:,0]>maxLim[0]-limit)[0]
    YA = np.where(node_coords[:,1]<limit)[0]
    YB = np.where(node_coords[:,1]>maxLim[1]-limit)[0]
    ZA = np.where(node_coords[:,2]<limit)[0]
    ZB = np.where(node_coords[:,2]>maxLim[2]-limit)[0]

    XAYA = XA[np.where(node_coords[XA,1]<limit)[0]]
    XAYB = XA[np.where(node_coords[XA,1]>maxLim[1]-limit)[0]]
    XBYA = XB[np.where(node_coords[XB,1]<limit)[0]]
    XBYB = XB[np.where(node_coords[XB,1]>maxLim[1]-limit)[0]]

    XAZA = XA[np.where(node_coords[XA,2]<limit)[0]]
    XAZB = XA[np.where(node_coords[XA,2]>maxLim[2]-limit)[0]]
    XBZA = XB[np.where(node_coords[XB,2]<limit)[0]]
    XBZB = XB[np.where(node_coords[XB,2]>maxLim[2]-limit)[0]]

    YAZA = YA[np.where(node_coords[YA,2]<limit)[0]]
    YAZB = YA[np.where(node_coords[YA,2]>maxLim[2]-limit)[0]]
    YBZA = YB[np.where(node_coords[YB,2]<limit)[0]]
    YBZB = YB[np.where(node_coords[YB,2]>maxLim[2]-limit)[0]]

    XAYAZA = XAYA[np.where(node_coords[XAYA,2]<limit)[0]]
    XAYBZA = XAYB[np.where(node_coords[XAYB,2]<limit)[0]]
    XBYAZA = XBYA[np.where(node_coords[XBYA,2]<limit)[0]]
    XBYBZA = XBYB[np.where(node_coords[XBYB,2]<limit)[0]]
    XAYAZB = XAYA[np.where(node_coords[XAYA,2]>maxLim[2]-limit)[0]]
    XAYBZB = XAYB[np.where(node_coords[XAYB,2]>maxLim[2]-limit)[0]]
    XBYAZB = XBYA[np.where(node_coords[XBYA,2]>maxLim[2]-limit)[0]]
    XBYBZB = XBYB[np.where(node_coords[XBYB,2]>maxLim[2]-limit)[0]]

    nNds = np.vstack((
    node_coords,
    node_coords[XA] + np.array([maxLim[0], 0, 0]),
    node_coords[YA] + np.array([0, maxLim[1], 0]),
    node_coords[ZA] + np.array([0, 0, maxLim[2]]),
    #
    node_coords[XB] + np.array([-maxLim[0], 0, 0]),
    node_coords[YB] + np.array([0, -maxLim[1], 0]),
    node_coords[ZB] + np.array([0, 0, -maxLim[2]]),
    #
    node_coords[XAYA] + np.array([maxLim[0], maxLim[1], 0]),
    node_coords[XBYA] + np.array([-maxLim[0], maxLim[1], 0]),
    node_coords[XAYB] + np.array([maxLim[0], -maxLim[1], 0]),
    node_coords[XBYB] + np.array([-maxLim[0], -maxLim[1], 0]),
    #
    node_coords[XAZA] + np.array([maxLim[0], 0, maxLim[2]]),
    node_coords[XBZA] + np.array([-maxLim[0], 0, maxLim[2]]),
    node_coords[XAZB] + np.array([maxLim[0], 0, -maxLim[2]]),
    node_coords[XBZB] + np.array([-maxLim[0], 0, -maxLim[2]]),
    #
    node_coords[YAZA] + np.array([0, maxLim[1], maxLim[2]]),
    node_coords[YBZA] + np.array([0, -maxLim[1], maxLim[2]]),
    node_coords[YAZB] + np.array([0, maxLim[1], -maxLim[2]]),
    node_coords[YBZB] + np.array([0, -maxLim[1], -maxLim[2]]),
    #
    node_coords[XAYAZA] + np.array([maxLim[0], maxLim[1], maxLim[2]]),
    node_coords[XBYAZA] + np.array([-maxLim[0], maxLim[1], maxLim[2]]),
    node_coords[XAYBZA] + np.array([maxLim[0], -maxLim[1], maxLim[2]]),
    node_coords[XAYAZB] + np.array([maxLim[0], maxLim[1], -maxLim[2]]),
    node_coords[XAYBZB] + np.array([maxLim[0], -maxLim[1], -maxLim[2]]),
    node_coords[XBYAZB] + np.array([-maxLim[0], maxLim[1], -maxLim[2]]),
    node_coords[XBYBZA] + np.array([-maxLim[0], -maxLim[1], maxLim[2]]),
    node_coords[XBYBZB] + np.array([-maxLim[0], -maxLim[1], -maxLim[2]])
    ))

    #masters = np.hstack(( masters,XA,YA,ZA,XB,YB,ZB,XAYA,XBYA,XAYB,XBYB,XAZA,XBZA,XAZB,XBZB,YAZA,YBZA,YAZB,YBZB,XAYAZA,XBYAZA,XAYBZA,XAYAZB,XAYBZB,XBYAZB,XBYBZA,XBYBZB ))
    if powerTes == True:
        radii = np.hstack(( radii, radii[np.hstack((XA,YA,ZA,XB,YB,ZB,XAYA,XBYA,XAYB,XBYB,XAZA,XBZA,XAZB,XBZB,YAZA,YBZA,YAZB,YBZB,XAYAZA,XBYAZA,XAYBZA,XAYAZB,XAYBZB,XBYAZB,XBYBZA,XBYBZB ))]))
    else:
        radii = []
    """
    nNds = np.asarray(nNds)
    if SHOW_PLOT:
        fig = plt.figure()
        ax = Axes3D(fig)
        #ax.auto_scale_xyz([-maxLim[0], 2*maxLim[0]], [-maxLim[1], 2*maxLim[1]], [-maxLim[2], 2*maxLim[2]])
        #ax.scatter(node_coords[:,0], node_coords[:,1], node_coords[:,2], color='r')
        ax.scatter(nNds[:,0], nNds[:,1], nNds[:,2], color='r')
        plt.show()
    """
    return nNds, mechBC_merged, mechInitC_merged, radii


def asssemble2dCircRVE(maxLim, minDist, trials, powerTes):
    dim = 2
    # lists for the model
    node_coords = []
    mechBC_merged = []
    mechInitC_merged = []

    print('Assembling 2d periodic circular RVE.')

    ########### getting nodes and radii


    #generate surface nodes
    nnodes = int(np.pi*maxLim/(minDist)/2+0.5)
    step = np.pi/nnodes
    node_coords_polar = np.zeros( (4*nnodes, dim) )    
    lim = 0.0001
    for i in range(nnodes):
        node_coords_polar[4*i,0] = i*step
        node_coords_polar[4*i+1,0] = (nnodes+i)*step
        node_coords_polar[4*i+2,0] = i*step
        node_coords_polar[4*i+3,0] = (nnodes+i)*step
        node_coords_polar[4*i+2,1] = (1.+lim)*maxLim/2.
        node_coords_polar[4*i+3,1] = (1.+lim)*maxLim/2.
        node_coords_polar[4*i,1] = (1.-lim)*maxLim/2.
        node_coords_polar[4*i+1,1] = (1.-lim)*maxLim/2.
    node_coords = np.column_stack((np.cos(node_coords_polar[:,0])*node_coords_polar[:,1],np.sin(node_coords_polar[:,0])*node_coords_polar[:,1]))
    radii = np.ones( len(node_coords_polar) )*0.4*minDist


    node_coords, polar_node_coords, radii = pointGenerators.generateParticlesSphere(maxLim, minDist*0.4, minDist, 0.8, dim, trials, node_coords, radii, allow_domain_overlap=False, periodic_distance=False)

    node_coords = np.asarray(node_coords)

    return node_coords, polar_node_coords, mechBC_merged, mechInitC_merged, radii#, masters



def assemble3DSSBeamBending (maxLim, minDist, trials, notch, loadWidth,  fracZoneWidth = 0.15, orthogonalFracZone=False, notchWidth = -1, coupled=False, node_coords_init=None, specifiedNodes=[], roughMinDistCoef=1, supportDivision=10,gapWidth=0,blank=0,powerTes=False,supportWidth=None,span=None,gradientZoneWidth=0):
    if span>maxLim[0]:
        print("Span has to be larger than XSIZE or none.")
        return

    if maxLim[0]/2 - fracZoneWidth/2 < 0:
        fracZoneWidth = maxLim[0]

    if maxLim[0]/2 - fracZoneWidth/2 - gradientZoneWidth < 0:
        fullgradientZoneWidth = gradientZoneWidth
        gradientZoneWidth = maxLim[0]/2 - fracZoneWidth/2
        roughMinDistCoef = roughMinDistCoef/fullgradientZoneWidth*gradientZoneWidth

    if span is None:
        span = maxLim[0]

    if supportWidth is None:
        supportWidth = span/40

    span -= 1e-7
    Xoverhang = maxLim[0]/2-span/2-supportWidth/2

    span_A = maxLim[0]/2-span/2
    span_B = maxLim[0]/2+span/2

    print ("span %s %s" %(span_A,span_B))

    L_sup_A = span_A-supportWidth/2
    if L_sup_A <0:
        L_sup_A = 1e-7
    L_sup_B = span_A+supportWidth/2

    R_sup_A = span_B-supportWidth/2
    R_sup_B = span_B+supportWidth/2
    if R_sup_B > maxLim[0]:
        R_sup_B = maxLim[0]-1e-7

    print("LS %s %s" %(L_sup_A,L_sup_B))
    print("RS %s %s" %(R_sup_A,R_sup_B))



    dim = 3

    if node_coords_init is None:
        node_coords = []
    else:
        node_coords = node_coords_init
    mechBC_merged = []
    mechInitC_merged = []
    govNodes = []
    govNodesMechBC = []
    rigidPlates = []

    #exact notch
    exactNotch = False

    indent = 1e-7
    notches=[]

    #notch heaight
    if not exactNotch:
        nHeight = maxLim[1]*notch - minDist/2
        #notch = nHeight / maxLim[1]


    if node_coords_init is None:
        if (len(specifiedNodes)>0):
            print ('appending specified nodes...')
            for node in specifiedNodes:
                node_coords.append((node))
            print (node_coords)

        node_coords.append( np.array([maxLim[0]/4, maxLim[1]/2, maxLim[2]/2]))
        #lineBC = np.array([0,0,0, 0,0,0,  -1,-1,-1,-1,-1,-1])
        #mBC = utilitiesMech.mechanicalBC(dim, 0, lineBC)
        #mechBC_merged.append(mBC)


        if notchWidth == -1:
            notchWidth = minDist /4
        else:
            notchWidth /= 2.0

        #generating notch points
        if (notch > 0):
            notchSide0 = []
            notchSide1 = []

            nodeA = np.array([maxLim[0]/2, maxLim[1]*notch, indent])
            nodeB = np.array([maxLim[0]/2, maxLim[1]*notch, maxLim[2]-indent])
            #pointGenerators.generateNodesLine3dRand(nodeA, nodeB, minDist, dim, node_coords, trials, catchCorners=False, equidist=False)


            if gapWidth >0:
                nx = 12
                nz = 6
                xind = gapWidth*0.2
                zind = maxLim[2]*0.1


                #pointGenerators.generate_flat_orthogonal_grid(nx*2, nz, node_coords,constants=[maxLim[0]/2-notchWidth-gapWidth*2  ,maxLim[1]-indent,zind],sizes=[gapWidth*4+notchWidth,1,maxLim[2]-indent*2-zind*2])

                oldl = len(node_coords)
                pointGenerators.generate_flat_orthogonal_grid(nx, nz, node_coords,constants=[maxLim[0]/2-notchWidth-gapWidth*1.5  ,indent,zind],sizes=[gapWidth*1.5-xind,1,maxLim[2]-indent*2-zind*2])
                newl = len(node_coords)
                gridlen=newl-oldl
                print('generated grid %d'%gridlen)
                """
                for i in range(gridlen):
                    if node_coords[oldl+i][0]>=maxLim[0]/2-notchWidth-gapWidth-1e-3:
                        mechBC = np.array([-1,-1,0,-1,-1,-1,  -1,2,-1,-1,-1,-1  ])
                        mBC = utilitiesMech.mechanicalBC(dim, oldl+i, mechBC)
                        mechBC_merged.append(mBC)
                        notchSide0.append(i+oldl)
                """

                oldl = len(node_coords)
                pointGenerators.generate_flat_orthogonal_grid(nx, nz, node_coords,constants=[maxLim[0]/2+notchWidth+xind  ,indent,zind],sizes=[gapWidth*1.5-xind,1,maxLim[2]-indent*2-zind*2])
                newl = len(node_coords)
                gridlen=newl-oldl
                print('generated grid %d'%gridlen)
                """
                for i in range(gridlen):
                    if node_coords[oldl+i][0]<=maxLim[0]/2+notchWidth+gapWidth+1e-3:
                        mechBC = np.array([-1,-1,0,-1,-1,-1,  -1,2,-1,-1,-1,-1  ])
                        mBC = utilitiesMech.mechanicalBC(dim, oldl+i, mechBC)
                        mechBC_merged.append(mBC)
                        notchSide1.append(i+oldl)
                """

            oldLen = len(node_coords)

            nodeA = np.array([maxLim[0]/2-notchWidth, maxLim[1]*notch-minDist/8, indent])
            nodeB = np.array([maxLim[0]/2-notchWidth, maxLim[1]*notch-minDist/8, maxLim[2]-indent])
            pointGenerators.generateNodesLine3dRand(nodeA, nodeB, minDist/4, dim, node_coords, trials, catchCorners=True, equidist=True)

            nTop = maxLim[1]*notch+indent-minDist/2
            if nTop < minDist:
                nTop = maxLim[1]*notch+indent

            # line from bottom to crack tip - front face
            nodeA = np.array([maxLim[0]/2-notchWidth, indent, indent])
            nodeB = np.array([maxLim[0]/2-notchWidth, nTop, indent])
            # NOTE JK here is the lmin on notch
            for i in node_coords:
                if i[1]<0:
                    print(i)
                    return

            pointGenerators.generateNodesLine3dRand(nodeA, nodeB, minDist, dim, node_coords, trials, catchCorners=True, equidist=True)
            # pointGenerators.generateNodesLine3dRand(nodeA, nodeB, minDist/4, dim, node_coords, trials, catchCorners=False, equidist=True)
            # line from bottom to crack tip - rear face
            for i in node_coords:
                if i[1]<0:
                    print(i)
                    return
            nodeA = np.array([maxLim[0]/2-notchWidth, indent, maxLim[2]-indent])
            nodeB = np.array([maxLim[0]/2-notchWidth, nTop, maxLim[2]-indent])
            # NOTE JK here is the lmin on notch
            pointGenerators.generateNodesLine3dRand(nodeA, nodeB, minDist, dim, node_coords, trials, catchCorners=True, equidist=True)
            # pointGenerators.generateNodesLine3dRand(nodeA, nodeB, minDist/4, dim, node_coords, trials, catchCorners=False, equidist=True)
            # line at the bottom of a specimen
            nodeA = np.array([maxLim[0]/2-notchWidth, indent, indent])
            nodeB = np.array([maxLim[0]/2-notchWidth, indent, maxLim[2]+indent-indent])
            # NOTE JK here is the lmin on notch
            if gapWidth ==0:
                pointGenerators.generateNodesLine3dRand(nodeA, nodeB, minDist, dim, node_coords, trials, catchCorners=False, equidist=True)

            # notch face
            nodeA = np.array([maxLim[0]/2-notchWidth, indent, indent])
            nodeB = np.array([maxLim[0]/2-notchWidth, nTop, maxLim[2]-indent])
            # NOTE JK here is the lmin on notch
            if powerTes == False:
                pointGenerators.generateNodesOrtoSurface3dRand(nodeA, nodeB, minDist/3*2, dim, node_coords, trials*2,minDistAmongNewPoints=True)
            else:
                pointGenerators.generateNodesOrtoSurface3dRand(nodeA, nodeB, minDist/3, dim, node_coords, trials*2,minDistAmongNewPoints=True)



            for i in range (oldLen, len(node_coords), 1):
                notchSide0.append(i)
                node_coords.append(np.array([maxLim[0]/2+notchWidth,
                                             node_coords[i][1],
                                             node_coords[i][2]]))
                notchSide1.append(len(node_coords)-1)



            notchL = []
            notchL.append(notchSide0)
            notchL.append(notchSide1)
            notches.append(notchL)

            leftTop = []
            oldLen = len(node_coords)
            if orthogonalFracZone:
                nodeA = np.array([maxLim[0]/2-notchWidth, maxLim[1]*notch-minDist/2, indent])
                nodeB = np.array([maxLim[0]/2-notchWidth, maxLim[1]*notch-minDist/2, maxLim[2]-indent])
            else:
                nodeA = np.array([maxLim[0]/2-notchWidth, maxLim[1]*notch, indent])
                nodeB = np.array([maxLim[0]/2-notchWidth, maxLim[1]*notch, maxLim[2]-indent])
            #pointGenerators.generateNodesLine3dRand(nodeA, nodeB, minDist, dim, node_coords, trials, catchCorners=False, equidist=True)
            for i in range (oldLen, len(node_coords), 1):
                if exactNotch:
                    leftTop.append(i)
                else:
                    notchSide0.append(i)

            rightTop = []
            oldLen = len(node_coords)
            if orthogonalFracZone:
                nodeA = np.array([maxLim[0]/2+notchWidth, maxLim[1]*notch-minDist/2, indent])
                nodeB = np.array([maxLim[0]/2+notchWidth, maxLim[1]*notch-minDist/2, maxLim[2]-indent])
            else:
                nodeA = np.array([maxLim[0]/2+notchWidth, maxLim[1]*notch, indent])
                nodeB = np.array([maxLim[0]/2+notchWidth, maxLim[1]*notch, maxLim[2]-indent])
            #pointGenerators.generateNodesLine3dRand(nodeA, nodeB, minDist, dim, node_coords, trials, catchCorners=False, equidist=True)
            for i in range (oldLen, len(node_coords), 1):
                if exactNotch:
                    rightTop.append(i)
                else:
                    notchSide1.append(i)

            notchTop = []
            notchTop.append(rightTop)
            notchTop.append(leftTop)
            if orthogonalFracZone: notches.append(notchTop)

            notch11 = []
            notch11.append(notchSide1)
            notch11.append(leftTop)
            notches.append(notch11)

            notch22 = []
            notch22.append(notchSide0)
            notch22.append(rightTop)
            notches.append(notch22)

    else:  # mark nodes on the side of the notch for remesher
        if notch > 0:
            notchWidth /= 2.0  # done the same way as in previous part just to correspond, but it's horrible!!
            tol = minDist/40
            # the following is done only for remesher performed in adaptivity
            # this can be don more elegant, but as a quick fix sufficient
            node_coords_nparray = np.asarray(node_coords)
            xs = node_coords_nparray[:, 0]
            ys = node_coords_nparray[:, 1]
            zs = node_coords_nparray[:, 2]
            xs_bounds = np.array([[maxLim[0] / 2 - notchWidth - tol, maxLim[0] / 2 - notchWidth + tol],
                                  [maxLim[0] / 2 + notchWidth - tol, maxLim[0] / 2 + notchWidth + tol]
                                  ])
            ys_bound = maxLim[1] * notch

            ymask = np.ma.masked_where(ys <= ys_bound, ys).mask

            xmask0 = np.ma.masked_inside(xs, xs_bounds[0, 0], xs_bounds[0, 1]).mask
            xmask1 = np.ma.masked_inside(xs, xs_bounds[1, 0], xs_bounds[1, 1]).mask

            ntch01 = []
            ntch01.append(np.argwhere(ymask & xmask0)[:,0].tolist())  # argwhere adds dimension
            ntch01.append(np.argwhere(ymask & xmask1)[:,0].tolist())
            notches.append(ntch01)


    if blank == 0:
        print("Generating support lines")
        if L_sup_A>2e-7:
            nodeA =  np.array([L_sup_A, indent, indent])
            nodeB =  np.array([L_sup_A, indent,  maxLim[2]-indent])
            pointGenerators.generateNodesLine3dRand(nodeA, nodeB, minDist, dim, node_coords, trials, True, True)
        nodeA =  np.array([L_sup_B, indent, indent])
        nodeB =  np.array([L_sup_B, indent,  maxLim[2]-indent])
        pointGenerators.generateNodesLine3dRand(nodeA, nodeB, minDist, dim, node_coords, trials, True, True)

        nodeA =  np.array([R_sup_A, indent, indent])
        nodeB =  np.array([R_sup_A, indent,  maxLim[2]-indent])
        pointGenerators.generateNodesLine3dRand(nodeA, nodeB, minDist, dim, node_coords, trials, True, True)
        if R_sup_B<maxLim[0]-2e-7:
            nodeA =  np.array([R_sup_B, indent, indent])
            nodeB =  np.array([R_sup_B, indent,  maxLim[2]-indent])
            pointGenerators.generateNodesLine3dRand(nodeA, nodeB, minDist, dim, node_coords, trials, True, True)


    ##################### CONSTRAINTS AND RIGID PLATES
    #rigid plate left support
    indentRP = 1e-4
    leftRigidPlateMechBC = np.array([0,0,0, 0,0,-1,  -1,-1,-1,-1,-1,-1])
    leftRigidPlate = utilitiesMech.RigidPlate(-1, 3, np.array([ L_sup_A-indentRP, L_sup_B+indentRP, -indentRP, 2*indentRP, -indentRP, maxLim[2]+indentRP  ]))
    rigidPlates.append(leftRigidPlate)
    leftSupportGN = np.array([ span_A, indent, maxLim[2]/2 ])
    govNodes.append(leftSupportGN)
    govNodesMechBC.append(utilitiesMech.mechanicalBC(dim, -1, leftRigidPlateMechBC))

    #rigid plate left support
    rightRigidPlateMechBC = np.array([-1,0,0, 0,0,-1,  -1,-1,-1,-1,-1,-1])
    rightRigidPlate = utilitiesMech.RigidPlate(-2, 3, np.array([ R_sup_A-indentRP, R_sup_B + indentRP, -indentRP, 2*indentRP, -indentRP, maxLim[2]+indentRP  ]))
    rigidPlates.append(rightRigidPlate)
    rightSupportGN = np.array([ span_B, indent, maxLim[2]/2 ])
    govNodes.append(rightSupportGN)
    govNodesMechBC.append(utilitiesMech.mechanicalBC(dim, -2, rightRigidPlateMechBC))

    if gapWidth >0 and blank == 0:
        #left rigid plate top load
        ltopRigidPlateMechBC = np.array([-1,1,0, -1,-1,-1,  -1,-1,-1,-1,-1,-1])
        ltopRigidPlate = utilitiesMech.RigidPlate(-3, 3,
        np.array([
         0.5*maxLim[0]-notchWidth-gapWidth,
        0.5*maxLim[0]-notchWidth,
        maxLim[1] - 2*indentRP,
        maxLim[1] + 2*indentRP,
        -indentRP,
        maxLim[2] + indentRP  ]))

        rigidPlates.append(ltopRigidPlate)
        govNodes.append(np.array([ 0.5*maxLim[0]-notchWidth-gapWidth/2, maxLim[1]-indent, maxLim[2]/2 ]))
        govNodesMechBC.append(utilitiesMech.mechanicalBC(dim, -3, ltopRigidPlateMechBC))

        #right rigid plate top load
        rtopRigidPlateMechBC = np.array([-1,1,0, -1,-1,-1,  -1,-1,-1,-1,-1,-1])
        rtopRigidPlate = utilitiesMech.RigidPlate(-4, 3,
        np.array([
         0.5*maxLim[0]+notchWidth,
        0.5*maxLim[0]+notchWidth+gapWidth,
        maxLim[1] - 2*indentRP,
        maxLim[1] + 2*indentRP,
        -indentRP,
        maxLim[2] + indentRP  ]))

        rigidPlates.append(rtopRigidPlate)
        govNodes.append(np.array([ 0.5*maxLim[0]+notchWidth+gapWidth/2, maxLim[1]-indent, maxLim[2]/2 ]))
        govNodesMechBC.append(utilitiesMech.mechanicalBC(dim, -4, rtopRigidPlateMechBC))

        #rigid plate left pad
        indentRP = 1e-4
        leftRigidPlateMechBC = np.array([-1,-1,0, -1,-1,-1,  -1,2,-1,-1,-1,-1])
        leftRigidPlate = utilitiesMech.RigidPlate(-5, 3, np.array([ maxLim[0]/2-notchWidth-gapWidth-1e-3, maxLim[0]/2-notchWidth+1e-5, -indentRP, indentRP, -indentRP, maxLim[2]+indentRP  ]))
        rigidPlates.append(leftRigidPlate)
        leftSupportGN = np.array([ maxLim[0]/2-notchWidth-gapWidth/2, indent, maxLim[2]/2 ])
        govNodes.append(leftSupportGN)
        govNodesMechBC.append(utilitiesMech.mechanicalBC(dim, -5, leftRigidPlateMechBC))

        #rigid plate right pad
        indentRP = 1e-4
        leftRigidPlateMechBC = np.array([-1,-1,0, -1,-1,-1,  -1,2,-1,-1,-1,-1])
        leftRigidPlate = utilitiesMech.RigidPlate(-6, 3, np.array([ maxLim[0]/2+notchWidth-1e-5, maxLim[0]/2+notchWidth+gapWidth+1e-3, -indentRP, indentRP, -indentRP, maxLim[2]+indentRP  ]))
        rigidPlates.append(leftRigidPlate)
        leftSupportGN = np.array([ maxLim[0]/2+notchWidth+gapWidth/2, indent, maxLim[2]/2 ])
        govNodes.append(leftSupportGN)
        govNodesMechBC.append(utilitiesMech.mechanicalBC(dim, -6, leftRigidPlateMechBC))

        """
        #rigid plate left pad side
        indentRP = 1e-4
        leftRigidPlateMechBC = np.array([-1,-1,-1, 0,0,-1,  -1,-1,3,-1,-1,-1])
        leftRigidPlate = utilitiesMech.RigidPlate(-7, 3, np.array([ maxLim[0]/2-notchWidth-gapWidth-1e-3, maxLim[0]/2+notchWidth+gapWidth+1e-3, indentRP, maxLim[1]-indentRP*2, -indentRP*2, indentRP/2  ]))
        rigidPlates.append(leftRigidPlate)
        leftSupportGN = np.array([ maxLim[0]/2-notchWidth-gapWidth/2, maxLim[1]/2, indent ])
        govNodes.append(leftSupportGN)
        govNodesMechBC.append(utilitiesMech.mechanicalBC(dim, -7, leftRigidPlateMechBC))

        #rigid plate right pad side
        indentRP = 1e-4
        leftRigidPlateMechBC = np.array([-1,-1,-1, 0,0,-1,  -1,-1,4,-1,-1,-1])
        leftRigidPlate = utilitiesMech.RigidPlate(-8, 3, np.array([ maxLim[0]/2-notchWidth-gapWidth-1e-3, maxLim[0]/2+notchWidth+gapWidth+1e-3, indentRP, maxLim[1]-indentRP*2, maxLim[2]-indentRP/2, maxLim[2]+indentRP/2  ]))
        rigidPlates.append(leftRigidPlate)
        leftSupportGN = np.array([ maxLim[0]/2+notchWidth+gapWidth/2, maxLim[1]/2, maxLim[2]])
        govNodes.append(leftSupportGN)
        govNodesMechBC.append(utilitiesMech.mechanicalBC(dim, -8, leftRigidPlateMechBC))
        """

        #top surf
        nodeA =  np.array([maxLim[0]/2-notchWidth-gapWidth , maxLim[1] - indent, indent])
        nodeB =  np.array([maxLim[0]/2+notchWidth+gapWidth, maxLim[1] - indent, maxLim[2] - indent])
        pointGenerators.generateNodesOrtoSurface3dRand(nodeA, nodeB, minDist, dim, node_coords, trials)
    elif blank == 0:
        #rigid plate top load
        topRigidPlateMechBC = np.array([-1,1,-1, -1,-1,-1,  -1,-1,-1,-1,-1,-1])
        topRigidPlate = utilitiesMech.RigidPlate(-3, 3,
        np.array([
         0.5*maxLim[0]-indentRP-notchWidth,
        0.5*maxLim[0]+indentRP+notchWidth,
        maxLim[1] - 2*indentRP,
        maxLim[1] + 2*indentRP,
        -indentRP,
        maxLim[2] + indentRP  ]))

        rigidPlates.append(topRigidPlate)
        govNodes.append(np.array([ maxLim[0]/2, maxLim[1]-indent, maxLim[2]/2 ]))
        govNodesMechBC.append(utilitiesMech.mechanicalBC(dim, -3, topRigidPlateMechBC))

        ############### loaded top face ###############
        lineBC = np.array([-1,-1,-1,-1,-1,-1,  -1, 1,-1,-1,-1,-1])
        nodeA =  np.array([0.5*maxLim[0], maxLim[1]-indent, indent])
        nodeB =  np.array([0.5*maxLim[0], maxLim[1]-indent, maxLim[2]-indent])
        pointGenerators.generateNodesLine3dRand(nodeA, nodeB, maxLim[2]/supportDivision, dim, node_coords, trials, False, True)


    if node_coords_init is None:
        """
        ###############generating of nodes, left horizontal support ###############
        nodeA = np.array([indent, indent, indent])
        nodeB = np.array([indent, indent, maxLim[2]-indent])
        pointGenerators.generateNodesLine3dRand(nodeA, nodeB, minDist, dim, node_coords, trials, True, True)
        """
        print('generating supports in span')
        nodeA = np.array([span_A, indent, indent])
        nodeB = np.array([span_A, indent, maxLim[2]-indent])
        pointGenerators.generateNodesLine3dRand(nodeA, nodeB, minDist, dim, node_coords, trials, True, True)



        ###############generating of nodes, right horizontal support ###############
        """
        nodeA = np.array([maxLim[0] - indent, indent, indent])
        nodeB = np.array([maxLim[0] - indent, indent, maxLim[2]-indent])
        pointGenerators.generateNodesLine3dRand(nodeA, nodeB, minDist, dim, node_coords, trials, True, True)
        """
        nodeA = np.array([span_B, indent, indent])
        nodeB = np.array([span_B, indent, maxLim[2]-indent])
        pointGenerators.generateNodesLine3dRand(nodeA, nodeB, minDist, dim, node_coords, trials, True, True)

        if gapWidth>10:
            #front surf
            nodeA =  np.array([maxLim[0]/2-notchWidth-gapWidth*2 ,  indent, indent])
            nodeB =  np.array([maxLim[0]/2+notchWidth+gapWidth*2, maxLim[1] - indent, indent])
            pointGenerators.generateNodesOrtoSurface3dRand(nodeA, nodeB, minDist, dim, node_coords, trials)

            #back surf
            nodeA =  np.array([maxLim[0]/2-notchWidth-gapWidth*2 ,  indent, maxLim[2]-indent])
            nodeB =  np.array([maxLim[0]/2+notchWidth+gapWidth*2, maxLim[1] - indent, maxLim[2]-indent])
            pointGenerators.generateNodesOrtoSurface3dRand(nodeA, nodeB, minDist, dim, node_coords, trials)

        ###############generating of nodes, front bottom line ###############
        nodeA = np.array([indent, indent, indent])
        nodeB = np.array([maxLim[0]-indent, indent, indent])
        #pointGenerators.generateNodesLine3dRand(nodeA, nodeB, minDist*.8, dim, node_coords, trials, False, False)
        ###############generating of nodes, rear bottom line ###############
        nodeA = np.array([indent, indent,  maxLim[2]-indent])
        nodeB = np.array([maxLim[0]-indent, indent,  maxLim[2]-indent])
        #pointGenerators.generateNodesLine3dRand(nodeA, nodeB, minDist*.8, dim, node_coords, trials, False, False)
        ###############generating of nodes, front top line ###############
        nodeA = np.array([indent, maxLim[1]-indent, indent])
        nodeB = np.array([maxLim[0]-indent, maxLim[1]-indent, indent])
        #pointGenerators.generateNodesLine3dRand(nodeA, nodeB, minDist*.8, dim, node_coords, trials, False, False)
        ###############generating of nodes, rear top line ###############
        nodeA = np.array([indent, maxLim[1]-indent,  maxLim[2]-indent])
        nodeB = np.array([maxLim[0]-indent, maxLim[1]-indent,  maxLim[2]-indent])
        #pointGenerators.generateNodesLine3dRand(nodeA, nodeB, minDist*.8, dim, node_coords, trials, False, False)


        #catchCorners, equidist

        oldLen = len(node_coords)
        if gapWidth>0:
            ##########################################generating of points, fracture zone
            maxLimF = np.array([
            maxLim[0]/2-notchWidth-gapWidth,#*2,
            maxLim[1] - indent,
            maxLim[2] - indent,
            maxLim[0]/2+notchWidth+gapWidth,#*2,
            maxLim[1]*notch-indent,
            indent])
        else:
            ##########################################generating of points, fracture zone
            maxLimF = np.array([
            maxLim[0]/2 - fracZoneWidth/2,
            maxLim[1] - indent,
            maxLim[2] - indent,
            maxLim[0]/2 + fracZoneWidth/2,
            maxLim[1]*notch-indent,
            indent])
        if not orthogonalFracZone and not powerTes:
            pointGenerators.generateNodesRect(maxLimF, minDist, dim, trials, node_coords, useLowBound=True)
        elif not powerTes:
            #pointGenerators.generateOrtogrid(maxLimF, minDist/2, dim, node_coords, maxLim[0]*fracZoneWidth)  #
            pointGenerators.generateOrtogrid_variable(maxLimF, minDist, node_coords, np.array([maxLim[0]*fracZoneWidth, maxLim[1]*(1-notch)-indent*2, maxLim[2]-indent*2]) )
            #(maxLim, minDist, dim, node_coords, size)
        fracZone = []
        for i in range (oldLen, len(node_coords), 1):
            fracZone.append(i)



        if (notch > 0):
            notchFrac = []
            notchFrac.append(notchSide0)
            notchFrac.append(fracZone)
            if exactNotch:
                 notches.append(notchFrac)

            notchFrac1 = []
            notchFrac1.append(notchSide1)
            notchFrac1.append(fracZone)
            if exactNotch:
                notches.append(notchFrac1)

        """
        #front surf
        nodeA =  np.array([indent , maxLim[1] - indent, indent])
        nodeB =  np.array([maxLim[0] - indent, maxLim[1] - indent, indent])
        pointGenerators.generateNodesOrtoSurface3dRand(nodeA, nodeB, minDist*roughMinDistCoef, dim, node_coords, trials)

        #back surf
        nodeA =  np.array([indent , maxLim[1] - indent, maxLim[2]-indent])
        nodeB =  np.array([maxLim[0] - indent, maxLim[1] - indent, maxLim[2]-indent])
        pointGenerators.generateNodesOrtoSurface3dRand(nodeA, nodeB, minDist*roughMinDistCoef, dim, node_coords, trials)
        """
        #top surf
        nodeA =  np.array([indent , maxLim[1] - indent, indent])
        nodeB =  np.array([maxLim[0] - indent, maxLim[1] - indent, maxLim[2] - indent])
        pointGenerators.generateNodesOrtoSurface3dRand(nodeA, nodeB, minDist*roughMinDistCoef, dim, node_coords, trials)


        if gapWidth==0:
            #bot surf
            nodeA =  np.array([indent , indent, indent])
            nodeB =  np.array([maxLim[0] - indent, indent,  maxLim[2] - indent])
            pointGenerators.generateNodesOrtoSurface3dRand(nodeA, nodeB, minDist*roughMinDistCoef, dim, node_coords, trials)

        #left face surf
        nodeA =  np.array([indent , indent, indent])
        nodeB =  np.array([indent, maxLim[1] - indent, maxLim[2] - indent])
        pointGenerators.generateNodesOrtoSurface3dRand(nodeA, nodeB, minDist*roughMinDistCoef, dim, node_coords, trials)

        #right face surf
        nodeA =  np.array([maxLim[0]-indent , indent, indent])
        nodeB =  np.array([maxLim[0]-indent, maxLim[1] - indent, maxLim[2] - indent])
        pointGenerators.generateNodesOrtoSurface3dRand(nodeA, nodeB, minDist*roughMinDistCoef, dim, node_coords, trials)





        if gapWidth>0:
            maxLimF = np.array([
            maxLim[0]/2-notchWidth-gapWidth*2,
            indent,#* notch,
            indent,
            maxLim[0]/2+notchWidth+gapWidth*2,
            maxLim[1] -indent,#* notch,
            maxLim[2] - indent
            ])
            pointGenerators.generateNodesRect(maxLimF, minDist, dim, trials, node_coords, useLowBound=True)


        else:
            ##########################################generating of points, fracture zone
            if powerTes == False:
                maxLimF = np.array([
                maxLim[0]/2 - fracZoneWidth/2,
                maxLim[1] -indent,#* notch,
                maxLim[2] - indent,
                maxLim[0]/2 + fracZoneWidth/2,
                indent,
                indent])
                pointGenerators.generateNodesRect(maxLimF, minDist, dim, trials, node_coords, useLowBound=True)
            else:
                maxLimF = np.array([
                maxLim[0]/2 + fracZoneWidth/2,
                maxLim[1] -indent,#* notch,
                maxLim[2] - indent,
                maxLim[0]/2 - fracZoneWidth/2,
                maxLim[1]*notch,
                indent])
                radii = np.zeros(len(node_coords))
                oldlen = len(node_coords)
                node_coords, radii = pointGenerators.generateParticlesRect(maxLimF, minDist*0.5, minDist, 0.8, dim, trials, np.asarray(node_coords), radii,useLowBound=True,allow_domain_overlap=True)

                maxLimF = np.array([
                maxLim[0]/2 + fracZoneWidth/2,
                maxLim[1]*notch,#* notch,
                maxLim[2] - indent,
                maxLim[0]/2 + notchWidth+minDist*0.3,
                indent,
                indent])
                oldlen = len(node_coords)
                node_coords, radii = pointGenerators.generateParticlesRect(maxLimF, minDist*0.5, minDist, 0.8, dim, trials, np.asarray(node_coords), radii,useLowBound=True,allow_domain_overlap=True)

                leftS = []
                for i in range (oldLen, len(node_coords), 1):
                    leftS.append( i)


                maxLimF = np.array([
                maxLim[0]/2 -notchWidth-minDist*0.3,
                maxLim[1]*notch,#* notch,
                maxLim[2] - indent,
                maxLim[0]/2 - fracZoneWidth/2,
                indent,
                indent])
                oldlen = len(node_coords)
                node_coords, radii = pointGenerators.generateParticlesRect(maxLimF, minDist*0.5, minDist, 0.8, dim, trials, np.asarray(node_coords), radii,useLowBound=True,allow_domain_overlap=True)

                rightS = []
                for i in range (oldLen, len(node_coords), 1):
                    rightS.append(i)

                node_coords = node_coords.tolist()



        if gapWidth>0:
            maxLimF = np.array([
            maxLim[0]/2-notchWidth-gapWidth*2,
            indent,
            indent,
            maxLim[0]/2-notchWidth-gapWidth,#*2,
            maxLim[1] -indent,#* notch,
            maxLim[2] - indent
            ])

        else:
            maxLimF = np.array([
            maxLim[0]/2-fracZoneWidth/2-gradientZoneWidth,
            indent,#* notch,
            indent,
            maxLim[0]/2-fracZoneWidth/2,
            maxLim[1] -indent,#* notch,
            maxLim[2] - indent
            ])
        print("gap ",maxLimF)
        pointGenerators.generateNodesRect(maxLimF, minDist, dim, trials, node_coords, useLowBound=True, bottomMinDist = minDist, topMinDist = minDist*roughMinDistCoef, gradienDirection=0)

        if gapWidth>0:
            maxLimF = np.array([
            maxLim[0]/2+notchWidth+gapWidth,
            indent,
            indent,
            maxLim[0]/2+notchWidth+gapWidth*2,
            maxLim[1] -indent,#* notch,
            maxLim[2] - indent
            ])

        else:
            maxLimF = np.array([
            maxLim[0]/2+fracZoneWidth/2,
            maxLim[1] -indent,#* notch,
            maxLim[2] - indent,
            maxLim[0]/2+fracZoneWidth/2+gradientZoneWidth,
            indent,#* notch,
            indent
            ])
        print("gap ",maxLimF)
        pointGenerators.generateNodesRect(maxLimF, minDist, dim, trials, node_coords, useLowBound=True, bottomMinDist = minDist*roughMinDistCoef, topMinDist = minDist, gradienDirection=0)

        if blank == 0:

            ##########################################generating of points, left support
            maxLimF = np.array([
            L_sup_B+minDist,
            supportWidth*3,
            maxLim[2]-indent,
            L_sup_A,
            indent,
            indent])
            pointGenerators.generateNodesRect(maxLimF, minDist, dim, trials, node_coords, useLowBound=True,
            topMinDist=minDist*roughMinDistCoef, bottomMinDist=minDist, minDistCenter=3)




            ##########################################generating of points, right support
            maxLimF = np.array([
            R_sup_B,
            supportWidth*3,
            maxLim[2]-indent,
            R_sup_A-minDist,
            indent,
            indent])
            pointGenerators.generateNodesRect(maxLimF, minDist, dim, trials, node_coords, useLowBound=True,
            topMinDist=minDist*roughMinDistCoef, bottomMinDist=minDist, minDistCenter=4)

            nodeA =  np.array([0.5*maxLim[0], maxLim[1]-indent, indent])
            nodeB =  np.array([0.5*maxLim[0], maxLim[1]-indent, maxLim[2]-indent])
            #pointGenerators.generateNodesLine3dRand(nodeA, nodeB, maxLim[2]/supportDivision, dim, node_coords, trials, False, True)




        print('last rect')
        ##########################################generating of points, homogeneous volume
        pointGenerators.generateNodesRect(maxLim, minDist*roughMinDistCoef, dim, trials, node_coords)

    #if coupled:
    #    notches = []

    if powerTes == True:
        newnodes = len(node_coords) - len(radii)
        for i in range(newnodes):
            radii = np.hstack((radii, 0));
    else:
        radii=None


    node_coords = np.asarray(node_coords)
    if False:
        fig = plt.figure()
        ax = Axes3D(fig)
        ax.scatter(node_coords[:,0], node_coords[:,1], node_coords[:,2])
        plt.show()


    return node_coords, mechBC_merged, mechInitC_merged, notches, govNodes, govNodesMechBC, rigidPlates,radii






def assemble3DSSBeamBendingBlank (maxLim, minDist, trials, notch, loadWidth,  fracZoneWidth = 0.15, orthogonalFracZone=False, notchWidth = -1, coupled=False, node_coords_init=None, specifiedNodes=[], roughMinDistCoef=1, supportDivision=10,gapWidth=0,blank=0,powerTes=False,supportWidth=None,span=None,gradientZoneWidth=0):
    dim = 3
    if node_coords_init is None:
        node_coords = []
    else:
        node_coords = node_coords_init
    mechBC_merged = []
    mechInitC_merged = []
    govNodes = []
    govNodesMechBC = []
    rigidPlates = []
    radii=[]
    #exact notch
    exactNotch = False

    indent = 1e-7
    notches=[]

    if powerTes == True:
        node_coords = np.zeros((0,dim))
        radii = np.zeros(len(node_coords))
        maxLimF= np.copy(maxLim)
        maxLimF[0]=maxLim[0]/4
        #node_coords, radii = pointGenerators.generateParticlesRect(maxLimF, 0.004, 0.010, 0.8, dim, trials, node_coords, radii, allow_domain_overlap = True, periodic_distance=True)
        #node_coords[:,0]+=maxLim[0]/4*1.5

        maxLimF= np.zeros((6))
        maxLimF[0] = maxLim[0]/3*2
        maxLimF[1] = maxLim[1]
        maxLimF[2] = maxLim[2]
        maxLimF[3] = maxLim[0]/3
        maxLimF[4] = 0
        maxLimF[5] = 0
        #TPB
        node_coords, radii = pointGenerators.generateParticlesRect(maxLimF, 0.012, 0.016, 0.9, dim, trials, 
                                                                   node_coords, radii, allow_domain_overlap = True, 
                                                                   periodic_distance=False,useLowBound=True)
        maxLimF= np.zeros((6))
        maxLimF[0] = maxLim[0]/3
        maxLimF[1] = maxLim[1]
        maxLimF[2] = maxLim[2]
        maxLimF[3] = 0
        maxLimF[4] = 0
        maxLimF[5] = 0
        #TPB
        node_coords, radii = pointGenerators.generateParticlesRect(maxLimF, 0.014, 0.016, 0.9, dim, trials, 
                                                                   node_coords, radii, allow_domain_overlap = True, 
                                                                   periodic_distance=False,useLowBound=True)
        maxLimF= np.zeros((6))
        maxLimF[0] = maxLim[0]
        maxLimF[1] = maxLim[1]
        maxLimF[2] = maxLim[2]
        maxLimF[3] = maxLim[0]/3*2
        maxLimF[4] = 0
        maxLimF[5] = 0
        #TPB
        node_coords, radii = pointGenerators.generateParticlesRect(maxLimF, 0.014, 0.016, 0.9, dim, trials, 
                                                                   node_coords, radii, allow_domain_overlap = True, 
                                                                   periodic_distance=False,useLowBound=True)
        #CYLINDER
        #node_coords, radii = pointGenerators.generateParticlesRect(maxLim, 0.012, 0.024, 0.95, dim, trials, node_coords, radii, allow_domain_overlap = True, periodic_distance=False)
        np.savetxt('radii.txt', np.asarray(radii), fmt='%f') 
    else:
        pointGenerators.generateNodesRect(maxLim, minDist, dim, trials, node_coords)



    return node_coords, mechBC_merged, mechInitC_merged, notches, govNodes, govNodesMechBC, rigidPlates,radii







def assemble3DBalbet(maxLim, minDist, trials, powerTes, shotRadius=0.02, shotGradientRadius=0.01, roughMinDistCoef=2,coupled=False, node_coords_init=None):
    dim = 3
    #lists for the model
    if node_coords_init is None:
        node_coords = []
    else:
        node_coords = node_coords_init
    mechBC_merged = []
    mechInitC_merged = []
    govNodes = []
    govNodesMechBC = []
    rigidPlates = []


    indent = 1e-5

    minDist *= roughMinDistCoef



    indent = 1e-7



    if node_coords_init is None:
        if powerTes == False:
            radii = []

            try:
                shots = []
                with open('shots.txt', 'r') as file:
                    for line in file:
                        x, y, z = map(float, line.split())
                        shots.append(np.array([x, y, z]))

            except Exception as e:
                print(f"No shots.txt file. Usin a single shot in the middle")
                shots = [maxLim/2]
                shots[0][2] = indent
                shots=np.asarray(shots)
            for shot in shots:
                shot[2]=indent
                node_coords.append(shot)
                mechBC_merged.append(utilitiesMech.mechanicalBC(dim, len(node_coords)-1, np.array([-1,-1, 1,   -1,-1,-1,    -1,-1,-1,  -1,-1,-1]) ))
                pointGenerators.generateNodesOrtoCilinder3dRand(shot, shotRadius, maxLim[2], 2, minDist/roughMinDistCoef, node_coords, trials,gradient_radius=shotGradientRadius,maxMinDist=minDist)

            oldLen = len(node_coords)
            nodeA = np.array([indent, indent,  maxLim[2]-indent])
            nodeB = np.array([maxLim[0]-indent, indent,  maxLim[2]-indent])
            pointGenerators.generateNodesLine3dRand(nodeA, nodeB, minDist, dim, node_coords, trials, False, False)
            nodeA = np.array([indent, maxLim[1]-indent,  maxLim[2]-indent])
            nodeB = np.array([maxLim[0]-indent, maxLim[1]-indent,  maxLim[2]-indent])
            pointGenerators.generateNodesLine3dRand(nodeA, nodeB, minDist, dim, node_coords, trials, False, False)
            nodeA = np.array([indent, indent, maxLim[2]-indent])
            nodeB = np.array([indent, maxLim[1]-indent, maxLim[2]-indent])
            pointGenerators.generateNodesLine3dRand(nodeA, nodeB, minDist, dim, node_coords, trials, False, False)
            nodeA = np.array([maxLim[0]-indent, indent, maxLim[2]-indent])
            nodeB = np.array([maxLim[0]-indent, maxLim[1]-indent, maxLim[2]-indent])
            pointGenerators.generateNodesLine3dRand(nodeA, nodeB, minDist, dim, node_coords, trials, False, False)

            newNodes = len(node_coords)-oldLen
            rpIdcs = []
            for i in range (newNodes):
                rpIdcs.append(oldLen + i)

            # boundary conditions
            topRigidPlate = utilitiesMech.RigidPlate(-1, 3, None, directIdcs = True)
            topRigidPlate.setDirectNodes(rpIdcs)
            topRigidPlateMechBC = np.array([0,0,0, 0,0,0,  -1,-1,-1,-1,-1,-1])
            rigidPlates.append(topRigidPlate)
            govNodes.append(maxLim/2)
            govNodesMechBC.append(utilitiesMech.mechanicalBC(dim, -1, topRigidPlateMechBC))



            ###############generating of nodes, front bottom line ###############
            nodeA = np.array([indent, indent, indent])
            nodeB = np.array([maxLim[0]-indent, indent, indent])
            pointGenerators.generateNodesLine3dRand(nodeA, nodeB, minDist, dim, node_coords, trials, True, False)
            ###############generating of nodes, front top line ###############
            nodeA = np.array([indent, maxLim[1]-indent, indent])
            nodeB = np.array([maxLim[0]-indent, maxLim[1]-indent, indent])
            pointGenerators.generateNodesLine3dRand(nodeA, nodeB, minDist, dim, node_coords, trials, True, False)


            ###############generating of nodes, left face lines ###############
            nodeA = np.array([indent, indent, indent])
            nodeB = np.array([indent, indent, maxLim[2]-indent])
            pointGenerators.generateNodesLine3dRand(nodeA, nodeB, minDist, dim, node_coords, trials, False, False)
            nodeA = np.array([indent, maxLim[1]-indent, indent])
            nodeB = np.array([indent, maxLim[1]-indent, maxLim[2]-indent])
            pointGenerators.generateNodesLine3dRand(nodeA, nodeB, minDist, dim, node_coords, trials, False, False)
            nodeA = np.array([indent, indent, indent])
            nodeB = np.array([indent, maxLim[1]-indent, indent])
            pointGenerators.generateNodesLine3dRand(nodeA, nodeB, minDist, dim, node_coords, trials, False, False)


            ###############generating of nodes, right face lines ###############
            nodeA = np.array([maxLim[0]-indent, indent, indent])
            nodeB = np.array([maxLim[0]-indent, indent, maxLim[2]-indent])
            pointGenerators.generateNodesLine3dRand(nodeA, nodeB, minDist, dim, node_coords, trials, False, False)
            nodeA = np.array([maxLim[0]-indent, maxLim[1]-indent, indent])
            nodeB = np.array([maxLim[0]-indent, maxLim[1]-indent, maxLim[2]-indent])
            pointGenerators.generateNodesLine3dRand(nodeA, nodeB, minDist, dim, node_coords, trials, False, False)
            nodeA = np.array([maxLim[0]-indent, indent, indent])
            nodeB = np.array([maxLim[0]-indent, maxLim[1]-indent, indent])
            pointGenerators.generateNodesLine3dRand(nodeA, nodeB, minDist, dim, node_coords, trials, False, False)



            ############### loaded top surf ###############
            nodeA =  np.array([indent, maxLim[1]-indent, indent])
            nodeB =  np.array([maxLim[0]-indent, maxLim[1]-indent, maxLim[2]-indent])
            pointGenerators.generateNodesOrtoSurface3dRand(nodeA, nodeB, minDist, dim, node_coords, trials)

            #front surf
            nodeA =  np.array([indent ,  indent, indent])
            nodeB =  np.array([maxLim[0] - indent, maxLim[1] - indent, indent])
            pointGenerators.generateNodesOrtoSurface3dRand(nodeA, nodeB, minDist, dim, node_coords, trials)

            #back surf
            nodeA =  np.array([indent , maxLim[1] - indent, maxLim[2]-indent])
            nodeB =  np.array([maxLim[0] - indent, maxLim[1] - indent, maxLim[2]-indent])
            pointGenerators.generateNodesOrtoSurface3dRand(nodeA, nodeB, minDist, dim, node_coords, trials)

            #top surf
            nodeA =  np.array([indent , maxLim[1] - indent, indent])
            nodeB =  np.array([maxLim[0] - indent, maxLim[1] - indent, maxLim[2] - indent])
            pointGenerators.generateNodesOrtoSurface3dRand(nodeA, nodeB, minDist, dim, node_coords, trials)

            #bot surf
            nodeA =  np.array([indent , indent, indent])
            nodeB =  np.array([maxLim[0] - indent, indent,  maxLim[2] - indent])
            pointGenerators.generateNodesOrtoSurface3dRand(nodeA, nodeB, minDist, dim, node_coords, trials)

            #left face surf
            nodeA =  np.array([indent , indent, indent])
            nodeB =  np.array([indent, maxLim[1] - indent, maxLim[2] - indent])
            pointGenerators.generateNodesOrtoSurface3dRand(nodeA, nodeB, minDist, dim, node_coords, trials)

            #right face surf
            nodeA =  np.array([maxLim[0]-indent , indent, indent])
            nodeB =  np.array([maxLim[0]-indent, maxLim[1] - indent, maxLim[2] - indent])
            pointGenerators.generateNodesOrtoSurface3dRand(nodeA, nodeB, minDist, dim, node_coords, trials)





            ##########################################generating of points, homogeneous volume
            pointGenerators.generateNodesRect(maxLim, minDist, dim, trials, node_coords)

        if powerTes == True:
            node_coords = np.zeros((0,dim))
            radii = np.zeros(len(node_coords))

            """
            mechBC = np.array([-1,0,0,-1,-1,-1,    -1,-1,-1,-1,-1,-1])
            mBC = utilitiesMech.mechanicalBC(dim, 0, mechBC)
            mechBC_merged.append(mBC)
            """

            node_coords = np.vstack((node_coords,   np.array([maxLim[0]/2, maxLim[1]/2, maxLim[2]/2])  ))
            radii = np.hstack((radii, minDist*0.4));


            #front surf
            nodeA =  np.array([indent ,  indent, indent])
            nodeB =  np.array([maxLim[0] - indent, maxLim[1] - indent, indent])
            node_coords, radii = pointGenerators.generateParticlesOrtoSurface(nodeA, nodeB, minDist*0.4, minDist, 0.8, dim, trials, node_coords, radii, allow_domain_overlap=True)

            #back surf
            nodeA =  np.array([indent ,  indent, maxLim[2] - indent])
            nodeB =  np.array([maxLim[0] - indent, maxLim[1] - indent, maxLim[2] -indent])
            node_coords, radii = pointGenerators.generateParticlesOrtoSurface(nodeA, nodeB, minDist*0.4, minDist, 0.8, dim, trials, node_coords, radii, allow_domain_overlap=True)

            #top surf
            nodeA =  np.array([indent , maxLim[1] - indent, indent])
            nodeB =  np.array([maxLim[0] - indent, maxLim[1] - indent, maxLim[2] - indent])
            node_coords, radii = pointGenerators.generateParticlesOrtoSurface(nodeA, nodeB, minDist*0.4, minDist, 0.8, dim, trials, node_coords, radii, allow_domain_overlap=True)

            #bot surf
            nodeA =  np.array([indent , indent, indent])
            nodeB =  np.array([maxLim[0] - indent, indent,  maxLim[2] - indent])
            node_coords, radii = pointGenerators.generateParticlesOrtoSurface(nodeA, nodeB, minDist*0.4, minDist, 0.8, dim, trials, node_coords, radii, allow_domain_overlap=True)

            #left face surf
            nodeA =  np.array([indent , indent, indent])
            nodeB =  np.array([indent, maxLim[1] - indent, maxLim[2] - indent])
            node_coords, radii = pointGenerators.generateParticlesOrtoSurface(nodeA, nodeB, minDist*0.4, minDist, 0.8, dim, trials, node_coords, radii, allow_domain_overlap=True)

            #right face surf
            nodeA =  np.array([maxLim[0]-indent , indent, indent])
            nodeB =  np.array([maxLim[0]-indent, maxLim[1] - indent, maxLim[2] - indent])
            node_coords, radii = pointGenerators.generateParticlesOrtoSurface(nodeA, nodeB, minDist*0.4, minDist, 0.8, dim, trials, node_coords, radii, allow_domain_overlap=True)

            # volume
            node_coords, radii = pointGenerators.generateParticlesRect(maxLim, minDist*0.4, minDist, 0.8, dim, trials, node_coords, radii, allow_domain_overlap = False, periodic_distance=True)



        node_coords = np.asarray(node_coords)
        """
        fig = plt.figure()
        ax = Axes3D(fig)
        ax.scatter(node_coords[:,0],node_coords[:,1],node_coords[:,2])
        plt.show()
        """


    #if coupled:
    #    notches = []
    return node_coords, mechBC_merged, mechInitC_merged, govNodes, govNodesMechBC, rigidPlates, radii











def assemble3Dcube(maxLim, minDist, trials, powerTes, coupled=False, node_coords_init=None, periodic = False):
    dim = 3
    #lists for the model
    if node_coords_init is None:
        node_coords = []
    else:
        node_coords = node_coords_init
    mechBC_merged = []
    mechInitC_merged = []
    govNodes = []
    govNodesMechBC = []
    rigidPlates = []


    indent = 1e-6


    ##################### CONSTRAINTS AND RIGID PLATES
    #rigid plate left support
    indentRP = indent
    leftRigidPlateMechBC = np.array([0,-1,-1, 0,0,0,  -1,-1,-1,-1,-1,-1])
    leftRigidPlate = utilitiesMech.RigidPlate(-1, 3,
    np.array([ -indent,  indent,
     -indent, maxLim[1],
     -indent, maxLim[2]  ]))
    rigidPlates.append(leftRigidPlate)
    govNodes.append(np.array([ 0, maxLim[1]/2, maxLim[2]/2]))
    govNodesMechBC.append(utilitiesMech.mechanicalBC(dim, -1, leftRigidPlateMechBC))
    #rigid plate right support
    rightRigidPlateMechBC = np.array([1,-1,-1, 0,0,0,  -1,-1,-1, -1,-1,-1])
    rightRigidPlate = utilitiesMech.RigidPlate(-2, 3,np.array([
    maxLim[0]-indent,  maxLim[0]+indent,
     -indent, maxLim[1],
     -indent, maxLim[2] ]))
    rigidPlates.append(rightRigidPlate)
    govNodes.append(np.array([ maxLim[0], maxLim[1]/2, maxLim[2]/2 ]))
    govNodesMechBC.append(utilitiesMech.mechanicalBC(dim, -2, rightRigidPlateMechBC))

    indent = 1e-7



    if node_coords_init is None:
        if powerTes == False:
            node_coords.append((  np.array([maxLim[0]/2, maxLim[1]/2, maxLim[2]/2])  ))
            radii = []
            """
            mechBC = np.array([-1,0,0,   -1,-1,-1,    -1,-1,-1,-1,-1,-1])
            mBC = utilitiesMech.mechanicalBC(dim, 0, mechBC)
            mechBC_merged.append(mBC)
            """

            ###############generating of nodes, front bottom line ###############
            nodeA = np.array([indent, indent, indent])
            nodeB = np.array([maxLim[0]-indent, indent, indent])
            pointGenerators.generateNodesLine3dRand(nodeA, nodeB, minDist, dim, node_coords, trials, True, False)
            ###############generating of nodes, rear bottom line ###############
            nodeA = np.array([indent, indent,  maxLim[2]-indent])
            nodeB = np.array([maxLim[0]-indent, indent,  maxLim[2]-indent])
            pointGenerators.generateNodesLine3dRand(nodeA, nodeB, minDist, dim, node_coords, trials, True, False)
            ###############generating of nodes, front top line ###############
            nodeA = np.array([indent, maxLim[1]-indent, indent])
            nodeB = np.array([maxLim[0]-indent, maxLim[1]-indent, indent])
            pointGenerators.generateNodesLine3dRand(nodeA, nodeB, minDist, dim, node_coords, trials, True, False)
            ###############generating of nodes, rear top line ###############
            nodeA = np.array([indent, maxLim[1]-indent,  maxLim[2]-indent])
            nodeB = np.array([maxLim[0]-indent, maxLim[1]-indent,  maxLim[2]-indent])
            pointGenerators.generateNodesLine3dRand(nodeA, nodeB, minDist, dim, node_coords, trials, True, False)

            ###############generating of nodes, left face lines ###############
            nodeA = np.array([indent, indent, indent])
            nodeB = np.array([indent, indent, maxLim[2]-indent])
            pointGenerators.generateNodesLine3dRand(nodeA, nodeB, minDist, dim, node_coords, trials, False, False)
            nodeA = np.array([indent, maxLim[1]-indent, indent])
            nodeB = np.array([indent, maxLim[1]-indent, maxLim[2]-indent])
            pointGenerators.generateNodesLine3dRand(nodeA, nodeB, minDist, dim, node_coords, trials, False, False)
            nodeA = np.array([indent, indent, indent])
            nodeB = np.array([indent, maxLim[1]-indent, indent])
            pointGenerators.generateNodesLine3dRand(nodeA, nodeB, minDist, dim, node_coords, trials, False, False)
            nodeA = np.array([indent, indent, maxLim[2]-indent])
            nodeB = np.array([indent, maxLim[1]-indent, maxLim[2]-indent])
            pointGenerators.generateNodesLine3dRand(nodeA, nodeB, minDist, dim, node_coords, trials, False, False)

            ###############generating of nodes, right face lines ###############
            nodeA = np.array([maxLim[0]-indent, indent, indent])
            nodeB = np.array([maxLim[0]-indent, indent, maxLim[2]-indent])
            pointGenerators.generateNodesLine3dRand(nodeA, nodeB, minDist, dim, node_coords, trials, False, False)
            nodeA = np.array([maxLim[0]-indent, maxLim[1]-indent, indent])
            nodeB = np.array([maxLim[0]-indent, maxLim[1]-indent, maxLim[2]-indent])
            pointGenerators.generateNodesLine3dRand(nodeA, nodeB, minDist, dim, node_coords, trials, False, False)
            nodeA = np.array([maxLim[0]-indent, indent, indent])
            nodeB = np.array([maxLim[0]-indent, maxLim[1]-indent, indent])
            pointGenerators.generateNodesLine3dRand(nodeA, nodeB, minDist, dim, node_coords, trials, False, False)
            nodeA = np.array([maxLim[0]-indent, indent, maxLim[2]-indent])
            nodeB = np.array([maxLim[0]-indent, maxLim[1]-indent, maxLim[2]-indent])
            pointGenerators.generateNodesLine3dRand(nodeA, nodeB, minDist, dim, node_coords, trials, False, False)

            ############### loaded top surf ###############
            nodeA =  np.array([indent, maxLim[1]-indent, indent])
            nodeB =  np.array([maxLim[0]-indent, maxLim[1]-indent, maxLim[2]-indent])
            pointGenerators.generateNodesOrtoSurface3dRand(nodeA, nodeB, minDist, dim, node_coords, trials)

            #front surf
            nodeA =  np.array([indent ,  indent, indent])
            nodeB =  np.array([maxLim[0] - indent, maxLim[1] - indent, indent])
            pointGenerators.generateNodesOrtoSurface3dRand(nodeA, nodeB, minDist, dim, node_coords, trials)

            #back surf
            nodeA =  np.array([indent , maxLim[1] - indent, maxLim[2]-indent])
            nodeB =  np.array([maxLim[0] - indent, maxLim[1] - indent, maxLim[2]-indent])
            pointGenerators.generateNodesOrtoSurface3dRand(nodeA, nodeB, minDist, dim, node_coords, trials)

            #top surf
            nodeA =  np.array([indent , maxLim[1] - indent, indent])
            nodeB =  np.array([maxLim[0] - indent, maxLim[1] - indent, maxLim[2] - indent])
            pointGenerators.generateNodesOrtoSurface3dRand(nodeA, nodeB, minDist, dim, node_coords, trials)

            #bot surf
            nodeA =  np.array([indent , indent, indent])
            nodeB =  np.array([maxLim[0] - indent, indent,  maxLim[2] - indent])
            pointGenerators.generateNodesOrtoSurface3dRand(nodeA, nodeB, minDist, dim, node_coords, trials)

            #left face surf
            nodeA =  np.array([indent , indent, indent])
            nodeB =  np.array([indent, maxLim[1] - indent, maxLim[2] - indent])
            pointGenerators.generateNodesOrtoSurface3dRand(nodeA, nodeB, minDist, dim, node_coords, trials)

            #right face surf
            nodeA =  np.array([maxLim[0]-indent , indent, indent])
            nodeB =  np.array([maxLim[0]-indent, maxLim[1] - indent, maxLim[2] - indent])
            pointGenerators.generateNodesOrtoSurface3dRand(nodeA, nodeB, minDist, dim, node_coords, trials)





            ##########################################generating of points, homogeneous volume
            pointGenerators.generateNodesRect(maxLim, minDist, dim, trials, node_coords)

        if powerTes == True:

            minDistX = 0.6*minDist
            indentX = 0.25*minDist

            node_coords = []

            #front surf
            nodeA =  np.array([indentX ,  indentX, indent])
            nodeB =  np.array([maxLim[0] - indentX, maxLim[1] - indentX, indent])
            na = len(node_coords)
            pointGenerators.generateNodesOrtoSurface3dRand(nodeA, nodeB, minDistX, dim, node_coords, trials)
            na = len(node_coords)-na

            #back surf
            if periodic:
                node_coords2 = np.copy(node_coords[-na:])
                node_coords2[:,2] += maxLim[2] - 2*indent
                node_coords =  node_coords + node_coords2.tolist()
            else:
                nodeA =  np.array([indentX ,  indentX, maxLim[2] - indent])
                nodeB =  np.array([maxLim[0] - indentX, maxLim[1] - indentX, maxLim[2] -indent])
                pointGenerators.generateNodesOrtoSurface3dRand(nodeA, nodeB, minDistX, dim, node_coords, trials)

            #bot surf
            nodeA =  np.array([indentX , indent, indentX])
            nodeB =  np.array([maxLim[0] - indentX, indent,  maxLim[2] - indentX])
            na = len(node_coords)
            pointGenerators.generateNodesOrtoSurface3dRand(nodeA, nodeB, minDistX, dim, node_coords, trials)
            na = len(node_coords)-na

            #top surf
            if periodic:
                node_coords2 = np.copy(node_coords[-na:])
                node_coords2[:,1] += maxLim[1] - 2*indent
                node_coords =  node_coords + node_coords2.tolist()
            else:
                nodeA =  np.array([indentX , maxLim[1] - indent, indentX])
                nodeB =  np.array([maxLim[0] - indentX, maxLim[1] - indent, maxLim[2] - indentX])
                pointGenerators.generateNodesOrtoSurface3dRand(nodeA, nodeB, minDistX, dim, node_coords, trials)

            #left face surf
            nodeA =  np.array([indent , indentX, indentX])
            nodeB =  np.array([indent, maxLim[1] - indentX, maxLim[2] - indentX])
            na = len(node_coords)
            pointGenerators.generateNodesOrtoSurface3dRand(nodeA, nodeB, minDistX, dim, node_coords, trials)
            na = len(node_coords)-na

            #right face surf
            if periodic:
                node_coords2 = np.copy(node_coords[-na:])
                node_coords2[:,0] += maxLim[0] - 2*indent
                node_coords =  node_coords + node_coords2.tolist()
            else:
                nodeA =  np.array([maxLim[0]-indent , indentX, indentX])
                nodeB =  np.array([maxLim[0]-indent, maxLim[1] - indentX, maxLim[2] - indentX])
                pointGenerators.generateNodesOrtoSurface3dRand(nodeA, nodeB, minDistX, dim, node_coords, trials)

            radii = np.ones(len(node_coords))*minDistX/2.
            node_coords = np.asarray(node_coords)
            # volume
            node_coords, radii = pointGenerators.generateParticlesRect(maxLim, minDist*0.4, minDist, 0.8, dim, trials, node_coords, radii, allow_domain_overlap = False, periodic_distance=True)



        node_coords = np.asarray(node_coords)
        """
        fig = plt.figure()
        ax = Axes3D(fig)
        ax.scatter(node_coords[:,0],node_coords[:,1],node_coords[:,2])
        plt.show()
        """


    #if coupled:
    #    notches = []
    return node_coords, mechBC_merged, mechInitC_merged, govNodes, govNodesMechBC, rigidPlates, radii




def assemble3DReinhardtTension (maxLim, minDist, trials, fracZoneWidth = 0.15,roughCoef=1,gapWidth=0):
    notch = 1
    dim = 3
    #lists for the model
    node_coords = []
    mechBC_merged = []
    mechInitC_merged = []
    govNodes = []
    govNodesMechBC = []
    rigidPlates = []
    notches=[]

    indent = 1e-5
    #### nodes for gauges
    node_coords.append( np.array([indent, maxLim[1]/2, maxLim[2]/2]))
    node_coords.append( np.array([maxLim[0]-indent, maxLim[1]/2, maxLim[2]/2]))

    #print (maxlim)

    notchWidth = 0.005/2
    notchDepth = 0.005
    #generating bottom notch
    if (notch > 0):
        notchSide0 = []
        oldLen = len(node_coords)
        nodeA = np.array([maxLim[0]/2-notchWidth, indent, indent])
        nodeB = np.array([maxLim[0]/2-notchWidth, indent, maxLim[2]-indent])
        pointGenerators.generateNodesLine3dRand(nodeA, nodeB, minDist*0.8, dim, node_coords, trials, catchCorners=True, equidist=True)
        nodeA = np.array([maxLim[0]/2-notchWidth, notchDepth, indent])
        nodeB = np.array([maxLim[0]/2-notchWidth, notchDepth, maxLim[2]-indent])
        pointGenerators.generateNodesLine3dRand(nodeA, nodeB, minDist*0.8, dim, node_coords, trials, catchCorners=True, equidist=True)
        nodeA = np.array([maxLim[0]/2-notchWidth, indent, indent])
        nodeB = np.array([maxLim[0]/2-notchWidth, notchDepth, indent])
        pointGenerators.generateNodesLine3dRand(nodeA, nodeB, minDist*0.8, dim, node_coords, trials, catchCorners=False, equidist=True)
        nodeA = np.array([maxLim[0]/2-notchWidth, indent, maxLim[2]-indent])
        nodeB = np.array([maxLim[0]/2-notchWidth, notchDepth, maxLim[2]-indent])
        pointGenerators.generateNodesLine3dRand(nodeA, nodeB, minDist*0.8, dim, node_coords, trials, catchCorners=False, equidist=True)

        nodeA = np.array([maxLim[0]/2-notchWidth, indent, indent])
        nodeB = np.array([maxLim[0]/2-notchWidth, notchDepth, maxLim[2]-indent])
        pointGenerators.generateNodesOrtoSurface3dRand(nodeA, nodeB, minDist*0.8, dim, node_coords, trials,minDistAmongNewPoints=True)
        for i in range (oldLen, len(node_coords), 1):
            notchSide0.append(i)


        notchSide1 = []
        oldLen = len(node_coords)
        nodeA = np.array([maxLim[0]/2+notchWidth, indent, indent])
        nodeB = np.array([maxLim[0]/2+notchWidth, indent, maxLim[2]-indent])
        pointGenerators.generateNodesLine3dRand(nodeA, nodeB, minDist*0.8, dim, node_coords, trials, catchCorners=True, equidist=True)
        nodeA = np.array([maxLim[0]/2+notchWidth, notchDepth, indent])
        nodeB = np.array([maxLim[0]/2+notchWidth, notchDepth, maxLim[2]-indent])
        pointGenerators.generateNodesLine3dRand(nodeA, nodeB, minDist*0.8, dim, node_coords, trials, catchCorners=True, equidist=True)
        nodeA = np.array([maxLim[0]/2+notchWidth, indent, indent])
        nodeB = np.array([maxLim[0]/2+notchWidth, notchDepth, indent])
        pointGenerators.generateNodesLine3dRand(nodeA, nodeB, minDist*0.8, dim, node_coords, trials, catchCorners=False, equidist=True)
        nodeA = np.array([maxLim[0]/2+notchWidth, indent, maxLim[2]-indent])
        nodeB = np.array([maxLim[0]/2+notchWidth, notchDepth, maxLim[2]-indent])
        pointGenerators.generateNodesLine3dRand(nodeA, nodeB, minDist*0.8, dim, node_coords, trials, catchCorners=False, equidist=True)

        nodeA = np.array([maxLim[0]/2+notchWidth, indent, indent])
        nodeB = np.array([maxLim[0]/2+notchWidth, notchDepth, maxLim[2]-indent])
        pointGenerators.generateNodesOrtoSurface3dRand(nodeA, nodeB, minDist*0.8, dim, node_coords, trials, minDistAmongNewPoints= True)
        for i in range (oldLen, len(node_coords), 1):
            notchSide1.append(i)





        notchSide2 = []
        oldLen = len(node_coords)
        nodeA = np.array([maxLim[0]/2-notchWidth, maxLim[1]-indent, indent])
        nodeB = np.array([maxLim[0]/2-notchWidth, maxLim[1]-indent, maxLim[2]-indent])
        pointGenerators.generateNodesLine3dRand(nodeA, nodeB, minDist*0.8, dim, node_coords, trials, catchCorners=True, equidist=True)
        nodeA = np.array([maxLim[0]/2-notchWidth, maxLim[1]-notchDepth, indent])
        nodeB = np.array([maxLim[0]/2-notchWidth, maxLim[1]-notchDepth, maxLim[2]-indent])
        pointGenerators.generateNodesLine3dRand(nodeA, nodeB, minDist*0.8, dim, node_coords, trials, catchCorners=True, equidist=True)
        nodeA = np.array([maxLim[0]/2-notchWidth, maxLim[1]-indent, indent])
        nodeB = np.array([maxLim[0]/2-notchWidth, maxLim[1]-notchDepth, indent])
        pointGenerators.generateNodesLine3dRand(nodeA, nodeB, minDist*0.8, dim, node_coords, trials, catchCorners=False, equidist=True)
        nodeA = np.array([maxLim[0]/2-notchWidth, maxLim[1]-indent, maxLim[2]-indent])
        nodeB = np.array([maxLim[0]/2-notchWidth, maxLim[1]-notchDepth, maxLim[2]-indent])
        pointGenerators.generateNodesLine3dRand(nodeA, nodeB, minDist*0.8, dim, node_coords, trials, catchCorners=False, equidist=True)

        nodeA = np.array([maxLim[0]/2-notchWidth, maxLim[1]-indent, indent])
        nodeB = np.array([maxLim[0]/2-notchWidth, maxLim[1]-notchDepth, maxLim[2]-indent])
        pointGenerators.generateNodesOrtoSurface3dRand(nodeA, nodeB, minDist*0.8, dim, node_coords, trials,minDistAmongNewPoints=True)
        for i in range (oldLen, len(node_coords), 1):
            notchSide2.append(i)


        notchSide3 = []
        oldLen = len(node_coords)
        nodeA = np.array([maxLim[0]/2+notchWidth, maxLim[1]-indent, indent])
        nodeB = np.array([maxLim[0]/2+notchWidth, maxLim[1]-indent, maxLim[2]-indent])
        pointGenerators.generateNodesLine3dRand(nodeA, nodeB, minDist*0.8, dim, node_coords, trials, catchCorners=True, equidist=True)
        nodeA = np.array([maxLim[0]/2+notchWidth, maxLim[1]-notchDepth, indent])
        nodeB = np.array([maxLim[0]/2+notchWidth, maxLim[1]-notchDepth, maxLim[2]-indent])
        pointGenerators.generateNodesLine3dRand(nodeA, nodeB, minDist*0.8, dim, node_coords, trials, catchCorners=True, equidist=True)
        nodeA = np.array([maxLim[0]/2+notchWidth, maxLim[1]-indent, indent])
        nodeB = np.array([maxLim[0]/2+notchWidth, maxLim[1]-notchDepth, indent])
        pointGenerators.generateNodesLine3dRand(nodeA, nodeB, minDist*0.8, dim, node_coords, trials, catchCorners=False, equidist=True)
        nodeA = np.array([maxLim[0]/2+notchWidth, maxLim[1]-indent, maxLim[2]-indent])
        nodeB = np.array([maxLim[0]/2+notchWidth, maxLim[1]-notchDepth, maxLim[2]-indent])
        pointGenerators.generateNodesLine3dRand(nodeA, nodeB, minDist*0.8, dim, node_coords, trials, catchCorners=False, equidist=True)

        nodeA = np.array([maxLim[0]/2+notchWidth, maxLim[1]-indent, indent])
        nodeB = np.array([maxLim[0]/2+notchWidth, maxLim[1]-notchDepth, maxLim[2]-indent])
        pointGenerators.generateNodesOrtoSurface3dRand(nodeA, nodeB, minDist*0.8, dim, node_coords, trials, minDistAmongNewPoints= True)
        for i in range (oldLen, len(node_coords), 1):
            notchSide3.append(i)



        notchL = []
        notchL.append(notchSide0)
        notchL.append(notchSide1)
        notches.append(notchL)
        notchP = []
        notchP.append(notchSide2)
        notchP.append(notchSide3)
        notches.append(notchP)

    if gapWidth >0:
        nx = 12
        nz = 6
        xind = gapWidth*0.2
        zind = maxLim[2]*0.1
        indentRP=indent*5


        pointGenerators.generate_flat_orthogonal_grid(nx, nz, node_coords,constants=[maxLim[0]/2-notchWidth-gapWidth*2  ,indent,zind],sizes=[gapWidth*2-xind,1,maxLim[2]-indent*2-zind*2])

        pointGenerators.generate_flat_orthogonal_grid(nx, nz, node_coords,constants=[maxLim[0]/2+notchWidth+xind  ,indent,zind],sizes=[gapWidth*2-xind,1,maxLim[2]-indent*2-zind*2])

        pointGenerators.generate_flat_orthogonal_grid(nx, nz, node_coords,constants=[maxLim[0]/2-notchWidth-gapWidth*2  ,maxLim[1]-indent,zind],sizes=[gapWidth*2-xind,1,maxLim[2]-indent*2-zind*2])

        pointGenerators.generate_flat_orthogonal_grid(nx, nz, node_coords,constants=[maxLim[0]/2+notchWidth+xind  ,maxLim[1]-indent,zind],sizes=[gapWidth*2-xind,1,maxLim[2]-indent*2-zind*2])


        #left rigid plate top load
        ltopRigidPlateMechBC = np.array([-1,-1,0, -1,-1,-1,  -1,2,-1,-1,-1,-1])
        ltopRigidPlate = utilitiesMech.RigidPlate(-3, 3,
        np.array([
         0.5*maxLim[0]-notchWidth-gapWidth,
        0.5*maxLim[0]-notchWidth,
        maxLim[1] - 2*indentRP,
        maxLim[1] + 2*indentRP,
        -indentRP,
        maxLim[2] + indentRP  ]))

        rigidPlates.append(ltopRigidPlate)
        govNodes.append(np.array([ 0.5*maxLim[0]-notchWidth-gapWidth/2, maxLim[1]-indent, maxLim[2]/2 ]))
        govNodesMechBC.append(utilitiesMech.mechanicalBC(dim, -3, ltopRigidPlateMechBC))

        #right rigid plate top load
        rtopRigidPlateMechBC = np.array([-1,-1,0, -1,-1,-1,  -1,2,-1,-1,-1,-1])
        rtopRigidPlate = utilitiesMech.RigidPlate(-4, 3,
        np.array([
         0.5*maxLim[0]+notchWidth,
        0.5*maxLim[0]+notchWidth+gapWidth,
        maxLim[1] - 2*indentRP,
        maxLim[1] + 2*indentRP,
        -indentRP,
        maxLim[2] + indentRP  ]))

        rigidPlates.append(rtopRigidPlate)
        govNodes.append(np.array([ 0.5*maxLim[0]+notchWidth+gapWidth/2, maxLim[1]-indent, maxLim[2]/2 ]))
        govNodesMechBC.append(utilitiesMech.mechanicalBC(dim, -4, rtopRigidPlateMechBC))

        #rigid plate left pad
        indentRP = 1e-4
        leftRigidPlateMechBC = np.array([-1,-1,0, -1,-1,-1,  -1,3,-1,-1,-1,-1])
        leftRigidPlate = utilitiesMech.RigidPlate(-5, 3, np.array([ maxLim[0]/2-notchWidth-gapWidth-1e-3, maxLim[0]/2-notchWidth+1e-5, -indentRP, indentRP, -indentRP, maxLim[2]+indentRP  ]))
        rigidPlates.append(leftRigidPlate)
        leftSupportGN = np.array([ maxLim[0]/2-notchWidth-gapWidth/2, indent, maxLim[2]/2 ])
        govNodes.append(leftSupportGN)
        govNodesMechBC.append(utilitiesMech.mechanicalBC(dim, -5, leftRigidPlateMechBC))

        #rigid plate right pad
        indentRP = 1e-4
        leftRigidPlateMechBC = np.array([-1,-1,0, -1,-1,-1,  -1,3,-1,-1,-1,-1])
        leftRigidPlate = utilitiesMech.RigidPlate(-6, 3, np.array([ maxLim[0]/2+notchWidth-1e-5, maxLim[0]/2+notchWidth+gapWidth+1e-3, -indentRP, indentRP, -indentRP, maxLim[2]+indentRP  ]))
        rigidPlates.append(leftRigidPlate)
        leftSupportGN = np.array([ maxLim[0]/2+notchWidth+gapWidth/2, indent, maxLim[2]/2 ])
        govNodes.append(leftSupportGN)
        govNodesMechBC.append(utilitiesMech.mechanicalBC(dim, -6, leftRigidPlateMechBC))

    ##################### CONSTRAINTS AND RIGID PLATES
    #rigid plate left support
    indentRP = 1e-3
    leftRigidPlateMechBC = np.array([0,0,0, 0,0,0,  -1,-1,-1,-1,-1,-1])
    leftRigidPlate = utilitiesMech.RigidPlate(-1, 3, np.array([
    -indentRP,
    +indentRP,
    -indentRP,
    maxLim[1]+indentRP,
    -indentRP,
    maxLim[2]+indentRP  ]))
    rigidPlates.append(leftRigidPlate)
    govNodes.append(np.array([ indent, maxLim[1]/2 , maxLim[2]/2 ]))
    govNodesMechBC.append(utilitiesMech.mechanicalBC(dim, -1, leftRigidPlateMechBC))
    #rigid plate left support
    rightRigidPlateMechBC = np.array([1,0,0, 0,0,0,  -1,-1,-1,-1,-1,-1])
    rightRigidPlate = utilitiesMech.RigidPlate(-1, 3, np.array([
    -indentRP+maxLim[0],
    +indentRP+maxLim[0],
    -indentRP,
    maxLim[1]+indentRP,
    -indentRP,
    maxLim[2]+indentRP  ]))
    rigidPlates.append(rightRigidPlate)
    govNodes.append(np.array([ maxLim[0]-indent, maxLim[1]/2 , maxLim[2]/2 ]))
    govNodesMechBC.append(utilitiesMech.mechanicalBC(dim, -2, rightRigidPlateMechBC))
        ####################






    nodeA =  np.array([indent , indent, indent])
    nodeB =  np.array([indent, maxLim[1] - indent, maxLim[2] - indent])
    pointGenerators.generateNodesOrtoSurface3dRand(nodeA, nodeB, minDist, dim, node_coords, trials)
    nodeA =  np.array([maxLim[0]-indent , indent, indent])
    nodeB =  np.array([maxLim[0]-indent, maxLim[1] - indent, maxLim[2] - indent])
    pointGenerators.generateNodesOrtoSurface3dRand(nodeA, nodeB, minDist, dim, node_coords, trials)


    ##########################################generating of points, left support area
    maxLimF = np.array([
    L_sup_B+minDist*10,
    maxLim[1]- indent,
    maxLim[2]- indent,
    indent,
    indent,
    indent])
    #pointGenerators.generateNodesRect(maxLimF, minDist, dim, trials, node_coords, useLowBound=True)
    pointGenerators.generateNodesRect(maxLimF, minDist, dim, trials, node_coords, useLowBound=True, bottomMinDist = minDist, topMinDist = minDist*roughCoef, gradienDirection=0)
    ##########################################generating of points, right support area
    maxLimF = np.array([
    maxLim[0]- indent,
    maxLim[1]- indent,
    maxLim[2]- indent,
    R_sup_A-minDist*10,
    indent,
    indent])
    #pointGenerators.generateNodesRect(maxLimF, minDist, dim, trials, node_coords, useLowBound=True)
    pointGenerators.generateNodesRect(maxLimF, minDist, dim, trials, node_coords, useLowBound=True, bottomMinDist = minDist*roughCoef, topMinDist = minDist, gradienDirection=0)



    ##########################################generating of points, fracture zone

    maxLimF = np.array([
    maxLim[0] - indent - 0.5*maxLim[0]*(1-fracZoneWidth),
    maxLim[1] - indent,
    maxLim[2] - indent,
    indent + 0.5*maxLim[0]*(1-fracZoneWidth),
    maxLim[1]*notch/2*0+indent,
    indent])
    pointGenerators.generateNodesRect(maxLimF, minDist, dim, trials, node_coords, useLowBound=True)

    maxLimF = np.array([
    maxLim[0] - indent - 0.5*maxLim[0]*(1-fracZoneWidth)+4*minDist*roughCoef,
    maxLim[1] - indent,
    maxLim[2] - indent,
    maxLim[0] - indent - 0.5*maxLim[0]*(1-fracZoneWidth),
    maxLim[1]*notch/2*0+indent,
    indent])
    pointGenerators.generateNodesRect(maxLimF, minDist, dim, trials, node_coords, useLowBound=True, bottomMinDist = minDist, topMinDist = minDist*roughCoef, gradienDirection=0)



    maxLimF = np.array([
    indent + 0.5*maxLim[0]*(1-fracZoneWidth),
    maxLim[1] - indent,
    maxLim[2] - indent,
    indent + 0.5*maxLim[0]*(1-fracZoneWidth)-4*minDist*roughCoef,
    maxLim[1]*notch/2*0+indent,
    indent])
    pointGenerators.generateNodesRect(maxLimF, minDist, dim, trials, node_coords, useLowBound=True, bottomMinDist = minDist*roughCoef, topMinDist = minDist, gradienDirection=0)


    ##########################################generating of points, homogeneous volume
    pointGenerators.generateNodesRect(maxLim, minDist*roughCoef, dim, trials, node_coords)
    """
    node_coords=np.asarray(node_coords)
    fig = plt.figure()
    ax = Axes3D(fig)
    ax.scatter(node_coords[:,0], node_coords[:,1], node_coords[:,2])
    plt.show()
    """

    return node_coords, mechBC_merged, mechInitC_merged, notches, govNodes, govNodesMechBC, rigidPlates




######## FUNCTION FOR CREATING OF A 3D SUPPORTED RECTANGE
def assemble3dCantileverBending(maxLim, minDist, trials):
    indent = 1e-5
    dim=3
    #lists for the model
    node_coords = []
    mechBC_merged = []
    mechInitC_merged = []


    ###############generating of points loaded line top right ###############
    mechBC = np.array([-1,-1, 1, -1,-1,-1,    -1,-1,-1, -1,-1,-1])
    nodeA = np.array([maxLim[0] - indent , indent, maxLim[2] -indent])
    nodeB = np.array([maxLim[0] - indent , maxLim[1] - indent, maxLim[2] -indent])

    oldLen = len(node_coords)
    pointGenerators.generateNodesLine3dRand(nodeA, nodeB, minDist, dim, node_coords,  trials, True)
    nrOfPoints =  (len(node_coords)) - oldLen
    #print (nrOfPoints)

    for n in range ( nrOfPoints ):
        mBC = utilitiesMech.mechanicalBC(dim, oldLen + n, mechBC)
        mechBC_merged.append(mBC)

    ###############generating of points supported surface  left face ###############
    mechBC = np.array([0,0,0,0,0,0,-1,-1,-1,-1,-1,-1])

    nodeA = np.array([ indent , indent, indent])
    nodeB = np.array([ indent , maxLim[1] - indent, maxLim[2] -indent])

    oldLen = len(node_coords)
    pointGenerators.generateNodesOrtoSurface3dRand(nodeA, nodeB, minDist, dim, node_coords, trials)
    nrOfPoints =  (len(node_coords)) - oldLen
    for n in range ( nrOfPoints ):
        mBC = utilitiesMech.mechanicalBC(dim, oldLen + n, mechBC)
        mechBC_merged.append(mBC)


    ###############generating of points rectangular volume ###############
    mechBC = np.array([-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1])

    #kvadr
    oldLen = len(node_coords)
    pointGenerators.generateNodesRect(maxLim, minDist, dim, trials, node_coords)
    newLen = len(node_coords)-1
    nrOfPoints =  (len(node_coords)) - oldLen
   # for n in range ( nrOfPoints ):
       # mBC = utilitiesGeom.mechanicalBC(dim, oldLen + n, mechBC)
       # mechBC_merged.append(mBC)
       # print('adding')
    ####################################################################################################

    return node_coords, mechBC_merged, mechInitC_merged







######## FUNCTION FOR CREATING OF A 3D SUPPORTED RECTANGE
def assemble3dCantileverUniPressFree(maxLim, minDist, trials):
    indent = 1e-5
    dim=3
    #lists for the model
    node_coords = []
    mechBC_merged = []
    mechInitC_merged = []

    pointGenerators.generateNodesRect(maxLim, minDist*100, dim, trials, node_coords)

    ###############generating of points supported surface left face ###############
    mechBC = np.array([0,0,0,0,0,0,-1,-1,-1,-1,-1,-1])

    nodeA = np.array([ indent , indent, indent])
    nodeB = np.array([ indent , maxLim[1] - indent, maxLim[2] -indent])

    oldLen = len(node_coords)
    pointGenerators.generateNodesOrtoSurface3dRand(nodeA, nodeB, minDist, dim, node_coords, trials)
    nrOfPoints =  (len(node_coords)) - oldLen
    for n in range ( nrOfPoints ):
        mBC = utilitiesMech.mechanicalBC(dim, oldLen + n, mechBC)
        mechBC_merged.append(mBC)

    ###############generating of points loaded surface right face ###############
    mechBC = np.array([1, 0, 0, 0, 0, 0,   -1,-1,-1,-1,-1,-1])

    nodeA = np.array([ maxLim[0] - indent , indent, indent])
    nodeB = np.array([ maxLim[0] - indent , maxLim[1] - indent, maxLim[2] -indent])

    oldLen = len(node_coords)
    pointGenerators.generateNodesOrtoSurface3dRand(nodeA, nodeB, minDist, dim, node_coords, trials)
    nrOfPoints =  (len(node_coords)) - oldLen
    for n in range ( nrOfPoints ):
        mBC = utilitiesMech.mechanicalBC(dim, oldLen + n, mechBC)
        mechBC_merged.append(mBC)


    ###############generating of points rectangular volume ###############
    mechBC = np.array([-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1])

    oldLen = len(node_coords)
    pointGenerators.generateNodesRect(maxLim, minDist, dim, trials, node_coords)
    newLen = len(node_coords)-1
    nrOfPoints =  (len(node_coords)) - oldLen
   # for n in range ( nrOfPoints ):
       # mBC = utilitiesGeom.mechanicalBC(dim, oldLen + n, mechBC)
       # mechBC_merged.append(mBC)
       # print('adding')
    ####################################################################################################

    return node_coords, mechBC_merged, mechInitC_merged





######## FUNCTION FOR CREATING OF A 3D SUPPORTED RECTANGE
def assemble3dCantileverUniPressConfined(maxLim, minDist, trials):
    indent = 1e-5
    dim=3
    #lists for the model
    node_coords = []
    mechBC_merged = []
    mechInitC_merged = []

    pointGenerators.generateNodesRect(maxLim, minDist*100, dim, trials, node_coords)

    ###############generating of points supported surface left face ###############
    mechBC = np.array([0,0,0,0,0,0,-1,-1,-1,-1,-1,-1])

    nodeA = np.array([ indent , indent, indent])
    nodeB = np.array([ indent , maxLim[1] - indent, maxLim[2] -indent])

    oldLen = len(node_coords)
    pointGenerators.generateNodesOrtoSurface3dRand(nodeA, nodeB, minDist, dim, node_coords, trials)
    nrOfPoints =  (len(node_coords)) - oldLen
    for n in range ( nrOfPoints ):
        mBC = utilitiesMech.mechanicalBC(dim, oldLen + n, mechBC)
        mechBC_merged.append(mBC)




    ###############generating of points loaded surface right face ###############
    mechBC = np.array([1, 0, 0, 0, 0, 0,   -1,-1,-1,-1,-1,-1])

    nodeA = np.array([ maxLim[0] - indent , indent, indent])
    nodeB = np.array([ maxLim[0] - indent , maxLim[1] - indent, maxLim[2] -indent])

    oldLen = len(node_coords)
    pointGenerators.generateNodesOrtoSurface3dRand(nodeA, nodeB, minDist, dim, node_coords, trials)
    nrOfPoints =  (len(node_coords)) - oldLen
    for n in range ( nrOfPoints ):
        mBC = utilitiesMech.mechanicalBC(dim, oldLen + n, mechBC)
        mechBC_merged.append(mBC)



    ###############generating of points supported surface top face ###############
    mechBC = np.array([-1,0,0,  -1,-1,-1,    -1,-1,-1,-1,-1,-1])

    nodeA = np.array([ indent , indent, maxLim[2] -indent])
    nodeB = np.array([ maxLim[0] - indent , maxLim[1] - indent, maxLim[2] -indent])

    oldLen = len(node_coords)
    pointGenerators.generateNodesOrtoSurface3dRand(nodeA, nodeB, minDist, dim, node_coords, trials)
    nrOfPoints =  (len(node_coords)) - oldLen
    for n in range ( nrOfPoints ):
        mBC = utilitiesMech.mechanicalBC(dim, oldLen + n, mechBC)
        mechBC_merged.append(mBC)


    ###############generating of points supported surface bottom face ###############
    mechBC = np.array([-1,0,0,  -1,-1,-1,    -1,-1,-1,-1,-1,-1])

    nodeA = np.array([ indent , indent, indent])
    nodeB = np.array([ maxLim[0] - indent , maxLim[1] - indent, indent])

    oldLen = len(node_coords)
    pointGenerators.generateNodesOrtoSurface3dRand(nodeA, nodeB, minDist, dim, node_coords, trials)
    nrOfPoints =  (len(node_coords)) - oldLen
    for n in range ( nrOfPoints ):
        mBC = utilitiesMech.mechanicalBC(dim, oldLen + n, mechBC)
        mechBC_merged.append(mBC)

    ###############generating of points supported surface front face ###############
    mechBC = np.array([-1,0,0,  -1,-1,-1,    -1,-1,-1,-1,-1,-1])

    nodeA = np.array([ indent , indent, indent])
    nodeB = np.array([ maxLim[0] - indent , indent, maxLim[2]-indent])

    oldLen = len(node_coords)
    pointGenerators.generateNodesOrtoSurface3dRand(nodeA, nodeB, minDist, dim, node_coords, trials)
    nrOfPoints =  (len(node_coords)) - oldLen
    for n in range ( nrOfPoints ):
        mBC = utilitiesMech.mechanicalBC(dim, oldLen + n, mechBC)
        mechBC_merged.append(mBC)

    ###############generating of points supported surface rear face ###############
    mechBC = np.array([-1,0,0,  -1,-1,-1,    -1,-1,-1,-1,-1,-1])

    nodeA = np.array([ indent , maxLim[1]-indent, indent])
    nodeB = np.array([ maxLim[0] - indent , maxLim[1]-indent, maxLim[2]-indent])

    oldLen = len(node_coords)
    pointGenerators.generateNodesOrtoSurface3dRand(nodeA, nodeB, minDist, dim, node_coords, trials)
    nrOfPoints =  (len(node_coords)) - oldLen
    for n in range ( nrOfPoints ):
        mBC = utilitiesMech.mechanicalBC(dim, oldLen + n, mechBC)
        mechBC_merged.append(mBC)





    ###############generating of points rectangular volume ###############
    mechBC = np.array([-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1])

    oldLen = len(node_coords)
    pointGenerators.generateNodesRect(maxLim, minDist, dim, trials, node_coords)
    newLen = len(node_coords)-1
    nrOfPoints =  (len(node_coords)) - oldLen
   # for n in range ( nrOfPoints ):
       # mBC = utilitiesGeom.mechanicalBC(dim, oldLen + n, mechBC)
       # mechBC_merged.append(mBC)
       # print('adding')
    ####################################################################################################

    return node_coords, mechBC_merged, mechInitC_merged




















def assemble3dcylinderUniPressFree(center, radius, height, minDist, trials, directionDim):
    indent = 1e-5
    dim=3
    #lists for the model
    node_coords = []
    mechBC_merged = []
    mechInitC_merged = []

    mechBC = np.array([0,0,0,0,0,0,    -1,-1,-1,-1,-1,-1])
    node_coords.append( center+indent)
    mBC = utilitiesMech.mechanicalBC(dim, 0, mechBC)
    mechBC_merged.append(mBC)

    ###############generating of points supported surface left face ###############
    mechBC = np.array([0,-1,-1,-1,-1,-1,    -1,-1,-1,-1,-1,-1])

    oldLen = len(node_coords)
    pointGenerators.generateNodesOrtoCircleBorder3dRand(center, radius, directionDim, minDist, node_coords, trials)
    pointGenerators.generateNodesOrtoCircle3dRand(center, radius, directionDim, minDist, node_coords, trials)
    nrOfPoints =  (len(node_coords)) - oldLen
    for n in range ( nrOfPoints ):
        mBC = utilitiesMech.mechanicalBC(dim, oldLen + n, mechBC)
        mechBC_merged.append(mBC)
    #print('%d nodes generated so far' %len(node_coords))
    ###############generating of points loaded surface right face ###############
    mechBC = np.array([1,-1,-1,-1,-1,-1,    -1,-1,-1,-1,-1,-1])

    nodeA = center.copy()
    nodeA[directionDim] += height

    oldLen = len(node_coords)
    pointGenerators.generateNodesOrtoCircleBorder3dRand(nodeA, radius,  directionDim, minDist, node_coords, trials)
    pointGenerators.generateNodesOrtoCircle3dRand(nodeA, radius,  directionDim, minDist, node_coords, trials)
    nrOfPoints =  (len(node_coords)) - oldLen
    for n in range ( nrOfPoints ):
        mBC = utilitiesMech.mechanicalBC(dim, oldLen + n, mechBC)
        mechBC_merged.append(mBC)


    ###############generating of points rectangular volume ###############
    pointGenerators.generateNodesOrtoCilinderSurf3dRand(center, radius-1e-5, height, directionDim, minDist,  node_coords, trials)

    #######################################################################

    ###############generating of points rectangular volume ###############
    pointGenerators.generateNodesOrtoCilinder3dRand(center, radius, height, directionDim, minDist,  node_coords, trials)
    #######################################################################


    return node_coords, mechBC_merged, mechInitC_merged



def assemble3dcylinderTorsionFree(center, radius, height, minDist, trials, directionDim, functions, powerTes = False):
    indent = 1e-5
    dim=3
    #lists for the model
    node_coords = []
    mechBC_merged = []
    mechInitC_merged = []


    if powerTes == False:
        radii = None
    else:
        radii = []


    mechBC = np.array([0,0,0,0,0,0,    -1,-1,-1,-1,-1,-1])
    node_coords.append( center+indent)
    if powerTes == True:
        radii.append(0)

    mBC = utilitiesMech.mechanicalBC(dim, 0, mechBC)
    mechBC_merged.append(mBC)

    ###############generating of points supported surface left face ###############
    mechBC = np.array([0,0,0, 0 , 0, 0,    -1,-1,-1,-1,-1,-1])

    oldLen = len(node_coords)
    pointGenerators.generateNodesOrtoCircleBorder3dRand(center, radius, directionDim, minDist, node_coords, trials)
    pointGenerators.generateNodesOrtoCircle3dRand(center, radius, directionDim, minDist, node_coords, trials)
    nrOfPoints =  (len(node_coords)) - oldLen
    for n in range ( nrOfPoints ):
        mBC = utilitiesMech.mechanicalBC(dim, oldLen + n, mechBC)
        mechBC_merged.append(mBC)
    #print('%d nodes generated so far' %len(node_coords))
    ###############generating of points loaded surface right face ###############
    nodeA = center.copy()
    nodeA[directionDim] += height

    oldLen = len(node_coords)
    pointGenerators.generateNodesOrtoCircleBorder3dRand(nodeA, radius,  directionDim, minDist, node_coords, trials)
    pointGenerators.generateNodesOrtoCircle3dRand(nodeA, radius,  directionDim, minDist, node_coords, trials)
    nrOfPoints =  (len(node_coords)) - oldLen
    for n in range ( nrOfPoints ):
        lnfc = len(functions)
        mechBC = np.array([lnfc, lnfc+1, lnfc+2, -1 , -1, -1,    -1,-1,-1,-1,-1,-1])
        point = node_coords[oldLen + n]
        rotAngles = np.array([0.01, 0, 0])
        value = 1
        period = 1

        funcRot0 = utilitiesNumeric.constSawToothRotationFunction(rotAngles, point, value, 0, period=period, sym = 1)
        funcRot1 = utilitiesNumeric.constSawToothRotationFunction(rotAngles, point, value, 1, period=period, sym = 1)
        funcRot2 = utilitiesNumeric.constSawToothRotationFunction(rotAngles, point, value, 2, period=period, sym = 1)
        #//funcRot3 = utilitiesNumeric.constSawToothRotationFunction(rotAngles, point, value, 3, period=period)
        #//funcRot4 = utilitiesNumeric.constSawToothRotationFunction(rotAngles, point, value, 4, period=period)
        #//funcRot5 = utilitiesNumeric.constSawToothRotationFunction(rotAngles, point, value, 5, period=period)
        functions.append(funcRot0)
        functions.append(funcRot1)
        functions.append(funcRot2)
        #//functions.append(funcRot3)
        #//functions.append(funcRot4)
        #//functions.append(funcRot5)

        mBC = utilitiesMech.mechanicalBC(dim, oldLen + n, mechBC)
        mechBC_merged.append(mBC)


    ###############generating of points rectangular volume ###############
    pointGenerators.generateNodesOrtoCilinderSurf3dRand(center, radius-1e-5, height, directionDim, minDist,  node_coords, trials)

    #######################################################################

    ###############generating of points rectangular volume ###############
    pointGenerators.generateNodesOrtoCilinder3dRand(center, radius, height, directionDim, minDist,  node_coords, trials)
    #######################################################################


    return node_coords, mechBC_merged, mechInitC_merged



def assemble3dcylinderTorsionPressFree(center, radius, height, minDist, trials, directionDim, functions, powerTes=False):
    indent = 1e-5
    dim=3
    #lists for the model
    node_coords = []
    mechBC_merged = []
    mechInitC_merged = []
    govNodes = []
    govNodesMechBC = []
    rigidPlates = []


    if powerTes == False:
        radii = []
    else:
        radii = []

    ### fixed nodes
    mechBC = np.array([0,0,0,-1,-1,-1,    -1,-1,-1,-1,-1,-1])
    node_coords.append( np.array([indent*3,indent,indent]))
    mBC = utilitiesMech.mechanicalBC(dim, 0, mechBC)
    mechBC_merged.append(mBC)
    mechBC = np.array([-1,0,0,-1,-1,-1,    -1,-1,-1,-1,-1,-1])
    node_coords.append( np.array([height-indent*3,indent,indent]))
    mBC = utilitiesMech.mechanicalBC(dim, 1, mechBC)
    mechBC_merged.append(mBC)

    for i in range (len(node_coords)):
        radii.append(0)


    ### nodes for gauges
    #node_coords.append( center+indent)
    #node_coords.append( np.array([height-2*indent, 0, 0])  )

    ##################### CONSTRAINTS AND RIGID PLATES
    #rigid plate left support
    indentRP = indent
    leftRigidPlateMechBC = np.array([0,-1,-1, 0,0,0,  -1,-1,-1,-1,-1,-1])
    leftRigidPlate = utilitiesMech.RigidPlate(-1, 3,
    np.array([ -indentRP,
     2*indentRP,
     -2*radius,
      2*radius,
      -2*radius,
       2*radius ]))
    rigidPlates.append(leftRigidPlate)
    govNodes.append(np.array([ 0, 0, 0]))
    govNodesMechBC.append(utilitiesMech.mechanicalBC(dim, -1, leftRigidPlateMechBC))
    #rigid plate left support
    rightRigidPlateMechBC = np.array([1,-1,-1, 2,0,0,  -1,-1,-1, -1,-1,-1])
    rightRigidPlate = utilitiesMech.RigidPlate(-2, 3,np.array([
    height-2*indentRP,
    height+2*indentRP,
     -2*radius,
      2*radius,
      -2*radius,
       2*radius ]))
    rigidPlates.append(rightRigidPlate)
    govNodes.append(np.array([ height, 0, 0]))
    govNodesMechBC.append(utilitiesMech.mechanicalBC(dim, -2, rightRigidPlateMechBC))

    print ('Nodes so far: %d' %len(node_coords))

    ###############generating of points supported surface left face ###############
    if powerTes==False or 1:
        pointGenerators.generateNodesOrtoCircleBorder3dRand(center, radius-1e-5, directionDim, minDist*0.4, node_coords, trials)
        pointGenerators.generateNodesOrtoCircle3dRand(center, radius-1e-5, directionDim, minDist*0.4, node_coords, trials)
    else:
        node_coords, radii=pointGenerators.generateParticlesOrtoCircle3dRand(center, radius-1e-5, directionDim, minDist*0.4, minDist, 0.8,  trials, node_coords, radii)



    ###############generating of points loaded surface right face ###############
    nodeA = center.copy()
    nodeA[directionDim] += height-indent
    if powerTes==False or 1:
        pointGenerators.generateNodesOrtoCircleBorder3dRand(nodeA, radius-1e-5,  directionDim, minDist*0.4, node_coords, trials)
        pointGenerators.generateNodesOrtoCircle3dRand(nodeA, radius-1e-5,  directionDim, minDist*0.4, node_coords, trials)
    else:
        node_coords, radii=pointGenerators.generateParticlesOrtoCircle3dRand(nodeA, radius-1e-5, directionDim, minDist*0.4, minDist, 0.8,  trials, node_coords, radii)


    ###############generating of points cylinder surf###############
    #pointGenerators.generateNodesOrtoCilinderSurf3dRand(center, radius-1e-5, height, directionDim, minDist,  node_coords, trials)


    if powerTes==False:
    ###############generating of points cylinder volume ###############
        pointGenerators.generateNodesOrtoCilinder3dRand(center, radius-1e-5, height, directionDim, minDist,  node_coords, trials)
        radii = []
    else:
        radii = np.zeros(len(node_coords))
        node_coords = np.array(node_coords)
        node_coords, radii=pointGenerators.generateParticlesOrtoCilinder3dRand(center, radius-1e-5, height, directionDim, minDist*0.4, minDist, 0.8,  trials, node_coords, radii)



    """
    node_coords = np.asarray(node_coords)
    if SHOW_PLOT:
        fig = plt.figure()
        ax = Axes3D(fig)
        ax.scatter(node_coords[:,0], node_coords[:,1], node_coords[:,2])
        plt.show()
    #"""




    return node_coords, mechBC_merged,  govNodes, govNodesMechBC, rigidPlates, radii


def assemble3dRWTHShearCylinder(center, radiusInp, heightInp, minDist, trials, directionDim, functions, notchRadLeft, notchRadRight, notchWidthLeft, notchWidthRight, notchDepth, quarter=False):
    indent = 1e-5
    dim=3

    radius = radiusInp - 1e-5
    height = heightInp - 1e-5
    #lists for the model
    node_coords = []
    mechBC_merged = []
    mechInitC_merged = []
    govNodes = []
    govNodesMechBC = []
    rigidPlates = []
    notches = []

    notchOuterRadLeft = notchRadLeft + notchWidthLeft/2
    notchInnerRadLeft = notchRadLeft - notchWidthLeft/2

    notchOuterRadRight = notchRadRight + notchWidthRight/2
    notchInnerRadRight = notchRadRight - notchWidthRight/2

    ### fixed nodes
    mechBC = np.array([-1,0,0,-1,-1,-1,    -1,-1,-1,-1,-1,-1])
    node_coords.append( np.array([indent*3,indent*3,indent*3]))
    mBC = utilitiesMech.mechanicalBC(dim, 0, mechBC)
    #mechBC_merged.append(mBC)
    mechBC = np.array([-1,0,0,-1,-1,-1,    -1,-1,-1,-1,-1,-1])
    node_coords.append( np.array([height-indent*3,indent,indent]))
    mBC = utilitiesMech.mechanicalBC(dim, 1, mechBC)
    #mechBC_merged.append(mBC)


    ### nodes for gauges
    #node_coords.append( center+indent)
    #node_coords.append( np.array([height-2*indent, 0, 0])  )

    ##################### CONSTRAINTS AND RIGID PLATES
    #rigid plate left support
    indentRP = indent
    leftRigidPlateMechBC = np.array([0,-1,-1, 0,0,0,  -1,-1,-1,   -1,-1,-1])
    leftRigidPlate = utilitiesMech.RigidPlate(-1, 3, np.array([ center[0],0, 0,notchInnerRadLeft, 1e-5 ]), radial=0, innerRad=None)
    rigidPlates.append(leftRigidPlate)
    govNodes.append(np.array([ 0, 0, 0]))
    govNodesMechBC.append(utilitiesMech.mechanicalBC(dim, -1, leftRigidPlateMechBC))


    #rigid plate left support
    rightRigidPlateMechBC = np.array([1,-1,-1, 0,0,0,  -1,-1,-1, -1,-1,-1])
    rightRigidPlate = utilitiesMech.RigidPlate(-2, 3, np.array([ height,0, 0,radius, 1e-6 ]), radial=0, innerRad = notchOuterRadRight-1e-5)
    #rightRigidPlate = utilitiesMech.RigidPlate(-1, 3, np.array([ 0,0, 0,notchInnerRadLeft, 1e-5 ]), radial=0, innerRad=None)
    rigidPlates.append(rightRigidPlate)
    govNodes.append(np.array([ height, 0, 0]))
    govNodesMechBC.append(utilitiesMech.mechanicalBC(dim, -2, rightRigidPlateMechBC))


    boxIndent = np.array([0,0,0])
    ###############generating of points left face ###############
    notchSideOuter = []
    notchSideInner = []


    roughMinDist = 0.003

    pointGenerators.generateNodesOrtoCircleBorder3dRand(center+boxIndent, radius, directionDim, roughMinDist, node_coords, trials)
    print ('Nodes so far: %d' %len(node_coords))

    oldLen = len(node_coords)
    pointGenerators.generateNodesOrtoCircleBorder3dRand(center+boxIndent, notchOuterRadLeft, directionDim, minDist, node_coords, trials)
    print ('Nodes so far: %d' %len(node_coords))
    pointGenerators.generateNodesOrtoCilinderSurf3dRand(center+boxIndent, notchOuterRadLeft, notchDepth, directionDim, roughMinDist,  node_coords, trials)
    for i in range (oldLen, len(node_coords), 1):
        notchSideOuter.append(i)

    print ('Nodes so far: %d' %len(node_coords))

    oldLen = len(node_coords)
    pointGenerators.generateNodesOrtoCircleBorder3dRand(center, notchInnerRadLeft, directionDim, minDist, node_coords, trials)
    print ('Nodes so far: %d' %len(node_coords))
    pointGenerators.generateNodesOrtoCilinderSurf3dRand(center, notchInnerRadLeft, notchDepth, directionDim, roughMinDist,  node_coords, trials)
    for i in range (oldLen, len(node_coords), 1):
        notchSideInner.append(i)

    oldLen = len(node_coords)
    pointGenerators.generateNodesOrtoAnnulus3dRand(center, radius, (radius-notchOuterRadLeft), directionDim, roughMinDist, node_coords, trials)
    pointGenerators.generateNodesOrtoTube3dRand(center+boxIndent, radius, notchDepth, radius-notchOuterRadLeft, directionDim, roughMinDist,  node_coords, trials*2)
    for i in range (oldLen, len(node_coords), 1):
        notchSideOuter.append(i)

    print ('Nodes so far: %d' %len(node_coords))

    oldLen = len(node_coords)
    pointGenerators.generateNodesOrtoCircle3dRand(center, notchInnerRadLeft,  directionDim, roughMinDist, node_coords, trials)
    print ('Nodes so far: %d' %len(node_coords))
    pointGenerators.generateNodesOrtoCilinder3dRand(center, notchInnerRadLeft, notchDepth, directionDim, roughMinDist,  node_coords, trials*2)
    for i in range (oldLen, len(node_coords), 1):
        notchSideInner.append(i)

    print ('Nodes so far: %d' %len(node_coords))
    center = np.array([height, 0 ,0])
    ###############generating of points right face ###############
    pointGenerators.generateNodesOrtoCircleBorder3dRand(center, radius, directionDim, roughMinDist, node_coords, trials)

    oldLen = len(node_coords)
    pointGenerators.generateNodesOrtoCircleBorder3dRand(center, notchOuterRadRight, directionDim, minDist, node_coords, trials)
    print ('Nodes so far: %d' %len(node_coords))
    pointGenerators.generateNodesOrtoCilinderSurf3dRand(center, notchOuterRadRight, -notchDepth, directionDim, minDist,  node_coords, trials)
    print ('Nodes so far: %d' %len(node_coords))
    pointGenerators.generateNodesOrtoAnnulus3dRand(center-boxIndent, radius, (radius-notchOuterRadRight), directionDim, roughMinDist, node_coords, trials)
    print ('Nodes so far: %d' %len(node_coords))
    pointGenerators.generateNodesOrtoTube3dRand(center, radius, -notchDepth, radius-notchOuterRadRight, directionDim, roughMinDist,  node_coords, trials*2)
    for i in range (oldLen, len(node_coords), 1):
        notchSideOuter.append(i)

    oldLen = len(node_coords)
    pointGenerators.generateNodesOrtoCircle3dRand(center-boxIndent, notchInnerRadRight,  directionDim, roughMinDist, node_coords, trials)
    print ('Nodes so far: %d' %len(node_coords))
    pointGenerators.generateNodesOrtoCilinderSurf3dRand(center-boxIndent, notchInnerRadRight, -notchDepth, directionDim, minDist,  node_coords, trials)
    print ('Nodes so far: %d' %len(node_coords))
    pointGenerators.generateNodesOrtoCircleBorder3dRand(center-boxIndent, notchInnerRadRight, directionDim, minDist, node_coords, trials)
    print ('Nodes so far: %d' %len(node_coords))
    pointGenerators.generateNodesOrtoCilinder3dRand(center-boxIndent, notchInnerRadRight, -notchDepth, directionDim, roughMinDist,  node_coords, trials*2)
    print ('Nodes so far: %d' %len(node_coords))
    for i in range (oldLen, len(node_coords), 1):
        notchSideInner.append(i)


    notchInnerVolume = []
    #"""
    #inner rough cylinder

    center = np.array([notchDepth, 0 ,0])
    notchRad = (notchRadLeft + notchRadRight)/2
    fineTubeInnerDim =  notchRad * 0.6
    fineTubeOuterDim =  notchRad * 1.4

    pointGenerators.generateNodesOrtoCilinder3dRand(center, fineTubeInnerDim, height-2*notchDepth, directionDim, roughMinDist,  node_coords, trials*2)
    pointGenerators.generateNodesOrtoTube3dRand(center, radius, height-2*notchDepth, (radius-fineTubeOuterDim), directionDim, roughMinDist,  node_coords, trials*2)

    oldLen = len(node_coords)
    pointGenerators.generateNodesOrtoTube3dRand(center, fineTubeOuterDim, height-2*notchDepth, (fineTubeOuterDim-fineTubeInnerDim), directionDim, minDist,  node_coords, trials)
    for i in range (oldLen, len(node_coords), 1):
        notchInnerVolume.append(i)
    print ('Nodes so far: %d' %len(node_coords))



    notchSideInnerNew = []
    notchSideOuterNew = []
    notchInnerVolumeNew = []

    if not quarter:
        notchSideInnerNew = notchSideInner
        notchSideOuterNew = notchSideOuter
        notchInnerVolumeNew = notchInnerVolume

    if quarter == True:
        print ('Dumping points outside quarter...')
        quarter_nodes = []
        for n in range (len(node_coords)):
            no = node_coords[n]
            if (no[1]>0 and no[2]>0):
                if n in notchSideInner:
                    notchSideInnerNew.append(len(quarter_nodes))
                if n in notchSideOuter:
                    notchSideOuterNew.append(len(quarter_nodes))
                if n in notchInnerVolume:
                    notchInnerVolumeNew.append(len(quarter_nodes))
                #
                quarter_nodes.append(no)

        indent = 1e-3
        node_coords = []
        print('old len %d' %len(node_coords))
        for n in range (len(quarter_nodes)):
            node_coords.append(quarter_nodes[n])
        print('new  len %d' %len(node_coords))

        print('Generating border surfaces...')
        mechBCsides = np.array([-1,0,0, 0,0,0,    -1,-1,-1,-1,-1,-1])
        veryOldLen = len (node_coords)

        surfMinDist = minDist

        oldLen = len(node_coords)
        nodeA =  np.array([indent, indent ,indent])
        nodeB =  np.array([notchDepth, indent , notchInnerRadLeft])
        pointGenerators.generateNodesOrtoSurface3dRand(nodeA, nodeB, surfMinDist, dim, node_coords, trials)
        nodeA =  np.array([height-notchDepth, indent ,indent])
        nodeB =  np.array([height, indent , notchInnerRadRight])
        pointGenerators.generateNodesOrtoSurface3dRand(nodeA, nodeB, surfMinDist, dim, node_coords, trials)
        for i in range (oldLen, len(node_coords), 1):
            notchSideInnerNew.append(i)

        oldLen = len(node_coords)
        nodeA =  np.array([indent, indent , notchOuterRadLeft])
        nodeB =  np.array([notchDepth, indent , radius])
        pointGenerators.generateNodesOrtoSurface3dRand(nodeA, nodeB, surfMinDist, dim, node_coords, trials)
        nodeA =  np.array([height-notchDepth, indent , notchOuterRadRight])
        nodeB =  np.array([height, indent , radius])
        pointGenerators.generateNodesOrtoSurface3dRand(nodeA, nodeB, surfMinDist, dim, node_coords, trials)
        for i in range (oldLen, len(node_coords), 1):
            notchSideOuterNew.append(i)


        oldLen = len(node_coords)
        nodeA =  np.array([notchDepth, indent , indent])
        nodeB =  np.array([height-notchDepth, indent , radius])
        pointGenerators.generateNodesOrtoSurface3dRand(nodeA, nodeB, surfMinDist, dim, node_coords, trials)
        #for i in range (oldLen, len(node_coords), 1):
        #    notchInnerVolume.append(i)


        oldLen = len(node_coords)
        nodeA =  np.array([indent, indent ,indent])
        nodeB =  np.array([notchDepth, notchInnerRadLeft , indent])
        pointGenerators.generateNodesOrtoSurface3dRand(nodeA, nodeB, surfMinDist, dim, node_coords, trials)
        nodeA =  np.array([height-notchDepth, indent ,indent])
        nodeB =  np.array([height, notchInnerRadRight, indent ])
        pointGenerators.generateNodesOrtoSurface3dRand(nodeA, nodeB, surfMinDist, dim, node_coords, trials)
        for i in range (oldLen, len(node_coords), 1):
            notchSideInnerNew.append(i)

        oldLen = len(node_coords)
        nodeA =  np.array([indent, notchOuterRadLeft, indent ])
        nodeB =  np.array([notchDepth, radius, indent ])
        pointGenerators.generateNodesOrtoSurface3dRand(nodeA, nodeB, surfMinDist, dim, node_coords, trials)
        nodeA =  np.array([height-notchDepth, notchOuterRadRight, indent])
        nodeB =  np.array([height, radius, indent])
        pointGenerators.generateNodesOrtoSurface3dRand(nodeA, nodeB, surfMinDist, dim, node_coords, trials)
        for i in range (oldLen, len(node_coords), 1):
            notchSideOuterNew.append(i)

        nodeA =  np.array([notchDepth, indent , indent])
        nodeB =  np.array([height-notchDepth, radius, indent  ])
        pointGenerators.generateNodesOrtoSurface3dRand(nodeA, nodeB, surfMinDist, dim, node_coords, trials)
        #for i in range (oldLen, len(node_coords), 1):
        #    notchInnerVolume.append(i)

        for i in range (veryOldLen, len(node_coords), 1):
            mBC = utilitiesMech.mechanicalBC(dim, i, mechBC)
            mechBC_merged.append(mBC)





    notchL = []
    notchL.append(notchSideOuterNew)
    notchL.append(notchSideInnerNew)
    notches.append(notchL)

    notchA = []
    notchA.append(notchSideOuterNew)
    notchA.append(notchInnerVolume)
    notches.append(notchA)

    notchB = []
    notchB.append(notchSideInnerNew)
    notchB.append(notchInnerVolume)
    notches.append(notchB)


    minX = 1000
    maxX = -1000
    maxRad = 0
    for n in node_coords:
        if n[0] > maxX: maxX = n[0]
        if n[0] < minX: minX = n[0]
        if np.linalg.norm(n[1:3]) > maxRad: maxRad = np.linalg.norm(n[1:3])
        if math.isnan(n[0]): print('nan')
        if math.isnan(n[1]): print('nan')
        if math.isnan(n[2]): print('nan')


    print ('minX %f maxX %f maxRad %f' %(minX, maxX,maxRad))

    """
    node_coords = np.asarray(node_coords)
    if SHOW_PLOT:
        fig = plt.figure()
        ax = Axes3D(fig)
        ax.scatter(node_coords[:,0], node_coords[:,1], node_coords[:,2])
        plt.show()
    #
    """

    #"""

    return node_coords, mechBC_merged,  govNodes, govNodesMechBC, rigidPlates, notches




def assemble2dRWTHShearCylinder (maxLim, minDist, trials, innerRadTop, innerRadBottom, notchWidth, notchDepth):
    dim = 2
    #lists for the model
    node_coords = []
    mechBC_merged = []
    mechInitC_merged = []
    govNodes = []
    govNodesMechBC = []
    rigidPlates = []

    notch=1

    #an indent due to mirroring of the data for voronoi tess.
    notches=[]
    indent = 1e-8
    #notchWidth = 1.5e-3 /2

    #radiusy nejsou kotovany na osu, ale na vnitrni hranu notche !!!!
    #generating notch points
    if (notch > 0):
        notchSide0 = []
        notchSide1 = []
        notchVolume = []


        # HORNI STRANA PRAVY NOTCH
        nodeA = np.array([maxLim[0]/2-notchWidth-innerRadTop, maxLim[1]-indent])
        nodeB = np.array([maxLim[0]/2-notchWidth-innerRadTop, maxLim[1]-notchDepth+minDist])
        oldLen = len(node_coords)
        pointGenerators.generateNodesLine2dRand(nodeA, nodeB, minDist, dim, node_coords, trials, catchCorners=True, equidist=True)
        for i in range (oldLen, len(node_coords), 1):
            notchSide0.append(i)


        nodeA = np.array([maxLim[0]/2-innerRadTop, maxLim[1]-indent])
        nodeB = np.array([maxLim[0]/2-innerRadTop, maxLim[1]-notchDepth+minDist])
        oldLen = len(node_coords)
        pointGenerators.generateNodesLine2dRand(nodeA, nodeB, minDist, dim, node_coords, trials, catchCorners=True, equidist=True)
        for i in range (oldLen, len(node_coords), 1):
            notchSide1.append(i)

        node_coords.append(np.array([maxLim[0]/2-notchWidth-innerRadTop, maxLim[1]-notchDepth]))
        node_coords.append(np.array([maxLim[0]/2-innerRadTop, maxLim[1]-notchDepth]))

        # HORNI STRANA LEVY NOTCH
        nodeA = np.array([maxLim[0]/2+innerRadTop, maxLim[1]-indent])
        nodeB = np.array([maxLim[0]/2+innerRadTop, maxLim[1]-notchDepth+minDist])
        oldLen = len(node_coords)
        pointGenerators.generateNodesLine2dRand(nodeA, nodeB, minDist, dim, node_coords, trials, catchCorners=True, equidist=True)
        for i in range (oldLen, len(node_coords), 1):
            notchSide0.append(i)

        nodeA = np.array([maxLim[0]/2+notchWidth+innerRadTop, maxLim[1]-indent])
        nodeB = np.array([maxLim[0]/2+notchWidth+innerRadTop, maxLim[1]-notchDepth+minDist])
        oldLen = len(node_coords)
        pointGenerators.generateNodesLine2dRand(nodeA, nodeB, minDist, dim, node_coords, trials, catchCorners=True, equidist=True)
        for i in range (oldLen, len(node_coords), 1):
            notchSide1.append(i)

        node_coords.append(np.array([maxLim[0]/2+notchWidth+innerRadTop, maxLim[1]-notchDepth]))
        node_coords.append(np.array([maxLim[0]/2+innerRadTop, maxLim[1]-notchDepth]))



        # DOLNI STRANA PRAVY NOTCH
        nodeA = np.array([maxLim[0]/2-notchWidth-innerRadBottom, indent])
        nodeB = np.array([maxLim[0]/2-notchWidth-innerRadBottom, notchDepth-minDist])
        oldLen = len(node_coords)
        pointGenerators.generateNodesLine2dRand(nodeA, nodeB, minDist, dim, node_coords, trials, catchCorners=True, equidist=True)
        for i in range (oldLen, len(node_coords), 1):
            notchSide0.append(i)

        nodeA = np.array([maxLim[0]/2-innerRadBottom, indent])
        nodeB = np.array([maxLim[0]/2-innerRadBottom, notchDepth-minDist])
        oldLen = len(node_coords)
        pointGenerators.generateNodesLine2dRand(nodeA, nodeB, minDist, dim, node_coords, trials, catchCorners=True, equidist=True)
        for i in range (oldLen, len(node_coords), 1):
            notchSide1.append(i)

        node_coords.append(np.array([maxLim[0]/2-notchWidth-innerRadBottom, indent+notchDepth]))
        node_coords.append(np.array([maxLim[0]/2-innerRadBottom, indent+notchDepth]))

        # DOLNI STRANA LEVY NOTCH
        nodeA = np.array([maxLim[0]/2+notchWidth+innerRadBottom, indent])
        nodeB = np.array([maxLim[0]/2+notchWidth+innerRadBottom, notchDepth-minDist])
        oldLen = len(node_coords)
        pointGenerators.generateNodesLine2dRand(nodeA, nodeB, minDist, dim, node_coords, trials, catchCorners=True, equidist=True)
        for i in range (oldLen, len(node_coords), 1):
            notchSide0.append(i)

        nodeA = np.array([maxLim[0]/2+innerRadBottom, indent])
        nodeB = np.array([maxLim[0]/2+innerRadBottom, notchDepth-minDist])
        oldLen = len(node_coords)
        pointGenerators.generateNodesLine2dRand(nodeA, nodeB, minDist, dim, node_coords, trials, catchCorners=True, equidist=True)
        for i in range (oldLen, len(node_coords), 1):
            notchSide1.append(i)

        node_coords.append(np.array([maxLim[0]/2+notchWidth+innerRadBottom, indent+notchDepth]))
        node_coords.append(np.array([maxLim[0]/2+innerRadBottom, indent+notchDepth]))




        notchA = []
        notchA.append(notchSide0)
        notchA.append(notchSide1)
        notches.append(notchA)


    ##################### CONSTRAINTS AND RIGID PLATES
    #rigid plate left support
    indentRP = 1e-3
    leftRigidPlateMechBC = np.array([0,0,0, -1,-1,-1])
    leftRigidPlate = utilitiesMech.RigidPlate(-1, 2, np.array([
    -indentRP, maxLim[0]/2-innerRadBottom-notchWidth+1e-3, -indentRP, 2*indentRP ]))
    rigidPlates.append(leftRigidPlate)
    govNodes.append(np.array([ indent, indent ]))
    govNodesMechBC.append(utilitiesMech.mechanicalBC(dim, -1, leftRigidPlateMechBC))
    #rigid plate left support
    rightRigidPlateMechBC = np.array([0,0,0, -1,-1,-1])
    rightRigidPlate = utilitiesMech.RigidPlate(-2, 2, np.array([
    maxLim[0]/2+innerRadBottom+notchWidth-1e-3, 2*maxLim[0]/2+1e-4,  -indentRP, 2*indentRP ]))
    rigidPlates.append(rightRigidPlate)
    govNodes.append(np.array([ maxLim[0], indent ]))
    govNodesMechBC.append(utilitiesMech.mechanicalBC(dim, -2, rightRigidPlateMechBC))


    #rigid plate top load
    topRigidPlateMechBC = np.array([0,1,0, -1,-1,-1])
    topRigidPlate = utilitiesMech.RigidPlate(-3, 2,
    np.array([
    maxLim[0]/2-innerRadTop-1e-4,
    maxLim[0]/2+innerRadTop+1e-4,
    maxLim[1] - 2*indentRP,
    maxLim[1] + 2*indentRP  ]))
    rigidPlates.append(topRigidPlate)
    govNodes.append(np.array([ maxLim[0]/2, maxLim[1]-indent ]))
    govNodesMechBC.append(utilitiesMech.mechanicalBC(dim, -3, topRigidPlateMechBC))
    ####################


    ###############generating of nodes, left horizontal support ###############
    #defining points of the line
    nodeA = np.array([indent, indent])
    nodeB = np.array([indent + maxLim[0]/2-innerRadBottom-notchWidth, indent])
    pointGenerators.generateNodesLine2dRand(nodeA, nodeB, minDist, dim, node_coords, trials, False, False)
    ###############generating of nodes, right horizontal support ###############
    #defining points of the line
    nodeA = np.array([maxLim[0]/2-innerRadBottom, indent])
    nodeB = np.array([maxLim[0]/2+innerRadBottom, indent])
    #pointGenerators.generateNodesLine2dRand(nodeA, nodeB, minDist, dim, node_coords, trials, True, False)

    #defining points of the line
    nodeA = np.array([maxLim[0]/2+innerRadBottom+notchWidth, indent])
    nodeB = np.array([2*maxLim[0]/2-indent, indent])
    pointGenerators.generateNodesLine2dRand(nodeA, nodeB, minDist, dim, node_coords, trials, False, False)

    ############### loaded top face ###############
    nodeA = np.array([maxLim[0]/2-innerRadTop, maxLim[1]-indent])
    nodeB = np.array([maxLim[0]/2+innerRadTop, maxLim[1]-indent])
    pointGenerators.generateNodesLine2dRand(nodeA, nodeB, minDist, dim, node_coords,  trials, False, True)



    #top left rect
    maxLimF = np.array([
    indent,
    maxLim[1]-notchDepth,
    maxLim[0]/2-innerRadTop-notchWidth,
    maxLim[1]-indent])
    pointGenerators.generateNodesRect(maxLimF, minDist*2, dim, trials*2, node_coords, useLowBound=True)

    #top mid rect
    maxLimF = np.array([
    maxLim[0]/2-innerRadTop,
    maxLim[1]-notchDepth,
    maxLim[0]/2+innerRadTop,
    maxLim[1]-indent])
    pointGenerators.generateNodesRect(maxLimF, minDist, dim, trials, node_coords, useLowBound=True)

    #top right rect
    maxLimF = np.array([
    maxLim[0]/2+innerRadTop+notchWidth,
    maxLim[1]-notchDepth,
    maxLim[0]/2*2 - indent,
    maxLim[1]-indent])
    pointGenerators.generateNodesRect(maxLimF, minDist*2, dim, trials*2, node_coords, useLowBound=True)

    notchVolume = []
    oldLen = len(node_coords)
    #mid mid rect
    maxLimF = np.array([
    indent,
    notchDepth,
    2*maxLim[0]/2-indent,
    maxLim[1]-notchDepth])
    pointGenerators.generateNodesRect(maxLimF, minDist/1, dim, trials, node_coords, useLowBound=True)
    for i in range (oldLen, len(node_coords), 1):
        notchVolume.append(i)

    notchB = []
    notchB.append(notchSide0)
    notchB.append(notchVolume)
    notches.append(notchB)

    notchC = []
    notchC.append(notchSide1)
    notchC.append(notchVolume)
    notches.append(notchC)

    #bot left rect
    maxLimF = np.array([
    indent,
    notchDepth,
    maxLim[0]/2-innerRadBottom-notchWidth,
    indent])
    pointGenerators.generateNodesRect(maxLimF, minDist, dim, trials, node_coords, useLowBound=True)

    #bot mid rect
    maxLimF = np.array([
    maxLim[0]/2-innerRadBottom,
    notchDepth,
    maxLim[0]/2+innerRadBottom,
    indent])
    pointGenerators.generateNodesRect(maxLimF, minDist, dim, trials, node_coords, useLowBound=True)

    #bot right rect
    maxLimF = np.array([
    maxLim[0]/2+innerRadBottom+notchWidth,
    notchDepth,
    maxLim[0]/2*2 - indent,
    indent])
    pointGenerators.generateNodesRect(maxLimF, minDist, dim, trials, node_coords, useLowBound=True)




    ##"""

    return node_coords, mechBC_merged, mechInitC_merged, notches, govNodes, govNodesMechBC, rigidPlates






def assembleCoupledBrazilianDisc(center, radius, height, minDist, trials, directionDim, powerTes=False ):
    indent = 1e-6
    dim=3
    #lists for the model
    node_coords = []
    mechBC_merged = []
    mechInitC_merged = []
    govNodes = []
    govNodesMechBC = []
    rigidPlates = []

    if powerTes == False:
        radii = None
    else:
        radii = []

    lineSupported = True

    node_coords.append( np.array([height/2, 0, radius/2]))
    node_coords.append( np.array([height/2, 0, -radius/2]))

    node_coords.append( np.array([height/2, 0, 0.0508/2]))
    node_coords.append( np.array([height/2, 0, -0.0508/2]))

    # tady se natvrdo sampluji nody ve vetknuti
    if lineSupported:
        #node_coords.append( np.array([height/2, radius/4, radius/4]))

        oldLen = len(node_coords)
        width = radius*0.24
        step = minDist*0.35
        da = int(width/step)
        pos = np.linspace(-width/2., width/2., da)
        for p in pos:
            nodeA = np.array([indent, radius*0.99, p])
            nodeB = np.array([height-indent, radius*0.99, p])
            pointGenerators.generateNodesLine3dRand(nodeA, nodeB, step, dim, node_coords, trials, catchCorners=True, equidist = True)
        newNodes = len(node_coords)-oldLen


        rpIdcs = []
        for i in range (newNodes):
            rpIdcs.append(oldLen + i)

        topRigidPlate = utilitiesMech.RigidPlate(-1, 3, None, directIdcs = True)
        topRigidPlate.setDirectNodes(rpIdcs)
        topRigidPlateMechBC = np.array([0,1,0, 0,0,0,  -1,-1,-1,-1,-1,-1])
        rigidPlates.append(topRigidPlate)
        govNodes.append(np.array([ height/2, radius, 0]))
        govNodesMechBC.append(utilitiesMech.mechanicalBC(dim, -1, topRigidPlateMechBC))

        oldLen = len(node_coords)
        for p in pos:
            nodeA = np.array([indent, -radius*0.99, p])
            nodeB = np.array([height-indent, -radius*0.99, p])
            pointGenerators.generateNodesLine3dRand(nodeA, nodeB, step, dim, node_coords, trials, catchCorners=True, equidist = True)
        newNodes = len(node_coords)-oldLen

        print()
        rpIdcs = []
        for i in range (newNodes):
            rpIdcs.append(oldLen + i)


        bottomRigidPlate = utilitiesMech.RigidPlate(-2, 3, None, directIdcs = True)
        bottomRigidPlate.setDirectNodes(rpIdcs)
        bottomRigidPlateMechBC = np.array([0,0,0, 0,0,0,  -1,-1,-1,-1,-1,-1])
        rigidPlates.append(bottomRigidPlate)
        govNodes.append(np.array([ height/2, -radius, 0]))
        govNodesMechBC.append(utilitiesMech.mechanicalBC(dim, -2, bottomRigidPlateMechBC))

        if powerTes == True:
            for i in range (len(node_coords)):
                radii.append(0)


    if not lineSupported:
        mechBC = np.array([0,-1,-1,-1,-1,-1,    -1,-1,-1,-1,-1,-1])

        node_coords.append( np.array([height/2, +radius-2*indent, 0]))
        mBC = utilitiesMech.mechanicalBC(dim, 0, mechBC)
        mechBC_merged.append(mBC)

        mechBC = np.array([0,0,0,-1,-1,-1,    -1,-1,-1,-1,-1,-1])
        node_coords.append( np.array([height/2, -radius+2*indent, 0]))
        mBC = utilitiesMech.mechanicalBC(dim, 1, mechBC)
        mechBC_merged.append(mBC)

        contactWidth = radius/3
        ## top support surf
        nodeA = np.array([ indent , radius-indent, -contactWidth/2])
        nodeB = np.array([ height-indent , radius-indent, contactWidth/2])
        #pointGenerators.generateNodesOrtoSurface3dRand(nodeA, nodeB, minDist, dim, node_coords, trials)

        #rigid plate top support
        indentRP = indent
        topRigidPlateMechBC = np.array([0,1,-1, 0,0,0,  -1,-1,-1,-1,-1,-1])
        topRigidPlate = utilitiesMech.RigidPlate(-1, 3,
        np.array([ -indentRP,
         height+indentRP,
         radius*0.95,
         radius-3*indent,
         -contactWidth/2-indentRP,
          +contactWidth/2+indentRP
                ]))
        rigidPlates.append(topRigidPlate)
        govNodes.append(np.array([ height/2, radius, 0]))
        govNodesMechBC.append(utilitiesMech.mechanicalBC(dim, -1, topRigidPlateMechBC))

        ## bottom support surf
        nodeA = np.array([ indent , -radius+indent, -contactWidth/2])
        nodeB = np.array([ height-indent , -radius+indent, contactWidth/2])
        #pointGenerators.generateNodesOrtoSurface3dRand(nodeA, nodeB, minDist, dim, node_coords, trials)

        #rigid plate bottom  support
        indentRP = indent
        bottomRigidPlateMechBC = np.array([0,0,-1, 0,0,0,  -1,-1,-1,-1,-1,-1])
        bottomRigidPlate = utilitiesMech.RigidPlate(-2, 3,
        np.array([ -indentRP,
         height+indentRP,
         -radius+3*indent,
         -radius*0.95,
         -contactWidth/2-indentRP,
          +contactWidth/2+indentRP
                ]))
        rigidPlates.append(bottomRigidPlate)
        govNodes.append(np.array([ height/2, -radius, 0]))
        govNodesMechBC.append(utilitiesMech.mechanicalBC(dim, -2, bottomRigidPlateMechBC))

    print (len(node_coords))

    ###############generating of points left face ###############
    #pointGenerators.generateNodesOrtoCircleBorder3dRand(center, radius-1e-5, directionDim, minDist, node_coords, trials)
    if powerTes == False:
        pointGenerators.generateNodesOrtoCircle3dRand(center, radius-1e-5, directionDim, minDist, node_coords, trials)
    else:
        pointGenerators.generateNodesOrtoCircle3dRand(center, radius-1e-5, directionDim, minDist*0.4, node_coords, trials)
    print (len(node_coords))

    ###############generating of points  right face ###############
    nodeA = center.copy()
    nodeA[directionDim] += height-indent
    #pointGenerators.generateNodesOrtoCircleBorder3dRand(nodeA, radius-1e-5,  directionDim, minDist, node_coords, trials)
    if powerTes == False:
        pointGenerators.generateNodesOrtoCircle3dRand(center, radius-1e-5, directionDim, minDist, node_coords, trials)
    else:
        pointGenerators.generateNodesOrtoCircle3dRand(center, radius-1e-5, directionDim, minDist*0.4, node_coords, trials)
    print (len(node_coords))
    ###############generating of points cylinder surf###############
    #pointGenerators.generateNodesOrtoCilinderSurf3dRand(center, radius-1e-5, height, directionDim, minDist,  node_coords, trials)

    ###############generating of points cylinder volume ###############
    if powerTes==False:
        pointGenerators.generateNodesOrtoCilinder3dRand(center, radius-1e-5, height, directionDim, minDist,  node_coords, trials)
        radii = np.zeros(len(node_coords))
    else:
        radii = np.zeros(len(node_coords))
        node_coords = np.array(node_coords)
        node_coords, radii=pointGenerators.generateParticlesOrtoCilinder3dRand(center, radius-1e-5, height, directionDim, minDist*0.4, minDist, 0.8,  trials, node_coords, radii)
    #######################################################################
    print (len(node_coords))

    """
    node_coords = np.asarray(node_coords)
    if SHOW_PLOT:
        fig = plt.figure()
        ax = Axes3D(fig)
        ax.scatter(node_coords[:,0], node_coords[:,1], node_coords[:,2])
        plt.show()
    #"""


    return node_coords, mechBC_merged, mechInitC_merged, govNodes, govNodesMechBC, rigidPlates, radii

def assemble3dDam(maxLim, minDist, trials, topsize):
    node_coords = np.zeros((0,3))
    radii = np.zeros(0)
    mechBC_merged = []
    mechIC_merged = []

    ##########################################generating of points, homogeneous volume
    rectBC = np.array([-1,-1,-1,-1,-1,-1])
    #rect
    oldLen = len(node_coords)
    #pointGenerators.generateNodesRect(maxLim, minDist, dim, trials, node_coords)
    node_coords, radii = pointGenerators.generateParticlesDam(maxLim, topsize, minDist/4., minDist, 0.8, 3, trials, node_coords, radii)
    #
    newLen = len(node_coords)-1

   # print (nrOfPoints)
  #  mBC = utilitiesGeom.mechanicalBC(dim, kvadrBC, oldLen, newLen)
   # mechBC_merged.append(mBC)
    ####################################################################################################

    return node_coords, radii, mechBC_merged, mechIC_merged

def assemble3dcylinderUniPressConfined(center, radius, height, minDist, trials, directionDim):
    indent = 1e-5
    dim=3
    #lists for the model
    node_coords = []
    mechBC_merged = []
    mechInitC_merged = []

    node_coords.append( center+indent)

    ###############generating of points supported surface bottom face ###############
    mechBC = np.array([0,0,0,0,0,0,    -1,-1,-1,-1,-1,-1])



    oldLen = len(node_coords)
    pointGenerators.generateNodesOrtoCircleBorder3dRand(center, radius, directionDim, minDist, node_coords, trials)
    pointGenerators.generateNodesOrtoCircle3dRand(center, radius, directionDim, minDist, node_coords, trials)
    nrOfPoints =  (len(node_coords)) - oldLen
    for n in range ( nrOfPoints ):
        mBC = utilitiesMech.mechanicalBC(dim, oldLen + n, mechBC)
        mechBC_merged.append(mBC)
    #print('%d nodes generated so far' %len(node_coords))
    ###############generating of points loaded surface top face ###############
    mechBC = np.array([1,0,0,0,0,0,    -1,-1,-1,-1,-1,-1])

    nodeA = center.copy()
    nodeA[directionDim] += height

    oldLen = len(node_coords)
    pointGenerators.generateNodesOrtoCircleBorder3dRand(nodeA, radius,  directionDim, minDist, node_coords, trials)
    pointGenerators.generateNodesOrtoCircle3dRand(nodeA, radius,  directionDim, minDist, node_coords, trials)
    nrOfPoints =  (len(node_coords)) - oldLen
    for n in range ( nrOfPoints ):
        mBC = utilitiesMech.mechanicalBC(dim, oldLen + n, mechBC)
        mechBC_merged.append(mBC)


    ###############generating of points rectangular volume ###############
    mechBC = np.array([-1,0,0,0,0,0,    -1,-1,-1,-1,-1,-1])
    oldLen = len(node_coords)
    pointGenerators.generateNodesOrtoCilinderSurf3dRand(center, radius-1e-5, height, directionDim, minDist,  node_coords, trials)
    nrOfPoints =  (len(node_coords)) - oldLen
    for n in range ( nrOfPoints ):
        mBC = utilitiesMech.mechanicalBC(dim, oldLen + n, mechBC)
        mechBC_merged.append(mBC)
    #######################################################################

    ###############generating of points rectangular volume ###############
    pointGenerators.generateNodesOrtoCilinder3dRand(center, radius, height, directionDim, minDist,  node_coords, trials)
    #######################################################################


    return node_coords, mechBC_merged, mechInitC_merged





def assemble3dtubeTorsionFree(center, radius, height, thickness, minDist, trials, directionDim, functions, rotationAngle):
    indent = 1e-5
    dim=3
    #lists for the model
    node_coords = []
    mechBC_merged = []
    mechInitC_merged = []

    mechBC = np.array([0,0,0,0,0,0,    -1,-1,-1,-1,-1,-1])
    node_coords.append( np.array([  0,  radius-thickness/2,  0 ]))
    mBC = utilitiesMech.mechanicalBC(dim, 0, mechBC)
    mechBC_merged.append(mBC)

    ###############generating of points supported surface left face ###############
    mechBC = np.array([0,0,0, 0,0,0,    -1,-1,-1,-1,-1,-1])

    oldLen = len(node_coords)
    pointGenerators.generateNodesOrtoCircleBorder3dRand(center, radius, directionDim, minDist, node_coords, trials)
    pointGenerators.generateNodesOrtoCircleBorder3dRand(center, radius-thickness, directionDim, minDist, node_coords, trials)
    pointGenerators.generateNodesOrtoAnnulus3dRand(center, radius, thickness, directionDim, minDist, node_coords, trials)
    nrOfPoints =  (len(node_coords)) - oldLen
    for n in range ( nrOfPoints ):
        mBC = utilitiesMech.mechanicalBC(dim, oldLen + n, mechBC)
        mechBC_merged.append(mBC)
    #print('%d nodes generated so far' %len(node_coords))
    ###############generating of points loaded surface right face ###############

    nodeA = center.copy()
    nodeA[directionDim] += float(height)

    oldLen = len(node_coords)
    pointGenerators.generateNodesOrtoCircleBorder3dRand(nodeA, radius, directionDim, minDist, node_coords, trials)
    pointGenerators.generateNodesOrtoCircleBorder3dRand(nodeA, radius-thickness, directionDim, minDist, node_coords, trials)
    pointGenerators.generateNodesOrtoAnnulus3dRand(nodeA, radius, thickness, directionDim, minDist, node_coords, trials)
    nrOfPoints =  (len(node_coords)) - oldLen
    for n in range ( nrOfPoints ):
        lnfc = len(functions)
        mechBC = np.array([lnfc, lnfc+1, lnfc+2, -1 , -1, -1,    -1,-1,-1,-1,-1,-1])
        point = node_coords[oldLen + n]
        rotAngles = np.array([rotationAngle, 0, 0])
        value = 1
        period = 1

        funcRot0 = utilitiesNumeric.constSawToothRotationFunction(rotAngles, point, value, 0, period=period, sym = 1)
        funcRot1 = utilitiesNumeric.constSawToothRotationFunction(rotAngles, point, value, 1, period=period, sym = 1)
        funcRot2 = utilitiesNumeric.constSawToothRotationFunction(rotAngles, point, value, 2, period=period, sym = 1)
        #//funcRot3 = utilitiesNumeric.constSawToothRotationFunction(rotAngles, point, value, 3, period=period)
        #//funcRot4 = utilitiesNumeric.constSawToothRotationFunction(rotAngles, point, value, 4, period=period)
        #//funcRot5 = utilitiesNumeric.constSawToothRotationFunction(rotAngles, point, value, 5, period=period)
        functions.append(funcRot0)
        functions.append(funcRot1)
        functions.append(funcRot2)
        #//functions.append(funcRot3)
        #//functions.append(funcRot4)
        #//functions.append(funcRot5)

        mBC = utilitiesMech.mechanicalBC(dim, oldLen + n, mechBC)
        mechBC_merged.append(mBC)


    ###############generating of points rectangular volume ###############
    pointGenerators.generateNodesOrtoCilinderSurf3dRand(center, radius-1e-5, height, directionDim, minDist,  node_coords, trials)
    pointGenerators.generateNodesOrtoCilinderSurf3dRand(center, radius-thickness+1e-5, height, directionDim, minDist,  node_coords, trials)
    #######################################################################

    ###############generating of points rectangular volume ###############
    pointGenerators.generateNodesOrtoTube3dRand(center, radius, height, thickness, directionDim, minDist,  node_coords, trials)
    #######################################################################


    return node_coords, mechBC_merged, mechInitC_merged




def assemble3dBiparvaTubeTransport(center, radius, height, thickness, minDist, trials):
    print ('Assembling Biparva tube...')
    directionDim = 0
    indent = 1e-5
    dim=3
    #lists for the model
    node_coords = []
    mechBC_merged = []
    mechInitC_merged = []
    govNodes = []
    govNodesMechBC = []
    rigidPlates = []

    center[0] += 1e-7
    #mechBC = np.array([0,0,0,0,0,0,    -1,-1,-1,-1,-1,-1])
    node_coords.append( np.array([  indent,  radius-thickness/2,  0 ]))
    #mBC = utilitiesMech.mechanicalBC(dim, 0, mechBC)
    #mechBC_merged.append(mBC)

    ###############generating of points supported surface left face ###############
    #pointGenerators.generateNodesOrtoCircleBorder3dRand(center, radius, directionDim, minDist, node_coords, trials)
    #pointGenerators.generateNodesOrtoCircleBorder3dRand(center, radius-thickness, directionDim, minDist, node_coords, trials)
    pointGenerators.generateNodesOrtoAnnulus3dRand(center, radius, thickness, directionDim, minDist, node_coords, trials)

    indentRP = 1e-6
    leftRigidPlateMechBC = np.array([0, 0,0, 0,-1,-1,  -1,-1,-1,   -1,-1,-1])
    leftRigidPlate = utilitiesMech.RigidPlate(-1, 3, np.array([
    -indentRP, indentRP,
    -radius-indentRP, radius+indentRP,
    -radius-indentRP, radius+indentRP, ]))
    rigidPlates.append(leftRigidPlate)
    govNodes.append(np.array([ 0,0,0 ]))
    govNodesMechBC.append(utilitiesMech.mechanicalBC(dim, -1, leftRigidPlateMechBC))

    center[0] -= 2e-7
    nodeA = center.copy()
    nodeA[directionDim] += float(height)

    #pointGenerators.generateNodesOrtoCircleBorder3dRand(nodeA, radius, directionDim, minDist, node_coords, trials)
    #pointGenerators.generateNodesOrtoCircleBorder3dRand(nodeA, radius-thickness, directionDim, minDist, node_coords, trials)
    pointGenerators.generateNodesOrtoAnnulus3dRand(nodeA, radius, thickness, directionDim, minDist, node_coords, trials)

    rightRigidPlateMechBC = np.array([1, 0,0, -1,-1,-1,  -1,-1,-1,   -1,-1,-1])
    rightRigidPlate = utilitiesMech.RigidPlate(-1, 3, np.array([
    -indentRP+height, indentRP+height,
    -radius-indentRP, radius+indentRP,
    -radius-indentRP, radius+indentRP, ]))
    rigidPlates.append(rightRigidPlate)
    govNodes.append(np.array([ height ,0,0 ]))
    govNodesMechBC.append(utilitiesMech.mechanicalBC(dim, -2, rightRigidPlateMechBC))

    ###############generating of points rectangular volume ###############
    #pointGenerators.generateNodesOrtoCilinderSurf3dRand(center, radius-1e-5, height, directionDim, minDist,  node_coords, trials)

    #pointGenerators.generateNodesOrtoCilinderSurf3dRand(center, radius-thickness+1e-5, height, directionDim, minDist,  node_coords, trials)
    #######################################################################

    ###############generating of points rectangular volume ###############
    pointGenerators.generateNodesOrtoTube3dRand(np.zeros((3)), radius-1e-5, height, thickness, directionDim, minDist,  node_coords, trials)
    #######################################################################


    return node_coords, mechBC_merged,  govNodes, govNodesMechBC, rigidPlates




def assemble3dTubeSplit(center, radius, height, thickness, minDist, trials,notchH=0, notchWidth=0):
    print ('Assembling tubesplit...')
    directionDim = 0
    indent = 1e-5
    dim=3
    #lists for the model
    node_coords = []
    mechBC_merged = []
    mechInitC_merged = []
    govNodes = []
    govNodesMechBC = []
    rigidPlates = []

    center[0] += 1e-7
    #mechBC = np.array([0,0,0,0,0,0,    -1,-1,-1,-1,-1,-1])
    node_coords.append( np.array([  indent,  radius-thickness/2,  0 ]))
    #mBC = utilitiesMech.mechanicalBC(dim, 0, mechBC)
    #mechBC_merged.append(mBC)

    ###############generating of points supported surface left face ###############
    #pointGenerators.generateNodesOrtoAnnulus3dRand(center, radius, thickness, directionDim, minDist, node_coords, trials)

    center[0] -= 2e-7
    nodeA = center.copy()
    nodeA[directionDim] += float(height)

    #pointGenerators.generateNodesOrtoAnnulus3dRand(nodeA, radius, thickness, directionDim, minDist, node_coords, trials)
    notches = []
    if notchH>0:
        print('notch',notchH,notchWidth)
        oldLen = len(node_coords)

        notchMinDist=minDist*0.9

        notchSide0 = []
        nodeA = np.array([indent, notchWidth/2, radius-thickness])
        nodeB = np.array([indent, notchWidth/2, radius-thickness+thickness*notchH])
        pointGenerators.generateNodesLine3dRand(nodeA, nodeB, notchMinDist, dim, node_coords, trials, catchCorners=True, equidist=True)

        nodeA = np.array([height-indent, notchWidth/2, radius-thickness])
        nodeB = np.array([height-indent, notchWidth/2, radius-thickness+thickness*notchH])
        pointGenerators.generateNodesLine3dRand(nodeA, nodeB, notchMinDist, dim, node_coords, trials, catchCorners=True, equidist=True)

        nodeA = np.array([indent, notchWidth/2, radius-thickness])
        nodeB = np.array([height-indent, notchWidth/2, radius-thickness])
        pointGenerators.generateNodesLine3dRand(nodeA, nodeB, notchMinDist, dim, node_coords, trials, catchCorners=False, equidist=True)

        nodeA = np.array([indent, notchWidth/2, radius-thickness+thickness*notchH])
        nodeB = np.array([height-indent, notchWidth/2, radius-thickness+thickness*notchH])
        pointGenerators.generateNodesLine3dRand(nodeA, nodeB, notchMinDist, dim, node_coords, trials, catchCorners=False, equidist=True)

        nodeA = np.array([indent, notchWidth/2, radius-thickness+thickness*notchH])
        nodeB = np.array([height-indent, notchWidth/2, radius-thickness])
        pointGenerators.generateNodesOrtoSurface3dRand(nodeA, nodeB, notchMinDist, dim, node_coords, trials)

        for i in range (oldLen, len(node_coords), 1):
            newNode = np.copy(node_coords[i])
            newNode[2] *=-1
            node_coords.append(newNode)

        for i in range (oldLen, len(node_coords), 1):
            notchSide0.append(i)

        notchSide1 = []

        for i in range (oldLen, len(node_coords), 1):
            newNode = np.copy(node_coords[i])
            newNode[1] -= notchWidth
            node_coords.append(newNode)
            notchSide1.append(len(node_coords)-1)


        notchA = []
        notchA.append(notchSide0)
        notchA.append(notchSide1)
        notches.append(notchA)


    ###############generating of points rectangular volume ###############
    pointGenerators.generateNodesOrtoCilinderSurf3dRand(center, radius-1e-5, height, directionDim, minDist,  node_coords, trials)

    oldlen = len(node_coords)
    pointGenerators.generateNodesOrtoCilinderSurf3dRand(center, radius-thickness+1e-5, height, directionDim, minDist,  node_coords, trials)
    #######################################################################
    newNodes = len(node_coords)-oldlen

    #select nodes in top third
    idcs_top = []
    for n in range(newNodes):
        h = oldlen+n
        if (node_coords[h][1] > (center[1]+0.7071*(radius-thickness))) :
            idcs_top.append(h)

    leftRigidPlateMechBC = np.array([-1, 1,-1,   -1,-1,-1, -1,-1,-1 ,-1,-1,-1])
    leftRigidPlate = utilitiesMech.RigidPlate(-1, 3, None, directIdcs=True )
    leftRigidPlate.setDirectNodes(idcs_top)
    rigidPlates.append(leftRigidPlate)
    govNodes.append(np.array([ height/2, center[1]+(radius-thickness),0 ]))
    govNodesMechBC.append(utilitiesMech.mechanicalBC(dim, -1, leftRigidPlateMechBC))

    print('Idcs top ' , len(idcs_top))
    #select nodes in top third
    idcs_top = []
    for n in range(newNodes):
        h = oldlen+n
        if (node_coords[h][1] < (center[1]-0.7071*(radius-thickness))) :
            idcs_top.append(h)
    print('Idcs bot ' , len(idcs_top))
    leftRigidPlateMechBC = np.array([0, 0,0,   0,0,0 , -1,-1,-1  ,-1,-1,-1])
    leftRigidPlate = utilitiesMech.RigidPlate(-2, 3, None, directIdcs=True )
    leftRigidPlate.setDirectNodes(idcs_top)
    rigidPlates.append(leftRigidPlate)
    govNodes.append(np.array([ height/2, center[1]-(radius-thickness),0]))
    govNodesMechBC.append(utilitiesMech.mechanicalBC(dim, -2, leftRigidPlateMechBC))

    ###############generating of points rectangular volume ###############
    pointGenerators.generateNodesOrtoTube3dRand(np.zeros((3)), radius-1e-5, height, thickness, directionDim, minDist,  node_coords, trials)
    #######################################################################


    return node_coords, mechBC_merged,  govNodes, govNodesMechBC, rigidPlates,notches






def assemble3dTubeInnerPressure(center, radius, height, thickness, minDist, trials, powerTes):
    print ('Assembling tube with inner pressure...')
    directionDim = 0
    indent = 1e-6
    indentRP = 1e-5
    dim=3
    #lists for the model
    node_coords = []
    mechBC_merged = []
    mechInitC_merged = []
    govNodes = []
    govNodesMechBC = []
    rigidPlates = []

    radii = []

    center[0] += indent
    node_coords.append( np.array([ height/2.,  radius,  1e-6 ]))
    node_coords.append( np.array([ height/2.,  -radius,  1e-6 ]))
    mechBC = np.array([0,-1,-1,  0,0,0,    -1,-1,-1,-1,-1,-1])
    mBC = utilitiesMech.mechanicalBC(dim, 0, mechBC)
    mechBC_merged.append(mBC)
    mBC = utilitiesMech.mechanicalBC(dim, 1, mechBC)
    mechBC_merged.append(mBC)

    circleLength = 2*np.pi*(radius-thickness)
    if powerTes:
        nrNodes = int ( circleLength / (minDist*0.4) )
        nrNodes = (int (nrNodes / 4) +1 ) * 4
        pointGenerators.generateNodesOrtoCilinderSurf3dRand(center, radius-thickness, height, directionDim, minDist*0.4, node_coords, trials, equiAngNodes=nrNodes )
    else:
        nrNodes = int ( circleLength / minDist )
        nrNodes = (int (nrNodes / 4) +1 ) * 4
        pointGenerators.generateNodesOrtoCilinderSurf3dRand(center, radius-thickness, height, directionDim, minDist,  node_coords, trials, equiAngNodes=nrNodes )

    borderRing_1 = []
    borderRing_2 = []
    for i in range (len(node_coords)):
        if (np.linalg.norm(node_coords[i][1:3]  -center[1:3] ) < radius-thickness+1e-3  ):
            if ( node_coords[i][0] < minDist*2/3  ):
                borderRing_1.append(i)
            if ( node_coords[i][0] > height-minDist*2/3   ):
                borderRing_2.append(i)



    leftRigidPlateMechBC = np.array([0, -1,-1,   -1,-1,-1 , -1,-1,-1,   -1,-1,-1])
    borderRing_RP1 = utilitiesMech.RigidPlate(-1, 1, None, directIdcs = True, dirs='6 1 0 0 0 0 0')
    borderRing_RP1.setDirectNodes(borderRing_1)
    rigidPlates.append(borderRing_RP1)
    govNodes.append(np.array([ 1e-4,0,0 ]))
    govNodesMechBC.append(utilitiesMech.mechanicalBC(dim, -1, np.copy(leftRigidPlateMechBC)) )

    borderRing_RP2 = utilitiesMech.RigidPlate(-2, 1, None, directIdcs = True, dirs='6 1 0 0 0 0 0')
    borderRing_RP2.setDirectNodes(borderRing_2)
    rigidPlates.append(borderRing_RP2)
    govNodes.append(np.array([ height-1e-4,0,0 ]))
    govNodesMechBC.append(utilitiesMech.mechanicalBC(dim, -2, np.copy(leftRigidPlateMechBC)) )


    #pointGenerators.generateNodesOrtoAnnulus3dRand(center, radius, thickness, directionDim, minDist, node_coords, trials)


    leftRigidPlateMechBC = np.array([0, -1,-1,   0,0,0,  -1,-1,-1,   -1,-1,-1])
    leftRigidPlate = utilitiesMech.RigidPlate(-3, 3, np.array([
    -indentRP, indentRP,
    -radius-indentRP, radius+indentRP,
    -radius-indentRP, radius+indentRP, ]), dirs='6 1 0 0 0 0 0')
    rigidPlates.append(leftRigidPlate)
    govNodes.append(np.array([ 0,0,0 ]))
    govNodesMechBC.append(utilitiesMech.mechanicalBC(dim, -3, np.copy(leftRigidPlateMechBC)) )


    center[0] -= 2*indent
    nodeA = center.copy()
    nodeA[directionDim] += float(height)

    #pointGenerators.generateNodesOrtoAnnulus3dRand(nodeA, radius, thickness, directionDim, minDist, node_coords, trials)

    rightRigidPlateMechBC = np.array([0, -1, -1, 0,0,0, -1,-1,-1,   -1,-1,-1])
    rightRigidPlate = utilitiesMech.RigidPlate(-4, 3, np.array([
    -indentRP+height, indentRP+height,
    -radius-indentRP, radius+indentRP,
    -radius-indentRP, radius+indentRP, ]), dirs='6 1 0 0 0 0 0')
    rigidPlates.append(rightRigidPlate)
    govNodes.append(np.array([ height ,0,0 ]))
    govNodesMechBC.append(utilitiesMech.mechanicalBC(dim, -4, rightRigidPlateMechBC))

    ###############generating of points rectangular volume ###############
    if powerTes:
        radii = np.zeros(len(node_coords))
        node_coords = np.array(node_coords)
        node_coords, radii = pointGenerators.generateParticlesOrtoTube3dRand(np.zeros((3)), radius-1e-5, height, thickness, directionDim, minDist*0.4, minDist, 0.8, node_coords, radii, trials)
    else:
        pointGenerators.generateNodesOrtoTube3dRand(np.zeros((3)), radius-1e-5, height, thickness, directionDim, minDist,  node_coords, trials)
        radii = np.zeros(len(node_coords))


    rebarBC = np.array([2, 0, -1,    -1, -1, -1,  -1,-1,-1,   -1,-1,-1])
    govNodesMechBC.append(utilitiesMech.mechanicalBC(dim, (-5), rebarBC))

    return node_coords, mechBC_merged,  govNodes, govNodesMechBC, rigidPlates, radii


def assemble3dslimTubeTorsionFree(center, radius, height, thickness, minDist, trials, directionDim, functions):
    indent = 1e-5
    dim=3
    #lists for the model
    node_coords = []
    mechBC_merged = []
    mechInitC_merged = []

    mechBC = np.array([0,0,0,0,0,0,    -1,-1,-1,-1,-1,-1])
    node_coords.append( np.array([  0,  radius-thickness/2,  0 ]))
    mBC = utilitiesMech.mechanicalBC(dim, 0, mechBC)
    mechBC_merged.append(mBC)


    ###############generating of points supported surface left face ###############
    mechBC = np.array([0,0,0, 0 , 0, 0,    -1,-1,-1,-1,-1,-1])
    oldLen = len(node_coords)
    pointGenerators.generateNodesOrtoCircleBorder3dRand(center, radius, directionDim, minDist, node_coords, trials)
    pointGenerators.generateNodesOrtoCircleBorder3dRand(center, radius-thickness, directionDim, minDist, node_coords, trials)
    nrOfPoints =  (len(node_coords)) - oldLen
    for n in range ( nrOfPoints ):
        mBC = utilitiesMech.mechanicalBC(dim, oldLen + n, mechBC)
        mechBC_merged.append(mBC)


    ###############generating of points loaded surface right face ###############
    nodeA = center.copy()
    nodeA[directionDim] += height
    oldLen = len(node_coords)
    pointGenerators.generateNodesOrtoCircleBorder3dRand(nodeA, radius, directionDim, minDist, node_coords, trials)
    pointGenerators.generateNodesOrtoCircleBorder3dRand(nodeA, radius-thickness, directionDim, minDist, node_coords, trials)
    nrOfPoints =  (len(node_coords)) - oldLen
    for n in range ( nrOfPoints ):
        lnfc = len(functions)
        mechBC = np.array([lnfc, lnfc+1, lnfc+2, -1 , -1, -1,    -1,-1,-1,-1,-1,-1])
        point = node_coords[oldLen + n]
        rotAngles = np.array([0.01, 0, 0])
        value = 1
        period = 1

        funcRot0 = utilitiesNumeric.constSawToothRotationFunction(rotAngles, point, value, 0, period=period, sym = 1)
        funcRot1 = utilitiesNumeric.constSawToothRotationFunction(rotAngles, point, value, 1, period=period, sym = 1)
        funcRot2 = utilitiesNumeric.constSawToothRotationFunction(rotAngles, point, value, 2, period=period, sym = 1)
        #//funcRot3 = utilitiesNumeric.constSawToothRotationFunction(rotAngles, point, value, 3, period=period)
        #//funcRot4 = utilitiesNumeric.constSawToothRotationFunction(rotAngles, point, value, 4, period=period)
        #//funcRot5 = utilitiesNumeric.constSawToothRotationFunction(rotAngles, point, value, 5, period=period)
        functions.append(funcRot0)
        functions.append(funcRot1)
        functions.append(funcRot2)
        #//functions.append(funcRot3)
        #//functions.append(funcRot4)
        #//functions.append(funcRot5)
        mBC = utilitiesMech.mechanicalBC(dim, oldLen + n, mechBC)
        mechBC_merged.append(mBC)


    ###############generating of points rectangular volume ###############
    pointGenerators.generateNodesOrtoCilinderSurf3dRand(center, radius, height, directionDim, minDist,  node_coords, trials)
    pointGenerators.generateNodesOrtoCilinderSurf3dRand(center, radius-thickness, height, directionDim, minDist,  node_coords, trials)
    #######################################################################

    return node_coords, mechBC_merged, mechInitC_merged
