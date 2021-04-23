import numpy as np
import random
import utilitiesGeom
import utilitiesMech
import utilitiesNumeric
import pointGenerators
import voronoi
import math
import matplotlib.pyplot as plt
#import voronoi_viewer
from mpl_toolkits.mplot3d import Axes3D

from scipy.spatial import Voronoi
from scipy.spatial import voronoi_plot_2d
from scipy.spatial import Delaunay
import tkinter
SHOW_PLOT = False


def assembleMeasuringGauges(type, D=-1, maxLim = None):
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

    if(type=='3pb2d'):
        coordsA = np.array([ maxLim[0]/2, maxLim[1]])
        coordsB = np.array([ 0, 0 ])
        measuringGauges.append(utilitiesMech.MeasuringGauge(coordsA, coordsB, 'totalLS', False))
        coordsA = np.array([ maxLim[0]/2-1e3, 0])
        coordsA = np.array([ maxLim[0]/2+1e3, 0])
        measuringGauges.append(utilitiesMech.MeasuringGauge(coordsA, coordsB, 'notchOpening', False))

    if(type=='3pb3d'):
        coordsA = np.array([ maxLim[0]/2, maxLim[1], maxLim[2]/2])
        coordsB = np.array([ 0, 0, 0 ])
        measuringGauges.append(utilitiesMech.MeasuringGauge(coordsA, coordsB, 'totalTopMid', False))
        coordsA = np.array([ maxLim[0]/2, 0, maxLim[2]/2])
        coordsB = np.array([ 0, 0, 0 ])
        measuringGauges.append(utilitiesMech.MeasuringGauge(coordsA, coordsB, 'totalBotMid', False))
        coordsA = np.array([ maxLim[0]/2-1e3, 0, maxLim[2]/2])
        coordsA = np.array([ maxLim[0]/2+1e3, 0, maxLim[2]/2])
        measuringGauges.append(utilitiesMech.MeasuringGauge(coordsA, coordsB, 'notchOpening', False))

    if(type=='3d_brazilianDisc'):
        #total length
        coordsA = np.array([ maxLim[0]/2, maxLim[1]/2, 0])
        coordsB = np.array([ maxLim[0]/2, -maxLim[1]/2, 0 ])
        measuringGauges.append(utilitiesMech.MeasuringGauge(coordsA, coordsB, 'totalLoadDir', False))
        coordsA = np.array([ maxLim[0]/2, 0, maxLim[2]/2])
        coordsB = np.array([ maxLim[0]/2, 0, -maxLim[2]/2 ])
        measuringGauges.append(utilitiesMech.MeasuringGauge(coordsA, coordsB, 'totalAcrossDir', False))

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


def assembleMaterialZones (elaX, dim, model='box', maxLim=None, D=None, thickness=None, limits=None, limits1=None, rebarDepth=None, rebarDiameter=None, rebarCount=None):
    #limits = xmin, ymin, zmin, xmax, ymax, zmax
    materialZones = []
    #matZone 1
    matZ = []
    #matZone 2
    matZ1 = []
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

    if (model =='3pb3d'):
        boundA = np.array(  [ limits[0]    ,  limits[1],  limits[2]  ] )
        matZ.append (boundA)
        boundB = np.array(  [  limits[3]   ,   limits[4],  limits[5] ] )
        matZ.append (boundB)
        materialZones.append(matZ)

        boundA1 = np.array(  [ limits1[0]    ,  limits1[1],  limits1[2]  ] )
        boundB1 = np.array(  [ limits1[3]    ,  limits1[4],  limits1[5]  ] )
        matZ1.append(boundA1)
        matZ1.append(boundB1)
        materialZones.append(matZ1)

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
            radius = rebarDiameter/2 * 1.01
            center = np.array([ (maxLim[0]/rebarCount)*(r+0.5), maxLim[1]-rebarDepth  ])
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

    print(node_coords)
    node_coords = np.asarray(node_coords)
    # if SHOW_PLOT:
    #     fig, ax = plt.subplots()
    #     ax.scatter(node_coords[:,0], node_coords[:,1])
    #     plt.show()

    return node_coords, mechBC_merged,  vor, areas, functions,govNodes, govNodesMechBC, rigidPlates



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

    print(node_coords)

    print('Conducting Voronoi tesselation...', end = '')
    vor, regions, vertices, polygons, areas, centroids, points = utilitiesNumeric.runMirroredVoronoi (node_coords, 2, maxLim, shifts=shifts)
    print('done.')

    #print(vor.points)
    print(areas)

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

    print(allVrtcs)
    for i in range (len(allVrtcs)):
        trsBC = utilitiesMech.transportBC(allVrtcs[i], noTrsprtBC)
        transportBC_merged.append(trsBC)


    return node_coords, mechBC_merged, transportBC_merged, vor, areas, functions





def create2dSSBeamUnifLoad(maxLim, minDist, trials, notch = -1,
                           loadWidth = 1, fracZoneWidth = 0.15,
                           orthogonalFracZone=False, notchWidth =-1,
                           node_coords_init=None,
                           activeTransport=False,
                           coupled = False):
    print('Creating 2d simply supported beam, uniform load.')
    #
    node_coords, mechBC_merged, mechInitC_merged, notches, govNodes, govNodesMechBC, rigidPlates  = assemble2DSSBeamBending(maxLim, minDist, trials, notch, loadWidth, fracZoneWidth, orthogonalFracZone=orthogonalFracZone, notchWidth=notchWidth, node_coords_init=node_coords_init,  coupled=coupled);

    print('Conducting Voronoi tesselation...', end = '')
    vor, regions, vertices, polygons, areas, centroids, points = utilitiesNumeric.runMirroredVoronoi (node_coords, 2, maxLim)
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
    func1.append( np.array([1, -1e-2]) )
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


    #return node_coords, mechBC_merged, transportBC_merged,  notches, govNodes, govNodesMechBC, rigidPlates, vor, areas, functions


    return node_coords, mechBC_merged, [], vor, areas, functions, notches, govNodes,    govNodesMechBC, rigidPlates, transportBC_merged




def create2dCantileverBending(maxLim, minDist, trials ):
    print('Creating 2d cantilever, bending.')
    ### sampling of nodes
    ### direct setting of mechanicalBCs
    node_coords, mechBC_merged, mechIC_merged  = assemble2DCantileverBending(maxLim, minDist, trials );

    print('Conducting Voronoi tesselation...', end = '')
    vor, regions, vertices, polygons, areas, centroids, points = utilitiesNumeric.runMirroredVoronoi (node_coords, 2, maxLim)
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

    return node_coords, mechBC_merged, mechIC_merged, transportBC_merged, transportIC_merged, vor, areas, functions


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







def create2dCorrosionRebar(maxLim, minDist, trials, rebarMinDist, rebarDiameter, rebarCount, rebarDepth, node_coords_init=None):
    print('Creating corrosion rebar model...')
    dim=2

    ### sampling of nodes
    ### direct setting of mechanicalBCs
    sampleBorders = True
    node_coords, mechBC_merged, mechInitC_merged,  govNodes, govNodesMechBC, rigidPlates  = assemble2dCorrosionRebar(maxLim, minDist, trials, rebarMinDist, rebarDiameter, rebarCount, rebarDepth, sampleBorders, node_coords_init=node_coords_init)




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

    fn1 = utilitiesNumeric.constantFunc(100)
    functions.append (fn1)

    trsptLeftRigidPlate = utilitiesMech.RigidPlate(-len(govNodesTrspt)-1, 3, None, directIdcs = True)
    trsptLeftRigidPlate.setDirectNodes(faces1)
    trsptLeftRigidPlateMechBC = np.array([0,-1])
    rigidPlatesTrspt.append(trsptLeftRigidPlate)
    govNodesTrspt.append(np.array([ 0, 0, 0]))
    govNodesTrsptBC.append(utilitiesMech.transportBC(govNodesTrspt[-1], trsptLeftRigidPlateMechBC))


    #"""
    ### selecting vertices on the top surface
    boundA = np.array([-1e-8, maxLim[1]-1e-8])
    boundB = np.array([maxLim[0], maxLim[1]+1e-8])
    faces2 = utilitiesGeom.returnSelectedPts(boundA, boundB, vor.vertices)
    vert = vor.vertices[faces1,:]

    trsptRightRigidPlate = utilitiesMech.RigidPlate(-len(govNodesTrspt)-1, 2, None, directIdcs = True)
    trsptRightRigidPlate.setDirectNodes(faces2)
    trsptRightRigidPlateMechBC = np.array([3,-1])
    rigidPlatesTrspt.append(trsptRightRigidPlate)
    govNodesTrspt.append(np.array([ -1, -1, -1]))
    govNodesTrsptBC.append(utilitiesMech.transportBC(govNodesTrspt[-1], trsptRightRigidPlateMechBC))




    return node_coords, mechBC_merged, transportBC_merged, govNodes, govNodesMechBC, rigidPlates, vor, areas, functions, rigidPlatesTrspt, govNodesTrspt, govNodesTrsptBC





def createCoupledBrazilianDisc(center, cylinderRad, cylinderHeight,  minDist, trials):
    dim = 3
    print('Creating coupled brazilian disc...')
    ### sampling of nodes
    ### direct setting of mechanicalBCs
    node_coords, mechBC_merged, mechInitC_merged, govNodes, govNodesMechBC, rigidPlates  = assembleCoupledBrazilianDisc(center, cylinderRad, cylinderHeight, minDist, trials, 0, [] )
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
    vor, volumes = utilitiesNumeric.runCylinderMirroredVoronoi (node_coords, center, cylinderRad, cylinderHeight, 0, quarter=False)
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


    fn3 = utilitiesNumeric.constantFunc(100)
    functions.append (fn3)


    """
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
    """

    #"""
    ### indirect setting of transportBCs by spatial selection of vertices
    transportBC_merged = []
    transportIC_merged = []
    ### selecting vertices on the left surface
    leftFaceBC = np.array([0,-1])
    boundA = np.array(  [-1e-4 , -cylinderRad*2, -cylinderRad*2] )
    boundB = np.array(  [ 1e-4 , cylinderRad*2, cylinderRad*2]  )
    leftFace = utilitiesGeom.returnSelectedPts(boundA, boundB, vor.vertices)
    """
    for i in range (len(leftFace)):
        trsBC = utilitiesMech.transportBC(leftFace[i], leftFaceBC)
        transportBC_merged.append(trsBC)
    """

    ### selecting vertices on the right surface
    rightFaceBC = np.array([2,-1])
    boundA = np.array(  [cylinderHeight-1e-4 , -cylinderRad*2, -cylinderRad*2] )
    boundB = np.array(  [cylinderHeight+1e-4 , cylinderRad*2, cylinderRad*2]  )
    rightFace = utilitiesGeom.returnSelectedPts(boundA, boundB, vor.vertices)
    """
    for i in range (len(rightFace)):
        trsBC = utilitiesMech.transportBC(rightFace[i], rightFaceBC)
        transportBC_merged.append(trsBC)
    """

    #setting of transport bc by using rigid plate
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

    return node_coords, mechBC_merged, transportBC_merged,  govNodes, govNodesMechBC, rigidPlates, vor, [], functions, rigidPlatesTrspt, govNodesTrspt, govNodesTrsptBC

def create2dPeriodicShear(maxLim, minDist, trials, powerTes ):
    print('Creating 2d periodic rectangle, shear loaded.')
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






def create2dDogBone(minDist, trials, D=1.0, excentricity = 50, symmetric=False, edgeMinDistCoef=1.0, roughDogBone=0, roughEdgeDogbone = 0, roughMinDistCoef=1, interLayerThickness=2):
    print('Creating 2d dog bone....')
    #


    node_coords_all, node_indices_dogbone, mechBC_merged, mechInitC_merged, node_count, govNodes, govNodesMechBC, rigidPlates  = assemble2dDogBone(D, minDist, trials, excentricity = excentricity, symmetric = symmetric, edgeMinDistCoef=edgeMinDistCoef, roughDogBone=roughDogBone, roughEdgeDogbone=roughEdgeDogbone, roughMinDistCoef=roughMinDistCoef, interLayerThickness=interLayerThickness);

    node_coords_all = np.asarray(node_coords_all)

    """
    fig, ax = plt.subplots()
    ax.scatter(node_coords_all[:,0], node_coords_all[:,1])
    ax.scatter(node_coords_all[node_indices_dogbone,0], node_coords_all[node_indices_dogbone,1])
    plt.show()
    #"""

    print('Conducting Voronoi tesselation...', end = '')
    vor = utilitiesNumeric.runMirroredVoronoiDogBone(node_coords_all, 2, D)
    print('done.')

    node_coords = np.copy(node_coords_all)
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


    return node_coords, mechBC_merged, mechInitC_merged, transportBC_merged, transportIC_merged, vor, areas, functions,  govNodes, govNodesMechBC, rigidPlates, node_indices_dogbone


def create3dDogBone(minDist, trials, D=1.0, excentricity = 20 ):
    print('Creating 3sd dog bone....')
    #
    node_coords, mechBC_merged, mechInitC_merged, node_count,govNodes, govNodesMechBC, rigidPlates  = assemble3dDogBone(D, minDist, trials, excentricity = excentricity);

    node_coords = np.asarray(node_coords)
    """
    if SHOW_PLOT:
        fig = plt.figure()
        ax = fig.add_subplot(111, projection='3d')
        ax.scatter(node_coords[:,0], node_coords[:,1], node_coords[:,2], c = 'b', marker='o')
        plt.show()
    """
    print('Conducting Voronoi tesselation...', end = '')
    vor = utilitiesNumeric.runMirroredVoronoiDogBone(node_coords, 3, D, thickness = 0.1)
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


    return node_coords, mechBC_merged, mechInitC_merged, transportBC_merged, transportIC_merged, vor, areas, functions, govNodes, govNodesMechBC, rigidPlates





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






def create3dSSBeamUnifLoad(maxLim, minDist, trials, notch = -1, loadWidth = 1, fracZoneWidth = 0.15, orthogonalFracZone = False, notchWidth = -1, coupled=False, node_coords_init=None ):
    print('Creating 3d simply supported beam, uniform load.')
    #govNodes, rigidPlates
    node_coords, mechBC_merged, mechInitC_merged, notches, govNodes, govNodesMechBC, rigidPlates  = assemble3DSSBeamBending(maxLim, minDist, trials, notch, loadWidth, fracZoneWidth=fracZoneWidth, orthogonalFracZone=orthogonalFracZone, notchWidth = notchWidth, coupled=coupled, node_coords_init=node_coords_init);
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

    return node_coords, mechBC_merged, mechInitC_merged,  vor, volumes, functions, notches, govNodes, govNodesMechBC, rigidPlates, transportBC_merged, transportIC_merged


def create3dCube(maxLim, minDist, trials, powerTes, coupled=False, node_coords_init=None ):
    print('Creating 3d cube. Power tesselation: %s' %powerTes)
    #govNodes, rigidPlates
    dim = 3
    node_coords, mechBC_merged, mechInitC_merged,  govNodes, govNodesMechBC, rigidPlates  = assemble3Dcube(maxLim, minDist, trials, powerTes, coupled=coupled, node_coords_init=node_coords_init);
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



def create3dReinhardtTension(maxLim, minDist, trials, fracZoneWidth = 0.15 ):
    print('Creating 3d simply supported beam, uniform load.')
    #govNodes, rigidPlates
    node_coords, mechBC_merged, mechInitC_merged, notchNodes, govNodes, govNodesMechBC, rigidPlates  = assemble3DReinhardtTension(maxLim, minDist, trials, fracZoneWidth=fracZoneWidth);
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





def create3dcylinderTorsionFree(center, radius, height, minDist, trials, directionDim ):

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


def create3dcylinderTorsionPressFree(center, radius, height, minDist, trials, directionDim, activeTransport):

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
    node_coords, mechBC_merged, govNodes, govNodesMechBC, rigidPlates  = assemble3dcylinderTorsionPressFree(center, radius, height, minDist, trials, directionDim, functions )

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

    maxLim = np.array([height, 2*radius, 2*radius])

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

    transportBC_merged = []
    transportIC_merged = []

    govNodesTrspt = []
    govNodesTrsptBC = []
    rigidPlatesTrspt = []

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


    return node_coords, mechBC_merged,  vor, volumes, functions,  govNodes, govNodesMechBC, rigidPlates, transportBC_merged, govNodesTrspt, rigidPlatesTrspt, govNodesTrsptBC

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

    if SHOW_PLOT:
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

    modelVertices = utilitiesGeom.returnSelectedPtsRadial (radius-thickness-1e-3 , radius+1e-3 , vor.vertices)
    ### selecting vertices on the outer surface

    outerFaceBC = np.array([2,-1])
    outerFace = utilitiesGeom.returnSelectedPtsRadial (radius-minDist/2 , radius+minDist/2 , vor.vertices)
    for i in range (len(outerFace)):
        trsBC = utilitiesMech.transportBC(outerFace[i], outerFaceBC)
        transportBC_merged.append(trsBC)

    if SHOW_PLOT:
        fig = plt.figure()
        ax = Axes3D(fig)
        ax.scatter(vor.vertices[modelVertices,0], vor.vertices[modelVertices,1], vor.vertices[modelVertices,2])
        ax.scatter(vor.vertices[outerFace,0], vor.vertices[outerFace,1], vor.vertices[outerFace,2])
        plt.show()

    innerFaceBC = np.array([-1,-1])
    innerFace = utilitiesGeom.returnSelectedPtsRadial ((radius-thickness)-minDist/2 , (radius-thickness)+minDist/2, vor.vertices)
    for i in range (len(innerFace)):
        trsBC = utilitiesMech.transportBC(innerFace[i], innerFaceBC)
        transportBC_merged.append(trsBC)

    if SHOW_PLOT:
        fig = plt.figure()
        ax = Axes3D(fig)
        ax.scatter(vor.vertices[modelVertices,0], vor.vertices[modelVertices,1], vor.vertices[modelVertices,2])
        ax.scatter(vor.vertices[innerFace,0], vor.vertices[innerFace,1], vor.vertices[innerFace,2])
        plt.show()



    radii = np.zeros((len(node_coords))) + minDist
    return node_coords, mechBC_merged, govNodes, govNodesMechBC, rigidPlates, transportBC_merged, vor, volumes, functions, radii










#
######## METHOD FOR CREATING OF A 2D SUPPORTED CANTILEVER MODEL
def assemble2DCantileverBending (maxLim, minDist, trials):
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

    #adding mech boundary conditions
    for n in range ( nrOfPoints ):
        mBC = utilitiesMech.mechanicalBC(dim, oldLen + n, lineBC)
        mechBC_merged.append(mBC)
        #print('adding')

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
                             notchWidth = -1, node_coords_init=None, coupled=False):
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


    #an indent due to mirroring of the data for voronoi tess.
    notches=[]
    indent = 1e-8
    # the following is for remesher that works so far just for the bema WITHOUT notch
    if node_coords_init is None:
        #notchWidth = 1.5e-3 /2
        if notchWidth == -1:
            notchWidth = minDist/2
        else:
            notchWidth /= 2

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
    supportWidth = maxLim[0] / 20


    ##################### CONSTRAINTS AND RIGID PLATES
    #rigid plate left support
    indentRP = 1e-3
    leftRigidPlateMechBC = np.array([0,0,-1, -1,-1,-1])
    leftRigidPlate = utilitiesMech.RigidPlate(-1, 2, np.array([ -indentRP, supportWidth+indentRP, -indentRP, 2*indentRP ]))
    rigidPlates.append(leftRigidPlate)
    govNodes.append(np.array([ indent+supportWidth/2, indent ]))
    govNodesMechBC.append(utilitiesMech.mechanicalBC(dim, -1, leftRigidPlateMechBC))
    #rigid plate left support
    rightRigidPlateMechBC = np.array([-1,0,-1, -1,-1,-1])
    rightRigidPlate = utilitiesMech.RigidPlate(-2, 2, np.array([ maxLim[0] - 2*indentRP-supportWidth, maxLim[0] + indentRP, -indentRP, 2*indentRP ]))
    rigidPlates.append(rightRigidPlate)
    govNodes.append(np.array([ maxLim[0]-indent-supportWidth/2, indent ]))
    govNodesMechBC.append(utilitiesMech.mechanicalBC(dim, -2, rightRigidPlateMechBC))
    #rigid plate top load
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
    ####################


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
        nodeA =  np.array([indent + 0.5*maxLim[0]*(1-loadWidth) , maxLim[1] - indent])
        nodeB =  np.array([maxLim[0] - indent - 0.5*maxLim[0]*(1-loadWidth) , maxLim[1] - indent])
        pointGenerators.generateNodesLine2dRand(nodeA, nodeB, minDist/2, dim, node_coords,  trials, True, True)
        #######top of frac zone
        nodeA =  np.array([0.5*maxLim[0] - 0.5*maxLim[1]*(1-notch) , maxLim[1] - indent])
        nodeB =  np.array([0.5*maxLim[0] + 0.5*maxLim[1]*(1-notch), maxLim[1] - indent])
        pointGenerators.generateNodesLine2dRand(nodeA, nodeB, minDist, dim, node_coords,  trials, False, False)


        ########################################## rest of  faces
        nodeA =  np.array([indent + 0.5*maxLim[0]*(1-loadWidth) , maxLim[1] - indent])
        nodeB =  np.array([indent, maxLim[1] - indent])
        pointGenerators.generateNodesLine2dRand(nodeA, nodeB, minDist*2, dim, node_coords,  trials, False, False)
        nodeA =  np.array([maxLim[0] - indent  , maxLim[1] - indent])
        nodeB =  np.array([maxLim[0] - indent - 0.5*maxLim[0]*(1-loadWidth) , maxLim[1] - indent])
        pointGenerators.generateNodesLine2dRand(nodeA, nodeB, minDist*2, dim, node_coords,  trials, False, False)
        nodeA =  np.array([indent  ,  indent])
        nodeB =  np.array([maxLim[0]/2-notchWidth ,  indent])
        pointGenerators.generateNodesLine2dRand(nodeA, nodeB, minDist*2, dim, node_coords,  trials, False, True)
        nodeA =  np.array([maxLim[0]-indent  ,  indent])
        nodeB =  np.array([maxLim[0]/2+notchWidth ,  indent])
        pointGenerators.generateNodesLine2dRand(nodeA, nodeB, minDist*2, dim, node_coords,  trials, False, True)
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
        #"""
        """
        maxLimF = np.array([
        maxLim[0]*0.5 - maxLim[0]*(1-notch)/2,
        maxLim[1] - indent,
        maxLim[0]*0.5 + maxLim[0]*(1-notch)/2,
        maxLim[1]*notch])
        """

        if not orthogonalFracZone:
            pointGenerators.generateNodesRect(maxLimF, minDist, dim, trials, node_coords, useLowBound=True)
        else:
            maxLimF = np.array([
            maxLim[0],
            maxLim[1] - indent,
            maxLim[0]*0.5 + maxLim[1]*(1-notch)/2,
            maxLim[1]*notch])

            pointGenerators.generateOrtogrid(maxLimF, minDist, dim, node_coords, maxLim[1]*(1-notch))
            #fracZoneWidth*maxLim[0])


        ## notch faces

        maxLimF = np.array([
        maxLim[0]/2 - 0.5*maxLim[1]*(1-notch)*1.5,
        maxLim[1],
        maxLim[0]/2 + 0.5*maxLim[1]*(1-notch)*1.5,
        indent+maxLim[1]*notch/2])
        pointGenerators.generateNodesRect(maxLimF, minDist, dim, trials, node_coords, useLowBound=True)


        """
        maxLimF = np.array([
        maxLim[0]/2 - 0.5*maxLim[1]*(1-notch)*2.5,
        maxLim[1],
        maxLim[0]/2 + 0.5*maxLim[1]*(1-notch)*2.5,
        indent+maxLim[1]*notch/2])
        pointGenerators.generateNodesRect(maxLimF, minDist*1.5, dim, trials, node_coords, useLowBound=True)
        """

        """
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
        pointGenerators.generateNodesRect(maxLim, minDist*2, dim, trials, node_coords)


    """
    node_coords = np.asarray(node_coords)
    plt.plot(node_coords[:,0], node_coords[:,1], 'o', color='black');
    plt.show()
    """

    return node_coords, mechBC_merged, mechInitC_merged, notches, govNodes, govNodesMechBC, rigidPlates




def assemble2dCorrosionRebar(maxLim, minDist, trials, rebarMinDist, rebarDiameter, rebarCount, rebarDepth, sampleBorders, node_coords_init=None):
    dim = 2

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

    #an indent due to mirroring of the data for voronoi tess.
    indent = 1e-8
    indentRP = 1e-7

    if node_coords_init is None:

        node_coords.append(np.array([indent, indent]))

        if sampleBorders:
            """
            #top
            nodeA = np.array([indent, maxLim[1]-indent])
            nodeB = np.array([maxLim[0]-indent, maxLim[1]-indent])
            pointGenerators.generateNodesLine2dRand(nodeA, nodeB, minDist, dim, node_coords, trials, catchCorners=False, equidist=False)

            topRigidPlateMechBC = np.array([-1, -1,-1,   -1,-1,-1])
            topRigidPlate = utilitiesMech.RigidPlate(-1, 2, np.array([-indentRP, maxLim[0]+indentRP, maxLim[1]-indentRP, maxLim[1]+indentRP ]))
            rigidPlates.append(topRigidPlate)
            govNodes.append(np.array([ maxLim[0]/2, maxLim[1]-indent ]))
            govNodesMechBC.append(utilitiesMech.mechanicalBC(dim, -1, topRigidPlateMechBC))
            """

            #bottom
            nodeA = np.array([indent, indent])
            nodeB = np.array([maxLim[0]-indent, indent])
            pointGenerators.generateNodesLine2dRand(nodeA, nodeB, minDist*4, dim, node_coords, trials, catchCorners=False, equidist=False)

            bottomRigidPlateMechBC = np.array([0, 0,0,   -1,-1,-1])
            bottomRigidPlate = utilitiesMech.RigidPlate(-1, 2, np.array([-indentRP, maxLim[0]+indentRP, -indentRP, indentRP ]))
            rigidPlates.append(bottomRigidPlate)
            govNodes.append(np.array([ maxLim[0]/2, indentRP ]))
            govNodesMechBC.append(utilitiesMech.mechanicalBC(dim, -2, bottomRigidPlateMechBC))

            """
            #left
            nodeA = np.array([indent, indent])
            nodeB = np.array([indent, maxLim[1]-indent])
            pointGenerators.generateNodesLine2dRand(nodeA, nodeB, minDist*3, dim, node_coords, trials, catchCorners=False, equidist=False)

            leftRigidPlateMechBC = np.array([0, -1,-1,   -1,-1,-1])
            leftRigidPlate = utilitiesMech.RigidPlate(-1, 2, np.array([-indentRP, indentRP, minDist/2, maxLim[1]-minDist/2 ]))
            rigidPlates.append(leftRigidPlate)
            govNodes.append(np.array([ indentRP, maxLim[1]/2 ]))
            govNodesMechBC.append(utilitiesMech.mechanicalBC(dim, -3, leftRigidPlateMechBC))


            #right
            nodeA = np.array([maxLim[0]-indent, indent])
            nodeB = np.array([maxLim[0]-indent, maxLim[1]-indent])
            pointGenerators.generateNodesLine2dRand(nodeA, nodeB, minDist*3, dim, node_coords, trials, catchCorners=False, equidist=False)

            rightRigidPlateMechBC = np.array([0, -1,-1,   -1,-1,-1])
            rightRigidPlate = utilitiesMech.RigidPlate(-1, 2, np.array([maxLim[0]-indentRP, maxLim[0]+indentRP, minDist/2, maxLim[1]-minDist/2 ]))
            rigidPlates.append(rightRigidPlate)
            govNodes.append(np.array([ maxLim[0], maxLim[1]/2 ]))
            govNodesMechBC.append(utilitiesMech.mechanicalBC(dim, -4, rightRigidPlateMechBC))
            """


        #rebar
        for r in range (rebarCount):
            #rebar edge
            centre = np.array([ (maxLim[0]/rebarCount)*(r+0.5), maxLim[1]-rebarDepth  ])
            pointGenerators.generateNodesCircle2dRand(centre, rebarDiameter/2, rebarMinDist, node_coords, trials )
            #rebar crossection
            pointGenerators.generateNodesOrtoCircle2dRand(centre, rebarDiameter/2, rebarMinDist, node_coords, trials)

            govNodes.append(np.array( np.copy(centre) ))


        fineRegDepth = 3*rebarDepth

        #fine top half rect
        fineTopBounds = np.array([     indent,          maxLim[1] -fineRegDepth,      maxLim[0],         maxLim[1]   ])
        pointGenerators.generateNodesRect(fineTopBounds, minDist, dim, trials, node_coords, useLowBound=True)

        #intermediate rect
        interHeight = maxLim[1]  / 3
        interBounds = np.array([     indent,      maxLim[1] -fineRegDepth - interHeight , maxLim[0],             maxLim[1] *0.8  ])
        pointGenerators.generateNodesRect(interBounds, minDist, dim, trials, node_coords, useLowBound=True, topMinDist = 5*minDist, bottomMinDist = minDist, gradienDirection=1)

        #bottom rough rect
        roughBottomBounds =  np.array([     indent,      indent,  maxLim[0],              maxLim[1] -fineRegDepth - interHeight ])
        pointGenerators.generateNodesRect(roughBottomBounds, minDist*5, dim, trials, node_coords, useLowBound=True)




    node_coords = np.asarray(node_coords)
    if SHOW_PLOT:
        plt.plot(node_coords[:,0], node_coords[:,1], 'o', color='black');
        plt.show()


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



def assemble2dDogBone(D, minDist, trials, excentricity = 50, symmetric=False, edgeMinDistCoef = 1.0, roughDogBone=0, roughEdgeDogbone=0, roughMinDistCoef=1, interLayerThickness=2 ):

    if roughDogBone >0 :
        sampleCircularBorders = False
    else:
        sampleCircularBorders = True

    dim = 2
    #lists for the model
    node_coords = []
    mechBC_merged = []
    mechInitC_merged = []
    govNodes = []
    govNodesMechBC = []
    rigidPlates = []

    indent = 1e-6

    oldLen = len(node_coords)

    if roughDogBone > 0:
        altMinDist = roughMinDistCoef * minDist
    else:
        altMinDist = minDist

    if symmetric == True:
        nodeA = np.array([0.2*D+2*indent,  3/4 * D - indent -minDist/2])
        nodeB = np.array([D*0.8-2*indent, 3/4 * D - indent -minDist/2])
        pointGenerators.generateNodesLine2dRand(nodeA, nodeB, minDist, dim, node_coords, trials, True, True)

        nodeA = np.array([0.2*D+2*indent,  3/4 * D - indent +minDist/2])
        nodeB = np.array([D*0.8-2*indent, 3/4 * D - indent +minDist/2])
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
    nodeB = np.array([indent+D, indent])
    pointGenerators.generateNodesLine2dRand(nodeA, nodeB, altMinDist*0.9, dim, node_coords, trials, True, True)

    #bottom line of dogbone
    nodeA = np.array([indent,  6/4 * D - indent])
    nodeB = np.array([D-indent, 6/4 * D - indent])
    pointGenerators.generateNodesLine2dRand(nodeA, nodeB, altMinDist*0.9, dim, node_coords, trials, True, True)

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
            angleLimitA =   -np.arcsin(  2* roughDogBone * minDist / radius)
            angleLimitB =   np.arcsin( 2* roughDogBone * minDist / radius)
        #if symmetric == True:
        #    angleLimitB =   0
        mirroredPointsA =  pointGenerators.generateNodesCircle2dRand(centreA, radius+indent, minDist*edgeMinDistCoef, node_coords, trials, angleLimitA=angleLimitA, angleLimitB=angleLimitB, mirrorIndent = indent*10 )
        angleLimitA +=   np.pi
        angleLimitB +=   np.pi
        #if symmetric == True:
        #    angleLimitA =    np.pi
        mirroredPointsB =  pointGenerators.generateNodesCircle2dRand(centreB, radius+indent, minDist*edgeMinDistCoef, node_coords, trials, angleLimitA=angleLimitA, angleLimitB=angleLimitB, mirrorIndent = indent*10)
    else:
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
        pointGenerators.generateNodesRect(maxLim, minDist, 2, trials, node_coords)
        nrOfPoints =  (len(node_coords)) - oldLen


    node_coords_all = np.copy ( node_coords )
    node_coords_dogbone = []
    node_indices_dogbone = []

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
    print('done.')
    node_coords_dogbone = np.asarray(node_coords_dogbone)






    if roughEdgeDogbone == 2 or roughEdgeDogbone == 3:
        mirrored_coords = []
        #mirroring rough edge dogbone circular borders
        dogboneRadius = 0.725*D
        leftCenter = np.array( [-0.525 * D, 3/4 * D] )
        rightCenter = np.array( [ 1.525 * D, 3/4 * D] )
        for node in node_coords_dogbone:
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

        mirrored_coords = np.asarray(mirrored_coords)
        """
        plt.plot(mirrored_coords[:,0], mirrored_coords[:,1], 'o', color='black');
        plt.plot(node_coords_dogbone[:,0], node_coords_dogbone[:,1], 'x', color='red');
        plt.show()
        #"""
        node_coords_dogbone = np.vstack((node_coords_dogbone, mirrored_coords))
        node_coords_all = np.copy(node_coords_dogbone)

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


    """
    plt.plot(node_coords_all[:,0], node_coords_all[:,1], 'o', color='green');
    plt.plot(node_coords_all[node_indices_dogbone,0], node_coords_all[node_indices_dogbone,1], 'x', color='red');
    plt.show()
    #"""

    node_count = len (node_coords_all)
    return node_coords_all, node_indices_dogbone, mechBC_merged, mechInitC_merged, node_count, govNodes, govNodesMechBC, rigidPlates




def assemble3dDogBone(D, minDist, trials, excentricity = 20):
    dim = 3
    thickness = 0.1
    #lists for the model
    node_coords = []
    mechBC_merged = []
    mechInitC_merged = []
    node_coords_out = []
    govNodes = []
    govNodesMechBC = []
    rigidPlates = []

    indent = 1e-5

    node_coords.append(np.array([(D-indent)/2, indent, thickness/2]))
    node_coords.append(np.array([(D-indent)/2, 6/4 * D - indent, thickness/2]))


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
    govNodes.append(np.array([ D/2+D/excentricity, indent, thickness/2 ]))
    govNodesMechBC.append(utilitiesMech.mechanicalBC(dim, -1, topRigidPlateMechBC))
    #bottom rigid plate
    bottomRigidPlateMechBC =  np.array([0,0,0,  -1,-1,-1, -1,-1,-1, -1,-1,-1,])
    bottomRigidPlate = utilitiesMech.RigidPlate(-1, 3, np.array([
    -indentRP,
    indentRP+D,
    -indentRP+6/4 * D,
     +indentRP+6/4 * D,
     -indentRP,
     thickness+indentRP  ]))
    rigidPlates.append(bottomRigidPlate)
    govNodes.append(np.array([ D/2+D/excentricity, 6/4 * D-indent, thickness/2 ]))
    govNodesMechBC.append(utilitiesMech.mechanicalBC(dim, -1, bottomRigidPlateMechBC))
    #####################

    #top surface of dogbone
    nodeA = np.array([indent, indent, indent])
    nodeB = np.array([D-indent, indent, thickness-indent ])
    pointGenerators.generateNodesOrtoSurface3dRand(nodeA, nodeB, minDist, dim, node_coords, trials)

    #bottom surface of dogbone
    nodeA = np.array([indent,  6/4 * D - indent, indent])
    nodeB = np.array([D-indent, 6/4 * D - indent, thickness-indent])
    pointGenerators.generateNodesOrtoSurface3dRand(nodeA, nodeB, minDist, dim, node_coords, trials)


    #sampling on circular borders
    centreA = np.array( [-0.525 * D, 3/4 * D, indent] )
    centreB = np.array( [ 1.525 * D, 3/4 * D, indent] )
    radius = 0.725*D
    angleLimit =   np.arcsin( 0.5*D / radius)
    mirroredPointsA = pointGenerators.generateNodesOrtoCilinderSurf3dRand(centreA, radius+indent, thickness, 2, minDist*0.9, node_coords, trials, angleLimitA=-angleLimit, angleLimitB=angleLimit, mirrorIndent = indent )

    mirroredPointsB = pointGenerators.generateNodesOrtoCilinderSurf3dRand(centreB, radius+indent, thickness, 2, minDist*0.9, node_coords, trials, angleLimitA=np.pi-angleLimit, angleLimitB=np.pi+angleLimit, mirrorIndent = indent)


    #rectangle of dogbone
    # xrectBC = np.array([-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1])
    oldLen = len(node_coords)
    maxLim = np.array([  D    ,  6/4*D, thickness ])
    pointGenerators.generateNodesRect(maxLim, minDist, 3, trials, node_coords)
    nrOfPoints =  (len(node_coords)) - oldLen

    #dumping points outside bone

    radius = np.linalg.norm( centreB[0:2] - np.array([D, 1/4*D]))
    print('Dumping points within bordering circles...')
    node_coords_out = []
    for node in node_coords:
        distA = np.linalg.norm( node[0:2] - centreA[0:2])
        distB = np.linalg.norm( node[0:2] - centreB[0:2])
        if (distA > radius and distB > radius):
            node_coords_out.append(node)

    print('done.')
    #mirrored circles. Not to be in the model at the end
    node_count = len(node_coords_out)
    node_coords_out = np.vstack( (node_coords_out, mirroredPointsA, mirroredPointsB) )


    #snode_coords_out = node_coords
    node_count = len(node_coords_out)

    return node_coords_out, mechBC_merged, mechInitC_merged, node_count,govNodes, govNodesMechBC, rigidPlates





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



def assemble3DSSBeamBending (maxLim, minDist, trials, notch, loadWidth,  fracZoneWidth = 0.15, orthogonalFracZone=False, notchWidth = -1, coupled=False, node_coords_init=None):
    minDist *=2
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

    #exact notch
    exactNotch = False

    indent = 1e-7
    notches=[]

    #notch heaight
    if not exactNotch:
        nHeight = maxLim[1]*notch - minDist/4
        #notch = nHeight / maxLim[1]


    if node_coords_init is None:
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
            oldLen = len(node_coords)
            nodeA = np.array([maxLim[0]/2-notchWidth, indent, indent])
            nodeB = np.array([maxLim[0]/2-notchWidth, indent, maxLim[2]-indent])
            pointGenerators.generateNodesLine3dRand(nodeA, nodeB, minDist/4, dim, node_coords, trials, catchCorners=True, equidist=True)
            """
            nodeA = np.array([maxLim[0]/2-notchWidth, maxLim[1]*notch-minDist/8, indent])
            nodeB = np.array([maxLim[0]/2-notchWidth, maxLim[1]*notch-minDist/8, maxLim[2]-indent])
            pointGenerators.generateNodesLine3dRand(nodeA, nodeB, minDist/4, dim, node_coords, trials, catchCorners=True, equidist=True)
            """
            nodeA = np.array([maxLim[0]/2-notchWidth, indent, indent])
            nodeB = np.array([maxLim[0]/2-notchWidth, maxLim[1]*notch-minDist/2, indent])
            pointGenerators.generateNodesLine3dRand(nodeA, nodeB, minDist/4, dim, node_coords, trials, catchCorners=False, equidist=True)
            nodeA = np.array([maxLim[0]/2-notchWidth, indent, maxLim[2]-indent])
            nodeB = np.array([maxLim[0]/2-notchWidth, maxLim[1]*notch-minDist/2, maxLim[2]-indent])
            pointGenerators.generateNodesLine3dRand(nodeA, nodeB, minDist/4, dim, node_coords, trials, catchCorners=False, equidist=True)

            nodeA = np.array([maxLim[0]/2-notchWidth, indent, indent])
            nodeB = np.array([maxLim[0]/2-notchWidth, maxLim[1]*notch-minDist/2, maxLim[2]-indent])
            pointGenerators.generateNodesOrtoSurface3dRand(nodeA, nodeB, minDist/3, dim, node_coords, 50000,minDistAmongNewPoints=True)


            for i in range (oldLen, len(node_coords), 1):
                notchSide0.append(i)


            notchSide1 = []
            oldLen = len(node_coords)
            nodeA = np.array([maxLim[0]/2+notchWidth, indent, indent])
            nodeB = np.array([maxLim[0]/2+notchWidth, indent, maxLim[2]-indent])
            pointGenerators.generateNodesLine3dRand(nodeA, nodeB, minDist/4, dim, node_coords, trials, catchCorners=True, equidist=True)
            """
            nodeA = np.array([maxLim[0]/2+notchWidth, maxLim[1]*notch-minDist/8, indent])
            nodeB = np.array([maxLim[0]/2+notchWidth, maxLim[1]*notch-minDist/8, maxLim[2]-indent])
            pointGenerators.generateNodesLine3dRand(nodeA, nodeB, minDist/4, dim, node_coords, trials, catchCorners=True, equidist=True)
            """
            nodeA = np.array([maxLim[0]/2+notchWidth, indent, indent])
            nodeB = np.array([maxLim[0]/2+notchWidth, maxLim[1]*notch-minDist/2, indent])
            pointGenerators.generateNodesLine3dRand(nodeA, nodeB, minDist/4, dim, node_coords, trials, catchCorners=False, equidist=True)
            nodeA = np.array([maxLim[0]/2+notchWidth, indent, maxLim[2]-indent])
            nodeB = np.array([maxLim[0]/2+notchWidth, maxLim[1]*notch-minDist/2, maxLim[2]-indent])
            pointGenerators.generateNodesLine3dRand(nodeA, nodeB, minDist/4, dim, node_coords, trials, catchCorners=False, equidist=True)

            nodeA = np.array([maxLim[0]/2+notchWidth, indent, indent])
            nodeB = np.array([maxLim[0]/2+notchWidth, maxLim[1]*notch-minDist/2, maxLim[2]-indent])
            pointGenerators.generateNodesOrtoSurface3dRand(nodeA, nodeB, minDist/3, dim, node_coords, 50000, minDistAmongNewPoints= True)



            for i in range (oldLen, len(node_coords), 1):
                notchSide1.append(i)

            notchL = []
            notchL.append(notchSide0)
            notchL.append(notchSide1)
            notches.append(notchL)

            leftTop = []
            oldLen = len(node_coords)
            if orthogonalFracZone:
                nodeA = np.array([maxLim[0]/2-notchWidth, maxLim[1]*notch-minDist/4, indent])
                nodeB = np.array([maxLim[0]/2-notchWidth, maxLim[1]*notch-minDist/4, maxLim[2]-indent])
            else:
                nodeA = np.array([maxLim[0]/2-notchWidth, maxLim[1]*notch, indent])
                nodeB = np.array([maxLim[0]/2-notchWidth, maxLim[1]*notch, maxLim[2]-indent])
            pointGenerators.generateNodesLine3dRand(nodeA, nodeB, minDist/4, dim, node_coords, trials, catchCorners=False, equidist=True)
            for i in range (oldLen, len(node_coords), 1):
                if exactNotch:
                    leftTop.append(i)
                else:
                    notchSide0.append(i)

            rightTop = []
            oldLen = len(node_coords)
            if orthogonalFracZone:
                nodeA = np.array([maxLim[0]/2+notchWidth, maxLim[1]*notch-minDist/4, indent])
                nodeB = np.array([maxLim[0]/2+notchWidth, maxLim[1]*notch-minDist/4, maxLim[2]-indent])
            else:
                nodeA = np.array([maxLim[0]/2+notchWidth, maxLim[1]*notch, indent])
                nodeB = np.array([maxLim[0]/2+notchWidth, maxLim[1]*notch, maxLim[2]-indent])
            pointGenerators.generateNodesLine3dRand(nodeA, nodeB, minDist/4, dim, node_coords, trials, catchCorners=False, equidist=True)
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

    """
    node_coords = np.asarray(node_coords)
    if SHOW_PLOT:
        fig = plt.figure()
        ax = Axes3D(fig)
        ax.scatter(node_coords[:,0], node_coords[:,1], node_coords[:,2])
        plt.show()
    """

    #width of the supports
    supportWidth = maxLim[0] / 20

    ##################### CONSTRAINTS AND RIGID PLATES
    #rigid plate left support
    indentRP = 1e-3
    leftRigidPlateMechBC = np.array([0,0,0, 0,0,-1,  -1,-1,-1,-1,-1,-1])
    leftRigidPlate = utilitiesMech.RigidPlate(-1, 3, np.array([ -indentRP, supportWidth+indentRP, -indentRP, 2*indentRP, -indentRP, maxLim[2]+indentRP  ]))
    rigidPlates.append(leftRigidPlate)
    govNodes.append(np.array([ indent+supportWidth/2, indent, maxLim[2]/2 ]))
    govNodesMechBC.append(utilitiesMech.mechanicalBC(dim, -1, leftRigidPlateMechBC))
    #rigid plate left support
    rightRigidPlateMechBC = np.array([-1,0,0, 0,0,-1,  -1,-1,-1,-1,-1,-1])
    rightRigidPlate = utilitiesMech.RigidPlate(-2, 3, np.array([ maxLim[0] - 2*indentRP-supportWidth, maxLim[0] + indentRP, -indentRP, 2*indentRP, -indentRP, maxLim[2]+indentRP  ]))
    rigidPlates.append(rightRigidPlate)
    govNodes.append(np.array([ maxLim[0]-indent-supportWidth/2, indent, maxLim[2]/2 ]))
    govNodesMechBC.append(utilitiesMech.mechanicalBC(dim, -2, rightRigidPlateMechBC))
    #rigid plate top load
    topRigidPlateMechBC = np.array([-1,1,-1, -1,-1,-1,  -1,-1,-1,-1,-1,-1])
    topRigidPlate = utilitiesMech.RigidPlate(-3, 3,
    np.array([
     0.5*maxLim[0]*(1-loadWidth)-indentRP,
    maxLim[0]  - 0.5*maxLim[0]*(1-loadWidth)+indentRP,
    maxLim[1] - 2*indentRP,
    maxLim[1] + 2*indentRP,
    -indentRP,
    maxLim[2] + indentRP  ]))
    rigidPlates.append(topRigidPlate)
    govNodes.append(np.array([ maxLim[0]/2, maxLim[1]-indent, maxLim[2]/2 ]))
    govNodesMechBC.append(utilitiesMech.mechanicalBC(dim, -3, topRigidPlateMechBC))
    ####################

    if node_coords_init is None:

        ###############generating of nodes, left horizontal support ###############
        nodeA = np.array([indent, indent, indent])
        nodeB = np.array([indent, indent, maxLim[2]-indent])
        pointGenerators.generateNodesLine3dRand(nodeA, nodeB, minDist/3, dim, node_coords, trials, True, True)
        nodeA = np.array([indent+supportWidth, indent, indent])
        nodeB = np.array([indent+supportWidth, indent, maxLim[2]-indent])
        pointGenerators.generateNodesLine3dRand(nodeA, nodeB, minDist/3, dim, node_coords, trials, True, True)


        ###############generating of nodes, right horizontal support ###############
        nodeA = np.array([maxLim[0] - indent, indent, indent])
        nodeB = np.array([maxLim[0] - indent, indent, maxLim[2]-indent])
        pointGenerators.generateNodesLine3dRand(nodeA, nodeB, minDist/3, dim, node_coords, trials, True, True)
        nodeA = np.array([maxLim[0] - indent-supportWidth, indent, indent])
        nodeB = np.array([maxLim[0] - indent-supportWidth, indent, maxLim[2]-indent])
        pointGenerators.generateNodesLine3dRand(nodeA, nodeB, minDist/3, dim, node_coords, trials, True, True)


        oldLen = len(node_coords)
        ##########################################generating of points, fracture zone
        maxLimF = np.array([
        maxLim[0]  - 0.5*maxLim[0]*(1-fracZoneWidth),
        maxLim[1] - indent,
        maxLim[2] - indent,
        0.5*maxLim[0]*(1-fracZoneWidth),
        maxLim[1]*notch-indent,
        indent])
        if not orthogonalFracZone:
            pointGenerators.generateNodesRect(maxLimF, minDist/2, dim, trials, node_coords, useLowBound=True)
        else:
            #pointGenerators.generateOrtogrid(maxLimF, minDist/2, dim, node_coords, maxLim[0]*fracZoneWidth)  #
            pointGenerators.generateOrtogrid_variable(maxLimF, minDist/2, node_coords, np.array([maxLim[0]*fracZoneWidth, maxLim[1]*(1-notch)-indent*2, maxLim[2]-indent*2]) )
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

        ##########################################generating of points, fracture zone
        maxLimF = np.array([
        maxLim[0] - indent - 0.5*maxLim[0]*(1-fracZoneWidth*2),
        maxLim[1] -indent,#* notch,
        maxLim[2] - indent,
        indent + 0.5*maxLim[0]*(1-fracZoneWidth*2),
        indent,
        indent])
        pointGenerators.generateNodesRect(maxLimF, minDist/2, dim, 50000, node_coords, useLowBound=True)


        ###############generating of nodes, front bottom line ###############
        nodeA = np.array([indent, indent, indent])
        nodeB = np.array([maxLim[0]-indent, indent, indent])
        pointGenerators.generateNodesLine3dRand(nodeA, nodeB, minDist, dim, node_coords, trials, False, False)
        ###############generating of nodes, rear bottom line ###############
        nodeA = np.array([indent, indent,  maxLim[2]-indent])
        nodeB = np.array([maxLim[0]-indent, indent,  maxLim[2]-indent])
        pointGenerators.generateNodesLine3dRand(nodeA, nodeB, minDist, dim, node_coords, trials, False, False)
        ###############generating of nodes, front top line ###############
        nodeA = np.array([indent, maxLim[1]-indent, indent])
        nodeB = np.array([maxLim[0]-indent, maxLim[1]-indent, indent])
        pointGenerators.generateNodesLine3dRand(nodeA, nodeB, minDist, dim, node_coords, trials, False, False)
        ###############generating of nodes, rear top line ###############
        nodeA = np.array([indent, maxLim[1]-indent,  maxLim[2]-indent])
        nodeB = np.array([maxLim[0]-indent, maxLim[1]-indent,  maxLim[2]-indent])
        pointGenerators.generateNodesLine3dRand(nodeA, nodeB, minDist, dim, node_coords, trials, False, False)

        ############### loaded top face ###############
        lineBC = np.array([-1,-1,-1,-1,-1,-1,  -1, 1,-1,-1,-1,-1])
        nodeA =  np.array([indent + 0.5*maxLim[0]*(1-loadWidth), maxLim[1] - indent, indent])
        nodeB =  np.array([maxLim[0] - indent - 0.5*maxLim[0]*(1-loadWidth), maxLim[1] - indent, maxLim[2] - indent])
        pointGenerators.generateNodesOrtoSurface3dRand(nodeA, nodeB, minDist/2, dim, node_coords, 10000)

        #front surf
        nodeA =  np.array([indent , maxLim[1] - indent, indent])
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


        """
        node_coords = np.asarray(node_coords)
        if SHOW_PLOT:
            fig = plt.figure()
            ax = Axes3D(fig)
            ax.scatter(node_coords[:,0],node_coords[:,1],node_coords[:,2])
            plt.show()
        """


        ##########################################generating of points, left support
        maxLimF = np.array([
        supportWidth*2,
        supportWidth*2,
        maxLim[2]-indent,
        indent,
        indent,
        indent])
        pointGenerators.generateNodesRect(maxLimF, minDist, dim, trials, node_coords, useLowBound=True)
        ##########################################generating of points, right support
        maxLimF = np.array([
        maxLim[0],
        supportWidth*2,
        maxLim[2]-indent,
        maxLim[0]-supportWidth*2,
        indent,
        indent])
        pointGenerators.generateNodesRect(maxLimF, minDist, dim, trials, node_coords, useLowBound=True)



        ##########################################generating of points, homogeneous volume
        pointGenerators.generateNodesRect(maxLim, minDist, dim, trials, node_coords)

    #if coupled:
    #    notches = []
    return node_coords, mechBC_merged, mechInitC_merged, notches, govNodes, govNodesMechBC, rigidPlates





















def assemble3Dcube(maxLim, minDist, trials, powerTes, coupled=False, node_coords_init=None):
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
            """
            mechBC = np.array([-1,0,0,   -1,-1,-1,    -1,-1,-1,-1,-1,-1])
            mBC = utilitiesMech.mechanicalBC(dim, 0, mechBC)
            mechBC_merged.append(mBC)
            """

            ###############generating of nodes, front bottom line ###############
            nodeA = np.array([indent, indent, indent])
            nodeB = np.array([maxLim[0]-indent, indent, indent])
            pointGenerators.generateNodesLine3dRand(nodeA, nodeB, minDist, dim, node_coords, trials, False, False)
            ###############generating of nodes, rear bottom line ###############
            nodeA = np.array([indent, indent,  maxLim[2]-indent])
            nodeB = np.array([maxLim[0]-indent, indent,  maxLim[2]-indent])
            pointGenerators.generateNodesLine3dRand(nodeA, nodeB, minDist, dim, node_coords, trials, False, False)
            ###############generating of nodes, front top line ###############
            nodeA = np.array([indent, maxLim[1]-indent, indent])
            nodeB = np.array([maxLim[0]-indent, maxLim[1]-indent, indent])
            pointGenerators.generateNodesLine3dRand(nodeA, nodeB, minDist, dim, node_coords, trials, False, False)
            ###############generating of nodes, rear top line ###############
            nodeA = np.array([indent, maxLim[1]-indent,  maxLim[2]-indent])
            nodeB = np.array([maxLim[0]-indent, maxLim[1]-indent,  maxLim[2]-indent])
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
            node_coords, radii = pointGenerators.generateParticlesRect(maxLim, minDist*0.4, minDist, 0.8, dim, trials, node_coords, radii, allow_domain_overlap = True, periodic_distance=True)



        node_coords = np.asarray(node_coords)
        """
        fig = plt.figure()
        ax = Axes3D(fig)
        ax.scatter(node_coords[:,0],node_coords[:,1],node_coords[:,2])
        plt.show()
        """


    #if coupled:
    #    notches = []
    return node_coords, mechBC_merged, mechInitC_merged, govNodes, govNodesMechBC, rigidPlates




def assemble3DReinhardtTension (maxLim, minDist, trials, fracZoneWidth = 0.15):
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
        pointGenerators.generateNodesLine3dRand(nodeA, nodeB, minDist/3, dim, node_coords, trials, catchCorners=True, equidist=True)
        nodeA = np.array([maxLim[0]/2-notchWidth, notchDepth, indent])
        nodeB = np.array([maxLim[0]/2-notchWidth, notchDepth, maxLim[2]-indent])
        pointGenerators.generateNodesLine3dRand(nodeA, nodeB, minDist/3, dim, node_coords, trials, catchCorners=True, equidist=True)
        nodeA = np.array([maxLim[0]/2-notchWidth, indent, indent])
        nodeB = np.array([maxLim[0]/2-notchWidth, notchDepth, indent])
        pointGenerators.generateNodesLine3dRand(nodeA, nodeB, minDist/3, dim, node_coords, trials, catchCorners=False, equidist=True)
        nodeA = np.array([maxLim[0]/2-notchWidth, indent, maxLim[2]-indent])
        nodeB = np.array([maxLim[0]/2-notchWidth, notchDepth, maxLim[2]-indent])
        pointGenerators.generateNodesLine3dRand(nodeA, nodeB, minDist/3, dim, node_coords, trials, catchCorners=False, equidist=True)

        nodeA = np.array([maxLim[0]/2-notchWidth, indent, indent])
        nodeB = np.array([maxLim[0]/2-notchWidth, notchDepth, maxLim[2]-indent])
        pointGenerators.generateNodesOrtoSurface3dRand(nodeA, nodeB, minDist/3, dim, node_coords, trials,minDistAmongNewPoints=True)
        for i in range (oldLen, len(node_coords), 1):
            notchSide0.append(i)


        notchSide1 = []
        oldLen = len(node_coords)
        nodeA = np.array([maxLim[0]/2+notchWidth, indent, indent])
        nodeB = np.array([maxLim[0]/2+notchWidth, indent, maxLim[2]-indent])
        pointGenerators.generateNodesLine3dRand(nodeA, nodeB, minDist/3, dim, node_coords, trials, catchCorners=True, equidist=True)
        nodeA = np.array([maxLim[0]/2+notchWidth, notchDepth, indent])
        nodeB = np.array([maxLim[0]/2+notchWidth, notchDepth, maxLim[2]-indent])
        pointGenerators.generateNodesLine3dRand(nodeA, nodeB, minDist/3, dim, node_coords, trials, catchCorners=True, equidist=True)
        nodeA = np.array([maxLim[0]/2+notchWidth, indent, indent])
        nodeB = np.array([maxLim[0]/2+notchWidth, notchDepth, indent])
        pointGenerators.generateNodesLine3dRand(nodeA, nodeB, minDist/3, dim, node_coords, trials, catchCorners=False, equidist=True)
        nodeA = np.array([maxLim[0]/2+notchWidth, indent, maxLim[2]-indent])
        nodeB = np.array([maxLim[0]/2+notchWidth, notchDepth, maxLim[2]-indent])
        pointGenerators.generateNodesLine3dRand(nodeA, nodeB, minDist/3, dim, node_coords, trials, catchCorners=False, equidist=True)

        nodeA = np.array([maxLim[0]/2+notchWidth, indent, indent])
        nodeB = np.array([maxLim[0]/2+notchWidth, notchDepth, maxLim[2]-indent])
        pointGenerators.generateNodesOrtoSurface3dRand(nodeA, nodeB, minDist/3, dim, node_coords, trials, minDistAmongNewPoints= True)
        for i in range (oldLen, len(node_coords), 1):
            notchSide1.append(i)





        notchSide2 = []
        oldLen = len(node_coords)
        nodeA = np.array([maxLim[0]/2-notchWidth, maxLim[1]-indent, indent])
        nodeB = np.array([maxLim[0]/2-notchWidth, maxLim[1]-indent, maxLim[2]-indent])
        pointGenerators.generateNodesLine3dRand(nodeA, nodeB, minDist/3, dim, node_coords, trials, catchCorners=True, equidist=True)
        nodeA = np.array([maxLim[0]/2-notchWidth, maxLim[1]-notchDepth, indent])
        nodeB = np.array([maxLim[0]/2-notchWidth, maxLim[1]-notchDepth, maxLim[2]-indent])
        pointGenerators.generateNodesLine3dRand(nodeA, nodeB, minDist/3, dim, node_coords, trials, catchCorners=True, equidist=True)
        nodeA = np.array([maxLim[0]/2-notchWidth, maxLim[1]-indent, indent])
        nodeB = np.array([maxLim[0]/2-notchWidth, maxLim[1]-notchDepth, indent])
        pointGenerators.generateNodesLine3dRand(nodeA, nodeB, minDist/3, dim, node_coords, trials, catchCorners=False, equidist=True)
        nodeA = np.array([maxLim[0]/2-notchWidth, maxLim[1]-indent, maxLim[2]-indent])
        nodeB = np.array([maxLim[0]/2-notchWidth, maxLim[1]-notchDepth, maxLim[2]-indent])
        pointGenerators.generateNodesLine3dRand(nodeA, nodeB, minDist/3, dim, node_coords, trials, catchCorners=False, equidist=True)

        nodeA = np.array([maxLim[0]/2-notchWidth, maxLim[1]-indent, indent])
        nodeB = np.array([maxLim[0]/2-notchWidth, maxLim[1]-notchDepth, maxLim[2]-indent])
        pointGenerators.generateNodesOrtoSurface3dRand(nodeA, nodeB, minDist/4, dim, node_coords, trials,minDistAmongNewPoints=True)
        for i in range (oldLen, len(node_coords), 1):
            notchSide2.append(i)


        notchSide3 = []
        oldLen = len(node_coords)
        nodeA = np.array([maxLim[0]/2+notchWidth, maxLim[1]-indent, indent])
        nodeB = np.array([maxLim[0]/2+notchWidth, maxLim[1]-indent, maxLim[2]-indent])
        pointGenerators.generateNodesLine3dRand(nodeA, nodeB, minDist/3, dim, node_coords, trials, catchCorners=True, equidist=True)
        nodeA = np.array([maxLim[0]/2+notchWidth, maxLim[1]-notchDepth, indent])
        nodeB = np.array([maxLim[0]/2+notchWidth, maxLim[1]-notchDepth, maxLim[2]-indent])
        pointGenerators.generateNodesLine3dRand(nodeA, nodeB, minDist/3, dim, node_coords, trials, catchCorners=True, equidist=True)
        nodeA = np.array([maxLim[0]/2+notchWidth, maxLim[1]-indent, indent])
        nodeB = np.array([maxLim[0]/2+notchWidth, maxLim[1]-notchDepth, indent])
        pointGenerators.generateNodesLine3dRand(nodeA, nodeB, minDist/3, dim, node_coords, trials, catchCorners=False, equidist=True)
        nodeA = np.array([maxLim[0]/2+notchWidth, maxLim[1]-indent, maxLim[2]-indent])
        nodeB = np.array([maxLim[0]/2+notchWidth, maxLim[1]-notchDepth, maxLim[2]-indent])
        pointGenerators.generateNodesLine3dRand(nodeA, nodeB, minDist/3, dim, node_coords, trials, catchCorners=False, equidist=True)

        nodeA = np.array([maxLim[0]/2+notchWidth, maxLim[1]-indent, indent])
        nodeB = np.array([maxLim[0]/2+notchWidth, maxLim[1]-notchDepth, maxLim[2]-indent])
        pointGenerators.generateNodesOrtoSurface3dRand(nodeA, nodeB, minDist/3, dim, node_coords, trials, minDistAmongNewPoints= True)
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
    rightRigidPlateMechBC = np.array([-1,0,0, 0,0,0,  1,-1,-1,-1,-1,-1])
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




    ###############generating of nodes, front bottom line ###############
    nodeA = np.array([indent, indent, indent])
    nodeB = np.array([maxLim[0]-indent, indent, indent])
    pointGenerators.generateNodesLine3dRand(nodeA, nodeB, minDist/1.5, dim, node_coords, trials, False, False)
    ###############generating of nodes, rear bottom line ###############
    nodeA = np.array([indent, indent,  maxLim[2]-indent])
    nodeB = np.array([maxLim[0]-indent, indent,  maxLim[2]-indent])
    pointGenerators.generateNodesLine3dRand(nodeA, nodeB, minDist/1.5, dim, node_coords, trials, False, False)
    ###############generating of nodes, front top line ###############
    nodeA = np.array([indent, maxLim[1]-indent, indent])
    nodeB = np.array([maxLim[0]-indent, maxLim[1]-indent, indent])
    pointGenerators.generateNodesLine3dRand(nodeA, nodeB, minDist/1.5, dim, node_coords, trials, False, False)
    ###############generating of nodes, rear top line ###############
    nodeA = np.array([indent, maxLim[1]-indent,  maxLim[2]-indent])
    nodeB = np.array([maxLim[0]-indent, maxLim[1]-indent,  maxLim[2]-indent])
    pointGenerators.generateNodesLine3dRand(nodeA, nodeB, minDist/1.5, dim, node_coords, trials, False, False)



    #front surf
    nodeA =  np.array([indent , maxLim[1] - indent, indent])
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
    #edges
    nodeA = np.array([indent, indent, indent])
    nodeB = np.array([indent, maxLim[1]-indent, indent])
    pointGenerators.generateNodesLine3dRand(nodeA, nodeB, minDist/2, dim, node_coords, trials, True, False)
    nodeA = np.array([indent, indent, maxLim[2]-indent])
    nodeB = np.array([indent, maxLim[1]-indent, maxLim[2]-indent])
    pointGenerators.generateNodesLine3dRand(nodeA, nodeB, minDist/2, dim, node_coords, trials, True, False)
    nodeA = np.array([indent, indent, indent])
    nodeB = np.array([indent, indent, maxLim[2]-indent])
    pointGenerators.generateNodesLine3dRand(nodeA, nodeB, minDist/2, dim, node_coords, trials, True, False)
    nodeA = np.array([indent, maxLim[1]-indent, indent])
    nodeB = np.array([indent, maxLim[1]-indent, maxLim[2]-indent])
    pointGenerators.generateNodesLine3dRand(nodeA, nodeB, minDist/2, dim, node_coords, trials, True, False)
    #surf
    nodeA =  np.array([indent , indent, indent])
    nodeB =  np.array([indent, maxLim[1] - indent, maxLim[2] - indent])
    pointGenerators.generateNodesOrtoSurface3dRand(nodeA, nodeB, minDist/2, dim, node_coords, trials)

    #right face surf
    #edges
    nodeA = np.array([maxLim[0]-indent, indent, indent])
    nodeB = np.array([maxLim[0]-indent, maxLim[1]-indent, indent])
    pointGenerators.generateNodesLine3dRand(nodeA, nodeB, minDist/2, dim, node_coords, trials, True, False)
    nodeA = np.array([maxLim[0]-indent, indent, maxLim[2]-indent])
    nodeB = np.array([maxLim[0]-indent, maxLim[1]-indent, maxLim[2]-indent])
    pointGenerators.generateNodesLine3dRand(nodeA, nodeB, minDist/2, dim, node_coords, trials, True, False)
    nodeA = np.array([maxLim[0]-indent, indent, indent])
    nodeB = np.array([maxLim[0]-indent, indent, maxLim[2]-indent])
    pointGenerators.generateNodesLine3dRand(nodeA, nodeB, minDist/2, dim, node_coords, trials, True, False)
    nodeA = np.array([maxLim[0]-indent, maxLim[1]-indent, indent])
    nodeB = np.array([maxLim[0]-indent, maxLim[1]-indent, maxLim[2]-indent])
    pointGenerators.generateNodesLine3dRand(nodeA, nodeB, minDist/2, dim, node_coords, trials, True, False)
    #surf
    nodeA =  np.array([maxLim[0]-indent , indent, indent])
    nodeB =  np.array([maxLim[0]-indent, maxLim[1] - indent, maxLim[2] - indent])
    pointGenerators.generateNodesOrtoSurface3dRand(nodeA, nodeB, minDist/2, dim, node_coords, trials)

    #fracZoneWidth = 0.15
    """
    node_coords = np.asarray(node_coords)
    if SHOW_PLOT:
        fig = plt.figure()
        ax = Axes3D(fig)
        ax.scatter(node_coords[:,0], node_coords[:,1], node_coords[:,2])
        plt.show()
    """
    ##########################################generating of points, left support area
    maxLimF = np.array([
    indent + maxLim[0]*fracZoneWidth/3,
    maxLim[1]- indent,
    maxLim[2]- indent,
    indent,
    indent,
    indent])
    pointGenerators.generateNodesRect(maxLimF, minDist/1.5, dim, trials, node_coords, useLowBound=True)

    ##########################################generating of points, right support area
    maxLimF = np.array([
    maxLim[0]- indent,
    maxLim[1]- indent,
    maxLim[2]- indent,
    maxLim[0]*(1-fracZoneWidth/3),
    indent,
    indent])
    pointGenerators.generateNodesRect(maxLimF, minDist/1.5, dim, trials, node_coords, useLowBound=True)


    ##########################################generating of points, fracture zone
    maxLimF = np.array([
    maxLim[0] - indent - 0.5*maxLim[0]*(1-fracZoneWidth),
    maxLim[1] - indent,
    maxLim[2] - indent,
    indent + 0.5*maxLim[0]*(1-fracZoneWidth),
    maxLim[1]*notch/2*0+indent,
    indent])
    pointGenerators.generateNodesRect(maxLimF, minDist/2, dim, trials, node_coords, useLowBound=True)



    ##########################################generating of points, homogeneous volume
    pointGenerators.generateNodesRect(maxLim, minDist/1, dim, trials, node_coords)

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



def assemble3dcylinderTorsionFree(center, radius, height, minDist, trials, directionDim, functions):
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



def assemble3dcylinderTorsionPressFree(center, radius, height, minDist, trials, directionDim, functions):
    indent = 1e-5
    dim=3
    #lists for the model
    node_coords = []
    mechBC_merged = []
    mechInitC_merged = []
    govNodes = []
    govNodesMechBC = []
    rigidPlates = []

    ### fixed nodes
    mechBC = np.array([0,0,0,-1,-1,-1,    -1,-1,-1,-1,-1,-1])
    node_coords.append( np.array([indent*3,indent,indent]))
    mBC = utilitiesMech.mechanicalBC(dim, 0, mechBC)
    mechBC_merged.append(mBC)
    mechBC = np.array([-1,0,0,-1,-1,-1,    -1,-1,-1,-1,-1,-1])
    node_coords.append( np.array([height-indent*3,indent,indent]))
    mBC = utilitiesMech.mechanicalBC(dim, 1, mechBC)
    mechBC_merged.append(mBC)


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
    pointGenerators.generateNodesOrtoCircleBorder3dRand(center, radius-1e-5, directionDim, minDist, node_coords, trials)
    pointGenerators.generateNodesOrtoCircle3dRand(center, radius-1e-5, directionDim, minDist, node_coords, trials)

    print ('Nodes so far: %d' %len(node_coords))

    ###############generating of points loaded surface right face ###############
    nodeA = center.copy()
    nodeA[directionDim] += height-indent
    pointGenerators.generateNodesOrtoCircleBorder3dRand(nodeA, radius-1e-5,  directionDim, minDist, node_coords, trials)
    pointGenerators.generateNodesOrtoCircle3dRand(nodeA, radius-1e-5,  directionDim, minDist, node_coords, trials)

    print ('Nodes so far: %d' %len(node_coords))

    ###############generating of points cylinder surf###############
    pointGenerators.generateNodesOrtoCilinderSurf3dRand(center, radius-1e-5, height, directionDim, minDist,  node_coords, trials)

    print ('Nodes so far: %d' %len(node_coords))

    ###############generating of points cylinder volume ###############
    pointGenerators.generateNodesOrtoCilinder3dRand(center, radius-1e-5, height, directionDim, minDist,  node_coords, trials)
    #######################################################################

    print ('Nodes so far: %d' %len(node_coords))

    """
    node_coords = np.asarray(node_coords)
    if SHOW_PLOT:
        fig = plt.figure()
        ax = Axes3D(fig)
        ax.scatter(node_coords[:,0], node_coords[:,1], node_coords[:,2])
        plt.show()
    #"""




    return node_coords, mechBC_merged,  govNodes, govNodesMechBC, rigidPlates


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
    leftRigidPlateMechBC = np.array([0,-1,-1, 0,0,0,  -1,-1,-1,-1,-1,-1])
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






def assembleCoupledBrazilianDisc(center, radius, height, minDist, trials, directionDim, functions ):
    indent = 1e-6
    dim=3
    #lists for the model
    node_coords = []
    mechBC_merged = []
    mechInitC_merged = []
    govNodes = []
    govNodesMechBC = []
    rigidPlates = []

    lineSupported = True

    node_coords.append( np.array([height/2, 0, radius/2]))
    node_coords.append( np.array([height/2, 0, -radius/2]))


    if lineSupported:
        node_coords.append( np.array([height/2, radius/4, radius/4]))


        oldLen = len(node_coords)
        nodeA = np.array([indent, radius*0.99, radius*0.14])
        nodeB = np.array([height-indent, radius*0.99, radius*0.14])
        #pointGenerators.generateNodesLine3dRand(nodeA, nodeB, minDist/2, dim, node_coords, trials, catchCorners=True, equidist = True)
        nodeA = np.array([indent, radius*0.99, radius*0.07])
        nodeB = np.array([height-indent, radius*0.99, radius*0.07])
        pointGenerators.generateNodesLine3dRand(nodeA, nodeB, minDist/2, dim, node_coords, trials, catchCorners=True, equidist = True)
        nodeA = np.array([indent, radius*0.99, radius*0.0])
        nodeB = np.array([height-indent, radius*0.99, radius*0.00])
        pointGenerators.generateNodesLine3dRand(nodeA, nodeB, minDist/2, dim, node_coords, trials, catchCorners=True, equidist = True)
        nodeA = np.array([indent, radius*0.99, -radius*0.14])
        nodeB = np.array([height-indent, radius*0.99, -radius*0.14])
        #pointGenerators.generateNodesLine3dRand(nodeA, nodeB, minDist/2, dim, node_coords, trials, catchCorners=True, equidist = True)
        nodeA = np.array([indent, radius*0.99, -radius*0.07])
        nodeB = np.array([height-indent, radius*0.99, -radius*0.07])
        pointGenerators.generateNodesLine3dRand(nodeA, nodeB, minDist/2, dim, node_coords, trials, catchCorners=True, equidist = True)
        nodeA = np.array([indent, radius*0.99, -radius*0.0])
        nodeB = np.array([height-indent, radius*0.99, -radius*0.00])
        pointGenerators.generateNodesLine3dRand(nodeA, nodeB, minDist/2, dim, node_coords, trials, catchCorners=True, equidist = True)
        newNodes = len(node_coords)-oldLen

        """
        oldLen = len(node_coords)
        nodeA = np.array([indent, radius*0.99, radius*0.14])
        nodeB = np.array([height-indent, radius*0.99, radius*0.14])
        pointGenerators.generateNodesLine3dRand(nodeA, nodeB, minDist/2, dim, node_coords, trials, catchCorners=True, equidist = True)
        nodeA = np.array([indent, radius*0.99, radius*0.10])
        nodeB = np.array([height-indent, radius*0.99, radius*0.10])
        pointGenerators.generateNodesLine3dRand(nodeA, nodeB, minDist/2, dim, node_coords, trials, catchCorners=True, equidist = True)
        nodeA = np.array([indent, radius*0.99, radius*0.06])
        nodeB = np.array([height-indent, radius*0.99, radius*0.06])
        pointGenerators.generateNodesLine3dRand(nodeA, nodeB, minDist/2, dim, node_coords, trials, catchCorners=True, equidist = True)
        nodeA = np.array([indent, radius*0.99, radius*0.02])
        nodeB = np.array([height-indent, radius*0.99, radius*0.02])
        pointGenerators.generateNodesLine3dRand(nodeA, nodeB, minDist/2, dim, node_coords, trials, catchCorners=True, equidist = True)
        nodeA = np.array([indent, radius*0.99, -radius*0.02])
        nodeB = np.array([height-indent, radius*0.99, -radius*0.02])
        pointGenerators.generateNodesLine3dRand(nodeA, nodeB, minDist/2, dim, node_coords, trials, catchCorners=True, equidist = True)
        nodeA = np.array([indent, radius*0.99, -radius*0.06])
        nodeB = np.array([height-indent, radius*0.99, -radius*0.06])
        pointGenerators.generateNodesLine3dRand(nodeA, nodeB, minDist/2, dim, node_coords, trials, catchCorners=True, equidist = True)
        nodeA = np.array([indent, radius*0.99, -radius*0.10])
        nodeB = np.array([height-indent, radius*0.99, -radius*0.10])
        pointGenerators.generateNodesLine3dRand(nodeA, nodeB, minDist/2, dim, node_coords, trials, catchCorners=True, equidist = True)
        nodeA = np.array([indent, radius*0.99, -radius*0.14])
        nodeB = np.array([height-indent, radius*0.99, -radius*0.14])
        pointGenerators.generateNodesLine3dRand(nodeA, nodeB, minDist/2, dim, node_coords, trials, catchCorners=True, equidist = True)
        newNodes = len(node_coords)-oldLen
        """

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
        nodeA = np.array([indent, -radius*0.99, radius*0.14])
        nodeB = np.array([height-indent, -radius*0.99, radius*0.14])
        #pointGenerators.generateNodesLine3dRand(nodeA, nodeB, minDist/2, dim, node_coords, trials, catchCorners=True, equidist = True)
        nodeA = np.array([indent, -radius*0.99, radius*0.07])
        nodeB = np.array([height-indent, -radius*0.99, radius*0.07])
        pointGenerators.generateNodesLine3dRand(nodeA, nodeB, minDist/2, dim, node_coords, trials, catchCorners=True, equidist = True)
        nodeA = np.array([indent, -radius*0.99, radius*0.0])
        nodeB = np.array([height-indent, -radius*0.99, radius*0.00])
        pointGenerators.generateNodesLine3dRand(nodeA, nodeB, minDist/2, dim, node_coords, trials, catchCorners=True, equidist = True)
        nodeA = np.array([indent, -radius*0.99, -radius*0.14])
        nodeB = np.array([height-indent, -radius*0.99, -radius*0.14])
        #pointGenerators.generateNodesLine3dRand(nodeA, nodeB, minDist/2, dim, node_coords, trials, catchCorners=True, equidist = True)
        nodeA = np.array([indent, -radius*0.99, -radius*0.07])
        nodeB = np.array([height-indent, -radius*0.99, -radius*0.07])
        pointGenerators.generateNodesLine3dRand(nodeA, nodeB, minDist/2, dim, node_coords, trials, catchCorners=True, equidist = True)
        nodeA = np.array([indent, -radius*0.99, -radius*0.0])
        nodeB = np.array([height-indent, -radius*0.99, -radius*0.00])
        pointGenerators.generateNodesLine3dRand(nodeA, nodeB, minDist/2, dim, node_coords, trials, catchCorners=True, equidist = True)
        newNodes = len(node_coords)-oldLen
        """
        oldLen = len(node_coords)
        nodeA = np.array([indent, -radius*0.99, radius*0.14])
        nodeB = np.array([height-indent, -radius*0.99, radius*0.14])
        pointGenerators.generateNodesLine3dRand(nodeA, nodeB, minDist/2, dim, node_coords, trials, catchCorners=True, equidist = True)
        nodeA = np.array([indent, -radius*0.99, radius*0.10])
        nodeB = np.array([height-indent, -radius*0.99, radius*0.10])
        pointGenerators.generateNodesLine3dRand(nodeA, nodeB, minDist/2, dim, node_coords, trials, catchCorners=True, equidist = True)
        nodeA = np.array([indent, -radius*0.99, radius*0.06])
        nodeB = np.array([height-indent, -radius*0.99, radius*0.06])
        pointGenerators.generateNodesLine3dRand(nodeA, nodeB, minDist/2, dim, node_coords, trials, catchCorners=True, equidist = True)
        nodeA = np.array([indent, -radius*0.99, radius*0.02])
        nodeB = np.array([height-indent, -radius*0.99, radius*0.02])
        pointGenerators.generateNodesLine3dRand(nodeA, nodeB, minDist/2, dim, node_coords, trials, catchCorners=True, equidist = True)
        nodeA = np.array([indent, -radius*0.99, -radius*0.02])
        nodeB = np.array([height-indent, -radius*0.99, -radius*0.02])
        pointGenerators.generateNodesLine3dRand(nodeA, nodeB, minDist/2, dim, node_coords, trials, catchCorners=True, equidist = True)
        nodeA = np.array([indent, -radius*0.99, -radius*0.06])
        nodeB = np.array([height-indent, -radius*0.99, -radius*0.06])
        pointGenerators.generateNodesLine3dRand(nodeA, nodeB, minDist/2, dim, node_coords, trials, catchCorners=True, equidist = True)
        nodeA = np.array([indent, -radius*0.99, -radius*0.10])
        nodeB = np.array([height-indent, -radius*0.99, -radius*0.10])
        pointGenerators.generateNodesLine3dRand(nodeA, nodeB, minDist/2, dim, node_coords, trials, catchCorners=True, equidist = True)
        nodeA = np.array([indent, -radius*0.99, -radius*0.14])
        nodeB = np.array([height-indent, -radius*0.99, -radius*0.14])
        pointGenerators.generateNodesLine3dRand(nodeA, nodeB, minDist/2, dim, node_coords, trials, catchCorners=True, equidist = True)
        newNodes = len(node_coords)-oldLen
        """

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



    ###############generating of points supported surface left face ###############
    pointGenerators.generateNodesOrtoCircleBorder3dRand(center, radius-1e-5, directionDim, minDist, node_coords, trials)
    pointGenerators.generateNodesOrtoCircle3dRand(center, radius-1e-5, directionDim, minDist, node_coords, trials)

    ###############generating of points loaded surface right face ###############
    nodeA = center.copy()
    nodeA[directionDim] += height-indent
    pointGenerators.generateNodesOrtoCircleBorder3dRand(nodeA, radius-1e-5,  directionDim, minDist, node_coords, trials)
    pointGenerators.generateNodesOrtoCircle3dRand(nodeA, radius-1e-5,  directionDim, minDist, node_coords, trials)

    ###############generating of points cylinder surf###############
    pointGenerators.generateNodesOrtoCilinderSurf3dRand(center, radius-1e-5, height, directionDim, minDist,  node_coords, trials)

    ###############generating of points cylinder volume ###############
    pointGenerators.generateNodesOrtoCilinder3dRand(center, radius-1e-5, height, directionDim, minDist,  node_coords, trials)
    #######################################################################


    """
    node_coords = np.asarray(node_coords)
    if SHOW_PLOT:
        fig = plt.figure()
        ax = Axes3D(fig)
        ax.scatter(node_coords[:,0], node_coords[:,1], node_coords[:,2])
        plt.show()
    #"""




    return node_coords, mechBC_merged, mechInitC_merged, govNodes, govNodesMechBC, rigidPlates

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
    print ('Assembling Biparva tube...', end='')
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


    mechBC = np.array([0,0,0,0,0,0,    -1,-1,-1,-1,-1,-1])
    node_coords.append( np.array([  indent,  radius-thickness/2,  0 ]))
    mBC = utilitiesMech.mechanicalBC(dim, 0, mechBC)
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
    pointGenerators.generateNodesOrtoTube3dRand(center, radius-1e-5, height, thickness, directionDim, minDist,  node_coords, trials)
    #######################################################################


    return node_coords, mechBC_merged,  govNodes, govNodesMechBC, rigidPlates




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
