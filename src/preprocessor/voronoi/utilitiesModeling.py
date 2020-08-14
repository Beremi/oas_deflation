import numpy as np
import random
import utilitiesGeom
import utilitiesMech
import utilitiesNumeric
import pointGenerators
import voronoi
import matplotlib.pyplot as plt
#import voronoi_viewer
from mpl_toolkits.mplot3d import Axes3D

from scipy.spatial import Voronoi
from scipy.spatial import voronoi_plot_2d
from scipy.spatial import Delaunay

def assembleMeasuringGauges(type, D=-1, maxLim = None):
    measuringGauges = []
    if (type == 'dogbone2d'):
        if (D==0.1):
            #total length LS
            coordsA = np.array([ D/2, 0])
            coordsB = np.array([ D/2, 6/4*D ])
            measuringGauges.append(utilitiesMech.MeasuringGauge(coordsA, coordsB, 'totalLS', False))

            #mid LS
            coordsA = np.array([ D/2, 3/4*D-0.075/2 ])
            coordsB = np.array([ D/2, 3/4*D+0.075/2 ])
            measuringGauges.append(utilitiesMech.MeasuringGauge(coordsA, coordsB, 'midLS', False))
            #left LC
            coordsA = np.array([ 0.2*D, 3/4*D-0.075/2 ])
            coordsB = np.array([ 0.2*D, 3/4*D+0.075/2 ])
            measuringGauges.append(utilitiesMech.MeasuringGauge(coordsA, coordsB, 'leftLS', False))
            #right LC
            coordsA = np.array([ D-0.2*D, 3/4*D-0.075/2 ])
            coordsB = np.array([ D-0.2*D, 3/4*D+0.075/2 ])
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

    if(type=='cylinder3d'):
        #total length
        coordsA = np.array([ 0, 0, 0])
        coordsB = np.array([ maxLim[0], 0, 0 ])
        measuringGauges.append(utilitiesMech.MeasuringGauge(coordsA, coordsB, 'total', True))

    if(type=='3pb2d'):
        coordsA = np.array([ maxLim[0]/2, maxLim[1]])
        coordsB = np.array([ 0, 0 ])
        measuringGauges.append(utilitiesMech.MeasuringGauge(coordsA, coordsB, 'totalLS', False))

    if(type=='3pb3d'):
        coordsA = np.array([ maxLim[0]/2, maxLim[1], maxLim[2]])
        coordsB = np.array([ 0, 0, 0 ])
        measuringGauges.append(utilitiesMech.MeasuringGauge(coordsA, coordsB, 'totalLS', False))

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


def assembleMaterialZones (elaX, dim, model='box', maxLim=None, D=None, thickness=None):
    materialZones = []
    #matZone 1
    matZ = []
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
            boundA1 = np.array(  [ -1e-8, 6/4*D+1e-8] )
            matZ.append (boundA1)
            boundB1 = np.array(  [ D  , 6/4*D - elaX]  )
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

    #fig = voronoi.voronoi_plot_2d(vor, show_vertices = True)
    #plt.show()

    print(node_coords)
    node_coords = np.asarray(node_coords)
    #fig, ax = plt.subplots()
    #ax.scatter(node_coords[:,0], node_coords[:,1])
    #plt.show()

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





def create2dSSBeamUnifLoad(maxLim, minDist, trials, notch = -1, loadWidth = 1, fracZoneWidth = 0.15,  orthogonalFracZone=False):
    print('Creating 2d simply supported beam, uniform load.')
    #
    node_coords, mechBC_merged, mechInitC_merged, notches, govNodes, govNodesMechBC, rigidPlates  = assemble2DSSBeamBending(maxLim, minDist, trials, notch, loadWidth, fracZoneWidth, orthogonalFracZone=orthogonalFracZone);

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



    return node_coords, mechBC_merged, mechInitC_merged,  vor, areas, functions, notches, govNodes, govNodesMechBC, rigidPlates


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

    #fig = voronoi_plot_2d(vor, show_vertices=True, line_colors='orange',  line_width=2, line_alpha=0.6, point_size=2)
    #plt.show()

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



def create2dPeriodicShear(maxLim, minDist, trials ):
    print('Creating 2d periodic rectangle, shear loaded.')
    ### sampling of nodes
    ### direct setting of mechanicalBCs
    node_coords, mechBC_merged, mechInitC_merged, nodePositions, coupledNodes, mirtype = asssemble2dPeriodicShear(maxLim, minDist, trials );

    print ('Conducting Voronoi tesselation...', end ='')
    vor = Voronoi(node_coords)
    regions, vertices, polygons, areas, centroids, points = voronoi.voronoi_2d(vor, maxLim)
    print('done.')


    #fig = voronoi_plot_2d(vor, show_vertices=False, line_colors='orange',  line_width=2, line_alpha=0.6, point_size=2)
    #plt.show()

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

    return node_coords, mechBC_merged, mechIC_merged, transportBC_merged, transportIC_merged, vor, areas, functions, nodePositions, coupledNodes, mirtype

    ########################################################################
    ### indirect setting of transportBCs by spatial selection of vertices
    transportBC_merged = []
    functions = []

    ### selecting vertices on the left surface
    boundA = np.array(  [-1e-8 , -1e-8] )
    boundB = np.array(  [maxLim[0]+1e-8, maxLim[1]+1e-8]  )
    faces1 = utilitiesGeom.returnSelectedPts(boundA, boundB, vor.vertices)
    vert = vor.vertices[faces1,:]
    boundA = np.array(  [1e-8 , 1e-8] )
    boundB = np.array(  [maxLim[0]-1e-8, maxLim[1]-1e-8]  )
    faces0 = utilitiesGeom.excludeSelectedPts(boundA, boundB, vert)
    faces = faces1[faces0]

    for i,k in enumerate(faces):
        fn1 = utilitiesNumeric.constantFunc(np.sin(vor.vertices[k,0])*np.exp(vor.vertices[k,1]))
        functions.append (fn1)
        trsBC = utilitiesMech.transportBC(k,[i,-1])
        transportBC_merged.append(trsBC)

    return node_coords, [], transportBC_merged, vor, areas, functions


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






def create2dDogBone(minDist, trials, D=1.0, excentricity = 20 ):
    print('Creating 2d dog bone....')
    #
    node_coords, mechBC_merged, mechInitC_merged, node_count, govNodes, govNodesMechBC, rigidPlates  = assemble2dDogBone(D, minDist, trials, excentricity = excentricity);

    node_coords = np.asarray(node_coords)
    """
    fig, ax = plt.subplots()
    ax.scatter(node_coords[:,0], node_coords[:,1])
    plt.show()
    """

    print('Conducting Voronoi tesselation...', end = '')
    vor = utilitiesNumeric.runMirroredVoronoiDogBone(node_coords, 2, D)
    print('done.')

    node_coords = node_coords[0:node_count]
    areas = []
    for i in range (node_count): areas.append(0)
    areas = np.asarray(areas)

    """
    fig = voronoi_plot_2d(vor, show_vertices=False, line_colors='orange',line_width=2, line_alpha=0.6, point_size=2)
    plt.show()
    """
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


def create3dDogBone(minDist, trials, D=1.0, excentricity = 20 ):
    print('Creating 3sd dog bone....')
    #
    node_coords, mechBC_merged, mechInitC_merged, node_count,govNodes, govNodesMechBC, rigidPlates  = assemble3dDogBone(D, minDist, trials, excentricity = excentricity);

    node_coords = np.asarray(node_coords)
    """
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

    #fig = voronoi_plot_2d(vor, show_vertices=False, line_colors='orange',line_width=2, line_alpha=0.6, point_size=2)
    #splt.show()

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






def create3dSSBeamUnifLoad(maxLim, minDist, trials, notch = -1, loadWidth = 1, fracZoneWidth = 0.15, orthogonalFracZone = False ):
    print('Creating 3d simply supported beam, uniform load.')
    #govNodes, rigidPlates
    node_coords, mechBC_merged, mechInitC_merged, notches, govNodes, govNodesMechBC, rigidPlates  = assemble3DSSBeamBending(maxLim, minDist, trials, notch, loadWidth, fracZoneWidth=fracZoneWidth, orthogonalFracZone=orthogonalFracZone);
    node_coords = np.asarray(node_coords)
    """
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




    return node_coords, mechBC_merged, mechInitC_merged,  vor, volumes, functions, notches, govNodes, govNodesMechBC, rigidPlates



def create3dReinhardtTension(maxLim, minDist, trials, fracZoneWidth = 0.15 ):
    print('Creating 3d simply supported beam, uniform load.')
    #govNodes, rigidPlates
    node_coords, mechBC_merged, mechInitC_merged, notchNodes, govNodes, govNodesMechBC, rigidPlates  = assemble3DReinhardtTension(maxLim, minDist, trials, fracZoneWidth=fracZoneWidth);
    node_coords = np.asarray(node_coords)
    """
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


def create3dcylinderTorsionPressFree(center, radius, height, minDist, trials, directionDim ):

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



    return node_coords, mechBC_merged,  vor, volumes, functions,govNodes, govNodesMechBC, rigidPlates




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

    oldLen = len(node_coords)
    #utilitiesGeom.generateNodesLine2dRand(nodeA, nodeB, minDist, dim, node_coords, trials, False, False)
    pointGenerators.generateSingleNode (nodeA, dim, node_coords)
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



def assemble2DSSBeamBending (maxLim, minDist, trials, notch, loadWidth, fracZoneWidth,  orthogonalFracZone=False):
    dim = 2
    #lists for the model
    node_coords = []
    mechBC_merged = []
    mechInitC_merged = []
    govNodes = []
    govNodesMechBC = []
    rigidPlates = []

    #an indent due to mirroring of the data for voronoi tess.
    notches=[]
    indent = 1e-8
    #notchWidth = 1.5e-3 /2
    notchWidth = minDist/2
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
    supportWidth = maxLim[0] / 40


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
    #nodeA =  np.array([indent + 0.5*maxLim[0]*(1-loadWidth) , maxLim[1] - indent])
    #nodeB =  np.array([maxLim[0] - indent - 0.5*maxLim[0]*(1-loadWidth) , maxLim[1] - indent])
    #pointGenerators.generateNodesLine2dRand(nodeA, nodeB, minDist/2, dim, node_coords,  trials, True, True)
    #######top of frac zone
    nodeA =  np.array([0.5*maxLim[0] - 0.5*maxLim[1]*(1-notch) , maxLim[1] - indent])
    nodeB =  np.array([0.5*maxLim[0] + 0.5*maxLim[1]*(1-notch), maxLim[1] - indent])
    pointGenerators.generateNodesLine2dRand(nodeA, nodeB, minDist/2, dim, node_coords,  trials, True, True)


    ########################################## rest of  faces
    nodeA =  np.array([indent + 0.5*maxLim[0]*(1-loadWidth) , maxLim[1] - indent])
    nodeB =  np.array([indent, maxLim[1] - indent])
    pointGenerators.generateNodesLine2dRand(nodeA, nodeB, minDist*1.5, dim, node_coords,  trials, True, True)
    nodeA =  np.array([maxLim[0] - indent  , maxLim[1] - indent])
    nodeB =  np.array([maxLim[0] - indent - 0.5*maxLim[0]*(1-loadWidth) , maxLim[1] - indent])
    pointGenerators.generateNodesLine2dRand(nodeA, nodeB, minDist*1.5, dim, node_coords,  trials, True, True)
    nodeA =  np.array([indent  ,  indent])
    nodeB =  np.array([maxLim[0]/2-notchWidth ,  indent])
    pointGenerators.generateNodesLine2dRand(nodeA, nodeB, minDist*1.5, dim, node_coords,  trials, False, True)
    nodeA =  np.array([maxLim[0]-indent  ,  indent])
    nodeB =  np.array([maxLim[0]/2+notchWidth ,  indent])
    pointGenerators.generateNodesLine2dRand(nodeA, nodeB, minDist*1.5, dim, node_coords,  trials, False, True)
    nodeA =  np.array([indent, indent])
    nodeB =  np.array([indent, maxLim[1] - indent])
    pointGenerators.generateNodesLine2dRand(nodeA, nodeB, minDist*1.5, dim, node_coords,  trials, False, False)
    nodeA =  np.array([maxLim[0]-indent, indent])
    nodeB =  np.array([maxLim[0]-indent, maxLim[1] - indent])
    pointGenerators.generateNodesLine2dRand(nodeA, nodeB, minDist*1.5, dim, node_coords,  trials, False, False)

    """
    nodeA =  np.array([indent  ,  indent])
    nodeB =  np.array([indent ,  maxLim[1]-indent])
    pointGenerators.generateNodesLine2dRand(nodeA, nodeB, minDist/1, dim, node_coords,  trials, False, False)
    nodeA =  np.array([maxLim[0] - indent   ,  indent])
    nodeB =  np.array([maxLim[0] - indent  ,  maxLim[1]-indent])
    pointGenerators.generateNodesLine2dRand(nodeA, nodeB, minDist/1, dim, node_coords,  trials, False, False)
    #"""
    ##########################################generating of points, fracture zone
    """
    maxLimF = np.array([
    maxLim[0] - indent - 0.5*maxLim[0]*(1-fracZoneWidth),
    maxLim[1] - indent,
    indent + 0.5*maxLim[0]*(1-fracZoneWidth),
    maxLim[1]*notch*0.8])
    #"""

    maxLimF = np.array([
    maxLim[0]*0.5 - maxLim[1]*(1-notch)/2,
    maxLim[1] - indent,
    maxLim[0]*0.5 + maxLim[1]*(1-notch)/2,
    maxLim[1]*notch])


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
    pointGenerators.generateNodesRect(maxLimF, minDist*1.5, dim, trials, node_coords, useLowBound=True)

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
    pointGenerators.generateNodesRect(maxLim, minDist*3, dim, trials, node_coords)


    return node_coords, mechBC_merged, mechInitC_merged, notches, govNodes, govNodesMechBC, rigidPlates





def assemble2dDogBone(D, minDist, trials, excentricity = 20):
    dim = 2
    #lists for the model
    node_coords = []
    mechBC_merged = []
    mechInitC_merged = []
    govNodes = []
    govNodesMechBC = []
    rigidPlates = []

    indent = 1e-5

    #####################nodes of interest
    node_coords.append(np.array([  D/2,  indent   ])) #top mid
    node_coords.append(np.array([  D/2,  6/4*D - indent  ]))  #bottom mid
    #gauges B
    if (D==0.1):
        node_coords.append( np.array([ D/2,        3/4*D-0.075/2 ])  )#mid LS
        node_coords.append( np.array([ D/2,        3/4*D+0.075/2 ])  )
        node_coords.append( np.array([ 0.2*D,  3/4*D-0.075/2 ])  ) #left LC
        node_coords.append( np.array([ 0.2*D,  3/4*D+0.075/2 ])  )
        node_coords.append( np.array([ D-0.2*D,  3/4*D-0.075/2 ])  )#right LC
        node_coords.append( np.array([ D-0.2*D,  3/4*D+0.075/2 ])  )
    #gauges D
    if (D==0.4):
        node_coords.append( np.array([ D/2,        3/4*D-0.240/2 ])  )#mid LS
        node_coords.append( np.array([ D/2,        3/4*D+0.240/2 ])  )
        node_coords.append( np.array([ D/2-0.2*D,  3/4*D-0.240/2 ])  )#left LC
        node_coords.append( np.array([ D/2-0.2*D,  3/4*D+0.240/2 ])  )
        node_coords.append( np.array([ D/2+0.2*D,  3/4*D-0.240/2 ])  )#right LC
        node_coords.append( np.array([ D/2+0.2*D,  3/4*D+0.240/2 ])  )

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


    #top line of dogbone
    nodeA = np.array([indent, indent])
    nodeB = np.array([indent+D, indent])
    pointGenerators.generateNodesLine2dRand(nodeA, nodeB, minDist*0.9, dim, node_coords, trials, True, False)


    #bottom line of dogbone
    nodeA = np.array([indent,  6/4 * D - indent])
    nodeB = np.array([D-indent, 6/4 * D - indent])
    pointGenerators.generateNodesLine2dRand(nodeA, nodeB, minDist*0.9, dim, node_coords, trials, True, False)


    #sampling on circular borders
    centreA = np.array( [-0.525 * D, 3/4 * D] )
    centreB = np.array( [ 1.525 * D, 3/4 * D] )
    radius = 0.725*D
    angleLimit =   np.arcsin( 0.5*D / radius)
    mirroredPointsA = pointGenerators.generateNodesCircle2dRand(centreA, radius+indent, minDist*0.9, node_coords, trials, angleLimitA=-angleLimit, angleLimitB=angleLimit, mirrorIndent = indent )
    mirroredPointsB = pointGenerators.generateNodesCircle2dRand(centreB, radius+indent, minDist*0.9, node_coords, trials, angleLimitA=np.pi-angleLimit, angleLimitB=np.pi+angleLimit, mirrorIndent = indent)


    #rectangle of dogbone
    ectBC = np.array([-1,-1,-1,-1,-1,-1])
    oldLen = len(node_coords)
    maxLim = np.array([  D    ,  6/4*D ])
    pointGenerators.generateNodesRect(maxLim, minDist, 2, trials, node_coords)
    nrOfPoints =  (len(node_coords)) - oldLen

    #dumping points outside bone

    radius = np.linalg.norm( centreB - np.array([D, 1/4*D]))
    print('Dumping points within bordering circles...', end='')
    node_coords_out = []
    for node in node_coords:
        distA = np.linalg.norm( node - centreA)
        distB = np.linalg.norm( node - centreB)
        if (distA > radius and distB > radius):
            node_coords_out.append(node)
    print('done.')

    #mirrored circles. Not to be in the model at the end
    node_count = len(node_coords_out)
    node_coords_out = np.vstack( (node_coords_out, mirroredPointsA, mirroredPointsB) )



    return node_coords_out, mechBC_merged, mechInitC_merged, node_count, govNodes, govNodesMechBC, rigidPlates




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





def asssemble2dPeriodicShear (maxLim, minDist, trials):
    dim = 2
    #lists for the model
    node_coords = []
    mechBC_merged = []
    mechInitC_merged = []

    periodicBand = 3 * minDist

    ###########generating of points in rectangle
    pointGenerators.generateNodesRectPeriodic(maxLim, minDist, dim, trials, node_coords)

    #np.savetxt('test.out', np.asarray(node_coords), delimiter='\t')
    #node_coords = np.loadtxt('test.out')

    mirtype = []
    coupledNodes = []
    nodePositions = []
    for i in range (len(node_coords)):
        nodePositions.append(i+1)
        mirtype.append(0)
    node_coords = np.asarray(node_coords)
    print(len(node_coords))

    #plt.plot(node_coords[:,0]+1, node_coords[:,1], 'x', color='red');
    ########### adding periodic points
    print('Adding periodic points')
    for i in range (len(node_coords)):
        point = np.asarray(node_coords[i,:])
        #
        xplus = False
        xminus = False
        yplus = False
        yminus = False
        #
        if(point[0]<periodicBand):
            xplus = True
        if(point[0]>maxLim[0]-periodicBand):
            xminus = True
        if(point[1]<periodicBand):
            yplus = True
        if(point[1]>maxLim[1]-periodicBand):
            yminus = True

        k = i+1# len(node_coords)+1

        #images due to one dimension
        if xplus:
            #k = len(node_coords)+1
            newPoint = np.copy(point)
            newPoint[0] += maxLim[0]
            nodePositions.append(-k)
            node_coords = np.vstack(( node_coords, newPoint ))
            coupledNodes.append( np.array([i, len(node_coords)-1]) )
            mirtype.append(1)
        if xminus:
            newPoint = np.copy(point)
            newPoint[0] -= maxLim[0]
            nodePositions.append(0)
            node_coords = np.vstack((node_coords, newPoint))
            mirtype.append(-1)
        if yplus:
            #k = len(node_coords)+1
            newPoint = np.copy(point)
            newPoint[1] += maxLim[1]
            nodePositions.append(-k)
            #plt.plot(  newPoint[0] , newPoint[1] ,'o', color='red')
            node_coords = np.vstack((node_coords, newPoint))
            coupledNodes.append( np.array([i, len(node_coords)-1]) )
            mirtype.append(2)
        if yminus:
            newPoint = np.copy(point)
            newPoint[1] -= maxLim[1]
            nodePositions.append(0)
            node_coords = np.vstack((node_coords, newPoint))
            mirtype.append(-1)

        if xplus and yplus:
            #k = len(node_coords)+1
            newPoint = np.copy(point)
            newPoint[0] += maxLim[0]
            newPoint[1] += maxLim[1]
            nodePositions.append(-k)
            node_coords = np.vstack((node_coords, newPoint))
            coupledNodes.append( np.array([i, len(node_coords)-1]) )
            mirtype.append(3)

        if xminus and yplus:
            newPoint = np.copy(point)
            newPoint[0] -= maxLim[0]
            newPoint[1] += maxLim[1]
            nodePositions.append(0)
            node_coords = np.vstack((node_coords, newPoint))
            mirtype.append(-1)
        if xplus and yminus:
            newPoint = np.copy(point)
            newPoint[0] += maxLim[0]
            newPoint[1] -= maxLim[1]
            nodePositions.append(0)
            node_coords = np.vstack((node_coords, newPoint))
            mirtype.append(-1)

        if xminus and yminus:
            newPoint = np.copy(point)
            newPoint[0] -= maxLim[0]
            newPoint[1] -= maxLim[1]
            nodePositions.append(0)
            node_coords = np.vstack((node_coords, newPoint))
            mirtype.append(-1)

    #plt.plot(node_coords[:,0], node_coords[:,1], 'o', color='black');

    #for i in range (len(coupledNodes)):
    #    plt.plot( node_coords[ abs(coupledNodes[i][0]),0 ] , node_coords[ abs(coupledNodes[i][0]),1 ] ,'o', color='red')
    #    plt.plot( node_coords[ abs(coupledNodes[i][1]),0 ] , node_coords[ abs(coupledNodes[i][1]),1 ] ,'o', color='green')

    #plt.show()

    #print (len(node_coords))

    #np.savetxt('test.out', np.asarray(node_coords), delimiter='\t')
    #node_coords = np.loadtxt('test.out')



    return node_coords, mechBC_merged, mechInitC_merged, nodePositions, coupledNodes, mirtype







def assemble3DSSBeamBending (maxLim, minDist, trials, notch, loadWidth,  fracZoneWidth = 0.15, orthogonalFracZone=False):
    dim = 3
    #lists for the model
    node_coords = []
    mechBC_merged = []
    mechInitC_merged = []
    govNodes = []
    govNodesMechBC = []
    rigidPlates = []

    node_coords.append( np.array([maxLim[0]/4, maxLim[1]/2, maxLim[2]/2]))
    #lineBC = np.array([0,0,0, 0,0,0,  -1,-1,-1,-1,-1,-1])
    #mBC = utilitiesMech.mechanicalBC(dim, 0, lineBC)
    #mechBC_merged.append(mBC)

    #an indent due to mirroring of the data for voronoi tess.
    notches=[]
    indent = 1e-8
    notchWidth = 5e-3/2
    #generating notch points
    if (notch > 0):
        notchSide0 = []
        oldLen = len(node_coords)
        nodeA = np.array([maxLim[0]/2-notchWidth, indent, indent])
        nodeB = np.array([maxLim[0]/2-notchWidth, indent, maxLim[2]-indent])
        pointGenerators.generateNodesLine3dRand(nodeA, nodeB, minDist/4, dim, node_coords, trials, catchCorners=True, equidist=True)
        nodeA = np.array([maxLim[0]/2-notchWidth, maxLim[1]*notch, indent])
        nodeB = np.array([maxLim[0]/2-notchWidth, maxLim[1]*notch, maxLim[2]-indent])
        pointGenerators.generateNodesLine3dRand(nodeA, nodeB, minDist/4, dim, node_coords, trials, catchCorners=True, equidist=True)
        nodeA = np.array([maxLim[0]/2-notchWidth, indent, indent])
        nodeB = np.array([maxLim[0]/2-notchWidth, maxLim[1]*notch, indent])
        pointGenerators.generateNodesLine3dRand(nodeA, nodeB, minDist/4, dim, node_coords, trials, catchCorners=False, equidist=True)
        nodeA = np.array([maxLim[0]/2-notchWidth, indent, maxLim[2]-indent])
        nodeB = np.array([maxLim[0]/2-notchWidth, maxLim[1]*notch, maxLim[2]-indent])
        pointGenerators.generateNodesLine3dRand(nodeA, nodeB, minDist/4, dim, node_coords, trials, catchCorners=False, equidist=True)

        nodeA = np.array([maxLim[0]/2-notchWidth, indent, indent])
        nodeB = np.array([maxLim[0]/2-notchWidth, maxLim[1]*notch, maxLim[2]-indent])
        pointGenerators.generateNodesOrtoSurface3dRand(nodeA, nodeB, minDist/4, dim, node_coords, trials,minDistAmongNewPoints=True)
        for i in range (oldLen, len(node_coords), 1):
            notchSide0.append(i)


        notchSide1 = []
        oldLen = len(node_coords)
        nodeA = np.array([maxLim[0]/2+notchWidth, indent, indent])
        nodeB = np.array([maxLim[0]/2+notchWidth, indent, maxLim[2]-indent])
        pointGenerators.generateNodesLine3dRand(nodeA, nodeB, minDist/4, dim, node_coords, trials, catchCorners=True, equidist=True)
        nodeA = np.array([maxLim[0]/2+notchWidth, maxLim[1]*notch, indent])
        nodeB = np.array([maxLim[0]/2+notchWidth, maxLim[1]*notch, maxLim[2]-indent])
        pointGenerators.generateNodesLine3dRand(nodeA, nodeB, minDist/4, dim, node_coords, trials, catchCorners=True, equidist=True)
        nodeA = np.array([maxLim[0]/2+notchWidth, indent, indent])
        nodeB = np.array([maxLim[0]/2+notchWidth, maxLim[1]*notch, indent])
        pointGenerators.generateNodesLine3dRand(nodeA, nodeB, minDist/4, dim, node_coords, trials, catchCorners=False, equidist=True)
        nodeA = np.array([maxLim[0]/2+notchWidth, indent, maxLim[2]-indent])
        nodeB = np.array([maxLim[0]/2+notchWidth, maxLim[1]*notch, maxLim[2]-indent])
        pointGenerators.generateNodesLine3dRand(nodeA, nodeB, minDist/4, dim, node_coords, trials, catchCorners=False, equidist=True)

        nodeA = np.array([maxLim[0]/2+notchWidth, indent, indent])
        nodeB = np.array([maxLim[0]/2+notchWidth, maxLim[1]*notch, maxLim[2]-indent])
        pointGenerators.generateNodesOrtoSurface3dRand(nodeA, nodeB, minDist/4, dim, node_coords, trials, minDistAmongNewPoints= True)
        for i in range (oldLen, len(node_coords), 1):
            notchSide1.append(i)

        notchL = []
        notchL.append(notchSide0)
        notchL.append(notchSide1)
        notches.append(notchL)

    """
    node_coords = np.asarray(node_coords)
    fig = plt.figure()
    ax = Axes3D(fig)
    ax.scatter(node_coords[:,0], node_coords[:,1], node_coords[:,2])
    plt.show()
    """

    #width of the supports
    supportWidth = maxLim[0] / 40

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

    ############### loaded top face ###############
    lineBC = np.array([-1,-1,-1,-1,-1,-1,  -1, 1,-1,-1,-1,-1])
    nodeA =  np.array([indent + 0.5*maxLim[0]*(1-loadWidth), maxLim[1] - indent, indent])
    nodeB =  np.array([maxLim[0] - indent - 0.5*maxLim[0]*(1-loadWidth), maxLim[1] - indent, maxLim[2] - indent])
    pointGenerators.generateNodesOrtoSurface3dRand(nodeA, nodeB, minDist/2, dim, node_coords, trials)

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


    ##########################################generating of points, fracture zone
    maxLimF = np.array([
    maxLim[0] - indent - 0.5*maxLim[0]*(1-fracZoneWidth),
    maxLim[1] - indent,
    maxLim[2] - indent,
    indent + 0.5*maxLim[0]*(1-fracZoneWidth),
    maxLim[1]*notch*0.8+indent,
    indent])
    if not orthogonalFracZone:
        pointGenerators.generateNodesRect(maxLimF, minDist/2, dim, trials, node_coords, useLowBound=True)
    else:
        pointGenerators.generateOrtogrid(maxLimF, minDist/2, dim, node_coords)


    ##########################################generating of points, fracture zone
    maxLimF = np.array([
    maxLim[0] - indent - 0.5*maxLim[0]*(1-fracZoneWidth),
    maxLim[1]*notch*0.8 -indent,
    maxLim[2] - indent,
    indent + 0.5*maxLim[0]*(1-fracZoneWidth),
    indent,
    indent])
    pointGenerators.generateNodesRect(maxLimF, minDist/2, dim, trials, node_coords, useLowBound=True)





    ##########################################generating of points, left support
    maxLimF = np.array([
    supportWidth*2,
    supportWidth*2,
    maxLim[2]-indent,
    indent,
    indent,
    indent])
    pointGenerators.generateNodesRect(maxLimF, minDist/2, dim, trials, node_coords, useLowBound=True)
    ##########################################generating of points, right support
    maxLimF = np.array([
    maxLim[0],
    supportWidth*2,
    maxLim[2]-indent,
    maxLim[0]-supportWidth*2,
    indent,
    indent])
    pointGenerators.generateNodesRect(maxLimF, minDist/2, dim, trials, node_coords, useLowBound=True)



    ##########################################generating of points, homogeneous volume
    pointGenerators.generateNodesRect(maxLim, minDist, dim, trials, node_coords)

    return node_coords, mechBC_merged, mechInitC_merged, notches, govNodes, govNodesMechBC, rigidPlates





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
    rightRigidPlateMechBC = np.array([1,-1,-1, -1,0,0,  -1,-1,-1, -1,-1,-1])
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
    fig = plt.figure()
    ax = Axes3D(fig)
    ax.scatter(node_coords[:,0], node_coords[:,1], node_coords[:,2])
    plt.show()
    #"""




    return node_coords, mechBC_merged,  govNodes, govNodesMechBC, rigidPlates


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
