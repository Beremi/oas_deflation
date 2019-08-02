import numpy as np
import random
import utilitiesGeom
import utilitiesMech
import utilitiesNumeric
import pointGenerators
import voronoi
import matplotlib.pyplot as plt
import voronoi_viewer
from mpl_toolkits.mplot3d import Axes3D

def createSingleSpringTestModel(length):
    print ('Creating single spring test model.')
    #defining functions
    functions = []
    #### Defining functions
    #0 sine func
    #fn = utilitiesNumeric.sineFunc(10,22)

    #constant zero
    fn = utilitiesNumeric.constantFunc(0)
    functions.append (fn)

    fn1 = utilitiesNumeric.constantFunc(1e-3)
    functions.append (fn1)
    """
    fn2 = utilitiesNumeric.sawToothConstFunc(value = 20, period = 11, sym =1)
    functions.append (fn2)

    fn3 = utilitiesNumeric.sineFunc(value = 22, period = 33)
    functions.append (fn3)

    fn4 = utilitiesNumeric.varyingSawToothFunction(fn1, fn2)
    functions.append(fn4)
    """
    dim = 2
    idt = length /2

    maxLim = np.array([  length*2    ,   2  ])

    node_coords,  mechBC_merged = assembleTwoNodeSpringTest(maxLim, idt)

    #print(node_coords)
    print('Conducting Voronoi tesselation...', end = '')
    vor, regions, vertices, polygons, areas, centroids, points = utilitiesNumeric.runMirroredVoronoi (node_coords, 2, maxLim)
    print('done.')

    #fig = voronoi.voronoi_plot_2d(vor, show_vertices = True)
    #plt.show()

    idt = 1e-1
    ### indirect setting of transportBCs by spatial selection of vertices
    transportBC_merged = []
    ### selecting all vertices in model
    noTrsprtBC = np.array([ -1 , -1 ])
    boundA = np.array(  [ -idt , -idt ] )
    boundB = np.array(  [ maxLim[0] + idt , maxLim[1] + idt  ]  )
    allVrtcs = utilitiesGeom.returnSelectedPts(boundA, boundB, vor.vertices)

    for i in range (len(allVrtcs)):
        trsBC = utilitiesMech.transportBC(allVrtcs[i], noTrsprtBC)
        transportBC_merged.append(trsBC)




    return node_coords, mechBC_merged, transportBC_merged, vor, areas, functions



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





def create2dSSBeamUnifLoad(maxLim, minDist, trials ):
    print('Creating 2d simply supported beam, uniform load.')
    #
    node_coords, mechBC_merged, mechInitC_merged  = assemble2DSSBeamBending(maxLim, minDist, trials );

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

    #function from txt data_table
    #funcFromTxt = utilitiesNumeric.PWLFuncFromTxt('data_table.txt')
    #functions.append(funcFromTxt)

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


    return node_coords, mechBC_merged, mechInitC_merged, transportBC_merged, transportIC_merged, vor, areas, functions


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
    fn1 = utilitiesNumeric.constantFunc(1e-3)#utilitiesNumeric.sawToothConstFunc(value = -0.5e-4, period = 10)
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
















def assembleTwoNodeSpringTest (maxLim, idt):
    node_coords = []
    mechBC_merged = []
    dim = 2

    nodeA = np.array ( [ 0 + idt , maxLim[1]/2 ] )
    nodeAmechBC = np.array([0, 0 , -1 ,      -1 , -1 , -1])
    pointGenerators.generateSingleNode(nodeA, dim, node_coords)
    mBC = utilitiesMech.mechanicalBC(dim, 0, nodeAmechBC)
    mechBC_merged.append(mBC)

    nodeB = np.array ( [ maxLim[0] - idt, maxLim[1]/2 ] )
    nodeBmechBC = np.array([-1, 0 , 0 ,      1, -1 , -1])
    pointGenerators.generateSingleNode(nodeB, dim, node_coords)
    mBC = utilitiesMech.mechanicalBC(dim, 1, nodeBmechBC)
    mechBC_merged.append(mBC)


    return node_coords,  mechBC_merged



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





def create3dcylinderUniPress(center, radius, height, minDist, trials, directionDim ):
    ### sampling of nodes
    ### direct setting of mechanicalBCs
    node_coords, mechBC_merged, mechIC_merged  = assemble3dcylinderUniPress(center, radius, height, minDist, trials, directionDim )

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



def assemble2DSSBeamBending (maxLim, minDist, trials):
    dim = 2
    #lists for the model
    node_coords = []
    mechBC_merged = []
    mechInitC_merged = []

    #an indent due to mirroring of the data for voronoi tess.
    indent = 1e-8

    #width of the supports
    supportWidth = maxLim[0] / 80

    ###############generating of nodes, left horizontal support ###############
    #mech bc
    lineBC = np.array([0,0,-1,-1,-1,-1])

    #defining points of the line
    nodeA = np.array([indent, indent])
    nodeB = np.array([indent + supportWidth, indent])

    oldLen = len(node_coords)
    pointGenerators.generateNodesLine2dRand(nodeA, nodeB, minDist, dim, node_coords, trials, True, False)
    nrOfPoints =  (len(node_coords)) - oldLen
    #print (nrOfPoints)

    #adding mech boundary conditions
    for n in range ( nrOfPoints ):
        mBC = utilitiesMech.mechanicalBC(dim, oldLen + n, lineBC)
        mechBC_merged.append(mBC)

        #print('adding')

    ###############generating of nodes, right horizontal support ###############
    #mech bc
    lineBC = np.array([-1,0,-1,-1,-1,-1])

    #defining points of the line
    nodeA = np.array([maxLim[0] - supportWidth -indent, indent])
    nodeB = np.array([maxLim[0] - indent, indent])

    oldLen = len(node_coords)
    pointGenerators.generateNodesLine2dRand(nodeA, nodeB, minDist, dim, node_coords, trials, True, False)
    nrOfPoints =  (len(node_coords)) - oldLen
    #print (nrOfPoints)

    #adding mech boundary conditions
    for n in range ( nrOfPoints ):
        mBC = utilitiesMech.mechanicalBC(dim, oldLen + n, lineBC)
        mechBC_merged.append(mBC)
        #print('adding')

    ############### loaded top face ###############
    lineBC = np.array([-1,-1,-1,-1, 1,-1])
    #lineIC = np.array([12.1 , 24.2  , 36.3])

    #defining points of the line
    nodeA =  np.array([indent , maxLim[1] - indent])
    nodeB =  np.array([maxLim[0] - indent , maxLim[1] - indent])

    #nodeA =  np.array([maxLim[0]/2 - maxLim[0]/1000  , maxLim[1] - indent])
    #nodeB =  np.array([maxLim[0]/2 + maxLim[0]/1000  , maxLim[1] - indent])

    oldLen = len(node_coords)
    pointGenerators.generateNodesLine2dRand(nodeA, nodeB, minDist, dim, node_coords,  trials, True, True)
    nrOfPoints =  (len(node_coords)) - oldLen
    #print (nrOfPoints)

    #adding mech boundary conditions
    for n in range ( nrOfPoints ):
        mBC = utilitiesMech.mechanicalBC(dim, oldLen + n, lineBC)
        mechBC_merged.append(mBC)
        #print('adding')
        #mIC = utilitiesGeom.mechanicalIC(dim, oldLen + n, lineIC)
        #mechInitC_merged.append(mIC)


    ##########################################generating of points, homogeneous volume
    rectBC = np.array([-1,-1,-1,-1,-1,-1])

    #rect
    oldLen = len(node_coords)
    pointGenerators.generateNodesRect(maxLim, minDist, dim, trials, node_coords)
    #
    newLen = len(node_coords)-1


    return node_coords, mechBC_merged, mechInitC_merged







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




















def assemble3dcylinderUniPress(center, radius, height, minDist, trials, directionDim):
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
    pointGenerators.generateNodesOrtoCilinderSurf3dRand(center, radius-1e-5, height, directionDim, minDist,  node_coords, trials)

    #######################################################################

    ###############generating of points rectangular volume ###############
    pointGenerators.generateNodesOrtoCilinder3dRand(center, radius, height, directionDim, minDist,  node_coords, trials)
    #######################################################################


    return node_coords, mechBC_merged, mechInitC_merged
