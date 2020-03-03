# uncompyle6 version 3.6.4
# Python bytecode 3.7 (3394)
# Decompiled from: Python 3.7.3 (default, Mar 27 2019, 22:11:17)
# [GCC 7.3.0]
# Embedded file name: /home/jm/GitWorkspace/partmod1/partmod/src/preprocessor/voronoi/Preprocessor.py
# Size of source mod 2**32: 11535 bytes
import Preprocessor, sys, time, numpy as np, random
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import scipy
from IPython.display import clear_output
import sys, os, math
from sklearn import preprocessing
from scipy.ndimage import rotate
from scipy.spatial import Voronoi
from scipy.spatial import voronoi_plot_2d
from scipy.spatial import Delaunay
from scipy.sparse.csgraph import reverse_cuthill_mckee
from scipy.sparse import csr_matrix
from scipy.sparse import csc_matrix
import utilitiesGeom, utilitiesMech, utilitiesModeling, utilitiesNumeric, voronoi


if __name__ == '__main__':
    print('\n%%%%%%%%% LATTICE PREPROCESSOR STARTED %%%%%%%%%')
    start = time.time()

    if len(sys.argv) > 1:
        seed = int(sys.argv[1])
    else:
        seed = np.random.randint(1000.0)
    np.random.seed(seed=seed)

    if len(sys.argv) > 1:
        minDist = float(sys.argv[2])
    else:
        minDist = 0.2

    periodicModel = 0
    nodePositions = []
    coupledNodes = []
    mirtype = []

    #type of solver. does not matter now
    #solver = "SteadyStateNonLinearSolver"
    solver = 'SteadyStateNonLinearSolver'

    #power tesselation on/off  does not matter now
    powerTes = False

    #dimension
    dim = 3

    print('Creating a %dd lattice model...' % dim)
    #coupled problem?
    activeTransport = 0
    activeMechanics = 1

    #dimensions of rectangle model
    Xdim = 2.0
    Ydim = 0.2
    Zdim = 0.2

    #size of grains (minimum distance between nodes)
    #be cautious with small grains!
    minDist = 0.08

    #trials of random node positioning
    trials = 40000

    #dimensions of cylinder model
    cylinderRad = 0.2
    cylinderHeight = 1
    tubeThickness = 0.05

    #dimensions of rectangle model
    if dim == 2:
        maxLim = np.array([Xdim, Ydim])
    if dim == 3:
        maxLim = np.array([Xdim, Ydim, Zdim])

    #volume of the model (later for check)
    volume = np.sum(maxLim)

    radius = minDist / 2
    elaX = minDist / Xdim * 2
    if dim == 2:
        dV = 3.141592 * radius ** 2
    if dim == 3:
        dV = 4.188789333333333 * radius ** 3
    expNodes = volume / dV * 0.5
    print('Expecting about %d nodes' % expNodes)

     #lists for the model
    node_coords = []
    mechBC_merged = []
    mechIC_merged = []
    trsprtBC_merged = []
    trsprtIC_merged = []
    functions = []
    radii = []
    coupledNodes = None
    totalNodeCount = 0

    master_folder = 'power_%.4f_%02d' % (minDist, seed)
    try:
        if not os.path.exists(master_folder):
            os.makedirs(master_folder)
    except:
        print('Please create directory %s! Code Exited.' % master_folder)
        sys.exit()

    notches = None

    #creating the model. Select the prepared models.
    if (dim == 2):

        #patch test for Transport
        #node_coords, mechBC_merged, trsprtBC_merged, vor, areas, functions, radii  = utilitiesModeling.createPatchTestTransport(maxLim, minDist, trials, dim, powerTes)

        #cantilever bending
        #node_coords, mechBC_merged, mechIC_merged, trsprtBC_merged, trsprtIC_merged, vor, areas, functions   = utilitiesModeling.create2dCantileverBending(maxLim, minDist, trials )
        #materialZones=None

        #cantilever  pressure free contraction
        #node_coords, mechBC_merged, mechIC_merged, trsprtBC_merged, trsprtIC_merged, vor, areas, functions   = utilitiesModeling.create2dCantileverUniTens(maxLim, minDist, trials)

        #confined  pressure
        #node_coords, mechBC_merged, mechIC_merged, trsprtBC_merged, trsprtIC_merged, vor, areas, functions  = utilitiesModeling.create2dbeamConfinedPress(maxLim, minDist, trials )

        #simply supported beam, uniform load
        #node_coords, mechBC_merged, mechIC_merged, trsprtBC_merged, trsprtIC_merged, vor, areas, functions   = utilitiesModeling.create2dSSBeamUnifLoad(maxLim, minDist, trials )

        #single spring test
        #node_coords, mechBC_merged, trsprtBC_merged, vor, areas, functions = utilitiesModeling.createSingleSpringTestModel( 2 )

        #diamond test
        #node_coords, mechBC_merged, trsprtBC_merged, vor, areas, functions = utilitiesModeling.createDiamondTestModel(1, 2)

        #periodic shear test
        #node_coords, mechBC_merged, mechIC_merged, trsprtBC_merged, trsprtIC_merged, vor, areas, functions, nodePositions, coupledNodes, mirtype   = utilitiesModeling.create2dPeriodicShear(maxLim, minDist, trials )
        #materialZones=None
        #periodicModel = 1
        #"""

        #simply supported NOTCHED beam, uniform load
        """
        notchH = 0.1 #notch height in percentage of total beam height
        node_coords, mechBC_merged, mechIC_merged, trsprtBC_merged, trsprtIC_merged, vor, areas, functions, notches  = utilitiesModeling.create2dSSBeamUnifLoad(maxLim, minDist, trials, notch=notchH, loadWidth=0.1)
        materialZones=None
        print(notches)
        """

        #2d dogbone
        #node_coords, mechBC_merged, mechIC_merged, trsprtBC_merged, trsprtIC_merged, vor, areas, functions   = utilitiesModeling.create2dDogBone(minDist, trials, D=0.2 )
        #materialZones=None

    if (dim == 3):

        #patch test for Transport
        #node_coords, mechBC_merged, trsprtBC_merged, vor, areas, functions, radii  = utilitiesModeling.createPatchTestTransport(maxLim, minDist, trials, dim, powerTes)
        #materialZones=None

        #cantilever bending
        """
        node_coords, mechBC_merged, mechIC_merged, trsprtBC_merged, trsprtIC_merged, vor, areas, functions   = utilitiesModeling.create3dCantileverBending(maxLim, minDist, trials )
        materialZones=None
        """

        #cantilever uniform pressure, free contraction
        #node_coords, mechBC_merged, mechIC_merged, trsprtBC_merged, trsprtIC_merged, vor, areas, functions   = utilitiesModeling.create3dCantileverUniPressFree(maxLim, minDist, trials )

        #cantilever uniform pressure, confined
        #node_coords, mechBC_merged, mechIC_merged, trsprtBC_merged, trsprtIC_merged, vor, areas, functions   = utilitiesModeling.create3dCantileverUniPressConfined(maxLim, minDist, trials )

        #cylinder uniform pressure free
        #node_coords, mechBC_merged, mechIC_merged, trsprtBC_merged, trsprtIC_merged, vor, areas, functions   = utilitiesModeling.create3dcylinderUniPressFree(np.zeros(3), cylinderRad, cylinderHeight,  minDist, trials, 0 )
        #materialZones=None

        #cylinder uniform pressure confined
        #node_coords, mechBC_merged, mechIC_merged, trsprtBC_merged, trsprtIC_merged, vor, areas, functions   = utilitiesModeling.create3dcylinderUniPressConfined(np.zeros(3), cylinderRad, cylinderHeight,  minDist, trials, 0 )
        #materialZones=None

        #cylinder torsion free
        #node_coords, mechBC_merged, mechIC_merged, trsprtBC_merged, trsprtIC_merged, vor, areas, functions   = utilitiesModeling.create3dcylinderTorsionFree(np.zeros(3), cylinderRad, cylinderHeight,  minDist, trials, 0 )

        #tube torsion free
        #node_coords, mechBC_merged, mechIC_merged, trsprtBC_merged, trsprtIC_merged, vor, areas, functions   = utilitiesModeling.create3dtubeTorsionFree(np.zeros(3), cylinderRad, cylinderHeight, tubeThickness, minDist, trials, 0 )

        #3d dogbone
        """
        D=0.2
        materialZones = utilitiesModeling.assembleMaterialZones (minDist*2.5, 3, model='dogbone',  D=D, thickness=0.1)
        node_coords, mechBC_merged, mechIC_merged, trsprtBC_merged, trsprtIC_merged, vor, areas, functions   = utilitiesModeling.create3dDogBone(minDist, trials, D=D )
        """

        #3d ss 3PB
        notchH = 0.15
        node_coords, mechBC_merged, mechIC_merged, trsprtBC_merged, trsprtIC_merged, vor, areas, functions, notches, govNodes, govNodesMechBC, rigidPlates = utilitiesModeling.create3dSSBeamUnifLoad(maxLim, minDist, trials, notch=notchH, loadWidth=0.1,fracZoneWidth = 0.2)
        materialZones = None

    node_coords = np.asarray(node_coords)
    node_count = len(node_coords)
    print('Model containing %d nodes successfuly generated.' % node_count)
    end = time.time() - start
    print('Model done in %.3f secs.' % end)
    sys.stdout.flush()

    #reordering nodes due to their connectivity
    #order = utilitiesNumeric.reorderToDiagonal(node_count, node_coords, vor)

    materials = []
    young = 30000000000.0
    poisson = 0.3
    density = 2200
    ft = 2000000.0
    Gt = 500
    marsMaterial = utilitiesMech.MarsMaterial(young, poisson, density, ft, Gt)
    materials.append(marsMaterial)
    fatigueMaterial = utilitiesMech.FatigueMaterial(35000000000.0, 0.3, 2200, 200000000.0, 35000000.0, 4000000000.0, 20000000000.0, -2e-07, 0.004, 4000000.0, 0.0, 10000000.0, 250.0, 0)


    transpC = 11
    transpS = 22
    transportMaterial = utilitiesMech.TransportMaterial(transpC, transpS)
    materials.append(transportMaterial)

    linElMaterial = utilitiesMech.linearElasticMaterial(young, poisson, density)
    materials.append(linElMaterial)

    print('')
    #Deconstructing Voronoi diagram and saving the geometry
    vert_count, verticesIdxDict, vertIdxStart, totalNodeCount = utilitiesGeom.extractGeometry(master_folder, dim, node_count, maxLim, vor, node_coords, areas, activeTransport, activeMechanics, mZ=materialZones, periodicModel=periodicModel, nodePositions=nodePositions, coupledNodes=coupledNodes, mirtype=mirtype, notches=notches)


    # saving rest of input
    utilitiesGeom.saveMaterials(master_folder, materials)
    utilitiesGeom.saveFunctions(master_folder, functions)
    if activeMechanics:
        utilitiesGeom.saveMechBC(master_folder, dim, mechBC_merged)
    if activeMechanics:
        if len(mechIC_merged) > 0:
            utilitiesGeom.saveMechIC(master_folder, dim, mechIC_merged)
    if activeTransport:
        utilitiesGeom.saveTransportBC(master_folder, trsprtBC_merged, verticesIdxDict, vertIdxStart)
    if activeTransport:
        if len(trsprtIC_merged) > 0:
            utilitiesGeom.saveTransportIC(master_folder, trsprtIC_merged)

    utilitiesGeom.saveExporters(master_folder, activeTransport, activeMechanics)


    if govNodes != None:
        if rigidPlates != None:
            if govNodesMechBC != None:
                #print(totalNodeCount)
                utilitiesGeom.saveConstraint(master_folder, dim, govNodes, govNodesMechBC, rigidPlates, totalNodeCount, node_coords)
                constraint = True



    if len(radii) > 0:
        if powerTes:
            utilitiesGeom.saveRadii(master_folder, radii)


    solStep = 0.01
    simTime = 1
    utilitiesGeom.saveMasterInput(master_folder, dim, solver, solStep, 0.000001, 0.1, simTime, activeTransport, activeMechanics, periodic=periodicModel, constraint=constraint)




    end = time.time() - end
    print('Saving done in %.3f secs.' % end)
    print('\nThe model contains:')
    print('Mech nodes: %d' % node_count)
    print('Aux nodes: %d' % (vertIdxStart - node_count))
    print('Vertices: %d' % vert_count)
    utilitiesGeom.checkSavedModel(master_folder, dim, activeMechanics, activeTransport)
    end = time.time() - start
    print('\nAll done in %.3f secs.' % end)
    print('%%%%%%%%% LATTICE PREPROCESSOR FINISHED %%%%%%%%%\n')
