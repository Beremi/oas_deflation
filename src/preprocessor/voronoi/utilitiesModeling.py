import numpy as np
import random
import utilitiesGeom
import utilitiesMech
import utilitiesNumeric
import voronoi

def create2dSSBeamUnifLoad(maxLim, minDist, trials ):
    node_coords,node_mechBC, mechBC_merged  = assemble2DSSBeamBending(maxLim, minDist, trials );
    
    print('Conducting Voronoi tesselation...')
    vor = utilitiesNumeric.runMirroredVoronoi (node_coords, 2, maxLim)

    regions, vertices, polygons, areas, centroids, points = voronoi.voronoi_2d(vor, maxLim)

    ### indirect setting of transportBCs by spatial selection of vertices
    transportBC_merged = []
    ### selecting vertices on the left surface
    leftFaceBC = np.array([2,-1])
    boundA = np.array(  [-1e-8 , 0] )
    boundB = np.array(  [ 1e-8 , maxLim[1]]  )
    leftFace = utilitiesGeom.returnSelectedPts(boundA, boundB, vor.vertices)
    #print(leftFace)
    for i in range (len(leftFace)):
        trsBC = utilitiesGeom.transportBC(leftFace[i], leftFaceBC)
        transportBC_merged.append(trsBC)

    ### selecting vertices on the right surface
    rightFaceBC = np.array([3,-1])
    boundA = np.array(  [maxLim[0] - 1e-8, 0] )
    boundB = np.array(  [maxLim[0] + 1e-8 , maxLim[1]]  )
    rightFace = utilitiesGeom.returnSelectedPts(boundA, boundB, vor.vertices)
    #print(rightFace)
    for i in range (len(rightFace)):
        trsBC = utilitiesGeom.transportBC(rightFace[i], rightFaceBC)
        transportBC_merged.append(trsBC)

    return node_coords,node_mechBC, mechBC_merged, transportBC_merged, vor, areas


def create2dCantileverBending(maxLim, minDist, trials ):
    ### sampling of nodes
    ### direct setting of mechanicalBCs
    node_coords,node_mechBC, mechBC_merged  = assemble2DCantileverBending(maxLim, minDist, trials );

    print('Conducting Voronoi tesselation...')
    ### conducting Voronoi tesselation
    vor = utilitiesNumeric.runMirroredVoronoi (node_coords, 2, maxLim)
    ### extracting characteristics of the Vor diagram
    regions, vertices, polygons, areas, centroids, points = voronoi.voronoi_2d(vor, maxLim)

    # fig = voronoi_plot_2d(vor, show_vertices=False, line_colors='orange', line_width=1, line_alpha=1, point_size=0.5)
    # fig.set_figwidth(maxLim[0]*3)
    # fig.set_figheight(maxLim[1]*3)
    #
    # for i in range(points.shape[0]):
    #     plt.text(points[i,0], points[i,1], str(i))
    #     #print (i)
    #
    # plt.xlim(-maxLim[0]*0.2, maxLim[0]*1.2)
    # plt.ylim(-maxLim[1]*0.2, maxLim[1]*1.2)
    # plt.show()

    ### indirect setting of transportBCs by spatial selection of vertices
    transportBC_merged = []
    ### selecting vertices on the left surface
    leftFaceBC = np.array([2,-1])
    boundA = np.array(  [-1e-8 , 0] )
    boundB = np.array(  [ 1e-8 , maxLim[1]]  )
    leftFace = utilitiesGeom.returnSelectedPts(boundA, boundB, vor.vertices)
    #print(leftFace)
    for i in range (len(leftFace)):
        trsBC = utilitiesGeom.transportBC(leftFace[i], leftFaceBC)
        transportBC_merged.append(trsBC)

    ### selecting vertices on the right surface
    rightFaceBC = np.array([3,-1])
    boundA = np.array(  [maxLim[0] - 1e-8, 0] )
    boundB = np.array(  [maxLim[0] + 1e-8 , maxLim[1]]  )
    rightFace = utilitiesGeom.returnSelectedPts(boundA, boundB, vor.vertices)
    #print(rightFace)
    for i in range (len(rightFace)):
        trsBC = utilitiesGeom.transportBC(rightFace[i], rightFaceBC)
        transportBC_merged.append(trsBC)

    return node_coords,node_mechBC, mechBC_merged, transportBC_merged, vor, areas





#
######## METHOD FOR CREATING OF A 2D SUPPORTED CANTILEVER MODEL
def assemble2DCantileverBending (maxLim, minDist, trials):
    dim = 2
    #lists for the model
    node_coords = []
    node_mechBC = []
    mechBC_merged = []

    #an indent due to mirroring of the data for voronoi tess.
    indent = 1e-8

    ###############generating of nodes, supported line left vertical ###############
    #mech bc
    lineBC = np.array([0,0,0,-1,-1,-1])

    #defining points of the line
    nodeA = np.array([indent, indent])
    nodeB = np.array([indent, maxLim[1]-indent])

    oldLen = len(node_coords)
    utilitiesGeom.generateNodesLine2dRand(nodeA, nodeB, minDist, dim, node_coords, node_mechBC, lineBC, trials, True)
    nrOfPoints =  (len(node_coords)) - oldLen
    #print (nrOfPoints)

    #adding mech boundary conditions
    for n in range ( nrOfPoints ):
        mBC = utilitiesGeom.mechanicalBC(dim, oldLen + n, lineBC)
        mechBC_merged.append(mBC)
        #print('adding')

    ###############generating a single point top right (a line of zero length) ###############
    lineBC = np.array([-1,-1,-1,0,1,0])

    #defining points of the line
    nodeA = np.array([maxLim[0] - indent, maxLim[1] - indent])
    nodeB = np.array([maxLim[0] - indent, maxLim[1] - indent])

    oldLen = len(node_coords)
    utilitiesGeom.generateNodesLine2dRand(nodeA, nodeB, minDist, dim, node_coords, node_mechBC, lineBC, trials, False)
    nrOfPoints =  (len(node_coords)) - oldLen
    #print (nrOfPoints)

    #adding mech boundary conditions
    for n in range ( nrOfPoints ):
        mBC = utilitiesGeom.mechanicalBC(dim, oldLen + n, lineBC)
        mechBC_merged.append(mBC)
        #print('adding')

    ##########################################generating of points, homogeneous volume
    rectBC = np.array([-1,-1,-1,-1,-1,-1])

    #rect
    oldLen = len(node_coords)
    utilitiesGeom.generateNodesRect(maxLim, minDist, dim, trials, node_coords, node_mechBC, rectBC)
    #
    newLen = len(node_coords)-1

   # print (nrOfPoints)
  #  mBC = utilitiesGeom.mechanicalBC(dim, kvadrBC, oldLen, newLen)
   # mechBC_merged.append(mBC)
    ####################################################################################################

    return node_coords, node_mechBC, mechBC_merged


def assemble2DSSBeamBending (maxLim, minDist, trials):
    dim = 2
    #lists for the model
    node_coords = []
    node_mechBC = []
    mechBC_merged = []

    #an indent due to mirroring of the data for voronoi tess.
    indent = 1e-8

    #width of the supports
    supportWidth = maxLim[0] / 40

    ###############generating of nodes, left horizontal support ###############
    #mech bc
    lineBC = np.array([0,0,-1,-1,-1,-1])

    #defining points of the line
    nodeA = np.array([indent, indent])
    nodeB = np.array([indent + supportWidth, indent])

    oldLen = len(node_coords)
    utilitiesGeom.generateNodesLine2dRand(nodeA, nodeB, minDist, dim, node_coords, node_mechBC, lineBC, trials, True)
    nrOfPoints =  (len(node_coords)) - oldLen
    #print (nrOfPoints)

    #adding mech boundary conditions
    for n in range ( nrOfPoints ):
        mBC = utilitiesGeom.mechanicalBC(dim, oldLen + n, lineBC)
        mechBC_merged.append(mBC)
        #print('adding')

    ###############generating of nodes, right horizontal support ###############
    #mech bc
    lineBC = np.array([-1,0,0,-1,-1,-1])

    #defining points of the line
    nodeA = np.array([maxLim[0] - supportWidth -indent, indent])
    nodeB = np.array([maxLim[0] - indent, indent])

    oldLen = len(node_coords)
    utilitiesGeom.generateNodesLine2dRand(nodeA, nodeB, minDist, dim, node_coords, node_mechBC, lineBC, trials, True)
    nrOfPoints =  (len(node_coords)) - oldLen
    #print (nrOfPoints)

    #adding mech boundary conditions
    for n in range ( nrOfPoints ):
        mBC = utilitiesGeom.mechanicalBC(dim, oldLen + n, lineBC)
        mechBC_merged.append(mBC)
        #print('adding')

    ############### loaded top face ###############
    lineBC = np.array([-1,-1,-1,-1, 1,-1])

    #defining points of the line
    nodeA =  np.array([indent, maxLim[1] - indent])
    nodeB =  np.array([maxLim[0] - indent , maxLim[1] - indent])

    oldLen = len(node_coords)
    utilitiesGeom.generateNodesLine2dRand(nodeA, nodeB, minDist, dim, node_coords, node_mechBC, lineBC, trials, False)
    nrOfPoints =  (len(node_coords)) - oldLen
    #print (nrOfPoints)

    #adding mech boundary conditions
    for n in range ( nrOfPoints ):
        mBC = utilitiesGeom.mechanicalBC(dim, oldLen + n, lineBC)
        mechBC_merged.append(mBC)
        #print('adding')



    ##########################################generating of points, homogeneous volume
    rectBC = np.array([-1,-1,-1,-1,-1,-1])

    #rect
    oldLen = len(node_coords)
    utilitiesGeom.generateNodesRect(maxLim, minDist, dim, trials, node_coords, node_mechBC, rectBC)
    #
    newLen = len(node_coords)-1

    #print (nrOfPoints)
  #  mBC = utilitiesGeom.mechanicalBC(dim, kvadrBC, oldLen, newLen)
   # mechBC_merged.append(mBC)
    ####################################################################################################

    return node_coords, node_mechBC, mechBC_merged





def create3dCantileverUniTens(maxLim, minDist, trials ):
    node_coords,node_mechBC, mechBC_merged  = assemble2DRectangle(maxLim, minDist, trials )

    areas = voronoi.voronoi_3d(vor, maxLim)

    return node_coords,node_mechBC, mechBC_merged

######## FUNCTION FOR CREATING OF A 3D SUPPORTED RECTANGE
def assemble3Dblock():
    indent = 1e-5

    transtBC = np.array([0,1])
    tsptBC = utilitiesGeom.transportBC(20, transtBC)
    transportBC_merged.append(tsptBC)

    ###############generating of points supported line bottom left ###############
    mechBC = np.array([0,0,0,0,0,0,0,0,0,0,0,0])

    nodeA = np.array([indent , indent, indent])
    nodeB = np.array([indent , maxLim[1] - indent, indent])

    oldLen = len(node_coords)
   # utilitiesGeom.generateNodesLine3dRand(nodeA, nodeB, minDist, dim, node_coords, node_mechBC, mechBC, trials, True)
    #
    nrOfPoints =  (len(node_coords)) - oldLen

    for n in range ( nrOfPoints ):
        mBC = utilitiesGeom.mechanicalBC(dim, oldLen + n, mechBC)
        mechBC_merged.append(mBC)
        #print('adding')

    ###############generating of points supported line top right ###############
    mechBC = np.array([-1,-1,-1, -1,-1,-1,    0,1,0, 0,0,0])

    nodeA = np.array([maxLim[0] - indent , indent, maxLim[2] -indent])
    nodeB = np.array([maxLim[0] - indent , maxLim[1] - indent, maxLim[2] -indent])

    oldLen = len(node_coords)
    utilitiesGeom.generateNodesLine3dRand(nodeA, nodeB, minDist, dim, node_coords, node_mechBC, mechBC, trials, True)
    #
    nrOfPoints =  (len(node_coords)) - oldLen
     #print (nrOfPoints)

    for n in range ( nrOfPoints ):
        mBC = utilitiesGeom.mechanicalBC(dim, oldLen + n, mechBC)
        mechBC_merged.append(mBC)
        #print('adding')


    ###############generating of points supported surface  left face ###############
    mechBC = np.array([0,0,0,0,0,0,-1,-1,-1,-1,-1,-1])

    nodeA = np.array([ indent , indent, indent])
    nodeB = np.array([ indent , maxLim[1] - indent, maxLim[2] -indent])

    oldLen = len(node_coords)
    utilitiesGeom.generateNodesOrtoSurface3dRand(nodeA, nodeB, minDist, dim, node_coords, node_mechBC, mechBC, trials)
    nrOfPoints =  (len(node_coords)) - oldLen
    for n in range ( nrOfPoints ):
        mBC = utilitiesGeom.mechanicalBC(dim, oldLen + n, mechBC)
        mechBC_merged.append(mBC)
        #print('adding')


    ###############generating of points rectangular volume ###############
    mechBC = np.array([-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1])

    #kvadr
    oldLen = len(node_coords)
    utilitiesGeom.generateNodesRect(maxLim, minDist, dim, trials, node_coords, node_mechBC, mechBC)
    newLen = len(node_coords)-1
    nrOfPoints =  (len(node_coords)) - oldLen
   # for n in range ( nrOfPoints ):
       # mBC = utilitiesGeom.mechanicalBC(dim, oldLen + n, mechBC)
       # mechBC_merged.append(mBC)
       # print('adding')
    ####################################################################################################
