import utilitiesModeling
import utilitiesGeom
import numpy as np

def createElasticZone(dim, maxLim, min_dist):
    elaX = min_dist / maxLim[0] * 2

    matZ = []
    if (dim==2):
        boundA = np.array(  [ -1e-8             , -1e-8          ] )
        matZ.append (boundA)
        boundB = np.array(  [ maxLim[0]*elaX    , maxLim[1] + 1e-8] )
        matZ.append (boundB)
        boundA1 = np.array(  [ maxLim[0]-maxLim[0]*elaX , - 1e-8] )
        matZ.append (boundA1)
        boundB1 = np.array(  [ maxLim[0] + 1e-8 , maxLim[1] + 1e8]  )
        matZ.append (boundB1)

    if (dim==3):
        boundA = np.array(  [ -1e-8   -maxLim[0]          , -1e-8    -maxLim[1]         , -1e8 -maxLim[2]] )
        matZ.append (boundA)
        boundB = np.array(  [ maxLim[0]*elaX    , maxLim[1] + 1e8   , maxLim[2] + 1e8  ] )
        matZ.append (boundB)
        boundA1 = np.array(  [ maxLim[0]-maxLim[0]*elaX , - 1e-8  -maxLim[1]      , -1e8-maxLim[2]] )
        matZ.append (boundA1)
        boundB1 = np.array(  [ maxLim[0] + 1e-8 , maxLim[1] + 1e8   , maxLim[2] + 1e8 ]  )
        matZ.append (boundB1)


    return matZ

def prepare_lists ():
    #lists for the model
    node_coords = []
    mechBC_merged = []
    mechIC_merged = []
    trsprtBC_merged = []
    trsprtIC_merged = []
    functions = []
    materials = []
    material_zones = []

    return node_coords, mechBC_merged, mechIC_merged, trsprtBC_merged, trsprtIC_merged, functions, materials, material_zones



def extract_geometry(node_coords, dim, maxLim, vor, areas, material_zones, materials, withoutTransport = False):
     node_coords = np.asarray(node_coords)
     node_count = len(node_coords)

     vert_count, verticesIdxDict, vertIdxStart = utilitiesGeom.extractGeometry(dim, node_count,  maxLim, vor, node_coords, areas, mZ=material_zones, withoutTransport = withoutTransport)

     return vert_count, verticesIdxDict, vertIdxStart

def save_input(materials, functions, dim, mechBC_merged, mechIC_merged, trsprtBC_merged, trsprtIC_merged, verticesIdxDict, vertIdxStart, withoutTransport):
    utilitiesGeom.saveMaterials(materials)
    utilitiesGeom.saveFunctions(functions)
    utilitiesGeom.saveMechBC(dim, mechBC_merged)
    if (len(mechIC_merged)>0):  utilitiesGeom.saveMechIC(dim, mechIC_merged)
    utilitiesGeom.saveTransportBC(trsprtBC_merged, verticesIdxDict, vertIdxStart)
    if (len(trsprtIC_merged)>0):utilitiesGeom.saveTransportIC(trsprtIC_merged)
    utilitiesGeom.saveExporters()

    solStep = 1e-2
    simTime = 100
    utilitiesGeom.saveMasterInput(dim, 0, solStep, 1e-4, 1e-1, simTime, withoutTransport=withoutTransport)



def GenerateModelType1 (trials, xDim, yDim, min_dist, yDispl, mechMaterial, trsprtMaterial):
     #2D: Cantilever bending
     maxLim = np.array([  float(xDim)   ,  float(yDim) ])

     node_coords, mechBC_merged, mechIC_merged, trsprtBC_merged, trsprtIC_merged, functions, materials, material_zones = prepare_lists()

     node_coords, mechBC_merged, mechIC_merged, trsprtBC_merged, trsprtIC_merged, vor, areas, functions = utilitiesModeling.create2dCantileverBending(maxLim, min_dist, trials )

     materials.append(mechMaterial)
     materials.append(trsprtMaterial)

     material_zones.append(createElasticZone(2, maxLim, min_dist))
     material_zones = None

     vert_count, verticesIdxDict, vertIdxStart = extract_geometry(node_coords, 2, maxLim, vor, areas, material_zones, materials)

     save_input(materials, functions, 2, mechBC_merged, mechIC_merged, trsprtBC_merged, trsprtIC_merged, verticesIdxDict, vertIdxStart)




#tube torsion
def GenerateModelType37 (trials, length, outerRadius, thickness, min_dist, rotAng, elaZones, withoutTransport, mechMaterial, mechMaterialEla, trsprtMaterial, rotationAngle):
     length = float(length)
     #tube torsion
     maxLim = np.array([  float(length)   ,  float(outerRadius), float(outerRadius) ])

     node_coords, mechBC_merged, mechIC_merged, trsprtBC_merged, trsprtIC_merged, functions, materials, material_zones = prepare_lists()

     node_coords, mechBC_merged, mechIC_merged, trsprtBC_merged, trsprtIC_merged, vor, areas, functions   = utilitiesModeling.create3dtubeTorsionFree(np.zeros(3), outerRadius, length, thickness, min_dist, trials, int(0), rotationAngle )

     materials.append(mechMaterial)
     materials.append(trsprtMaterial)

     if (elaZones):
         matZ = createElasticZone(3, maxLim, min_dist)
         material_zones.append(matZ)
         #
         materials.append(mechMaterialEla)
     else:
         material_zones = None



     vert_count, verticesIdxDict, vertIdxStart = extract_geometry(node_coords, 3, maxLim, vor, areas, material_zones, materials, withoutTransport=withoutTransport)

     save_input(materials, functions, 3, mechBC_merged, mechIC_merged, trsprtBC_merged, trsprtIC_merged, verticesIdxDict, vertIdxStart, withoutTransport)
