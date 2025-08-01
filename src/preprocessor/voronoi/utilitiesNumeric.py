import numpy as np
import matplotlib.pyplot as plt
import scipy
import math
import matplotlib.pylab as plt
from scipy.sparse.csgraph import reverse_cuthill_mckee
from scipy.sparse import csr_matrix
from scipy.sparse import csc_matrix
import voronoi, power
from power_tesselation import PowerTesselation
import matplotlib
#import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

#2d voronoi a teselace
from scipy.spatial import Voronoi
from scipy.spatial import voronoi_plot_2d
from scipy.spatial import Delaunay

SHOW_PLOT = False
AXIS_ASPECT_EQUAL = False  # True may cause error using newer matplotlib versions
##run voronoi, mirrored data
def runMirroredVoronoi (node_coords, dim, maxLim, shifts=0, notch=None):
    print(maxLim)
    vor = Voronoi(voronoi.mirror_dataBeam(node_coords, dim, maxLim, shifts,notch=notch)[:,:dim]) #the last column might be present representing radii

    if (dim == 2):
        regions, vertices, polygons, areas, centroids, points = voronoi.voronoi_2d(vor, maxLim, shifts = shifts)
        return vor, regions, vertices, polygons, areas, centroids, points
    if (dim == 3):
        volumes = voronoi.voronoi_3d(vor, maxLim);
        return vor, volumes

# run voronoi, rebars 2d
def runMirroredVoronoiRebars (data, dim, sizes, rebarDiameter, rebarDepth, rebarCount):
    vor = Voronoi(voronoi.mirror_data_rebars(data, dim, sizes, rebarDiameter, rebarDepth, rebarCount))

    if (dim == 2):
        regions, vertices, polygons, areas, centroids, points = voronoi.voronoi_2d(vor, sizes)
        return vor, regions, vertices, polygons, areas, centroids, points

    if (dim==3):
        volumes = voronoi.voronoi_3d(vor, sizes);
        return vor, volumes


# run voronoi, rebars 2d
def runMirroredVoronoiClover (data, dim, sizes, holeDiameter):
    vor = Voronoi(voronoi.mirror_data_clover(data, dim, sizes, holeDiameter))

    if (dim == 2):
        regions, vertices, polygons, areas, centroids, points = voronoi.voronoi_2d(vor, sizes)
        return vor, regions, vertices, polygons, areas, centroids, points

    if (dim==3):
        volumes = voronoi.voronoi_3d(vor, sizes);
        return vor, volumes

# run voronoi, rebars 2d
def runMirroredVoronoiTDCB (data, dim, sizes, holeDiameter):
    vor = Voronoi(voronoi.mirror_data_TDCB(data, dim, sizes, holeDiameter))

    if (dim == 2):
        regions, vertices, polygons, areas, centroids, points = voronoi.voronoi_2d(vor, sizes)
        return vor, regions, vertices, polygons, areas, centroids, points

    if (dim==3):
        volumes = voronoi.voronoi_3d(vor, sizes);
        return vor, volumes


##run voronoi, mirrored data
def mirror_circ(data, dist, S, D):
    dist = np.expand_dims(dist, 1)
    v = data - S
    v_normed = v / dist
    return (v - v_normed * (dist - 0.725 * D) * 2) + S

def is_inside_dogbone(points, radii, D, SL, SR):
    ''' return boolean array (True == inside) and
    indices array of positions of points inside of dogbone
    '''
    mask = (((points[:, 0] > radii) &
            (points[:, 0] < (D - radii))) &
            ((points[:, 1] > radii) &
            (points[:, 1] < (1.5 * D - radii))))

    mask_tmp = mask.copy()
    mask_tmp[((points[:, 1] < D/4) | (points[:, 1] > (D/4 + D)))] = False
    mask_circ = (
                ((np.sum((points[mask_tmp] - SL)**2, axis=1)**0.5 - 0.725 * D) - radii[mask_tmp] > 0) &
                 ((np.sum((points[mask_tmp] - SR)**2, axis=1)**0.5 - 0.725 * D) - radii[mask_tmp] > 0))
    mask[mask_tmp] = mask_circ
    return mask

def runMirroredVoronoiDogBone (node_coords, dim, D, shifts=0, thickness = None):
    points = node_coords
    bound_lim = D/2
    SL = np.array([-0.725*D+D/5, 0.75*D]) # left circle centeroid
    SR = np.array([D+0.725*D-D/5, 0.75*D]) # right circle centeroid
    # Find vertices outside of dogbone and find its ridge_point
    pt = PowerTesselation(points, None, limits='auto')#Voronoi(points)#PowerTesselation(points, radii, limits='auto')
    mask_outside = ~is_inside_dogbone(pt.vertices, np.zeros(pt.vertices.shape[0]), D, SL, SR)
    outside_idx = np.where(mask_outside)[0]
    out_region = []
    for region in pt.regions[1:]:
        out_region.append(np.any(np.array(region)[:, None] == outside_idx[None, :]))
    points_to_mirror_idx = (pt.point_region[out_region]-1).tolist()
    points_to_mirror = points[points_to_mirror_idx]
    #points_to_mirror = points

    CLB = np.array([0, 0]) # np.array([0, D/4])
    CLT = np.array([0, 2* D/4 + D]) # np.array([0, D/4 + D])
    CRB = np.array([D, 0]) # np.array([D, D/4])
    CRT = np.array([D, 2 * D/4 + D]) # np.array([D, D/4 + D])
    vLB = SL - CLB
    vLT = CLT - SL
    mask1 = ~(np.cross(points_to_mirror-CLB, vLB) < 0)
    mask2 = ~(np.cross(points_to_mirror-CLT, vLT) < 0)
    vRB = SR - CRB
    vRT = CRT - SR
    mask3 = ~(np.cross(points_to_mirror-CRB, vRB) > 0)
    mask4 = ~(np.cross(points_to_mirror-CRT, vRT) > 0)
    mask_circle = mask1 & mask2 & mask3 & mask4
    # print(mask_circle.shape, points_to_mirror.shape)
    points_to_mirror_circle = points_to_mirror[mask_circle]


    mask_top = points_to_mirror[:, 1] > (1.5 * D - bound_lim)
    mask_bottom = points_to_mirror[:, 1] < (bound_lim)
    # mask_left = (points_to_mirror[:, 0] < (bound_lim)) & ((points_to_mirror[:, 1] < (D/4)) | (points_to_mirror[:, 1] > (D/4 + D)))
    # mask_right = (points_to_mirror[:, 0] > (D - bound_lim)) & ((points_to_mirror[:, 1] < (D/4)) | (points_to_mirror[:, 1] > (D/4 + D)))
    mask_left = (points_to_mirror[:, 0] < (bound_lim)) & ((points_to_mirror[:, 1] < (D/2)) | (points_to_mirror[:, 1] > (D)))
    mask_right = (points_to_mirror[:, 0] > (D - bound_lim)) & ((points_to_mirror[:, 1] < (D/2)) | (points_to_mirror[:, 1] > (D)))
    dist2_L = np.sum((points_to_mirror_circle - SL)**2, axis=1)
    dist_L = dist2_L**0.5
    mask_left_circ = ((dist_L < (0.725 * D + bound_lim)) &
                    ((points_to_mirror_circle[:, 1] > (D/8)) &
                    (points_to_mirror_circle[:, 1] < (D/4 + D/8 + D))))
    # mask_left_circ = ((dist_L < (0.725 * D + bound_lim)) &
    #                   ((points_to_mirror_circle[:, 1] > (D/4)) &
    #                    (points_to_mirror_circle[:, 1] < (D/4 + D))))
    dist2_R = np.sum((points_to_mirror_circle - SR)**2, axis=1)
    dist_R = dist2_R**0.5
    mask_right_circ = ((dist_R < (0.725 * D + bound_lim)) &
                    ((points_to_mirror_circle[:, 1] > (D/8)) &
                        (points_to_mirror_circle[:, 1] < (D/4 + D/8 + D))))
    # mask_right_circ = ((dist_R < (0.725 * D + bound_lim)) &
    #                    ((points_to_mirror_circle[:, 1] > (D/4)) &
    #                     (points_to_mirror_circle[:, 1] < (D/4 + D))))


    pointsOut= np.vstack((
    points,
    np.array([0,3*D]) + points_to_mirror[mask_top, :] * np.array([1,-1]), #nahoru
    np.array([0,0]) + points_to_mirror[mask_bottom, :] * np.array([1,-1]), #dolu
    np.array([0,0]) + points_to_mirror[mask_left] * np.array([-1,1]), #doleva
    np.array([D*2,0]) + points_to_mirror[mask_right] * np.array([-1,1]), #doprava
    mirror_circ(points_to_mirror_circle[mask_left_circ], dist_L[mask_left_circ], SL, D),
    mirror_circ(points_to_mirror_circle[mask_right_circ], dist_R[mask_right_circ], SR, D),

    # np.array([0,0]) + points * np.array([-1,-1]), #nahoru doleva
    # np.array([2*D,0]) + points * np.array([-1,-1]), #nahoru doprava
    # np.array([0,6/4*2*D]) + points * np.array([-1,-1]), #dolu doleva
    # np.array([2*D,6/4*2*D]) + points * np.array([-1,-1]), #dolu doprava
    ))

    points = pointsOut

    #vor = Voronoi(voronoi.mirror_dataDogBone(node_coords, dim, D, thickness=thickness)[0])
    vor = PowerTesselation(points, None, limits='auto')#Voronoi(points)
    if SHOW_PLOT:
        fig, ax = plt.subplots()
        voronoi_plot_2d(vor, ax=ax)
        ax.set_aspect('equal')
        plt.show()
    return vor


##run power, mirrored data
def runMirroredPowerDogBone (node_coords, dim, D, shifts=0, thickness = None, radii = []):
    points = node_coords
    bound_lim = D/4
    SL = np.array([-0.725*D+D/5, 0.75*D]) # left circle centeroid
    SR = np.array([D+0.725*D-D/5, 0.75*D]) # right circle centeroid
    # Find vertices outside of dogbone and find its ridge_point
    pt = PowerTesselation(points, radii, limits='auto')
    mask_outside = ~is_inside_dogbone(pt.vertices, np.zeros(pt.vertices.shape[0]), D, SL, SR)
    outside_idx = np.where(mask_outside)[0]
    out_region = []
    for region in pt.regions[1:]:
        out_region.append(np.any(np.array(region)[:, None] == outside_idx[None, :]))
    points_to_mirror_idx = (pt.point_region[out_region]-1).tolist()
    points_to_mirror = points[points_to_mirror_idx]
    #points_to_mirror = points
    radii_to_mirror = radii[points_to_mirror_idx]
    #radii_to_mirror = radii

    CLB = np.array([0, 0]) # np.array([0, D/4])
    CLT = np.array([0, 2* D/4 + D]) # np.array([0, D/4 + D])
    CRB = np.array([D, 0]) # np.array([D, D/4])
    CRT = np.array([D, 2 * D/4 + D]) # np.array([D, D/4 + D])
    vLB = SL - CLB
    vLT = CLT - SL
    mask1 = ~(np.cross(points_to_mirror-CLB, vLB) < 0)
    mask2 = ~(np.cross(points_to_mirror-CLT, vLT) < 0)
    vRB = SR - CRB
    vRT = CRT - SR
    mask3 = ~(np.cross(points_to_mirror-CRB, vRB) > 0)
    mask4 = ~(np.cross(points_to_mirror-CRT, vRT) > 0)
    mask_circle = mask1 & mask2 & mask3 & mask4
    # print(mask_circle.shape, points_to_mirror.shape)
    points_to_mirror_circle = points_to_mirror[mask_circle]
    radii_to_mirror_circle = radii_to_mirror[mask_circle]


    mask_top = points_to_mirror[:, 1] > (1.5 * D - bound_lim)
    mask_bottom = points_to_mirror[:, 1] < (bound_lim)
    mask_left = (points_to_mirror[:, 0] < (bound_lim)) & ((points_to_mirror[:, 1] < (D/4)) | (points_to_mirror[:, 1] > (D/4 + D)))
    mask_right = (points_to_mirror[:, 0] > (D - bound_lim)) & ((points_to_mirror[:, 1] < (D/4)) | (points_to_mirror[:, 1] > (D/4 + D)))
    dist2_L = np.sum((points_to_mirror_circle - SL)**2, axis=1)
    dist_L = dist2_L**0.5
    mask_left_circ = ((dist_L < (0.725 * D + bound_lim)) &
                    ((points_to_mirror_circle[:, 1] > (D/8)) &
                    (points_to_mirror_circle[:, 1] < (D/4 + D/8 + D))))
    # mask_left_circ = ((dist_L < (0.725 * D + bound_lim)) &
    #                   ((points_to_mirror_circle[:, 1] > (D/4)) &
    #                    (points_to_mirror_circle[:, 1] < (D/4 + D))))
    dist2_R = np.sum((points_to_mirror_circle - SR)**2, axis=1)
    dist_R = dist2_R**0.5
    mask_right_circ = ((dist_R < (0.725 * D + bound_lim)) &
                    ((points_to_mirror_circle[:, 1] > (D/8)) &
                        (points_to_mirror_circle[:, 1] < (D/4 + D/8 + D))))
    # mask_right_circ = ((dist_R < (0.725 * D + bound_lim)) &
    #                    ((points_to_mirror_circle[:, 1] > (D/4)) &
    #                     (points_to_mirror_circle[:, 1] < (D/4 + D))))


    pointsOut= np.vstack((
    points,
    np.array([0,3*D]) + points_to_mirror[mask_top, :] * np.array([1,-1]), #nahoru
    np.array([0,0]) + points_to_mirror[mask_bottom, :] * np.array([1,-1]), #dolu
    np.array([0,0]) + points_to_mirror[mask_left] * np.array([-1,1]), #doleva
    np.array([D*2,0]) + points_to_mirror[mask_right] * np.array([-1,1]), #doprava
    mirror_circ(points_to_mirror_circle[mask_left_circ], dist_L[mask_left_circ], SL, D),
    mirror_circ(points_to_mirror_circle[mask_right_circ], dist_R[mask_right_circ], SR, D),

    # np.array([0,0]) + points * np.array([-1,-1]), #nahoru doleva
    # np.array([2*D,0]) + points * np.array([-1,-1]), #nahoru doprava
    # np.array([0,6/4*2*D]) + points * np.array([-1,-1]), #dolu doleva
    # np.array([2*D,6/4*2*D]) + points * np.array([-1,-1]), #dolu doprava
    ))

    if len(radii) > 0:
        #radii = np.tile(radii, 9) #hstack radii 9x
        radii = np.hstack((radii,
                        radii_to_mirror[mask_top],
                        radii_to_mirror[mask_bottom],
                        radii_to_mirror[mask_left],
                        radii_to_mirror[mask_right],
                        radii_to_mirror_circle[mask_left_circ],
                        radii_to_mirror_circle[mask_right_circ]
                        ))

    points = pointsOut

    vor = PowerTesselation(points, radii, limits='auto')
    #
    # # PLOT
    # fig, ax = plt.subplots()
    #
    # for i, p in enumerate(points):
    #     circle = plt.Circle(p, radii[i], color='grey', fill=True)
    #     ax.add_artist(circle)
    # ax.scatter(points.T[0], points.T[1])
    # voronoi_plot_2d(vor, ax=ax, show_points=False, show_vertices=False)
    # for r in vor.ridge_points:
    #     if (r[0] >= 0) and (r[1] >= 0):
    #         p1 = vor.points[r[0]]
    #         p2 = vor.points[r[1]]
    #         ax.plot([p1[0], p2[0]], [p1[1], p2[1]], color='red')
    #
    # ax.set_aspect('equal')
    # plt.show()
    #
    return vor

def runMirroredPowerDogBone_old (node_coords, dim, D, shifts=0, thickness = None, radii = []):
    # vor = Voronoi(voronoi.mirror_dataDogBone(node_coords, dim, D, thickness=thickness)[0])

    node_coords, radii = voronoi.mirror_dataDogBone(node_coords, dim, D, thickness=thickness, radii = radii)

    # vor = Voronoi(node_coords) #the last column might be present representing radii
    vor = PowerTesselation(node_coords, weights=radii, limits='auto') #(points.min(axis=0)-.5).tolist()+(points.max(axis=0)+.5).tolist())

    # PLOT
    # fig, ax = plt.subplots()
    #
    # for i, p in enumerate(node_coords):
    #     circle = plt.Circle(p, radii[i], color='grey', fill=True)
    #     ax.add_artist(circle)
    # plt.scatter(node_coords.T[0], node_coords.T[1])
    # ax.set_aspect('equal')
    # plt.show()

    return vor


##run power, no mirroring data
def runPowerPlain (node_coords, radii, dim, maxLim, Xtop=0, shifts=0):
    vor = PowerTesselation(node_coords, weights=radii, limits='auto') #(points.min(axis=0)-.5).tolist()+(points.max(axis=0)+.5).tolist())

    """
    fig, ax = plt.subplots()
    voronoi_plot_2d(vor, ax=ax)
    ax.scatter(vor.vertices[:, 0], vor.vertices[:, 1], color='r', zorder=100)
    for (x, y), r in zip(points, radii):
        circle = plt.Circle((x, y), r, color='r', fill=False)
        ax.add_artist(circle)
    ax.set_xlim(0, maxLim[0])
    ax.set_ylim(0, maxLim[1])
    if AXIS_ASPECT_EQUAL:
        ax.set_aspect('equal')
    if SHOW_PLOT:
        plt.show()
    """

    if (dim == 2):
        regions, vertices, polygons, areas, centroids, points = voronoi.voronoi_2d(vor, maxLim, shifts = shifts) #power.power_2d(node_coords, radii, dim, maxLim)
        return vor, regions, vertices, polygons, areas, centroids, points
    if (dim == 3):
        volumes = voronoi.voronoi_3d(vor, maxLim);
        return vor, volumes


##run power, mirrored data
def runMirroredPower (node_coords, radii, dim, maxLim, Xtop=0, shifts=0, nomirror = False, notch=None):

    if nomirror == False:
        points, radii = voronoi.mirror_dataBeam(node_coords, dim, maxLim, shifts, weights=radii,notch=notch)#[:,:dim]
    else:
        points = np.copy(node_coords)

    vor = PowerTesselation(points, weights=radii, limits='auto') #(points.min(axis=0)-.5).tolist()+(points.max(axis=0)+.5).tolist())

    """
    fig, ax = plt.subplots()
    voronoi_plot_2d(vor, ax=ax)
    ax.scatter(vor.vertices[:, 0], vor.vertices[:, 1], color='r', zorder=100)
    for (x, y), r in zip(points, radii):
        circle = plt.Circle((x, y), r, color='r', fill=False)
        ax.add_artist(circle)
    ax.set_xlim(0, maxLim[0])
    ax.set_ylim(0, maxLim[1])
    if AXIS_ASPECT_EQUAL:
        ax.set_aspect('equal')
    if SHOW_PLOT:
        plt.show()
    """

    if (dim == 2):
        regions, vertices, polygons, areas, centroids, points = voronoi.voronoi_2d(vor, maxLim, shifts = shifts) #power.power_2d(node_coords, radii, dim, maxLim)
        return vor, regions, vertices, polygons, areas, centroids, points
    if (dim == 3):
        volumes = voronoi.voronoi_3d(vor, maxLim);
        return vor, volumes

def runMirroredPowerDam (node_coords, radii, dim, maxLim, topsize, shifts=0):
    points, radii = voronoi.mirror_dataDam(node_coords, topsize, dim, maxLim, shifts, weights=radii)#[:,:dim]
    vor = PowerTesselation(points, weights=radii, limits='auto') #(points.min(axis=0)-.5).tolist()+(points.max(axis=0)+.5).tolist())
    return vor

def runCylinderMirroredVoronoi  (node_coords, center, radius, height, directionDim, quarter = False, weights=[]):

    if len(weights) == 0:
        vor = Voronoi(voronoi.mirror_dataCylinder(node_coords, center, radius, height, directionDim, quarter = quarter))
    else:
        points, radii = voronoi.mirror_dataCylinder(node_coords, center, radius, height, directionDim, quarter = quarter, weights=weights)
        vor = PowerTesselation(points, weights=radii, limits='auto')

    volumes = voronoi.volumesCylinder3d (vor, center, radius, height, directionDim )
    return vor, volumes

def runTubeMirroredVoronoi  (node_coords, center, radius, height,thickness, directionDim):
    points, radii = voronoi.mirror_dataTube(node_coords, center, radius, height,  thickness, directionDim)
    vor = Voronoi(points)
    volumes = voronoi.volumesCylinder3d (vor, center, radius, height, directionDim )
    return vor, volumes

def runTubeMirroredPower  (node_coords, center, radius, height,thickness, directionDim, radii):
    points, radii = voronoi.mirror_dataTube(node_coords, center, radius, height,  thickness, directionDim, radii)
    vor = PowerTesselation(points, weights=radii, limits='auto')
    volumes = voronoi.volumesCylinder3d (vor, center, radius, height, directionDim )
    return vor, volumes

##################################################
#### General function set by table ####
class generalFunc:
    def __init__(self, table):
        self.table = table

    def getString(self):
        print (self.table)
        line = 'PWLFunction\t%d'%(len(self.table))

        for i in range (len(self.table)):
            line += '\t%f'%(self.table[i][0])

        for i in range (len(self.table)):
            line += '\t%f'%(self.table[i][1])

        return line
####################################################
class PWLFuncFromTxt:
    def __init__(self, filename):
        self.table = np.loadtxt(filename)

    def getString(self):
        line = 'PWLFunction\t%d'%(len(self.table))

        for i in range (len(self.table)):
            line += '\t%f'%(self.table[i][0])

        for i in range (len(self.table)):
            line += '\t%f'%(self.table[i][1])

        return line


##################################################
#### General function set by table ####
class constantFunc:
    def __init__(self, val):
        self.val = val

    def getString(self):
        line = 'PWLFunction\t1'
        line += '\t%f'%(0)
        line += '\t%f'%(self.val)

        return line
####################################################
# displType 0-dx, 1-dy, 2-dz, 3-rx, 4-ry, 5-rz
class constSawToothRotationFunction:
    def __init__ (self, angles, nodeCoords, value, displType, period =None, sym=None, lower=None, time=None, num_cycles=None):
        self.angles = angles.copy()
        self.nodeCoords = nodeCoords.copy()
        self.value = value
        self.period = period
        self.sym = sym
        self.lower = lower
        self.time = time
        self.num_cycles = num_cycles
        self.displType = displType

    def getString (self):
        line = 'ConstSawToothRotationFunction'
        line += '\trotAngles'
        line += '\t%f' %(self.angles[0])
        line += '\t%f' %(self.angles[1])
        line += '\t%f' %(self.angles[2])
        line += '\tinitNodeCrds'
        line += '\t%f' %(self.nodeCoords[0])
        line += '\t%f' %(self.nodeCoords[1])
        line += '\t%f' %(self.nodeCoords[2])
        line += '\tdisplType'
        line += '\t%d' %(self.displType)


        line += '\tvalue\t%f' %(self.value)
        if (self.period!=None):
            line += '\tperiod\t%f' %(self.period)
        if (self.sym!=None):
            line += '\tsym\t%f' %(self.sym)
        if (self.lower!=None):
            line += '\tlower\t%f' %(self.lower)
        if (self.time!=None):
            line += '\ttime\t%f' %(self.time)
        if (self.num_cycles!=None):
            line += '\tnum_cycles\t%f' %(self.num_cycles)

        return line

##################################################
#### Sine function ####
class sineFunc:
    def __init__ (self, value, period, shift = None):
        self.value = value
        self.period = period
        self.shift = shift
    def getString(self):
        line = 'SinusFn'
        line += '\tvalue\t%f' %(self.value)
        line += '\tperiod\t%f' %(self.period)
        if (self.shift!=None):
            line += '\tshift\t%f' %(self.shift)

        return line
##################################################

##################################################
#### Saw tooth constant function ####
class sawToothConstFunc:
    def __init__(self, value, period =None, sym=None, lower=None, time=None, num_cycles=None):
        self.value = value
        self.period = period
        self.sym = sym
        self.lower = lower
        self.time = time
        self.num_cycles = num_cycles
    def getString (self):
        line = 'ConstSawToothFn'
        line += '\tvalue\t%f' %(self.value)
        if (self.period!=None):
            line += '\tperiod\t%f' %(self.period)
        if (self.sym!=None):
            line += '\tsym\t%f' %(self.sym)
        if (self.lower!=None):
            line += '\tlower\t%f' %(self.lower)
        if (self.time!=None):
            line += '\ttime\t%f' %(self.time)
        if (self.num_cycles!=None):
            line += '\tnum_cycles\t%f' %(self.num_cycles)

        return line
##################################################

##################################################
#### Linear saw tooth constant function ####
class linearSawToothFunc:
    def __init__(self, value, period=None, sym=None, lower=None, time=None, num_cycles=None, multiplier = None):
        self.value = value
        self.period = period
        self.sym = sym
        self.lower = lower
        self.time = time
        self.num_cycles = num_cycles
        self.multiplier = mutiplier
    def getString (self):
        line = 'LinSawToothFn'
        line += '\tvalue\t%f' %(self.value)
        if (self.period!=None):
            line += '\tperiod\t%f' %(self.period)
        if (self.sym!=None):
            line += '\tsym\t%f' %(self.sym)
        if (self.lower!=None):
            line += '\tlower\t%f' %(self.lower)
        if (self.time!=None):
            line += '\ttime\t%f' %(self.time)
        if (self.num_cycles!=None):
            line += '\tnum_cycles\t%f' %(self.num_cycles)
        if (self.multiplier!=None):
            line += '\tmultiplier\t%f' %(self.multiplier)
        return line

##################################################

##################################################
#### varying saw tooth function
class varyingSawToothFunction:
     def __init__(self, pwlFn, constSawToothFn ):
         self.constSawToothFn = constSawToothFn
         self.pwlFn = pwlFn
     def getString (self):
         line = 'VaryingSawToothFn\t'
         line += self.pwlFn.getString()[12:]
         line += '\t'
         line += self.constSawToothFn.getString()[16:]
         return line




#reordering of indices
def reorderToDiagonal (node_count, node_coords, vor):
    A = np.zeros( (node_count,node_count) )

    validRidgeIdxs = []
    for i in range (vor.ridge_points.shape[0]):
            pr = False
            if (vor.ridge_points[i][0] < node_count and vor.ridge_points[i][1] < node_count):
                    pr=True
            if (pr):
                #print(vor.ridge_points[i,:])
                validRidgeIdxs.append(i)

    validRidgeIdxs = np.asarray(validRidgeIdxs)

    for i in range (validRidgeIdxs.size):
        pointA = vor.ridge_points[validRidgeIdxs[i]][0]
        pointB = vor.ridge_points[validRidgeIdxs[i]][1]

        coordsA = node_coords[pointA,:]
        coordsB = node_coords[pointB,:]

        dist = np.linalg.norm( coordsB - coordsA )

        A [pointA][pointB] = 1
        A [pointB][pointA] = 1
        #print(dist)


    #print('original connectivity matrix')
    if SHOW_PLOT:
        fig = plt.figure(figsize=(10, 10))

        ax = fig.add_subplot(1,1,1)
        if AXIS_ASPECT_EQUAL:
            ax.set_aspect('equal')
            #plt.imshow(A)
            #plt.colorbar()
        plt.show()

    C = np.zeros( (node_count,node_count) )
    C = csr_matrix(A)
    order = reverse_cuthill_mckee(C, symmetric_mode=True)
    order = np.asarray(order)

    B = A[order][:,order]
    #print(B)


    if SHOW_PLOT:
        #print('reordered connectivity matrix')
        fig = plt.figure(figsize=(10, 10))

        ax = fig.add_subplot(1,1,1)
        if AXIS_ASPECT_EQUAL:
            ax.set_aspect('equal')
            #plt.imshow(B)
            #plt.colorbar()
        plt.show()

    return order
