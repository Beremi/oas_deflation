# distutils: language=c++
cimport cython
import numpy as np
cimport numpy as np
import random
from libc.stdio cimport printf
from libc.stdlib cimport rand, RAND_MAX
from libcpp.vector cimport vector
#from cpprandom cimport mt19937_64, uniform_real_distribution
from libcpp cimport bool
from libc.math cimport sin, cos, abs, fmin, sqrt
import sys

@cython.boundscheck(False)  # Deactivate bounds checking
@cython.wraparound(False)   # Deactivate negative indexing.
@cython.cdivision(True)
def generateNodesRect_cython(double[:] maxLim,
                      double minDist,
                      int dim,
                      int trials,
                      list node_coords,
                      bool useLowBound=False, double topMinDist = -1, double bottomMinDist=-1, int gradienDirection = -1,                        int minDistCenter=-1):
    print('Generating {:d}d block segment of size: {}'.format(dim, ' / '.join('{:f}'.format(i) for i in maxLim)))

    cdef:
        int generatedPoints = 0
        int p, d, i = 0
        int tr = 0
        double distInt
        vector[double] coords
        int node_coords_input_len = 0
        double expectedNodeCount = 0
        vector[double] node_coords_temp
        double[:] topBound, lowBound
        bint distIsGood
        vector[double]  mdC
        #mt19937_64 gen = mt19937_64()
        #uniform_real_distribution[double] dist = uniform_real_distribution[double](0.0, 1.0)
    for d in range(dim):
        coords.push_back(0.0)
        mdC.push_back(0.0)

    minDistGradient = False
    if topMinDist > 0 and bottomMinDist > 0:
        minDistGradient = True

    if node_coords:
        node_coords_input_len = len(node_coords)
        for i in range(node_coords_input_len):
            for d in range(dim):
                node_coords_temp.push_back(node_coords[i][d])


    if dim == 2:
        expectedNodeCount = (maxLim[2]-maxLim[0])*(maxLim[3]-maxLim[1])/(np.pi*(minDist/2)**2) /2
    if dim == 3:
        expectedNodeCount = (maxLim[3]-maxLim[0])*(maxLim[4]-maxLim[1])*(maxLim[4]-maxLim[2])/(4/3*np.pi*(minDist/2)**3) /3


    while (tr < trials):
        tr = 0
        distIsGood = False
        while (not distIsGood) and (tr < trials):
            if (useLowBound==True):
                topBound = maxLim[0:dim]
                lowBound = maxLim[dim:2*dim]
                for d in range(dim):
                    coords[d] = lowBound[d] + np.random.rand() * (topBound[d] - lowBound[d])

            else:
                for d in range(dim):
                    coords[d] = np.random.rand() * maxLim[d]

            distIsGood = True

            for p in range(node_coords_input_len + generatedPoints):
                distInt = 0
                for d in range(dim):
                    dx = node_coords_temp[p*dim+d]
                    dx -= coords[d]
                    distInt += dx * dx
                distInt = distInt**0.5

                currentMinDist = minDist
                if minDistGradient == True:
                    relativePosition = 0
                    if minDistCenter==-1:
                        if useLowBound==True:
                            lowBound = maxLim[dim:2*dim]
                            relativePosition = (coords[gradienDirection]-lowBound[gradienDirection])/(maxLim[gradienDirection]-lowBound[gradienDirection])
                        else:
                            relativePosition = (coords[gradienDirection])/(maxLim[gradienDirection])

                        minDistDiff = topMinDist - bottomMinDist

                        currentMinDist = bottomMinDist + relativePosition * minDistDiff

                    #minDistCenter 1 - top left
                    #minDistCenter 2 - top right
                    #minDistCenter 3 - bot left
                    #minDistCenter 4 - bot right
                    if minDistCenter>0:
                        if minDistCenter == 1 :
                            if useLowBound:
                                mdC[0] = min(min(0, maxLim[0]), maxLim[3])
                                mdC[1] = max(max(0, maxLim[1]), maxLim[4])
                            else:
                                mdC[0] = min(0, maxLim[0])
                                mdC[1] = max(0, maxLim[1])
                        if minDistCenter == 2 :
                            if useLowBound:
                                mdC[0] = max(max(0, maxLim[0]), maxLim[3])
                                mdC[1] = max(max(0, maxLim[1]), maxLim[4])
                            else:
                                mdC[0] = max(0, maxLim[0])
                                mdC[1] = max(0, maxLim[1])
                        if minDistCenter == 3 :
                            if useLowBound:
                                mdC[0] = min(min(0, maxLim[0]), maxLim[3])
                                mdC[1] = min(min(0, maxLim[1]), maxLim[4])
                            else:
                                mdC[0] = min(0, maxLim[0])
                                mdC[1] = min(0, maxLim[1])
                        if minDistCenter == 4 :
                            if useLowBound:
                                mdC[0] = max(max(0, maxLim[0]), maxLim[3])
                                mdC[1] = min(min(0, maxLim[1]), maxLim[4])
                            else:
                                mdC[0] = max(0, maxLim[0])
                                mdC[1] = min(0, maxLim[1])

                          #  print('MDC 0 %f' %mdC[0])
                          #  print('MDC 1 %f' %mdC[1])





                        relativePosition = 0
                        for d in range (2):
                            relativePosition += (mdC[d]-coords[d])*(mdC[d]-coords[d])

                        if not useLowBound:
                            relativePosition = sqrt(relativePosition) / maxLim[0]
                        else:
                            relativePosition = sqrt(relativePosition) / abs(maxLim[0]-maxLim[3])

                        minDistDiff = topMinDist - bottomMinDist
                        currentMinDist = bottomMinDist + relativePosition * minDistDiff

                if (distInt < currentMinDist):
                    distIsGood = False
                    tr += 1
                    break

        # Adding node coords
        if distIsGood and (tr < trials):
            for d in range(dim):
                node_coords_temp.push_back(coords[d])
            generatedPoints += 1
            sys.stdout.write("\033[F") #back to previous line
            sys.stdout.write("\033[K") #clear line
            if gradienDirection ==-1:
                print('Expected %s, generated points %s, trials %s.' %(expectedNodeCount,generatedPoints, tr))
            else:
                print('Generated points %s, trials %s.' %(generatedPoints, tr))

    # Copy back to lists
    for p in range(node_coords_input_len, node_coords_input_len + generatedPoints):
        if (dim==2):
          node_coords.append([node_coords_temp[p*dim], node_coords_temp[p*dim+1]])
        if (dim==3):
          node_coords.append([node_coords_temp[p*dim], node_coords_temp[p*dim+1], node_coords_temp[p*dim+2]])


@cython.boundscheck(False)  # Deactivate bounds checking
@cython.wraparound(False)   # Deactivate negative indexing.
@cython.cdivision(True)
def generateNodesOrtoSurface3dRand_cython(double[:] nodeA,double[:] nodeB,
                      double minDist,
                      int dim,
                      list node_coords,
                      int trials    ,bool minDistAmongNewPoints=False   ):
    print('Generating 3d surface')
    cdef:
        int generatedPoints = 0
        int p, d, i = 0
        int tr = 0
        double distInt
        vector[double] coords
        int node_coords_input_len = 0
        vector[double] node_coords_temp
        bint distIsGood
        vector[double]  mdC
        #mt19937_64 gen = mt19937_64()
        #uniform_real_distribution[double] dist = uniform_real_distribution[double](0.0, 1.0)
    for d in range(dim):
        coords.push_back(0.0)
        mdC.push_back(0.0)

    if node_coords:
        node_coords_input_len = len(node_coords)
        for i in range(node_coords_input_len):
            for d in range(dim):
                node_coords_temp.push_back(node_coords[i][d])

    while (tr < trials):
        tr = 0
        distIsGood = False
        while (not distIsGood) and (tr < trials):
            for c in range (dim):
                if (nodeA[c] == nodeB[c]):
                    coords[c] = nodeA[c]
                else:
                    coords[c] = (nodeB[c] - nodeA[c])*np.random.uniform()  + nodeA[c]

            distIsGood = True

            for p in range(node_coords_input_len + generatedPoints):
                distInt = 0
                for d in range(dim):
                    dx = node_coords_temp[p*dim+d]
                    dx -= coords[d]
                    distInt += dx * dx
                distInt = distInt**0.5

                currentMinDist = minDist

                if (distInt < currentMinDist):
                    distIsGood = False
                    tr += 1
                    break

        # Adding node coords
        if distIsGood and (tr < trials):
            for d in range(dim):
                node_coords_temp.push_back(coords[d])
            generatedPoints += 1
            #sys.stdout.write("\033[F") #back to previous line
            #sys.stdout.write("\033[K") #clear line
            #print('points %s, trials %s' %(generatedPoints, tr))





    # Copy back to lists
    for p in range(node_coords_input_len, node_coords_input_len + generatedPoints):
        if (dim==2):
          node_coords.append([node_coords_temp[p*dim], node_coords_temp[p*dim+1]])
        if (dim==3):
          node_coords.append([node_coords_temp[p*dim], node_coords_temp[p*dim+1], node_coords_temp[p*dim+2]])


@cython.boundscheck(False)  # Deactivate bounds checking
@cython.wraparound(False)   # Deactivate negative indexing.
@cython.cdivision(True)
def generateParticlesRect_cython(np.ndarray[np.float64_t, ndim=1] maxLim, double minDiam, double maxDiam, double volumeRatio, int dim, int trials,
                          np.ndarray[np.float64_t, ndim=2] node_coords, np.ndarray[np.float64_t, ndim=1] radii, bool allow_domain_overlap=False, bool periodic_distance=False,useLowBound=False):
        cdef:
            double gap = 0.1
            vector[double] diam, freq, point
            vector[double] node_coords_temp, delta, dist2, radii_temp
            double Volume = 1
            double saturation = 0
            int iters = 0
            int d = 0
            int di = 0
            bool approved
            double radius, dist, dist_tmp, rnd
            double min_dist2 = 1.0e30
            int node_coords_input_len = node_coords.shape[0]
            int generatedPoints = 0, p = 0, i = 0
            vector[double] topBound, lowBound

        if useLowBound == True:
            topBound = maxLim[0:dim]
            lowBound = maxLim[dim:2*dim]
            Volume = 1
            for d in range(dim):
                Volume*= (topBound[d]-lowBound[d])
        else:
            for d in range(dim):
                Volume *= maxLim[d]
        print(Volume)

        if node_coords_input_len > 0:
            for i in range(node_coords_input_len):
                for d in range(dim):
                    node_coords_temp.push_back(node_coords[i, d])

        if node_coords_input_len > 0:
            for i in range(node_coords_input_len):
                radii_temp.push_back(radii[i])

        for d in range(dim):
            point.push_back(0.0)



        cdef int diam_len = 100
        for i in range(diam_len):
            dist_tmp = (minDiam * 0.5 - maxDiam)/(diam_len-1.)
            diam.push_back(maxDiam + dist_tmp*i )
        #diam = np.linspace(maxDiam, minDiam * 0.5, 30)

        for i in range(diam_len):
            freq.push_back(diam[i] / maxDiam)
        if dim == 2:
            for i in range(diam_len):
                freq[i] = (1.065 * sqrt(freq[i]) - 0.053*pow(freq[i], 4) -
                        0.012 * pow(freq[i], 6) - 0.0045*pow(freq[i], 8) -
                        0.0025 * pow(freq[i], 10))
        elif dim == 3:
            for i in range(diam_len):
                freq[i] = sqrt(freq[i])

        radius = maxDiam / 2.
        while ((2 * radius > minDiam) and (iters < trials)):

            #rnd = np.random.rand(dim)
            if (allow_domain_overlap):
                for d in range(dim):
                    point[d] = random.random() * maxLim[d]
                    #point[d] = rnd[d] * maxLim[d]
                else:
                    for d in range(dim):
                        point[d] = maxLim[3+d] + random.random()*(maxLim[d] - maxLim[3+d])
            else:
                for d in range(dim):
                    #point[d] = rnd[d] * (maxLim[d] - radius * 2) + radius
                    point[d] = random.random() * (maxLim[d] - radius * 2) + radius


            approved = False
            if ( (node_coords_input_len + generatedPoints)==0):
                approved = True
            else:
                dist2.clear()
                delta.clear()
                for i in range(node_coords_input_len + generatedPoints):
                    for d in range(dim):
                        delta.push_back(abs(node_coords_temp[i*dim+d] - point[d]))
                if (periodic_distance):
                    for i in range(node_coords_input_len + generatedPoints):
                        dist = 0
                        dist_tmp = 0
                        for d in range(dim):
                            dist_tmp = fmin(delta[i*dim+d], maxLim[d] - delta[i*dim+d])
                            dist_tmp *= dist_tmp
                            dist += dist_tmp
                        dist2.push_back(dist)
                else:
                    for i in range(node_coords_input_len + generatedPoints):
                        dist = 0
                        for d in range(dim):
                            dist += delta[i*dim+d] * delta[i*dim+d]
                        dist2.push_back(dist)
                #printf('XXX %d - %d\n', node_coords_input_len + generatedPoints, len(delta))
                #sys.stdout.flush()
                min_dist2 = 1e30 # reset min
                for i in range(node_coords_input_len + generatedPoints):
                    dist_tmp = ((1 + gap) * (radii_temp[i] + radius))
                    dist_tmp = dist2[i] - dist_tmp * dist_tmp
                    if min_dist2 > dist_tmp:
                        min_dist2 = dist_tmp
                if min_dist2>0:
                    approved = True
            if ( approved) :
                for d in range(dim):
                    node_coords_temp.push_back(point[d])
                radii_temp.push_back(radius)
                iters = 0
                if (dim==2):  saturation += np.pi * radius*radius / Volume
                elif (dim==3): saturation += np.pi/6 * (radius*2.)* (radius*2.)* (radius*2.) / Volume

                while ((1. - saturation / volumeRatio ) < freq[di]):
                    di+=1
                radius = (diam[di] + (diam[di - 1] - diam[di]) / (freq[di - 1] - freq[di])*(1. - saturation / volumeRatio - freq[di])) / 2.

                generatedPoints += 1
                sys.stdout.write('\r'+'Particles:' +  str(node_coords_input_len + generatedPoints) + ' \t Radius ' +str(radius/2))
                sys.stdout.flush()
            else: iters += 1

        print("Saturation of the volume: ", saturation)

        # Copy back to lists
        cdef np.ndarray node_coords_res = np.zeros([generatedPoints, dim], dtype=np.float64)
        cdef np.ndarray radii_res = np.zeros([generatedPoints], dtype=np.float64)
        node_coords = np.vstack((node_coords, node_coords_res))
        radii = np.hstack((radii, radii_res))
        for p in range(node_coords_input_len, node_coords_input_len + generatedPoints):
            if (dim==2):
                node_coords[p, 0] = node_coords_temp[p*dim]
                node_coords[p, 1] = node_coords_temp[p*dim+1]
            if (dim==3):
                node_coords[p, 0] = node_coords_temp[p*dim]
                node_coords[p, 1] = node_coords_temp[p*dim+1]
                node_coords[p, 2] = node_coords_temp[p*dim+2]
            radii[p] = radii_temp[p]

        return node_coords, radii


@cython.boundscheck(False)  # Deactivate bounds checking
@cython.wraparound(False)   # Deactivate negative indexing.
@cython.cdivision(True)
def generateParticlesOrtoCilinder3dRand_cython(double[:] center,
                                        double radiusCyl,
                                        double height,
                                        int directionDim,
                                        double minDiam,
                                        double maxDiam,
                                        double volumeRatio,
                                        int trials,
                                        np.ndarray[np.float64_t, ndim=2] node_coords,
                                        np.ndarray[np.float64_t, ndim=1] radii):
                                        #np.ndarray[np.float64_t, ndim=2] node_coords, np.ndarray[np.float64_t, ndim=1] radii
    print ('POWER Generating a 3d cylinder segment... (cython)')
    cdef:
        double gap = 0.1
        double PI = np.pi
        double Volume = PI * radiusCyl**2 * height
        vector[double] diam, freq, point
        vector[double] node_coords_temp, delta, dist2, radii_temp
        double saturation = 0
        int iters = 0
        int d = 0
        int di = 0
        bool approved
        double radius, dist, dist_tmp, rnd, angle
        double min_dist2 = 1.0e30
        int node_coords_input_len = node_coords.shape[0]
        int generatedPoints = 0, p = 0, i = 0
        int dim = 3
        double rn

    if node_coords_input_len > 0:
        for i in range(node_coords_input_len):
            for d in range(dim):
                node_coords_temp.push_back(node_coords[i, d])

    if node_coords_input_len > 0:
        for i in range(node_coords_input_len):
            radii_temp.push_back(radii[i])

    cdef int diam_len = 30
    for i in range(diam_len):
        dist_tmp = (minDiam * 0.5 - maxDiam)/(diam_len-1.)
        diam.push_back(maxDiam + dist_tmp*i )
    #diam = np.linspace(maxDiam, minDiam * 0.5, 30)

    for i in range(diam_len):
        freq.push_back(diam[i] / maxDiam)
    if dim == 2:
        for i in range(diam_len):
            freq[i] = (1.065 * sqrt(freq[i]) - 0.053*pow(freq[i], 4) -
                    0.012 * pow(freq[i], 6) - 0.0045*pow(freq[i], 8) -
                    0.0025 * pow(freq[i], 10))
    elif dim == 3:
        for i in range(diam_len):
            freq[i] = sqrt(freq[i])

    for d in range(dim):
        point.push_back(0.0)

    radius = maxDiam / 2.
    while ((2 * radius > minDiam) and (iters < trials)):
        angle = random.random() * PI * 2
        rn = random.random() * radiusCyl
        for d in range(dim):
            point[d] = 0
            # point[d] += center[d]
        if (directionDim == 0 ):
            point[0] = height * random.random()
            point[1] = cos(angle) * rn
            point[2] = sin(angle) * rn
        if (directionDim == 1):
            point[0] = cos(angle) * rn
            point[1] = height * random.random()
            point[2] = sin(angle) * rn
        if (directionDim == 2):
            point[0] = cos(angle) * rn
            point[1] = sin(angle) * rn
            point[2] = height * random.random()
        #for d in range(dim):
        #    point[d] = 0
            # point[d] += center[d]
        #    if d == directionDim:
        #        point[d] = random.random() * height
        #    else:
        #        if directionDim > 0:
        #            point[0] = rn * cos(angle)

        approved = False
        if ( (node_coords_input_len + generatedPoints)==0):
            approved = True
        else:
            dist2.clear()
            delta.clear()
            for i in range(node_coords_input_len + generatedPoints):
                for d in range(dim):
                    delta.push_back(abs(node_coords_temp[i*dim+d] - point[d]))
            for i in range(node_coords_input_len + generatedPoints):
                dist = 0
                for d in range(dim):
                    dist += delta[i*dim+d] * delta[i*dim+d]
                dist2.push_back(dist)
            #printf('XXX %d - %d\n', node_coords_input_len + generatedPoints, len(delta))
            #sys.stdout.flush()
            min_dist2 = 1e30 # reset min
            for i in range(node_coords_input_len + generatedPoints):
                dist_tmp = ((1 + gap) * (radii_temp[i] + radius))
                dist_tmp = dist2[i] - dist_tmp * dist_tmp
                if min_dist2 > dist_tmp:
                    min_dist2 = dist_tmp
            if min_dist2>0:
                approved = True
        if ( approved) :
            for d in range(dim):
                node_coords_temp.push_back(point[d])
            radii_temp.push_back(radius)
            iters = 0
            if (dim==2):  saturation += PI * radius*radius / Volume
            elif (dim==3): saturation += PI/6 * (radius*2.)* (radius*2.)* (radius*2.) / Volume

            while ((1. - saturation / volumeRatio ) < freq[di]):
                di+=1
            radius = (diam[di] + (diam[di - 1] - diam[di]) / (freq[di - 1] - freq[di])*(1. - saturation / volumeRatio - freq[di])) / 2.
            generatedPoints += 1
            sys.stdout.write('\r'+'Particles:' +  str(node_coords_input_len + generatedPoints))
            sys.stdout.flush()
        else: iters += 1

    print("Saturation of the volume: ", saturation)

    # Copy back to lists
    cdef np.ndarray node_coords_res = np.zeros([generatedPoints, dim], dtype=np.float64)
    cdef np.ndarray radii_res = np.zeros([generatedPoints], dtype=np.float64)
    node_coords = np.vstack((node_coords, node_coords_res))
    radii = np.hstack((radii, radii_res))
    for p in range(node_coords_input_len, node_coords_input_len + generatedPoints):
        if (dim==2):
            node_coords[p, 0] = node_coords_temp[p*dim]
            node_coords[p, 1] = node_coords_temp[p*dim+1]
        if (dim==3):
            node_coords[p, 0] = node_coords_temp[p*dim]
            node_coords[p, 1] = node_coords_temp[p*dim+1]
            node_coords[p, 2] = node_coords_temp[p*dim+2]
        radii[p] = radii_temp[p]

    return node_coords, radii



@cython.boundscheck(False)  # Deactivate bounds checking
@cython.wraparound(False)   # Deactivate negative indexing.
@cython.cdivision(True)
def generateNodesOrtoCilinder3dRand_cython(
                      double[:] center,
                      double radius,
                      double height,
                      int directionDim,
                      double minDist,
                      list node_coords,
                      int trials
                      ):
    print ('Generating a 3d cylinder segment. ')
    cdef:
        int generatedPoints = 0
        int p, d
        int tr = 0
        int dim = 3
        double distInt
        double angle
        vector[double] coords
        int node_coords_input_len = len(node_coords)
        vector[double] node_coords_temp
        double[:] rn
        bint distIsGood
        #mt19937_64 gen = mt19937_64()
        #uniform_real_distribution[double] dist = uniform_real_distribution[double](0.0, 1.0)
    for d in range(3):
        coords.push_back(0.0)

    if node_coords:
        for i in range(node_coords_input_len):
            for d in range(dim):
                node_coords_temp.push_back(node_coords[i][d])

    while (tr < trials):
        tr = 0
        distIsGood = False
        while (not distIsGood) and (tr < trials):
            rn = np.random.rand(5)
            angle = rn[0]  * np.pi * 2
            if (directionDim == 0 ):
                coords[0] = center[0] + height * rn[1]
                coords[1] = center[1] +radius * cos(angle) * rn[4]
                coords[2] = center[2] +radius * sin(angle) * rn[4]
            if (directionDim == 1):
                coords[0] = center[0] +radius * cos(angle) * rn[4]
                coords[1] = center[1] +height * rn[2]
                coords[2] = center[2] +radius * sin(angle) * rn[4]
            if (directionDim == 2):
                coords[0] = center[0] +radius * cos(angle) * rn[4]
                coords[1] = center[1] +radius * sin(angle) * rn[4]
                coords[2] = center[2] +height * rn[3]



            distIsGood = True

            for p in range(node_coords_input_len + generatedPoints):
                distInt = 0
                for d in range(3):
                    dx = node_coords_temp[p*dim+d]
                    dx -= coords[d]
                    distInt += dx * dx
                distInt = distInt**0.5
                if (distInt < minDist):
                    distIsGood = False
                    tr += 1
                    break

        # Adding node coords
        if distIsGood and (tr < trials):
            for d in range(3):
                node_coords_temp.push_back(coords[d])
            generatedPoints += 1

    # Copy back to lists
    for p in range(node_coords_input_len, node_coords_input_len + generatedPoints):
        if (dim==2):
          node_coords.append([node_coords_temp[p*dim], node_coords_temp[p*dim+1]])
        if (dim==3):
          node_coords.append([node_coords_temp[p*dim], node_coords_temp[p*dim+1], node_coords_temp[p*dim+2]])





@cython.boundscheck(False)  # Deactivate bounds checking
@cython.wraparound(False)   # Deactivate negative indexing.
@cython.cdivision(True)
def generateNodesOrtoTube3dRand_cython(
                      double[:] center,
                      double outerRad,
                      double height,
                      double thickness,
                      int directionDim,
                      double minDist,
                      list node_coords,
                      int trials
                      ):
    print ('Generating a 3d tube segment. cython ')
    cdef:
        int generatedPoints = 0
        int p, d
        int tr = 0
        int dim = 3
        double distInt
        double angle
        vector[double] coords
        int node_coords_input_len = len(node_coords)
        vector[double] node_coords_temp
        bint distIsGood
        #mt19937_64 gen = mt19937_64()
        #uniform_real_distribution[double] dist = uniform_real_distribution[double](0.0, 1.0)
    for d in range(3):
        coords.push_back(0.0)

    if node_coords:
        for node in node_coords:
            for d in range(3):
                node_coords_temp.push_back(node[d])

    innerRad = outerRad - thickness
    while (tr < trials):
        tr = 0
        distIsGood = False
        while (not distIsGood) and (tr < trials):
            angle = random.random()  * np.pi * 2
            rn = random.random()
            if (directionDim == 0 ):
                coords[0] = center[0] + height * random.random()
                coords[1] = center[1] + (innerRad + thickness * rn) * np.cos(angle)
                coords[2] = center[2] + (innerRad + thickness * rn) * np.sin(angle)
            if (directionDim == 1 ):
                coords[0] = center[0] + (innerRad + thickness * rn) * np.cos(angle)
                coords[1] = center[1] + height * random.random()
                coords[2] = center[2] + (innerRad + thickness * rn) * np.sin(angle)
            if (directionDim == 2 ):
                coords[0] = center[0] + (outerRad + thickness * rn) * np.cos(angle)
                coords[1] = center[1] + (outerRad + thickness * rn) * np.sin(angle)
                coords[2] = center[2] + height * random.random()



            distIsGood = True

            for p in range(node_coords_input_len + generatedPoints):
                distInt = 0
                for d in range(3):
                    dx = node_coords_temp[p*dim+d]
                    dx -= coords[d]
                    distInt += dx * dx
                distInt = distInt**0.5
                if (distInt < minDist):
                    distIsGood = False
                    tr += 1
                    break

        # Adding node coords
        if distIsGood and (tr < trials):
            for d in range(3):
                node_coords_temp.push_back(coords[d])
            generatedPoints += 1

    # Copy back to lists
    for p in range(node_coords_input_len, node_coords_input_len + generatedPoints):
        if (dim==2):
          node_coords.append([node_coords_temp[p*dim], node_coords_temp[p*dim+1]])
        if (dim==3):
          node_coords.append([node_coords_temp[p*dim], node_coords_temp[p*dim+1], node_coords_temp[p*dim+2]])



@cython.boundscheck(False)  # Deactivate bounds checking
@cython.wraparound(False)   # Deactivate negative indexing.
@cython.cdivision(True)
def generateParticlesOrtoTube3dRand_cython(
                      double[:] center,
                      double outerRad,
                      double height,
                      double thickness,
                      int directionDim,
                      double minDiam,
                      double maxDiam,
                      double volumeRatio,
                      np.ndarray[np.float64_t, ndim=2] node_coords,
                      np.ndarray[np.float64_t, ndim=1] radii,
                      int trials
                      ):

    print ('POWER Generating a 3d tube segment. cython ')
    cdef:
        double gap = 0.1
        double PI = np.pi
        double Volume = PI * (outerRad**2 - (outerRad - thickness)**2) * height
        vector[double] diam, freq, point
        vector[double] node_coords_temp, delta, dist2, radii_temp
        double saturation = 0
        int iters = 0
        int d = 0
        int di = 0
        bool approved
        double radius, dist, dist_tmp, rnd, angle
        double min_dist2 = 1.0e30
        int node_coords_input_len = node_coords.shape[0]
        int generatedPoints = 0, p = 0, i = 0
        int dim = 3
        double rn

    if node_coords_input_len > 0:
        for i in range(node_coords_input_len):
            for d in range(dim):
                node_coords_temp.push_back(node_coords[i, d])

    if node_coords_input_len > 0:
        for i in range(node_coords_input_len):
            radii_temp.push_back(radii[i])

    cdef int diam_len = 30
    for i in range(diam_len):
        dist_tmp = (minDiam * 0.5 - maxDiam)/(diam_len-1.)
        diam.push_back(maxDiam + dist_tmp*i )
    #diam = np.linspace(maxDiam, minDiam * 0.5, 30)

    for i in range(diam_len):
        freq.push_back(diam[i] / maxDiam)
    if dim == 2:
        for i in range(diam_len):
            freq[i] = (1.065 * sqrt(freq[i]) - 0.053*pow(freq[i], 4) -
                    0.012 * pow(freq[i], 6) - 0.0045*pow(freq[i], 8) -
                    0.0025 * pow(freq[i], 10))
    elif dim == 3:
        for i in range(diam_len):
            freq[i] = sqrt(freq[i])

    for d in range(dim):
        point.push_back(0.0)

    radius = maxDiam / 2.
    while ((2 * radius > minDiam) and (iters < trials)):
        angle = random.random() * PI * 2
        rn = outerRad - radius - random.random() * (thickness - 2*radius)
        for d in range(dim):
            point[d] = 0
            # point[d] += center[d]
        if (directionDim == 0 ):
            point[0] = (height-2*radius) * random.random() + radius
            point[1] = cos(angle) * rn
            point[2] = sin(angle) * rn
        if (directionDim == 1):
            point[0] = cos(angle) * rn
            point[1] = (height-2*radius)  * random.random() + radius
            point[2] = sin(angle) * rn
        if (directionDim == 2):
            point[0] = cos(angle) * rn
            point[1] = sin(angle) * rn
            point[2] = (height-2*radius)  * random.random() + radius
        #for d in range(dim):
        #    point[d] = 0
            # point[d] += center[d]
        #    if d == directionDim:
        #        point[d] = random.random() * height
        #    else:
        #        if directionDim > 0:
        #            point[0] = rn * cos(angle)

        approved = False
        if ( (node_coords_input_len + generatedPoints)==0):
            approved = True
        else:
            dist2.clear()
            delta.clear()
            for i in range(node_coords_input_len + generatedPoints):
                for d in range(dim):
                    delta.push_back(abs(node_coords_temp[i*dim+d] - point[d]))
            for i in range(node_coords_input_len + generatedPoints):
                dist = 0
                for d in range(dim):
                    dist += delta[i*dim+d] * delta[i*dim+d]
                dist2.push_back(dist)
            #printf('XXX %d - %d\n', node_coords_input_len + generatedPoints, len(delta))
            #sys.stdout.flush()
            min_dist2 = 1e30 # reset min
            for i in range(node_coords_input_len + generatedPoints):
                dist_tmp = ((1 + gap) * (radii_temp[i] + radius))
                dist_tmp = dist2[i] - dist_tmp * dist_tmp
                if min_dist2 > dist_tmp:
                    min_dist2 = dist_tmp
            if min_dist2>0:
                approved = True
        if ( approved) :
            for d in range(dim):
                node_coords_temp.push_back(point[d])
            radii_temp.push_back(radius)
            iters = 0
            if (dim==2):  saturation += PI * radius*radius / Volume
            elif (dim==3): saturation += PI/6 * (radius*2.)* (radius*2.)* (radius*2.) / Volume

            while ((1. - saturation / volumeRatio ) < freq[di]):
                di+=1
            radius = (diam[di] + (diam[di - 1] - diam[di]) / (freq[di - 1] - freq[di])*(1. - saturation / volumeRatio - freq[di])) / 2.
            generatedPoints += 1
            sys.stdout.write('\r'+'Particles:' +  str(node_coords_input_len + generatedPoints))
            sys.stdout.flush()
        else: iters += 1

    print("Saturation of the volume: ", saturation)

    # Copy back to lists
    cdef np.ndarray node_coords_res = np.zeros([generatedPoints, dim], dtype=np.float64)
    cdef np.ndarray radii_res = np.zeros([generatedPoints], dtype=np.float64)
    node_coords = np.vstack((node_coords, node_coords_res))
    radii = np.hstack((radii, radii_res))
    for p in range(node_coords_input_len, node_coords_input_len + generatedPoints):
        if (dim==2):
            node_coords[p, 0] = node_coords_temp[p*dim]
            node_coords[p, 1] = node_coords_temp[p*dim+1]
        if (dim==3):
            node_coords[p, 0] = node_coords_temp[p*dim]
            node_coords[p, 1] = node_coords_temp[p*dim+1]
            node_coords[p, 2] = node_coords_temp[p*dim+2]
        radii[p] = radii_temp[p]

    return node_coords, radii


@cython.boundscheck(False)  # Deactivate bounds checking
@cython.wraparound(False)   # Deactivate negative indexing.
@cython.cdivision(True)
def generateNodesOrtoCylinderSurf3dRand_cython(
                      double[:] center,
                      double radius,
                      double height,
                      int directionDim,
                      double minDist,
                      list node_coords,
                      int trials,
                       angleLimitA=None, angleLimitB=None, mirrorIndent=None, equiAngNodes=0
                      ):
    print ('Generating a 3d cylinder surf cython. ')
    cdef:#
        int generatedPoints = 0
        int p, d
        int tr = 0
        int dim = 3
        double distInt
        double angle
        vector[double] coords
        int node_coords_input_len = len(node_coords)
        vector[double] node_coords_temp
        bint distIsGood
        #mt19937_64 gen = mt19937_64()
        #uniform_real_distribution[double] dist = uniform_real_distribution[double](0.0, 1.0)
    for d in range(3):
        coords.push_back(0.0)

    if node_coords:
        for node in node_coords:
            for d in range(3):
                node_coords_temp.push_back(node[d])

    if equiAngNodes >0:

      for i in range (equiAngNodes):
        if (directionDim==0):
          coords[1] = center[1] + radius * np.cos(2*np.pi / equiAngNodes  *i)
          coords[2] = center[2] + radius * np.sin(2*np.pi / equiAngNodes  *i)
          for z in range (int(height/minDist)):
            coords[0] =  center[0] + (  height / int(height/minDist) ) * (z+0.5)
            for d in range(3):
                node_coords_temp.push_back(coords[d])
            generatedPoints += 1

        if (directionDim==2):
          coords[0] = center[0] + radius * np.cos(2*np.pi / equiAngNodes  *i)
          coords[1] = center[1] + radius * np.sin(2*np.pi / equiAngNodes  *i)
          for z in range (int(height/minDist)):
            coords[2] =  center[2] + (  height / int(height/minDist) ) * (z+0.5)
            for d in range(3):
                node_coords_temp.push_back(coords[d])
            generatedPoints += 1


    if equiAngNodes <=0:
      while (tr < trials):
          tr = 0
          distIsGood = False
          while (not distIsGood) and (tr < trials):
              angle = random.random()  * np.pi * 2
              rn = random.random()
              if (directionDim == 0 ):
                  coords[0] = center[0] + height * random.random()
                  coords[1] = center[1] + radius * np.cos(angle)
                  coords[2] = center[2] + radius * np.sin(angle)
              if (directionDim == 1 ):
                  coords[0] = center[0] + radius * np.cos(angle)
                  coords[1] = center[1] + height * random.random()
                  coords[2] = center[2] + radius * np.sin(angle)
              if (directionDim == 2 ):
                  coords[0] = center[0] + radius * np.cos(angle)
                  coords[1] = center[1] + radius * np.sin(angle)
                  coords[2] = center[2] + height * random.random()



              distIsGood = True

              for p in range(node_coords_input_len + generatedPoints):
                  distInt = 0
                  for d in range(3):
                      dx = node_coords_temp[p*dim+d]
                      dx -= coords[d]
                      distInt += dx * dx
                  distInt = distInt**0.5
                  if (distInt < minDist):
                      distIsGood = False
                      tr += 1
                      break

          # Adding node coords
          if distIsGood and (tr < trials):
              for d in range(3):
                  node_coords_temp.push_back(coords[d])
              generatedPoints += 1

    # Copy back to lists
    for p in range(node_coords_input_len, node_coords_input_len + generatedPoints):
        if (dim==2):
          node_coords.append([node_coords_temp[p*dim], node_coords_temp[p*dim+1]])
        if (dim==3):
          node_coords.append([node_coords_temp[p*dim], node_coords_temp[p*dim+1], node_coords_temp[p*dim+2]])



@cython.boundscheck(False)  # Deactivate bounds checking
@cython.wraparound(False)   # Deactivate negative indexing.
@cython.cdivision(True)
def checkMutDistancesLoops_cython(int dim,
                                  double minDist,
                                  list currentNodes,
                                  list newNode):
    cdef:
        bool distIsGood = True
        int p, d
        int currentNodes_len = len(currentNodes)
        double distInt, dx

    # NOTE JK: currentNodes list must be faltten, but maybe in some older versions it was already sent here as a flatt array (but then does not make sense (len() and range(p*dim + d)))
    if currentNodes_len == 0:
        return distIsGood
    #if type(currentNodes[0]) is not float:
    #  currentNodes = np.concatenate(currentNodes).ravel().tolist()

    minDist *= minDist
    for p in range (currentNodes_len):
        distInt = 0
        for d in range(dim):
            #dx = currentNodes[p*dim+d] - newNode[d]
            dx = currentNodes[p][d] - newNode[d]
            distInt += dx * dx
        #distInt = distInt**0.5
        if (distInt < minDist):
            distIsGood = False
            break


    return distIsGood






@cython.boundscheck(False)  # Deactivate bounds checking
@cython.wraparound(False)   # Deactivate negative indexing.
@cython.cdivision(True)
def generateNodesRectPeriodic_cython(double[:] maxLim,
                      double minDist,
                      int dim,
                      int trials,
                      list node_coords):
    print('Generating {:d}d block segment of size: {}'.format(dim, ' / '.join('{:f}'.format(i) for i in maxLim)))
    cdef:
        int generatedPoints = 0
        int p, d
        int tr = 0
        double distInt, dx
        vector[double] coords
        int node_coords_input_len = len(node_coords)
        vector[double] node_coords_temp
        bint distIsGood
        #mt19937_64 gen = mt19937_64()
        #uniform_real_distribution[double] dist = uniform_real_distribution[double](0.0, 1.0)
    for d in range(dim):
        coords.push_back(0.0)

    if node_coords:
        for node in node_coords:
            for d in range(dim):
                node_coords_temp.push_back(node[d])


    while (tr < trials):
        tr = 0
        distIsGood = False
        while (not distIsGood) and (tr < trials):
            #print(tr)
            #for d in range(dim):
            #    coords[d] = dist(gen) * maxLim[d]
            for d in range(dim):
                 coords[d] = np.random.rand() * maxLim[d]
            distIsGood = True

            for p in range(node_coords_input_len + generatedPoints):
                distInt = 0
                for d in range(dim):
                    dx = node_coords_temp[p*dim+d]
                    dx -= coords[d]

                    if (abs(dx+maxLim[d]) < abs(dx) ):
                        dx += maxLim[d]

                    elif (abs(dx-maxLim[d]) < abs(dx) ):
                        dx -= maxLim[d]


                    distInt += dx * dx
                distInt = distInt**0.5

                if (distInt < minDist):
                    #print(distInt)
                    distIsGood = False
                    tr += 1
                    break

        # Adding node coords
        if distIsGood and (tr < trials):
            for d in range(dim):
                node_coords_temp.push_back(coords[d])
            generatedPoints += 1


    # Copy back to lists
    for p in range(node_coords_input_len, node_coords_input_len + generatedPoints):
        if (dim==2):
          node_coords.append([node_coords_temp[p*dim], node_coords_temp[p*dim+1]])
        if (dim==3):
          node_coords.append([node_coords_temp[p*dim], node_coords_temp[p*dim+1], node_coords_temp[p*dim+2]])






@cython.boundscheck(False)  # Deactivate bounds checking
@cython.wraparound(False)   # Deactivate negative indexing.
@cython.cdivision(True)
def checkMutDistancesLoopsPeriodic_cython(int dim, double minDist, list currentNodes, list newNode, double[:] maxLim):
    cdef:
        bool distIsGood = True
        int p, d
        int currentNodes_len = len(currentNodes)
        double distInt, dx

    for p in range (currentNodes_len):
        distInt = 0
        for d in range(dim):
            dx = currentNodes[p*dim+d] - newNode[d]

            if (abs(dx+maxLim[d]) < abs(dx) ):
                dx += maxLim[d]
            elif (abs(dx-maxLim[d]) < abs(dx) ):
                dx -= maxLim[d]

            distInt += dx * dx

        distInt = distInt**0.5

        if (distInt < minDist):
            distIsGood = False
            break

    return distIsGood
