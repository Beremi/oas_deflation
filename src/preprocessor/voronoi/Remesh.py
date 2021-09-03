import matplotlib as mpl
mpl.use('Agg')
import errno
import os
from ModelGenerator import *
# from regions import *
from pointGenerators import generateNodesRemesh
from regions import *
import numpy as np
import matplotlib.pyplot as plt
import scipy.spatial as ss
import utilitiesGeom

SHOW_PLOT = False

def loadNodes(filename, dim):
    print("loading %s" % (filename))
    if not os.path.isfile(filename):
        raise FileNotFoundError(
                errno.ENOENT, os.strerror(errno.ENOENT), filename)
    try:
        node_coords = np.genfromtxt(filename, skip_header=1,
                                    usecols=( ( dim == 3 ) and (1, 2, 3) or (1, 2) )
                                    )
        if node_coords.size == 0:
            return list()
    except IOError as e:
        print(e, end=' ... probably no entry in the file %s\n'  % (filename))
        return list()

    try:
        len(node_coords[1])
    except Exception as e:
        print(e, end=', probably only one entry in the file %s\n'  % (filename))
        return [ node_coords,  ]
    if (len(node_coords[0]) < dim ):
        return list()
    return [np.array(coor) for coor in node_coords.tolist()]


def loadRegionsToSkip(region_file):
    rtSkip = list()
    if os.path.isfile(region_file):
        with open(region_file) as file:
            lines = file.readlines()

            for line in lines:
                if line[0] == "#":
                    continue
                # print(line)
                ln = line.split()
                if line.startswith("block"):
                    obj = Block(Point(float(ln[1]), float(ln[2]), float(ln[3])),
                                Point(float(ln[4]), float(ln[5]), float(ln[6]))
                                )

                elif line.startswith("sphere"):
                    obj = Sphere(Point(float(ln[1]), float(ln[2]), float(ln[3])),
                                 float(ln[4])
                                 )
                else:
                    continue
                rtSkip.append(obj)
    else:
        print("no such file: %s" % region_file)
    return rtSkip


if __name__ == '__main__':
    print('\n%%%%%%%%% LATTICE REMESH STARTED %%%%%%%%%')
    start = time.time()

    # print("-----------------------------------------------------")
    # print(len(sys.argv))
    # print("-----------------------------------------------------")

    if len(sys.argv) < 9:
        print("not enough information provided for Remesher")
        sys.exit(1)

    prep_input_file = os.path.join(sys.argv[3], sys.argv[1])
    # print(prep_input_file)

    remeshDir = os.path.join(sys.argv[3], sys.argv[2])
    # print(remeshDir)
    # sys.exit()

    if not os.path.isfile(prep_input_file):
        print('No such input file: \'%s\' Exiting.' % prep_input_file)
        sys.exit()

    model = None
    solver = None
    exporters=[]

    # skipLines = ['#', ' ', '\n', '\t']

    f = open (prep_input_file, 'r')
    for row in f:
        if row and row.strip() and not row.startswith('#'):
        # if not (row[0] in skipLines):
            r = row.split()

            if (r[0]=='Model'):
                model = Model(row)
            if (r[0]=='Solver'):
                solver = Solver(row)
            if (r[0]=='Material'):
                if model != None:
                    model.addMaterial(row)
            if (r[0]=='Exporter'):
                exporters.append(r[1:])
            if (r[0]=='SpecifiedNode'):
                node = np.zeros (len(r)-1)
                for c in range (len(r)-1):
                    node[c] = float (r[1+c])
                model.specifiedNodes.append(node)
    f.close()

    if model == None:
        print ('Missing model!! Exiting...')
        sys.exit()
    if solver == None:
        print ('Missing solver!! Exiting...')
        sys.exit()
    if len(model.materials)==0:
        print ('Missing some material!! Exiting...')
        sys.exit()

    ##########################################################################
    ##########################################################################
    # oldDir = "/home/jose/Soft/ParticleModel/TESTS/adaptivity_pokus/TPB_no_notch"
    # node_coords_old = loadNodes(os.path.join(oldDir, "nodes.inp"), model.dimension)

    node_coords_ini = loadNodes(os.path.join(remeshDir, "nodes.out"), model.dimension)

    ctr = loadNodes(os.path.join(remeshDir, "centersToRemesh.out"), model.dimension)
    cpr = loadNodes(os.path.join(remeshDir, "centersFine.out"), model.dimension)
    # print("--------------------------------------------------")
    # print(ctr)
    # print(cpr)
    # print("--------------------------------------------------")
    # where to load them from?
    # should be already in input (specified in master prep file)
    radiusRemesh = float(sys.argv[4])
    radiusTransitional = float(sys.argv[5])
    useExistingFineNodes = bool(int(sys.argv[6]))

    if useExistingFineNodes:
        # print("loading existing fine geoemtry <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<")
        # print("len nodes before: %d" % (len(node_coords_ini)))
        node_coords_ini.extend(loadNodes(os.path.join(remeshDir, "nodesFine.out"), model.dimension))
        # print("len nodes after: %d" % (len(node_coords_ini)))

    remesherSeed = int(sys.argv[7])
    # print(node_coords_ini)

    rtSkip = loadRegionsToSkip(str(sys.argv[8]))

    if len(sys.argv) > 9:
        minDistRemesh = float(sys.argv[9])
    else:
        minDistRemesh = model.minDist / 3.

    # print("radiusRemesh %lg \nradiusTransitional %lg \nminDistRemesh %lg" % (
    #     radiusRemesh, radiusTransitional, minDistRemesh))

    # print("--- before ------------------------------------")
    # print(len(node_coords_ini))
    # print("--------------------------------------------------")
    # print("minDistRemesh = %lg" % minDistRemesh)
    # print("--------------------------------------------------")

    # print(model.trials)

    # print(rect_lims)
    # exit(0)

    # add specified nodes if in regions to remesh
    for node in model.specifiedNodes:
        # print("reading specified node")
        # distance is True if point is outside
        # is outside regions already remeshed?
        if utilitiesGeom.checkMutDistancesLoops(model.dimension,
                                                    radiusRemesh, cpr,
                                                    list(node)):
            # print("is not in already remeshed")
            # is inside regions to remeshed?
            if not utilitiesGeom.checkMutDistancesLoops(model.dimension,
                                                    radiusTransitional,
                                                    ctr, list(node)):
                # if using existing fine nodes, append only is in transitional area
                if not useExistingFineNodes or utilitiesGeom.checkMutDistancesLoops(model.dimension,
                                                        radiusRemesh,
                                                        ctr, list(node)):
                    # print("is in regions to remesh")
                    node_coords_ini.append(node)
                    print("appending specified node")
                    print(node)
    # empty specified nodes, since they are already appended if necessary
    model.specifiedNodes = []

    node_coords = generateNodesRemesh(node_coords_ini.copy(),
                        trials=model.trials, maxLim=model.maxLim,
                        minDistRemesh=minDistRemesh, minDist=model.minDist,
                        centersToRemesh=ctr, centersPreviouslyRemeshed=cpr,
                        regionsToSkip=rtSkip,
                        radiusRemesh=radiusRemesh,
                        radiusTransitional=radiusTransitional,
                        dim=model.dimension,
                        useExistingFineNodes=useExistingFineNodes,
                        remesherSeed=remesherSeed)

    # exit(1)
    # print("--- after --------------------------------------")
    # print(len(node_coords_ini))
    # print(len(node_coords))

    # fig, ax = plt.subplots(figsize=(15, 8))
    #
    # ax.scatter( np.array(node_coords)[:len(node_coords_ini), 0],
    #             np.array(node_coords)[:len(node_coords_ini), 1], color='b')
    #
    # ax.scatter( np.array(node_coords)[len(node_coords_ini):, 0],
    #             np.array(node_coords)[len(node_coords_ini):, 1], color='r' )



    # sys.exit(1)

    # dirNam = '/home/jose/Soft/ParticleModel/TESTS/adaptivity_pokus' + '/TPB_test_adaptive_II'
    dirNam = remeshDir

    if model != None:
        for i in range (model.nr_models):
            print('\nCreating model #%d' %i)
            model.setDirectory(dirNam)
            model.createModel(node_coords_init=node_coords)  # TODO with existing nodes
            model.saveGeometry()
            model.saveRest(solver, prep_input_file, exporters)


    # ax.scatter( np.array(model.node_coords)[:, 0],
    #             np.array(model.node_coords)[:, 1],
    #             color='g', marker='x' )
    #
    # ax.scatter( np.array(node_coords_old)[:, 0],
    #             np.array(node_coords_old)[:, 1], color='m', marker='+' )

    # print(model.vert_count)
    # print("-----------------------------------------")
    # print(model.verticesIdxDict)
    # print("-----------------------------------------")
    # print(model.vertIdxStart)
    # print("-----------------------------------------")
    # print(model.totalNodeCount)

    # ss.voronoi_plot_2d(model.vor, ax=ax)
    #
    # if SHOW_PLOT:
    #     plt.show()
    # plt.close()
    # sys.exit(1)

    print('\n%%%%%%%%% LATTICE REMESH DONE %%%%%%%%%')
    #print('\n%%%%%%%%% %d NODES MODEL %%%%%%%%%' %len(model.node_coords))
    print('All done in %f seconds. (%d minutes).' %((time.time()-start), (time.time()-start)/60.0))
