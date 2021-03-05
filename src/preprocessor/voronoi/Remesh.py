import errno
import os
from ModelGenerator import *
# from regions import *
from pointGenerators import generateNodesRemesh
import numpy as np
import matplotlib.pyplot as plt
import scipy.spatial as ss

SHOW_PLOT = False

def loadNodes(filename, dim):
    if not os.path.isfile(filename):
        raise FileNotFoundError(
                errno.ENOENT, os.strerror(errno.ENOENT), filename)
    nodes = list()
    node_coords = np.genfromtxt(filename, skip_header=1,
                                usecols=( ( dim == 3 ) and (1, 2, 3) or (1, 2) )
                                )
    # print(node_coords)
    # print("--------------------------------------------------")
    # print([np.array(coor) for coor in node_coords.tolist()])
    # sys.exit(1)
    try:
        len(node_coords[0])
    except Exception as e:
        print(e)
        return [ node_coords,  ]
    return [np.array(coor) for coor in node_coords.tolist()]


def loadRegionsToKeep(region_file):

    return


if __name__ == '__main__':
    print('\n%%%%%%%%% LATTICE REMESH STARTED %%%%%%%%%')
    start = time.time()

    # print("-----------------------------------------------------")
    # print(len(sys.argv))
    # print("-----------------------------------------------------")

    if len(sys.argv) < 6:
        print("not enough information provided for Remesher")
        sys.exit(1)

    prep_input_file = os.path.join(sys.argv[3], sys.argv[1])
    print(prep_input_file)

    remeshDir = os.path.join(sys.argv[3], sys.argv[2])
    print(remeshDir)
    # sys.exit()

    if not os.path.isfile(prep_input_file):
        print('Missing input master file in argument! Exiting.')
        sys.exit()

    model = None
    solver = None

    skipLines = ['#', ' ', '\n', '\t']

    f = open (prep_input_file, 'r')
    for row in f:
        if not (row[0] in skipLines):
            r = row.split()

            if (r[0]=='Model'):
                model = Model(row)
            if (r[0]=='Solver'):
                solver = Solver(row)
            if (r[0]=='Material'):
                if model != None:
                    model.addMaterial(row)
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

    node_coords_ini = loadNodes(os.path.join(remeshDir, "nodes.inp"), model.dimension)

    ctr = loadNodes(os.path.join(remeshDir, "centersToRemesh.out"), model.dimension)
    # print("--------------------------------------------------")
    # print(ctr)
    # print("--------------------------------------------------")
    # where to load them from?
    # should be already in input (specified in master prep file)
    radiusRemesh = float(sys.argv[4])
    radiusTransitional = float(sys.argv[5])

    if len(sys.argv) > 6:
        minDistRemesh = float(sys.argv[6])
    else:
        minDistRemesh = model.minDist / 3.

    # print("--- before ------------------------------------")
    # print(len(node_coords_ini))
    # print("--------------------------------------------------")
    # print("minDistRemesh = %lg" % minDistRemesh)
    # print("--------------------------------------------------")

    cpr = list()

    # print(model.trials)



    rect_lims = [[0., 0.], [model.maxLim[0], model.maxLim[1]]]

    # model.maxLim[0] += 1e-8
    # model.maxLim[1] += 1e-8

    # print(rect_lims)
    # exit(0)

    node_coords = generateNodesRemesh(node_coords_ini.copy(),
                        trials=model.trials, maxLim=model.maxLim,
                        minDistRemesh=minDistRemesh, minDist=model.minDist,
                        centersToRemesh=ctr, centersPreviouslyRemeshed=cpr,
                        radiusRemesh=radiusRemesh,
                        radiusTransitional=radiusTransitional,
                        dim=2, rectLims=rect_lims)

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
            model.saveRest(solver, prep_input_file)


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
