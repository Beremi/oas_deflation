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




class Model:
    def __init__ (self, row):
        print('Setting model...', end='')
        self.seed = np.random.randint(1000.0)
        self.trials = 3e4

        self.powerTes = False
        self.activeMechanics = False
        self.activeTransport = False
        self.periodicModel = 0

        self.nodePositions = []
        self.coupledNodes = []
        self.mirtype = []

        self.govNodes = []
        self.rigidPlates = []
        self.constraint = False

        self.node_coords = []

        self.mechBC_merged = []
        self.mechIC_merged = []
        self.trsprtBC_merged = []
        self.trsprtIC_merged = []
        self.govNodesMechBC = []

        self.functions = []
        self.radii = []
        self.totalNodeCount = -1
        self.notches = []
        self.materialZones= []
        self.measuringGauges = []

        self.vor = None
        self.areas = None

        self.materials = []
        self.vert_count = []
        self.verticesIdxDict = None
        self.vertIdxStart = -1
        #"""
        r = row.split()
        for i in range (len(r)):
            self.modelType = r[1]

            if (r[i]=='seed'):
                self.seed = r[i+1]
            if (r[i]=='dimension'):
                self.dimension = int(r[i+1])
                self.maxLim = np.zeros((self.dimension))
            if (r[i]=='Xsize'):
                self.maxLim[0] = float(r[i+1])
            if (r[i]=='Ysize'):
                self.maxLim[1] = float(r[i+1])
            if (r[i]=='Zsize'):
                if (self.dimension==3):
                    self.maxLim[2] = float(r[i+1])
            if (r[i]=='minDist'):
                self.minDist = float(r[i+1])
            if (r[i]=='trials'):
                self.trials = int(r[i+1])

            if (r[i]=='powerTes'):
                if (int(r[i+1])==1): self.powerTes = True
                if (int(r[i+1])==0): self.powerTes = False
            if (r[i]=='activeMechanics'):
                if (int(r[i+1])==1): self.activeMechanics = True
                if (int(r[i+1])==0): self.activeMechanics = False
            if (r[i]=='activeTransport'):
                if (int(r[i+1])==1): self.activeTransport = True
                if (int(r[i+1])==0): self.activeTransport = False
            if (r[i]=='periodicModel'):
                if (int(r[i+1])==1): self.periodicModel = True
                if (int(r[i+1])==0): self.periodicModel = False

            if (r[i]=='cylinderRad'):
                self.cylinderRad = float(r[i+1])
            if (r[i]=='cylinderHeight'):
                self.cylinderHeight = float(r[i+1])
            if (r[i]=='tubeThickness'):
                self.tubeThickness = float(r[i+1])
            if (r[i]=='dogboneD'):
                self.tubeThickness = float(r[i+1])
            if (r[i]=='notchH'):
                self.notchH = float(r[i+1])
            if (r[i]=='loadWidth'):
                self.loadWidth = float(r[i+1])

        np.random.seed(seed=self.seed)
        self.master_folder = 'power_%.4f_%02d' % (self.minDist, self.seed)
        try:
            if not os.path.exists(self.master_folder):
                os.makedirs(self.master_folder)
        except:
            print('Please create directory %s! Code Exited.' % self.master_folder)
            sys.exit()

        print('done.')



    def createModel(self):
        print ('Creating model %s...' %self.modelType)
        #
        if self.modelType == '2d_notched3pb':
            self.run_2d_notched3pb()


    def run_2d_notched3pb(self):
        self.node_coords,
        self.mechBC_merged,
        self.mechIC_merged,
        self.vor,
        self.areas,
        self.functions,
        self.notches,
        self.govNodes,
        self.govNodesMechBC,
        self.rigidPlates  = utilitiesModeling.create2dSSBeamUnifLoad(self.maxLim, self.minDist, self.trials, notch=self.notchH, loadWidth=self.loadWidth)

        self.measuringGauges = utilitiesModeling.assembleMeasuringGauges('3pb2d', maxLim=self.maxLim)

        if (len(self.node_coords))==0:
            print('Nody se neprenesly!!')
            sys.exit()
        else:
            print(self.node_coords)
            sys.exit()






    def extractGeometry(self):
        self.vert_count,
        self.verticesIdxDict,
        self.vertIdxStart,
        self.totalNodeCount = utilitiesGeom.extractGeometry(
            self.master_folder,
            self.dimension,
            len(self.node_coords),
            self.maxLim, self.vor,
            self.node_coords,
            self.areas,
            self.activeTransport, self.activeMechanics,
            mZ=self.materialZones, periodicModel=self.periodicModel,
            nodePositions=self.nodePositions, coupledNodes=self.coupledNodes,
            mirtype=self.mirtype, notches=self.notches)




    def addMaterial(self, row):
        r = row.split()
        print ('Adding %s...' %r[1], end='')
        #LinearElasticMaterial
        if (r[1]=='DisMechMaterial'):
            #
            young = None
            alpha = None
            density = None
            #
            for i in range (len(r)):
                if (r[i]=='young'):
                    young = float(r[i+1])
                if (r[i]=='alpha'):
                    alpha = float(r[i+1])
                if (r[i]=='density'):
                    density = float(r[i+1])
            #
            if (young == None or alpha == None or density == None):
                print ('!! DisMechMaterial incomplete. Exiting. !!')
                sys.exit()

            linElMaterial = utilitiesMech.linearElasticMaterial(young, alpha, density)
            self.materials.append(linElMaterial)
            print('done.')

        #TransportMaterial
        if (r[1]=='TrsprtMaterial'):
            #
            transpC = None
            transpS = None
            #
            for i in range (len(r)):
                if (r[i]=='transpC'):
                    transpC = float(r[i+1])
                if (r[i]=='transpS'):
                    transpS = float(r[i+1])
            #
            if (transpC == None or transpS == None):
                print ('!! TrsprtMaterial incomplete. Exiting. !!')
                sys.exit()

            transportMaterial = utilitiesMech.TransportMaterial(transpC, transpS)
            self.materials.append(transportMaterial)
            print('done.')

        #MarsMaterial
        if (r[1]=='MarsMaterial'):
            #
            young = None
            alpha = None
            density = None
            ft = None
            Gt = None
            #
            for i in range (len(r)):
                if (r[i]=='young'):
                    young = float(r[i+1])
                if (r[i]=='alpha'):
                    alpha = float(r[i+1])
                if (r[i]=='density'):
                    density = float(r[i+1])
                if (r[i]=='ft'):
                    ft = float(r[i+1])
                if (r[i]=='Gt'):
                    Gt = float(r[i+1])
            #
            if (young == None or alpha == None or density == None or ft == None or  Gt == None):
                print ('!! MarsMaterial incomplete. Exiting. !!')
                sys.exit()

            marsMaterial = utilitiesMech.MarsMaterial(young, alpha, density, ft, Gt)
            self.materials.append(marsMaterial)
            print('done.')

        #FatigueMaterial
        #E0 35e9 alpha 0.300000 density 2200.0 fc 20e6 ft 3.5e6 KinN 4e9 gammaN 20e9 m -0.2e-6 Ad 4000e-6 tauBar 4.0e6 Kin 0 gamma 10.0e9 S 0.0025e3 a 0
        if (r[1]=='FatigueMaterial'):
            #
            young = None
            alpha = None
            density = None
            fc = None
            ft = None
            KinN = None
            gammaN = None
            m = None
            Ad = None
            tauBar = None
            Kin = None
            gamma = None
            S = None
            a = None
            #
            for i in range (len(r)):
                if (r[i]=='E0'):
                    young = float(r[i+1])
                if (r[i]=='alpha'):
                    alpha = float(r[i+1])
                if (r[i]=='density'):
                    density = float(r[i+1])
                if (r[i]=='fc'):
                    fc = float(r[i+1])
                if (r[i]=='ft'):
                    ft = float(r[i+1])
                if (r[i]=='KinN'):
                    KinN = float(r[i+1])
                if (r[i]=='gammaN'):
                    gammaN = float(r[i+1])
                if (r[i]=='m'):
                    m = float(r[i+1])
                if (r[i]=='Ad'):
                    Ad = float(r[i+1])
                if (r[i]=='tauBar'):
                    tauBar = float(r[i+1])
                if (r[i]=='Kin'):
                    Kin = float(r[i+1])
                if (r[i]=='gamma'):
                    gamma = float(r[i+1])
                if (r[i]=='S'):
                    S = float(r[i+1])
                if (r[i]=='a'):
                    a = float(r[i+1])
            #
            if (young == None or
            alpha == None or
            density == None or
            fc == None or
            ft == None or
            KinN == None or
            gammaN == None or
            m == None or
            Ad == None or
            tauBar == None or
            Kin == None or
            gamma == None or
            S == None or
            a == None):
                print ('!! FatigueMaterial incomplete. Exiting. !!')
                sys.exit()

            fatigueMaterial = utilitiesMech.FatigueMaterial(young, alpha, density, fc, ft, KinN, gammaN, m, Ad, tauBar, Kin, gamma, S, a)
            self.materials.append(fatigueMaterial)
            print('done.')

    """
    def run_2d_notched3pb(self):
        self.node_coords,
        self.mechBC_merged,
        self.mechIC_merged,
        self.vor,
        self.areas,
        self.functions,
        self.notches,
        self.govNodes,
        self.govNodesMechBC,
        self.rigidPlates  = utilitiesModeling.create2dSSBeamUnifLoad(self.maxLim, self.minDist, self.trials, notch=self.notchH, loadWidth=self.loadWidth)

        self.materialZones=None
        self.measuringGauges = utilitiesModeling.assembleMeasuringGauges('3pb2d', maxLim=self.maxLim)

        print(node_coords)
    #"""

















class Solver:
    def __init__(self, row):
        print('Setting solver...', end='')
        #setting default solver
        self.solverType = 'SteadyStateNonLinearSolver'
        self.time_step = 1e-6
        self.max_time_step = 1e-1
        self.min_time_step = 1e-8
        self.total_time = 1
        self.limit_tolerance = 1e-5
        self.maxIt = 20

        r = row.split()
        for i in range (len(r)):
            self.solverType = r[1]

            if (r[i]=='time_step'):
                self.time_step = float(r[i+1])
            if (r[i]=='max_time_step'):
                self.max_time_step = float(r[i+1])
            if (r[i]=='min_time_step'):
                self.min_time_step = float(r[i+1])
            if (r[i]=='total_time'):
                self.total_time = float(r[i+1])
            if (r[i]=='limit_tolerance'):
                self.limit_tolerance = float(r[i+1])
            if (r[i]=='maxIt'):
                self.maxIt = int(r[i+1])

        print('done.')
        #return solverType, time_step, max_time_step, min_time_step, total_time, limit_tolerance, maxIt





















if __name__ == '__main__':
    print('\n%%%%%%%%% LATTICE PREPROCESSOR STARTED %%%%%%%%%')
    start = time.time()

    if len(sys.argv) == 2:
        file = sys.argv[1]
        print('Loading master file %s' %file)
        """
        f = open (file, 'r')
        for row in f:
            print(row)
        #"""
    else:
        print('Missing input master file in argument! Exiting.')
        sys.exit()

    #print('Extracting data from master file...')
    f = open (file, 'r')
    for row in f:
        if not (row[0]=='#'):
            r = row.split()
            #
            if (r[0]=='Model'):
                model = Model(row)
            if (r[0]=='Solver'):
                solver = Solver(row)
            if (r[0]=='Material'):
                model.addMaterial(row)


    if Model != None:
        model.createModel()
        model.extractGeometry()


    print('\n%%%%%%%%% LATTICE PREPROCESSOR DONE %%%%%%%%%')
    print('All done in %f seconds.' %(time.time()-start))




























#
