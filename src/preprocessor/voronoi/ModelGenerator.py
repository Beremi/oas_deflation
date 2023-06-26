import matplotlib as mpl
#mpl.use('Agg')
import Preprocessor, sys, time, numpy as np, random
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import scipy
from IPython.display import clear_output
import sys, os, math
from shutil import copyfile
from sklearn import preprocessing
from scipy.ndimage import rotate
from scipy.spatial import Voronoi
from scipy.spatial import voronoi_plot_2d
from scipy.spatial import Delaunay
from scipy.sparse.csgraph import reverse_cuthill_mckee
from scipy.sparse import csr_matrix
from scipy.sparse import csc_matrix
import utilitiesGeom, utilitiesMech, utilitiesModeling, utilitiesNumeric, voronoi
from utilitiesGeom import mechBCFile
from pathlib import Path
import shutil


# Disable
def blockPrint():
    sys.stdout = open(os.devnull, 'w')

# Restore
def enablePrint():
    sys.stdout = sys.__stdout__



class Model:
    def __init__ (self, row):
        #print('Setting model...', end='')

        self.nr_models = 1
        self.printout = False
        self.userSeed = -1
        self.seed = -1
        self.master_folder = ''

        self.specifiedNodes = []

        self.auxmechelements = False

        self.trials = 3e4

        self.powerTes = self.activeMechanics = self.activeTransport = self.coupled = False
        self.periodicModel = 0

        self.nodePositions = []
        self.coupledNodes = []
        self.mirtype = []

        self.govNodes = []
        self.govNodesTrspt = []
        self.rigidPlates = []
        self.rigidPlatesTrspt = []
        self.constraint = False
        self.constraintTrspt = False

        self.edgeMinDistCoef=1.0

        self.node_coords = []

        self.node_coords_polar = []

        self.mechBC_merged = self.mechIC_merged = self.trsprtBC_merged = self.trsprtIC_merged = self.govNodesMechBC = self.govNodesTrsptBC =self.expansionRings= None

        self.tpbwedgeheight = 0.2

        self.loading = '3pb'

        self.functions = []
        self.radii = []
        self.totalNodeCount = -1
        self.notches = []
        self.materialZones= []
        self.measuringGauges = []
        self.symmetric = 0
        self.weakboundary = 0

        self.ansysOut = False

        self.vor = None
        self.areas = None

        self.materials = []
        self.vert_count = []
        self.verticesIdxDict = None
        self.vertIdxStart = -1
        #"""
        r = row.split()
        self.modelType = r[1]
        self.dimension = int(r[1][0])
        self.maxLim = np.zeros((self.dimension))    # float for circRVE

        self.masters = []

        self.defaultFilesFolder = None

        #RWTH cylinderRad
        #cylinderHeight 0.1 cylinderRad 0.1 minDist 0.02 notchRadLeft 0.08 notchRadRight 0.05 notchWidth 0.01
        self.notchRadLeft = None
        self.notchRadRight = None
        self.innerRadTop = None
        self.innerRadBottom = None
        self.notchWidth = None
        self.notchDepth = None
        self.RWTHQuarter = False
        self.elasticZone = 0

        self.node_indices_dogbone = []
        self.interfaceVertexIndices = []

        self.rebarMinDist = None
        self.rebarDiameter = None
        self.rebarCount = None
        self.rebarDepth = None
        self.interfaceMinDist = None

        self.fracZoneHeight = None
        self.fracZoneOverhang = None
        self.fem = False


        self.masterSolver = self.masterMaterials = self.masterFunctions = False

        self.supportDivision = None

        for i in range (len(r)):
            if (r[i]=='supportDivision'):
                self.supportDivision = int(r[i+1])
            if (r[i]=='holeMinDist'):
                self.holeMinDist = float(r[i+1])
            if (r[i]=='fracZoneHeight'):
                self.fracZoneHeight = float(r[i+1])
            if (r[i]=='fracZoneOverhang'):
                self.fracZoneOverhang = float(r[i+1])
            if (r[i]=='holeDiameter'):
                self.holeDiameter = float(r[i+1])
            if (r[i]=='rebarMinDist'):
                self.rebarMinDist = float(r[i+1])
            if (r[i]=='rebarDiameter'):
                self.rebarDiameter = float(r[i+1])
            if (r[i]=='rebarCount'):
                self.rebarCount = int(r[i+1])
            if (r[i]=='rebarDepth'):
                self.rebarDepth = float(r[i+1])
            if (r[i]=='interfaceMinDist'):
                self.interfaceMinDist = float(r[i+1])
            if (r[i]=='fineRingThickness'):
                self.fineRingThickness = float(r[i+1])
            if (r[i]=='fineRegDepth'):
                self.fineRegDepth = float(r[i+1])
            if (r[i]=='fineWidth'):
                self.fineWidth = float(r[i+1])
            if (r[i]=='gradientRegDepth'):
                self.gradientRegDepth = float(r[i+1])

            if (r[i]=='defaultFilesFolder'):
                self.defaultFilesFolder = str(r[i+1])


            if (r[i]=='printout'):
                if r[i+1] == 1:
                    self.printout = True
            if (r[i]=='fem'):
                if int(r[i+1]) == 1:
                    self.fem = True

            if (r[i]=='ansysOut'):
                if int(r[i+1]) == 1:
                    self.ansysOut = True



            if (r[i]=='nr_models'):
                self.nr_models = int(r[i+1])

            if (r[i]=='seed'):
                self.userSeed = int(r[i+1])
            if (r[i]=='Xsize'):
                self.maxLim[0] = float(r[i+1])
            if (r[i]=='Ysize'):
                self.maxLim[1] = float(r[i+1])
            if (r[i]=='Zsize'):
                if (self.dimension==3):
                    self.maxLim[2] = float(r[i+1])
            if (r[i] == 'RVEdiameter'):
                self.maxLim = float(r[i+1])
            if (r[i]=='minDist'):
                self.minDist = float(r[i+1])
            if (r[i]=='trials'):
                self.trials = int(r[i+1])

            if (r[i]=='tpbwedgeheight'):
                self.tpbwedgeheight = float(r[i+1])

            self.randomizeMaterial = False
            if (r[i]=='randomizeMaterial'):
                if (int(r[i+1])==1): self.randomizeMaterial = True
                if (int(r[i+1])==0): self.randomizeMaterial = False

            if (r[i]=='auxmechelements'):
                if (int(r[i+1])==1): self.auxmechelements = True
                if (int(r[i+1])==0): self.auxmechelements = False

            if (r[i]=='masterSolver'):
                if (int(r[i+1])==1): self.masterSolver = True
                if (int(r[i+1])==0): self.masterSolver = False

            if (r[i]=='masterMaterials'):
                if (int(r[i+1])==1): self.masterMaterials = True
                if (int(r[i+1])==0): self.masterMaterials = False

            if (r[i]=='masterFunctions'):
                if (int(r[i+1])==1): self.masterFunctions = True
                if (int(r[i+1])==0): self.masterFunctions = False


            if (r[i]=='powerTes'):
                if (int(r[i+1])==1): self.powerTes = True
                if (int(r[i+1])==0): self.powerTes = False
            if (r[i]=='activeMechanics'):
                if (int(r[i+1])==1): self.activeMechanics = True
                if (int(r[i+1])==0): self.activeMechanics = False
            if (r[i]=='activeTransport'):
                if (int(r[i+1])==1): self.activeTransport = True
                if (int(r[i+1])==0): self.activeTransport = False
            if (r[i]=='coupled'):
                if (int(r[i+1])==1): self.activeMechanics = True
                if (int(r[i+1])==1): self.activeTransport = True
                if (int(r[i+1])==1): self.coupled = True
            if (r[i]=='periodicModel'):
                if (int(r[i+1])==1): self.periodicModel = True
                if (int(r[i+1])==0): self.periodicModel = False
            if (r[i]=='elasticZone'):
                self.elasticZone = float(r[i+1])
            #"""
            """
            keys = ['powerTes', 'activeMechanics', 'activeTransport', 'periodicModel']
            if r[i] in keys:
                setattr(self, r[i], bool(r[i+1]))
            """

            if (r[i]=='roughDogBone'):
                self.roughDogBone = int(r[i+1])
            if (r[i]=='roughEdgeDogbone'):
                self.roughEdgeDogbone = int(r[i+1])

            if (r[i]=='roughMinDistCoef'):
                self.roughMinDistCoef = float(r[i+1])

            if (r[i]=='interLayerThickness'):
                self.interLayerThickness = float(r[i+1])

            if (r[i]=='loading'):
                self.loading = (r[i+1])

            if (r[i]=='edgeMinDistCoef'):
                self.edgeMinDistCoef = float(r[i+1])
            if (r[i]=='elasticHeightCoef'):
                self.elasticHeightCoef = float(r[i+1])

            if (r[i]=='dmgBand'):
                self.dmgBand = float(r[i+1])

            if (r[i]=='cylinderRad'):
                self.cylinderRad = float(r[i+1])
            if (r[i]=='cylinderHeight'):
                self.cylinderHeight = float(r[i+1])
            if (r[i]=='notchRadLeft'):
                self.notchRadLeft = float(r[i+1])
            if (r[i]=='notchRadRight'):
                self.notchRadRight = float(r[i+1])
            if (r[i]=='innerRadTop'):
                self.innerRadTop = float(r[i+1])
            if (r[i]=='innerRadBottom'):
                self.innerRadBottom = float(r[i+1])
            if (r[i]=='notchWidth'):
                self.notchWidth = float(r[i+1])
            if (r[i]=='notchDepth'):
                self.notchDepth = float(r[i+1])
            if (r[i]=='RWTHQuarter'):
                if (int(r[i+1])==1): self.RWTHQuarter = True
                if (int(r[i+1])==0): self.RWTHQuarter = False

            if (r[i]=='tubeThickness'):
                self.tubeThickness = float(r[i+1])
            if (r[i]=='dogboneD'):
                self.dogboneD = float(r[i+1])
            if (r[i]=='dogboneExcentricityFrac'):
                self.dogboneExcentricityFrac = float(r[i+1])
            if (r[i] == 'dogboneExcentricity_Z'):
                self.dogboneExcentricity_Z = float(r[i + 1])
            if (r[i]=='notchH'):
                self.notchH = float(r[i+1])
            if (r[i]=='loadWidth'):
                self.loadWidth = float(r[i+1])
            if (r[i]=='fracZoneWidth'):
                self.fracZoneWidth = float(r[i+1])
            if (r[i]=='orthogonalFracZone'):
                if (int(r[i+1])==1): self.orthogonalFracZone = True
                if (int(r[i+1])==0): self.orthogonalFracZone = False
            if (r[i]=='Xtopsize'):
                self.Xtopsize = float(r[i+1])
            if (r[i]=='symmetric'):
                self.symmetric =  int(r[i+1])
            if (r[i]=='weakboundary'):
                self.weakboundary =  float(r[i+1])

            if (r[i]=='adaptivityReady'):
                if (int(r[i+1])==1): self.adaptivityReady = True
                if (int(r[i+1])==0): self.adaptivityReady = False

        print('done.')

    def setDirectory(self, dirNam=None):
        if self.userSeed == -1:
            self.seed = np.random.randint(100000.0)
            np.random.seed(seed=self.seed)
        else:
            self.seed = self.userSeed
            np.random.seed(seed=self.seed)

        if dirNam is None:
            self.master_folder = '%s_minDist%.4f_seed%02d' % (self.modelType,self.minDist, self.seed)
        else:
            self.master_folder = dirNam

        try:
            if not os.path.exists(self.master_folder):
                os.makedirs(self.master_folder)
        except:
            print('Please create directory %s! Code Exited.' % self.master_folder)
            sys.exit(1)


    def createModel(self, node_coords_init=None, modelnr=-1):
        print ('Creating model %s' %self.modelType)
        #
        #if (self.printout == False): blockPrint()
        if self.modelType == '2d_transportPatchTest':
            self.run_2d_transportPatchTest()
        if self.modelType == '3d_transportPatchTest':
            self.run_3d_transportPatchTest()
        if self.modelType == '3d_BiparvaTubeTransport':
            self.run_3d_BiparvaTubeTransport()

        if self.modelType == '2d_notched3pb':
            self.run_2d_notched3pb(node_coords_init=node_coords_init)
        if self.modelType == '3d_notched3pb':
            self.run_3d_notched3pb(node_coords_init=node_coords_init, modelnr=modelnr)

        if self.modelType == '2d_coupledCompression':
            self.run_2d_coupledCompression()
        if self.modelType == '3d_coupledCompression':
            self.run_2d_coupledCompression()

        if self.modelType == '2d_dogbone':
            self.run_2d_dogbone()
        if self.modelType == '2d_dogboneStrip':
            self.run_2d_dogboneStrip()
        if self.modelType == '2d_dogboneBand':
            self.run_2d_dogboneBand()

        if self.modelType == '3d_dogbone':
            self.run_3d_dogbone()

        if self.modelType == '3d_cylinderTorsionPress':
            self.run_3d_torsionPress()

        if self.modelType == '3d_ReinhardtTension':
            self.run_3d_ReinhardtTension()

        if self.modelType == '3d_RWTHShearCylinder':
            self.run_3d_RWTHShearCylinder()
        if self.modelType == '2d_RWTHShearCylinder':
            self.run_2d_RWTHShearCylinder()

        if self.modelType == '2d_periodicShear':
            self.run_2d_periodicShear()

        if self.modelType == '3d_periodicShear':
            self.run_3d_periodicShear()

        if self.modelType == '3d_cube':
            self.run_3d_cube()

        if self.modelType == '3d_consolidation':
            self.run_3d_consolidation()

        if self.modelType == '2d_singleSpring':
            self.run_2d_singleSpring()

        if self.modelType == '2d_singleContact':
            self.run_2d_singleContact()
        if self.modelType == '3d_singleContact':
            self.run_3d_singleContact()

        if self.modelType == '2d_coupledArtificialCrack':
            self.run_2d_coupledArtificialCrack()

        if self.modelType == '3d_coupledArtificialCrack':
            self.run_3d_coupledArtificialCrack()

        if self.modelType == '2d_CFRAC_Clover':
            self.run_2d_CFRAC_Clover()

        if self.modelType == '3d_CFRAC_Clover':
            self.run_3d_CFRAC_Clover()

        if self.modelType == '2d_CFRAC_TDCB':
            self.run_2d_CFRAC_TDCB()
        if self.modelType == '3d_CFRAC_TDCB':
            self.run_3d_CFRAC_TDCB()
        if self.modelType == '2d_Hanging_FracZone':
            self.run_2d_Hanging_FracZone()
        if self.modelType == '3d_Hanging_FracZone':
            self.run_3d_Hanging_FracZone()

        if self.modelType == '2d_corrosionRebar':
            self.run_2d_corrosionRebar(node_coords_init=node_coords_init)
        if self.modelType == '3d_corrosionRebar':
            self.run_3d_corrosionRebar(node_coords_init=node_coords_init)

        if self.modelType == '3d_TubeInnerPressure':
            self.run_3d_TubeInnerPressure()

        if self.modelType == '2d_coupledRVE':
            self.run_2d_coupledRVE()

        if self.modelType == '2d_circRVE':
            self.run_2d_circRVE()

        if self.modelType == '2d_coupledPress':
            self.run_2d_coupledPress()

        if self.modelType =='3d_coupledBrazilianDisc':
            self.run_3d_coupledBrazilianDisc()

        if self.modelType =='2d_cantileverBending':
            self.run_2d_cantileverBending()

        if self.modelType == '3d_dam':
            self.run_3d_dam()


        #if (self.printout == False): enablePrint()

    def run_2d_singleSpring(self):
        (self.node_coords, self.mechBC_merged,  self.vor, self.areas, self.functions,self.govNodes, self.govNodesMechBC, self.rigidPlates) = utilitiesModeling.createSingleSpringTestModel( self.minDist, self.master_folder )
        self.measuringGauges = utilitiesModeling.assembleMeasuringGauges('2d_singleSpring', maxLim=np.array([self.minDist, 1]) )
        materialZones = None

    def run_2d_singleContact(self):
        self.maxLim=np.array([self.minDist*2,self.minDist])
        (self.node_coords, self.mechBC_merged,  self.vor, self.areas, self.functions,self.govNodes, self.govNodesMechBC, self.rigidPlates) = utilitiesModeling.create_2d_SingleContat( self.minDist, self.master_folder,self.maxLim )
        self.measuringGauges = utilitiesModeling.assembleMeasuringGauges('2d_singleSpring', maxLim=self.maxLim )
        materialZones = None
    def run_3d_singleContact(self):
        print(self.dimension)
        self.maxLim=np.array([self.minDist*2,self.minDist,self.minDist])
        (self.node_coords, self.mechBC_merged,  self.vor, self.areas, self.functions,self.govNodes, self.govNodesMechBC, self.rigidPlates) = utilitiesModeling.create_3d_SingleContat( self.minDist, self.master_folder,self.maxLim )
        self.measuringGauges = utilitiesModeling.assembleMeasuringGauges('3d_singleContact', maxLim=self.maxLim )
        materialZones = None

    def run_2d_notched3pb(self, node_coords_init=None):
        supportWidth = self.maxLim[0] / 20

        # hned tady na zacatku se model rozsiri o sirku podpor, aby skutecne rozpeti podpor se rovnalo zadanemu X lim
        #na kazdou stranu se prida pulka podpory
        self.maxLim[0] = self.maxLim[0] + 2 * 0.5 * supportWidth
        #print('rozsireni modelu')

        (self.node_coords, self.mechBC_merged, self.mechIC_merged, self.vor, self.areas, self.functions, self.notches, self.govNodes,
        self.govNodesMechBC, self.rigidPlates, self.trsprtBC_merged)     = utilitiesModeling.create2dSSBeamUnifLoad(self.maxLim, self.minDist, self.trials, notch=self.notchH, loadWidth=self.loadWidth, fracZoneWidth = self.fracZoneWidth, orthogonalFracZone=self.orthogonalFracZone, notchWidth=self.notchWidth, node_coords_init=node_coords_init, activeTransport=self.activeTransport, coupled=self.coupled, specifiedNodes=self.specifiedNodes, loading=self.loading)
        self.measuringGauges = utilitiesModeling.assembleMeasuringGauges('3pb2d', maxLim=self.maxLim)
        if self.loading == "4pb":
            self.materialZones = utilitiesModeling.assembleMaterialZones(0, 2, model='4pb2d', maxLim=self.maxLim, minDist=self.minDist)

    def run_3d_notched3pb(self, node_coords_init=None,modelnr=-1):
        print('running notched')
        #width of the supports, nemenit, je i v assemble3DSSBeamBending !!!!
        supportWidth = self.maxLim[0] / 20

        if modelnr ==0 :
            # hned tady na zacatku se model rozsiri o sirku podpor, aby skutecne rozpeti podpor se rovnalo zadanemu X lim
            #na kazdou stranu se prida pulka podpory
            self.maxLim[0] = self.maxLim[0] + 2 * 0.5 * supportWidth
            #print('rozsireni modelu')

        #(self.node_coords, self.mechBC_merged, self.mechIC_merged, self.vor, self.areas, self.functions, self.notches, self.govNodes, self.govNodesMechBC, self.rigidPlates, self.trsprtBC_merged, self.trsprtIC_merged)
        #(self.node_coords, self.mechBC_merged, self.mechInitC_merged,  self.vor, self.areas, self.functions, self.notches, self.govNodes, self.govNodesMechBC, self.rigidPlates, self.rigidPlatesTrspt, self.govNodesTrspt, self.govNodesTrsptBC, self.trsprtBC_merged)
        (self.node_coords, self.mechBC_merged, self.trsprtBC_merged, self.govNodes, self.govNodesMechBC, self.rigidPlates, self.vor, self.areas, self.functions, self.rigidPlatesTrspt, self.govNodesTrspt, self.govNodesTrsptBC, self.radii, self.notches)  = utilitiesModeling.create3dSSBeamUnifLoad(self.maxLim, self.minDist, self.trials, notch=self.notchH, loadWidth=self.loadWidth, fracZoneWidth = self.fracZoneWidth, orthogonalFracZone=self.orthogonalFracZone, notchWidth=self.notchWidth, coupled=self.coupled, node_coords_init=node_coords_init, specifiedNodes=self.specifiedNodes, roughMinDistCoef=self.roughMinDistCoef, supportDivision=self.supportDivision)
        self.measuringGauges = utilitiesModeling.assembleMeasuringGauges('3pb3d', maxLim=self.maxLim)



        self.materialZones = utilitiesModeling.assembleMaterialZones(0, 2, model='3pb3d', maxLim=self.maxLim, minDist=self.minDist)

        if self.elasticZone ==1:
            lim = np.array([
            self.maxLim[0]*0.5  + self.maxLim[1] * 0.2,#+ self.maxLim[0]*self.fracZoneWidth*1,
            self.maxLim[1] * (1-self.tpbwedgeheight),
            -self.maxLim[2]* 0.1,
            self.maxLim[0]*0.5  - self.maxLim[1] * 0.2,#- self.maxLim[0]*self.fracZoneWidth*1,
            self.maxLim[1] * 1.1,
            self.maxLim[2] * 1.1
            ])
            #print (lim)
            lim1=np.array([
            self.maxLim[0]*0.5  - self.maxLim[0] * self.fracZoneWidth * 0.4,
            self.maxLim[1] * 0.0,
            -self.maxLim[2]* 0.1,
            self.maxLim[0]*0.5  + self.maxLim[0] * self.fracZoneWidth * 0.4,
            self.maxLim[1] * self.notchH*1.8,
            self.maxLim[2]* 1.1,
            ])
            # print (lim1)
            self.materialZones = utilitiesModeling.assembleMaterialZones(0,3, model='3pb3d', limits=lim, limits1=lim1, maxLim=self.maxLim, minDist=self.minDist)




    def run_3d_cube(self, node_coords_init=None):
        (self.node_coords, self.mechBC_merged, self.mechIC_merged, self.vor, self.areas, self.functions, self.govNodes, self.govNodesMechBC, self.rigidPlates, self.trsprtBC_merged, self.trsprtIC_merged) = utilitiesModeling.create3dCube(self.maxLim, self.minDist, self.trials, self.powerTes, coupled=self.coupled, node_coords_init=node_coords_init )
        self.materialZones=None
        self.measuringGauges = utilitiesModeling.assembleMeasuringGauges('box3d', maxLim=self.maxLim)


    def run_3d_consolidation(self, node_coords_init=None):
        (self.node_coords, self.mechBC_merged, self.mechIC_merged, self.vor, self.areas, self.functions, self.govNodes, self.govNodesMechBC, self.rigidPlates, self.trsprtBC_merged, self.trsprtIC_merged, self.radii) = utilitiesModeling.create3dConsolidation(self.maxLim, self.minDist, self.trials, self.powerTes, coupled=self.coupled, node_coords_init=node_coords_init )
        self.materialZones=None

    def run_2d_dogbone(self):
        (self.node_coords,self.mechBC_merged,self.mechIC_merged,self.trsprtBC_merged,self.trsprtIC_merged,self.vor,self.areas,self.functions,self.govNodes,self.govNodesMechBC,self.rigidPlates, self.node_indices_dogbone)   = utilitiesModeling.create2dDogBone(self.minDist, self.trials, D=self.dogboneD, excentricity=self.dogboneExcentricityFrac, symmetric=self.symmetric, edgeMinDistCoef=self.edgeMinDistCoef, roughDogBone=self.roughDogBone, roughEdgeDogbone = self.roughEdgeDogbone, roughMinDistCoef=self.roughMinDistCoef, interLayerThickness=self.interLayerThickness, powerTes = self.powerTes, weakboundary = self.weakboundary)
        self.materialZones=None
        if self.elasticZone > 0:
            elaHeight = 1/4*self.dogboneD
            if self.roughDogBone >1:
                elaHeight = 3/4*self.dogboneD - (self.roughDogBone-1)*self.minDist
            if self.symmetric:
                if self.dmgBand == 1:
                    elaHeight = 3 / 4 * self.dogboneD - self.minDist / 4

                self.materialZones = utilitiesModeling.assembleMaterialZones(elaHeight, 2, model='dogboneStrip', D=self.dogboneD)
            else:
                self.materialZones = utilitiesModeling.assembleMaterialZones(elaHeight, 2, model='dogbone', D=self.dogboneD)
        if self.weakboundary > 0:
            print("Weak boundary activated - material at second row is taken as boundary material...\n")
            self.materialZones = utilitiesModeling.assembleMaterialZones(self.weakboundary, 2, model='dogbone', D=self.dogboneD, weakboundary = True)
        self.measuringGauges = utilitiesModeling.assembleMeasuringGauges('dogbone2d', D=self.dogboneD)


    def run_2d_dogboneStrip(self):
        (self.node_coords,self.mechBC_merged,self.mechIC_merged,self.trsprtBC_merged,self.trsprtIC_merged,self.vor,self.areas,self.functions,self.govNodes,self.govNodesMechBC,self.rigidPlates)   = utilitiesModeling.create2dDogBoneStrip(self.minDist, D=self.dogboneD, excentricity=self.dogboneExcentricityFrac, randomStrip=self.roughDogBone)
        self.materialZones=None

        self.measuringGauges = utilitiesModeling.assembleMeasuringGauges('dogbone2dStrip', D=self.dogboneD)



    def run_2d_dogboneBand(self):


        (self.node_coords,self.mechBC_merged,self.mechIC_merged,self.trsprtBC_merged,self.trsprtIC_merged,self.vor,self.areas,self.functions,self.govNodes,self.govNodesMechBC,self.rigidPlates)   = utilitiesModeling.create2dDogBoneBand(self.maxLim, self.minDist, roughMinDistCoef=self.roughMinDistCoef, elasticHeightCoef=self.elasticHeightCoef)

        elasticHeight = self.maxLim[1] *  self.elasticHeightCoef / (2* self.elasticHeightCoef + 1)
        dmgHeight =  self.maxLim[1] * 1 / (2* self.elasticHeightCoef + 1)
        width =  self.maxLim[0]
        self.materialZones= utilitiesModeling.assembleMaterialZones(0, 2, model='dogboneBand', maxLim=np.array([elasticHeight, dmgHeight, width]) )



        self.measuringGauges = utilitiesModeling.assembleMeasuringGauges('dogbone2dBand', maxLim= self.maxLim)



    def run_3d_dogbone(self):
        self.maxLim = np.array([self.dogboneD, 6/4*self.dogboneD, 0.1])
        #self.materialZones = utilitiesModeling.assembleMaterialZones (self.minDist*2, 3, model='dogbone',  D=self.dogboneD, thickness=0.1)
        self.materialZones = None
        (self.node_coords, self.mechBC_merged, self.mechIC_merged, self.trsprtBC_merged, self.trsprtIC_merged, self.vor, self.areas, self.functions, self.govNodes, self.govNodesMechBC, self.rigidPlates, self.node_indices_dogbone)   = utilitiesModeling.create3dDogBone(self.minDist, self.trials, D=self.dogboneD, excentricity_X=self.dogboneExcentricityFrac, excentricity_Z=self.dogboneExcentricity_Z, symmetric=self.symmetric)
        self.measuringGauges = utilitiesModeling.assembleMeasuringGauges('dogbone3d', D=self.dogboneD)

    def run_3d_torsionPress(self):
        self.maxLim = np.array([self.cylinderHeight, 2*self.cylinderRad, 2*self.cylinderRad])
        #self.materialZones = utilitiesModeling.assembleMaterialZones (self.minDist*2, 3, model='box', maxLim=self.maxLim)
        (self.node_coords, self.mechBC_merged,  self.vor, self.areas, self.functions, self.govNodes, self.govNodesMechBC, self.rigidPlates,
        self.trsprtBC_merged, self.govNodesTrspt, self.rigidPlatesTrspt, self.govNodesTrsptBC, self.radii)    = utilitiesModeling.create3dcylinderTorsionPressFree(np.zeros(3), self.cylinderRad, self.cylinderHeight,  self.minDist, self.trials, 0 , self.activeTransport, powerTes=self.powerTes)
        self.measuringGauges = utilitiesModeling.assembleMeasuringGauges('cylinder3d', maxLim=self.maxLim)
        self.materialZones =  None
        if self.elasticZone >0:
            self.materialZones = utilitiesModeling.assembleMaterialZones(self.elasticZone,3, model='3dcylinder', maxLim=self.maxLim)

    def run_3d_ReinhardtTension(self):
        #DURHAM - prism tension 250x60x50
        #Reinhardt "Tensile tests and failure analysis of concrete 1986"
        #maxLim = np.array([0.25,0.06,0.05])
        print(self.maxLim)
        self.materialZones = utilitiesModeling.assembleMaterialZones (self.minDist*2, self.dimension, model='box', maxLim=self.maxLim)
        (self.node_coords, self.mechBC_merged, self.mechIC_merged, self.vor, self.areas, self.functions, self.notches, self.govNodes, self.govNodesMechBC, self.rigidPlates) = utilitiesModeling.create3dReinhardtTension(self.maxLim, self.minDist, self.trials, fracZoneWidth=self.fracZoneWidth)
        self.measuringGauges = utilitiesModeling.assembleMeasuringGauges('reinhardt3d', maxLim=self.maxLim)

    def run_3d_RWTHShearCylinder(self):
        self.maxLim = np.array([self.cylinderHeight, 2*self.cylinderRad, 2*self.cylinderRad])
        #self.materialZones = utilitiesModeling.assembleMaterialZones (self.minDist*2, 3, model='box', maxLim=self.maxLim)
        (self.node_coords, self.mechBC_merged,  self.vor, self.areas, self.functions, self.govNodes, self.govNodesMechBC, self.rigidPlates, self.notches)    = utilitiesModeling.create3dRWTHShearCylinder(np.array([1e-6,0,0]), self.cylinderRad, self.cylinderHeight, self.minDist, self.trials, self.notchRadLeft, self.notchRadRight, self.notchWidth, self.notchWidth, self.notchDepth, quarter = self.RWTHQuarter )
        self.measuringGauges = utilitiesModeling.assembleMeasuringGauges('cylinder3d', maxLim=self.maxLim)
        self.materialZones =  None


    def run_2d_RWTHShearCylinder(self):
        self.maxLim = np.array([2*self.cylinderRad, self.cylinderHeight])
        #self.materialZones = utilitiesModeling.assembleMaterialZones (self.minDist*2, 3, model='box', maxLim=self.maxLim)
        (self.node_coords, self.mechBC_merged,  self.vor, self.volumes, self.functions,self.govNodes, self.govNodesMechBC, self.rigidPlates, self.notches)    = utilitiesModeling.create2dRWTHShearCylinder(self.cylinderRad, self.cylinderHeight, self.minDist, self.trials, self.innerRadTop, self.innerRadBottom, self.notchWidth, self.notchDepth)
        self.measuringGauges = utilitiesModeling.assembleMeasuringGauges('2d_RWTH', maxLim=self.maxLim)
        self.materialZones =  None

    def run_2d_periodicShear(self):
        (self.node_coords, self.mechBC_merged, self.mechIC_merged, self.trsprtBC_merged, self.trsprtIC_merged, self.vor, self.areas, self.functions, self.radii)   = utilitiesModeling.create2dPeriodicShear(self.maxLim, self.minDist, self.trials, self.powerTes )
        self.materialZones=None
        self.periodicModel = 1

    def run_3d_periodicShear(self):
        (self.node_coords, self.mechBC_merged, self.mechIC_merged, self.trsprtBC_merged, self.trsprtIC_merged, self.vor, self.areas, self.functions, self.radii)   = utilitiesModeling.create3dPeriodicShear(self.maxLim, self.minDist, self.trials, self.powerTes )
        self.materialZones=None
        self.periodicModel = 1

    def run_2d_transportPatchTest(self):
        (self.node_coords, self.mechBC_merged, self.trsprtBC_merged, self.vor, self.areas, self.functions, self.radii)  = utilitiesModeling.createPatchTestTransport(self.maxLim, self.minDist, self.trials, self.dimension, self.powerTes)

    def run_2d_coupledCompression(self):
        (self.node_coords, self.mechBC_merged, self.trsprtBC_merged, self.vor, self.areas, self.functions, self.radii)  = utilitiesModeling.createCoupledCompression(self.maxLim, self.minDist, self.trials, self.dimension, self.powerTes)

    def run_3d_transportPatchTest(self):
        (self.node_coords, self.mechBC_merged, self.trsprtBC_merged, self.vor, self.areas, self.functions, self.radii)  = utilitiesModeling.createPatchTestTransport(self.maxLim, self.minDist, self.trials, self.dimension, self.powerTes)


    def run_3d_BiparvaTubeTransport(self):
        #node_coords, mechBC_merged, govNodes, govNodesMechBC, rigidPlates, transportBC_merged, vor, volumes, functions, radii
        self.maxLim = np.array([self.cylinderHeight, self.cylinderRad, self.cylinderRad])
        (self.node_coords, self.mechBC_merged,  self.govNodes, self.govNodesMechBC, self.rigidPlates, self.trsprtBC_merged, self.vor, self.areas, self.functions, self.radii, self.rigidPlatesTrspt, self.govNodesTrspt, self.govNodesTrsptBC)  = utilitiesModeling.create3dBiparvaTubeTransport(self.cylinderRad, self.cylinderHeight, self.tubeThickness, self.minDist, self.trials, self.maxLim)
        self.measuringGauges = utilitiesModeling.assembleMeasuringGauges('cylinder3d', maxLim=self.maxLim)


    def run_3d_TubeInnerPressure(self):
        self.maxLim = np.array([self.cylinderHeight, self.cylinderRad, self.cylinderRad])
        (self.node_coords, self.mechBC_merged,  self.govNodes, self.govNodesMechBC, self.rigidPlates, self.trsprtBC_merged, self.vor, self.areas, self.functions, self.radii, self.rigidPlatesTrspt, self.govNodesTrspt, self.govNodesTrsptBC, self.interfaceVertexIndices)  = utilitiesModeling.create3dTubeInnerPressure(self.cylinderRad, self.cylinderHeight, self.tubeThickness, self.minDist, self.trials, self.maxLim, self.activeTransport, self.powerTes)
        self.measuringGauges = utilitiesModeling.assembleMeasuringGauges('cylinder3d', maxLim=self.maxLim)



    def run_2d_coupledArtificialCrack(self):
        (self.node_coords, self.mechBC_merged, self.trsprtBC_merged, self.notches, self.govNodes, self.govNodesMechBC, self.rigidPlates, self.vor, self.areas, self.functions)  = utilitiesModeling.create2dCoupledArtificialCrack(self.maxLim, self.minDist, self.trials, self.notchH)

    def run_2d_coupledPress(self):
        (self.node_coords, self.mechBC_merged, self.trsprtBC_merged,  self.govNodes, self.govNodesMechBC, self.rigidPlates, self.vor, self.areas, self.functions, self.rigidPlatesTrspt, self.govNodesTrspt, self.govNodesTrsptBC)  = utilitiesModeling.create2dCoupledPress(self.maxLim, self.minDist, self.trials)
        self.measuringGauges = utilitiesModeling.assembleMeasuringGauges('coupledPress2d', maxLim=self.maxLim)


    def run_3d_coupledArtificialCrack(self):
        (self.node_coords, self.mechBC_merged, self.trsprtBC_merged, self.govNodes, self.govNodesMechBC, self.rigidPlates, self.vor, self.areas, self.functions, self.rigidPlatesTrspt, self.govNodesTrspt, self.govNodesTrsptBC)  = utilitiesModeling.create3dCoupledArtificialCrack(self.maxLim, self.minDist, self.trials, self.notchH)



    def run_3d_coupledBrazilianDisc(self):
        (self.node_coords, self.mechBC_merged, self.trsprtBC_merged, self.govNodes, self.govNodesMechBC, self.rigidPlates, self.vor, self.areas, self.functions, self.rigidPlatesTrspt, self.govNodesTrspt, self.govNodesTrsptBC, self.radii)  = utilitiesModeling.createCoupledBrazilianDisc(np.zeros(3), self.cylinderRad, self.cylinderHeight,  self.minDist, self.trials, powerTes = self.powerTes)
        self.maxLim = np.array([self.cylinderHeight, 2*self.cylinderRad, 2*self.cylinderRad])
        self.measuringGauges = utilitiesModeling.assembleMeasuringGauges('3d_brazilianDisc', maxLim=self.maxLim)

    def run_2d_corrosionRebar(self, node_coords_init=None):
        (self.node_coords, self.mechBC_merged, self.trsprtBC_merged, self.govNodes, self.govNodesMechBC, self.rigidPlates, self.vor, self.areas, self.functions, self.rigidPlatesTrspt, self.govNodesTrspt, self.govNodesTrsptBC, self.interfaceVertexIndices)  = utilitiesModeling.create2dCorrosionRebar(self.maxLim, self.minDist, self.trials, self.rebarMinDist, self.rebarDiameter, self.rebarCount, self.rebarDepth, node_coords_init=node_coords_init, interfaceMinDist=self.interfaceMinDist, roughMinDistCoef=self.roughMinDistCoef, adaptivityReady=self.adaptivityReady,
        fineRingThickness=self.fineRingThickness, fineRegDepth=self.fineRegDepth, gradientRegDepth=self.gradientRegDepth)
        self.materialZones= utilitiesModeling.assembleMaterialZones(0,  2, model='2d_corrosionRebar', maxLim=self.maxLim, rebarDepth=self.rebarDepth, rebarDiameter=self.rebarDiameter, rebarCount=self.rebarCount)
        if self.activeTransport == False:
            self.rigidPlatesTrspt = []
            self.govNodesTrspt=[]
            self.govNodesTrsptBC = []

        expansionRingsProps = []
        expansionRingsProps.append(self.rebarCount)
        expansionRingsProps.append(self.rebarDepth)
        expansionRingsProps.append(self.rebarDiameter)
        expansionRingsProps.append(self.maxLim)

        self.measuringGauges = utilitiesModeling.assembleMeasuringGauges('2d_corrosionRebar', expansionRingsProps=expansionRingsProps)

    def run_2d_CFRAC_Clover(self, node_coords_init=None):
        self.activeTransport = False

        (self.node_coords, self.mechBC_merged, self.trsprtBC_merged, self.govNodes, self.govNodesMechBC, self.rigidPlates, self.vor, self.areas, self.functions)  = utilitiesModeling.create2d_CFRAC_Clover(self.maxLim, self.minDist, self.trials, self.holeMinDist, self.holeDiameter, -1, 1)

        self.materialZones= []
        if self.activeTransport == False:
            self.rigidPlatesTrspt = []
            self.govNodesTrspt=[]
            self.govNodesTrsptBC = []

        self.measuringGauges = utilitiesModeling.assembleMeasuringGauges('2d_CFRAC_Clover')

    def run_3d_CFRAC_Clover(self, node_coords_init=None):
        self.activeTransport = False

        (self.node_coords, self.mechBC_merged, self.trsprtBC_merged, self.govNodes, self.govNodesMechBC, self.rigidPlates, self.vor, self.areas, self.functions)  = utilitiesModeling.create3d_CFRAC_Clover(self.maxLim, self.minDist, self.trials, self.holeMinDist, self.holeDiameter, -1, 1)

        self.materialZones= []
        if self.activeTransport == False:
            self.rigidPlatesTrspt = []
            self.govNodesTrspt=[]
            self.govNodesTrsptBC = []

        self.measuringGauges = utilitiesModeling.assembleMeasuringGauges('3d_CFRAC_Clover')

    def run_2d_CFRAC_TDCB(self, node_coords_init=None):
        self.activeTransport = False
        self.fineWidth *= self.minDist
        elazonewidth = self.fineWidth/2-self.minDist

        (self.node_coords, self.mechBC_merged, self.trsprtBC_merged, self.govNodes, self.govNodesMechBC, self.rigidPlates, self.vor, self.areas, self.functions,self.notches,self.node_indices_dogbone)  = utilitiesModeling.create2d_CFRAC_TDCB(self.maxLim, self.minDist, self.trials, self.holeMinDist, self.holeDiameter, -1, self.roughMinDistCoef,elazonewidth=self.fineWidth/2)

        self.materialZones= utilitiesModeling.assembleMaterialZones (self.maxLim[0]/2-elazonewidth, 2, model='tdcb', maxLim=self.maxLim)
        if self.activeTransport == False:
            self.rigidPlatesTrspt = []
            self.govNodesTrspt=[]
            self.govNodesTrsptBC = []

        self.measuringGauges = utilitiesModeling.assembleMeasuringGauges('2d_CFRAC_TDCB')

    def run_2d_Hanging_FracZone(self, node_coords_init=None):
        self.activeTransport = False

        (self.node_coords, self.mechBC_merged, self.trsprtBC_merged, self.govNodes, self.govNodesMechBC, self.rigidPlates, self.vor, self.areas, self.functions)  = utilitiesModeling.create2d_Hanging_FracZone(self.maxLim, self.minDist, self.trials)

        coords = np.asarray(self.node_coords)
        lims = [np.amin(coords[:,0]),np.amax(coords[:,0]),np.amin(coords[:,1]),np.amax(coords[:,1])]

        self.materialZones= utilitiesModeling.assembleMaterialZones (0, 2, model='hangingfraczone', maxLim=lims)

        if self.activeTransport == False:
            self.rigidPlatesTrspt = []
            self.govNodesTrspt=[]
            self.govNodesTrsptBC = []

        self.measuringGauges = utilitiesModeling.assembleMeasuringGauges('hangingfraczone')

    def run_3d_Hanging_FracZone(self, node_coords_init=None):
        self.activeTransport = False

        (self.node_coords, self.mechBC_merged, self.trsprtBC_merged, self.govNodes, self.govNodesMechBC, self.rigidPlates, self.vor, self.areas, self.functions)  = utilitiesModeling.create3d_Hanging_FracZone(self.maxLim, self.minDist, self.trials)

        coords = np.asarray(self.node_coords)
        lims = [np.amin(coords[:,0]),np.amax(coords[:,0]),np.amin(coords[:,1]),
        np.amax(coords[:,1]),np.amin(coords[:,2]),np.amax(coords[:,2]) ]

        self.materialZones= utilitiesModeling.assembleMaterialZones (0, 3, model='hangingfraczone', maxLim=lims)

        if self.activeTransport == False:
            self.rigidPlatesTrspt = []
            self.govNodesTrspt=[]
            self.govNodesTrsptBC = []

        self.measuringGauges = utilitiesModeling.assembleMeasuringGauges('hangingfraczone')

    def run_3d_CFRAC_TDCB(self, node_coords_init=None):
        self.activeTransport = False
        self.fineWidth *= self.minDist
        elazonewidth = self.fineWidth /2

        (self.node_coords, self.mechBC_merged, self.trsprtBC_merged, self.govNodes, self.govNodesMechBC, self.rigidPlates, self.vor, self.areas, self.functions,self.notches,self.node_indices_dogbone)  = utilitiesModeling.create3d_CFRAC_TDCB(self.maxLim, self.minDist, self.trials, self.holeMinDist, self.holeDiameter, -1, self.roughMinDistCoef,elazonewidth=self.fineWidth/2,notchWidth=self.notchWidth,fracZoneHeight=self.fracZoneHeight,fracZoneOverhang=self.fracZoneOverhang, fem=self.fem)


        if self.fem == False:
            self.materialZones= utilitiesModeling.assembleMaterialZones (elazonewidth, 3, model='tdcb', maxLim=self.maxLim)
        else:
            coords = np.asarray(self.node_coords)
            lims = [np.amin(coords[:,0]),np.amax(coords[:,0]),np.amin(coords[:,1]),np.amax(coords[:,1]),np.amin(coords[:,2]),np.amax(coords[:,2]) ]
            self.materialZones= utilitiesModeling.assembleMaterialZones (elazonewidth, 3, model='tdcbfem', maxLim=lims)

        if self.activeTransport == False:
            self.rigidPlatesTrspt = []
            self.govNodesTrspt=[]
            self.govNodesTrsptBC = []

        self.measuringGauges = utilitiesModeling.assembleMeasuringGauges('3d_CFRAC_TDCB')

    def run_3d_corrosionRebar(self, node_coords_init=None):
        (self.node_coords, self.mechBC_merged, self.trsprtBC_merged, self.govNodes, self.govNodesMechBC, self.rigidPlates, self.vor, self.areas, self.functions, self.rigidPlatesTrspt, self.govNodesTrspt, self.govNodesTrsptBC, self.interfaceVertexIndices)  = utilitiesModeling.create3dCorrosionRebar(self.maxLim, self.minDist, self.trials, self.rebarMinDist, self.rebarDiameter, self.rebarCount, self.rebarDepth, node_coords_init=node_coords_init, interfaceMinDist=self.interfaceMinDist, roughMinDistCoef=self.roughMinDistCoef, adaptivityReady=self.adaptivityReady,
        fineRingThickness=self.fineRingThickness, fineRegDepth=self.fineRegDepth, gradientRegDepth=self.gradientRegDepth)
        self.materialZones= utilitiesModeling.assembleMaterialZones(0,  2, model='3d_corrosionRebar', maxLim=self.maxLim, rebarDepth=self.rebarDepth, rebarDiameter=self.rebarDiameter, rebarCount=self.rebarCount)
        if self.activeTransport == False:
            self.rigidPlatesTrspt = []
            self.govNodesTrspt=[]
            self.govNodesTrsptBC = []

        expansionRingsProps = []
        expansionRingsProps.append(self.rebarCount)
        expansionRingsProps.append(self.rebarDepth)
        expansionRingsProps.append(self.rebarDiameter)
        expansionRingsProps.append(self.maxLim)

        self.measuringGauges = utilitiesModeling.assembleMeasuringGauges('3d_corrosionRebar', expansionRingsProps=expansionRingsProps)

    def run_2d_coupledRVE(self):
        print('2d_coupledRVE')
        (self.node_coords, self.mechBC_merged, self.mechIC_merged, self.trsprtBC_merged, self.trsprtIC_merged, self.vor, self.areas, self.functions, self.radii)   = utilitiesModeling.create2dCoupledRVE(self.maxLim, self.minDist, self.trials, self.powerTes )
        self.materialZones=None
        self.periodicModel = 1

    def run_2d_circRVE(self):
        (self.node_coords, self.node_coords_polar, self.mechBC_merged, self.mechIC_merged, self.trsprtBC_merged, self.trsprtIC_merged, self.vor, self.areas, self.functions, self.radii)   = utilitiesModeling.create2dCircRVE(self.maxLim, self.minDist, self.trials, self.powerTes )
        self.materialZones = None
        self.periodicModel = 2

    def run_2d_cantileverBending(self):
        print('2d_cantilever bending')
        (self.node_coords, self.mechBC_merged, self.mechIC_merged, self.trsprtBC_merged, self.trsprtIC_merged, self.vor, self.areas, self.functions)   = utilitiesModeling.create2dCantileverBending(self.maxLim, self.minDist, self.trials )
        self.materialZones=None

    def run_3d_dam(self):
        self.materialZones = None
        (self.node_coords, self.mechBC_merged, self.mechIC_merged, self.trsprtBC_merged, self.trsprtIC_merged, self.vor, self.areas, self.functions)   = utilitiesModeling.create3dDam(self.maxLim, self.minDist, self.trials, self.Xtopsize)

    def saveGeometry(self):
        print('SAVING GEOMETRY')
        tube = False
        if self.modelType == '3d_BiparvaTubeTransport' or self.modelType=='3d_TubeInnerPressure': tube = True

        actvTrsprt = self.activeTransport
        if self.ansysOut == True:
            actvTrsprt = True
        #print('Extracting geometry...', end='')
        #if (self.printout == False): blockPrint()
        self.node_coords = np.asarray(self.node_coords)
        print('nodes before %d' %len(self.node_coords))
        (new_node_coords,self.vert_count,
        self.verticesIdxDict,
        self.vertIdxStart,
        self.totalNodeCount) = utilitiesGeom.extractGeometry(
            self.master_folder,
            self.dimension,
            len(self.node_coords),
            self.maxLim, self.vor,
            self.node_coords,
            self.node_coords_polar,
            self.areas,
            actvTrsprt, self.activeMechanics,
            mZ=self.materialZones, periodicModel=self.periodicModel,
            notches=self.notches, isTube=tube, coupled=self.coupled, minDist=self.minDist, node_indices_dogbone=self.node_indices_dogbone, randomizeMaterial=self.randomizeMaterial, auxmechelements=self.auxmechelements)

        print('NODES NODES after %d' %len(new_node_coords))

        self.node_coords = np.copy(new_node_coords)
        self.node_count = len(self.node_coords)
        #if (self.printout == False): enablePrint()

        if len(self.interfaceVertexIndices)>0:
            for i in range (len(self.interfaceVertexIndices)):
                for j in range (len(self.interfaceVertexIndices[i])):
                    self.interfaceVertexIndices[i][j] += int(self.vertIdxStart)

        print ('saving dome mdome.')

    def saveRest(self, solver, master_file, exporters, generateFineNodes=-1):
        print('saving rest')
        # NOTE JK: folder and bc_file already exist, then it is only appended, which results in error while bc are loaded to solver (two bc applied on the same dof). This cannot be done in saveMechBC, because it can be used by save constraints
        bc_path = os.path.join(self.master_folder,
                               utilitiesGeom.mechBCFile)
        if os.path.isfile(bc_path):
            os.remove(bc_path)
        print('Saving files...', end='')
        utilitiesGeom.saveMaterials(self.master_folder, self.materials)
        utilitiesGeom.saveFunctions(self.master_folder, self.functions)
        if self.activeMechanics:
            utilitiesGeom.saveMechBC(self.master_folder, self.dimension, self.mechBC_merged)
            if (self.mechIC_merged != None and len(self.mechIC_merged)>0):
                utilitiesGeom.saveMechIC(self.master_folder, self.dimension, self.mechIC_merged)
        if self.activeTransport:
            utilitiesGeom.saveTransportBC(self.master_folder, self.trsprtBC_merged, self.verticesIdxDict, self.vertIdxStart)
            if (self.trsprtIC_merged != None and len(self.trsprtIC_merged)>0):
                utilitiesGeom.saveTransportIC(self.master_folder, self.trsprtIC_merged)

        utilitiesGeom.saveExporters(self.master_folder, self.activeTransport, self.activeMechanics, exporters=exporters)



        if self.govNodes != None:
            if self.rigidPlates != None:
                if self.govNodesMechBC != None:
                    if self.modelType == '2d_corrosionRebar' or self.modelType == '3d_corrosionRebar' or self.modelType == '3d_TubeInnerPressure':
                        if self.modelType == '3d_TubeInnerPressure':

                            expansionRingsProps = []
                            expansionRingsProps.append(1)
                            expansionRingsProps.append(0)
                            expansionRingsProps.append(self.cylinderRad*2-self.tubeThickness*2)
                            expansionRingsProps.append(self.maxLim)
                            expansionRingsProps.append(0) #dir
                        if self.modelType == '2d_corrosionRebar' or self.modelType == '3d_corrosionRebar':
                            expansionRingsProps = []
                            expansionRingsProps.append(self.rebarCount)
                            expansionRingsProps.append(self.rebarDepth)
                            expansionRingsProps.append(self.rebarDiameter)
                            expansionRingsProps.append(self.maxLim)
                            expansionRingsProps.append(2) #dir

                        mechDofIndices = utilitiesGeom.saveConstraint(self.master_folder, self.dimension, self.govNodes, self.govNodesMechBC, self.rigidPlates, self.totalNodeCount, self.node_coords, expansionRingsProps=expansionRingsProps, virtualDoF=expansionRingsProps[0], nodesMechBC=self.mechBC_merged)
                        self.constraint = True
                    else:
                        utilitiesGeom.saveConstraint(self.master_folder, self.dimension, self.govNodes, self.govNodesMechBC, self.rigidPlates, self.totalNodeCount, self.node_coords,nodesMechBC=self.mechBC_merged)
                        self.constraint = True

        self.totalNodeCount += len(self.govNodes)

        if self.modelType == '2d_corrosionRebar' or self.modelType == '3d_corrosionRebar'  :
            self.totalNodeCount += self.rebarCount

        if self.modelType == '3d_TubeInnerPressure':
            self.totalNodeCount += 1

        #saving transport rigid plates
        if self.govNodesTrspt != None:
            if self.rigidPlatesTrspt != None:
                if self.govNodesTrsptBC != None:
                    utilitiesGeom.saveConstraintTransport(self.master_folder, self.dimension, self.govNodesTrspt, self.govNodesTrsptBC, self.rigidPlatesTrspt, self.totalNodeCount, self.node_coords, self.vert_count, self.verticesIdxDict, self.vertIdxStart)
                    self.constraintTrspt = True
                    #
                    if len(self.interfaceVertexIndices)>0:
                        #print(self.interfaceVertexIndices)

                        if (self.modelType == '2d_corrosionRebar' or self.modelType == '3d_corrosionRebar' ):
                            direction = 0
                            materialId = 1

                            circleLength = 2*np.pi*self.rebarDiameter/2
                            nrNodes = int ( circleLength / self.interfaceMinDist )
                            nrNodes = (int (nrNodes / 4) +1 ) * 4

                            interfaceElementLength = circleLength / nrNodes
                            rebarArea = np.pi * (self.rebarDiameter/2)**2
                            if self.dimension ==3:
                                rebarArea *=self.maxLim[2]

                            coeff = interfaceElementLength / rebarArea

                        if (self.modelType == '3d_TubeInnerPressure'):
                            direction = 0
                            materialId = 1


                            circleLength = 2*np.pi*(self.cylinderRad-self.tubeThickness)

                            nrNodes = int ( circleLength / self.minDist )
                            nrNodes = (int (nrNodes / 4) +1 ) * 4

                            rebarArea = np.pi * (self.cylinderRad)**2 * self.cylinderHeight

                            interfaceElementLength = circleLength / nrNodes
                            coeff = interfaceElementLength / rebarArea

                        utilitiesGeom.saveCoupledConstraint(self.master_folder, self.interfaceVertexIndices, mechDofIndices, direction, materialId, coeff)


        if (self.measuringGauges!=None):
            utilitiesGeom.saveMeasuringGauges(self.master_folder, self.dimension, self.measuringGauges)

        if len(self.radii) > 0:
            if self.powerTes:
                utilitiesGeom.saveRadii(self.master_folder, self.radii)

        utilitiesGeom.saveMasterInput(self.master_folder, self.dimension, solver.solverType, solver.time_step, solver.min_time_step, solver.max_time_step, solver.total_time, self.activeTransport, self.activeMechanics, periodic=self.periodicModel, constraint=self.constraint, constraintTrspt=self.constraintTrspt,
        limitTolerance= solver.limit_tolerance, maxIt=solver.maxIt, tolerance=solver.tolerance, auxMechElements=self.auxmechelements, masterSolver=self.masterSolver, masterMaterials=self.masterMaterials, masterFunctions=self.masterFunctions)

        # if src and dest are same, copyfile raises SameFileError Exception https://docs.python.org/3/library/shutil.html#shutil.SameFileError
        # get only the filename from master file string https://docs.python.org/3/library/os.path.html#os.path.basename

        dst_file = os.path.join(self.master_folder, os.path.basename(master_file))
        if not os.path.isfile(dst_file):
            print ('Copying prep_master used...', end='')
            copyfile(master_file, dst_file)


        print('master %s' %self.master_folder)
        print('cwd %s' %os.getcwd())

        df = self.defaultFilesFolder
        print('df %s' %df)

        if df != None:
            files=os.listdir(df)
            for fname in files:
                print(os.path.join(self.defaultFilesFolder ,fname))
                print(os.path.isfile(os.path.join(self.defaultFilesFolder ,fname)))
                if os.path.isfile(os.path.join(self.defaultFilesFolder ,fname)):
                   shutil.copy2(os.path.join(self.defaultFilesFolder ,fname), self.master_folder)


        """
        if generateFineNodes:
            fine_nodes, mechBC_merged, mechInitC_merged, notches, govNodes, govNodesMechBC, rigidPlates  = utilitiesModeling.assemble3DSSBeamBending(self.maxLim, generateFineNodes, self.trials, self.notchH, self.loadWidth, fracZoneWidth=self.fracZoneWidth, orthogonalFracZone=False, notchWidth = self.notchWidth, coupled=False, node_coords_init=None, specifiedNodes=self.specifiedNodes, roughMinDistCoef=1, supportDivision=self.supportDivision);
            fine_nodes = np.asarray(fine_nodes)
            fine_nodes = np.hstack(( fine_nodes, np.zeros((len(fine_nodes),1))   ))
            utilitiesGeom.saveNodes(self.master_folder, fine_nodes, "Particle",3, 'nodesFine.inp')
        """
        print ('done.')

        utilitiesGeom.checkSavedModel(self.master_folder, self.dimension, self.activeMechanics, self.activeTransport)

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
                sys.exit(1)

            linElMaterial = utilitiesMech.linearElasticMaterial(young, alpha, density)
            self.materials.append(linElMaterial)
            print('done.')

        #TransportMaterial
        if (r[1]=='TrsprtMaterial'):
            #
            viscosity = None
            permeability = None
            density = None
            capacity = None
            crack_turtuosity = 1.
            biot_coeff = 0.
            #
            for i in range (len(r)):
                if (r[i]=='capacity'):
                    capacity = float(r[i+1])
                if (r[i]=='density'):
                    density = float(r[i+1])
                if (r[i]=='permeability'):
                    permeability = float(r[i+1])
                if (r[i]=='viscosity'):
                    viscosity = float(r[i+1])
                if (r[i]=='crack_turtuosity'):
                    crack_turtuosity = float(r[i+1])
                if (r[i]=='biot_coeff'):
                    biot_coeff = float(r[i+1])

            #
            if (viscosity == None or permeability == None or density == None or capacity == None):# or crack_turtuosity==None):
                print ('!! TrsprtMaterial incomplete. Exiting. !!')
                sys.exit()

            transportMaterial = utilitiesMech.TransportMaterial(viscosity, permeability, density, capacity, crack_turtuosity, biot_coeff=biot_coeff, coupled=self.coupled)
            self.materials.append(transportMaterial)
            print('done.')

        #CSLMaterial
        if (r[1]=='CSLMaterial'):

            young = None
            alpha = None
            density = None
            ft = None
            Gt = None
            #"""
            """
            params_list=['young', 'alpha', 'density', 'ft', 'Gt']
            params = {param: None for param in params_list}
            for i in range (len(r)):
                if r[i] in params:
                    setattr(params, r[i], bool(r[i+1]))

            if None in params.values():
                print ('!! CSLMaterial incomplete. Exiting. !!')
                exit()

            #"""
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
            #"""
            if (young == None or alpha == None or density == None or ft == None or  Gt == None):
                print ('!! CSLMaterial incomplete. Exiting. !!')
                sys.exit()

            CSLMaterial = utilitiesMech.CSLMaterial(young, alpha, density, ft, Gt, coupled = self.coupled)
            self.materials.append(CSLMaterial)
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
        self.tolerance = 1e-3

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
            if (r[i]=='tolerance'):
                self.tolerance = float(r[i+1])

        print('done.')
        #return solverType, time_step, max_time_step, min_time_step, total_time, limit_tolerance, maxIt





















if __name__ == '__main__':
    print('\n%%%%%%%%% LATTICE PREPROCESSOR STARTED %%%%%%%%%')
    start = time.time()

    if len(sys.argv) > 1:
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

    model = None
    solver = None
    exporters=[]

    f = open (file, 'r')
    for row in f:
        if row and row.strip() and not row.startswith('#'):
            r = row.split()
            #print(r)
            #
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

    if len(model.specifiedNodes)>0:
        print ("Specified Nodes:")
        print (model.specifiedNodes)
    if model == None:
        print ('Missing model!! Exiting...')
        sys.exit(1)

    if solver == None:
        print ('Missing solver!! Exiting...')
        sys.exit(1)

    if len(model.materials)==0:
        print ('Missing some material!! Exiting...')
        #sys.exit(1)

    if model != None:
        for i in range (model.nr_models):

            if (model.modelType=='2d_dogboneBand') and i==0:
                model.maxLim [0]= float(model.maxLim [0]*model.minDist)
                model.maxLim [1]= float(model.maxLim [1]*model.minDist*(2*model.elasticHeightCoef+1))



            print('\nCreating model #%d' %i)
            if len(sys.argv) > 3:
                model.userSeed = int(sys.argv[3])
            model.setDirectory(len(sys.argv) > 2 and sys.argv[2] or None) # directory to generate in can be specified in the input string (after prep_master.inp)
            dst_file = os.path.join(model.master_folder, os.path.basename(file))
            #print ('prep_master used... %s' %dst_file, end='')
            if not os.path.isfile(dst_file):
                print ('Copying prep_master used... %s' %dst_file)
                copyfile(file, dst_file)

            model.createModel(modelnr = i)
            model.saveGeometry()
            model.saveRest(solver, file, exporters, generateFineNodes=0.02)

    print('\n%%%%%%%%% LATTICE PREPROCESSOR DONE %%%%%%%%%')
    #print('\n%%%%%%%%% %d NODES MODEL %%%%%%%%%' %len(model.node_coords))
    print('All done in %f seconds. (%d minutes).' %((time.time()-start), (time.time()-start)/60.0))




























#
