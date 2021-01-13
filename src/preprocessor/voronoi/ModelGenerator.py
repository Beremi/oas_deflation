
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

        self.node_coords = []

        self.mechBC_merged = self.mechIC_merged = self.trsprtBC_merged = self.trsprtIC_merged = self.govNodesMechBC = self.govNodesTrsptBC = None

        self.functions = []
        self.radii = []
        self.totalNodeCount = -1
        self.notches = []
        self.materialZones= []
        self.measuringGauges = []
        self.symmetric = False

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
        self.maxLim = np.zeros((self.dimension))

        #RWTH cylinderRad
        #cylinderHeight 0.1 cylinderRad 0.1 minDist 0.02 notchRadLeft 0.08 notchRadRight 0.05 notchWidth 0.01
        self.notchRadLeft = None
        self.notchRadRight = None
        self.innerRadTop = None
        self.innerRadBottom = None
        self.notchWidth = None
        self.notchDepth = None
        self.RWTHQuarter = False

        for i in range (len(r)):

            if (r[i]=='printout'):
                if r[i+1] == 1:
                    print('print')
                    self.printout = True

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
            if (r[i]=='coupled'):
                if (int(r[i+1])==1): self.activeMechanics = True
                if (int(r[i+1])==1): self.activeTransport = True
                if (int(r[i+1])==1): self.coupled = True
            if (r[i]=='periodicModel'):
                if (int(r[i+1])==1): self.periodicModel = True
                if (int(r[i+1])==0): self.periodicModel = False
            #"""
            """
            keys = ['powerTes', 'activeMechanics', 'activeTransport', 'periodicModel']
            if r[i] in keys:
                setattr(self, r[i], bool(r[i+1]))
            """
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
            if (r[i]=='notchH'):
                self.notchH = float(r[i+1])
            if (r[i]=='loadWidth'):
                self.loadWidth = float(r[i+1])
            if (r[i]=='fracZoneWidth'):
                self.fracZoneWidth = float(r[i+1])
            if (r[i]=='orthogonalFracZone'):
                if (int(r[i+1])==1): self.orthogonalFracZone = True
                if (int(r[i+1])==0): self.orthogonalFracZone = False

            if (r[i]=='symmetric'):
                if (int(r[i+1])==1): self.symmetric = True



        print('done.')

    def setDirectory(self, dirNam=None):
        if dirNam is None:
            if self.userSeed == -1:
                self.seed = np.random.randint(1000.0)
                np.random.seed(seed=self.seed)
            else:
                self.seed = self.userSeed
                np.random.seed(seed=self.seed)

            self.master_folder = 'power_%.4f_%02d' % (self.minDist, self.seed)
        else:
            self.master_folder = dirNam

        try:
            if not os.path.exists(self.master_folder):
                os.makedirs(self.master_folder)
        except:
            print('Please create directory %s! Code Exited.' % self.master_folder)
            sys.exit()


    def createModel(self, node_coords_init=None):
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
            self.run_3d_notched3pb(node_coords_init=node_coords_init)

        if self.modelType == '2d_dogbone':
            self.run_2d_dogbone()
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


        if self.modelType == '2d_singleSpring':
            self.run_2d_singleSpring()

        if self.modelType == '2d_coupledArtificialCrack':
            self.run_2d_coupledArtificialCrack()

        if self.modelType == '2d_coupledRVE':
            self.run_2d_coupledRVE()

        if self.modelType =='3d_coupledBrazilianDisc':
            self.run_3d_coupledBrazilianDisc()

        if self.modelType =='2d_cantileverBending':
            self.run_2d_cantileverBending()


        #if (self.printout == False): enablePrint()

    def run_2d_singleSpring(self):
        (self.node_coords, self.mechBC_merged,  self.vor, self.areas, self.functions,self.govNodes, self.govNodesMechBC, self.rigidPlates) = utilitiesModeling.createSingleSpringTestModel( self.minDist, self.master_folder )
        self.measuringGauges = utilitiesModeling.assembleMeasuringGauges('2d_singleSpring', maxLim=np.array([self.minDist, 1]) )
        materialZones = None

    def run_2d_notched3pb(self, node_coords_init=None):
        (self.node_coords, self.mechBC_merged, self.mechIC_merged, self.vor, self.areas, self.functions, self.notches, self.govNodes,
        self.govNodesMechBC, self.rigidPlates)  = utilitiesModeling.create2dSSBeamUnifLoad(self.maxLim, self.minDist, self.trials, notch=self.notchH, loadWidth=self.loadWidth, fracZoneWidth = self.fracZoneWidth, orthogonalFracZone=self.orthogonalFracZone, notchWidth=self.notchWidth, node_coords_init=node_coords_init)
        self.measuringGauges = utilitiesModeling.assembleMeasuringGauges('3pb2d', maxLim=self.maxLim)

    def run_3d_notched3pb(self, node_coords_init=None):
        (self.node_coords, self.mechBC_merged, self.mechIC_merged, self.vor, self.areas, self.functions, self.notches, self.govNodes, self.govNodesMechBC, self.rigidPlates, self.trsprtBC_merged, self.trsprtIC_merged) = utilitiesModeling.create3dSSBeamUnifLoad(self.maxLim, self.minDist, self.trials, notch=self.notchH, loadWidth=self.loadWidth, fracZoneWidth = self.fracZoneWidth, orthogonalFracZone=self.orthogonalFracZone, notchWidth=self.notchWidth, coupled=self.coupled, node_coords_init=node_coords_init)
        self.measuringGauges = utilitiesModeling.assembleMeasuringGauges('3pb3d', maxLim=self.maxLim)

        materialZones = None

    def run_2d_dogbone(self):
        (self.node_coords,self.mechBC_merged,self.mechIC_merged,self.trsprtBC_merged,self.trsprtIC_merged,self.vor,self.areas,self.functions,self.govNodes,self.govNodesMechBC,self.rigidPlates)   = utilitiesModeling.create2dDogBone(self.minDist, self.trials, D=self.dogboneD, excentricity=self.dogboneExcentricityFrac, symmetric=self.symmetric )
        self.materialZones=None
        self.measuringGauges = utilitiesModeling.assembleMeasuringGauges('dogbone2d', D=self.dogboneD)

    def run_3d_dogbone(self):
        self.maxLim = np.array([self.dogboneD, 6/4*self.dogboneD, 0.1])
        #self.materialZones = utilitiesModeling.assembleMaterialZones (self.minDist*2, 3, model='dogbone',  D=self.dogboneD, thickness=0.1)
        self.materialZones = None
        (self.node_coords, self.mechBC_merged, self.mechIC_merged, self.trsprtBC_merged, self.trsprtIC_merged, self.vor, self.areas, self.functions, self.govNodes, self.govNodesMechBC, self.rigidPlates)   = utilitiesModeling.create3dDogBone(self.minDist, self.trials, D=self.dogboneD, excentricity=self.dogboneExcentricityFrac )


    def run_3d_torsionPress(self):
        self.maxLim = np.array([self.cylinderHeight, 2*self.cylinderRad, 2*self.cylinderRad])
        #self.materialZones = utilitiesModeling.assembleMaterialZones (self.minDist*2, 3, model='box', maxLim=self.maxLim)
        (self.node_coords, self.mechBC_merged,  self.vor, self.areas, self.functions, self.govNodes, self.govNodesMechBC, self.rigidPlates)    = utilitiesModeling.create3dcylinderTorsionPressFree(np.zeros(3), self.cylinderRad, self.cylinderHeight,  self.minDist, self.trials, 0 )
        self.measuringGauges = utilitiesModeling.assembleMeasuringGauges('cylinder3d', maxLim=self.maxLim)
        self.materialZones =  None

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
        (self.node_coords, self.mechBC_merged, self.mechIC_merged, self.trsprtBC_merged, self.trsprtIC_merged, self.vor, self.areas, self.functions, self.nodePositions, self.coupledNodes, self.mirtype)   = utilitiesModeling.create2dPeriodicShear(self.maxLim, self.minDist, self.trials )
        self.materialZones=None
        self.periodicModel = 1

    def run_3d_periodicShear(self):
        (self.node_coords, self.mechBC_merged, self.mechIC_merged, self.trsprtBC_merged, self.trsprtIC_merged, self.vor, self.areas, self.functions)   = utilitiesModeling.create3dPeriodicShear(self.maxLim, self.minDist, self.trials )
        self.materialZones=None
        self.periodicModel = 1

    def run_2d_transportPatchTest(self):
        (self.node_coords, self.mechBC_merged, self.trsprtBC_merged, self.vor, self.areas, self.functions, self.radii)  = utilitiesModeling.createPatchTestTransport(self.maxLim, self.minDist, self.trials, self.dimension, self.powerTes)

    def run_3d_transportPatchTest(self):
        (self.node_coords, self.mechBC_merged, self.trsprtBC_merged, self.vor, self.areas, self.functions, self.radii)  = utilitiesModeling.createPatchTestTransport(self.maxLim, self.minDist, self.trials, self.dimension, self.powerTes)

    def run_3d_BiparvaTubeTransport(self):
        self.maxLim = np.array([self.cylinderHeight, self.cylinderRad, self.cylinderRad])
        (self.node_coords, self.mechBC_merged, self.trsprtBC_merged, self.vor, self.areas, self.functions, self.radii)  = utilitiesModeling.create3dBiparvaTubeTransport(self.cylinderRad, self.cylinderHeight, self.tubeThickness, self.minDist, self.trials, self.maxLim)

    def run_2d_coupledArtificialCrack(self):
        (self.node_coords, self.mechBC_merged, self.trsprtBC_merged, self.notches, self.govNodes, self.govNodesMechBC, self.rigidPlates, self.vor, self.areas, self.functions)  = utilitiesModeling.createCoupledArtificialCrack(self.maxLim, self.minDist, self.trials, self.notchH)


    def run_3d_coupledBrazilianDisc(self):
        (self.node_coords, self.mechBC_merged, self.trsprtBC_merged, self.govNodes, self.govNodesMechBC, self.rigidPlates, self.vor, self.areas, self.functions, self.rigidPlatesTrspt, self.govNodesTrspt, self.govNodesTrsptBC)  = utilitiesModeling.createCoupledBrazilianDisc(np.zeros(3), self.cylinderRad, self.cylinderHeight,  self.minDist, self.trials)
        self.maxLim = np.array([self.cylinderHeight, 2*self.cylinderRad, 2*self.cylinderRad])
        self.measuringGauges = utilitiesModeling.assembleMeasuringGauges('3d_brazilianDisc', maxLim=self.maxLim)


    def run_2d_coupledRVE(self):
        print('2d_coupledRVE')
        (self.node_coords, self.mechBC_merged, self.mechIC_merged, self.trsprtBC_merged, self.trsprtIC_merged, self.vor, self.areas, self.functions, self.nodePositions, self.coupledNodes, self.mirtype)   = utilitiesModeling.create2dCoupledRVE(self.maxLim, self.minDist, self.trials )
        self.materialZones=None
        self.periodicModel = 1

    def run_2d_cantileverBending(self):
        print('2d_cantilever bending')
        (self.node_coords, self.mechBC_merged, self.mechIC_merged, self.trsprtBC_merged, self.trsprtIC_merged, self.vor, self.areas, self.functions)   = utilitiesModeling.create2dCantileverBending(self.maxLim, self.minDist, self.trials )
        self.materialZones=None


    def saveGeometry(self):
        tube = False
        if self.modelType == '3d_BiparvaTubeTransport': tube = True

        #print('Extracting geometry...', end='')
        #if (self.printout == False): blockPrint()
        self.node_coords = np.asarray(self.node_coords)
        (self.vert_count,
        self.verticesIdxDict,
        self.vertIdxStart,
        self.totalNodeCount) = utilitiesGeom.extractGeometry(
            self.master_folder,
            self.dimension,
            len(self.node_coords),
            self.maxLim, self.vor,
            self.node_coords,
            self.areas,
            self.activeTransport, self.activeMechanics,
            mZ=self.materialZones, periodicModel=self.periodicModel,
            nodePositions=self.nodePositions, coupledNodes=self.coupledNodes,
            mirtype=self.mirtype, notches=self.notches, isTube=tube, coupled=self.coupled)

        #if (self.printout == False): enablePrint()
        print ('done.')

    def saveRest(self, solver, master_file):
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

        utilitiesGeom.saveExporters(self.master_folder, self.activeTransport, self.activeMechanics)

        if self.govNodes != None:
            if self.rigidPlates != None:
                if self.govNodesMechBC != None:
                    if not (self.activeTransport==True and self.coupled==False):
                        utilitiesGeom.saveConstraint(self.master_folder, self.dimension, self.govNodes, self.govNodesMechBC, self.rigidPlates, self.totalNodeCount, self.node_coords)
                        self.constraint = True

        self.totalNodeCount += len(self.govNodes)


        #saving transport rigid plates
        if self.govNodesTrspt != None:
            if self.rigidPlatesTrspt != None:
                if self.govNodesTrsptBC != None:
                    utilitiesGeom.saveConstraintTransport(self.master_folder, self.dimension, self.govNodesTrspt, self.govNodesTrsptBC, self.rigidPlatesTrspt, self.totalNodeCount, self.node_coords, self.vert_count, self.verticesIdxDict, self.vertIdxStart)
                    self.constraintTrspt = True



        if (self.measuringGauges!=None):
            utilitiesGeom.saveMeasuringGauges(self.master_folder, self.dimension, self.measuringGauges)

        if len(self.radii) > 0:
            if self.powerTes:
                utilitiesGeom.saveRadii(self.master_folder, self.radii)

        utilitiesGeom.saveMasterInput(self.master_folder, self.dimension, solver.solverType, solver.time_step, solver.min_time_step, solver.max_time_step, solver.total_time, self.activeTransport, self.activeMechanics, periodic=self.periodicModel, constraint=self.constraint, constraintTrspt=self.constraintTrspt,
        limitTolerance= solver.limit_tolerance, maxIt=solver.maxIt, tolerance=solver.tolerance)

        # if src and dest are same, copyfile raises SameFileError Exception https://docs.python.org/3/library/shutil.html#shutil.SameFileError
        dst_file = os.path.join(self.master_folder,master_file)
        if not os.path.isfile(dst_file):
            print ('Copying prep_master used...', end='')
            copyfile(master_file, dst_file)

        print ('done.')

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
            viscosity = None
            permeability = None
            density = None
            capacity = None
            crack_turtuosity = 1.
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
            #
            if (viscosity == None or permeability == None or density == None or capacity == None):# or crack_turtuosity==None):
                print ('!! TrsprtMaterial incomplete. Exiting. !!')
                sys.exit()

            transportMaterial = utilitiesMech.TransportMaterial(viscosity, permeability, density, capacity, crack_turtuosity, coupled=self.coupled)
            self.materials.append(transportMaterial)
            print('done.')

        #MarsMaterial
        if (r[1]=='MarsMaterial'):

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
                print ('!! MarsMaterial incomplete. Exiting. !!')
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

    model = None
    solver = None

    f = open (file, 'r')
    for row in f:
        if not (row[0]=='#'):
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

    if model == None:
        print ('Missing model!! Exiting...')
        sys.exit()

    if solver == None:
        print ('Missing solver!! Exiting...')
        sys.exit()

    if len(model.materials)==0:
        print ('Missing some material!! Exiting...')
        sys.exit()

    if model != None:
        for i in range (model.nr_models):
            print('\nCreating model #%d' %i)
            model.setDirectory()
            model.createModel()
            model.saveGeometry()
            model.saveRest(solver, file)




    print('\n%%%%%%%%% LATTICE PREPROCESSOR DONE %%%%%%%%%')
    #print('\n%%%%%%%%% %d NODES MODEL %%%%%%%%%' %len(model.node_coords))
    print('All done in %f seconds. (%d minutes).' %((time.time()-start), (time.time()-start)/60.0))




























#
