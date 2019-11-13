from traits.api import *
from traitsui.api import *
from traitsui.menu import *
import ModelGenerators
from MaterialTraits import *

class MechMaterialSelector(HasTraits):
    material = Enum ( mechMaterials_keys )
    selected = Property( Instance (MechMaterial), depends_on = 'material' )

    @cached_property
    def _get_selected(self):
        return mechMaterials_dict[self.material]

    trait_view = View(
                    UItem('material', style='simple'),
                    UItem('selected', style='custom' )
                    )


class TrsprtMaterialSelector(HasTraits):
    material = Enum ( trsprtMaterials_keys )
    selected = Property( Instance (TrsprtMaterial), depends_on = 'material' )

    @cached_property
    def _get_selected(self):
        return trsprtMaterials_dict[self.material]

    trait_view = View(
                    UItem('material', style='simple'),
                    UItem('selected', style='custom' )
                    )

# Model abstract class
class ModelType(HasTraits):
    name = Str
    dimension = Int
    trials = Int
    min_dist = Float
    runButton = Button('Generate model!')


class ModelType1(ModelType):
    name = '2D: Cantilever bending'
    dimension = Int(2)
    trials = Int(50000)
    xDim = Float(1.0)
    yDim = Float(0.2)
    min_dist = Float(0.03)
    yDispl = Float(0.001)

    mechMaterialSelector = MechMaterialSelector()
    mechMaterial = mechMaterialSelector.selected

    trsprtMaterialSelector = TrsprtMaterialSelector()
    trsprtMaterial = trsprtMaterialSelector.selected

    traits_view = View(Group(
        HSplit(
                (
                    Item('trials', label = 'Trials'),
                    Item('xDim'),
                    Item('yDim'),
                    Item('min_dist'),
                    Item('yDispl'),
                    UItem('runButton')
                ),
                (
                    Item('mechMaterial', style='simple', label = 'Mech Material')
                ),
                (
                    Item('trsprtMaterialSelector', style='simple', label = 'Trsprt Material')
                )
            )
        )
    )

    def _runButton_fired(self):
        print ('\nGenerating %s ...' % self.name)
        ModelGenerators.GenerateModelType1 (self.trials, self.xDim, self.yDim, self.min_dist, self.yDispl, self.mechMaterial.get_material_object(), self.trsprtMaterial.get_material_object())


###############################################
# 3D models
class ModelType37(ModelType):
    name = '3D: Tube cyclic torsion, displacement load'
    dimension = Int(3)
    trials = Int(50000)
    length = Float(1.0)
    outerRadius = Float(0.25)
    thickness = Float(0.1)
    min_dist = Float(0.05)
    rotAng = Float(0.005)

    elasticZones = Bool(True)
    withoutTransport = Bool(True)

    mechMaterialSelector = Instance ( MechMaterialSelector )
    mechMaterialElaSelector = Instance ( MechMaterialSelector )
    trsprtMaterialSelector = Instance ( TrsprtMaterialSelector )

    mechMaterial = Instance ( MechMaterial )
    mechMaterialEla = Instance ( MechMaterial )
    trsprtMaterial = Instance ( TrsprtMaterial )

    mechGroup = Group(
        UItem('mechMaterialSelector', style='custom'),
        show_border = True, show_labels=True, label = 'Main mechanical material'
    )
    mechElaGroup = Group(
        UItem('mechMaterialElaSelector', style='custom'),
        show_border = True, show_labels=True, label = 'Ela zone mechanical material'
    )
    trsprtGroup = Group(
        UItem('trsprtMaterialSelector', style='custom'),
        show_border = True, show_labels=True, label = 'Transport material'
    )


    traits_view = View(
        Group(
                (Group(
                    Item('trials', label = 'Trials'),
                    Item('length'),
                    Item('outerRadius'),
                    Item('thickness'),
                    Item('min_dist'),
                    Item('rotAng'),
                    Item('elasticZones', label = 'ElaZones'),
                    Item('withoutTransport', label = 'NO Trsprt'),
                    UItem('runButton'),
                    show_border = True, show_labels=True, label = 'Model properties'
                    )
                ),
                ( mechGroup   ),
                ( mechElaGroup ),
                ( trsprtGroup ), orientation='horizontal', columns = 4
                )
            )

    def _runButton_fired(self):
        mechMaterial = self.mechMaterialSelector.selected
        mechMaterialEla = self.mechMaterialElaSelector.selected
        trsprtMaterial = self.trsprtMaterialSelector.selected

        print ('\nGenerating %s ...' % self.name)
        ModelGenerators.GenerateModelType37 (self.trials, self.length, self.outerRadius, self.thickness, self.min_dist, self.rotAng, self.elasticZones, self.withoutTransport, mechMaterial.get_material_object(), mechMaterialEla.get_material_object(), trsprtMaterial.get_material_object(), self.rotAng)
