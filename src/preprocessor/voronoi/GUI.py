from traits.api import *
from traitsui.api import *
from traitsui.menu import *

from ModelTraits import *
from MaterialTraits import *

models = [
    #2D: Cantilever bending
     #ModelType1(),
    #3D:Tube torsion
    ModelType37(mechMaterialSelector = MechMaterialSelector(), mechMaterialElaSelector = MechMaterialSelector(),  trsprtMaterialSelector = TrsprtMaterialSelector() )
]
model_dict = { m.name: m for m in models }
model_keys = list(model_dict.keys())


class ModelSelector(HasTraits):
    model = Enum ( model_keys )
    selected = Property( Instance (ModelType), depends_on = 'model' )
    @cached_property
    def _get_selected(self):
        return model_dict[self.model]

    trait_view = View(
                    Item('model', style='simple'),
                    UItem('selected', style='custom' ),
                    resizable = True,  title = 'Preprocessor GUI'
                    )


if __name__ == '__main__':
    print ('WORK IN PROGRESS GUI')
    modelSelector = ModelSelector()
    modelSelector.configure_traits()
