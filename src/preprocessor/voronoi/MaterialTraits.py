from traits.api import *
from traitsui.api import *
#from traitsui.menu import *

import utilitiesMech


class Material(HasTraits):
    name = Str

class MechMaterial(Material):
    pass

class TrsprtMaterial(Material):
    pass

class LinearElasticMaterial(MechMaterial):
    name = 'Linear Elastic Material'
    young_modulus = Float(30e9)
    poisson = Float(0.3)
    density = Float(2200)

    traits_view = View(
        Item('young_modulus', label='E [Pa]'),
        Item('poisson', label = 'Poisson [-]'),
        Item('density', label = 'Density [kg/m3]')
    )

    def get_material_object(self):
        return utilitiesMech.linearElasticMaterial(self.young_modulus, self.poisson, self.density)


class MarsMaterial(MechMaterial):
    name = 'Mars Material'
    young_modulus = Float(30e9)
    poisson = Float(0.3)
    density = Float(2200)
    ft = Float(2e6)
    Gt = Float(500)

    traits_view = View(
        Item('young_modulus', label='E [Pa]'),
        Item('poisson', label = 'Poisson [-]'),
        Item('density', label = 'Density [kg/m3]'),
        Item('ft', label = 'ft'),
        Item('Gt', label = 'Gt')
    )

    def get_material_object(self):
        return utilitiesMech.MarsMaterial(self.young_modulus, self.poisson, self.density, self.ft, self.Gt)


#E0	35e9	alpha	0.300000    density 2200.0 fc 200e6 ft 35e6 KinN 4e9 gammaN 20e9 m -0.2e-6 Ad 4000e-6 tauBar 4.0e6 Kin 0.0 gamma 10.0e6 S 0.00025e6 a 0
class FatigueMaterial(MechMaterial):
    name = 'Fatigue Material'
    young_modulus = Float(35e9)
    poisson = Float(0.3)
    density = Float(2200)
    fc = Float(200e6)
    ft = Float(35e6)
    KinN = Float(4e9)
    gammaN = Float(20e9)
    m = Float(-0.2e-6)
    Ad = Float(4000e-6)
    tau_bar = Float(4e6)
    Kin = Float(0.0)
    gamma = Float(10e6)
    S = Float(0.00025e6)
    a = Float(0)

    traits_view = View(
        Item('young_modulus', label='E [Pa]'),
        Item('poisson', label = 'Poisson [-]'),
        Item('density', label = 'Density [kg/m3]'),
        Item('fc', label = 'fc'),
        Item('ft', label = 'ft'),
        Item('KinN', label = 'KinN'),
        Item('gammaN', label = 'gammaN'),
        Item('m', label = 'm'),
        Item('Ad', label = 'Ad'),
        Item('tau_bar', label = 'tau_bar'),
        Item('Kin', label = 'Kin'),
        Item('gamma', label = 'gamma'),
        Item('S', label = 'S'),
        Item('a', label = 'a')
        )

    def get_material_object (self):
        return utilitiesMech.FatigueMaterial( self.young_modulus, self.poisson, self.density, self.fc, self.ft, self.KinN, self.gammaN, self.m, self.Ad, self.tau_bar, self.Kin, self.gamma, self.S, self.a)




class LinearTrsprtMaterial(TrsprtMaterial):
    name = 'Linear Transport Material'
    capacity = 1
    conductivity = 2

    traits_view = View(
        Item('capacity'),
        Item('conductivity')
        )

    def get_material_object (self):
        return utilitiesMech.TransportMaterial( self.capacity, self.conductivity)







mechMaterials = [
    LinearElasticMaterial(),
    MarsMaterial(),
    FatigueMaterial()
]
mechMaterials_dict = { m.name: m for m in mechMaterials }
mechMaterials_keys = list(mechMaterials_dict.keys())

trsprtMaterials = [
    LinearTrsprtMaterial()
]
trsprtMaterials_dict = { m.name: m for m in trsprtMaterials }
trsprtMaterials_keys = list(trsprtMaterials_dict.keys())
