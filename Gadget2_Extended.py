from amuse.community.gadget2.interface import Gadget2
import numpy as np
from amuse.units import constants, units
from amuse.couple.bridge import CalculateFieldForParticles

class Gadget2_Extended(Gadget2):
    def __init__(self, radius=None, unit_converter=None, mode='normal', **options):
        Gadget2.__init__(self, unit_converter=unit_converter, mode=mode, **options)
        self.radius_for_gravity_calc = radius

    def get_gravity_at_point(self, radius, x, y, z):
        field_code = CalculateFieldForParticles(particles=self.gas_particles, G=constants.G)
        return field_code.get_gravity_at_point(radius, x, y, z)

    def get_potential_at_point(self, radius, x, y, z):
        field_code = CalculateFieldForParticles(particles=self.gas_particles, G=constants.G)
        return field_code.get_potential_at_point(radius, x, y, z)

    def add_particles(self, particles):
        self.gas_particles.add_particles(particles)

