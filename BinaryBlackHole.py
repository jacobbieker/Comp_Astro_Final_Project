from amuse.lab import *
import numpy as np
from amuse.ext.solarsystem import get_position
from amuse.units import units, constants, nbody_system
from amuse.units.quantities import zero
from amuse.datamodel import Particle, Particles
from amuse.ext.orbital_elements import new_binary_from_orbital_elements
from amuse.ext.orbital_elements import orbital_elements_from_binary
from amuse.community.huayno.interface import Huayno
from amuse.community.smalln.interface import SmallN
from amuse.community.hermite0.interface import Hermite


class BinaryBlackHole(object):

    def __init__(self):
        """
        This model is to generate a binary black hole system when given an initial Particle to split

        Create binaries, and then distribute them according to Plummer sphere

        Have it set center of mass velocity and venter of mass location

        Keep binaries by themselves


        """
        raise NotImplementedError

    def get_orbital_period(self, orbital_separation, total_mass):
        return 2 * np.pi * (orbital_separation ** 3 / (constants.G * total_mass)).sqrt()

    def get_semi_major_axis(self, orbital_period, total_mass):
        return (constants.G * total_mass * orbital_period ** 2 / (4 * np.pi ** 2)) ** (1. / 3)

    def center_of_mass(self):
        raise NotImplementedError

    def center_of_mass_velocity(self):
        raise NotImplementedError

    def set_binary_location_and_velocity(self, center_of_mass, center_of_mass_velocity):
        raise NotImplementedError

    def set_merge_conditions(self):
        raise NotImplementedError

    def merge_particles(self):
        raise NotImplementedError