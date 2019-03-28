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

    def __init__(self, mass_one, mass_two, orbital_period, eccentricity=0.0, inclincation=0.0,
                 orbital_fraction_timestep=0.5):
        """
        This model is to generate a binary black hole system when given an initial Particle to split

        Create binaries, and then distribute them according to Plummer sphere

        Have it set center of mass velocity and venter of mass location

        Keep binaries by themselves


        """
        self.orbital_period = orbital_period
        self.eccentricity = eccentricity
        self.inclincation = inclincation
        self.blackholes = Particles(2)
        self.blackholes.mass = [mass_one, mass_two] | units.MSun
        self.total_mass = self.blackholes.mass.sum()
        self.semi_major_axis = self.get_semi_major_axis()
        self.timestep = self.orbital_period * orbital_fraction_timestep

        binary_position, binary_velocity = get_position(self.blackholes[0].mass, self.blackholes[1].mass,
                                                        self.eccentricity, self.semi_major_axis,
                                                        180, self.inclincation, 180, 0, self.timestep)
        self.blackholes[1].position = binary_position
        self.blackholes[1].velocity = binary_velocity
        self.blackholes.move_to_center()

    def particles(self):
        return self.blackholes

    def get_orbital_period(self, orbital_separation, total_mass):
        return 2 * np.pi * (orbital_separation ** 3 / (constants.G * total_mass)).sqrt()

    def get_semi_major_axis(self):
        return (constants.G * self.total_mass * self.orbital_period ** 2 / (4 * np.pi ** 2)) ** (1. / 3)

    def set_center_of_mass(self, new_center_of_mass):
        """
        Sets the center of mass of the binary system
        :param new_center_of_mass: center of mass needed
        :return:
        """

        try:
            self.blackholes[0].position -= new_center_of_mass
            self.blackholes[1].position -= new_center_of_mass
        except:
            # Need to convert to a list with units, instead of a list of elements with units
            for index, element in enumerate(new_center_of_mass):
                new_center_of_mass[index] = element.value_in(units.m)
            new_center_of_mass = new_center_of_mass | units.m
            self.blackholes[0].position -= new_center_of_mass
            self.blackholes[1].position -= new_center_of_mass

    def set_center_of_mass_velocity(self, new_center_of_mass_velocity):
        """
        Sets the center of mass velocity of the binary system
        :param new_center_of_mass_velocity:
        :return:
        """
        try:
            self.blackholes[0].position -= new_center_of_mass_velocity
            self.blackholes[1].position -= new_center_of_mass_velocity
        except:
            # Need to convert to a list with units, instead of a list of elements with units
            for index, element in enumerate(new_center_of_mass_velocity):
                new_center_of_mass_velocity[index] = element.value_in(units.m / units.s)
            new_center_of_mass = new_center_of_mass_velocity | units.m / units.s
            self.blackholes[0].velocity -= new_center_of_mass
            self.blackholes[1].velocity -= new_center_of_mass

    def set_binary_location_and_velocity(self, center_of_mass, center_of_mass_velocity):
        self.set_center_of_mass(center_of_mass)
        self.set_center_of_mass_velocity(center_of_mass_velocity)

    def set_merge_conditions(self):
        raise NotImplementedError

    def merge_particles(self):
        """
        Merges the binary particles into a single particle, with the
        :return:
        """
        raise NotImplementedError
