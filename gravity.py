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


def get_orbital_elements_of_triple(particles):
    """
    Returns the eccentricity and semimajor axis of the inner and outer binaries
    :return:
    """
    inner_binary = particles[0] + particles[1]
    outer_binary = Particles(1)
    outer_binary[0].mass = inner_binary.mass.sum()
    outer_binary[0].position = inner_binary.center_of_mass()
    outer_binary[0].velocity = inner_binary.center_of_mass_velocity()
    outer_binary.add_particle(particles[2])
    _, _, semimajor_axis_in, eccentricity_in, _, _, _, _ \
        = orbital_elements_from_binary(inner_binary, G=constants.G)
    _, _, semimajor_axis_out, eccentricity_out, _, _, _, _ \
        = orbital_elements_from_binary(outer_binary, G=constants.G)
    return semimajor_axis_in, eccentricity_in, semimajor_axis_out, eccentricity_out


def get_orbital_period(orbital_separation, total_mass):
    return 2 * np.pi * (orbital_separation ** 3 / (constants.G * total_mass)).sqrt()


def get_semi_major_axis(orbital_period, total_mass):
    return (constants.G * total_mass * orbital_period ** 2 / (4 * np.pi ** 2)) ** (1. / 3)

