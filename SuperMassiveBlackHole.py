from __future__ import division, print_function
from AccretionDisk import AccretionDisk
import numpy as np
from amuse.units import units, constants
from amuse.datamodel import Particle


class SuperMassiveBlackHole(object):

    def __init__(self, mass=1E8 | units.MSun):
        self.mass = mass
        self.super_massive_black_hole = Particle()
        self.super_massive_black_hole.mass = self.mass
        self.radius = (2*constants.G*self.super_massive_black_hole.mass)/(constants.c**2)
        self.super_massive_black_hole.position = (0, 0, 0) | units.AU
        self.super_massive_black_hole.velocity = (0, 0, 0) | units.kms
        self.super_massive_black_hole.radius = self.radius
