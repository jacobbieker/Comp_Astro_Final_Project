from AccretionDisk import AccretionDisk
import numpy as np
from amuse.units import units, constants
from amuse.datamodel import Particle


class SuperMassiveBlackHole(object):

    def __init__(self, mass=1E8 | units.MSun, has_accretion_disk=False, fraction_in_disk=0.1):
        self.mass = mass
        self.has_accretion_disk = has_accretion_disk
        self.fraction_in_disk = fraction_in_disk

        self.super_massive_black_hole = Particle()
        self.super_massive_black_hole.mass = self.mass
        self.radius = (2*constants.G*self.super_massive_black_hole.mass)/(constants.c**2)
        self.super_massive_black_hole.position = (0,0,0) | units.AU
        self.super_massive_black_hole.velocity = (0,0,0) | units.kms
        self.super_massive_black_hole.radius = self.radius
        self.disk = None

    def create_accretion_disk(self):
        """
        Needs the accretion disk to start at 2 times the scwarzchild radius, and goes out to 10 radii

        Binaries start outside of 2 radii, just select the ones outside this to add them

        100 times schwarzhild radius


        :return:
        """

        self.disk = AccretionDisk(disk_min=2*(self.radius.value_in(units.parsec)),
                                  disk_max=100*(self.radius.value_in(units.parsec)),
                                  fraction_of_central_blackhole_mass=self.fraction_in_disk,
                                  )

        self.has_accretion_disk = True

        raise NotImplementedError