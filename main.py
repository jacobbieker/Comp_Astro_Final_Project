from __future__ import division, print_function
import numpy as np
import matplotlib.pyplot as plt
from SuperMassiveBlackHole import SuperMassiveBlackHole
from BinaryBlackHole import BinaryBlackHole
from AccretionDisk import AccretionDisk
from amuse.lab import Particle, units, nbody_system, constants, Particles
from amuse.community.huayno.interface import Huayno
# from amuse.units import units, constants
# from amuse.datamodel import Particles, Particle


SMBH = SuperMassiveBlackHole(mass=1e6 | units.MSun)
SMBH_mass = SMBH.mass

binaries = BinaryBlackHole(10, 10, SMBH_mass, orbital_period=10 | units.yr)

# print(SMBH.mass)
# print(binaries.total_mass)
# print(binaries.eccentricity)
# print(binaries.orbital_period)
print(SMBH.super_massive_black_hole)

converter = nbody_system.nbody_to_si(1 | units.MSun, 1 | units.AU)
gravity = Huayno(converter)

gravity.particles.add_particle(SMBH.super_massive_black_hole)
gravity.particles.add_particle(binaries.blackholes)
#
# channel_from_grav_to_binaries = gravity.particles.new_channel_to(all_particles)


# gravity.stop()