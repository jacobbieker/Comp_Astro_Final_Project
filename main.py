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


smbh = SuperMassiveBlackHole(mass=1e6 | units.MSun)
smbh_mass = smbh.mass

binaries = BinaryBlackHole(10, 10, smbh_mass, orbital_period=10 | units.yr)

all_particles = Particles()
all_particles.add_particle(smbh.super_massive_black_hole)
all_particles.add_particles(binaries.blackholes)

converter = nbody_system.nbody_to_si(smbh_mass, 100000 * smbh.radius)
print (smbh.radius.in_(units.parsec))
gravity = Huayno(converter)


gravity.particles.add_particle(all_particles)
print (gravity.particles)
# gravity.particles.add_particle(smbh.super_massive_black_hole)
# gravity.particles.add_particle(binaries.blackholes)
#
# channel_from_grav_to_binaries = gravity.particles.new_channel_to(all_particles)


# gravity.stop()
