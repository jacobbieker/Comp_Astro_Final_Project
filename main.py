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


smbh = SuperMassiveBlackHole(mass=1e7 | units.MSun)
smbh_mass = smbh.mass

blackhole_masses = np.random.uniform(low=10, high=15, size=2)

inner_boundary = (smbh.radius.value_in(units.parsec)) * 100
outer_boundary = (smbh.radius.value_in(units.parsec)) * 100000

binaries = BinaryBlackHole(blackhole_masses[0], blackhole_masses[1], smbh_mass,
                         initial_outer_semi_major_axis=np.random.uniform(inner_boundary, outer_boundary, size=1)[0] | units.parsec,
                         initial_outer_eccentricity=np.random.uniform(0, 0.99, size=1)[0],
                         eccentricity=np.random.uniform(0.0, 0.99, size=1),
                         inclincation=np.random.uniform(0.0, 180.0, size=1),
                         )


all_particles = Particles()
all_particles.add_particle(smbh.super_massive_black_hole)
all_particles.add_particles(binaries.blackholes)

converter = nbody_system.nbody_to_si(smbh_mass, 100000 * smbh.radius)
print ('smbh radius: ', smbh.radius.in_(units.parsec))
gravity = Huayno(converter)


gravity.particles.add_particle(all_particles)
print(gravity.particles)
print('semi major axis: ', binaries.initial_outer_semi_major_axis.in_(units.AU), \
      '\nbinaries hill radius: ', binaries.hill_radius.in_(units.AU), \
      '\neccentricity: ', binaries.initial_outer_eccentricity, \
      '\ntotal mass: ', binaries.total_mass.in_(units.MSun), \
      '\nsmbh mass: ', binaries.central_blackhole.mass.in_(units.MSun), \
      '\nbinaries max orbital period: ', binaries.binary_max_orbital_period.in_(units.yr), \
      '\nbinaries min orbital period: ', binaries.binary_min_orbital_period.in_(units.yr), \
      '\nblackhole1 position: ', binaries.blackholes[0].position.in_(units.AU), \
      '\nblackhole2 position: ', binaries.blackholes[1].position.in_(units.AU), \
      '\nbinaries mass: ', binaries.blackholes.mass.in_(units.MSun), \
      '\nbinaries distance: ', binaries.blackholes_distance.in_(units.AU))

# channel_from_grav_to_binaries = gravity.particles.new_channel_to(all_particles)


# gravity.stop()




















#
#
# from AccretionDisk import AccretionDisk
# from SuperMassiveBlackHole import SuperMassiveBlackHole
# from BinaryBlackHole import BinaryBlackHole
# from amuse.ic.plummer import new_plummer_model
# from amuse.datamodel import Particle, Particles
# from amuse.couple.bridge import Bridge
# from amuse.community.huayno.interface import Huayno
# from amuse.units import units
# import numpy as np
# from BinaryBlackHolesWithAGN import BinaryBlackHolesWithAGN
#
# simulation = BinaryBlackHolesWithAGN(1e6, 2, 0.5, converter=converter)
