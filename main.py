from __future__ import division, print_function
import numpy as np
import matplotlib.pyplot as plt
from SuperMassiveBlackHole import SuperMassiveBlackHole
from BinaryBlackHole import BinaryBlackHole
from AccretionDisk import AccretionDisk
from amuse.lab import Particle, units, nbody_system, constants, Particles
from amuse.community.huayno.interface import Huayno
from amuse.community.ph4.interface import ph4
from mpl_toolkits.mplot3d import Axes3D


# from amuse.units import units, constants
# from amuse.datamodel import Particles, Particle


number_of_binaries = 1000
steps_of_inclination = 19

smbh = SuperMassiveBlackHole(mass=1e6 | units.MSun)
smbh_mass = smbh.mass

inner_boundary = smbh.radius*1e2
outer_boundary = smbh.radius*1e6

print(inner_boundary.in_(units.parsec), outer_boundary.in_(units.parsec))


all_gravity_particles = Particles()
initial_outer_semi_major_axis = np.linspace(inner_boundary.value_in(outer_boundary.unit), outer_boundary.value_in(outer_boundary.unit), number_of_binaries//steps_of_inclination) | (outer_boundary.unit)
initial_inclination = np.linspace(0,180, steps_of_inclination)

sma_incl_list = []

for i in range(len(initial_inclination)):
    for j in range(len(initial_outer_semi_major_axis)):
        sma_incl_list.append([initial_inclination[i], initial_outer_semi_major_axis[j]])
for i in range(len(sma_incl_list)):
    blackhole_masses = [30,30]

    binaries = BinaryBlackHole(blackhole_masses[0], blackhole_masses[1], smbh_mass,
                               initial_outer_semi_major_axis=sma_incl_list[i][1],
                               initial_outer_eccentricity=0.6,
                               inner_eccentricity=0.6,
                               inclination=sma_incl_list[i][0],
                               )


    # print (sma_incl_list[i][1].in_(units.parsec))
    all_gravity_particles.add_particles(binaries.blackholes)

all_gravity_particles.add_particle(smbh.super_massive_black_hole)

converter = nbody_system.nbody_to_si(all_gravity_particles.mass.sum(), all_gravity_particles.virial_radius())
'''
gravity = ph4(converter)
gravity.particles.add_particles(all_gravity_particles)
channel_from_grav_to_binaries = gravity.particles.new_channel_to(all_gravity_particles)
channel_from_binaries_to_grav = all_gravity_particles.new_channel_to(gravity.particles)

# ----------------------- must become parameter -----------------------#
end_time = 5. | units.Myr
timestep = 0.1 | end_time.unit
# ---------------------------------------------------------------------#
sim_time = 0. | end_time.unit

while sim_time < end_time:
    sim_time += timestep
    print('lego')
    gravity.evolve_model(sim_time)
    print('letsgo')

    channel_from_grav_to_binaries.copy()


# print(gravity.particles)
'''
print('initial outer semi major axis: ', binaries.initial_outer_semi_major_axis.in_(units.AU), \
      '\nbinaries hill radius: ', binaries.hill_radius.in_(units.AU), \
      '\ninitial outereccentricity: ', binaries.initial_outer_eccentricity, \
      '\nbinaries max orbital period: ', binaries.binary_max_orbital_period.in_(units.yr), \
      '\nbinaries min orbital period: ', binaries.binary_min_orbital_period.in_(units.yr), \
      '\ntotal binary mass: ', binaries.total_mass.in_(units.MSun), \
      '\nsmbh mass: ', binaries.central_blackhole.mass.in_(units.MSun), \
      '\nbinary blackhole mass: ', binaries.blackholes.mass.in_(units.MSun), \
      '\nsmbh radius: ', smbh.radius.in_(units.AU), \
      '\nbinary blackhole distance: ', binaries.blackholes_distance.in_(units.AU))


fig = plt.figure(figsize=(10, 10))
ax = fig.add_subplot(111, projection='3d')
graph = ax.scatter(all_gravity_particles.x.value_in(units.parsec), all_gravity_particles.y.value_in(units.parsec), all_gravity_particles.z.value_in(units.parsec))
ax.scatter(0,0,0, color='red')
plt.savefig('binaries_positions.pdf')
plt.show()
# gravity.stop()
