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
import matplotlib.animation
from amuse.lab import *
from amuse.units import units
from amuse.io import read_set_from_file


def get_orbital_elements_of_triple(particles):
    """
    Returns the eccentricity and semimajor axis of the inner and outer binaries
    :return:
    """
    inner_binary = particles[2] + particles[1]
    outer_binary = Particles(1)
    outer_binary[0].mass = inner_binary.mass.sum()
    outer_binary[0].position = inner_binary.center_of_mass()
    outer_binary[0].velocity = inner_binary.center_of_mass_velocity()
    outer_binary.add_particle(particles[0])
    _, _, semimajor_axis_in, eccentricity_in, _, _, _, _ \
        = orbital_elements_from_binary(inner_binary, G=constants.G)
    _, _, semimajor_axis_out, eccentricity_out, _, _, _, _ \
        = orbital_elements_from_binary(outer_binary, G=constants.G)
    return semimajor_axis_in, eccentricity_in, semimajor_axis_out, eccentricity_out


def get_orbital_period(orbital_separation, total_mass):
    return 2 * np.pi * (orbital_separation ** 3 / (constants.G * total_mass)).sqrt()


def get_semi_major_axis(orbital_period, total_mass):
    return (constants.G * total_mass * orbital_period ** 2 / (4 * np.pi ** 2)) ** (1. / 3)


def main():
    bodies = Particles(3)

    bodies[0].mass = 1e9 | units.MSun
    bodies[0].radius = 0.5 | units.RSun

    bodies[0].position = (0,0,0) | units.km
    bodies[0].velocity = (0,0,0) | (units.m/units.s)

    bodies.mass = [1e9, 1, 2] | units.MSun

    binary = Particles(2)
    binary[0].mass = bodies.mass[1]
    binary[1].mass = bodies.mass[2]

    ain_0 = get_semi_major_axis(500 | units.yr, binary.mass.sum())
    print(ain_0.value_in(units.AU))
    timestep = 250. | units.yr

    r,v = get_position(bodies[2].mass, bodies[1].mass, 0.2, ain_0, 180, 0, 180, 0, timestep)
    binary[1].position = r
    binary[1].velocity = v
    binary.move_to_center()

    ###### Now Need to make sure binary is on its own away from the SMBH ############

    r,v = get_position(bodies[0].mass, bodies[1].mass+bodies[2].mass, 0.0, ain_0*10, 0, 45, 0, 0, timestep)

    # Convert to center of mass of the binary system, but in nbody units as of right now
    print(binary[0].position)
    print(binary.center_of_mass())
    print(v)
    for i, e in enumerate(r):
        r[i] = e.value_in(units.m)
        v[i] = v[i].value_in(units.m/units.s)
    r = r | units.m
    v = v | units.m / units.s
    binary[1].position -= r
    binary[1].velocity -= v
    binary[0].position -= r
    binary[0].velocity -= v

    bodies[1].position = binary[0].position
    bodies[1].velocity = binary[0].velocity
    bodies[2].position = binary[1].position
    bodies[2].velocity = binary[1].velocity

    bodies.move_to_center()

    aout = 0.0002 | units.parsec

    converter = nbody_system.nbody_to_si(bodies.mass.sum(), aout)

    gravity = SmallN(converter)
    gravity.particles.add_particles(bodies)

    channel_from_grav = gravity.particles.new_channel_to(bodies)

    gravity.particles.move_to_center()
    channel_from_grav.copy()

    Pout = get_orbital_period(aout, bodies.mass.sum())

    time_start = 0. | units.Myr
    time = 1.0 | units.Myr

    ains = []
    eins = []
    from mpl_toolkits.mplot3d import Axes3D
    import matplotlib.pyplot as plt
    x = []
    y = []
    z = []
    cm_x = []
    bh_x = []
    cm_y = []
    bh_y = []

    outer_timing = 0 | time_start.unit
    ain_0 = 500 | units.yr
    while time_start < time:
        time_start += timestep
        print(time_start.value_in(units.Myr), flush=True)

        gravity.evolve_model(time_start)
        channel_from_grav.copy()

        ain, ein, aout, eout = get_orbital_elements_of_triple(bodies)

        ains.append(ain.value_in(units.AU))
        eins.append(ein)
        outer_timing += timestep
        if outer_timing >= ain_0*5:
            bodies[1].vx = 0.01*bodies[1].vx
            bodies[2].vx = 0.01*bodies[2].vx
            bodies[1].vy = 0.01*bodies[1].vy
            bodies[2].vy = 0.01*bodies[2].vy
            bodies[1].vz = 0.01*bodies[1].vz
            bodies[2].vz = 0.01*bodies[2].vz
            outer_timing = outer_timing % (ain_0*5)
        x.append((bodies[1].x.value_in(units.AU), bodies[2].x.value_in(units.AU)))
        y.append((bodies[1].y.value_in(units.AU), bodies[2].y.value_in(units.AU)))
        z.append((bodies[1].z.value_in(units.AU), bodies[2].z.value_in(units.AU)))

        cm_x.append(bodies[1:].center_of_mass()[0].value_in(units.AU))
        cm_y.append(bodies[1:].center_of_mass()[1].value_in(units.AU))

        bh_x.append(bodies[0].x.value_in(units.AU))
        bh_y.append(bodies[0].y.value_in(units.AU))
        #print(bh_x)


    gravity.stop()

    plt.plot([i[0] for i in x], [i[0] for i in y])
    plt.plot([i[1] for i in x], [i[1] for i in y])
    plt.show()
    cm_y = np.asarray(cm_y)
    cm_x = np.asarray(cm_x)
    distances = (cm_x - bodies[0].x.value_in(units.AU))**2 + (cm_y - bodies[0].y.value_in(units.AU))**2
    print(min(distances))
    plt.plot([i for i in range(len(distances))], distances)
    plt.show()
    print(distances[distances <= ((2*constants.G*bodies[0].mass)/(constants.c**2))])
    plt.plot(cm_x, cm_y)
    plt.scatter(bh_x, bh_y, c='g', s=500)
    plt.show()
    t = np.asarray([i for i in range(len(ains))])
    print(max(ains))
    print(min(ains))
    print(max(eins))
    print(min(eins))
    print(aout.value_in(units.AU))

main()
