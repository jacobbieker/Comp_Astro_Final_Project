from __future__ import print_function
from amuse.lab import *
from Gadget2_Extended import Gadget2_Extended
from amuse.ext.sph_to_grid import convert_SPH_to_grid
from matplotlib import pyplot as plt
import numpy  as np
from SuperMassiveBlackHole import SuperMassiveBlackHole
from amuse.couple.bridge import Bridge
from amuse.community.ph4.interface import ph4

def gas_sphere(N, Mtot, Rvir):
    converter = nbody_system.nbody_to_si(Mtot, Rvir)
    gas = new_plummer_gas_model(N, convert_nbody=converter)
    return gas, converter


#
#
def make_map(hydro, grid_points=100, L=1):
    x, y = np.indices((grid_points + 1, grid_points + 1))

    x = L * (x.flatten() - grid_points / 2.) / grid_points
    y = L * (y.flatten() - grid_points / 2.) / grid_points
    z = x * 0.
    vx = 0. * x
    vy = 0. * x
    vz = 0. * x

    x = units.parsec(x)
    y = units.parsec(y)
    z = units.parsec(z)
    vx = units.kms(vx)
    vy = units.kms(vy)
    vz = units.kms(vz)

    rho, rhovx, rhovy, rhovz, rhoe = hydro.get_hydro_state_at_point(x, y, z, vx, vy, vz)
    rho = rho.reshape((grid_points + 1, grid_points + 1))

    return rho

def main(N, Mtot, Rvir, t_end, dt):
    from amuse.ext.protodisk import ProtoPlanetaryDisk
    smbh = SuperMassiveBlackHole(mass=1e6 | units.MSun)
    inner_boundary = smbh.radius * 100
    outer_boundary = smbh.radius * 100000
    disk_mass_fraction = 0.1
    disk_convert = nbody_system.nbody_to_si(smbh.super_massive_black_hole.mass, inner_boundary)
    gadget_convert = nbody_system.nbody_to_si(disk_mass_fraction*smbh.super_massive_black_hole.mass, outer_boundary)
    gas_particles = ProtoPlanetaryDisk(100000,
                                       convert_nbody=disk_convert,
                                       densitypower=1.0,
                                       Rmin=1.,
                                       Rmax=1e4,
                                       q_out=1.0,
                                       discfraction=0.1).result

    gas_particles.move_to_center()

    hydro = Gadget2_Extended(radius=40 | units.AU, convert_nbody=gadget_convert,
                    number_of_workers=6)
    hydro.gas_particles.add_particles(gas_particles)

    grav_converter = nbody_system.nbody_to_si(smbh.super_massive_black_hole.mass, smbh.super_massive_black_hole.radius)

    grav = ph4(convert_nbody=grav_converter)
    grav.particles.add_particle(smbh.super_massive_black_hole)
    bridge = Bridge(verbose=True)
    bridge.timestep = 100 | units.yr
    bridge.add_system(hydro, (grav, ))
    bridge.add_system(grav, (hydro,))

    bridge.evolve_model(1000 | units.yr)
    bridge.evolve_model(500 | units.yr)
    grav.stop()
    hydro.stop()

    rho = make_map(hydro, 100)

    plt.imshow(rho.value_in(units.amu / units.cm ** 3))
    plt.savefig('density_map_disk_converter_nolog.pdf')
    plt.show()

    grid = convert_SPH_to_grid(hydro, (100,100,250), do_scale=True)


    #rho = make_map(hydro, 100)
    '''
    plt.imshow(np.log10(1.e-5 + rho.value_in(units.amu / units.cm ** 3)))
    plt.savefig('density_map_gridded_converter.pdf')
    plt.show()

    summed_density = np.sum(grid.rho.value_in(grid.rho.unit), axis=0)

    plt.imshow(summed_density)
    plt.title("X")
    plt.show()

    summed_density = np.sum(grid.rho.value_in(grid.rho.unit), axis=1)

    plt.imshow(summed_density)
    plt.title("Y")
    plt.show()

    summed_density = np.sum(grid.rho.value_in(grid.rho.unit), axis=2)

    plt.imshow(summed_density)
    plt.title("Z")
    plt.show()
    '''
    #write_set_to_file(grid, "Unitless_Disk_test.h5", "hdf5")

    print("Did Grid")

    #rho = make_map(hydro, 100)

    plt.imshow(np.log10(1.e-5 + rho.value_in(units.amu / units.cm ** 3)))
    plt.savefig('density_map.pdf')
    plt.show()
    Etot_init = hydro.kinetic_energy \
                + hydro.potential_energy + hydro.thermal_energy

    time = 0 | units.yr
    while time < t_end:
        time += dt
        hydro.evolve_model(time)
        print(time, flush=True)
    write_set_to_file(hydro.particles, "hydro.h5", "hdf5")

    Ekin = hydro.kinetic_energy
    Epot = hydro.potential_energy
    Eth = hydro.thermal_energy
    Etot = Ekin + Epot + Eth
    Q = (Ekin + Eth) / Epot
    dE = (Etot_init - Etot) / Etot
    com = hydro.gas_particles.center_of_mass()

    print("T=", hydro.get_time(), "M=", hydro.gas_particles.mass.sum(), end=' ')
    print("E= ", Etot, "Q= ", Q, "dE=", dE, "CoM=", com.in_(units.RSun))

    x_gas = gas.x.value_in(units.AU)
    y_gas = gas.y.value_in(units.AU)
    z_gas = gas.z.value_in(units.AU)

    z_0 = 0. * x_gas
    vx_gas = (0. | (units.Hz)) * x_gas
    vy_gas = (0. | (units.Hz)) * x_gas
    vz_gas = (0. | (units.Hz)) * x_gas

    rho = make_map(hydro)

    plt.imshow(np.log10(1.e-5 + rho.value_in(units.amu / units.cm ** 3)))
    plt.savefig('density_map.pdf')
    plt.show()

    hydro.stop()


# plt.scatter(x_gas, y_gas, z_gas)
# plt.savefig('disk_scatter.png')
# plt.show()

def new_option_parser():
    from amuse.units.optparse import OptionParser
    result = OptionParser()

    result.add_option("-N", dest="N", type="int",
                      default=10000,
                      help="no. of particles [%default]")
    result.add_option("--M_tot", unit=units.MSun, dest="Mtot", type="float",
                      default=5e6 | units.MSun,
                      help="Total disk mass [%default]")
    result.add_option("--Rvir", unit=units.AU, dest="Rvir", type="float",
                      default=1974626.290440254 | units.AU,
                      help="Radius of initial sphere [%default]")
    result.add_option("--t_end", unit=units.yr, dest="t_end", type="float",
                      default=5 | units.Myr,
                      help="End time of simulation [%default]")
    result.add_option("--dt", unit=units.yr, dest="dt", type="float",
                      default=10 | units.yr,
                      help="Hydro timestep [%default]")
    return result


if __name__ in ("__main__"):
    o, arguments = new_option_parser().parse_args()
    print(o.__dict__)
    main(**o.__dict__)
