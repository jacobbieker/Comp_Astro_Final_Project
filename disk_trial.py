from __future__ import print_function
from amuse.lab import *
from Gadget2_Gravity import Gadget2_Gravity
from amuse.ext.sph_to_grid import convert_SPH_to_grid
import numpy as np
from SuperMassiveBlackHole import SuperMassiveBlackHole
from amuse.couple.bridge import Bridge
from amuse.community.ph4.interface import ph4
from BinaryBlackHole import BinaryBlackHole
from amuse.units.generic_unit_converter import ConvertBetweenGenericAndSiUnits

def main(N, Mtot, Rvir, t_end, dt):
    from amuse.ext.protodisk import ProtoPlanetaryDisk
    smbh = SuperMassiveBlackHole(mass=1e6 | units.MSun)
    inner_boundary = smbh.radius * 100
    outer_boundary = smbh.radius * 1000000
    disk_mass_fraction = 0.1
    disk_convert = nbody_system.nbody_to_si(smbh.super_massive_black_hole.mass, inner_boundary)
    gadget_convert = nbody_system.nbody_to_si(disk_mass_fraction*smbh.super_massive_black_hole.mass, outer_boundary)
    gas_particles = ProtoPlanetaryDisk(100000,
                                       convert_nbody=disk_convert,
                                       densitypower=1.0,
                                       Rmin=1.,
                                       Rmax=inner_boundary/outer_boundary,
                                       q_out=1.0,
                                       discfraction=disk_mass_fraction).result

    gas_particles.move_to_center()

    blackhole_masses = [30,30]
    all_grav = Particles()
    all_grav.add_particle(smbh.super_massive_black_hole)
    binary = BinaryBlackHole(blackhole_masses[0], blackhole_masses[1], smbh.super_massive_black_hole.mass,
                             initial_outer_semi_major_axis=inner_boundary * 5,
                             inner_eccentricity=0.6,
                             inclination=45.,
                             )
    all_grav.add_particles(binary.blackholes)

    gen_convert = ConvertBetweenGenericAndSiUnits(constants.c, units.s)
    hydro = Gadget2_Gravity(convert_nbody=gadget_convert,
                            number_of_workers=10)
    hydro.parameters.time_max = 2*gen_convert.to_generic(5 | units.Myr)
    hydro.gas_particles.add_particles(gas_particles)

    grav_converter = nbody_system.nbody_to_si(smbh.super_massive_black_hole.mass, smbh.super_massive_black_hole.radius)

    grav = ph4(convert_nbody=grav_converter, number_of_workers=2)
    grav.particles.add_particles(all_grav)

    channl_from_grav_to_particles = grav.particles.new_channel_to(all_grav)
    channel_from_hydro_to_disk = hydro.gas_particles.new_channel_to(gas_particles)
    write_set_to_file(all_grav, "Test_Grav_Particles_Initial.hdf5", "amuse")
    write_set_to_file(gas_particles, "Test_Gas_Particles_Initial.hdf5", "amuse")
    bridge = Bridge(verbose=True, use_threading=True)
    bridge.timestep = 1000 | units.yr
    start_time = 0. | units.Myr
    timestep = 10000 | units.yr
    end_time = 10 | units.Myr
    bridge.add_system(hydro, (grav, ))
    bridge.add_system(grav, (hydro,))
    while start_time < end_time:
        bridge.evolve_model(5 | units.Myr)
        channel_from_hydro_to_disk.copy()
        channl_from_grav_to_particles.copy()
        write_set_to_file(all_grav, "Test_Grav_Particles_Final.hdf5", "amuse")
        write_set_to_file(gas_particles, "Test_Gas_Particles_Final.hdf5", "amuse")
    grav.stop()
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
    main(**o.__dict__)
