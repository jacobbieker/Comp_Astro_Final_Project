from amuse.units import units
from BinaryBlackHolesWithAGN import BinaryBlackHolesWithAGN
from amuse.units.optparse import OptionParser


def new_option_parser():
    result = OptionParser()
    result.add_option("--mass_of_central_black_hole", unit=units.MSun, dest="mass_of_central_black_hole", type="float",
                      default=1e6 | units.MSun, help="Mass of central SMBH [%default]")
    result.add_option("--number_of_binaries", dest="number_of_binaries", type="int", default=50,
                      help="No. of binaries [%default]")
    result.add_option("--number_of_gas_particles", dest="number_of_gas_particles", type="int", default=100000,
                      help="No. of gas particles [%default]")
    result.add_option("--end_time", unit=units.Myr, dest="end_time", type="float", default=10 | units.Myr,
                      help="End time of simulation [%default]")
    result.add_option("--blackhole_mass", unit=units.MSun, dest="blackhole_mass", type="float", default=30 | units.MSun,
                      help="Mass of the Blackholes in the binaries [%default]")
    result.add_option("--gravity_timestep", unit=units.yr, dest="gravity_timestep", type="float",
                      default=100 | units.yr,
                      help="Timestep for the gravity code [%default]")
    result.add_option("--bridge_timestep", unit=units.Myr, dest="bridge_timestep", type="float", default=0.1 | units.Myr,
                      help="Timestep for the Bridge [%default]")
    result.add_option("--smbh_as_potential", dest="smbh_as_potential",action="store_true", default=False,
                      help="Whether to simulate the SMBH as a potential vs a particle, 0 is false, 1 is true [%default]")
    result.add_option("--binaries_affect_disk", dest="binaries_affect_disk", action="store_true", default=False,
                      help="Whether the binaries affect the accretion disk [%default]")
    result.add_option("--disk_mass_fraction", dest="disk_mass_fraction", type="float", default=0.1,
                      help="Disk mass fraction [%default]")
    result.add_option("--number_of_hydro_workers", dest="number_of_hydro_workers", type="int", default=6,
                      help="Number of workers for hydro code [%default]")
    result.add_option("--number_of_grav_workers", dest="number_of_grav_workers", type="int", default=12,
                      help="Number of workers for gravity code [%default]")
    result.add_option("--filename", dest="filename", type="string", default="BinaryBlackHoles",
                      help="Filename [%default]")

    return result


def main(mass_of_central_black_hole,
         number_of_binaries,
         number_of_gas_particles,
         end_time,
         gravity_timestep,
         bridge_timestep,
         blackhole_mass,
         smbh_as_potential,
         binaries_affect_disk,
         disk_mass_fraction,
         number_of_hydro_workers,
         number_of_grav_workers,
         filename):
    BinaryBlackHolesWithAGN(mass_of_central_black_hole=mass_of_central_black_hole,
                            number_of_binaries=number_of_binaries,
                            number_of_gas_particles=number_of_gas_particles,
                            disk_mass_fraction=disk_mass_fraction,
                            binaries_affect_disk=binaries_affect_disk,
                            smbh_as_potential=smbh_as_potential,
                            blackhole_masses=blackhole_mass,
                            timestep=bridge_timestep,
                            gravity_timestep=gravity_timestep,
                            end_time=end_time,
                            number_of_hydro_workers=number_of_hydro_workers,
                            number_of_grav_workers=number_of_grav_workers,
                            filename=filename)


if __name__ == "__main__":
    o, arguments = new_option_parser().parse_args()
    print(o)
    main(**o.__dict__)
