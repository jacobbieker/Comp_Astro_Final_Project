from amuse.units import units
from BinaryBlackHolesWithAGN import BinaryBlackHolesWithAGN
import numpy as np

def new_option_parser():
    from amuse.units.optparse import OptionParser

    file_num = np.random.randint(0, 100000, size=1)

    result = OptionParser()
    result.add_option("--mass_of_central_black_hole", unit=units.MSun, dest="mass_of_central_black_hole", type="float",
                      default=1e6 | units.MSun, help="Mass of central SMBH [%default]")
    result.add_option("--number_of_binaries", dest="number_of_binaries", type="int", default=50,
                      help="No. of binaries [%default]")
    result.add_option("--number_of_gas_particles", dest="number_of_gas_particles", type="int", default=100000,
                      help="No. of gas particles [%default]")
    result.add_option("--end_time", unit=units.Myr, dest="end_time", type="float", default=10 | units.Myr,
                      help="End time of simulation [%default]")
    result.add_option("--disk_mass_fraction", dest="disk_mass_fraction", type="float", default=0.1,
                      help="Disk mass fraction [%default]")
    result.add_option("--number_of_hydro_workers", dest="number_of_hydro_workers", type="int", default=6,
                      help="Number of workers for hydro code [%default]")
    result.add_option("--number_of_grav_workers", dest="number_of_grav_workers", type="int", default=12,
                      help="Number of workers for gravity code [%default]")
    result.add_option("--filename", dest="filename", type="string", default=str(file_num),
                      help="Filename [%default]")
    return result


def main(mass_of_central_black_hole,
         number_of_binaries,
         number_of_gas_particles,
         end_time,
         disk_mass_fraction,
         number_of_hydro_workers,
         number_of_grav_workers,
         filename):

    simulation = BinaryBlackHolesWithAGN(mass_of_central_black_hole,
                                         number_of_binaries,
                                         number_of_gas_particles,
                                         disk_mass_fraction,
                                         number_of_hydro_workers,
                                         number_of_grav_workers,
                                         end_time,
                                         filename)


if __name__ == "__main__":
    o, arguments = new_option_parser().parse_args()
    main(**o.__dict__)
