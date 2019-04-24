from amuse.units import units
from BinaryBlackHolesWithAGN import BinaryBlackHolesWithAGN
import numpy as np

file_num = np.random.randint(0,100000, size=1)
simulation = BinaryBlackHolesWithAGN(mass_of_central_black_hole=1e6 | units.MSun,
                                     number_of_binaries=50,
                                     number_of_gas_particles=0,
                                     disk_mass_fraction=0.1,
                                     number_of_hydro_workers=6,
                                     number_of_grav_workers=12,
                                     end_time=10 | units.Myr,
                                     filename=str(file_num))
