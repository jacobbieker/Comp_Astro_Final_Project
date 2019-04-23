from amuse.units import units
from BinaryBlackHolesWithAGN import BinaryBlackHolesWithAGN



simulation = BinaryBlackHolesWithAGN(mass_of_central_black_hole=1e6 | units.MSun,
                                     number_of_binaries=50,
                                     number_of_gas_particles=1000000,
                                     disk_mass_fraction=0.1,
                                     number_of_hydro_workers=4,
                                     number_of_grav_workers=5,
                                     end_time=100 | units.Myr)
