from AccretionDisk import AccretionDisk
from SuperMassiveBlackHole import SuperMassiveBlackHole
from BinaryBlackHole import BinaryBlackHole
from amuse.ic.plummer import new_plummer_model
from amuse.datamodel import Particle, Particles
from amuse.couple.bridge import Bridge
from amuse.community.huayno.interface import Huayno
from amuse.units import units
import numpy as np
from BinaryBlackHolesWithAGN import BinaryBlackHolesWithAGN

simulation = BinaryBlackHolesWithAGN(mass_of_central_black_hole=1e7 | units.MSun,
                                     number_of_binaries=5,
                                     number_of_gas_particles=100,
                                     disk_mass_fraction=0.5,
                                     number_of_workers=4)