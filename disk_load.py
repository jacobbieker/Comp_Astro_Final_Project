from __future__ import print_function
from amuse.lab import *
from amuse.community.gadget2.interface import Gadget2
from amuse.ext.sph_to_grid import convert_SPH_to_grid
from matplotlib import pyplot as plt
import numpy  as np


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


gas, converter = gas_sphere(10000, 5e6 | units.MSun, 19.756 | units.AU)

grid = read_set_from_file("/home/jacob/Development/Comp_Astro_Final_Project/Unitless_Disk_medium.h5", "hdf5")

from mpl_toolkits.mplot3d import Axes3D
from matplotlib import pyplot as plt

fig = plt.figure(figsize=(10, 10))
ax = fig.add_subplot(111, projection='3d')
print(grid.x[:,0,0])
graph = ax.scatter(grid.x.value_in(grid.x.unit), grid.y.value_in(grid.y.unit), grid.rho.value_in(grid.rho.unit))
plt.show()


hydro = Gadget2(unit_converter=converter, mode='periodic')

hydro.gas_particles.add_particles(grid)