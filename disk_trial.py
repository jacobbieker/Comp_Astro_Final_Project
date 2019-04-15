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


#
#
# def plot_hydro(time, sph, i, L=10):
#     x_label = "x [pc]"
#     y_label = "y [pc]"
#     # fig = single_frame(x_label, y_label, logx=False, logy=False, xsize=12, ysize=12)
#
# 	fig=pyplot.figure(figsize=(12,12))
#
#     rho=make_map(sph,N=200,L=L)
#     pyplot.imshow(numpy.log10(1.e-5+rho.value_in(units.amu/units.cm**3)), extent=[-L/2,L/2,-L/2,L/2],vmin=1,vmax=5)
#     pyplot.savefig("GMC_"+str(i)+".png")
# #    subplot.set_title("GMC at zero age")
# #pyplot.title("Molecular cloud at time="+time.as_string_in(units.Myr))
# #pyplot.xlabel("x [pc]")
# #pyplot.ylabel("x [pc]")
# #pyplot.title("GMC at time="+time.as_string_in(units.Myr))
# #pyplot.savefig("GMC_"+str(i)+".png")


def main(N, Mtot, Rvir, t_end, dt):
    gas, converter = gas_sphere(N, Mtot, Rvir)

    hydro = Gadget2(converter, mode='periodic', number_of_workers=10)
    hydro.gas_particles.add_particles(gas)

    grid = convert_SPH_to_grid(hydro, (100,100,100), do_scale=True)

    print("Did Grid")
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
                      default=19.756 | units.AU,
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
