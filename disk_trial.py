from amuse.lab import *
from amuse.community.gadget2.interface import Gadget2
from matplotlib import pyplot as plt




def gas_sphere(N, Mtot, Rvir):
	converter=nbody_system.nbody_to_si(Mtot, Rvir)
	gas = new_plummer_gas_model(N, convert_nbody=converter)
	return gas, converter








def main(N, Mtot, Rvir, t_end, dt):
	gas, converter = gas_sphere(N, Mtot, Rvir)

	hydro = Gadget2(converter)
	hydro.gas_particles.add_particles(gas)

	Etot_init = hydro.kinetic_energy \
	+ hydro.potential_energy + hydro.thermal_energy

	time = 0 | units.yr
	while time < t_end:
		time += dt
		hydro.evolve_model(time)
		print time


	write_set_to_file(hydro.particles, "hydro.h5", "hdf5")

	Ekin = hydro.kinetic_energy
	Epot = hydro.potential_energy
	Eth = hydro.thermal_energy
	Etot = Ekin + Epot + Eth
	Q = (Ekin+Eth)/Epot
	dE = (Etot_init-Etot)/Etot
	com = hydro.gas_particles.center_of_mass()

	print "T=", hydro.get_time(), "M=", hydro.gas_particles.mass.sum(),
	print "E= ", Etot, "Q= ", Q, "dE=", dE, "CoM=", com.in_(units.RSun)

	hydro.stop()

	x_gas = gas.x.value_in(units.AU)
	y_gas = gas.y.value_in(units.AU)
	z_gas = gas.z.value_in(units.AU)

	plt.scatter(x_gas, y_gas, z_gas)
	plt.savefig('disk_scatter.png')
	plt.show()

def new_option_parser():
	from amuse.units.optparse import OptionParser
	result = OptionParser()

	result.add_option("-N", dest="N", type="int",
					  default = 10000,
					  help="no. of particles [%default]")
	result.add_option("--M_tot", unit = units.MSun, dest="Mtot", type="float",
					  default = 10 | units.MSun,
					  help="Total disk mass [%default]")
	result.add_option("--Rvir", unit = units.AU, dest="Rvir", type="float",
					  default = 100 | units.AU,
					  help="Radius of initial sphere [%default]")
	result.add_option("--t_end", unit = units.yr, dest="t_end", type="float",
					  default = 1000 | units.yr,
					  help="End time of simulation [%default]")
	result.add_option("--dt", unit = units.yr, dest="dt", type="float",
					  default = 10 | units.yr,
					  help="Hydro timestep [%default]")
	return result


if __name__ in ("__main__"):
    o, arguments = new_option_parser().parse_args()
    print o.__dict__
    main(**o.__dict__)
