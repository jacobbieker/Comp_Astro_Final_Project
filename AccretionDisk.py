from Gadget2_Gravity import Gadget2_Gravity
from amuse.units import units, constants
from amuse.units.generic_unit_converter import ConvertBetweenGenericAndSiUnits
import numpy as np
from amuse.ext.protodisk import ProtoPlanetaryDisk


class AccretionDisk(object):
    """
    This class creates an accretion disk around a central particle, ideally a black hole

    This class encapsulates the Hydrodynamics code, built around Gadget2, and enables some convience functions, such
    as getting the information needed for plotting, etc.

    Could also be extended with Radiative Transport if needd

    E_dot, I_dot, A_dot

    Start with different hill radii

    Make 30 solar mass, half hill radius apart, 0.6 eccentricity, set apocenter as half hill radius

    Set in different orbits around the black hole

    Prevent randomness

    Use smallN for all of the


    """

    def __init__(self, number_of_particles=100000, gadget_converter=None, disk_converter=None,
                 number_of_workers=1, disk_min=1, disk_max=1e4, fraction_of_central_blackhole_mass=0.1,
                 powerlaw=1e-2, end_time=5 | units.Myr):
        self.gadget_converter = gadget_converter
        self.disk_converter = disk_converter
        self.gen_convert = ConvertBetweenGenericAndSiUnits(constants.c, units.s)
        self.code = Gadget2_Gravity(self.gadget_converter, number_of_workers=number_of_workers)
        self.code.parameters.time_max = 2*self.gen_convert.to_generic(end_time)
        self.number_of_particles = number_of_particles
        self.disk_min = disk_min
        self.disk_max = disk_max
        self.fraction_of_central_blackhole_mass = fraction_of_central_blackhole_mass
        self.powerlaw = powerlaw
        self.gas_particles = self.make_disk(number_of_particles)
        self.code.gas_particles.add_particles(self.gas_particles)
        self.hydro_channel_to_particles = self.code.gas_particles.new_channel_to(self.gas_particles)
        self.particles_channel_to_hydro = self.gas_particles.new_channel_to(self.code.gas_particles)

    def make_disk(self, number_of_particles):
        """
        Makes the accretion disk around the center of disk
        :param center_of_disk: Center of disk coordinates, in (x,y,z) format
        :param number_of_particles: Number of particles to use
        :return: The particles in the disk, as a Particles set
        """
        gas_particles = ProtoPlanetaryDisk(number_of_particles,
                                           convert_nbody=self.disk_converter,
                                           densitypower=self.powerlaw,
                                           Rmin=self.disk_min,
                                           Rmax=self.disk_max,
                                           q_out=1.0,
                                           discfraction=self.fraction_of_central_blackhole_mass).result
        gas_particles.move_to_center()
        return gas_particles

    @property
    def hydro_code(self):
        """
        Returns the direct Hydro code for easier manipulation
        :return:
        """
        return self.code

    def get_density_map(self, num_points=1000, z_plane=None):
        """
        Returns the density map of the gas, mostly for visualization, sampled num_points times in
        the x,y,and z directions if z_plane = None, or in x and y at the specific z value if otherwise

        :return:
        """

        if z_plane is None:
            x, y, z = np.indices((num_points + 1, num_points + 1, num_points + 1))
            x = (x.flatten() - num_points / 2.) / num_points
            y = (y.flatten() - num_points / 2.) / num_points
            z = (z.flatten() - num_points / 2.) / num_points

            vx = 0. * x
            vy = 0. * x
            vz = 0. * x
        else:
            x, y = np.indices((num_points + 1, num_points + 1))
            x = (x.flatten() - num_points / 2.) / num_points
            y = (y.flatten() - num_points / 2.) / num_points
            z = z_plane * np.ones(x.shape)

            vx = 0. * x
            vy = 0. * x
            vz = 0. * x

        x = units.AU(x)
        y = units.AU(y)
        z = units.AU(z)
        vx = units.kms(vx)
        vy = units.kms(vy)
        vz = units.kms(vz)

        rho, rhovx, rhovy, rhovz, rhoe = self.code.get_hydro_state_at_point(x, y, z, vx, vy, vz)
        rho = rho.reshape((num_points + 1, num_points + 1, num_points + 1))
        return rho

    def get_total_energy(self):
        """
        Gets the total energy of the hydro system
        :return:
        """

        return self.code.kinetic_energy + self.code.potential_energy + self.code.thermal_energy
