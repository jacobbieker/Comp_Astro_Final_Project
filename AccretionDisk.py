from amuse.community.gadget2.interface import Gadget2
from amuse.datamodel import Particles
from amuse.units import units, constants
import numpy as np


class AccretionDisk(object):
    """
    This class creates an accretion disk around a central particle, ideally a black hole

    This class encapsulates the Hydrodynamics code, built around Gadget2, and enables some convience functions, such
    as getting the information needed for plotting, etc.

    Could also be extended with Radiative Transport if needd
    """

    def __init__(self, number_of_particles=10000, mode='normal', converter=None, number_of_workers=1):
        self.code = Gadget2(unit_converter=converter, mode=mode, number_of_workers=number_of_workers)
        self.number_of_particles = number_of_particles
        self.converter = converter
        self.gas_particles = self.make_disk(number_of_particles)

    def distribute_particles(self):
        raise NotImplementedError

    def make_disk(self, number_of_particles, center_of_disk=(0,0,0)):
        """
        Makes the accretion disk around the center of disk
        :param center_of_disk: Center of disk coordinates, in (x,y,z) format
        :param number_of_particles: Number of particles to use
        :return: The particles in the disk, as a Particles set
        """



        raise NotImplementedError

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
            x,y,z = np.indices((num_points+1,num_points+1,num_points+1))
            x=(x.flatten()-num_points/2.)/num_points
            y=(y.flatten()-num_points/2.)/num_points
            z=(z.flatten()-num_points/2.)/num_points

            vx=0.*x
            vy=0.*x
            vz=0.*x
        else:
            x,y = np.indices((num_points+1,num_points+1))
            x=(x.flatten()-num_points/2.)/num_points
            y=(y.flatten()-num_points/2.)/num_points
            z=z_plane*np.ones(x.shape)

            vx=0.*x
            vy=0.*x
            vz=0.*x

        x=units.AU(x)
        y=units.AU(y)
        z=units.AU(z)
        vx=units.kms(vx)
        vy=units.kms(vy)
        vz=units.kms(vz)

        rho,rhovx,rhovy,rhovz,rhoe=self.code.get_hydro_state_at_point(x,y,z,vx,vy,vz)
        rho=rho.reshape((num_points+1,num_points+1, num_points+1))
        return rho

    def get_total_energy(self):
        """
        Gets the total energy of the hydro system
        :return:
        """

        return self.code.kinetic_energy + self.code.potential_energy + self.code.thermal_energy
