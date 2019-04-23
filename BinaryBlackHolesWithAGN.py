from __future__ import print_function
from AccretionDisk import AccretionDisk
from SuperMassiveBlackHole import SuperMassiveBlackHole
from BinaryBlackHole import BinaryBlackHole
from amuse.couple.bridge import Bridge
from amuse.community.huayno.interface import Huayno
import numpy as np
from amuse.lab import units, nbody_system, constants, Particles
from amuse.io import write_set_to_file


class SuperMassiveBlackHolePotential(object):
    def __init__(self,R, M):
        self.radius=R
        self.mass=M
        self.alpha = 1.0

    def get_gravity_at_point(self,eps,x,y,z):
        r2=x**2+y**2+z**2
        r=r2.sqrt()
        m=self.mass*(r/self.radius)**self.alpha
        fr=constants.G*m/r2
        ax=-fr*x/r
        ay=-fr*y/r
        az=-fr*z/r
        return ax,ay,az

    def get_potential_at_point(self,eps,x,y,z):
        r=(x**2+y**2+z**2)**0.5
        c=constants.G*self.mass/self.radius**self.alpha
        phi=c/(self.alpha-1)*(r**(self.alpha-1)-self.radius**(self.alpha-1))
        return phi


class BinaryBlackHolesWithAGN(object):

    def __init__(self, mass_of_central_black_hole, number_of_binaries, number_of_gas_particles, disk_mass_fraction, binaries_affect_disk=False,
                 radiative_transfer=False, timestep=0.1 | units.Myr, end_time = 5 | units.Myr, number_of_hydro_workers=1, number_of_grav_workers=1, steps_of_inclination = 18,
                 disk_powerlaw=1):
        self.smbh = SuperMassiveBlackHole(mass=mass_of_central_black_hole)
        self.smbh_potential = SuperMassiveBlackHolePotential(M=self.smbh.super_massive_black_hole.mass, R=self.smbh.radius)
        self.inner_boundary = self.smbh.radius * 100
        self.outer_boundary = self.smbh.radius * 100000
        self.steps_of_inclination = steps_of_inclination
        self.end_time = end_time
        self.binary_codes = []
        self.binary_code_from_channels = []
        self.binary_code_to_channels = []
        self.number_of_gas_particles = number_of_gas_particles
        self.disk_converter = nbody_system.nbody_to_si(self.smbh.super_massive_black_hole.mass, self.inner_boundary)
        self.gadget_converter = nbody_system.nbody_to_si(disk_mass_fraction*self.smbh.super_massive_black_hole.mass, self.outer_boundary)
        self.disk = AccretionDisk(fraction_of_central_blackhole_mass=disk_mass_fraction,
                                  number_of_particles=self.number_of_gas_particles,
                                  disk_min=1.,
                                  disk_max=self.outer_boundary/self.inner_boundary,
                                  number_of_workers=number_of_hydro_workers,
                                  gadget_converter=self.gadget_converter,
                                  disk_converter=self.disk_converter,
                                  powerlaw=disk_powerlaw,
                                  end_time=self.end_time)

        self.binaries = Particles()
        self.merged_blackholes = Particles()
        self.minimum_distance = 100 * self.get_schwarzschild_radius(30 | units.MSun)
        self.binaries_affect_disk = binaries_affect_disk
        self.number_of_binaries = number_of_binaries
        self.hydro_code = self.disk.hydro_code
        # Generate the binary locations and masses
        self.all_grav_particles = Particles()
        self.generate_binaries()
        self.gravity_converter = nbody_system.nbody_to_si(self.all_grav_particles.mass.sum(), self.all_grav_particles.virial_radius())

        # Now add them to a combined gravity code
        self.grav_code = Huayno(self.gravity_converter, number_of_workers=number_of_grav_workers)
        self.grav_code.timestep = 100 | units.yr

        # Adding them gravity
        self.grav_code.particles.add_particles(self.all_grav_particles)

        # Channels to update the particles here
        self.channel_from_grav_to_binaries = self.grav_code.particles.new_channel_to(self.all_grav_particles)
        self.channel_from_binaries_to_grav = self.all_grav_particles.new_channel_to(self.grav_code.particles)

        self.timestep = timestep
        self.bridge = self.create_bridges(timestep)
        self.evolve_model(self.end_time)

    def evolve_model(self, end_time):

        sim_time = 0. | self.end_time.unit


        while sim_time < end_time:
            # Now extract information such as inclination to each other and the disk

            # New particle superset of all particles in the simulation
            all_sim_particles = self.bridge.particles
            print(len(all_sim_particles))
            # Now extract information
            write_set_to_file(all_sim_particles, "Particles_{}_Binaries_{}_Gas_AGN_sim.hdf5".format(self.number_of_binaries, self.number_of_gas_particles), "hdf5")
            # Now evolve the total model of hydro and gravity
            sim_time += self.timestep
            self.bridge.evolve_model(sim_time)
            print('Time: {}'.format(sim_time.value_in(units.yr)))

            self.channel_from_grav_to_binaries.copy()
            self.disk.hydro_channel_to_particles.copy()

        self.grav_code.stop()
        self.disk.hydro_code.stop()


    def generate_binaries(self):
        """
        Generate a number of blackhole binaries with random initial outer semi major axis and inclination within the boundaries
        """
        blackhole_masses = [30,30]
        for _ in range(self.number_of_binaries):
            # Make sure to generate the binaries randomly but within the disk radius and not too close the SMBH
            initial_outer_semi_major_axis = np.random.uniform(self.inner_boundary.value_in(self.outer_boundary.unit), self.outer_boundary.value_in(self.outer_boundary.unit), 1)[0]
            initial_outer_eccentricity = np.random.uniform(0, 180, 1)[0]
            binaries = BinaryBlackHole(blackhole_masses[0], blackhole_masses[1], self.smbh.super_massive_black_hole.mass,
                                       initial_outer_semi_major_axis= initial_outer_semi_major_axis | (self.outer_boundary.unit),
                                       initial_outer_eccentricity=0.6,
                                       inner_eccentricity=0.6,
                                       inclination=initial_outer_eccentricity,
                                       )

            # Add the particles in the gravity particles
            self.all_grav_particles.add_particles(binaries.blackholes)
            self.binaries.add_particles(binaries.blackholes)

    def create_bridges(self, timestep=0.1 | units.Myr):
        """
        Bridge between SMBH potential and disk one way (smbh affects disk)
        Bridge between SMBH potential and binaries one way (smbh affects binaries)
        Bridge between disk and binaries one way (disk affects binaries)
        :return:
        """

        self.bridge = Bridge(use_threading=True, verbose=True)
        self.bridge.timestep = timestep
        self.bridge.add_system(self.grav_code, (self.smbh_potential,))
        self.bridge.add_system(self.hydro_code, (self.grav_code, self.smbh_potential ))
        # self.bridge.add_system(self.grav_code, (self.hydro_code,))

        return self.bridge









    self.blackholes_distance = (self.blackholes[0].position - self.blackholes[1].position).length()
    merge_condition = self.set_merge_conditions(self.blackholes_distance, self.minimum_distance)

    if merge_condition:
        print('binaries merged')
        self.merge_blackholes()


    def set_merge_conditions(self, blackholes_distance, minimum_distance):
        merge_condition = blackholes_distance < minimum_distance
        return merge_condition

    def merge_blackholes(self, merging_blackholes, fraction_of_total_mass=0.95):
        """
        Merges the binary particles into a single particle after merging
        :param fraction_of_total_mass: sets the mass of the merged blackhole
        :return:
        """

        merged_blackhole_location = merging_blackholes.center_of_mass_position()
        merged_blackhole_velocity = merging_blackholes.center_of_mass_velocity()
        # Set the initial position and velocity of the merged_blackholes to be the same as the center of mass position and velocity
        merged_blackhole = Particle()
        merged_blackhole.mass = fraction_of_total_mass * self.total_mass
        merged_blackhole.position = merged_blackhole_location
        merged_blackhole.velocity = merged_blackhole_velocity

        self.grav_code.particles.remove_particles(self.blackholes[0], self.blackholes[1])
        self.grav_code.particles.add_particle(merged_blackhole)
