from AccretionDisk import AccretionDisk
from SuperMassiveBlackHole import SuperMassiveBlackHole
from BinaryBlackHole import BinaryBlackHole
from amuse.datamodel import Particle, Particles, ParticlesSuperset
from amuse.couple.bridge import Bridge
from amuse.community.ph4.interface import ph4
from amuse.community.hermite0.interface import Hermite
import numpy as np
from amuse.lab import units, nbody_system, constants, Particles
from amuse.io import write_set_to_file


class BinaryBlackHolesWithAGN(object):

    def __init__(self, mass_of_central_black_hole, number_of_binaries, number_of_gas_particles, disk_mass_fraction, binaries_affect_disk=False,
                 radiative_transfer=False, timestep=0.1 | units.Myr, end_time = 5 | units.Myr, number_of_hydro_workers=1, number_of_grav_workers=1, steps_of_inclination = 18,
                 disk_powerlaw=1):
        self.smbh = SuperMassiveBlackHole(mass=mass_of_central_black_hole)
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
        write_set_to_file(self.disk.gas_particles, "Initial_AccretionDisk_SMBH_Mass_{}_MSun.hdf5".format(self.smbh.super_massive_black_hole.mass.value_in(units.MSun)), "hdf5")

        self.binaries = Particles()
        self.binaries_affect_disk = binaries_affect_disk
        self.number_of_binaries = number_of_binaries
        self.hydro_code = self.disk.hydro_code
        # Generate the binary locations and masses
        self.all_grav_particles = self.smbh.super_massive_black_hole #Particles()
        #self.all_grav_particles.add_particle(self.smbh.super_massive_black_hole)
        #self.all_grav_particles.add_particles(self.binaries)
        self.gravity_converter = nbody_system.nbody_to_si(self.all_grav_particles.mass, self.all_grav_particles.radius)

        # Now add them to a combined gravity code
        self.grav_code = ph4(self.gravity_converter, number_of_workers=number_of_grav_workers)
        # Adding them gravity
        self.grav_code.particles.add_particle(self.smbh.super_massive_black_hole)
        # Adding them together because not sure how to split the channels into different particle groups


        # Channels to update the particles here
        # self.channel_from_grav_to_binaries = self.grav_code.particles.new_channel_to(self.smbh.super_massive_black_hole)
        # self.channel_from_binaries_to_grav = self.smbh.super_massive_black_hole.new_channel_to(self.grav_code.particles)

        self.timestep = timestep
        self.bridge = self.create_bridges(timestep)
        self.generate_binaries()
        self.add_bridge(self.binary_codes)
        self.evolve_model(self.end_time)

    def evolve_model(self, end_time):

        sim_time = 0. | self.end_time.unit

        # New particle superset of all particles in the sim
        # Initial Conditions
        all_sim_particles = ParticlesSuperset([self.grav_code.particles, self.disk.hydro_code.gas_particles])
        write_set_to_file(all_sim_particles, "{}_Binaries_{}_Gas_AGN_sim.hdf5".format(self.number_of_binaries, self.number_of_gas_particles), "amuse")

        while sim_time < end_time:
            sim_time += self.timestep
            self.bridge.evolve_model(sim_time)
            print('Time: {}'.format(sim_time.value_in(units.yr)), flush=True)

            #self.channel_from_grav_to_binaries.copy()
            self.disk.hydro_channel_to_particles.copy()
            for channel in self.binary_code_from_channels:
                channel.copy()

            # New particle superset of all particles in the sim
            all_sim_particles = ParticlesSuperset([self.grav_code.particles, self.disk.hydro_code.gas_particles])
            write_set_to_file(all_sim_particles, "{}_Binaries_{}_Gas_AGN_sim.hdf5".format(self.number_of_binaries, self.number_of_gas_particles), "amuse")

        self.grav_code.stop()
        self.disk.hydro_code.stop()


    def generate_binaries(self):
        # Now only use those outer particle positions to generate the binaries,
        # since nothing is within 2 radii of the black hole
        initial_outer_semi_major_axis = np.linspace(self.inner_boundary.value_in(self.outer_boundary.unit), self.outer_boundary.value_in(self.outer_boundary.unit), number_of_binaries//steps_of_inclination)
        inclination_list = np.linspace(0,180, self.steps_of_inclination)
        #
        #
        #
        # need to feix the iteration
        #
        #
        # 

        for i in range(self.number_of_binaries):

            blackhole_masses = [30,30]

            binary = BinaryBlackHole(blackhole_masses[0], blackhole_masses[1], self.smbh.super_massive_black_hole.mass,
                                     initial_outer_semi_major_axis=initial_outer_semi_major_axis[i] | self.outer_boundary.unit,
                                     eccentricity=0.6,
                                     inclincation=np.random.uniform(0.0, 180.0, size=1),
                                     )

            smbh_and_binary = Particles()
            smbh_and_binary.add_particle(self.smbh.super_massive_black_hole)
            smbh_and_binary.add_particles(binary.blackholes)
            binary_converter = nbody_system.nbody_to_si(smbh_and_binary.mass.sum(), smbh_and_binary.virial_radius())
            one_binary_code = Hermite(binary_converter)
            one_binary_code.particles.add_particles(binary.blackholes)

            self.binary_codes.append(one_binary_code)
            self.binaries.add_particles(binary.blackholes)

            channel_from_grav_to_binaries = one_binary_code.particles.new_channel_to(binary.blackholes)
            channel_from_binaries_to_grav = binary.blackholes.new_channel_to(one_binary_code.particles)
            self.binary_code_from_channels.append(channel_from_grav_to_binaries)
            self.binary_code_to_channels.append(channel_from_binaries_to_grav)

    def create_bridges(self, timestep=0.1 | units.Myr):
        """
        Bridge between SMBH and disk one way
        Bridge between disk and Binaries, one way
        SMBH and Binaries are in the same gravity, no bridge
        Possibly bridge from binaries to the gas

        :return:
        """

        self.bridge = Bridge(use_threading=True)
        self.bridge.timestep = timestep
        #self.bridge.add_system(self.grav_code, (self.hydro_code,))
        self.bridge.add_system(self.hydro_code, (self.grav_code, ))

        return self.bridge

    def add_bridge(self, code):
        self.bridge.add_system(self.grav_code, set(code))
        self.bridge.add_system(self.hydro_code, set(code))
