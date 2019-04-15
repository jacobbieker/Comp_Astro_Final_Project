from AccretionDisk import AccretionDisk
from SuperMassiveBlackHole import SuperMassiveBlackHole
from BinaryBlackHole import BinaryBlackHole
from amuse.datamodel import Particle, Particles, ParticlesSuperset
from amuse.couple.bridge import Bridge
from amuse.community.ph4.interface import ph4
import numpy as np
from amuse.lab import units, nbody_system, constants, Particles
from amuse.io import write_set_to_file


class BinaryBlackHolesWithAGN(object):

    def __init__(self, mass_of_central_black_hole, number_of_binaries, number_of_gas_particles, disk_mass_fraction, binaries_affect_disk=False,
                 radiative_transfer=False, timestep=0.1 | units.Myr, end_time = 5 | units.Myr, number_of_workers=1,
                 disk_powerlaw=1):
        self.smbh = SuperMassiveBlackHole(mass=mass_of_central_black_hole)
        self.inner_boundary = (self.smbh.radius) * 100
        self.outer_boundary = (self.smbh.radius) * 100000
        self.end_time = end_time
        self.converter = nbody_system.nbody_to_si(self.smbh.super_massive_black_hole.mass, self.outer_boundary) # Converter is wrong
        self.number_of_gas_particles = number_of_gas_particles
        self.disk = AccretionDisk(fraction_of_central_blackhole_mass=disk_mass_fraction,
                                  number_of_particles=self.number_of_gas_particles,
                                  disk_min=self.inner_boundary.value_in(self.outer_boundary.unit),
                                  disk_max=self.outer_boundary.value_in(self.outer_boundary.unit),
                                  end_of_disk=self.outer_boundary,
                                  number_of_workers=number_of_workers,
                                  converter=None,
                                  powerlaw=disk_powerlaw)

        print("Made Disk")
        exit()
        self.binaries = Particles()
        self.binaries_affect_disk = binaries_affect_disk
        self.number_of_binaries = number_of_binaries
        self.hydro_code = self.disk.hydro_code
        # Generate the binary locations and masses
        self.generate_binaries()

        # Now add them to a combined gravity code
        self.grav_code = ph4(self.converter, number_of_workers=number_of_workers)

        # Adding them together because not sure how to split the channels into different particle groups
        self.all_grav_particles = Particles() #ParticlesSuperset([self.smbh.super_massive_black_hole, self.binaries])
        self.all_grav_particles.add_particle(self.smbh.super_massive_black_hole)
        self.all_grav_particles.add_particles(self.binaries)
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

        # New particle superset of all particles in the sim
        # Initial Conditions
        self.disk.channel_from_hydro_to_grid.copy()
        all_sim_particles = ParticlesSuperset([self.grav_code.particles, self.disk.hydro_code.gas_particles])
        write_set_to_file(all_sim_particles, "{}_Binaries_{}_Gas_AGN_sim.hdf5".format(self.number_of_binaries, self.number_of_gas_particles), "amuse")

        while sim_time < end_time:
            sim_time += self.timestep

            self.grav_code.evolve_model(sim_time)
            self.bridge.evolve_model(sim_time)
            print('letsgo')

            self.channel_from_grav_to_binaries.copy()
            #self.disk.hydro_channel_to_particles.copy()
            self.disk.channel_from_hydro_to_grid.copy()

            # New particle superset of all particles in the sim
            all_sim_particles = ParticlesSuperset([self.grav_code.particles, self.disk.hydro_code.gas_particles])
            write_set_to_file(all_sim_particles, "{}_Binaries_{}_Gas_AGN_sim.hdf5".format(self.number_of_binaries, self.number_of_gas_particles), "amuse")

        self.grav_code.stop()
        self.disk.hydro_code.stop()


    def generate_binaries(self):
        # Now only use those outer particle positions to generate the binaries,
        # since nothing is within 2 radii of the black hole

        for _ in range(self.number_of_binaries):

            blackhole_masses = np.random.uniform(low=10, high=15, size=2)

            binary = BinaryBlackHole(blackhole_masses[0], blackhole_masses[1], self.smbh.super_massive_black_hole.mass,
                                     initial_outer_semi_major_axis=np.random.uniform(self.inner_boundary.value_in(self.outer_boundary.unit), self.outer_boundary.value_in(self.outer_boundary.unit), size=1)[0] | self.outer_boundary.unit,
                                     initial_outer_eccentricity=np.random.uniform(0, 0.99, size=1)[0],
                                     eccentricity=np.random.uniform(0.0, 0.99, size=1),
                                     inclincation=np.random.uniform(0.0, 180.0, size=1),
                                     )

            self.binaries.add_particles(binary.blackholes)

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
        self.bridge.add_system(self.grav_code, (self.hydro_code,))
        self.bridge.add_system(self.hydro_code, (self.grav_code,))

        return self.bridge
