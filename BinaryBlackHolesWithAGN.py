from AccretionDisk import AccretionDisk
from SuperMassiveBlackHole import SuperMassiveBlackHole
from BinaryBlackHole import BinaryBlackHole
from amuse.ic.plummer import new_plummer_model
from amuse.datamodel import Particle, Particles
from amuse.couple.bridge import Bridge
from amuse.community.huayno.interface import Huayno
from amuse.units import units
import numpy as np


class BinaryBlackHolesWithAGN(object):

    def __init__(self, mass_of_central_black_hole, number_of_binaries, disk_mass_fraction, binaries_affect_disk=False,
                 radiative_transfer=False, timestep=0.1 | units.Myr, converter=None, number_of_workers=1,
                 disk_powerlaw=1):
        self.smbh = SuperMassiveBlackHole(mass=mass_of_central_black_hole)
        self.inner_boundary = (self.smbh.radius.value_in(units.parsec)) * 2
        self.outer_buondary = (self.smbh.radius.value_in(units.parsec)) * 100
        self.disk = AccretionDisk(fraction_of_central_blackhole_mass=disk_mass_fraction,
                                  disk_min=self.inner_boundary,
                                  disk_max=self.outer_buondary,
                                  number_of_workers=number_of_workers,
                                  converter=converter,
                                  powerlaw=disk_powerlaw)
        self.binaries = Particles()
        self.binaries_affect_disk = binaries_affect_disk
        self.number_of_binaries = number_of_binaries
        self.hydro_code = self.disk.hydro_code
        self.converter = converter
        # Generate the binary locations and masses
        self.generate_binaries(new_plummer_model)

        # Now add them to a combined gravity code
        self.grav_code = Huayno(converter, number_of_workers=number_of_workers)

        # Adding them together because not sure how to split the channels into different particle groups
        self.all_grav_particles = self.smbh.super_massive_black_hole + self.binaries
        self.grav_code.add_particles(self.all_grav_particles)

        # Channels to update the particles here
        self.channel_from_grav_to_binaries = self.grav_code.particles.new_channel_to(self.all_grav_particles)
        self.channel_from_binaries_to_grav = self.all_grav_particles.new_channel_to(self.grav_code.particles)

        self.timestep = timestep
        self.bridge = self.create_bridges(timestep)

    def evolve_model(self, end_time):

        sim_time = 0. | end_time.units

        while sim_time < end_time:
            sim_time += self.timestep

            self.bridge.evolve_model(sim_time)

            self.channel_from_grav_to_binaries.copy()
            self.disk.hydro_channel_to_particles.copy()


    def generate_binaries(self, method):
        binary_locations = method(self.number_of_binaries, convert_nbody=self.converter)

        com = binary_locations.center_of_mass()
        # determine bump's local velocity
        outer_particles = binary_locations.select(lambda r: (com - r).length() > 2 * self.smbh.radius,
                                                  ["position"])

        # Now only use those outer particle positions to generate the binaries,
        # since nothing is within 2 radii of the black hole

        for i in self.number_of_binaries:
            blackhole_masses = np.random.uniform(low=10, high=15, size=2)
            binary = BinaryBlackHole(blackhole_masses[0], blackhole_masses[1],
                                     orbital_period=1 | units.yr,
                                     eccentricity=np.random.uniform(0.0, 0.99, size=1),
                                     inclincation=np.random.uniform(0.0, 90.0, size=1))

            binary.set_in_orbit_around_central_blackhole(central_blackhole=self.smbh.super_massive_black_hole,
                                                         eccentricity=np.random.uniform(0.0, 0.99, size=1),
                                                         inclination=np.random.uniform(0.0, 90.0, size=1),
                                                         semi_major_axis=np.random.uniform(self.inner_boundary, self.outer_buondary, size=1) | units.parsec,
                                                         )
            self.binaries.add_particles(binary)

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
