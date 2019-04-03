from AccretionDisk import AccretionDisk
from SuperMassiveBlackHole import SuperMassiveBlackHole
from BinaryBlackHole import BinaryBlackHole
from amuse.ic.plummer import new_plummer_model
from amuse.datamodel import Particle, Particles
from amuse.couple.bridge import Bridge
from amuse.units import units


class System(object):

    def __init__(self, mass_of_central_black_hole, number_of_binaries, disk_mass_fraction, binaries_affect_disk=False,
                 radiative_transfer=False, timestep=0.1 | units.Myr, converter = None):
        self.supermassive_blackhole = SuperMassiveBlackHole(mass=mass_of_central_black_hole)
        self.disk = AccretionDisk(fraction_of_central_blackhole_mass=disk_mass_fraction)
        self.binaries = Particles()
        self.binaries_affect_disk = binaries_affect_disk
        self.number_of_binaries = number_of_binaries
        self.grav_code = None
        self.hydro_code = None
        self.converter = converter
        # Generate the binary locations and masses



        self.bridge = self.create_bridges(timestep)

        raise NotImplementedError

    def evolve_model(self, end_time):

        raise NotImplementedError

    def generate_binaries(self, method):

        binary_locations = method(self.number_of_binaries, convert_nbody=self.converter)

        com = binary_locations.center_of_mass()
        # determine bump's local velocity
        inner_particles = binary_locations.select(lambda r: (com-r).length()<self.supermassive_blackhole.radius,["position"])



        raise NotImplementedError

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
