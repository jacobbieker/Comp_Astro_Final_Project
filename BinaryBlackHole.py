from __future__ import division, print_function
import numpy
from amuse.datamodel import Particles, Particle
from amuse.ext.solarsystem import get_position
from amuse.units import units, constants


class BinaryBlackHole(object):

    def __init__(self, central_blackhole_mass, initial_outer_semi_major_axis, mass_one=30. | units.MSun, mass_two=30. | units.MSun,
                 initial_outer_eccentricity=0.6, inner_eccentricity=0.6, inclination=0.0,
                 orbital_fraction_timestep=0.5):
        """
        This model is to generate a binary black hole system when given an initial Particle to split

        Create binaries, and then distribute them according to Plummer sphere

        Have it set center of mass velocity and center of mass location

        Keep binaries by themselves

        Put Nbody in bridge

        Have it lose energy as it goes through the disk

        Put that in the bridge so that it loses energy in conjunction with the hydro code

        Put in dynamical friction term in bridge for losing energy then, so it does it

        But need to figure out dynamical friction equation for that

        Make a protoplanetary disk -> AGN disk is the same, just larger

        """
        self.central_blackhole = Particle(1)
        self.central_blackhole.mass = central_blackhole_mass

        self.blackholes = Particles(2)
        self.blackholes[0].mass = mass_one
        self.blackholes[1].mass = mass_two
        self.total_mass = self.blackholes.mass.sum()

        self.initial_outer_semi_major_axis = initial_outer_semi_major_axis
        self.initial_outer_eccentricity = initial_outer_eccentricity

        # Calculating the hill radius to act as a maximum boundary for the initial inner semi major axis and therefore the inner orbital period
        self.hill_radius = self.get_hill_radius(self.initial_outer_semi_major_axis, self.initial_outer_eccentricity,
                                                self.total_mass, self.central_blackhole.mass)
        self.binary_max_orbital_period = self.get_orbital_period(0.5 * self.hill_radius, self.total_mass)
        self.binary_min_orbital_period = self.get_orbital_period(
            1000000 * self.get_schwarzschild_radius(self.blackholes[0].mass), self.total_mass)
        # Restricting the initial inner orbital period
        self.orbital_period = numpy.random.uniform(self.binary_min_orbital_period.value_in(units.yr),
                                                   self.binary_max_orbital_period.value_in(units.yr), size=1) | units.yr

        self.inner_semi_major_axis = self.get_semi_major_axis(self.total_mass, self.orbital_period)
        self.inner_eccentricity = inner_eccentricity
        self.inclination = inclination

        self.timestep = self.orbital_period * orbital_fraction_timestep

        # Create the binary blackholes and set them in orbit around each other
        binary_position, binary_velocity = get_position(self.blackholes[0].mass, self.blackholes[1].mass,
                                                        self.inner_eccentricity, self.inner_semi_major_axis,
                                                        180, self.inclination, 180, 0, self.timestep)
        self.blackholes[1].position = binary_position
        self.blackholes[1].velocity = binary_velocity

        self.blackholes_distance = (self.blackholes[0].position - self.blackholes[1].position).length()
        self.minimum_distance = 100 * self.get_schwarzschild_radius(self.blackholes[0].mass)
        # self.merged_blackhole = Particle()

        self.blackholes.move_to_center()
        # Set the (inner binary) center of mass and center of mass velocity in orbit around the supermassive black hole
        self.initial_outer_semi_major_axis, self.initial_outer_eccentricity = self.set_in_orbit_around_central_blackhole(
            central_blackhole_mass,
            self.initial_outer_semi_major_axis,
            self.initial_outer_eccentricity,
            time_to_advance=(self.orbital_period/2.))

    def get_orbital_period(self, orbital_separation, total_mass):
        return 2 * numpy.pi * (orbital_separation ** 3 / (constants.G * total_mass)).sqrt()

    def get_semi_major_axis(self, total_mass, orbital_period):
        return (constants.G * total_mass * orbital_period ** 2 / (4 * numpy.pi ** 2)) ** (1. / 3.)

    def get_hill_radius(self, semi_major_axis, eccentricity, total_binary_mass, central_blackhole_mass):
        return semi_major_axis * (1 - eccentricity) * (total_binary_mass / (3 * central_blackhole_mass)) ** (1. / 3.)

    def get_schwarzschild_radius(self, mass):
        return (2 * constants.G * mass) / (constants.c ** 2)

    def set_center_of_mass_position(self, new_center_of_mass):
        """
        Sets the center of mass position of the binary system
        :param new_center_of_mass: center of mass needed
        :return:
        """

        try:
            self.blackholes[0].position -= new_center_of_mass
            self.blackholes[1].position -= new_center_of_mass
        except:
            # Need to convert to a list with units, instead of a list of elements with units
            for index, element in enumerate(new_center_of_mass):
                new_center_of_mass[index] = element.value_in(units.m)
            new_center_of_mass = new_center_of_mass | units.m
            self.blackholes[0].position -= new_center_of_mass
            self.blackholes[1].position -= new_center_of_mass

    def set_center_of_mass_velocity(self, new_center_of_mass_velocity):
        """
        Sets the center of mass velocity of the binary system
        :param new_center_of_mass_velocity:
        :return:
        """
        try:
            self.blackholes[0].velocity -= new_center_of_mass_velocity
            self.blackholes[1].velocity -= new_center_of_mass_velocity
        except:
            # Need to convert to a list with units, instead of a list of elements with units
            for index, element in enumerate(new_center_of_mass_velocity):
                new_center_of_mass_velocity[index] = element.value_in(units.m / units.s)
            new_center_of_mass = new_center_of_mass_velocity | units.m / units.s
            self.blackholes[0].velocity -= new_center_of_mass
            self.blackholes[1].velocity -= new_center_of_mass

    def set_binary_location_and_velocity(self, center_of_mass, center_of_mass_velocity):
        self.set_center_of_mass_position(center_of_mass)
        self.set_center_of_mass_velocity(center_of_mass_velocity)

    def set_in_orbit_around_central_blackhole(self, central_blackhole_mass, semi_major_axis,
                                              eccentricity, mean_anomaly=0, inclination=0,
                                              argument_of_perhilion=0, longitude_of_ascending_node=0,
                                              time_to_advance=5 | units.day):
        """
        Sets the binary's center of mass in orbit around the supermassive central blackhole
        :param central_blackhole:
        :param eccentricity:
        :param semi_major_axis:
        :param mean_anomaly:
        :param inclination:
        :param argument_of_perhilion:
        :param longitude_of_ascending_node:
        :param time_to_advance:
        :return:
        """

        binary_orbital_position, binary_orbital_velocity = get_position(central_blackhole_mass,
                                                                        self.blackholes.mass.sum(),
                                                                        ecc=eccentricity,
                                                                        semi=semi_major_axis,
                                                                        mean_anomaly=mean_anomaly,
                                                                        incl=inclination,
                                                                        argument=argument_of_perhilion,
                                                                        longitude=longitude_of_ascending_node,
                                                                        delta_t=time_to_advance)
        # Get the center of mass position and velocity
        self.set_binary_location_and_velocity(binary_orbital_position, binary_orbital_velocity)

        return semi_major_axis, eccentricity
