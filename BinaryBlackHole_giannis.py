import numpy as np
from amuse.datamodel import Particles
from amuse.ext.solarsystem import get_position
from amuse.units import units, constants


class BinaryBlackHole(object):

    def __init__(self, mass_one, mass_two, orbital_period, eccentricity=0.0, inclincation=0.0,
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
        self.orbital_period = orbital_period
        self.eccentricity = eccentricity
        self.inclincation = inclincation
        self.blackholes = Particles(2)
        self.blackholes.mass = [mass_one, mass_two] | units.MSun
        self.total_mass = self.blackholes.mass.sum()
        self.semi_major_axis = self.get_semi_major_axis()
        self.timestep = self.orbital_period * orbital_fraction_timestep

        binary_position, binary_velocity = get_position(self.blackholes[0].mass, self.blackholes[1].mass,
                                                        self.eccentricity, self.semi_major_axis,
                                                        180, self.inclincation, 180, 0, self.timestep)
        self.blackholes[1].position = binary_position
        self.blackholes[1].velocity = binary_velocity
        self.blackholes.move_to_center()

    def particles(self):
        return self.blackholes

    def get_orbital_period(self, orbital_separation, total_mass):
        return 2 * np.pi * (orbital_separation ** 3 / (constants.G * total_mass)).sqrt()

    def get_semi_major_axis(self):
        return (constants.G * self.total_mass * self.orbital_period ** 2 / (4 * np.pi ** 2)) ** (1. / 3)

    def set_center_of_mass(self, new_center_of_mass):
        """
        Sets the center of mass of the binary system
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
            self.blackholes[0].position -= new_center_of_mass_velocity
            self.blackholes[1].position -= new_center_of_mass_velocity
        except:
            # Need to convert to a list with units, instead of a list of elements with units
            for index, element in enumerate(new_center_of_mass_velocity):
                new_center_of_mass_velocity[index] = element.value_in(units.m / units.s)
            new_center_of_mass = new_center_of_mass_velocity | units.m / units.s
            self.blackholes[0].velocity -= new_center_of_mass
            self.blackholes[1].velocity -= new_center_of_mass

    def set_binary_location_and_velocity(self, center_of_mass, center_of_mass_velocity):
        self.set_center_of_mass(center_of_mass)
        self.set_center_of_mass_velocity(center_of_mass_velocity)

    def set_merge_conditions(self, minimum_distance = 350 | units.km):
        blackholes_distance = (self.blackholes[0].position - self.blackholes[1].position).length()
        merge_condition = blackholes_distance < minimum_distance
        return merge_condition

    def set_in_orbit_around_central_blackhole(self, central_blackhole, eccentricity,
                                              semi_major_axis, mean_anomaly=0, inclination=0,
                                              argument_of_perhilion=0, longitude_of_ascending_node=0,
                                              time_to_advance=5 | units.day,
                                              merge_condition):
        """
        Sets the binary's orbit around the supermassive central blackhole
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

        binary_orbital_position, binary_orbital_velocity = get_position(central_blackhole.mass,
                                                                       self.blackholes.mass.sum(),
                                                                       ecc=eccentricity,
                                                                       semi=semi_major_axis,
                                                                       mean_anomaly=mean_anomaly,
                                                                       incl=inclination,
                                                                       argument=argument_of_perhilion,
                                                                       longitude=longitude_of_ascending_node,
                                                                       delta_t=time_to_advance)

        self.set_binary_location_and_velocity(binary_orbital_position, binary_orbital_velocity)



        self.set_merge_conditions()
        try:
            self.merge_particles(merge_condition)
        except:
            pass





    def merged_blackhole_attributes(self, fraction_of_total_mass = 0.95, total_mass):
        """
        Merges the binary particles into a single particle, with the
        :return:
        """
        # 1) Get set_binary_location_and_velocity value
        # 2) Add the merged particle to merged_blackholes
        # 3) Remove the binary blackholes from blackholes set
        # 4) Bridge the new particle set to SMBH gravity
        self.merged_blackhole = Particle()
        merged_blackhole_location, merged_blackhole_velocity = self.set_binary_location_and_velocity(binary_orbital_position, binary_orbital_velocity)
        # Set the initial position and velocity of the merged_blackholes to be the same as was the last values from --- set_binary_location_and_velocity ---
        self.merged_blackhole.mass = fraction_of_total_mass * self.total_mass
        self.merged_blackhole.radius =  (2*constants.G*self.merged_blackhole.mass)/(constants.c**2) # Could as well be zero cause we dont care anymore about it
                                                                                                    # I only left it in because we want to avoid the new particle
                                                                                                    # colliding with other particles
        self.merged_blackhole.position = merged_blackhole_location
        self.merged_blackhole.velocity = merged_blackhole_velocity

        self.blackholes.remove_particles(self.blackholes[0], self.blackholes[1])
        self.blackholes.add_particle(self.merged_blackhole)




        # Need to add these particles to gravity instead of the binaries
        raise NotImplementedError


    def merge_particles(self, merge_condition):
        if merge_condition:
            merged_blackhole_attributes(total_mass)
