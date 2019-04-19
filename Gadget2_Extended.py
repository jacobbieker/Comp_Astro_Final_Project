from amuse.community.gadget2.interface import Gadget2
import numpy as np
from amuse.units import constants, units
from amuse.couple.bridge import CalculateFieldForParticles

class Gadget2_Extended(Gadget2):
    def __init__(self, radius=None, unit_converter=None, mode='normal', **options):
        Gadget2.__init__(self, unit_converter=unit_converter, mode=mode, **options)
        self.radius_for_gravity_calc = radius

    def get_gravity_at_point(self, alpha, x, y, z):
        print("Start Field Code", flush=True)
        field_code = CalculateFieldForParticles(particles=self.gas_particles, G=constants.G)
        print("Made Field Code", flush=True)
        #print(field_code.get_gravity_at_point(alpha, x, y, z))
        return field_code.get_gravity_at_point(alpha, x, y, z)
        '''
        if self.radius_for_gravity_calc is None:
            self.radius_for_gravity_calc = 2. | x.unit
        center_point = np.asarray([x.value_in(x.unit), y.value_in(x.unit),z.value_in(x.unit)]) | x.unit
        nearby_particles = self.gas_particles.select(
            lambda r: np.abs((r - center_point).length()) < self.radius_for_gravity_calc, ["position"])

        if len(nearby_particles) <= 0:
            # No particles close enough to affect it
            return np.asarray([[0.0], [0.0], [0.0]]) | units.m / (units.s*units.s)
        print(len(nearby_particles), flush=True)

        # Now get center of mass of the particles
        # Then calculate acceleration from that "particle" to the point

        center_of_mass = nearby_particles.center_of_mass()
        total_mass = nearby_particles.mass.sum()

        # GMm/r^2
        x_y_z = [center_of_mass[0].value_in(center_of_mass.unit)-center_point[0].value_in(center_of_mass.unit),
                 center_of_mass[1].value_in(center_of_mass.unit)-center_point[1].value_in(center_of_mass.unit),
                 center_of_mass[2].value_in(center_of_mass.unit)-center_point[2].value_in(center_of_mass.unit)] | center_of_mass.unit
        print(x_y_z, flush=True)
        ax,ay,az = (constants.G*total_mass)/((x_y_z)**2)
        print(ax, flush=True)
        return ax,ay,az
        '''

    def add_particles(self, particles):
        self.gas_particles.add_particles(particles)

