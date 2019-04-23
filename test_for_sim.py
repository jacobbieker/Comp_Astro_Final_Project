from amuse.ext.orbital_elements import orbital_elements_from_binary

def get_orbital_elements_of_triple(self, smbh, stars):
    inner_binary = stars
    outer_binary = Particles(1)
    outer_binary[0].mass = inner_binary.mass.sum()
    outer_binary[0].position = inner_binary.center_of_mass()
    outer_binary[0].velocity = inner_binary.center_of_mass_velocity()
    outer_binary.add_particle(smbh)
    M1, M2, ain, ein, ta_in, inc_in, lan_in, aop_in \
        = orbital_elements_from_binary(inner_binary, G=constants.G)
    M12, M3, aout, eout, ta_out, inc_out, lan_out, aop_out \
        = orbital_elements_from_binary(outer_binary, G=constants.G)
    return ain, ein, aout, eout, inc_in, inc_out
