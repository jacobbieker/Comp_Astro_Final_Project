# Comp_Astro_Final_Project
## *"Dynamical Hardening of Black Hole Binaries in Active Galactic Nuclei"*

__Contributors__:   Giannis Politopoulos, 
		Jacob Bieker, 
		Zorry Belcheva, 
		Dylan Natoewal, 
		Benjamin Gilliam, 
		Lindsey Oberhelman

### Summary
The aim of this project, is to simulate an AGN around which, there is a massive disk and a number of binary stellar-mass blackholes. The AGN consists of a SMBH ~1e6 MSun and the binary blackholes ~30 MSun whilst the disk ~10% of the SMBH mass.
The idea behind the simulation, is that the binary blackholes will interact with the AGN disk while orbiting the SMBH, lose energy hence reducing their orbital period and semi major axis. This effect becomes greater the closer the binary blackholes get resulting to their merging. Such mergers are potential candidates for the gravitational waves detected by LIGO.

The total code consists of 7 files.

__SuperMAssiveBlackHole.py__ : Creates the SMBH at the center of the grid

__AccretionDisk.py__ : Creates a disk, using the ProtoPlanetaryDisk in AMUSE

__BinaryBlackHole.py__ : Creates a binary blackhole pair and sets it in orbit around the center of mass. It then sets the center of mass of the binary in orbit around the SMBH.

__BinaryBlackHolesWithAGN.py__ : Calls the aforementioned classes to create the system to be simulated.

__Gadget2_Gravity.py__ : Is an extended version of Gadget2 AMUSE package, to include the function: get_gravity_at_point for the gravitational interaction of the disk with the binaries

__main.py__ : Wraps up everything to be ran for the simulation

__plotting.py__ : Uses everything saved in the simulation for plotting and animation

## To run the simulation
The simulation is ran through AMUSE by running the __main.py__ file. The simulation parameters can be seen and tuned through the Option Parser implemented.

__The main parameters are:__ 

  **mass_of_central_black_hole**,	_default=1000000.0 MSun_
  
  **number_of_binaries**,	_default=50_
  
  **number_of_gas_particles**,	_default=100000_
  
**end_time**,	_default=10 Myr_
  
  **blackhole_mass**,	_default=30 MSun_
  
  **gravity_timestep**,	_default=100 yr_
  
  **bridge_timestep**,	_default=0.1 Myr_
  
  **smbh_as_potential**,	_default=False_
  
  **binaries_affect_disk**,	_default=False_
  
  **disk_mass_fraction**,	_default=0.1_
  
  **number_of_hydro_workers**,	_default=6_
  
  **number_of_grav_workers**,	_default=12_
  
**filename**,	_default=BinaryBlackHoles_

