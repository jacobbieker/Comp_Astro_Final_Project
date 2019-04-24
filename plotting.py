import numpy
from amuse.support.data import particle_attributes
from amuse.io import read_set_from_file
from amuse.lab import units
#from amuse.lab import *
import h5py
from matplotlib.pyplot import plot, xlabel, ylabel, title, show
from mpl_toolkits.mplot3d import Axes3D
from matplotlib import pyplot as plt
from Gadget2_Gravity import Gadget2_Gravity

many_snapshots = read_set_from_file("[56962]_Particles_50_Binaries_1000000_Gas_AGN.h5", 'hdf5')


x = []
y = []
z = []

def grouped(iterable, n):
    return zip(*[iter(iterable)]*n)


def make_map(hydro, grid_points=5000, L=1):
    x, y = numpy.indices((grid_points + 1, grid_points + 1))

    x = L * (x.flatten() - grid_points / 2.) / grid_points
    y = L * (y.flatten() - grid_points / 2.) / grid_points
    z = x * 0.
    vx = 0. * x
    vy = 0. * x
    vz = 0. * x

    x = units.parsec(x)
    y = units.parsec(y)
    z = units.parsec(z)
    vx = units.kms(vx)
    vy = units.kms(vy)
    vz = units.kms(vz)

    rho, rhovx, rhovy, rhovz, rhoe = hydro.get_hydro_state_at_point(x, y, z, vx, vy, vz)
    rho = rho.reshape((grid_points + 1, grid_points + 1))

    return rho

def make_density_map():
    from mpl_toolkits.mplot3d import Axes3D
    from matplotlib import pyplot as plt
    i = 0
    for grid in many_snapshots.history:
        i += 1

        fig = plt.figure(figsize=(10, 10))
        ax = fig.add_subplot(111, projection='3d')
        graph = ax.scatter(grid.x.value_in(grid.x.unit), grid.y.value_in(grid.y.unit), grid.z.value_in(grid.z.unit))
        plt.savefig("Rho_Map_3d_{}.png".format(i), dpi=300)

        hydro = Gadget2_Gravity(mode='normal', number_of_workers=12)

        hydro.gas_particles.add_particles(grid)

        rho = make_map(hydro)
        plt.imshow(rho)
        plt.savefig("Density_Map_Center_{}.png".format(i), dpi=300)



def sep_vs_time():
    """
    Makes to plot of the separation between each of the binaries as 
    a function of time
    """

    sep = []
    average_separation_time = []
    all_separations = []
    list_variances = []
    num_time_steps = 0

    for si in many_snapshots.history:
        #print(si[:20].mass.in_(units.MSun))
        #x.append(si[:number_of_particles].x.value_in(units.AU))
        #y.append(si[:number_of_particles].y.value_in(units.AU))
        #z.append(si[:number_of_particles].x.value_in(units.AU))
        for i in grouped(si, 2):
            sep.append((i[0].position-i[1].position).length().value_in(units.AU))

        separations = numpy.asarray(sep)
        all_separations.append(separations) #Stores all the seperations in arrays per timestep

        mean_sep_per_time_step = numpy.mean(separations) #Get the mean separation per timestep
        average_separation_time.append(mean_sep_per_time_step)

        #Finding the Variance for each time step
        stdv = numpy.std(separations)
        variance = stdv**2
        list_variances.append(variance)

        num_time_steps += 1

        sep = []


    #separations = np.asarray(separations)
    time_steps = numpy.linspace(0,0.1*num_time_steps, len(average_separation_time))
    print (time_steps)
    average_separation_time = numpy.asarray(average_separation_time)
    print(average_separation_time)

    #fig1 = plt.figure()
    #ax = fig1.add_subplot(111, projection='3d')
    #ax.scatter(x,y,z)
    #plt.show()

    slopes = []
    for binary in all_separations:
        slopes.append((binary[1]-binary[-1])/(time_steps[1]-time_steps[-1]))



    plt.clf()
    plt.plot(time_steps[1:], all_separations[1:])
    #plt.errorbar(average_separation_time, time_steps, yerr=list_variances)
    plt.xlabel('Time [Myr]')
    plt.ylabel('Separation [AU]')
    plt.savefig("All_Big_SepTime.png", dpi=300)

sep_vs_time()
#make_density_map()

