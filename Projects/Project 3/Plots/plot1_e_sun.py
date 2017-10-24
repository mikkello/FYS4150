import numpy as np
import matplotlib.pyplot as plt

#########################
## Importing variables ##
infile = open("variables1.txt", "r")
infile.readline()
line = infile.readline()
variables = line.split()
N_steps = int(variables[0]) # Number of steps
N_planets = int(variables[1]) # Number of planets
N_years = int(variables[2]) # Number of years
dt = float(variables[3]) # Time step
t = np.linspace(0,N_years,N_steps) # Time vector for plotting
infile.close()

######################################
## Importing data from Euler solver ##
infile = open("positions1.xyz", "r")

x = np.zeros((N_planets, N_steps))
y = np.zeros((N_planets, N_steps))
z = np.zeros((N_planets, N_steps))
E_tot = np.zeros((N_planets, N_steps))
Ang_m = np.zeros((N_planets, N_steps))

for dim in range(N_steps):
    for i in range(N_planets):
        line = infile.readline()
        words = line.split()
        x[i][dim] = words[0]
        y[i][dim] = words[1]
        z[i][dim] = words[2]
        E_tot[i][dim] = words[3]
        #Ang_m[i][dim] = words[4]
        
infile.close()

################################################
## Importing data from Velocity Verlet solver ##

infile = open("positions2.xyz", "r")
x2 = np.zeros((N_planets, N_steps))
y2 = np.zeros((N_planets, N_steps))
z2 = np.zeros((N_planets, N_steps))
E_tot2 = np.zeros((N_planets, N_steps))
#Ang_m2 = np.zeros((N_planets, N_steps))


for dim in range(N_steps):
    for i in range(N_planets):
        line = infile.readline()
        words = line.split()
        x2[i][dim] = words[0]
        y2[i][dim] = words[1]
        z2[i][dim] = words[2]
        E_tot2[i][dim] = words[3]
        #Ang_m2[i][dim] = words[4]

#######################################################
## Dividing total energy values by the initial value ##
E_in = E_tot[1][0]
E_tot = np.divide(E_tot[1],E_in)
E_in2 = E_tot2[1][0]
E_tot2 = np.divide(E_tot2[1],E_in2)

####################################################################
## Plotting relative total energies for Euler and Velocity Verlet ##
plt.figure(num=None, figsize=(12, 10), dpi=120, facecolor='w', edgecolor='k')
plt.rcParams.update({'font.size': 20})
plt.plot(t[:],E_tot, label='Euler')
plt.plot(t[:],E_tot2, label='Velocity Verlet')
plt.xlabel('Years')
plt.ylabel('$frac{Total energy}{Initial energy}$')

#####################################################################################################
## Plotting the calculated orbits for the Sun-Earth system using Euler and Velocity Verlet methods ##
"""
if N_planets == 2:
    plt.figure(num=None, figsize=(12, 10), dpi=120, facecolor='w', edgecolor='k')
    plt.plot(x[1], y[1], 'b', label='Earth (FE)')
    plt.plot(x[0], y[0],'ro', label='Sun')
    plt.plot(x2[1], y2[1], 'g', label='Earth (VV)')
    plt.xlabel('x [AU]')
    plt.ylabel('y [AU]')
    plt.legend()
    plt.title('Earth-Sun system using Forward Euler and Velocity Verlet (dt = %s)' %(dt) )
    #plt.savefig('earth_sun_dt_%s.png'%(dt))
    plt.show()
    """
