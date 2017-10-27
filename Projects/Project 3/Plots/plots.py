import numpy as np
import matplotlib.pyplot as plt

#########################
## Importing variables ##
infile = open("Variables_2_dt_0.000100.txt", "r")
infile.readline()
line = infile.readline()
variables = line.split()
N_steps = int(variables[0]) # Number of steps
N_planets = int(variables[1]) # Number of planets
N_years = int(variables[2]) # Number of years
dt = float(variables[3]) # Time step
t = np.linspace(0,N_years,N_steps) # Time vector for plotting
infile.close()




"""
######################################
## Importing data from Euler solver ##
infile = open("positions1.xyz", "r")

x = np.zeros((N_planets, N_steps))
y = np.zeros((N_planets, N_steps))
z = np.zeros((N_planets, N_steps))
E_tot = np.zeros((N_planets, N_steps))
E_pot = np.zeros((N_planets, N_steps))
E_kin = np.zeros((N_planets, N_steps))
Ang_m = np.zeros((N_planets, N_steps))

for dim in range(N_steps):
    for i in range(N_planets):
        line = infile.readline()
        words = line.split()
        x[i][dim] = words[0]
        y[i][dim] = words[1]
        z[i][dim] = words[2]
        E_tot[i][dim] = words[3]
        E_pot[i][dim] = words[4]
        E_kin[i][dim] = words[5]
        Ang_m[i][dim] = words[6]
        
infile.close()
"""
"""
###################################
## Estimation of escape velocity ##
k = [1.33, 1.35, 1.37, 1.39, 1.41, 1.43, 1.45, 1.47]
x = np.zeros((int(len(k)), N_planets, N_steps))
y = np.zeros((int(len(k)), N_planets, N_steps))

for a in range(len(k)):
    filename = "positions3_k_" + str(k[a]) + ".xyz"
    infile = open(filename, "r")
    for dim in range(N_steps):
        for i in range(N_planets):
            line = infile.readline()
            words = line.split()
            x[a][i][dim] = words[0]
            y[a][i][dim] = words[1]
infile.close()
plt.figure(num=None, figsize=(12, 10), dpi=120, facecolor='w', edgecolor='k')
plt.rcParams.update({'font.size': 16})
plt.grid()
for i in range(len(k)):
    plt.plot(x[i][1],y[i][1], label='k = %s' % (k[i]))
plt.legend()
plt.title(r'Estimation of the escape velocity of Earth. $V_{0} \mathrm{=} k\cdot \rm 2 \cdot \pi$ AU/yr')
plt.xlabel('x [AU]')
plt.ylabel('y [AU]')
plt.savefig('escape_velocity_dt_%s.png'%(dt))
"""


################################################
## Importing data from Velocity Verlet solver ##

infile = open("Positions_2_dt_0.000100_B4.xyz", "r")
x2 = np.zeros((N_planets, N_steps))
y2 = np.zeros((N_planets, N_steps))
z2 = np.zeros((N_planets, N_steps))
E_tot2 = np.zeros((N_planets, N_steps))
E_pot2 = np.zeros((N_planets, N_steps))
E_kin2= np.zeros((N_planets, N_steps))
Ang_m2 = np.zeros((N_planets, N_steps))


for dim in range(N_steps):
    for i in range(N_planets):
        line = infile.readline()
        words = line.split()
        x2[i][dim] = words[0]
        y2[i][dim] = words[1]
        z2[i][dim] = words[2]
       
 

#######################################################
## Dividing total energy values by the initial value ##
"""
E_in = E_tot[1][0]
E_tot_rel = np.divide(E_tot[1],E_in)
E_in = E_pot[1][0]
E_pot_rel = np.divide(E_pot[1],E_in)
E_in = E_kin[1][0]
E_kin_rel = np.divide(E_kin[1],E_in)

E_in = E_tot2[1][0]
E_tot2_rel = np.divide(E_tot2[1],E_in)
E_in = E_pot2[1][0]
E_pot2_rel = np.divide(E_pot2[1],E_in)
E_in = E_kin2[1][0]
E_kin2_rel = np.divide(E_kin2[1],E_in)
"""
##################################################################################
## Plotting relative total energies and ang. mom. for Euler and Velocity Verlet ##
"""
plt.figure(num=None, figsize=(12, 10), dpi=120, facecolor='w', edgecolor='k')
plt.rcParams.update({'font.size': 16})
plt.plot(t[:],E_tot_rel, 'b',  label='Total energy (FE)')
plt.plot(t[:],E_tot2_rel, 'g', label='Total energy (VV)')
plt.xlabel('Years')
plt.ylabel(r'$\frac{Total / energy}{Initial / energy}$')
plt.legend()
plt.title('Conservation of total energy for the Euler and Velocity Verlet methods.\n dt = %s' %(dt) )
#plt.savefig('cons_E_rel_dt_%s.png'%(dt))


plt.figure(num=None, figsize=(14, 10), dpi=120, facecolor='w', edgecolor='k')
plt.rcParams.update({'font.size': 16})
plt.subplot(2, 1, 1)
plt.title('Conservation of energies for the Euler and Velocity Verlet methods.\n dt = %s' %(dt) )
plt.plot(t[:],E_tot[1], 'b',  label='Total energy (FE)')
plt.plot(t[:],E_kin[1], 'g',  label='Kinetic energy (FE)')
plt.plot(t[:],E_pot[1], 'r',  label='Potential energy (FE)')
plt.ylabel('Energy [$M_{sun}AU^{2}/yr^{2}$]')
plt.legend()

plt.subplot(2, 1, 2)
plt.plot(t[:],E_tot2[1], 'b', label='Total energy (VV)')
plt.plot(t[:],E_kin2[1], 'g', label='Kinetic energy (VV)')
plt.plot(t[:],E_pot2[1], 'r', label='Potential energy (VV)')
plt.ylabel('Energy [$M_{sun}AU^{2}/yr^{2}$]')
plt.xlabel('Years')
plt.legend()

#plt.savefig('cons_E_dt_%s.png'%(dt))

plt.figure(num=None, figsize=(14, 10), dpi=120, facecolor='w', edgecolor='k')
plt.rcParams.update({'font.size': 16})
plt.plot(t[:],Ang_m[1], 'b',  label='Angular momentum (FE)')
plt.plot(t[:],Ang_m2[1], 'g', label='Angular momentum (VV)')
plt.xlabel('Years')
plt.ylabel(r'Angular momentum$[M_{Sun}AU^{2}/yr]$')
plt.legend()
plt.title('Conservation of angular momentum for the Euler and Velocity Verlet methods.\n dt = %s' %(dt) )
#plt.savefig('cons_Ang_dt_%s.png'%(dt))
"""

#####################################################################################################
## Plotting the calculated orbits for the Sun-Earth system using Euler and Velocity Verlet methods ##

if N_planets == 2:
    plt.figure(num=None, figsize=(12, 10), dpi=120, facecolor='w', edgecolor='k')
    plt.rcParams.update({'font.size': 16})
    #plt.plot(x[1], y[1], 'b', label='Earth (FE)')
    plt.plot(x2[0], y2[0],'ro', label='Sun')
    plt.plot(x2[1], y2[1], 'g', label='Earth (VV)')
    plt.xlabel('x [AU]')
    plt.ylabel('y [AU]')
    plt.legend()
    plt.title('Earth-Sun system using Forward Euler and Velocity Verlet (dt = %s)' %(dt) )
    plt.savefig('test3.png')
    #plt.savefig('earth_sun_dt_%s.png'%(dt))
    plt.show()

"""
if N_planets == 3:
    plt.figure(num=None, figsize=(12, 10), dpi=120, facecolor='w', edgecolor='k')
    plt.rcParams.update({'font.size': 16})
    plt.plot(x2[1], y2[1], 'g', label='Earth')
    plt.plot(x2[0], y2[0],'ro', label='Sun')
    plt.plot(x2[2], y2[2], 'b', label='Jupiter')
    plt.xlabel('x [AU]')
    plt.ylabel('y [AU]')
    plt.legend()
    plt.title('Earth-Jupiter-Sun system (dt = %s)' %(dt) )
    #plt.savefig('test3.png')
    plt.savefig('earth_jupiter_notfixed_dt_%s.png'%(dt))
    plt.show()"""