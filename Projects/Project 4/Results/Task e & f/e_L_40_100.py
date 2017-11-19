import numpy as np
import matplotlib.pyplot as plt
from scipy.interpolate import CubicSpline

#########################
## Declaring variables ##


lattice = [40,60,80,100]#[40,60,80,100] # Runs of program with different cycles
N_files = len(lattice)
temps = [2.15, 2.2, 2.25, 2.3, 2.35]
colors = ['-ro','-bo','-go','-yo']
N_temps = len(temps)

E_tot = np.zeros((N_files, N_temps))
M_tot = np.zeros((N_files, N_temps))
C_v = np.zeros((N_files, N_temps))
chi = np.zeros((N_files, N_temps))
Onsager = [2.269]


##########################
## Importing variables  ##
for i in range(N_files):
    lattices = str(lattice[i])
    filename = "L"+lattices+"_100000.dat"
    infile = open(filename, "r")
    infile.readline()
    for j in range(N_temps):
        line = infile.readline()
        variables = line.split()
        E_tot[i][j] = variables[0]
        M_tot[i][j] = variables[1]
        C_v[i][j] = variables[2]
        chi[i][j] = variables[3]
infile.close()


########################################
## Calculating critial temperature at ##
lattice2 = [40,60,80,100]
T_crit = [2.25, 2.3, 2.3, 2.25] # C_v max
T_crit2 = [2.35, 2.3, 2.3, 2.3] # Chi max
"""
a = []
for i in range(0, len(T_crit)-1):
	for j in range(i+1, len(T_crit)):
		a.append( float(T_crit[i] - T_crit[j])/(1./lattice2[i] - 1./lattice2[j]))

T_infty = -sum(a)/len(a) *1./lattice2[1] + T[1];
print(T_infty)
"""
a = 0 
count = 0 
for i in range(len(T_crit)): 
    for j in range(len(T_crit)): 
        if i != j: 
            a +=  ((lattice2[i]*T_crit[i]) - (lattice2[j]*T_crit[j])) / ( (lattice2[i]) - (lattice2[j]) )
            count += 1
print(a/count)


"""
####################################################
## Plotting properties for multiple lattice sizes ##


plt.figure(num=None, figsize=(16, 10), dpi=120, facecolor='w', edgecolor='k')
plt.rcParams.update({'font.size': 16})
    
plt.subplot(2, 2, 1)
plt.title('a) Mean energy')
for i in range(N_files):
    plt.plot(temps[:],E_tot[i][:], colors[i],  label="L = %s" % lattice[i])
plt.ylabel(r'$\langle E\rangle$')
plt.legend()
plt.grid()
 

plt.subplot(2, 2, 2)
plt.title('b) Mean magnetization')
for i in range(N_files):
    plt.plot(temps[:],M_tot[i][:], colors[i],  label="L = %s" % lattice[i])
plt.ylabel(r"$\langle |M|\rangle$")
plt.legend()
plt.grid()


plt.subplot(2, 2, 3)
plt.title('c) Heat capacity')
for i in range(N_files):
    plt.plot(temps[:],C_v[i][:], colors[i],  label="L = %s" % lattice[i])
plt.axvline(Onsager[:], color='m', label='Onsager exact')
plt.ylabel(r'$C_v$')
plt.xlabel('T')
plt.legend()
plt.grid()


plt.subplot(2, 2, 4)
plt.title('d) Susceptibility')
for i in range(N_files):
    plt.plot(temps[:],chi[i][:], colors[i],  label="L = %s" % lattice[i])
plt.axvline(Onsager[:], color='m', label='Onsager exact')
plt.ylabel('$\chi$')
plt.xlabel('T')
plt.legend()
plt.grid()

plt.savefig("e_L_40_100.png")
"""