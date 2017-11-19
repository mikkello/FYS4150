from numpy import *
import matplotlib.pyplot as plt

#########################
## Declaring variables ##
states = ["random", "ordered"]
temps = [1, 1.2, 1.4, 1.6, 1.8, 2.0, 2.2, 2.4]
N_states = len(states)
N_cycles = 100000 # Monte Carlo cycles
L = 20 # Lattice size LxL
N_temps = len(temps)
acc = zeros((N_states,N_temps)) # Accepted configurations

####################################
## Importing variables from files ##
for i in range(N_states):
    filename = "L20_100000_"+states[i]+".dat"
    infile = open(filename, "r")
    infile.readline()
    for j in range(N_temps):
        line = infile.readline()
        variables = line.split()
        acc[i][j]= float(variables[5])/(N_cycles*L*L) # Normalizing
infile.close()


###################################################################
## Plotting accepeted configurations (normalized) vs temperature ##

plt.figure(num=None, figsize=(16, 10), dpi=120, facecolor='w', edgecolor='k')
plt.rcParams.update({'font.size': 16})
    
plt.plot(temps[:],acc[0][:], 'b',  label="Random")
plt.title('Accepted configurations (normalized) as a function of temperature.')
plt.ylabel("Accepeted configurations (normalized)")
plt.xlabel("T")
plt.legend()
plt.grid()
plt.savefig("taskc3_acc_vs_T_L20_1e5")

