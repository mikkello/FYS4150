from numpy import *
import matplotlib.pyplot as plt

#########################
## Declaring variables ##
N = (100000-10000) # Monte Carlo cycles (after steady state)
temperatures = ["1.00", "1.00", "2.40", "2.40"]
states = ["ordered", "random", "ordered", "random"]
N_temp = len(temperatures)
L = 20

E = np.zeros((N_temp,N))

#########################################################
## Computing expectation values and standard deviation ##
## Based on datafiles from calculation                 ##
E_var_random_1_00 = [(-1.997*L*L), sqrt(0.0231573*L*L)]
E_var_random_2_40 = [(-1.23497*L*L), sqrt(1.41253*2.4*2.4*L*L)]
E_var_ordered_1_00 = [(-1.99717*L*L), sqrt(0.0467845*L*L)]
E_var_ordered_2_40 = [(-1.23639*L*L), sqrt(1.40413*2.4*2.4*L*L)]



##########################
## Importing variables  ##
for i in range(N_temp):
    filename = "E_M_"+states[i]+"_T"+temperatures[i]+".dat"
    infile = open(filename, "r")
    for j in range(N):
        infile.readline()
        infile.readline()
        infile.readline()
        E_M_line = infile.readline()
        variables = E_M_line.split()
        E[i][j] = variables[0]
    infile.close()
    
#######################################    
## Plotting probability distribution ##

plt.figure(num=None, figsize=(16, 10), dpi=120, facecolor='w', edgecolor='k')
plt.rcParams.update({'font.size': 16})


plt.subplot(2, 2, 1)
plt.title('a) Ordered. T = 1.0')
weights = np.ones_like(E[0][:])/float(len(E[0][:]))
plt.hist(E[0][:], weights=weights)
plt.ylabel(r'$P(E)$')
plt.legend()
plt.grid()

plt.subplot(2, 2, 2)
plt.title('a) T = 1.0')
weights = np.ones_like(E[1][:])/float(len(E[1][:]))
plt.hist(E[1][:], weights=weights)
plt.ylabel(r'$P(E)$')
plt.legend()
plt.grid()

plt.subplot(2, 2, 3)
weights = np.ones_like(E[2][:])/float(len(E[2][:]))
plt.hist(E[2][:], bins=500, weights=weights)
plt.axvline(E_var_ordered_2_40[0], color='r', label='Expectation value')
plt.axvline((E_var_ordered_2_40[0]-E_var_ordered_2_40[1]), color='g', label='Standard deviation')
plt.axvline((E_var_ordered_2_40[0]+E_var_ordered_2_40[1]), color='g')
plt.title('c) Ordered. T = 2.4')
plt.ylabel(r'$P(E)$')
plt.xlabel(r'$E$')
plt.legend(prop={'size': 11})
plt.grid()

plt.subplot(2, 2, 4)
weights = np.ones_like(E[3][:])/float(len(E[3][:]))
plt.hist(E[3][:], bins=500, weights=weights)
plt.axvline(E_var_random_2_40[0], color='r', label='Expectation value')
plt.axvline((E_var_random_2_40[0]-E_var_random_2_40[1]), color='g', label='Standard deviation')
plt.axvline((E_var_random_2_40[0]+E_var_random_2_40[1]), color='g')
plt.title('b) T = 2.4')
plt.ylabel(r'$P(E)$')
plt.xlabel(r'$E$')
plt.legend(prop={'size': 11})
plt.grid()

plt.savefig('taskd_P(E)_vs_E_1e5_L20_2')