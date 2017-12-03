import numpy as np
import matplotlib.pyplot as plt

################
## Parameters ##
N_agents = 500
N_MCC = 2000-2
N_transactions = 10000000
bin_width = 0.05
N_bins = int(N_agents/bin_width)
Lambda = "0.000"
Alpha = "0.000"
Gamma = "0.000"
filename = 'CYCLE_agents_%s_MCC_%s_N_trans_%s_lambda_%s_alpha_%s_gamma_%s.dat' % (N_agents, N_MCC+2, N_transactions, Lambda, Alpha, Gamma)

# Arrays for Monte Carlo cycle numbers and variance
MCC = np.zeros(N_MCC)
sigma = np.zeros(N_MCC)


#########################
## Importing variables ##
infile = open(filename, "r")
infile.readline()
for i in range(N_MCC):
    line = infile.readline()
    variables = line.split()
    MCC[i] = float(variables[0])
    sigma[i] = float(variables[1])
infile.close()

##############################
## Plotting variance vs MCC ##
plt.figure(num='none', figsize=(12, 10), dpi=120, facecolor='w', edgecolor='k')
plt.rcParams.update({'font.size': 16})
plt.plot(MCC, sigma, color='b')
plt.grid()
plt.xlabel('Monte Carlo cycles')
plt.ylabel(r'$\sigma^{2}_{m}$')
plt.title('Equilibrium analysis')
plt.savefig('task_b_equi')
plt.tight_layout()
plt.show()