import numpy as np
import matplotlib.pyplot as plt
import scipy.special

################
## Parameters ##
N_agents = 500
N_MCC = 1000
N_transactions = 10000000
bin_width = 0.05
N_bins = int(N_agents/bin_width)
Lambda = "0.000"
Alpha = "0.000"
Gamma = "0.000"
filename = 'MONEY_agents_%s_MCC_%s_N_trans_%s_lambda_%s_alpha_%s_gamma_%s.dat' % (N_agents, N_MCC, N_transactions, Lambda, Alpha, Gamma)

# Arrays for wealth dist. histogram
m = np.zeros(N_bins)
counts = np.zeros(N_bins)

#########################
## Importing variables ##
infile = open(filename, "r")
for i in range(N_bins):
    line = infile.readline()
    variables = line.split()
    m[i] = float(variables[0])
    counts[i] = float(variables[1])
infile.close()

counts = np.trim_zeros(counts, 'b')
m = m[:len(counts)]

gibbs =  ( 1**1/scipy.special.gamma(1) ) * ( np.exp(-m) )

##################################
## Plotting wealth distribution ##

# Histogram
plt.figure(num='none', figsize=(12, 10), dpi=120, facecolor='w', edgecolor='k')
plt.rcParams.update({'font.size': 16})
plt.subplot(2, 1, 1)
plt.plot(m, gibbs, label='Gibbs distribution', color='b')
plt.hist(m, weights=counts, bins=m, normed=True, label='Simulation', color='r')
plt.grid()
plt.title(r'a)')
plt.ylabel('Agent distribution (normalized)')
plt.legend()

# Log plot
plt.rcParams.update({'font.size': 16})
plt.subplot(2, 1, 2)
plt.plot(m, counts/int(N_agents),'ro',label='Simulation')
plt.plot(m, 0.05*gibbs, label='Gibbs distribution', color='b')
plt.grid()
plt.yscale('log')
plt.ylabel('Probability distribution')
plt.xlabel('Wealth')
plt.title(r'b)')
plt.legend()

plt.suptitle("Wealth distribution from simple model")
plt.savefig('task_b_histogram_agents_%s_MCC_%s_N_trans_%s' % (N_agents, N_MCC, N_transactions))
