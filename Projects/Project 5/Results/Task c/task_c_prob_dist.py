import numpy as np
import matplotlib.pyplot as plt
import scipy.special
import scipy.optimize

################
## Parameters ##
N_agents = 500
N_MCC = 1000
N_transactions = 10000000
bin_width = 0.05
N_bins = int(N_agents/bin_width)
Lambda = ["0.000", "0.250", "0.500", "0.900"]
Alpha = "0.000"
Gamma = "0.000"

# Arrays for wealth distribution/histogram and theoretical Gibbs dist.
m = np.zeros((len(Lambda),N_bins))
counts = np.zeros((len(Lambda),N_bins))
gibbs = np.zeros((len(Lambda),N_bins))
n = np.zeros(len(Lambda))

def filename(llambda):
    filename = 'MONEY_agents_%s_MCC_%s_N_trans_%s_lambda_%s_alpha_%s_gamma_%s.dat' % (N_agents, N_MCC, N_transactions, Lambda[llambda], Alpha, Gamma)
    return filename

#########################
## Importing variables ##
for i in range(len(Lambda)):
    file = filename(i)
    infile = open(file, "r")
    for j in range(N_bins):
        line = infile.readline()
        variables = line.split()
        m[i][j] = float(variables[0])
        counts[i][j] = float(variables[1])/int(N_agents)
    infile.close()



##################################
## Plotting wealth distribution ##
plt.figure(num=0, figsize=(12, 10), dpi=120, facecolor='w', edgecolor='k')
plt.rcParams.update({'font.size': 16})
for i in range(len(Lambda)):
    plt.plot(m[i], counts[i], label=r'$\lambda$ = %s' % (Lambda[i]))

plt.xlim((0,2.4))
plt.grid()
plt.title("Probability distributions for different savings parameters")
plt.ylabel('Probability distribution')
plt.xlabel('Wealth')
plt.legend()
plt.savefig('task_c_prob_dist_agents_%s_MCC_%s_N_trans_%s_tails' % (N_agents, N_MCC, N_transactions))
plt.show()
