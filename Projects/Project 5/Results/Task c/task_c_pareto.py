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

#############################################################
## Importing variables and writing pareto exponent to file ##
f = open('task_c_tails.txt', 'w')
f.write('N_agents   Lambda   Alpha   Gamma   nu   variance' + '\n')
for i in range(len(Lambda)):
    file = filename(i)
    infile = open(file, "r")
    for j in range(N_bins):
        line = infile.readline()
        variables = line.split()
        m[i][j] = float(variables[0])
        counts[i][j] = float(variables[1])/int(N_agents)
        
    countss = np.trim_zeros(counts[i], 'b')
    mm = m[i][:len(countss)]
    
    count = 0
    n = len(countss) - 1 
    
    while count < 0.10 :
        count += countss[n]
        n -= 1
    countss = countss
   
    pareto = lambda x, nu: (x**(-1 -nu))
    nuu, var = scipy.optimize.curve_fit(pareto, mm[-n:], countss[-n:])
    
    f.write(r'%s  %s  %s  %s  %.5f  %.5f' % (N_agents, Lambda[i], Alpha, Gamma, nuu[0], var[0])+'\n')  
    
    infile.close()
f.close()
