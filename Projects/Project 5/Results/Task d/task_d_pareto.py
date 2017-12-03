import numpy as np
import matplotlib.pyplot as plt

################
## Parameters ##
N_agents = ["500", "1000"]
N_MCC = 1000
N_transactions = 10000000
bin_width = 0.05
N_bins = [int(int(N_agents[0])/bin_width),int(int(N_agents[1])/bin_width)]
Lambda = ["0.000", "0.500"]
Alpha = ["0.000", "0.500", "1.00", "1.50", "2.00"]
Gamma = "0.000"

# Arrays for N_agents = 500 and N_agents = 1000
m500 = np.zeros((len(Alpha),len(Lambda),N_bins[0]))
m1000 = np.zeros((len(Alpha),len(Lambda),N_bins[1]))
counts500 = np.zeros((len(Alpha),len(Lambda),N_bins[0]))
counts1000 = np.zeros((len(Alpha),len(Lambda),N_bins[1]))

def filename(aagents, llambda, aalpha):
    filename = 'MONEY_agents_%s_MCC_%s_N_trans_%s_lambda_%s_alpha_%s_gamma_%s.dat' % (N_agents[aagents], N_MCC, N_transactions, Lambda[llambda], Alpha[aalpha], Gamma)
    return filename


#############################################################
## Importing variables and writing pareto exponent to file ##
f = open('task_d_tails.txt', 'w')
f.write('N_agents   Lambda   Alpha   Gamma   nu   variance' + '\n')

# N_agents = 500
for i in range(len(Alpha)):
    for j in range(len(Lambda)):
        file = filename(0,j,i)
        infile = open(file, "r")
        for k in range(N_bins[0]):
            line = infile.readline()
            variables = line.split()
            m500[i][j][k] = float(variables[0])
            counts500[i][j][k] = float(variables[1])/int(N_agents[0])
        
        countss = np.trim_zeros(counts500[i][j], 'b')
        mm = m500[i][j][:len(countss)]
        
        count = 0
        n = len(countss) - 1 
    
        while count < 0.10 :
            count += countss[n]
            n -= 1
        countss = countss
   
        pareto = lambda x, nu: (x**(-1 -nu))
        nuu, var = scipy.optimize.curve_fit(pareto, mm[-n:], countss[-n:])
    
        f.write(r'%s  %s  %s  %s  %.5f  %.5f' % (N_agents[0], Lambda[j], Alpha[i], Gamma, nuu[0], var[0])+'\n') 
        infile.close()
 

# N_agents = 1000
for i in range(len(Alpha)):
    for j in range(len(Lambda)):
        file = filename(1,j,i)
        infile = open(file, "r")
        for k in range(N_bins[1]):
            line = infile.readline()
            variables = line.split()
            m1000[i][j][k] = float(variables[0])
            counts1000[i][j][k] = float(variables[1])/int(N_agents[1])
        countss = np.trim_zeros(counts1000[i][j], 'b')
        mm = m1000[i][j][:len(countss)]
        
        count = 0
        n = len(countss) - 1 
    
        while count < 0.10 :
            count += countss[n]
            n -= 1
        countss = countss
   
        pareto = lambda x, nu: (x**(-1 -nu))
        nuu, var = scipy.optimize.curve_fit(pareto, mm[-n:], countss[-n:])
    
        f.write(r'%s  %s  %s  %s  %.5f  %.5f' % (N_agents[1], Lambda[j], Alpha[i], Gamma, nuu[0], var[0])+'\n')
        infile.close()
f.close