import numpy as np
import matplotlib.pyplot as plt

################
## Parameters ##
N_agents = ["1000"]
N_MCC = 1000
N_transactions = 10000000
bin_width = 0.05
N_bins = [int(int(N_agents[0])/bin_width)]
Lambda = ["0.000", "0.500"]
Alpha = ["1.00", "2.00"]
Gamma = ["0.000", "1.00", "2.00", "3.00", "4.00"]

# Arrays for Lambda = 0 and Lambda = 0.5
m0 = np.zeros((len(Alpha),len(Gamma),N_bins[0]))
m05 = np.zeros((len(Alpha),len(Gamma),N_bins[0]))
counts0 = np.zeros((len(Alpha),len(Gamma),N_bins[0]))
counts05 = np.zeros((len(Alpha),len(Gamma),N_bins[0]))

def filename(llambda, aalpha, ggamma):
    filename = 'MONEY_agents_%s_MCC_%s_N_trans_%s_lambda_%s_alpha_%s_gamma_%s.dat' % (N_agents[0], N_MCC, N_transactions, Lambda[llambda], Alpha[aalpha], Gamma[ggamma])
    return filename


#############################################################
## Importing variables and writing pareto exponent to file ##
f = open('task_e_tails.txt', 'w')
f.write('N_agents   Lambda   Alpha   Gamma   nu   variance' + '\n')

# Lambda = 0, Alpha = [1.0,2.0], Gamma = [0.0-4.0]
for i in range(len(Alpha)):
    for j in range(len(Gamma)):
        file = filename(0,i,j)
        infile = open(file, "r")
        for k in range(N_bins[0]):
            line = infile.readline()
            variables = line.split()
            m0[i][j][k] = float(variables[0])
            counts0[i][j][k] = float(variables[1])/int(N_agents[0])
        
        countss = np.trim_zeros(counts0[i][j], 'b')
        mm = m0[i][j][:len(countss)]
        
        count = 0
        n = len(countss) - 1 
        while count < 0.10 :
            count += countss[n]
            n -= 1
        countss = countss
   
        pareto = lambda x, nu: (x**(-1 -nu))
        nuu, var = scipy.optimize.curve_fit(pareto, mm[-n:], countss[-n:])
        f.write(r'%s  %s  %s  %s  %.5f  %.5f' % (N_agents[0], Lambda[0], Alpha[i], Gamma[j], nuu[0], var[0])+'\n') 
        
        infile.close()
        

# Lambda = 0.5, Alpha = [1.0,2.0], Gamma = [0.0-4.0]
for i in range(len(Alpha)):
    for j in range(len(Gamma)):
        file = filename(1,i,j)
        infile = open(file, "r")
        for k in range(N_bins[0]):
            line = infile.readline()
            variables = line.split()
            m05[i][j][k] = float(variables[0])
            counts05[i][j][k] = float(variables[1])/int(N_agents[0])
        
        countss = np.trim_zeros(counts05[i][j], 'b')
        mm = m05[i][j][:len(countss)]
        
        count = 0
        n = len(countss) - 1 
        while count < 0.10 :
            count += countss[n]
            n -= 1
        countss = countss
   
        pareto = lambda x, nu: (x**(-1 -nu))
        nuu, var = scipy.optimize.curve_fit(pareto, mm[-n:], countss[-n:])
        f.write(r'%s  %s  %s  %s  %.5f  %.5f' % (N_agents[0], Lambda[1], Alpha[i], Gamma[j], nuu[0], var[0])+'\n')
        infile.close()

f.close()