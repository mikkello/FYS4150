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


#########################
## Importing variables ##

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
        infile.close()


##################################
## Plotting wealth distribution ##

# Alpha = 1, Lambda = 0
plt.figure(num='none', figsize=(12, 10), dpi=120, facecolor='w', edgecolor='k')
plt.rcParams.update({'font.size': 16})
plt.subplot(2, 2, 1)
for i in range(len(Gamma)):
    plt.loglog(m0[0][i][:], counts0[0][i][:], label=r'$\gamma$ = %.2f' % float(Gamma[i]))
plt.grid()
plt.legend()
plt.title(r'a) $\alpha$ = 0, $\lambda$ = 0.0')
plt.ylabel('Probability distribution')



# Alpha = 1, Lambda = 0.5
plt.subplot(2, 2, 2)
for i in range(len(Gamma)):
    plt.loglog(m05[0][i][:], counts05[0][i][:], label=r'$\gamma$ = %.2f' % float(Gamma[i]))
plt.grid()
plt.legend()
plt.title(r'b) $\alpha$ = 1, $\lambda$ = 0.5')
#plt.ylabel('Probability distribution')



# Alpha = 2, Lambda = 0
plt.subplot(2, 2, 3)
for i in range(len(Gamma)):
    plt.loglog(m0[1][i][:], counts0[1][i][:], label=r'$\gamma$ = %.2f' % float(Gamma[i]))
plt.grid()
plt.legend()
plt.title(r'c) $\alpha$ = 2, $\lambda$ = 0.0')
plt.xlabel('Wealth')
plt.ylabel('Probability distribution')




# Alpha = 2, Lambda = 0.5
plt.subplot(2, 2, 4)
for i in range(len(Gamma)):
    plt.loglog(m05[1][i][:], counts05[1][i][:], label=r'$\gamma$ = %.2f' % float(Gamma[i]))
plt.grid()
plt.legend()
plt.title(r'a) $\alpha$ = 2, $\lambda$ = 0.5')
plt.xlabel('Wealth')
#plt.ylabel('Probability distribution')

plt.suptitle(r'Distribution of wealth with varying memory parameter $\gamma$')
plt.savefig('task_e_gamma')
plt.show()
