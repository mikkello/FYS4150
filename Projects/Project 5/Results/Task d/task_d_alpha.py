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


#########################
## Importing variables ##

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
        infile.close()



##################################
## Plotting wealth distribution ##

# N_agents = 500, Lambda = 0.0
plt.figure(num='none', figsize=(12, 10), dpi=120, facecolor='w', edgecolor='k')
plt.rcParams.update({'font.size': 16})
plt.subplot(2, 2, 1)
for i in range(len(Alpha)):
    plt.loglog(m500[i][0][:], counts500[i][0][:], label=r'$\alpha$ = %.2f' % float(Alpha[i]))
plt.grid()
plt.legend()
plt.title(r'a) $N_{agents}$ = 500, $\lambda$ = 0.0')
plt.ylabel(r'Probability distribution')

# N_agents = 500, Lambda = 0.5
plt.subplot(2, 2, 2)
for i in range(len(Alpha)):
    plt.loglog(m500[i][1][:], counts500[i][1][:], label=r'$\alpha$ = %.2f' % float(Alpha[i]))
plt.grid()
plt.legend()
plt.title(r'b) $N_{agents}$ = 500, $\lambda$ = 0.5')
#plt.ylabel(r'Probability distribution')

# N_agents = 1000, Lambda = 0.0
plt.subplot(2, 2, 3)
for i in range(len(Alpha)):
    plt.loglog(m1000[i][0][:], counts1000[i][0][:], label=r'$\alpha$ = %.2f' % float(Alpha[i]))
plt.grid()
plt.legend()
plt.ylabel(r'Probability distribution')
plt.title(r'c) $N_{agents}$ = 1000, $\lambda$ = 0.0')
plt.xlabel(r'Wealth')


# N_agents = 1000, Lambda = 0.5
plt.subplot(2, 2, 4)
for i in range(len(Alpha)):
    plt.loglog(m1000[i][1][:], counts1000[i][1][:], label=r'$\alpha$ = %.2f' % float(Alpha[i]))
plt.grid()
plt.legend()
#plt.ylabel(r'Probability distribution')
plt.title(r'a) $N_{agents}$ = 1000, $\lambda$ = 0.5')
plt.xlabel(r'Wealth')


plt.suptitle(r'Distribution of wealth with varying nearest neighbour parameter $\alpha$')
plt.savefig('task_d_alpha')
plt.show()
