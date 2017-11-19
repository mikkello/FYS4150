from numpy import *
import matplotlib.pyplot as plt

#########################
## Declaring variables ##
N = 100000 # Monte Carlo cycles
temperatures = ["1.00", "1.00", "2.40", "2.40"]
states = ["ordered", "random", "ordered", "random"]
N_temp = len(temperatures)

E = np.zeros((N_temp,N))
M = np.zeros((N_temp,N))
MCC = np.linspace(1,N,N)
accept = np.zeros((N_temp,N))


##########################
## Importing variables  ##
for i in range(N_temp):
    filename = "E_M_"+states[i]+"_T"+temperatures[i]+".dat"
    infile = open(filename, "r")
    for j in range(N):
        infile.readline()
        acceptN = infile.readline()
        norm = float(MCC[j]*20*20)
        accept[i][j] = float(acceptN)/norm
        E_M_line = infile.readline()
        variables = E_M_line.split()
        E[i][j] = variables[0]
        M[i][j] = variables[1]
        infile.readline()
        
    infile.close()

##################################################################
## Plotting accepted states (normalized) as a function of #MCC  ##
plt.figure(num=None, figsize=(16, 10), dpi=120, facecolor='w', edgecolor='k')
plt.rcParams.update({'font.size': 16})
plt.title('Accepted configurations (normalized) as a function of number of Monte Carlo cycles.')
plt.plot(MCC[:],accept[1][:], 'b',  label="T = 1.0")
plt.plot(MCC[:],accept[3][:], 'r',  label="T = 2.4")
plt.ylabel("Accepted configurations (normalized)")
plt.xlabel("Monte Carlo cycles")
plt.xscale('log')
plt.legend()
plt.grid()

plt.savefig("taskc_2")

"""
#######################################################################################################
## Plotting mean energy, mean magnetization, heat capacity and susceptibility as a function of #MCC  ##

plt.figure(num=None, figsize=(16, 10), dpi=120, facecolor='w', edgecolor='k')
plt.rcParams.update({'font.size': 16})
    
plt.subplot(2, 2, 1)
plt.title('a) Random')
plt.plot(MCC[:],E[1][:], 'b',  label="T = 1.0")
plt.plot(MCC[:],E[3][:], 'r',  label="T = 2.4")
plt.ylabel(r"$\langle E\rangle$")
plt.xscale('log')
plt.legend()
plt.grid()
 

plt.subplot(2, 2, 2)
plt.title('b) Ordered')
plt.plot(MCC[:],E[0][:], 'b',  label="T = 1.0")
plt.plot(MCC[:],E[2][:], 'r',  label="T = 2.4")
plt.ylabel(r"$\langle E\rangle$")
plt.xscale('log')
plt.legend()
plt.grid()


plt.subplot(2, 2, 3)
plt.title('c) Random')
plt.plot(MCC[:],M[1][:], 'b',  label="T = 1.0")
plt.plot(MCC[:],M[3][:], 'r',  label="T = 2.4")
plt.ylabel(r"$\langle |M|\rangle$")
plt.xlabel('Monte Carlo cycles')
plt.xscale('log')
plt.legend()
plt.grid()


plt.subplot(2, 2, 4)
plt.title('d) Ordered')
plt.plot(MCC[:],M[0][:], 'b',  label="T = 1.0")
plt.plot(MCC[:],M[2][:], 'r',  label="T = 2.4")
plt.ylabel(r"$\langle |M|\rangle$")
plt.xlabel('Monte Carlo cycles')
plt.xscale('log')
plt.legend()
plt.grid()

plt.savefig("taskc_1.png")
"""