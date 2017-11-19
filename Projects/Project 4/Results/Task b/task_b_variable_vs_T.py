import numpy as np
import matplotlib.pyplot as plt

#########################
## Declaring variables ##
dt = 0.1 # Temperature step
N = int(2.1/0.1)
T = np.linspace(1,3,N)
number = [100,1000,10000,100000,1000000,10000000] # Runs of program with different cycles

E_tot = np.zeros((6, N))
E_tot_a = np.zeros(N)
M_tot = np.zeros((6, N))
M_tot_a = np.zeros(N)
C_v = np.zeros((6, N))
C_v_a = np.zeros(N)
chi = np.zeros((6, N))
chi_a = np.zeros(N)


###########################
## Analytical solutionss ##
def analyticalE(T):
    eightJB = 8/T
    E = 8 * (np.exp(-eightJB) - np.exp(eightJB))/(6 + np.exp(eightJB) + np.exp(-eightJB))
    return E

def analyticalEsqr(T):
    eightJB = 8/T
    E_squared = ((8)*(8))*(np.exp(-eightJB) + np.exp(eightJB)) / (6 + np.exp(eightJB) + np.exp(-eightJB))
    return E_squared

def analyticalM(T):
    M = 4/(2* (6+ np.exp(-8/(T)) + np.exp(8/(T)) )) * (2 * np.exp(8/T) + 4)
    return M
    
def analyticalMsqr(T):
    eightJB = 8/T
    M_squared = 16*(np.exp(eightJB) + 1)/(6 + np.exp(eightJB) + np.exp(-eightJB))
    return M_squared

def analyticalC_v(T):
    C_v = (( analyticalEsqr(T) - ( np.power((analyticalE(T)),2) )) / (T*T) )
    return C_v

def analyticalChi(T):
    chi = (analyticalMsqr(T)- (np.power(analyticalM(T),2)))/T;
    return chi

##########################
## Importing variables  ##
for i in range(6):
    cycles = str(number[i])
    filename = "L2_"+cycles+".dat"
    infile = open(filename, "r")
    infile.readline()
    for j in range(N):
        line = infile.readline()
        variables = line.split()
        E_tot[i][j] = variables[0]
        M_tot[i][j] = variables[1]
        C_v[i][j] = variables[2]
        chi[i][j] = variables[3]
infile.close()

#########################
## Analytical arrays   ##
for i in range(N):
        E_tot_a[i] = analyticalE(T[i])/4
        M_tot_a[i] = analyticalM(T[i])/4
        C_v_a[i] = analyticalC_v(T[i])/4
        chi_a[i] = analyticalChi(T[i])/4

################
## Formatting ##
def format_e(n):
    a = '%E' % n
    return a.split('E')[0].rstrip('0').rstrip('.') + 'E' + a.split('E')[1]


#########################
## Relative error plot ##
a = np.zeros(6)
b = np.zeros(6)
c = np.zeros(6)
d = np.zeros(6)

for i in range(6):
    a[i]= abs(1-((E_tot[i][10])/(E_tot_a[10])))
    b[i]= abs( 1 - ( (M_tot[i][10]) / (M_tot_a[10]) ) )
    c[i] = abs( 1 - ( (C_v[i][10]) / (C_v_a[10]) ) )
    d[i] = abs( 1 - ( (chi[i][10]) / (chi_a[10]) ) )
    
print(analyticalChi(1)/4)
"""
plt.figure(num=None, figsize=(16, 10), dpi=120, facecolor='w', edgecolor='k')
plt.rcParams.update({'font.size': 16})
    
plt.subplot(2, 2, 1)
plt.title('a) Mean energy')
plt.plot(number[:],a[:], '-bo',  label=r"$\langle E\rangle$")
plt.ylabel('Relative error')
plt.xscale('log')
plt.legend()
plt.grid()
 

plt.subplot(2, 2, 2)
plt.title('b) Mean magnetization')
plt.plot(number[:],b[:], '-bo',  label=r"$\langle |M|\rangle$")
plt.ylabel('Relative error')
plt.xscale('log')
plt.legend()
plt.grid()


plt.subplot(2, 2, 3)
plt.title('c) Heat capacity')
plt.plot(number[:],c[:], '-bo',  label='$C_v$')
plt.ylabel('Relative error')
plt.xlabel('Monte Carlo cycles')
plt.xscale('log')
plt.legend()
plt.grid()


plt.subplot(2, 2, 4)
plt.title('d) Susceptibility')
plt.plot(number[:],d[:], '-bo',  label='$\chi$')
plt.ylabel('Relative error')
plt.xlabel('Monte Carlo cycles')
plt.xscale('log')
plt.legend()
plt.grid()

plt.savefig("relerr_T10.png")
"""
######################################
## Plotting numerical vs analytical ##
"""
for i in range(6):
    plt.figure(num=None, figsize=(16, 10), dpi=120, facecolor='w', edgecolor='k')
    plt.rcParams.update({'font.size': 16})
    
    cycles = str(format_e(number[i]))
    #legend = cycles + " cycles"
    
    
    plt.subplot(2, 2, 1)
    plt.title('a) Mean energy')
    plt.plot(T[:],E_tot[i], 'bo',  label="Numerical")
    plt.plot(T[:],E_tot_a[:],'r', label="Analytical")
    plt.ylabel(r'$\langle E\rangle$')
    plt.legend()
    plt.grid()

    plt.subplot(2, 2, 2)
    plt.title('b) Mean magnetization')
    plt.plot(T[:],M_tot[i], 'bo',  label="Numerical")
    plt.plot(T[:],M_tot_a[:],'r', label="Analytical")
    plt.ylabel(r'$\langle |M|\rangle$')
    plt.legend()
    plt.grid()

    plt.subplot(2, 2, 3)
    plt.title('c) Heat capacity')
    plt.plot(T[:],C_v[i], 'bo',  label="Numerical")
    plt.plot(T[:],C_v_a[:],'r', label="Analytical")
    plt.xlabel('T')
    plt.ylabel(r'$C_v$')
    plt.legend()
    plt.grid()

    plt.subplot(2, 2, 4)
    plt.title('d) Susceptibility')
    plt.plot(T[:],chi[i], 'bo',  label="Numerical")
    plt.plot(T[:],chi_a[:],'r', label="Analytical")
    plt.xlabel('T')
    plt.ylabel(r'$\chi_{abs}$')
    plt.legend()
    plt.grid()
    
    filename = "L2_"+cycles+".png"
    
    plt.savefig(filename)

"""