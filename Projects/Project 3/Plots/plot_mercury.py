import numpy as np
import matplotlib.pyplot as plt
#########################
## Importing variables ##
infile = open("Variables_7_dt_0.000001.txt", "r")
infile.readline()
line = infile.readline()
variables = line.split()
N_steps = int(variables[0]) # Number of steps
N_planets = int(variables[1]) # Number of planets
N_years = int(variables[2]) # Number of years
dt = float(variables[3]) # Time step
infile.close()
time = []
"""
# Newtonian force
infile = open("Perihelion_6_dt_0.000100_f.txt", "r")
Perihilion_6_dt_100k = []
infile.readline()
for line in infile:
    words = line.split()
    Perihilion_6_dt_100k.append((float(words[0])))
    """
    
infile.close()

# Relativistic correction
infile = open("Perihelion_7_dt_0.000001.txt", "r")
Perihilion_7_dt_100k = []
infile.readline()
for line in infile:
    words = line.split()
    Perihilion_7_dt_100k.append(float(words[0]))
    time.append((float(words[1]))/0.240846)
infile.close()

plt.figure(num=None, figsize=(14, 10), dpi=120, facecolor='w', edgecolor='k')
plt.rcParams.update({'font.size': 16})
plt.plot(time, Perihilion_7_dt_100k)
plt.xlabel("Time [Mercury years]")
plt.ylabel("Angle [Arcseconds]")
plt.title('Perihelion precession of Mercury.\n dt = %s' %(dt))
plt.legend(["Newtonian", "Relativistic"])
plt.grid()
plt.savefig('mercury_perihelion_comparison_dt_%s' %(dt))
plt.show()

