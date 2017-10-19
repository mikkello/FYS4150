import numpy as np
import matplotlib.pyplot as plt

infile = open("variables1.txt", "r")

infile.readline()
line = infile.readline()
variables = line.split()
steps = int(variables[0])
num_pos = int(steps/2)
num_planets = int(variables[1])
num_years = int(variables[2])
dt = float(variables[3])
infile.close()

infile = open("positions1.xyz", "r")

x = np.zeros((num_planets, steps))
y = np.zeros((num_planets, steps))
z = np.zeros((num_planets, steps))

for dim in range(steps):
    for i in range(num_planets):
        line = infile.readline()
        words = line.split()
        x[i][dim] = words[0]
        y[i][dim] = words[1]
        z[i][dim] = words[2]
        
infile.close()


plt.figure(num=None, figsize=(12, 10), dpi=120, facecolor='w', edgecolor='k')
plt.rcParams.update({'font.size': 20})

if num_planets == 2:
    plt.plot(x[1], y[1], label='earth')
    plt.plot(x[0], y[0],'ro', label='sun')
    plt.xlabel('x [AU]')
    plt.ylabel('y [AU]')
    plt.legend()
    plt.title('Earth-Sun system using Forward Euler (dt = %s)' %(dt) )
    plt.savefig('earth_sun_dt_%s.png'%(dt))
    plt.show()
