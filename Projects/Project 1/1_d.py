from numpy import *
from matplotlib.pyplot import *
import time

# Step size
n = 10
h = 1./(n+1)  


# Filling vectors for diagonals
b = zeros(n+1)
b.fill(-1)
a_t = zeros(n+1)
a_t.fill(2)
c = zeros(n+1)
c.fill(-1)


# General Gaussian elimination algorithm for tridiagonal matrix
def GeneralSolver(n, b, a_t, c, d):
# Forward substitution
	for i in range(n-1):  
		a_t[i+1] -= (b[i]*c[i]) / a_t[i]
		d[i+1] -= (b[i]*d[i]) / a_t[i]
# Backwards substitution
	for i in range(n-1,1,-1): 
		d[i-1] -= (c[i-1]*d[i]) / a_t[i]
# Calculating solution
	for i in range(1,n):
		d[i] = d[i] / a_t[i]
	return d


# Computing d vector from source function
x = linspace(0,1,n+1)
f = 100*exp(-10*x)
d = f*(h**2)
d[0] = d[n] = 0


print("Grid points: n = %d" %(n))
ans = GeneralSolver(n, b, a_t, c, d)

# Error computation
u = 1-x*(1-exp(-10))-exp(-10*x)
error = zeros(n+1)
for i in range(1,n-2):
	error[i] = (abs((ans[i]-u[i]) / u[i]))
max_error = log10(abs(1-max(error)))
print("Maximum error [log10]: %f" %(max_error))

epsilon=[-0.210019,-1.218700,-2.221513,-3.221815,-4.221843,-5.247077,-5.800997]
h=[1,2,3,4,5,6,7]


# Plotting
plot(h, epsilon,'b')
xlabel("log10(h)")
ylabel("$log10(\epsilon)$", rotation=0, position=(0,0.55))
legend(["Maximum error"])
title("Maximum error as a function of step size")
savefig('h_epsilon.png')
show()