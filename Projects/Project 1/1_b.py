from numpy import *
from matplotlib.pyplot import *
import time

# Step size
n = 100
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


# Computation time
t_start = time.clock()
ans = GeneralSolver(n, b, a_t, c, d)
t_stop = time.clock()
t_diff = t_stop-t_start
print("Grid points: n = %d" %(n))
print("Computation time: %f s" %(t_diff))


# Error computation
u = 1-x*(1-exp(-10))-exp(-10*x)
error = zeros(n+1)
for i in range(1,n-1):
	error[i] = abs((ans[i]-u[i]) / u[i])
max_error = abs(1-max(error))
print("Maximum error: %f" %(max_error))


# Plotting
plot(x, ans,'b', x, u, 'r')
xlabel("x")
ylabel("u(x)", rotation=0, position=(0,0.47))
legend(["Numerical, n=%d" %(n), "Analytical"])
title("Gen. alg. for solving lin. equations using Gaussian elimination")
#savefig('1b_1000.png')
show()