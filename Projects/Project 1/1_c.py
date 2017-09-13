from numpy import *
from matplotlib.pyplot import *
import time

#Step size
n = 1000000
h = 1./(n+1)


# Pre calculating a_tilde values
a_t = zeros(n)
for i in range(1,n,1):
	a_t[i] = (i+1)/(i)
a_t[0] = 2


# Generating vectors (d, d tilde and u)
x = linspace(0,1,n+1)
f = 100*exp(-10*x)
d = (h**2)*f
d_t = zeros(n)
d_t[0] = d[0]
u = zeros(n+1)


# Specialized Gaussian elimination algorithm for tridiagonal matrix
def SpecialSolver(a_t, d):
	
	for i in range(1,n,1):
		d_t[i] = d[i] + d_t[i-1] / a_t[i-1]
        
	
	u[0] = 0
	u[-1] = 0
	u[-2]= d_t[-1] / a_t[-1]

	for i in range(n-2, 0, -1):
		u[i] = (d_t[i]+u[i+1]) / a_t[i]
	return u


# Computation time
t_start = time.clock()
ans = SpecialSolver(a_t, d)
t_stop = time.clock()
t_diff = t_stop-t_start
print("Grid points: n = %d" %(n))
print("Computation time: %f s" %(t_diff))


# Error computation
u_ex = 1-x*(1-exp(-10))-exp(-10*x)
error = zeros(n+1)
for i in range(2,n-2):
	error[i] = abs((ans[i]-u_ex[i]) / u_ex[i])
max_error = abs(max(error))
print("Maximum error: %f" %(max_error))


# Plotting
plot(x, ans,'b', x, u_ex, 'r')
xlabel("x")
ylabel("u(x)", rotation=0, position=(0,0.47))
legend(["Numerical, n=%d" %(n), "Analytical"])
title("Special alg. for solving lin. equations using Gaussian elimination")
#savefig('1c_1000.png')
show()