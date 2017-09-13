import numpy as np
import time
import scipy.linalg

n = 100000
h = 1./(n)

# Filling diagonal vectors
a = zeros(n-1)
a.fill(-1)
b = zeros(n)
b.fill(2)
c = zeros(n-1)
c.fill(-1)


# Generating Toeplitz matrix
def tridiag(a, b, c, d1=-1, d2=0, d3=1):
    return np.diag(a, d1) + np.diag(b, d2) + np.diag(c, d3)
A = tridiag(a, b, c)


# Computing right hand side of equation and exact solution
x = np.linspace(0,1,n)
f = 100*np.exp(-10*x)
u = (h**2)*f
u_ex = 1-x*(1-exp(-10))-exp(-10*x)


# Computation time
t_start = time.clock()
LU = scipy.linalg.lu_factor(A) # lU decomposition of the Toeplitz matrix
ans = scipy.linalg.lu_solve(LU, u) # Solving equation using L and U matrices
t_stop = time.clock()
t_diff = t_stop-t_start


print("Grid points: n = %d" %(n))
print("Computation time: %f s" %(t_diff))


plot(x, ans,'b', x, u_ex, 'r')
