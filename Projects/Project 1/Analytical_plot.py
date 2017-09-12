from numpy import *
from matplotlib.pyplot import *

def u(x):
    return 1-(1-exp(-10))*x-exp(-10*x)

x = linspace(0, 1, 1000)    
y = zeros(len(x))  
       
for i in range(len(x)):
    y[i] = u(x[i])

plot(x, y, 'r')
xlabel('x')
ylabel('u(x)',rotation=0)
legend(['$u(x)=1-(1-e^{-10})x-e^{-10x}$'])
title('Analytical solution')
savefig('analytical.png')
show()
