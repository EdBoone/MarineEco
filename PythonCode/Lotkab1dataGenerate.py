#################################################################
#################################################################
#
#  Properties of this Lotka-Volterra model.
#
# Fixed points
# x = 0 and y = 0
#  and
# y = a1/b1  and x = g1/d1
#
# Jacobian
#
# Stability of fixed points
#
# J(x,y) = [ a1 - b1*y      -b1*x
#            d1*y          d1*x -g1]
#
#
#  Under first fixed point
#
# J(0,0) = [ a1    0
#            0    -g1]
#
#  Eigen values are a1 and -g1
#
#  Under second fixed point
#
# J(g1/d1, a1/b1) = [ 0        -b1*g1/d1
#                    a1*d1/b1      0     ]
#
#  Eigen values are i*sqrt(a1*g1) and -i*sqrt(a1*g1)
#  Oscillations around fixed points
#
#

# Simple Lotka-Volterra Model
import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import odeint
from scipy import integrate
import os

os.chdir('C:\Users\ed\Google Drive\MarineEco')
 
#=======================================================
def LV1(par,initial_cond,start_t,end_t,incr):
     #-time-grid-----------------------------------
     t  = np.linspace(start_t, end_t,incr)
     #differential-eq-system----------------------
     def funct(y,t):
        Fi=y[0]
        Ri=y[1]
        a1, b1, g1, d1 = par
        # the model equations (see Munz et al. 2009)
        #dx = x*(alpha - beta*y)
        #dy  = -y*(gamma - delta*x)
        f0 = Fi*(a1 - b1*Ri)
        f1 = -Ri*(g1 - d1*Fi)
        return [f0, f1]
     #integrate------------------------------------
     ds = integrate.odeint(funct,initial_cond,t)
     return (ds[:,0],ds[:,1],t)
#=======================================================
   
#parameters   
a1 = 1.0   # prey natural reproduction rate  a1 > 0
b1 = 0.01  # rate of predation  b1 > 0
g1 = 1.0   # Preditor death rate  g1 > 0
d1 = 0.02  # death rate of preditor  d1 > 0

# initial conditions
F0 = 10               #  initial Fox population
R0 = 100              #  initial Rabbit population
y0 = [F0, R0]         #  initial condition vector

#=======================================================
#  Create the gradient to influence b1
x1 = np.linspace(start=10, stop=20, num=11 )

#  Create the parameter gradient in logscale.
b1log = ( (np.log(0.02) - np.log(0.01))/10*(x1-10) + np.log(0.01) 
    + np.random.normal(0,0.01,11) )

# initial conditions
F0 = 10               #  initial Fox population
R0 = 100              #  initial Rabbit population
y0 = [F0, R0]         #  initial condition vector

# Create containers for the output
F0out = np.zeros(shape=(101,11))
F1out = np.zeros(shape=(101,11))
Tout = np.zeros(shape=(101,11))
P0out = np.zeros(shape=(101,11))
P1out = np.zeros(shape=(101,11))

# Fit the models across the      
for i in xrange(0,11):
    b1 = np.exp(b1log[i])    
    rates=(a1,b1,g1,d1)
    # Fit the model to create the means 
    F0out[:,i],F1out[:,i],Tout[:,i] = LV1(rates,y0,0,100.,101)
    # Create counts
    P0out[:,i] = np.random.poisson(lam=F0out[:,i])
    P1out[:,i] = np.random.poisson(lam=F1out[:,i])
    # Create the plots for the means
    plt.figure()
    plt.plot(Tout[:,i],F0out[:,i],'-b',Tout[:,i],F1out[:,i],'-r')
    plt.ylim( [0,350])
    plt.legend(('Foxes', 'Rabbits'),'upper center',ncol=2)
    title1 = 'b = ' + str(np.round(b1,4))
    plt.title(title1)
    pathname1 = 'LV' + str(i) +'mean.pdf'
    #plt.show()
    plt.savefig(pathname1, format='pdf')
    plt.close()
    plt.figure()
    plt.plot(Tout[:,i],P0out[:,i],'-b',Tout[:,i],P1out[:,i],'-r')
    plt.ylim( [0,350])
    plt.legend(('Foxes', 'Rabbits'),'upper center',ncol=2)
    title1 = 'b = ' + str(np.round(b1,4))
    plt.title(title1)
    pathname1 = 'LV' + str(i) +'value.pdf'
    #plt.show()
    plt.savefig(pathname1, format='pdf')
    plt.close()

# Save the datasets
data1out = np.concatenate( (Tout,P0out,P1out), axis= 1)
np.savetxt("Lotka1b1data.csv", data1out, delimiter=",")
np.savetxt("Lotka1x1data.csv", x1, delimiter=",")