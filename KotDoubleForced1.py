# Simple Kot Models

import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import odeint
from scipy import integrate
#import os

#os.chdir('C:\Users\ed\Google Drive\MarineEco')
 
#=======================================================
def KotDF1(par1,initial_cond,start_t,end_t,incr):
     #-time-grid-----------------------------------
     t  = np.linspace(start_t, end_t, num = incr)
     #differential-eq-system----------------------
     def funct(t,x):
        D1, si, mu1, mu2, y1, y2, k1, k2, epsilon1, T1 = par1
        omega1 = 2*pi/T1
        xp1 = D1*( si*(1 + epsilon1*sin(omega1*t))-x[0])-mu1*x[0]*x[1]/y1/(k1+x[0])
        xp2 = mu1*x[0]*x[1]/(k1+x[0]) - D1*x[1] - mu2*x[1]*x[2]/y2/(k2+x[1])
        xp3 = mu2*x[1]*x[2]/(k2+x[1]) - D1*x[2]
        return [xp1, xp2, xp3]
     #integrate------------------------------------
     ds = integrate.odeint(funct, t, initial_cond)
     return (t, ds[:,0], ds[:,1], ds[:,2])
#=======================================================


D1 = 0.1        # Dilution rate 1/hour
si = 115        # units of si are mg/l
# maximum specific growth rate of prey and preditor (1/hour)
mu1 = 0.5       # prey
mu2 = 0.2       # preditor
# yields
y1 = 0.4        # yield of prey per unit mass of substrate
y2 = 0.6        # biomass yield of predator per unit mass of prey
# half-saturation (Michaelis-Menton)
k1 = 8          # prey
k2 = 9          # predator
# Other variables
epsilon1 = 0.6
T1 = 15

# Create containers
x1 = np.zeros(shape=(101))
x2 = np.zeros(shape=(101))
x3 = np.zeros(shape=(101))
t1 = np.zeros(shape=(101))

# Initial Conditions
x0 = [0.25, 0.5, 0.5]

# Pack the parameters together
par1 = (D1, si, mu1, mu2, y1, y2, k1, k2, epsilon1, T1)
t1[:], x1[:], x2[:], x3[:] = KotDF1( par1, x0, 0, 15., 101 )
 
plt.figure()
plt.plot(Tout[:,i],F0out[:,i],'-b',Tout[:,i],F1out[:,i],'-r')
plt.ylim( [0,350])
plt.legend(('Foxes', 'Rabbits'),'upper center',ncol=2)
title1 = 'b = ' + str(np.round(b1,4))
plt.title(title1)
 