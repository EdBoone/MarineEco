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
     def funct(x,t):
        D1, si, mu1, mu2, y1, y2, k1, k2, epsilon1, T1 = par1
        omega1 = 2*np.pi/T1
        xp1 = D1*( si*(1 + epsilon1*np.sin(omega1*t))-x[0])-mu1*x[0]*x[1]/y1/(k1+x[0])
        xp2 = mu1*x[0]*x[1]/(k1+x[0]) - D1*x[1] - mu2*x[1]*x[2]/y2/(k2+x[1])
        xp3 = mu2*x[1]*x[2]/(k2+x[1]) - D1*x[2]
        return [xp1, xp2, xp3]
     #integrate------------------------------------
     ds = integrate.odeint(funct, initial_cond, t)
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
x1 = np.zeros(shape=(1001))
x2 = np.zeros(shape=(1001))
x3 = np.zeros(shape=(1001))
t1 = np.zeros(shape=(1001))

# Initial Conditions
x0 = [0.25, 0.5, 0.5]

# Pack the parameters together
par1 = (D1, si, mu1, mu2, y1, y2, k1, k2, epsilon1, T1)
t1[:], x1[:], x2[:], x3[:] = KotDF1( par1, x0, 0, 45., 1001 )
 
plt.figure()
plt.plot(t1, x1,'-b', t1 , x2,'-r', t1, x3, '-k')
plt.ylim( [0,100])
#plt.legend(('Foxes', 'Rabbits'),'upper center',ncol=2)
title1 = 'It works???'
plt.title(title1)
 
 
D1 = 0.1        # Dilution rate 1/hour
si = 115        # units of si are mg/l
# maximum specific growth rate of prey and preditor (1/hour)
mu1 = 0.51       # prey
mu2 = 0.19       # preditor
# yields
y1 = 0.1        # yield of prey per unit mass of substrate
y2 = 0.59        # biomass yield of predator per unit mass of prey
# half-saturation (Michaelis-Menton)
k1 = 8.1          # prey
k2 = 9.1          # predator
# Other variables
epsilon1 = 0.7
T1 = 15

# Create containers
x1a = np.zeros(shape=(1001))
x2a = np.zeros(shape=(1001))
x3a = np.zeros(shape=(1001))
t1a = np.zeros(shape=(1001))

# Initial Conditions
x0 = [0.25, 0.5, 0.5]

# Pack the parameters together
par1 = (D1, si, mu1, mu2, y1, y2, k1, k2, epsilon1, T1)
t1a[:], x1a[:], x2a[:], x3a[:] = KotDF1( par1, x0, 0, 50., 1001 )

s1s = range( 0, 75, 5 )
s2s = range( 0, 275, 5)
x1d = np.random.poisson(x1a)
x2d = np.random.poisson(x2a)

  
plt.figure()
plt.scatter(t1a[s1s], x1d[s1s])
plt.plot( t1 , x1, '-r' )
plt.ylim( [0,100] )
plt.xlim( [0,45] )
plt.vlines(max(t1a[s1s]),0,100, linestyles='dashed', colors ='r')
#plt.legend(('Foxes', 'Rabbits'),'upper center',ncol=2)
title1 = 'Data Assimilate: Step 1 \n x1'
plt.title(title1)
plt.savefig('DSStep1x1.pdf', format='pdf')

plt.figure()
plt.scatter(t1a[s1s], x2d[s1s])
plt.plot( t1 , x2, '-r' )
plt.ylim( [0,100] )
plt.xlim( [0,45] )
plt.vlines(max(t1a[s1s]),0,100, linestyles='dashed', colors ='r')
#plt.legend(('Foxes', 'Rabbits'),'upper center',ncol=2)
title1 = 'Data Assimilate: Step 1 \n x2'
plt.title(title1)
plt.savefig('DSStep1x2.pdf', format='pdf') 
 
plt.figure()
plt.scatter(t1a[s2s], x1d[s2s])
plt.plot(t1, x1, '-r', t1a , x1a, '-b' )
plt.ylim( [0,100] )
plt.xlim( [0,45] )
plt.vlines(max(t1a[s1s]),0,100, linestyles='dashed', colors ='r')
plt.vlines(max(t1a[s2s]),0,100, linestyles='dashed', colors ='b')
#plt.legend(('Foxes', 'Rabbits'),'upper center',ncol=2)
title1 = 'Data Assimilate: Step 2 \n x1'
plt.title(title1)
plt.savefig('DSStep2x1.pdf', format='pdf')

plt.figure()
plt.scatter(t1a[s2s], x2d[s2s])
plt.plot( t1, x2, '-r', t1a , x2a, '-b',  )
plt.ylim( [0,100] )
plt.xlim( [0,45] )
plt.vlines(max(t1a[s1s]),0,100, linestyles='dashed', colors ='r')
plt.vlines(max(t1a[s2s]),0,100, linestyles='dashed', colors ='b')
#plt.legend(('Foxes', 'Rabbits'),'upper center',ncol=2)
title1 = 'Data Assimilate: Step 2 \n x2'
plt.title(title1)
plt.savefig('DSStep2x2.pdf', format='pdf') 