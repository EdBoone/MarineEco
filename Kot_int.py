# Simple Kot Models

import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import odeint
from scipy import integrate
#import os

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
