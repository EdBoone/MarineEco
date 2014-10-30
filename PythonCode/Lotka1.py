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
#################################################################

# Simple Lotka-Volterra Model
import numpy as np
from scipy.integrate import odeint
from scipy import integrate
import scipy.stats as scs
#from numba import jit
 
#@jit
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

#@jit
#=======================================================
#  Poisson log likelihood
def LVLogLike1(P01,P11,par,initial_cond,start_t,end_t,incr):
    F0out,F1out,Tout = LV1(par,initial_cond,start_t,end_t,incr)
    out1 = 0.0   
    l1 = scs.poisson.logpmf(P01,F0out)
    l2 = scs.poisson.logpmf(P11,F1out)
    out1 = np.sum(l1)+np.sum(l2)
    return out1

#@jit    
#=======================================================    
def LVLogPost1(P0,P1,x1_data,a1,b0,b1,g1,d1,initial_cond,start_t,end_t,incr):
    r1 = P0.shape[1]
    lpost1 = 0
    for i in range(r1):
        b0x1 = b0 + b1*(x1_data[i]-10)        
        b1x = np.exp(b0x1)  
        P01 = P0[ : , i ]
        P11 = P1[ : , i ]  
        rates=(a1,b1x,g1,d1)
        lpost1 =  lpost1 + LVLogLike1(P01,P11,rates,initial_cond,start_t,end_t,incr)
    lpost1 = lpost1 + scs.norm.logpdf(b0, 0, 1)
    lpost1 = lpost1 + scs.norm.logpdf(b1, 0, 1)
    lpost1 = lpost1 + scs.norm.logpdf(a1, 0, 1)
    lpost1 = lpost1 + scs.norm.logpdf(g1, 0, 1)
    lpost1 = lpost1 + scs.norm.logpdf(d1, 0, 1)    
    return lpost1