#  Code for the MarineEco project

import numpy as np
import scipy as sc
from scipy.integrate import odeint
from scipy import integrate

# Set up the solver for the system of differential equations.
#=======================================================
def MarineDF1(par1,initial_cond,start_t,end_t,incr):
     #-time-grid-----------------------------------
     t  = np.linspace(start_t, end_t, num = incr)
     #differential-eq-system----------------------
     # This needs fixed to represent the system we are looking at
     def funct(x,t):
        r1, r2, r3, alpha12, alpha13, alpha21, alpha23, alpha31, alpha32, beta24, beta35, c4, c5, K1, K2, K3 = par1
        xp1 = r1*x[0] + r1*x[0]*( (x[0] + alpha12*x[1] + alpha13*x[2])/K1 )
        xp2 = r2*x[1] + r2*x[1]*( (x[1] + alpha21*x[0] + alpha23*x[2])/K2 ) - beta24*x[1]*x[3]
        xp3 = r3*x[2] + r3*x[2]*( (x[2] + alpha31*x[0] + alpha32*x[1])/K3 ) - beta35*x[2]*x[4]
        xp4 = beta24*x[1]*x[3] - c4*x[3]
        xp5 = beta35*x[2]*x[4] - c5*x[4]
        return [xp1, xp2, xp3, xp4, xp5]
     #integrate------------------------------------
     ds = integrate.odeint(funct, initial_cond, t)
     return (t, ds[:,0], ds[:,1], ds[:,2], ds[:,3], ds[:,4])
#=======================================================

# Create the log Normal density on the log scale.
def LogNormalLn( x1, mean1, sd1 ):
    out1 = -np.log(x1*sd1*np.sqrt(2.0*np.pi)) - 0.5*( (np.log(x1) - mean1)/sd1 )**2
    return out1

# Create the chi square distribution on the Log scale    
def chisqln( x1, df1 ):
    out1 = -df1/2*np.log( 2 ) - sc.gammaln( df1/2 )
    out1 = out1 + ( df1/2 - 1 )*np.log( x1 ) - x1/2
    return out1

# Create the log-likelihood for the marine eco data
def MarineLk1( X1, par1, variance1, initial_cond, start_t, end_t, incr ):
    sig1, sig2, sig3, sig4, sig5 = variance1
    XP1 = MarineDF1(par1, initial_cond, start_t, end_t, incr)
    # solve the system first with the current parameters
    # Measure the likelihood of the data under the current parameters
    
    #define XP2 as the array of observation times for the integrated model
    # NOTE: this assumes the first observation is at time 1 and the model
    #begins at time 0.  For the case the first observation is at time zero
    #and model begins time zero, use XP1[::1.0/incr,:] instead 
    XP2 = XP1[1.0/incr::1.0/incr,:]
    
    #Note: I changed the index of the arrays to begin with entry 0
    res1 = 0.0
    res1 = np.sum(LogNormalLn( X1[:,0], XP2[:,0], sig1 ))          # For X1
    res1 = res1 + np.sum(LogNormalLn( X1[:,1], XP2[:,1], sig2 ))   # For X2
    res1 = res1 + np.sum(LogNormalLn( X1[:,2], XP2[:,2], sig3 ))   # For X3
    res1 = res1 + np.sum(LogNormalLn( X1[:,3], XP2[:,3], sig4 ))   # For X4
    res1 = res1 + np.sum(LogNormalLn( X1[:,4], XP2[:,4], sig5 ))   # For X5
    return res1

# Create the log-posterior (up to proportionality constant) for the marine eco data
def MarinePst1( X1, par1, variance1, par2m, par2v, variance1df, initial_cond, start_t, end_t, incr):
    r1, r2, r3, alpha12, alpha13, alpha21, alpha23, alpha31, alpha32, beta24, beta35, c4, c5, K1, K2, K3 = par1
    r1m, r2m, r3m, alpha12m, alpha13m, alpha21m, alpha23m, alpha31m, alpha32m, beta24m, beta35m, c4m, c5m = par2m
    r1v, r2v, r3v, alpha12v, alpha13v, alpha21v, alpha23v, alpha31v, alpha32v, beta24v, beta35v, c4v, c5v = par2v
    sig1, sig2, sig3, sig4, sig5 = variance1
    df1, df2, df3, df4, df5 = variance1df
    # Obtain the likelihood of the data.
    res1 = MarineLk1( X1, par1, variance1, initial_cond, start_t, end_t, incr )
    # Start putting the priors on
    # In this case K1, K2 and K3 are fixed as I think bad things will happen otherwise
    res1 = res1 + LogNormalLn( r1, r1m, r1v )
    res1 = res1 + LogNormalLn( r2, r2m, r2v )
    res1 = res1 + LogNormalLn( r3, r3m, r3v )
    res1 = res1 + LogNormalLn( alpha12, alpha12m, alpha12v )
    res1 = res1 + LogNormalLn( alpha13, alpha13m, alpha13v )
    res1 = res1 + LogNormalLn( alpha21, alpha21m, alpha21v )
    res1 = res1 + LogNormalLn( alpha23, alpha23m, alpha23v )
    res1 = res1 + LogNormalLn( alpha31, alpha31m, alpha31v )
    res1 = res1 + LogNormalLn( alpha32, alpha32m, alpha32v )
    res1 = res1 + LogNormalLn( beta24, beta24m, beta24v )
    res1 = res1 + LogNormalLn( beta35, beta35m, beta35v )
    res1 = res1 + LogNormalLn( c4, c4m, c4v )
    res1 = res1 + LogNormalLn( c5, c5m, c5v )
    res1 = res1 + chisqln( sig1**2, df1 )
    res1 = res1 + chisqln( sig2**2, df2 )
    res1 = res1 + chisqln( sig3**2, df3 )
    res1 = res1 + chisqln( sig4**2, df4 )
    res1 = res1 + chisqln( sig5**2, df5 )
    return res1
    
