#  Code for the MarineEco project

from pylab import *
import matplotlib.pyplot as plt
from scipy.integrate import odeint
from scipy import integrate

###############################################################################
# Define the time derivative for the competition model

def Comp_D(initial_cond,t,A,B,C,K,R):
    """This is the time derivative for the competition predator-prey model
     
    The inputs are the array of initial conditions, a dummy time parameter
    for the ODE solver, and the model parameters.  Inputs are of the form

    x1  - coral
    x2  - macro algae
    x3  - turf
    x4  - browser
    x5  - grazer/ detritivore
    
    a12 - competition of x2 on x1 
    a13 - competition of x3 on x1
    a21 - competition of x1 on x2
    a23 - competition of x3 on x2
    a31 - competition of x1 on x3
    a32 - competition of x2 on x3
     
    b24 - predation success of x4 on x2
    b35 - predation sucesss of x5 on x3
     
    c4  - mortality rate of x4
    c5  - mortality rate of x5
     
    k1  - carrying capacity of x1, up to 65%
    k2  - carrying capacity of x2, up to 100%     
    k3  - carrying capacity of x3, up to 60%
     
    Comp_D returns the array [dx1,dx2,dx3,dx4,dx5], the time derivative at
    the initial condition."""

    #Unpack the initial conditions and the parameters
    [x1, x2, x3, x4, x5]      = initial_cond
    
    [a12,a13,a21,a23,a31,a32] = A
    
    [b24,b35]                 = B
    
    [c4,c5]                   = C
    
    [k1,k2,k3]                = K
     
    [r1,r2,r3]                = R
     
    #Calculate the derivatives
    dx1 = r1*x1 - r1*x1*(x1+a12*x2+a13*x3)/k1
    
    dx2 = r2*x2 - r2*x2*(x2+a21*x1+a23*x3)/k2 - b24*x2*x4
    
    dx3 = r3*x3 - r3*x3*(x3+a31*x1+a32*x2)/k3 - b35*x3*x5
    
    dx4 = b24*x2*x4 - c4*x4
    
    dx5 = b35*x3*x5 - c5*x5
    
    return array([dx1,dx2,dx3,dx4,dx5])
###############################################################################
# Create the log Normal density on the log scale.

def LogNormalLn( x1, mean1, sd1 ):
    res1 = 1
    return res1
    
###############################################################################
#Here we define the log scale log-normal 
#likelyhood function for the parameter estimates

def LN(mean,variance,sample):
    
    """This function evaluates the likelyhood of a given sample
    
    This is a vectorized parameter likelyhood function which takes an array of
    data (mean), the associated array of variances (variance), and the
    parameter dependent outcomes (sample)
    
    Inputs are of the form
    
    number of observed variables = shape(mean), shape(variance), shape(sample)

    Output is an array on lognormal scale of the likelyhood of each outcome
    (sample) on given the associated mean and variance of the data."""
    
    #define the log of the denominator
    # with multiplication across arrays with Schur product
    denom  = -log(sample*variance*sqrt(2.0*pi))
    #define the numerator in log scale
    #again operating on each array element individually
    numer  = -(log(sample) - mean)**2/(2.0*variance**2)
    
    #calclulate the likelyhood
    output = numer + denom
    
    return(output)

###############################################################################
# Create the log-scale likelihood for the marine eco data
#likelyhood distribution is by default log-normal 

def MarineLk1( X1, par1, initial_cond, start_t, end_t, incr,likely=LN ):
    # solve the system first with the current parameters
    # Measure the likelihood of the data under the current parameters
    
    
    
    return res1

    
###############################################################################
# Create the log-posterior (up to proportionality constant)
#for the marine eco data

def MarinePst1( X1, par1, par2, initial_cond, start_t, end_t, incr):
    res1 = 1
    return res1

###############################################################################