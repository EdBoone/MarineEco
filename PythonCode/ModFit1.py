#####################################################################
#####################################################################
###
###  Fit the Lotka-Volterra Model
###
#####################################################################
#####################################################################

# Simple Lotka-Volterra Model
import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import odeint
from scipy import integrate
import Lotka1 as lv1
import scipy.stats as scs
import time as ts

start_time = ts.time()

# Set the number of MCMC Samples
MCMCSamps1 = 100

# Read in the data
lv_data = np.genfromtxt("/home/ed/Dropbox/Ben and Ed/DynModel/Datasets/Lotka1b1data.csv", delimiter=",")
x1_data = np.genfromtxt("/home/ed/Dropbox/Ben and Ed/DynModel/Datasets/Lotka1x1data.csv", delimiter=",")

# Separate the data into the necessary peices.
T0 = lv_data[ : , 0 ]
P0 = lv_data[ : , [11,12,13,14,15,16,17,18,19,20,21] ]
P1 = lv_data[ : , [22,23,24,25,26,27,28,29,30,31,32] ]

# Starting values for a1, g1, d1... b1 is a function of x1
# These are at true values... at the moment
a1 = 1.0   # prey natural reproduction rate  a1 > 0
#b1 = 0.01  # rate of predation  b1 > 0
g1 = 1.0   # Preditor death rate  g1 > 0
d1 = 0.02  # death rate of preditor  d1 > 0

# initial conditions
F0 = 10               #  initial Fox population
R0 = 100              #  initial Rabbit population
y0 = [F0, R0]         #  initial condition vector

# Starting values for b0 and b1
b0 = np.log(0.01)
b1 = (np.log(0.02) - np.log(0.01))/10

# Create containers to hold the sampled values
b0hold = np.zeros( (MCMCSamps1,1) )
b1hold = np.zeros( (MCMCSamps1,1) )
a1hold = np.zeros( (MCMCSamps1,1) )
g1hold = np.zeros( (MCMCSamps1,1) )
d1hold = np.zeros( (MCMCSamps1,1) )

# Set the step values for the M-H sampler
b0step = 0.001
b1step = 0.00001
a1step = 0.00001
g1step = 0.00001
d1step = 0.00001

# Initialize the logposterior
lpost1 = lv1.LVLogPost1(P0,P1,x1_data,
                    a1,b0,b1,g1,d1,
                    y0,0,100.,101)

for iMCMC1 in range(MCMCSamps1):
    #=====================================================
    # Change b0
    b0c = b0 + np.random.normal(0,b0step,1)
    lpost1c =  lv1.LVLogPost1(P0,P1,x1_data,
                    a1,b0c,b1,g1,d1,
                    y0,0,100.,101)
    rand1u = np.random.uniform( 0.0 , 1.0 , 1 )
    diff1 = lpost1c - lpost1
    if np.log(rand1u) < diff1:
        b0 = b0c
        lpost1 = lpost1c
    #=====================================================
    # Change b1
    b1c = b1 + np.random.normal(0,b1step,1)
    lpost1c =  lv1.LVLogPost1(P0,P1,x1_data,
                    a1,b0,b1c,g1,d1,
                    y0,0,100.,101)
    rand1u = np.random.uniform( 0.0 , 1.0 , 1 )
    diff1 = lpost1c - lpost1
    if np.log(rand1u) < diff1:
        b1 = b1c
        lpost1 = lpost1c
    #=====================================================
    # Change a1
    a1c = a1 + np.random.normal(0,a1step,1)[0]
    lpost1c =  lv1.LVLogPost1(P0,P1,x1_data,
                    a1c,b0,b1,g1,d1,
                    y0,0,100.,101)
    rand1u = np.random.uniform( 0.0 , 1.0 , 1 )
    diff1 = lpost1c - lpost1
    if np.log(rand1u) < diff1:
        a1 = a1c
        lpost1 = lpost1c
    #=====================================================
    # Change g1
    g1c = g1 + np.random.normal(0,g1step,1)[0]
    lpost1c =  lv1.LVLogPost1(P0,P1,x1_data,
                         a1,b0,b1,g1c,d1,
                         y0,0,100.,101)
    rand1u = np.random.uniform( 0.0 , 1.0 , 1 )
    diff1 = lpost1c - lpost1
    if np.log(rand1u) < diff1:
        g1 = g1c
        lpost1 = lpost1c
    #=====================================================
    # Change d1
    d1c = d1 + np.random.normal(0,d1step,1)[0]
    lpost1c =  lv1.LVLogPost1(P0,P1,x1_data,
                    a1,b0,b1,g1,d1c,
                    y0,0,100.,101)
    rand1u = np.random.uniform( 0.0 , 1.0 , 1 )
    diff1 = lpost1c - lpost1
    if np.log(rand1u) < diff1:
        d1 = d1c
        lpost1 = lpost1c
    #=====================================================
    # Store the values 
    b0hold[ iMCMC1 ] = b0
    b1hold[ iMCMC1 ] = b1
    a1hold[ iMCMC1 ] = a1
    g1hold[ iMCMC1 ] = g1
    d1hold[ iMCMC1 ] = d1




#=========================================================
end_time = ts.time()
print(end_time - start_time)

