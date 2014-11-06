from pylab import *

###############################################################################
#Here we define the log scale likelyhood functions for the parameter estimates

###############################################################################

def LN_para_likelyhood(mean,variance,sample):
    
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