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
    
    #define the denominator, multiplication across arrays with Schur product
    denom  = sample*variance*sqrt(2*pi)
    #define the numerator, again operating on each array element individually
    numer  = exp((-(log(sample) - mean)**2)/(2*variance**2))
    #the log scale likelyhood is computed for each array element
    output = log(numer/denom)
    
    return(output)
    
###############################################################################