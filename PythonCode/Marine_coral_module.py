#  Code for the MarineEco project

from pylab import *
from scipy.integrate import odeint
from scipy.special import gammaln


###############################################################################
# Define the time derivative for the competition model

def comp_d(initial_cond,t,A,B,C,R,K):
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
     
    r1  - recruitment rate of coral
    r2  - recruitment rate of macro algae
    r3  - recruitment rate of turf

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
         
    [r1,r2,r3]                = R

    [k1,k2,k3]                = K
     
    #Calculate the derivatives
    dx1 = r1*x1 - r1*x1*(x1+a12*x2+a13*x3)/k1
    
    dx2 = r2*x2 - r2*x2*(x2+a21*x1+a23*x3)/k2 - b24*x2*x4
    
    dx3 = r3*x3 - r3*x3*(x3+a31*x1+a32*x2)/k3 - b35*x3*x5
    
    dx4 = b24*x2*x4 - c4*x4
    
    dx5 = b35*x3*x5 - c5*x5
    
    return array([dx1,dx2,dx3,dx4,dx5])
        
###############################################################################
#Here we define the log scale log-normal 
#likelyhood function for the parameter estimates

def ln_ls(mean, sd, outcomes):
    
    """This function evaluates the likelyhood of a given sample
    
    This is a vectorized parameter likelyhood function which takes an array of
    data (mean), the associated array of standard deviations (sd1), 
    and the parameter dependent outcomes (sample)
    
    Inputs are of the form
    
    [number of obs, number of variables]  = shape(mean)
    
    [number of variables]                 = shape(sd) 

    [number of obs, number of variables]  = shape(outcomes)

    Output is an array on lognormal scale of the likelyhood of each outcome
    given the associated mean and (fixed) standard deviation of the data."""
    
    #define the log of the denominator
    # with multiplication across arrays with Schur product
    denom  = -log(outcomes*sd*sqrt(2.0*pi))
    #define the numerator in log scale
    #again operating on each array element individually
    numer  = -0.5*((log(outcomes) - mean)/sd)**2
    
    #calclulate the likelyhood
    output = numer + denom
    output = sum(output)
    
    return(output)


###############################################################################
# Create the chi square likelyhood on the Log scale    
def chisq_ls(outcome,K):

    """Log scale Chi Square likelyhood distribution given integer std normals
    
    This returns the log scale likelyhood of the outcome, given that it is
    the sum of integer independent standard normal random variables.  This is
    vectorized to take an array of independent outcomes and associated array
    K of integers for the number of standard normals associated to the 
    array element of the outcomes.  It returns the likelyhood of all
    independent outcomes as their product (in log scale)."""
    
    #log of denominator
    denom  = -K/2.0*log(2.0) - gammaln(K/2.0)
    #log of numerator
    numer  =  (K/2.0 - 1.0)*log(outcome) - outcome/2.0
    #calculate the likelyhood of entire array in log scale
    output = numer + denom
    output = sum(output)
    
    return output

###############################################################################
# Create the log-scale likelihood for the marine eco data
#likelyhood distribution is by default log-normal 

def coral_lk(obs, var, par, init_cond, end_t, incr):
    
    """This solves the system first with given params and measures likelyhood
    
    Observations are assumed to all be of the full state dimension and given
    yearly wihtout gaps.
    
    The inputs are given as 
    
    obs          = array containing all observations of state
    var          = array of (fixed) variances of the observations
    par          = current guess at parameters for the likelyhood funciton
    initial_cond = the initial conditions for the state variables
    start_t      = begining time of the experiment
    end_t        = end time of the experiemnt
    incr         = time step for the ODE solver
    likely       = likelyhood funciton, defaul log scale log-normal
    
    The function returns the likelyhood function for the model run under
    the initial conditions and parameters given the observational data."""
    
    #unpack the parameters and put them in the input form for the derivative
    A = par[0:6]
    B = par[6:8]
    C = par[8:10]
    R = par[10:13]
    K = par[13:]
    
    #define the array of time steps to integrate on
    time = linspace(0,end_t, end_t/incr +1)

    #trajectory is integrated
    traj = odeint(comp_d, init_cond, time, args = (A,B,C,R,K))
    print(traj[traj<0])
    #extract the points to compare the trajectory to the observations
    outcomes = traj[::1/incr,:]  

    #trajectory at obsservation times is compared with the observations
    L_hood = ln_ls(obs,var,outcomes)
    
    return L_hood

    
###############################################################################
# Create the log-posterior (up to proportionality constant)
#for the marine eco data

def para_post(init_cond, par, obs, obs_sd, obs_df, 
              par_mean, par_sd, end_t, incr):
    
    """This calculates the probability of the given parameters in the posterior
    
    Inputs are the initial conditions for the model, including all state
    variables, the parameters to measure, the observational data to measure
    against, the standard deviations of the observational data, degrees of
    freedom for the observational variances,
    the means of the parameters for the prior distribution, the variance 
    of these parameters in the prior, the end time of the experiment 
    and the integration interval.
    
    Note: this is written in log scale, so the sum is taken over the likelyhood
    and the prior."""
    
    #define the parameters for the varience likelyhood function
    par_l = par[:13]
    
    #find the likelyhood of the parameters given the observed data
    post = coral_lk(obs, obs_sd, par, init_cond,end_t,incr) 
    #multiply (in log scale) the likelyhood against the prior for the variances
    post = post + chisq_ls(obs_sd**2,obs_df)
    #multiple (in log scale) with the prior probability for the parameters
    post = post + ln_ls(par_mean,par_sd,par_l)

    return post

###############################################################################