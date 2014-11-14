from pylab import *
from Marine_coral_module import para_post
import matplotlib

###############################################################################
#This module will calculate the empirical posterior density
#with the Metropolis-Hasting algorithm for the coral
#predator-prey-competition-model

###############################################################################
#Define the storage for the observations of state varaibles to be read from
#the files

x1=[]
x2=[]
x3=[]
x4=[]
x5=[]

#Read data and append squentially in time into the associated state varaible
f1 = open('test_c_a_data.csv')
for line in f1:
    a = line.strip().split(',')
    if a[3] == 'Coral':
        x1.append(float(a[2]))
    elif a[3] == 'Turf':
        x2.append(float(a[2]))
    elif a[3] == 'Macroalgae':
        x3.append(float(a[2]))
    
f2 = open('test_herb_data.csv')
for line in f2:
    a = line.strip().split(',')
    if a[1] == 'Browser':
        x4.append(float(a[3]))
    elif a[1] == 'Grazer/detritivore':
        x5.append(float(a[3]))        

#create observational array, each row representing a year, each column
#a specific state variable
#NOTE: this matches the experiment and the observations with
# the starting year to be 2006
obs = array([x1[1:],x2[1:],x3[1:],x4,x5]).T
init_cond = obs[0,:]
end_t = len(obs[:,0]) - 1
incr = .01
###############################################################################
# Set the prior parameters

# Prior parameter means
a12m = 0.1
a13m = 0.1
a21m = 0.1
a23m = 0.1
a31m = 0.1
a32m = 0.1
b24m = 0.0001
b35m = 0.0001
c4m  = 0.001
c5m  = 0.001
r1m  = 1.0
r2m  = 1.0
r3m  = 1.0
par_mean = [a12m, a13m, a21m, a23m, a31m, a32m, b24m, b35m,
            c4m, c5m,r1m, r2m, r3m]

# Prior parameter sandard deviations
a12sd = 0.001
a13sd = 0.001
a21sd = 0.001
a23sd = 0.001
a31sd = 0.001
a32sd = 0.001
b24sd = 0.000001
b35sd = 0.000001
c4sd  = 0.00001
c5sd  = 0.00001
r1sd  = 0.01
r2sd  = 0.01
r3sd  = 0.01
par_sd = [a12sd, a13sd, a21sd, a23sd, a31sd, a32sd,
          b24sd, b35sd, c4sd, c5sd, r1sd, r2sd, r3sd]

# Degrees of freedom parameters for the variances.
df1 = 1.0
df2 = 1.0
df3 = 1.0
df4 = 1.0
df5 = 1.0
obs_df = array([df1, df2, df3, df4, df5])

###############################################################################
#initial parameter vals for the Metropolis-Hasting algorithm

# Set the initial values for the parameters
a12 = 0.1
a13 = 0.1
a21 = 0.1
a23 = 0.1
a31 = 0.1
a32 = 0.1
b24 = 0.0001
b35 = 0.0001
c4 = 0.001
c5 = 0.001
r1 = 1.0
r2 = 1.0
r3 = 1.0
K1 = 65.0
K2 = 100.0
K3 = 60.0
par = array([a12, a13, a21, a23, a31, a32, b24, b35, 
             c4, c5, r1, r2, r3, K1, K2, K3])

# Standard deviation
sig1 = 1.0
sig2 = 1.0
sig3 = 1.0
sig4 = 1.0
sig5 = 1.0
obs_sd = array([sig1, sig2, sig3, sig4, sig5])

###############################################################################
# Set the step values for the M-H algorithm
a12s = 0.001
a13s = 0.001
a21s = 0.001
a23s = 0.001
a31s = 0.001
a32s = 0.001
b24s = 0.00001
b35s = 0.00001
c4s = 0.001
c5s = 0.001
r1s = 0.001
r2s = 0.001
r3s = 0.001

step = array([a12s,a13s,a21s,a23s,a31s,a32s,b24s,b35s,c4s,c5s,r1s,r2s,r3s])

###############################################################################
#Define the size of the ensemble and run the algorithm

#number of samples to be accepted, varying per parameter (total samples
#is equal to n_samples*number of parameters)
n_samples = 2
par_temp = zeros(13)
prob_0 = para_post(init_cond,par,obs,obs_sd,obs_df,par_mean,par_sd,end_t,incr)


#storage for the samples and the weights
post = zeros([13*n_samples, 13])
weights = zeros([13*n_samples])

#loop over M-H algorithm, varying one parameter at a time
for i in range(13):
    #initialize the number of samples taken for this parameter
    j = 0
    #initialize the search with the first parameters
    par_temp[:] = par[:13]
    #initialize the initial probability with that of the first parameters
    prob_temp = float(prob_0)
    #initialize the parameter value to be varied
    p_hold = par[i]
    while j < n_samples:
        #perturb the currently varying parameter value by a random normal
        #but make sure the parameters remain positive
        temp = p_hold + normal(0.0,step[i])
        while temp <= 0:
            temp = p_hold + normal(0.0,step[i])
        par_temp[i] = temp  
        #determine the posterior probability of this value
        prob = para_post(init_cond,par_temp,obs,obs_sd,obs_df,
                         par_mean,par_sd,end_t,incr)
        #generate uniform random for acceptance/ rejection criteria
        u = log(uniform())
        if (prob - prob_temp) > u:
            #if the ratios of probs of the new and old parameters is greater
            #than the acceptance criteria, store the parameters and the
            #weight according to the posterior distribution
            post[j+n_samples*i,:] = par_temp
            weights[j+n_samples*i] = prob
            
            #the new parameter value  and probability is stored for the
            #next iteration of the perturbation loop
            p_hold = par_temp[i]
            prob_temp = float(prob)
            
            #sample count for the parameter goes up by 1
            j = j+1
        
print(post)
print(weights)
