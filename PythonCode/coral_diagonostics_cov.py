from pylab import *
from Marine_coral_module import para_post_alt
import matplotlib

###############################################################################
#This module will calculate the empirical posterior density
#with the Metropolis-Hasting algorithm for the coral
#predator-prey-competition-model

###############################################################################
#Define the storage for the observations of state varaibles to be read from
#the files

#Import the data sets and process for the model
x1 = genfromtxt('coral_data_long')
x2 = genfromtxt('macro_data_long')
x3 = genfromtxt('turf_data_long')
x4=[]
x5=[]

f2 = open('test_herb_data.csv')
for line in f2:
    a = line.strip().split(',')
    if a[1] == 'Browser':
        x4.append(float(a[3]))
    elif a[1] == 'Grazer/detritivore':
        x5.append(float(a[3]))

f2.close()

x1_init = x1[5,1]
x2_init = mean(x2[:,1])
x3_init = mean(x3[:3,1])
x4_init = mean(x4[:3])
x5_init = mean(x5[:3])

c_obs = x1[5:16,1]
m_obs = x2[:3,1]
t_obs = x3[:3,1]
b_obs = x4[:3]
g_obs = x5[:3]

#create observational array, each row representing a year, each column
#a specific state variable
#NOTE: this matches the experiment and the observations with
# the starting year to be 2006
obs_1 = c_obs
obs_2 = array([m_obs,t_obs,b_obs,g_obs]).T
init_cond = array([x1_init,x2_init,x3_init,x4_init,x5_init])
end_t = 10
incr = .01
###############################################################################
# Set the prior parameters

# Prior parameter means
a12m = 2.21859143
a13m = 1.82923341
a21m = 2.63540104e-01
a23m = 1.54167405
a31m = 2.68116615e-02
a32m = 5.88265886e-01
b24m = 2.28485761
b35m = 1.32683022
c4m  = 1.81193764
c5m  = 2.23199371
r1m  = 1.27508204e+02
r2m  = 8.22849835e+01
r3m  = 1.58610308e+02
par_mean = [a12m, a13m, a21m, a23m, a31m, a32m, b24m, b35m,
            c4m, c5m,r1m, r2m, r3m]

# Prior parameter sandard deviations
a12sd = 0.1
a13sd = 0.1
a21sd = 0.1
a23sd = 0.1
a31sd = 0.1
a32sd = 0.1
b24sd = 0.01
b35sd = 0.01
c4sd  = 0.01
c5sd  = 0.01
r1sd  = 1.0
r2sd  = 1.0
r3sd  = 1.0
par_sd = [a12sd, a13sd, a21sd, a23sd, a31sd, a32sd,
          b24sd, b35sd, c4sd, c5sd, r1sd, r2sd, r3sd]

# Degrees of freedom parameters for the variances.
df1 = 1.0
df2 = 1.0
df3 = 1.0
df4 = 1.0
df5 = 1.0
obs_df_1 = df1
obs_df_2 = array([df2, df3, df4, df5])

###############################################################################
#initial parameter vals for the Metropolis-Hasting algorithm

# Set the initial values for the parameters
a12 = 2.21859143
a13 = 1.82923341
a21 = 2.63540104e-01
a23 = 1.54167405
a31 = 2.68116615e-02
a32 = 5.88265886e-01
b24 = 2.28485761
b35 = 1.32683022
c4  = 1.81193764
c5  = 2.23199371
r1  = 1.27508204e+02
r2  = 8.22849835e+01
r3  = 1.58610308e+02
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
obs_sd_1 = sig1
obs_sd_2 = array([sig2, sig3, sig4, sig5])

###############################################################################
# Set the step values for the M-H algorithm
a12s = 0.01
a13s = 0.01
a21s = 0.01
a23s = 0.01
a31s = 0.01
a32s = 0.01
b24s = 0.01
b35s = 0.01
c4s = 0.01
c5s = 0.01
r1s = 1.0
r2s = 1.0
r3s = 1.0

I = identity(13)
step = I*array([a12s,a13s,a21s,a23s,a31s,a32s,b24s,b35s,c4s,c5s,r1s,r2s,r3s])

###############################################################################
#Define the size of the ensemble and run the algorithm

#number of samples to be accepted, varying per parameter (total samples
#is equal to n_samples*number of parameters)
n_samples = 10000
par_temp = zeros(len(par))
prob_0 = para_post_alt(init_cond,par,obs_1,obs_2,obs_sd_1,obs_sd_2,obs_df_1,
                       obs_df_2,par_mean,par_sd,end_t,incr,pt=True)


#storage for the samples and the weights
f_post   = open('posterior_points','w+')
f_weight = open('posterior_weights','w+')

#loop over M-H algorithm, varying over all parameters:
j = 0
#initialize the search with the first parameters
par_temp[:] = par[:]
#initialize the initial probability with that of the first parameters
prob_temp = float(prob_0)
trial = zeros(len(par))
#initialize the parameter value to be varied
while j < n_samples:
    #perturb the parameter values by a random normal
    #but make sure the parameters remain positive
    trial[:] = par_temp[:]    
    trial[:13] = par_temp[:13] + multivariate_normal(zeros(13),step)
    while any(trial) <= 0:
        trial[:13] = par_temp[:13]+multivariate_normal(zeros(len(p_pert)),step)  

    #determine the posterior probability of this value
    prob = para_post_alt(init_cond,trial,obs_1,obs_2,obs_sd_1,obs_sd_2,
                         obs_df_1, obs_df_2, par_mean,par_sd,end_t,incr)
    #generate uniform random for acceptance/ rejection criteria
    u = log(uniform())
    if (prob - prob_temp) > u:
        #if the ratios of probs of the new and old parameters is greater
        #than the acceptance criteria, store the parameters and the
        #weight according to the posterior distribution
        f_post.write(str(trial)+'\n')
        f_weight.write(str(prob)+'\n')
            
        #the new parameter value  and probability is stored for the
        #next iteration of the perturbation loop
        par_temp[:] = trial[:]
        prob_temp = float(prob)
            
        #sample count for the parameter goes up by 1
        j = j+1
        print(j)
        
f_post.close()
f_weight.close()
