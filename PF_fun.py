#### PF-AUS Filter with Gassian resampling for the Lorenz 3 variable system ###
from pylab import *
from Kot_int import KotDF1
from EL63 import EL63V
from resamp_gauss import resamp_gauss
import mpl_toolkits.mplot3d

def PF(truth, obs, R,initial_cloud,h):
 
    """Returns [part, mean_state, mean_error,avg_mean_error] for
    mandatory inputs truth trajectory, observations, and intitial cloud.
    
    [state_dimension, tfin+1]          = shape(truth)    
    [state_dimension, Nanl]            = shape(obs)
    [state_dimension, state_dimension] = shape(R)
    [state_dimension, particle_number] = shape(initial_cloud)
    [state_dimension, Nanl]            = shape(bred_perts)
    [state_dimension, time_step]       = shape(mean_state)
    [state_dimension, time_step]       = shape(mean_error)
    [time_step]                        = shape(avg_mean_error)
    # of times resampler triggered     = rsn
    
    Optional arguments include sigma=a=10, beta=8/3, rho=29, tanl=20.
    The duration of the experiment tfin=2000, and analysis  
    time particle cloud is propagated between analysis times,
    at which point the observation is incorporated through the Bayesian 
    update of the weights.  Resampling is initiated according to
    resamp_gauss when Neff is less than the threshold."""
    
    #Recover the integration length
    tfin = int(len(truth[0,:])-1)
    # NUMBER OF ANALYSES
    Nanl = len(obs[0,:])
    #Recover the time between analyses
    tanl = int(tfin/Nanl)
    
    #Set initial trajectory matricies and weights
    [state_dimension, particle_number] = shape(initial_cloud)    
    part = zeros([state_dimension,tfin+1,particle_number])
    part[:,0,:] = initial_cloud
    part_weights = ones(particle_number)*(1.0/particle_number)
    Q = inv(R)
    mean_state = zeros([state_dimension,tfin+1])
    """
    figure()
    subplot(projection='3d')       
    plot(part[0,0,:],
         part[1,0,:],
         part[2,0,:],'y.')
    show()    
    """
    #corrections = zeros([3,Nanl])

    #loop over number of analyses
    for j in range(Nanl):

        #propagate each particle to next analysis step
        part[:,j*tanl:(j+1)*tanl+1,:] = EL63V(tanl,h,part[:,j*tanl,:],a,r,b)
        """
        for jj in range(tanl):
            figure()
            subplot(projection='3d')       
            plot(part[0,j*tanl + jj,:],
            part[1,j*tanl + jj,:],
            part[2,j*tanl + jj,:],'r.')
            plot([truth[0,j*tanl+jj]],
                 [truth[1,j*tanl+jj]],
                 [truth[2,j*tanl+jj]],'bo')
            show()
        print('analysis time' + str(j))
        """
        #Calculate Mean at each propagated state up to analysis time step
        weighted_cloud = part[:,j*tanl:(j+1)*tanl,:]*part_weights
        mean_state[:,j*tanl:(j+1)*tanl]=mean(weighted_cloud,2)*particle_number

        #Bayesian update is applied to particle wieghts    
        innov = obs[:,j] - part[:,(j+1)*tanl,:].T
        temp = array([sum(innov*Q[:,0],1),
                      sum(innov*Q[:,1],1),
                      sum(innov*Q[:,2],1)])
        w_temp =  part_weights*exp(-0.5*sum(temp*innov.T,0))**(1.0/3.0)
        if sum(w_temp)<exp(-9):
            w_temp = ones([particle_number])*(1.0/particle_number)
        else:
             w_temp = w_temp/sum(w_temp)
        part_weights = part_weights*w_temp
        part_weights = part_weights/sum(part_weights)
        
        #effective number of samples computed
        n_eff = 1/(part_weights.dot(part_weights))
        
        #Resampling step if number of effective samples falls below threshold 
        if (n_eff < particle_number/2):
            part[:,(j+1)*tanl,:] = resamp_gauss(part[:,(j+1)*tanl,:],
                                                    part_weights,.0004,
                                                    .1/particle_number)
                        
                                                    
            part_weights = ones(particle_number)*(1.0/particle_number)                       
        #Find mean state updated weights      
        weighted_cloud = part[:,(j+1)*tanl,:]*part_weights             
        mean_state[:,(j+1)*tanl] = mean(weighted_cloud,1)*particle_number
        """
        figure()
        subplot(projection='3d')       
        plot(part[0,(j+1)*tanl,:],
             part[1,(j+1)*tanl,:],
             part[2,(j+1)*tanl,:],'y.')
        plot([obs[0,j]],[obs[1,j]],[obs[2,j]],'bo')
        show()
        """
    #mean_error_av = mean_error*mean_error
    mean_error = abs(mean_state - truth)
    error_dist = sqrt(sum(mean_error*mean_error,0))    
    r_avg_mean_error = zeros(tfin+1)
    print('PF')
    for i in range(tfin+1):
        r_avg_mean_error[i] = mean(error_dist[0:i+1])
    return(part, mean_state, mean_error, r_avg_mean_error)