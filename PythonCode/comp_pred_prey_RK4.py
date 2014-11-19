from pylab import *
from Marine_coral_module import comp_d

def CPPRK4(init_cond,end_time,time_step,A,B,C,R,K):
    
    """This is a RK4 integrator for the comp_d function
    
    This will return the trajectories of the initial conditions based on the
    RK4 integration scheme.  Inputs include the initial conditions, the time
    steps to be integrated on, and the array parameters A,B,C,R,K."""
    
    #Define the storage
    state_dim = len(init_cond)
    exp_len   = int(end_time/time_step)
    #the trajectory starts at time zero and ends at the final time
    traj      = zeros([exp_len+1, state_dim])
    traj[0,:] = init_cond
    
    for i in range(exp_len):
        k1   = comp_d(traj[i,:],i,A,B,C,R,K)
        temp = traj[i,:] + (time_step/2.0)*k1
        k2   = comp_d(temp,i,A,B,C,R,K)
        temp = traj[i,:] + (time_step/2.0)*k2
        k3   = comp_d(temp,i,A,B,C,R,K)
        temp = traj[i,:] + time_step*k3
        k4   = comp_d(temp,i,A,B,C,R,K)
        
        traj[i+1,:] = traj[i,:] + (time_step/6.0)*(k1+2.0*k2+2.0*k3+k4)
        
        #create a break in case system parameters are singular
        if (any(traj[i+1,:]) <= exp(-20)) or (any(traj[i+1,:]) >= exp(20)):
            traj[:,:] = inf
            break
    return(traj)
    