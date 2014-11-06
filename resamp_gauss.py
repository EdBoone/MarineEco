import mpl_toolkits.mplot3d

from pylab import *    

def resamp_gauss(parts,weights,cov_fac,weight_min,obs_err=False,
                 plt=False,obs=False):
    
    """A generic resample algorithm, using Gaussian clouds about particles
    
    This resampler will reproduce samples not at the exact location of existing
    samples but in a Gaussian cloud of covariance which is a fraction cov_fac
    of the original covariance (useful for deterministic models, 
    or for EnKF etc.)

    cov_fac = 0.0004 ===> 2% std dev
    cov_fac = 0.0025 ===> 5% std dev

    [state_dimension, particle_number] = shape(parts)
    [particle_number]                  = shape(weights)

    weights less than weight_min are set equal to zero, particle is eliminated
    (useful in PF update step -- use weight_min ~ 0.1/particle_number)"""

    [s_dim, particle_number] = shape(parts)
    r_samp = zeros([s_dim, particle_number])
    if plt:
        figure(1)
        subplot(projection='3d')       
        plot(parts[0,:],parts[1,:],parts[2,:],'r.')
        plot([obs[0]],[obs[1]],[obs[2]],'bo')
    if type(obs_err) != bool:
        for ii in range(particle_number):
            r_samp[:,ii] = multivariate_normal(parts[:,ii],obs_err)
    else:
        #make all weights falling below the threshold equal to zero
        save_index = weights > weight_min
        n_out = around(weights*save_index*particle_number)
        sm = sum(n_out)    
        df = (particle_number - sm)  
        sn = sign(df)
        df = int(abs(df))
        #Randomly refills the weights up to the number of the sample size
        #for use in calculating the weighted covariance for resampling
        for ii in range(df):
            idx = randint(particle_number) 
            newnum = n_out[idx] + sn
            while ((newnum < 0) or (newnum > particle_number)):
                idx = randint(particle_number)
                newnum = n_out[idx] + sn
            n_out[idx] = newnum

        cov_weights = n_out/sum(n_out) 
        #Calculate the weighted mean and weighted covariance
        mean_state = mean((parts*cov_weights),1)*particle_number
        delta_part = (mean_state - parts.T).T*sqrt(cov_weights)
        covar = cov(delta_part)*particle_number
    
        idx=0 
        # select n_out(ii) samples from a Gaussian with mean xin(ii) 
        # and covariance covar*cov_fac
        for ii in range(particle_number):
            if (n_out[ii] > 0):
                r_samp[:,idx:idx+n_out[ii]]= multivariate_normal(parts[:,ii], 
                                                                 covar*cov_fac,
                                                                 [n_out[ii]]).T
            
                idx = idx + n_out[ii]
        if plt:
            figure(2)
            subplot(projection='3d')       
            plot(r_samp[0,:],r_samp[1,:],r_samp[2,:],'y.')
            plot([obs[0]],[obs[1]],[obs[2]],'bo')
        #Jitting step
        Hsq=0.1
        Aj=sqrt(1.0-Hsq) 
        m_r_samp = mean(r_samp,1)
        mi = Aj*r_samp+(1-Aj)*tile(m_r_samp,(particle_number,1)).T
        Vt = cov(r_samp)
        for ii in range(particle_number):
            r_samp[:,ii]=multivariate_normal(mi[:,ii],Hsq*Vt)
    if plt:
        figure(3)
        subplot(projection='3d')       
        plot(r_samp[0,:],r_samp[1,:],r_samp[2,:],'m.')
        plot([obs[0]],[obs[1]],[obs[2]],'bo')
    show()
    return(r_samp)
    