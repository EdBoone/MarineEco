from pylab import *
from scipy.optimize import minimize
from comp_pred_prey_RK4 import CPPRK4
from Marine_coral_module import comp_d
import matplotlib

###############################################################################
#This module will attempt to fit the predator/prey/competition parameters
#to the data from '97 to '07 where these dynamics dominate.
#Initial conditions for coral cover is drawn from the data but the
#rest of the initial conditions are guesses based upon the data in 2006/2007
###############################################################################
#Import the data sets and process for the model
x1 = genfromtxt('coral_data_long')
x2 = genfromtxt('macro_data_long')
x3 = genfromtxt('turf_data_long')
x4 = []
x5 = []

a12 = 1.0
a13 = 1.0
a21 = 1.0
a23 = 1.0
a31 = 0.5
a32 = 0.5
b24 = .5
b35 = .0001
c4 = 0.5
c5 = 0.01
K1 = 65.0
K2 = 100.0
K3 = 60.0
r1 = 130.07
r2 = 80.76
r3 = 156.0

f2 = open('test_herb_data.csv')
for line in f2:
    a = line.strip().split(',')
    if a[1] == 'Browser':
        x4.append(float(a[3]))
    elif a[1] == 'Grazer/detritivore':
        x5.append(float(a[3]))

x1_init = x1[5,1]
x2_init = mean(x2[:3,1])
x3_init = mean(x3[:3,1])
x4_init = mean(x4[:3])
x5_init = mean(x5[:3])

c_obs = x1[5:16,1]
m_obs = x2[:3,1]
t_obs = x3[:3,1]
b_obs = x4[:3]
g_obs = x5[:3]

X = array([x1_init,x2_init,x3_init,x4_init,x5_init])
R = array([r1,r2,r3])
K = array([K1,K2,K3])

#inputs defined for the minimization routine
par = array([a12, a13, a21, a23, a31, a32, b24, b35, c4, c5])
init_conds = (X,R,K,c_obs,m_obs,t_obs,b_obs,g_obs)

###############################################################################
#Here we will define a funciton which will return the misfit from the data
#for the least squares estimate

def misfit(parameters,X,R,K,c_obs,m_obs,t_obs,b_obs,g_obs):
    
    end_time  = 10.0
    time_step = 0.01
    time = linspace(0,end_time,int(end_time/time_step)+1)
    A = parameters[:6]
    B = parameters[6:8]
    C = parameters[8:]

    traj = CPPRK4(X,end_time,time_step,A,B,C,R,K)

    plot(time,traj)
    show()
    
    diff = sum(abs(traj[::1.0/time_step,0] - c_obs))
    diff = diff + sum(abs(m_obs - traj[800::1.0/time_step,1]))
    diff = diff + sum(abs(t_obs - traj[800::1.0/time_step,2]))
    diff = diff + sum(abs(b_obs - traj[800::1.0/time_step,3]))
    diff = diff + sum(abs(g_obs - traj[800::1.0/time_step,4]))
    
    return diff

###############################################################################
#optimization is run over the parameters
opt = minimize(misfit,par,args = init_conds)

print(opt)