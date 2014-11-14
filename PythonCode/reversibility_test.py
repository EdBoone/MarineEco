from pylab import *
from scipy.integrate import odeint
from Marine_coral_module import comp_d
import matplotlib

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
r1 = 2.0
r2 = 2.0
r3 = 2.0
K1 = 65.0
K2 = 100.0
K3 = 60.0

A = array([a12,a13,a21,a23,a31,a32])
B = array([b24,b35])
C = array([c4,c5])
R = array([r1,r2,r3])
K = array([K1,K2,K3])
             
init_cond = array([37.57998016,3.55818651,3.15670238,21.66666667,189.1666667])
time = linspace(0,30,3001)

#trajectory is integrated
traj = odeint(comp_d, init_cond, time, args = (A,B,C,R,K))
b_traj = odeint(comp_d,traj[-1,:],time[::-1],args = (A,B,C,R,K))
diff = abs(traj - b_traj[::-1,:])

figure()
plot(time,diff[:,0], label='Coral DIFF')
plot(time,diff[:,1], label='Macro DIFF')
plot(time,diff[:,2], label='Turf Diff')
plot(time,diff[:,3], label='Brow Diff')
plot(time,diff[:,4], label='Detr Diff')
legend()
show()