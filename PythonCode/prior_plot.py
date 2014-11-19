from pylab import *
from Marine_coral_module import CPPRK4
import matplotlib

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

init_cond = array([x1_init,x2_init,x3_init,x4_init,x5_init])
end_t = 10
incr = .01
time = linspace(0,end_t,end_t/incr + 1)

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
             
A = [a12, a13, a21, a23, a31, a32]
B = [b24, b35]
C = [c4, c5]
R = [r1, r2, r3]
K = [K1, K2, K3]

traj = CPPRK4(init_cond,end_t,incr,A,B,C,R,K)

figure()
plot(time,traj[:,0],'g-',label='Coral')
plot(linspace(0,10,11),c_obs,'yo',label='Observations')
title('Coral and Observations - Prior Mean')
xlabel('Year')
ylabel('% Coverage')
legend(bbox_to_anchor=(1.05, 1),loc=2, borderaxespad=0.)
show()

figure()
plot(time,traj[:,1],'g-',label='Macro Algae')
plot(linspace(8,10,3),m_obs,'yo',label='Observations')
title('Macro Algae and Observations - Prior Mean')
xlabel('Year')
ylabel('% Coverage')
legend(bbox_to_anchor=(1.05, 1),loc=2, borderaxespad=0.)
show()

figure()
plot(time,traj[:,2],'g-',label='Turf Algae')
plot(linspace(8,10,3),t_obs,'yo',label='Observations')
title('Turf Algae and Observations - Prior Mean')
xlabel('Year')
ylabel('% Coverage')
legend(bbox_to_anchor=(1.05, 1),loc=2, borderaxespad=0.)
show()

figure()
plot(time,traj[:,3],'g-',label='Browser')
plot(linspace(8,10,3),b_obs,'yo',label='Observations')
title('Browsers and Observations - Prior Mean')
xlabel('Year')
ylabel('Biomass')
legend(bbox_to_anchor=(1.05, 1),loc=2, borderaxespad=0.)
show()

figure()
plot(time,traj[:,4],'g-',label='Grazer/Detritivore')
plot(linspace(8,10,3),g_obs,'yo',label='Observations')
title('Grazer/Detrivivores and Observations - Prior Mean')
xlabel('Year')
ylabel('Biomass')
legend(bbox_to_anchor=(1.05, 1),loc=2, borderaxespad=0.)
show()

figure()
plot(time,traj[:,0],'r-',label='Coral')
plot(time,traj[:,1],'b-',label='Macro')
plot(time,traj[:,2],'g-',label='Turf')
plot(time,traj[:,3],'m-',label='Browser')
plot(time,traj[:,4],'c-',label='Graser/ Detritivore')
title('All - Prior')
xlabel('Year')
ylabel('Species')
legend(bbox_to_anchor=(1.05, 1),loc=2, borderaxespad=0.)
show()