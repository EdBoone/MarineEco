function lotkavolterra(r0,f0,T)
% Solve the Lotka-Volterra equations for a population
% of rabbits r(t) and foxes f(t) on time interval [0,T]
% and plot the population dynamics
%
% r0 - Initial rabbit population
% f0 - Initial fox population
% T - Final time
%  Try r0 = 1500
%      f0 = 30
%       T = 200
[t,y]=ode45(@lvsystem,[0 T],[r0 f0]);

subplot(1,2,1)
plot(t,y(:,1),'b-',t,y(:,2),'r-','LineWidth',2);
legend('Rabbits','Foxes');
xlabel('time','FontSize',16);
ylabel('population','FontSize',16);
subplot(1,2,2)

plot(y(:,1),y(:,2),'k-','LineWidth',2)
xlabel('Rabbits','FontSize',16);
ylabel('Foxes','FontSize',16);

function dydt=lvsystem(t,y)
dydt=[0.05*y(1)*(1-0.01*y(2)); 0.1*y(2)*(0.005*y(1)-2)];