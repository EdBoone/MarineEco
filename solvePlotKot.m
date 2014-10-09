% script for finding and plotting the solution to the ode system
clear all ;

% [t,y] = ode45('kot1',[0,15],[1.6e-7;8.2e-7;9.8e-5;2.3e-5]); % dimensional equations

[t,y] = ode45('kotDoubleForced2',[0,100],[0.25;0.2;0.1]); % dimensionless equations

% t = 0:0.01:15 ;
% y = ode4('kotDoubleForced2', t, [0.25;0.5;0.5]);

% The above ode23 command can be changed to any of the ode solvers
% in the MATLAB suite:  ode45, ode113, ode23s, etc.  Just type
% "help ode45" to get a list of the other temporal integration 
% routines available.


% The first column of y should be the solution for S (or x),
% the second column for H (or y), the third for P (or z).

figure(1)
% plot(t,y(:,1),'r',t,y(:,2),'b',t,y(:,3),'g')
% plot(t,y(:,2),'r');
subplot(3,1,1); plot(t,y(:,1),'r');
subplot(3,1,2); plot(t,y(:,2),'b');
subplot(3,1,3); plot(t,y(:,3),'g');

figure(2)
plot(y(:,1),y(:,2),'r',y(:,1),y(:,3),'b');
% plot(y(:,1),y(:,2),'b');

figure(3)
plot(y(:,2),y(:,3),'b');

figure(4)
plot3(y(:,1), y(:,2), y(:,3),'b',0.25,0.2,0.1,'g*');
grid on