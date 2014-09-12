%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  Basic Lotka-Volterra stuff
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Define initial conditions.
t0 = 0;
tfinal = 15;
y0 = [20 20]';
% Simulate the differential equation.
tfinal = tfinal*(1+eps);
[t,y] = ode23('lotka',[t0 tfinal],y0);

% Time plot
subplot(1,2,1)
plot(t,y)
title('Time history')

% Phase plane plot
subplot(1,2,2)
plot(y(:,1),y(:,2))
title('Phase plane plot')

% Use a different solver
[T,Y] = ode45('lotka',[t0 tfinal],y0);

% Phase plane plot
subplot(1,1,1)
title('Phase plane plot')
plot(y(:,1),y(:,2),'-',Y(:,1),Y(:,2),'-');
legend('ode23','ode45')