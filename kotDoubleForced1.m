function xp = kotDoubleForced1(t,x)

xp = zeros(3,1);

% These are the dimensional equations from the Kot paper

% dilution rate (1/hour)

D = 0.1 ;
%
% units of Si are mg/l
si = 115 ;
%
% maximum specific growth rate of prey and predator (1/hour)
mu1 = 0.5 ; % prey
mu2 = 0.2; % predator
%
% yield of prey per unit mass of substrate (dimensionless)
y1 = 0.4 ;
%
% biomass yield of predator per unit mass of prey (dimensionless)
%
y2 = 0.6 ;
% half-saturation (Michaelis-Menten) constants for prey and predator (mg/l)
k1 = 8; % prey
k2 = 9; % predator
% 
% other variables
epsilon = 0.6 ;
T = 15 ;
omega = 2*pi/T ;
% omega = 5*pi/6 ;

%
%
% Differential equations

xp(1,1) = D*( Si*(1+epsilon*sin(omega*t)) - x(1,1)) - mu1*x(1,1)*x(2,1)/y1/(k1+x(1,1)) ;

xp(2,1) = mu1*x(1,1)*x(2,1)/(k1+x(1,1)) - D*x(2,1) - mu2*x(2,1)*x(3,1)/y2/(k2+x(2,1)) ;

xp(3,1) = mu2*x(2,1)*x(3,1)/(k2+x(2,1)) - D*x(3,1) ;



