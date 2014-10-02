function xp = kotDoubleForced2(t,x)

xp = zeros(3,1);

% These are the dimensionless equations from the Kot paper

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
% constants in dimensionless equations
%
A = mu1/D ;
a = k1/si ;
B = mu2/D ;
b = k2/y1/si ;
epsilon = 0.6 ;
T = 15 ;
% omega = 2*pi/D/T ;
% omega = 4*pi ;
omega = 5*pi/6 ;

%
%
% Differential equations

xp(1,1) = 1 + epsilon*sin(omega*t) - x(1,1) - A*x(1,1)*x(2,1)/( a+x(1,1) ) ;

xp(2,1) = A*x(1,1)*x(2,1)/(a+x(1,1)) - x(2,1) - B*x(2,1)*x(3,1)/( b+x(2,1) );

xp(3,1) = B*x(2,1)*x(3,1)/(b+x(2,1)) - x(3,1) ;


