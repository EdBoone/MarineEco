function xp = kot(t,x)

xp = zeros(3,1);

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
%
% Differential equations

xp(1,1) = x(1,1)*(munr*x(4,1)/(x(4,1)+knr) - dr) - mupr/ypr*x(3,1)*x(1,1)/(x(1,1)+kpr) - D*x(1,1) ;

xp(2,1) = x(2,1)*(munc*x(4,1)/(x(4,1)+knc) - dc) - mupc/ypc*x(3,1)*x(2,1)/(x(2,1)+kpc) - D*x(2,1) ;

xp(3,1) = x(3,1)*(mupr*x(1,1)/(x(1,1)+kpr) + mupc*x(2,1)/(x(2,1)+kpc) - dp) - D*x(3,1) ;

xp(4,1) = D*N0 - x(1,1)*(munr/ynr*x(4,1)/(x(4,1)+knr)) - x(2,1)*munc/ync*x(4,1)/(x(4,1)+knc) - D*x(4,1);
