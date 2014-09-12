function [dF, dR] = LV1RkFn(Fi,Ri,a1,b1,e1,c1)
        % a1 = Natural Growth Rate of R in absence of F
        % b1 = Death rate per encounter or R due to predation
        % g1 = efficiency fo turning R into F
        % d1 = Natural death rate of F in absence of R.
        dR = a1*Ri - b1*Ri*Fi;      % R
        dF = e1*b1*Ri*Fi - c1*Fi;   % F
