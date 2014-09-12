function [y1, x1] = LV1FixedPoint(a1,b1,e1,c1)
        % a1 = Natural Growth Rate of R in absence of F
        % b1 = Death rate per encounter or R due to predation
        % g1 = efficiency fo turning R into F
        % d1 = Natural death rate of F in absence of R.
        y1 = a1/b1;      % R
        x1 = c1/(e1*b1);   % F
