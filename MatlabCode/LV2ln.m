%=======================================================
% Solve the Lotka-Volterra ODE using Runga-Kutta 4th order method
% It seems to be stable
function [t1, Fi, Ri, loglike1] = LV2ln( F0, R0, par, initial_cond, start_t1, end_t1 )
     %-time-grid-----------------------------------
     t1  = start_t1:end_t1;
     n1 = size(t1,1);
     Fi = zeros(n1,1);
     Ri = zeros(n1,1);
     Fi(1) = initial_cond(1);
     Ri(1) = initial_cond(2);
     a1 = exp(par(1));
     b1 = exp(par(2));
     e1 = exp(par(3));
     c1 = exp(par(4));
     h1 = 1;
     loglk1 = 0;
     %differential-eq-system----------------------
    for i = start_t1:(end_t1-1)
        % the model equations (see Munz et al. 2009)
        %dx = x*(alpha - beta*y)
        %dy  = -y*(gamma - delta*x)
        %f0 = Fi*(a1 - b1*Ri)
        %f1 = -Ri*(g1 - d1*Fi)
        [FiK1, RiK1] = LV1RkFn(Fi(i), Ri(i), a1,b1,e1,c1);
        [FiK2, RiK2] = LV1RkFn(Fi(i) + 0.5*h1*FiK1, Ri(i) + 0.5*h1*RiK1, a1,b1,e1,c1);
        [FiK3, RiK3] = LV1RkFn(Fi(i) + 0.5*h1*FiK2, Ri(i) + 0.5*h1*RiK2, a1,b1,e1,c1);
        [FiK4, RiK4] = LV1RkFn(Fi(i) + h1*FiK3, Ri(i) + h1*RiK3, a1,b1,e1,c1);
        Fi(i+1) = Fi(i) + h1/6*(FiK1 + 2*FiK2 + 2*FiK3 + FiK4);
        Ri(i+1) = Ri(i) + h1/6*(RiK1 + 2*RiK2 + 2*RiK3 + RiK4);
        loglk1 = loglk1 + lnpoisspdf(F0(i+1),Fi(i+1));
        loglk1 = loglk1 + lnpoisspdf(R0(i+1),Ri(i+1));        
    end
    loglike1 = loglk1;
end