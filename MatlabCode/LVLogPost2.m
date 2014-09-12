% Lotka-Volterra Log Posterior
% 
function res1 = LVLogPost2(F01,R11, priorm1, priorsd1, ...
       rates,initial_cond,start_t1,end_t1)
    lpost1 = 0;
    lpost1 = lpost1 + LVLogLike2(F01,R11,rates,initial_cond,start_t1,end_t1);
    lpost1 = lpost1 + sum(normln(rates(1), priorm1(1), priorsd1(1)));
    lpost1 = lpost1 + sum(normln(rates(2), priorm1(2), priorsd1(2)));
    lpost1 = lpost1 + sum(normln(rates(3), priorm1(3), priorsd1(3)));
    lpost1 = lpost1 + sum(normln(rates(4), priorm1(4), priorsd1(4)));   
    res1 = lpost1;