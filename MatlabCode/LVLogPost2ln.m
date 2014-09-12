% Lotka-Volterra Log Posterior
% 
function res1 = LVLogPost2ln(F01, R11, priorm1, priorsd1, ...
       rates, b01, b11, pmb01, pmb11, psdb01, psdb11, ...
       initial_cond, start_t1, end_t1)
    lpost1 = 0;
    lpost1 = lpost1 + LVLogLike2ln(F01,R11,rates,initial_cond,start_t1,end_t1);
    lpost1 = lpost1 + sum(normln(rates(1), priorm1(1), priorsd1(1)));
    lpost1 = lpost1 + sum(normln(rates(3), priorm1(3), priorsd1(3)));
    lpost1 = lpost1 + sum(normln(rates(4), priorm1(4), priorsd1(4))); 
    lpost1 = lpost1 + normln( b01, pmb01, psdb01 );  
    lpost1 = lpost1 + normln( b11, pmb11, psdb11 ); 
    res1 = lpost1;