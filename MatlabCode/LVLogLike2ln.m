%=======================================================
%  Poisson log likelihood
function res1 = LVLogLike2ln(F01,R11, rates,initial_cond,start_t1,end_t1)
    [ ~ , ~ , ~ , loglike1 ] = LV2ln(F01, R11, rates, ...
                 initial_cond,start_t1,end_t1);
    res1 = loglike1;  