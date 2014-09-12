%=======================================================
%  Poisson log likelihood
function res1 = LVLogLike1(F01,R11,rates,initial_cond,start_t1,end_t1)
    [~,F0out,F1out] = LV1(rates,initial_cond,start_t1,end_t1);
    l1 = 0;
    l2 = 0;
    for i = 1:size(F01,1)
        l1 = l1 + lnpoisspdf(F01(i),F0out(i));
        l2 = l2 + lnpoisspdf(R11(i),F1out(i));
    end
    out1 = l1 + l2;
    res1 = out1;
 %[t1, F0out, F1out] = LV1(rates,y0,1,101);    