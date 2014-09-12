function out1 = lnpoisspdf(x1,theta1)
    if (theta1 < 0 )
        res1 = -Inf;
    else
      res1 = -theta1 + x1*log(theta1)-sum(log(1:x1));
    end
  out1 = res1;