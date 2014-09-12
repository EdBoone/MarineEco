function [mean1, var1] = HierMeanVarSampUniv(x1,sig1,mu1p,sig1p, sigmean1, lambda1);
  n1 = size(x1,1);
  sig2post = 1/( 1/sig1p + n1/sig1 );
  mu1post = sig2post*( mu1p/sig1p + sum(x1)/sig1 );
  res1 = mu1post + randn(1,1)*sqrt(sig2post);
  resid1 = x1 - res1*ones(n1,1);
  rhold = sigmean1 + resid1'*resid1;
  df1 = (lambda1 + n1);
  ch1 = chisqr(1,1,df1);
  sig2 = rhold/ch1;
  mean1 = res1;
  var1 = sig2;
  
  
 
