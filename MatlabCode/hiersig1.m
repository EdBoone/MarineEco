function sig2 = hiersig1(Y1, X1, b1, sigmean1, lambda1);
  [n1, n2] = size(Y1);
  resid1 = Y1 - X1*b1;
  rhold = sigmean1 + resid1'*resid1;
  df1 = (lambda1 + (n1-n2))/2;
  ch1 = chisqr(1,1,df1);
  sig2 = rhold/ch1;
  