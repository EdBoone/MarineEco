function [mean1, var1] = HierMeanVarSamp(Y1,X1,V1,Z1,B,VB, sigmean1, lambda1);
  X2 = X1;
  V2i = inv(V1);
  Y2 = Y1;  
  VBi = inv(VB);
        XtXi = inv(X2'*V2i*X2+VBi);
        B1a = (XtXi)*(X2'*V2i*Y2 + VBi*Z1*B);
  Bout = B1a;
  Varout = XtXi;
  res1 = Bout + randn(1,1)*sqrt(Varout);
  [n1, n2] = size(Y1);
  resid1 = Y1 - X1*res1;
  rhold = sigmean1 + resid1'*resid1;
  df1 = (lambda1 + (n1-n2))/2;
  ch1 = chisqr(1,1,df1);
  sig2 = rhold/ch1;
  mean1 = res1;
  var1 = sig2;
  
  
 
