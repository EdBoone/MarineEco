%======================================================
%  Hierarchical Regression mean and variance
%
function [Bout, Varout] = hierbmean(Y1,X1,V1,Z1,B,VB);
  X2 = X1;
  V2i = inv(V1);
  Y2 = Y1;  
  VBi = inv(VB);
        XtXi = inv(X2'*V2i*X2+VBi);
        B1a = (XtXi)*(X2'*V2i*Y2 + VBi*Z1*B);
  Bout = B1a;
  Varout = XtXi;
