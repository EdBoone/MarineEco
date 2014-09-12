function res1 = normln( z1, m1, s1)
  out1 = -0.5*log( 2*pi*s1.^2 ) - 0.5*( (z1 - m1)./s1 ).^2;
  res1 = out1;