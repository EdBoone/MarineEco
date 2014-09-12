%====================================================================
%
% Analysis of MCMC code
%
%====================================================================

close all;
clear all;

  load('MCMC1/MCMC1Out1R4');
  a1hold2 = a1hold;
  b1hold2 = b1hold;
  d1hold2 = d1hold;
  e1hold2 = e1hold;
  
  for i = 2:5
      load(['MCMC1/MCMC1Out',num2str(i),'R4']);
      a1hold2 = [ a1hold2 ; a1hold ];
      b1hold2 = [ b1hold2 ; b1hold ];
      d1hold2 = [ d1hold2 ; d1hold ];
      e1hold2 = [ e1hold2 ; e1hold ];      
  end
  
  


    n1 = size( b1hold2, 1 );
    fix1holdx = zeros( n1, 11 );
    fix1holdy = zeros( n1, 11 );
  for j = 1:11
    % Determine fixed points.
    for iMCMC1 = 1:n1

        a1 = a1hold2( iMCMC1, j);
        b1 = b1hold2( iMCMC1, j);
        d1 = d1hold2( iMCMC1, j);
        e1 = d1hold2( iMCMC1, j);

        fix1x = exp(a1)/exp(b1);
        fix1y = exp(d1)/( exp(e1)*exp(b1));

        fix1holdx( iMCMC1, j ) = fix1x;
        fix1holdy( iMCMC1, j ) = fix1y;

    end
    
  

  end
  
  
  for j = 1:11
    % Make a plot of the data
    Figure1 = figure(1);clf;
    scatter(fix1holdx( :, j ),fix1holdy( :, j ),'b');
    xlabel('x');
    ylabel('y');
    % Save the plots
    print('-depsc',['Plots/FixedPoint',num2str(j),'.eps']);
  end

%[bandwidth,density,X,Y]=kde2d(fix1hold);
% plot the data and the density estimate
%contour(X,Y,density), hold on
%plot(fix1hold(:,1),fix1hold(:,2),'r.','MarkerSize',5)
%hold off