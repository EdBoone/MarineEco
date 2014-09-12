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
  b01hold2 = b01hold;
  b11hold2 = b11hold;
  
  for i = 2:5
      load(['MCMC1/MCMC1Out',num2str(i),'R4']);
      a1hold2 = [ a1hold2 ; a1hold ];
      b1hold2 = [ b1hold2 ; b1hold ];
      d1hold2 = [ d1hold2 ; d1hold ];
      e1hold2 = [ e1hold2 ; e1hold ];   
      b01hold2 = [ b01hold2 ; b01hold ];
      b11hold2 = [ b11hold2 ; b11hold ];
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

t1 = csvread('Datasets/t1.csv');  
F1 = csvread('Datasets/F1.csv');
R1 = csvread('Datasets/R1.csv');
  
n1size = size(t1,1);
initial_cond1 = [max([F1( n1size, 11),1]), max([R1(n1size, 11),1])];
Fout2 = zeros( n1, 100 );
Rout2 = zeros( n1, 100 );

  % Generate Predictions for the first site.
  for i = 1:n1
      initial_cond = [poissrnd(initial_cond1(1)), ...
                      poissrnd(initial_cond1(2))];
    % Separate the parameters for the model
    par = [exp(a1hold2(i,11)), ...
           exp(b1hold2(i,11)), ...
           exp(e1hold2(i,11)), ... 
           exp(d1hold2(i,11))];
    % Generate the means
    [t1a, F1a, R1a] = LV1( par, initial_cond, 1, 100 );
    % Add the poisson noise
    F1p = poissrnd(F1a);
    R1p = poissrnd(R1a);
    % Store the result
    Fout2(i,:) = F1p;
    Rout2(i,:) = R1p;

  end


  R1Q = quantile(Rout2, [0.025, 0.5, 0.975], 1 );
  F1Q = quantile(Fout2, [0.025, 0.5, 0.975], 1 );
  
  Fextinct1 = Fout2(:,100);
  1 - sum( Fextinct1 > 0 )/size(Fextinct1,1)
  
  figure;
  hold on;
    plot(t1, R1(:,11), 'b', t1, F1(:,11), 'r');
    plot(t1a+250, R1Q(2,:), 'b', ...
       t1a+250, R1Q(1,:), '-.b', ...
       t1a+250, R1Q(3,:), '-.b');
    plot(t1a+250, F1Q(2,:), 'r', ...
       t1a+250, F1Q(1,:), '-.r', ...
       t1a+250, F1Q(3,:), '-.r');
    plot( [250,250], [0,max(R1Q(3,:))], '--k');
    xlabel('Time');
    ylabel('Abundance');
    hleg1 = legend('Prey','Predator'); 
    print('-depsc',['Plots/Predfig1.eps']);
    
    
  figure;
  hold on;
    plot( R1Q(2,:), F1Q(2,:), 'b' );
    plot( R1Q(3,:), F1Q(3,:), '--b');
    plot( R1Q(1,:), F1Q(1,:), '--b');
    xlabel('Prey');
    ylabel('Predator');
    print('-depsc',['Plots/PredPhase11.eps']);
  
  figure;
  hold on;
    plot( R1(:,1), F1(:,1), 'Color',[ 0, 0, 0]);
    plot( R1(:,2), F1(:,2), 'Color',[ 1, 0, 1]);
    plot( R1(:,3), F1(:,3), 'Color',[ 1, 0, 0.5]);
    plot( R1(:,4), F1(:,4), 'Color',[ 0, 0.5, 1 ]);
    plot( R1(:,5), F1(:,5), 'Color',[ 0, 1, 1 ]);
    plot( R1(:,6), F1(:,6), 'Color',[ 0, 1, 0.5 ]);
    plot( R1(:,7), F1(:,7), 'Color',[ 1, 0, 0 ]);
    plot( R1(:,8), F1(:,8), 'Color',[ 0.5, 0, 1 ]);
    plot( R1(:,9), F1(:,9), 'Color',[ 0, 0, 1 ]);
    plot( R1(:,10), F1(:,10), 'Color',[ 0.5, 0, 1 ]);
    plot( R1(:,11), F1(:,11), 'Color',[ 0, 1, 0]);
    xlabel('Prey');
    ylabel('Predator');
    hleg1 = legend('Site 1','Site 2','Site 3', ...
                   'Site 4','Site 5','Site 6', ...
                   'Site 7','Site 8','Site 9', ...
                   'Site 10', 'Site 11');
    print('-depsc',['Plots/Phasefig1.eps']);

 figure;
 hold on;
    [f1,xi1] = ksdensity( a2hold );
    plot( xi1, f1 );
    xlabel( 'a' );
    ylabel( 'Density' );
    print('-depsc',['Plots/Adensity.eps']);

  figure;
  hold on;
    plot( 1:1000, a2hold);
    xlabel( 'Sample' );
    ylabel( 'a' );
    print('-depsc',['Plots/Atrace.eps']);
    
    
 figure;
 hold on;
    [f1,xi1] = ksdensity( d2hold );
    plot( xi1, f1 );
    xlabel( 'd' );
    ylabel( 'Density' );
    print('-depsc',['Plots/Ddensity.eps']);

    
  figure;
  hold on;
    plot( 1:1000, d2hold);
    xlabel( 'Sample' );
    ylabel( 'd' );
    print('-depsc',['Plots/Dtrace.eps']);

    
 figure;
 hold on;
    [f1,xi1] = ksdensity( b11hold2 );
    plot( xi1, f1 );
    xlabel( 'b' );
    ylabel( 'Density' );
    print('-depsc',['Plots/Bdensity.eps']);

    
  figure;
  hold on;
    plot( 1:5000, b11hold2 );
    xlabel( 'Sample' );
    ylabel( 'b' );
    print('-depsc',['Plots/Btrace.eps']);    
