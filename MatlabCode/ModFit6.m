%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%
%%%  Fit the Lotka-Volterra Model
%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

close all
clear all

% Simple Lotka-Volterra Model
tic

% Set the number of MCMC Samples
MCMCSamps1 = 200;
MCMCSamps1long = 1000;
MCMCtrain1 = 4;
ChainNo1 = 3;

% Set the flag to determine whether or not to read in previous starting
% values:  0 means do not read in starting values
%          1 means read in starting values
StartStepFlag1 = 1;

% Read in the data
x1_data = csvread('Datasets/x1.csv');
T0 = csvread('Datasets/t1.csv');
F1 = csvread('Datasets/F1.csv');
R1 = csvread('Datasets/R1.csv');

% Set up the directory paths to save start step values and results
StartStepName1 = 'MCMC1/MCMC1S';




% Starting values for a1, g1, d1... b1 is a function of x1
% These are at true values... at the moment
a1a = 0.04;    % prey natural reproduction of Prey 
               % in absence of Pred rate a1 > 0
b1a = 0.001;   % Death Rate of Prey due to Pred  b1 > 0
e1a = 0.1;     % Preditor death rate  g1 > 0
d1a = 0.2;     % Death rate of Pred in absence of Prey  d1 > 0

% Keep a vector of ones around.
r1 = size(x1_data,1);
ones1 = ones(size(x1_data));

% Determine the begin and end of the times
start_t1 = min(T0);
end_t1 = max(T0);

% initial conditions
F0 = 15;               %  initial Fox population
R0 = 1000;             %  initial Rabbit population

% Set the initial conditions
initial_cond = [F0, R0];

if ( StartStepFlag1 == 0 )
    % Starting values for b0 and b1
    b01 = log(0.0001);
    b11 = (log(0.001) - log(0.0001))/10;
    
    % Convert b1 to regular values
    b1 = b01 + b11.*x1_data;
    a1 = log(a1a).*ones1;
    e1 = log(e1a).*ones1;
    d1 = log(d1a).*ones1;
    
    % Set the step values for the M-H sampler
    b1step = 0.001.*ones1;
    a1step = 0.001.*ones1;
    e1step = 0.001.*ones1;
    d1step = 0.001.*ones1;
    b01step = 0.0001;
    b11step = 0.0001;
else
    load(StartStepName1);
end

% Create containers to hold the sampled values
b1hold = zeros( MCMCSamps1, r1 );
a1hold = zeros( MCMCSamps1, r1 );
e1hold = zeros( MCMCSamps1, r1 );
d1hold = zeros( MCMCSamps1, r1 );
a2hold = zeros( MCMCSamps1, 1 );
a2sig1hold = zeros( MCMCSamps1, 1 );
e2hold = zeros( MCMCSamps1, 1 );
e2sig1hold = zeros( MCMCSamps1, 1 );
d2hold = zeros( MCMCSamps1, 1 );
d2sig1hold = zeros( MCMCSamps1, 1 );
b01hold = zeros(MCMCSamps1, 1);
b11hold = zeros(MCMCSamps1, 1);
fixed1holdx = zeros( MCMCSamps1, 1);

% Set up prior parameters
priorm1 = zeros(r1,4);
priorsd1 = ones(r1,4);
pmb01 = b01;
pmb11 = b11;  
psdb01 = 0.1;  
psdb11 = 0.1;

% a1 internal variance
a2v1 = 0.01;
e2v1 = 0.01;
d2v1 = 0.01;
sigmean1 = 1;
lambda1 = 1;



lpost1a = zeros(r1,1);
lpost1ta = zeros(r1,1);
% Initialize the logposterior
for i = 1:r1
    F01 = F1( : , i );
    R11 = R1( : , i ); 
    rates = [a1(i), b1(i), e1(i), d1(i) ];
    lpost1 = LVLogPost2ln(F01,R11, priorm1(i,:), priorsd1(i,:), ...
                    rates, b01, b11, pmb01, pmb11, psdb01, psdb11, ...
                    initial_cond, start_t1, end_t1);
    lpost1a(i) = lpost1;
end

for tMCMC1 = 1:MCMCtrain1
tic
  for iMCMC1 = 1:MCMCSamps1
    %=====================================================
    for i = 1:r1
        % Grab the correct data
        F01 = F1( : , i );
        R11 = R1( : , i ); 
        %=====================================================
        % Change a1
        a1c = a1(i) + randn(1,1)*a1step(i);
        rates = [a1c, b1(i), e1(i), d1(i) ];
        lpost1c = LVLogPost2ln(F01, R11, priorm1(i,:), priorsd1(i,:), ...
                        rates, b01, b11, pmb01, pmb11, psdb01, psdb11, ...
                        initial_cond, start_t1, end_t1);
        rand1u = rand(1,1);
        diff1 = lpost1c - lpost1a(i);
        if ( log(rand1u) < diff1 )
            a1(i) = a1c;
            lpost1a(i) = lpost1c;
        end
        %=====================================================
        % Change e1
        e1c = e1(i) + randn(1,1)*e1step(i);
        rates = [a1(i), b1(i), e1c, d1(i) ];
        lpost1c = LVLogPost2ln(F01,R11, priorm1(i,:), priorsd1(i,:), ...
                        rates, b01, b11, pmb01, pmb11, psdb01, psdb11, ...
                        initial_cond,start_t1,end_t1);
        rand1u = rand(1,1);
        diff1 = lpost1c - lpost1a(i);
        if( log(rand1u) < diff1)
            e1(i) = e1c;
            lpost1a(i) = lpost1c;
        end
        %=====================================================
        % Change c1
        d1c = d1(i) + randn(1,1)*d1step(i);
        rates = [a1(i), b1(i), e1(i), d1c ];
        lpost1c = LVLogPost2ln(F01,R11, priorm1(i,:), priorsd1(i,:), ...
                        rates, b01, b11, pmb01, pmb11, psdb01, psdb11, ...
                        initial_cond,start_t1,end_t1);        
        rand1u = rand(1,1);
        diff1 = lpost1c - lpost1a(i);
        if (log(rand1u) < diff1)
            d1(i) = d1c;
            lpost1a(i) = lpost1c;
        end
        %=====================================================
    end
    %=====================================================
    % Sample the second level parameters
    %
    % Do a1 first
    [a2m1, a2v1] = HierMeanVarSampUniv(a1,a2v1,log(a1a),1,sigmean1, lambda1); 
    [e2m1, e2v1] = HierMeanVarSampUniv(e1,e2v1,log(e1a),1,sigmean1, lambda1); 
    [d2m1, d2v1] = HierMeanVarSampUniv(d1,d2v1,log(d1a),1,sigmean1, lambda1); 
    %a2meansamp1 = (a2m1 + randn(1,1)*sqrt(a2v1));
    %a2sigsamp1 = sqrt(hiersig1(a1, ones1, a2meansamp1, sigmean1, lambda1)); 
    priorm1(:,1) = a2m1.*ones1;
    priorsd1(:,1) = a2v1.*ones1;
    priorm1(:,3) = e2m1.*ones1;
    priorsd1(:,3) = e2v1.*ones1;
    priorm1(:,4) = d2m1.*ones1;
    priorsd1(:,4) = d2v1.*ones1;
    
    % Update the logposterior
    for i = 1:r1
        F01 = F1( : , i );
        R11 = R1( : , i ); 
        rates = [a1(i), b1(i), e1(i), d1(i) ];
        lpost1 = LVLogPost2ln(F01,R11, priorm1(i,:), priorsd1(i,:), ...
                        rates, b01, b11, pmb01, pmb11, psdb01, psdb11, ...
                        initial_cond, start_t1, end_t1);
        lpost1a(i) = lpost1;
    end
    lpost1tot = sum(lpost1a);
    
    % Sample the regression parameters (intercept)
    b01c = b01 + b01step*randn(1,1);
    b1c = b01c + b11*x1_data;
    % Determine the candidate logposterior
    for i = 1:r1
        F01 = F1( : , i );
        R11 = R1( : , i ); 
        rates = [a1(i), b1c(i), e1(i), d1(i) ];
        lpost1c = LVLogPost2ln(F01,R11, priorm1(i,:), priorsd1(i,:), ...
                        rates, b01, b11, pmb01, pmb11, psdb01, psdb11, ...
                        initial_cond, start_t1, end_t1);
        lpost1ta(i) = lpost1c;
    end
    lpost1totc = sum(lpost1ta);
    rand1u = rand(1,1);
    diff1 = lpost1totc - lpost1tot;
    if (log(rand1u) < diff1)
        b1 = b1c;
        b01 = b01c;
        lpost1a = lpost1ta;
    end
    lpost1tot = sum(lpost1a);
   
    % Sample the regression parameters (slope)
    b11c = b11 + b11step*randn(1,1);
    if (b11c > 0)
      b1c = b01 + b11c*x1_data;
      % Determine the candidate logposterior
      for i = 1:r1
          F01 = F1( : , i );
          R11 = R1( : , i ); 
          rates = [a1(i), b1c(i), e1(i), d1(i) ];
          lpost1c = LVLogPost2ln(F01,R11, priorm1(i,:), priorsd1(i,:), ...
                          rates, b01, b11, pmb01, pmb11, psdb01, psdb11, ...
                          initial_cond, start_t1, end_t1);
          lpost1ta(i) = lpost1c;
      end
      lpost1totc = sum(lpost1ta);
      rand1u = rand(1,1);
      diff1 = lpost1totc - lpost1tot;
      if (log(rand1u) < diff1)
          b1 = b1c;
          b11 = b11c;
          lpost1a = lpost1ta;
      end
      lpost1tot = sum(lpost1a);
    end 

    fixed1 = LV1FixedPoint( a1(1), b1(1), e1(1), d1(1) );
    %=====================================================
    % Store the values 

    b1hold( iMCMC1, : ) = b1;
    a1hold( iMCMC1, : ) = a1;
    e1hold( iMCMC1, : ) = e1;
    d1hold( iMCMC1, : ) = d1;
    a2hold( iMCMC1, : ) = a2m1;
    a2sig1hold( iMCMC1, : ) = a2v1;
    e2hold( iMCMC1, : ) = e2m1;
    e2sig1hold( iMCMC1, : ) = e2v1;
    d2hold( iMCMC1, : ) = d2m1;
    d2sig1hold( iMCMC1, : ) = d2v1;
    b01hold( iMCMC1, : ) = b01;
    b11hold( iMCMC1, : ) = b11;
    fixed1hold( iMCMC1, : ) = fixed1;
    
    %toc
    if (mod(iMCMC1,200) == 0)
       display(iMCMC1)
    end
  end
  if( tMCMC1 == (MCMCtrain1 - 1))
    MCMCSamps1 = MCMCSamps1long; 
  end
%=========================================================
toc

% Update step values
b1step = max( [std(b1hold,0,1)'./2, b1step/2], [], 2 );
a1step = max( [std(a1hold,0,1)'./2, a1step/2], [], 2 );
e1step = max( [std(e1hold,0,1)'./2, e1step/2], [], 2 );
d1step = max( [std(d1hold,0,1)'./2, d1step/2], [], 2 );
b01step = max( [std(b01hold,0,1)'./2, b01step/2], [], 2 );
b11step = max( [std(b11hold,0,1)'./2, b11step/2], [], 2 );

save(StartStepName1,'b01','b11','b1','a1','e1', 'd1',...
     'b01step', 'b11step', 'b1step', 'a1step', 'e1step','d1step');


figure
subplot(2,2,1)
plot(b1hold)
subplot(2,2,2)
plot(a1hold)
subplot(2,2,3)
plot(e1hold)
subplot(2,2,4)
plot(d1hold)

figure
subplot(2,2,1)
plot(b01hold)
subplot(2,2,2)
plot(b11hold)
subplot(2,2,3)
plot(a2hold)
subplot(2,2,4)
plot(e2hold)

SampleName1 = ['MCMC1/MCMC1Out',num2str(ChainNo1),'R',num2str(tMCMC1)];

save(SampleName1,'b1hold','a1hold','e1hold','d1hold',...
     'b01hold', 'b11hold', 'a2hold', 'e2hold', 'd2hold', ...
     'fixed1hold');

end



