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
MCMCSamps1 = 1000;
MCMCtrain1 = 10;

% Read in the data
x1_data = csvread('Datasets/x1.csv');
T0 = csvread('Datasets/t1.csv');
F1 = csvread('Datasets/F1.csv');
R1 = csvread('Datasets/R1.csv');

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

% Starting values for b0 and b1
b0c = log(0.0001);
b1c = (log(0.001) - log(0.0001))/10;

% Convert b1 to regular values
b1 = exp( b0c + b1c*x1_data );
a1 = a1a*ones1;
e1 = e1a*ones1;
d1 = d1a*ones1;

% Create containers to hold the sampled values
b1hold = zeros( MCMCSamps1, r1 );
a1hold = zeros( MCMCSamps1, r1 );
e1hold = zeros( MCMCSamps1, r1 );
d1hold = zeros( MCMCSamps1, r1 );

% Set the step values for the M-H sampler
b0step = 0.001;
b1step = 0.00000001*ones1;
a1step = 0.00000001*ones1;
e1step = 0.00000001*ones1;
d1step = 0.00000001*ones1;

% Set up prior parameters
priorm1 = zeros(r1,4);
priorsd1 = ones(r1,4);

lpost1a = zeros(r1,1);
% Initialize the logposterior
for i = 1:r1
    F01 = F1( : , i );
    R11 = R1( : , i ); 
    rates = [a1(i), b1(i), e1(i), d1(i) ];
    lpost1 = LVLogPost2(F01,R11, priorm1(i,:), priorsd1(i,:), ...
                    rates, initial_cond, start_t1, end_t1);
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
        % Change b1
        b1c = b1(i) + randn(1,1)*b1step(i);
        rates = [a1(i), b1c, e1(i), d1(i) ];
        lpost1c = LVLogPost2(F01,R11, priorm1(i,:), priorsd1(i,:), ...
                        rates, initial_cond,start_t1,end_t1);
        rand1u = rand(1,1);
        diff1 = lpost1c - lpost1a(i);
        if ( log(rand1u) < diff1 )
            b1(i) = b1c;
            lpost1a(i) = lpost1c;
        end
        %=====================================================
        % Change a1
        a1c = a1(i) + randn(1,1)*a1step(i);
        rates = [a1c, b1(i), e1(i), d1(i) ];
        lpost1c = LVLogPost2(F01,R11, priorm1(i,:), priorsd1(i,:), ...
                        rates, initial_cond,start_t1,end_t1);
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
        lpost1c = LVLogPost2(F01,R11, priorm1(i,:), priorsd1(i,:), ...
                        rates, initial_cond,start_t1,end_t1);
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
        lpost1c = LVLogPost2(F01,R11, priorm1(i,:), priorsd1(i,:), ...
                        rates, initial_cond,start_t1,end_t1);        
        rand1u = rand(1,1);
        diff1 = lpost1c - lpost1a(i);
        if (log(rand1u) < diff1)
            d1(i) = d1c;
            lpost1a(i) = lpost1c;
        end
        %=====================================================
    end
    [a2m1, a2v1] = hierbmean(a1,ones1,a2v2*eye(r1),1,0,1); 
    prior1m(:,1) = (a2m1 + randn(1,1)*sqrt(a2v1)).*ones1;
    prior1sd(:,1) = sqrt(a2v1)*ones1;  % Fix this...
    %=====================================================
    % Store the values 
    b1hold( iMCMC1, : ) = b1;
    a1hold( iMCMC1, : ) = a1;
    e1hold( iMCMC1, : ) = e1;
    d1hold( iMCMC1, : ) = d1;
    
    display(iMCMC1)
    toc
end



%=========================================================
toc

% Update step values
b1step = max( [std(b1hold,0,1)', b1step/2], [], 2 );
a1step = max( [std(a1hold,0,1)', a1step/2], [], 2 );
e1step = max( [std(e1hold,0,1)', e1step/2], [], 2 );
d1step = max( [std(d1hold,0,1)', d1step/2], [], 2 );

figure
subplot(2,2,1)
plot(b1hold)
subplot(2,2,2)
plot(a1hold)
subplot(2,2,3)
plot(e1hold)
subplot(2,2,4)
plot(d1hold)

end