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
MCMCSamps1 = 100;

% Read in the data
x1_data = csvread('Datasets/x1.csv');
T0 = csvread('Datasets/t1.csv');
F1 = csvread('Datasets/F1.csv');
R1 = csvread('Datasets/R1.csv');

% Starting values for a1, g1, d1... b1 is a function of x1
% These are at true values... at the moment
a1a = 0.04;    % prey natural reproduction of Prey 
              % in absence of Pred rate a1 > 0
b1a = 0.001;  % Death Rate of Prey due to Pred  b1 > 0
e1a = 0.1;     % Preditor death rate  g1 > 0
c1a = 0.2;     % Death rate of Pred in absence of Prey  d1 > 0

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
c1 = c1a*ones1;

% Create containers to hold the sampled values
b1hold = zeros( MCMCSamps1, r1 );
a1hold = zeros( MCMCSamps1, r1 );
g1hold = zeros( MCMCSamps1, r1 );
d1hold = zeros( MCMCSamps1, r1 );

% Set the step values for the M-H sampler
b0step = 0.001;
b1step = 0.00001;
a1step = 0.00001;
e1step = 0.00001;
c1step = 0.00001;

% Set up prior parameters
priorm1 = zeros(r1,4);
priorsd1 = ones(r1,4);


% Initialize the logposterior
lpost1 = 0;
for i = 1:r1
    F01 = F1( : , i );
    R11 = R1( : , i ); 
    rates = [a1(i), b1(i), e1(i), c1(i) ];
    lpost1 = lpost1 + LVLogPost2(F01,R11, priorm1(i,:), priorsd1(i,:), ...
                    rates, initial_cond, start_t1, end_t1);
end

for iMCMC1 = 1:MCMCSamps1
    %=====================================================
    for i = 1:r1
        % Grab the correct data
        F01 = F1( : , i );
        R11 = R1( : , i ); 
        %=====================================================
        % Change b1
        b1c = b1(i) + randn(1,1)*b1step;
        a1c = a1(i) + randn(1,1)*a1step;
        c1c = c1(i) + randn(1,1)*c1step;
        e1c = e1(i) + randn(1,1)*e1step;
        rates = [a1c, b1c, e1c, c1c ];
        lpost1c = LVLogPost2(F01,R11, priorm1, priorsd1, ...
                        rates, initial_cond,start_t1,end_t1);
        rand1u = rand(1,1);
        diff1 = lpost1c - lpost1;
        if ( log(rand1u) < diff1 )
            b1(i) = b1c;
            a1(i) = a1c;
            c1(i) = c1c;
            e1(i) = e1c;
            lpost1 = lpost1c;
        end
%         %=====================================================
%         % Change a1
%         a1c = a1(i) + randn(1,1)*a1step;
%         rates = [a1c, b1(i), e1(i), c1(i) ];
%         lpost1c = LVLogPost2(F01,R11, priorm1, priorsd1, ...
%                         rates, initial_cond,start_t1,end_t1);
%         rand1u = rand(1,1);
%         diff1 = lpost1c - lpost1;
%         if ( log(rand1u) < diff1 )
%             a1(i) = a1c;
%             lpost1 = lpost1c;
%         end
%         %=====================================================
%         % Change e1
%         e1c = e1(i) + randn(1,1)*e1step;
%         rates = [a1(i), b1(i), e1c, c1(i) ];
%         lpost1c = LVLogPost2(F01,R11, priorm1, priorsd1, ...
%                         rates, initial_cond,start_t1,end_t1);
%         rand1u = rand(1,1);
%         diff1 = lpost1c - lpost1;
%         if( log(rand1u) < diff1)
%             e1(i) = e1c;
%             lpost1 = lpost1c;
%         end
%         %=====================================================
%         % Change c1
%         c1c = c1(i) + randn(1,1)*c1step;
%         rates = [a1(i), b1(i), e1(i), c1c ];
%         lpost1c = LVLogPost2(F01,R11, priorm1, priorsd1, ...
%                         rates, initial_cond,start_t1,end_t1);        
%         rand1u = rand(1,1);
%         diff1 = lpost1c - lpost1;
%         if (log(rand1u) < diff1)
%             c1(i) = c1c;
%             lpost1 = lpost1c;
%         end
        %=====================================================
    end
        
    %=====================================================
    % Store the values 
    b1hold( iMCMC1, : ) = b1;
    a1hold( iMCMC1, : ) = a1;
    g1hold( iMCMC1, : ) = e1;
    d1hold( iMCMC1, : ) = c1;
    
    display(iMCMC1)
    toc
end



%=========================================================
toc

