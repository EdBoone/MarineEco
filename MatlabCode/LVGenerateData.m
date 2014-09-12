%%%  Lotka-Volterra Data Generation


% Starting values for a1, g1, d1... b1 is a function of x1
% These are at true values... at the moment
a1a = 0.04;    % prey natural reproduction of Prey 
               % in absence of Pred rate a1 > 0
b1a = 0.001;   % Death Rate of Prey due to Pred  b1 > 0
e1a = 0.1;     % Preditor death rate  g1 > 0
c1a = 0.2;     % Death rate of Pred in absence of Prey  d1 > 0

% initial conditions
F0 = 15;               %  initial Fox population
R0 = 1000;             %  initial Rabbit population

% Starting values for b0 and b1
b0c = log(0.0001);
b1c = (log(0.001) - log(0.0001))/10;

% Generate the x data
x1_data = 0:10;

% Add noise to the parameters
bx1 = b0c + b1c*x1_data + randn(1,11)*0.01;
b1 = exp(bx1);
a1 = a1a + randn(1,11)*0.001;
e1 = e1a + randn(1,11)*0.01;
c1 = c1a + randn(1,11)*0.01;

% Get parameters to run the model
r1 = size(b1,2);

% Set up model parameters
start_t = 1;
end_t = 250;

% Set the initial conditions
initial_cond = [F0, R0];

% Create containers for the outputs.
tout1 = start_t:end_t;
Rout2 = zeros(end_t,11);
Fout2 = zeros(end_t,11);

for i = 1:r1
    % Separate the parameters for the model
    par = [a1(i), b1(i), e1(i), c1(i)];
    % Generate the means
    [t1, F1, R1] = LV1( par, initial_cond, start_t, end_t );
    % Add the poisson noise
    F1p = poissrnd(F1);
    R1p = poissrnd(R1);
    % Store the result
    Fout2(:,i) = F1p;
    Rout2(:,i) = R1p;
    % Make a plot of the data
    Figure1=figure(1);clf;
    set(Figure1,'defaulttextinterpreter','latex');
    plot(t1,F1p,'-r',t1,R1p,'-b');
    xlabel('Time');
    ylabel('Abundance');
    title(['$\beta$= ',num2str(b1(i))]);
    hleg1 = legend('Predator','Prey');
    % Save the plots
    print('-depsc',['Plots/Rawfig',num2str(i),'.eps']);
end

% Write the generated data files out.
csvwrite('Datasets/t1.csv',tout1');
csvwrite('Datasets/F1.csv',Fout2);
csvwrite('Datasets/R1.csv',Rout2);
csvwrite('Datasets/x1.csv',x1_data');


% %%
% figure()
% plot(u1,'k');
% hold on
% upredq=quantile(upred1hold,[.025,.975],1);
% 
% % for i=1:nBOOTSamples
% %     plot(upred1hold(i,:),'g');
% % end
% plot(upredq(1,:),'g');
% plot(upredq(2,:),'g');
% plot(u1,'k');
% title('Displacement')
% hold off
% 
% figure()
% plot(v1,'k');
% hold on
% vpredq=quantile(vpred1hold,[.025,.975],1);
% 
% % for i=1:nBOOTSamples
% %     plot(upred1hold(i,:),'g');
% % end
% plot(vpredq(1,:),'g');
% plot(vpredq(2,:),'g');
% plot(v1,'k');
% title('Velocity')
% hold off
% 
% 
