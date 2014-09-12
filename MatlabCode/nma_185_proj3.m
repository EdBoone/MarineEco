function nma_185_proj3
%Solve Lotka-Volterra 2-ODE system

%
%
% Author: Nasser Abbasi
% May 26, 2003
%

initial_R = 1000;
initial_F = 15; 
nSteps    = 1000; 
h         = 1;

F=zeros(nSteps,1);
R=zeros(nSteps,1);

F(1)=initial_F;
R(1)=initial_R;

i=0;
time=0;
while(1)
    i=i+1;
    if(i>nSteps)
        break;
    end
    
    k_1R = f( time        , R(i)           , F(i));
    k_1F = g( time        , F(i)           , R(i));
    
    k_2R = f( time+0.5*h , R(i)+0.5*h*k_1R , F(i)+0.5*h*k_1F);
    k_2F = g( time+0.5*h , F(i)+0.5*h*k_1F , R(i)+0.5*h*k_1R);
    
    k_3R = f( time+0.5*h , R(i)+0.5*h*k_2R , F(i)+0.5*h*k_2F);
    k_3F = g( time+0.5*h , F(i)+0.5*h*k_2F , R(i)+0.5*h*k_2R);
    
    k_4R = f( time+h     , R(i)+h*k_3R     , F(i)+h*k_3F);
    k_4F = g( time+h     , F(i)+h*k_3F     , R(i)+h*k_3R);
    
    R(i+1) = R(i)+ (h/6)*(k_1R + 2*k_2R + 2*k_3R + k_4R);
    F(i+1) = F(i)+ (h/6)*(k_1F + 2*k_2F + 2*k_3F + k_4F);       
    
    time = time+h;    
end

figure;
plot(R,F);
title('Rabbit vs. Fox State Space');
xlabel('Rabbits population');
ylabel('Foxes population');

figure;
subplot(2,1,1); plot(R); xlabel('time'); ylabel('Rabbits');
title('population of rabbits vs time');

subplot(2,1,2); plot(F); xlabel('time'); ylabel('Foxes');
title('population of Foxes vs time');


figure;
plot(R); hold on; plot(F,'--');
title('Time domain representation of rabbits and foxes');
legend('Rabbits','Foxes');
xlabel('time'); ylabel('population');



%%%%%%%%%%%%%%
%
%
%%%%%%%%%%%%%%
function v=f(t,R,F)
a=0.04; b=0.0005; c=0.2; e=0.1;
v=a*R-b*R*F;
end
%%%%%%%%%%%%%%
%
%
%%%%%%%%%%%%%%
function v=g(t,F,R)
a=0.04; b=0.0005; c=0.2; e=0.1;
v=e*b*R*F-c*F;
end
end

