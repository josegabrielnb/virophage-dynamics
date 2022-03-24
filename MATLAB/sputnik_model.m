% ODE model 2: non-independent entry, integration and cell death
% Sputnik mechanism

% State variables

% y(1): uninfected cells
% y(2): cells infected by giant virus
% y(3): cells infected by giant virus and exogenous virophage
% y(4): cells with an integrated virophage
% y(5): cells with an integrated virophage infected by a giant virus
% y(6): y(5) cells infected by an exogenous virophage
% y(7): giant virus population
% y(8): virophage population
% y(9): giant virus-virophage complex

% %Parameters
global r1 r2 K b1 b2 a1 a2 f lambda phi1 phi2 gamma1 gamma2 death k
r1 = 1; %intrinsic growth rate of the uninfected cells
r2 = 0.8; %intrinsic growth rate of the cells with a provirophage
K = 1e6; %carrying capacity of cells
b1 = 1e-7; %rate of infection by giant viruses
b2 = 1e-7; %rate of infection by complex
a1 = 1; %death rate of cells infected by giant virus
a2 = 0.9; %death rate of cells infected by giant virus and virophage  
f = 1; %inhibition of giant virus replication by the virophage (0-1)   
lambda = 1e-7; %rate of virophage integration into uninfected cells USE 0.001
phi1 = 100; %conversion of dying cells to giant virus
phi2 = 1000; %conversion of dying cells to virophages
gamma1 = 1.2; %rate of decay of giant viruses
gamma2 = 2.6; %rate of decay of virophages    
death = 0; % programmed cell death parameter in the cell population
k = 8e-7; %complex formation rate

%time span
tspan = [0 100];

%initial conditions
y0 = [1e5;0;0;0;0;0;1e5;1e5;0];
% % 
% % %numerical integration
tic
[t,y] = ode45(@odemodel, tspan, y0);
toc
% 
%plot results
semilogy(t,y(:,1),'-')
ylim([1 1e9])
hold on;
%semilogy(t,y(:,2),'-')
%semilogy(t,y(:,3),'-')
semilogy(t,y(:,4),'-')
%semilogy(t,y(:,5),'-')
%semilogy(t,y(:,6),'-')
semilogy(t,y(:,7),'-')
semilogy(t,y(:,8),'-')
lgd = legend('Cx','Ci','G','V','FontSize',16,'location','southeast');
hold off;
title("Cell-virus-virophage dynamics")
xlabel("time")
ylabel("population size")

%ODE model function

function dydt = odemodel(t,y)
    
    dydt = zeros(9,1);
    global r1 r2 K b1 b2 a1 a2 f lambda phi1 phi2 gamma1 gamma2 death k
    
    dydt(1) = r1*y(1)*(1-(y(1)+y(2)+y(3)+y(4)+y(5)+y(6))/K)-b1*y(1)*y(7)...
        -b2*y(1)*y(9);
    dydt(2) = b1*y(1)*y(7)-a1*y(2)-death*y(2);
    dydt(3) = b2*y(1)*y(9)-lambda*y(3)-a2*y(3)-death*y(3);
    dydt(4) = r2*y(4)*(1-(y(1)+y(2)+y(3)+y(4)+y(5)+y(6))/K)+lambda*y(3)...
        -b1*y(4)*y(7)-b2*y(4)*y(9);
    dydt(5) = b1*y(4)*y(7)-a1*y(5)-death*y(5);
    dydt(6) = b2*y(4)*y(9)-a2*y(6)-death*y(6);
    dydt(7) = a1*phi1*y(2)+a2*phi1*f*y(3)+a1*phi1*y(5)+a2*phi1*f*y(6)...
        -b1*y(1)*y(7)-b1*y(4)*y(7)-k*y(7)*y(8)-gamma1*y(7)+gamma2*y(9);
    dydt(8) = a2*phi2*y(3)+a1*phi2*y(5)+a2*phi2*y(6)-k*y(7)*y(8)...
        -gamma2*y(8)+gamma1*y(9);
    dydt(9) = k*y(7)*y(8)-b2*y(1)*y(9)-b2*y(4)*y(9)-(gamma1+gamma2)*y(9);
end