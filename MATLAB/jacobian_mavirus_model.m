%% Jacobian of the Mavirus model

syms f y1 y2 y3 y4 y5 y6 y7 y8

r1 = 1; %intrinsic growth rate of the uninfected cells
r2 = 0.8; %intrinsic growth rate of the cells with a provirophage
K = 1e6; %carrying capacity of cells
b1 = 1e-7; %rate of infection by giant viruses
b2 = 1e-6; %rate of infection by virophages
a1 = 1; %death rate of cells infected by giant virus
a2 = 0.9; %death rate of cells infected by giant virus and virophage  
%f = 1; %inhibition of giant virus replication by the virophage (0-1)   
lambda = 1e-4; %rate of virophage integration into uninfected cells
phi1 = 100; %conversion of dying cells to giant virus
phi2 = 1000; %conversion of dying cells to virophages
gamma1 = 1.2; %rate of decay of giant viruses
gamma2 = 2.6; %rate of decay of virophages

eq1 = r1*y1*(1-(y1+y2+y3+y4+y5+y6)/K)-b1*y1*y7...
        -lambda*y1*y8;
eq2 = b1*y1*y7-b2*y2*y8-a1*y2;
eq3 = b2*y2*y8-a2*y3;
eq4 = r2*y4*(1-(y1+y2+y3+y4+y5+y6)/K)+lambda*y1*y8...
        -b1*y4*y7;
eq5 = b1*y4*y7-b2*y5*y8-a1*y5;
eq6 = b2*y5*y8-a2*y6;
eq7 = a1*phi1*y2+a2*phi1*f*y3+a1*phi1*y5+a2*phi1*f*y6...
        -b1*y1*y7-b1*y4*y7-gamma1*y7;
eq8 = a2*phi2*y3+a1*phi2*y5+a2*phi2*y6...
        -b2*y2*y8-b2*y5*y8-gamma2*y8;

eqns1 = [eq1,eq2,eq3,eq4,eq5,eq6,eq7,eq8];

J_mavirus = jacobian(eqns1,[y1,y2,y3,y4,y5,y6,y7,y8]);

%Case f = 0, point attractor

x0 = [0, 0, 0, 0, 960012.98008304027602952674885044, 2250.9514710663194935683420108415, 16025.524304024258978738076857282, 173684.35313495316398533465825146, 6407500.1433902041312891334281628];

J0_mavirus_x0 = subs(J_mavirus,{f, y1,y2,y3,y4,y5,y6,y7,y8},{x0});

Eig_mavirus_x0 = eig(J0_mavirus_x0);

vpa(Eig_mavirus_x0, 5)
% The dominant eigenvalue -640.75 is negative -> stable equilibrium
% complex eigenvalues with negative real parts -> damping wave

%Case f = 1, oscillatory dynamics

x1 = [1, 0, 0, 0, 121212.12121212121212121212121212, 2517.3793778921806897125931889197, 82704.926402477665053025086282137, 6348524.5840600715370884015945346, 29568222.579369172822649018563534];
J0_mavirus_x1 = subs(J_mavirus,{f, y1,y2,y3,y4,y5,y6,y7,y8},{x1});

Eig_mavirus_x1 = eig(J0_mavirus_x1);

vpa(Eig_mavirus_x1, 5)
% The dominant eigenvalue -2956.7 is negative -> stable equilibrium
% complex eigenvalues with positive real parts -> oscillations 
