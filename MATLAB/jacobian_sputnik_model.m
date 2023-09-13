%% Jacobian of the Sputnik model

syms r1 r2 K b1 b2 a1 a2 f lambda phi1 phi2 gamma1 gamma2 k...
    y1 y2 y3 y4 y5 y6 y7 y8 y9

r1 = 1; %intrinsic growth rate of the uninfected cells
r2 = 0.8; %intrinsic growth rate of the cells with a provirophage
K = 1e6; %carrying capacity of cells
b1 = 1e-7; %rate of infection by giant viruses
b2 = 1e-7; %rate of infection by virophages
a1 = 1; %death rate of cells infected by giant virus
a2 = 0.9; %death rate of cells infected by giant virus and virophage  
%f = 1; %inhibition of giant virus replication by the virophage (0-1)   
lambda = 1e-4; %rate of virophage integration into uninfected cells
phi1 = 100; %conversion of dying cells to giant virus
phi2 = 1000; %conversion of dying cells to virophages
gamma1 = 1.2; %rate of decay of giant viruses
gamma2 = 2.6; %rate of decay of virophages
k = 8e-7; %complex formation rate
    
eq1 = r1*y1*(1-(y1+y2+y3+y4+y5+y6)/K)-b1*y1*y7...
        -b2*y1*y9;
eq2 = b1*y1*y7-a1*y2;
eq3 = b2*y1*y9-lambda*y3-a2*y3;
eq4 = r2*y4*(1-(y1+y2+y3+y4+y5+y6)/K)+lambda*y3...
        -b1*y4*y7-b2*y4*y9;
eq5 = b1*y4*y7-a1*y5;
eq6 = b2*y4*y9-a2*y6;
eq7 = a1*phi1*y2+a2*phi1*f*y3+a1*phi1*y5+a2*phi1*f*y6...
        -b1*y1*y7-b1*y4*y7-k*y7*y8-gamma1*y7+gamma2*y9;
eq8 = a2*phi2*y3+a1*phi2*y5+a2*phi2*y6-k*y7*y8...
        -gamma2*y8+gamma1*y9;
eq9 = k*y7*y8-b2*y1*y9-b2*y4*y9-(gamma1+gamma2)*y9;

eqns2 = [eq1,eq2,eq3,eq4,eq5,eq6,eq7,eq8,eq9];

J_sputnik = jacobian(eqns2,[y1,y2,y3,y4,y5,y6,y7,y8,y9]);

%Case f = 0, point attractor

%Solution 1
%x0 = [0,807766.45022954746461081539624803, 13039.610526179262382417363569978, 76919.51708268502669947322195886, 377.59407513598045628167837567714, 6.0954248290045529092383137639948, 35.960370791486622711664774683817, 161427.97862517964272162423611324, 25757017.760062257823807325264127, 857119.74428313810403229013147686];

%Solution 2
x0 = [0,0, 0, 0, 817046.28866955540519494066480919, 10939.257615238426031057215178226, 65327.845029583338293765063158208, 133887.85637899981791735269161309, 26078673.095100480824452911415152, 719605.01310598282592454376322191];

J0_sputnik_x0 = subs(J_sputnik,{f, y1,y2,y3,y4,y5,y6,y7,y8,y9},{x0});

Eig_sputnik_x0 = eig(J0_sputnik_x0);

vpa(Eig_sputnik_x0,5)
%Solution 1: the dominant eigenvalue  is negative -24.897 -> stable equilibrium 
%the complex conjugate eigenvalues have negative real parts -> point equilibrium

%Solution2: the dominant eigenvalue  is negative -25.153 -> stable equilibrium 
%the complex conjugate eigenvalues have negative real parts -> point equilibrium

%Case f = 1, oscillatory dynamics

%Solution 1
%x1 = [1,121167.42783688037684683817604539, 15785.390565662340911611878065947, 86926.572314270179216306551466905, 56.008153782818888139500905150426, 7.2966027100426189778970121628085, 40.18520383010274514472918174516, 1302775.0813455542642906914937413, 23619180.655261826224061194902326, 6457396.1118830871234391211794857];

%Solution 2
x1 = [1,0, 0, 0, 121212.12121212121212121212121212, 12634.086457778848545685751116936, 71573.473627647957564081976343104, 1042312.1327667550050190744671472, 24295591.350992188070694550197686, 5314330.4168528608491330867434755];

J0_sputnik_x1 = subs(J_sputnik,{f, y1,y2,y3,y4,y5,y6,y7,y8,y9},{x1});

Eig_sputnik_x1 = eig(J0_sputnik_x1);

vpa(Eig_sputnik_x1,5)
%Solution 1: the dominant eigenvalue  is negative -23.73 -> stable equilibrium 
% positive real part of the complex conjugate eigenvalue -> oscillations

%Solution 2: the dominant eigenvalue  is negative -24.089 -> stable equilibrium
% positive real part of a complex conjugate eigenvalue -> oscillations

