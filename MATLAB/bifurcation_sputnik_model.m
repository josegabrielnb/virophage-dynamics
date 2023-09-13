%% Equilibria of the Sputnik model

% Sputnik model

syms r1 r2 K b1 b2 a1 a2 f lambda phi1 phi2 gamma1 gamma2 k...
    y1 y2 y3 y4 y5 y6 y7 y8 y9
    
eq1 = r1*y1*(1-(y1+y2+y3+y4+y5+y6)/K)-b1*y1*y7-b2*y1*y9 == 0;
eq2 = b1*y1*y7-a1*y2 == 0;
eq3 = b2*y1*y9-lambda*y3-a2*y3 == 0;
eq4 = r2*y4*(1-(y1+y2+y3+y4+y5+y6)/K)+lambda*y3-b1*y4*y7-b2*y4*y9 == 0;
eq5 = b1*y4*y7-a1*y5 == 0;
eq6 = b2*y4*y9-a2*y6 == 0;
eq7 = a1*phi1*y2+a2*phi1*f*y3+a1*phi1*y5+a2*phi1*f*y6...
        -b1*y1*y7-b1*y4*y7-k*y7*y8-gamma1*y7+gamma2*y9 == 0;
eq8 = a2*phi2*y3+a1*phi2*y5+a2*phi2*y6-k*y7*y8...
        -gamma2*y8+gamma1*y9 == 0;
eq9 = k*y7*y8-b2*y1*y9-b2*y4*y9-(gamma1+gamma2)*y9 == 0;

%Define system of equations

eqns1 = [eq1,eq2,eq3,eq4,eq5,eq6,eq7,eq8,eq9];

%Evaluate system with parameter values

f_test = 0.86; %set f_test

%Case f = f_test

a = subs(eqns1,{r1,r2,K,b1,b2,a1,a2,f,lambda,phi1,phi2,gamma1,gamma2,k},...
    {1,0.8,1e6,1e-7,1e-7,1,0.9,f_test,1e-4,100,1000,1.2,2.6,8e-7});

%Solve system of equations

sol0 = solve(a,y1,y2,y3,y4,y5,y6,y7,y8,y9, 'Real', true);

x0 = [vpa(sol0.y1),vpa(sol0.y2),vpa(sol0.y3),vpa(sol0.y4),vpa(sol0.y5),...
    vpa(sol0.y6),vpa(sol0.y7),vpa(sol0.y8),vpa(sol0.y9)];

%f=0.8
%x0(4,:) is a valid solution
%y0 = [0, 0, 0, 146640.55579883521313414847608888, 12564.709892728087087743777805631, 84606.235067917121390126398511217, 856837.30699743367802210675232608, 28897396.94472418894692072624352, 5192670.6869267229490817440284281];

%f=0.85
%x0(4,:) is a valid solution
%y0 = [0, 0, 0, 139195.99915031113152687566879497, 12584.677109983144125084726439454, 80950.968004845583161476435308023, 904097.61679957126827296286179977, 27599678.633666770117845145746485, 5234049.2290793098612195424938607];

%f=0.8305
%x0(4,:) is a valid solution
%y0 = [0, 0, 0, 141996.21479459212848211515281745, 12577.109863193166871314512998696, 82341.209261539189218610323222411, 885735.57269726328570493756288183, 28092625.347991494507576278478254, 5218948.1559481408377187425248097];

%f=0.82
%x0(4,:)
%y0 = [0, 0, 0, 143557.38238233239938924689438068, 12572.919514062749493463207541206, 83108.243271544174192401520908733, 875811.4215664396325827953374779, 28364927.875028460567126540239449, 5210280.2170900457828163116798771];

%f=0.8301
%x0(4,:)
%y0 = [0, 0, 0, 142054.98773773704710804146349639, 12576.951724779954586992931523365, 82370.189792547821933069208202919, 885357.98179783834512516573209776, 28102909.387861318967617477872729, 5218624.9841616430658500054421208];

%f=0.8302
%x0(4,:)
%y0 = [0, 0, 0, 142040.28936445723542302202013195, 12576.991270437058143352919854286, 82362.942891603083573939215883356, 885452.38303240184068347856203978, 28100337.719043377236529710599386, 5218705.8287556191421940081910035];

%f=0.8303
%x0(4,:)
%y0 = [0, 0, 0, 142025.594416965352771238335369, 12577.030808724444882627487806552, 82355.697169617245006326091129563, 885546.78192722167037411601961699, 28097766.489523258730046142234497, 5218786.6389103219883443486659421];

%f=0.86
%x0(4,:) is a valid solution
y0 = [0, 0, 0, 137807.72651317459793878964945502, 12588.450903712776978703713692826, 80254.8206426285282218396754726, 913479.32530541418481888580853488, 27353133.518820036295687306236822, 5241312.6902184585900664498825015];

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

%Case f = f_test

x0 = [f_test, y0];

J0_sputnik_x0 = subs(J_sputnik,{f, y1,y2,y3,y4,y5,y6,y7,y8,y9},{x0});

Eig_sputnik_x0 = eig(J0_sputnik_x0);

vpa(Eig_sputnik_x0,5)
