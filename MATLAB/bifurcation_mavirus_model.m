%% Equilibria of the Mavirus model

% Mavirus model

syms r1 r2 K b1 b2 a1 a2 f lambda phi1 phi2 gamma1 gamma2...
    y1 y2 y3 y4 y5 y6 y7 y8

eq1 = r1*y1*(1-(y1+y2+y3+y4+y5+y6)/K)-b1*y1*y7-lambda*y1*y8 == 0;
eq2 = b1*y1*y7-b2*y2*y8-a1*y2 == 0;
eq3 = b2*y2*y8-a2*y3 == 0;
eq4 = r2*y4*(1-(y1+y2+y3+y4+y5+y6)/K)+lambda*y1*y8-b1*y4*y7 == 0;
eq5 = b1*y4*y7-b2*y5*y8-a1*y5 == 0;
eq6 = b2*y5*y8-a2*y6 == 0;
eq7 = a1*phi1*y2+a2*phi1*f*y3+a1*phi1*y5+a2*phi1*f*y6-b1*y1*y7-b1*y4*y7-gamma1*y7 == 0;
eq8 = a2*phi2*y3+a1*phi2*y5+a2*phi2*y6-b2*y2*y8-b2*y5*y8-gamma2*y8 == 0;

%Define system of equations

eqns1 = [eq1,eq2,eq3,eq4,eq5,eq6,eq7,eq8];

%Evaluate system with parameter values

f_test = 0.85; %set f_test

%Case f = f_test

a = subs(eqns1,{r1,r2,K,b1,b2,a1,a2,f,lambda,phi1,phi2,gamma1,gamma2},...
    {1,0.8,1e6,1e-7,1e-6,1,0.9,f_test,1e-4,100,1000,1.2,2.6});

%Solve system of equations

sol0 = solve(a,y1,y2,y3,y4,y5,y6,y7,y8);

x0 = [vpa(sol0.y1),vpa(sol0.y2),vpa(sol0.y3),vpa(sol0.y4),vpa(sol0.y5),...
    vpa(sol0.y6),vpa(sol0.y7),vpa(sol0.y8)];

%f=0.8
%x0(12,:) is a valid solution
%y0 = [0, 0, 0, 150829.9707693383680536588717017, 2529.686872489732336573744306471, 97612.246422081372180086322180681, 5992224.7674887242194374484944892, 34728022.165608882332188214068893];

%f=0.85
%x0(5,:) is a valid solution
y0 = [0, 0, 0, 142116.93192112370449639660696195, 2526.6176857179064880438965972944, 93442.333015120279251561331017595, 6095312.9390243048781119853233853, 33284853.576774057524277957428928];

%f=0.82
%x0(3,:) is a valid solution
%y0 = [0, 0, 0, 147217.0813867513240284073298385, 2528.4604599919615851972161887813, 95904.548033115881953634811581197, 6034799.2809611266594620851391322, 34137015.229448635019119336646706];

%f=0.81
%x0(13,:) is a valid solution
%y0 = [0, 0, 0, 149001.19633542704545354783782678, 2529.0739073829159801461525020821, 96751.606480859260884829271999225, 6013744.9862106462214518139013753, 34430170.497816722092960058060225];

%f=0.805
%x0(12,:) is a valid solution
%y0 = [0, 0, 0, 149909.89711245959287276956592147, 2529.3804533007123189875991309728, 97180.212821835250161175094288727, 6003044.0768992355571765419252706, 34578503.769774148361115433233064];

%f=0.8075
%x0(8,:) is a valid solution
%y0 = [0, 0, 0, 149454.13831527848120792558352769, 2529.2271957909070853646474837979, 96965.483240681217461546097467204, 6008409.2099859951539613093721705, 34504189.683648996040400286044831];

%f=0.8073
%x0(3,:) is a valid solution
%y0 = [0, 0, 0, 149490.4950135486470582480811731, 2529.2394575390186874140398482291, 96982.630170497874292164825315837, 6007981.0828673156796973844293027, 34510123.939936022101345936265269];


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

%Case f=f_test

x0 = [f_test, y0];

J0_mavirus_x0 = subs(J_mavirus,{f, y1,y2,y3,y4,y5,y6,y7,y8},{x0});

Eig_mavirus_x0 = eig(J0_mavirus_x0);

vpa(Eig_mavirus_x0, 5)




