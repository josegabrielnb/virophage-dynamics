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

%Case f = 0 (total inhibition)

a = subs(eqns1,{r1,r2,K,b1,b2,a1,a2,f,lambda,phi1,phi2,gamma1,gamma2},...
    {1,0.8,1e6,1e-7,1e-6,1,0.9,0,1e-4,100,1000,1.2,2.6});

%Solve system of equations

sol0 = solve(a,y1,y2,y3,y4,y5,y6,y7,y8);

x0 = [vpa(sol0.y1),vpa(sol0.y2),vpa(sol0.y3),vpa(sol0.y4),vpa(sol0.y5),...
    vpa(sol0.y6),vpa(sol0.y7),vpa(sol0.y8)];

%x0(12,:) is a valid solution
%[0, 0, 0, 960012.98008304027602952674885044, 2250.9514710663194935683420108415, 16025.524304024258978738076857282, 173684.35313495316398533465825146, 6407500.1433902041312891334281628]

%Case f = 1

b = subs(eqns1,{r1,r2,K,b1,b2,a1,a2,f,lambda,phi1,phi2,gamma1,gamma2},...
    {1,0.8,1e6,1e-7,1e-6,1,0.9,1,1e-4,100,1000,1.2,2.6});

sol1 = solve(b,y1,y2,y3,y4,y5,y6,y7,y8);

x1 = [vpa(sol1.y1),vpa(sol1.y2),vpa(sol1.y3),vpa(sol1.y4),vpa(sol1.y5),...
    vpa(sol1.y6),vpa(sol1.y7),vpa(sol1.y8)];

%x1(6,:) is a valid solution
%[0, 0, 0, 121212.12121212121212121212121212, 2517.3793778921806897125931889197, 82704.926402477665053025086282137, 6348524.5840600715370884015945346, 29568222.579369172822649018563534]





