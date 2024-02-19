%{
16 paramers
kappa1 = params(1);
kappa2 = params(2);
kappa3 = params(3);
lam1 = params(4);
lam2 = params(5);
lam3 = params(6);
sigma1 = params(7);
sigma2 = params(8);
sigma3 = params(9);
h1 = prams(10);
h2 = params(11);
h3 = params(12);
h4 = params(13);
h5 = params(14);
h6 = params(15);
h7 = params(16);
mu = params(17);
ini is params, deltat,tau1,tau2,tau3,tau4,tau5,tau6,tau7
%}
rng(2022); %seed
mdl = ssm(@(params)vasicekCS3(params,30/360, 0.5,1,2,5,7,10,30)); % model1, specify different maturities 
mdl2 = ssm(@(params)vasicekCS3(params,30/360, 0.5,1,2,3,5,7,10)); % model2, specify different maturities

% starting guess
params0mdl = [0.65;0.1;0.1;0.15;0.09;0.02;0.02;0.02;0.02;0.4;0.4;0.01;0.01;0.01;0.01;0.01;0.1]; 

%%mean

%% estimation results;
%% put likelihood table into excel;

%% load your data before running
yield = readtable("xxxx.xlsx");
y2 = yield{1:369,["x6m", "x12m", "x24m","x36m","x60m","x84m","x120m"]};


Options = optimoptions(@fminunc,'MaxFunctionEvaluations',1.6e+04);
Options2 = optimoptions(@fmincon,'MaxFunctionEvaluations',1.6e+03);
[resmdl,estparmdl,estparcovmdl,loglmdl,outputmdl]= estimate(mdl2,y2,params0mdl,"lb",[-Inf;-Inf;-Inf;0;0;0;0;0;0;0;0;0;0;0;0;0;0],'Options',Options2, 'CovMethod','sandwich');

%{
%% filter the model implied state variable;
fitXmdl1 = filter(resmdl1, y);
fitXmdl2 = filter(resmdl2, y);

%% expected futures prices implied by the model;
fitYmdl1 = fitXmdl1 * transpose(resmdl1.C);
fitYmdl2 = fitXmdl2 * transpose(resmdl2.C);
%}



