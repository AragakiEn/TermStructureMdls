function [A,B,C,D,Mean0,Cov0,StateType ] = vasicekCS3(params,deltat,tau1,tau2,tau3,tau4,tau5,tau6,tau7)
% Detailed explanation goes here
% Three factor vasicek
rng(2022); %seed
kappa1 = params(1);
kappa2 = params(2);
kappa3 = params(3);
lam1 = params(4);
lam2 = params(5);
lam3 = params(6);
sigma1 = params(7);
sigma2 = params(8);
sigma3 = params(9);
h1 = params(10);
h2 = params(11);
h3 = params(12);
h4 = params(13);
h5 = params(14);
h6 = params(15);
h7 = params(16);
mu = params(17);

thetaq1 = sigma1*lam1/kappa1;
thetaq2 = sigma2*lam2/kappa2; 
thetaq3 = sigma3*lam3/kappa3;

b11 = (1-exp(-kappa1*tau1))/kappa1;
b12 = (1-exp(-kappa2*tau1))/kappa2;
b13 = (1-exp(-kappa3*tau1))/kappa3;

b21 = (1-exp(-kappa1*tau2))/kappa1;
b22 = (1-exp(-kappa2*tau2))/kappa2;
b23 = (1-exp(-kappa3*tau2))/kappa3;


b31 = (1-exp(-kappa1*tau3))/kappa1;
b32 = (1-exp(-kappa2*tau3))/kappa2;
b33 = (1-exp(-kappa3*tau3))/kappa3;

b41 = (1-exp(-kappa1*tau4))/kappa1;
b42 = (1-exp(-kappa2*tau4))/kappa2;
b43 = (1-exp(-kappa3*tau4))/kappa3;

b51 = (1-exp(-kappa1*tau5))/kappa1;
b52 = (1-exp(-kappa2*tau5))/kappa2;
b53 = (1-exp(-kappa3*tau5))/kappa3;

b61 = (1-exp(-kappa1*tau6))/kappa1;
b62 = (1-exp(-kappa2*tau6))/kappa2;
b63 = (1-exp(-kappa3*tau6))/kappa3;

b71 = (1-exp(-kappa1*tau7))/kappa1;
b72 = (1-exp(-kappa2*tau7))/kappa2;
b73 = (1-exp(-kappa3*tau7))/kappa3;

term1_1 = (thetaq1-(sigma1^2)/(2*kappa1^2));
term1_2 = (thetaq2-(sigma2^2)/(2*kappa2^2));
term1_3 = (thetaq3-(sigma3^2)/(2*kappa3^2));

term2_1 = sigma1^2/(4*kappa1);
term2_2 = sigma2^2/(4*kappa2);
term2_3 = sigma3^2/(4*kappa3);

a11 = term1_1*(b11-tau1)-(b11^2)*term2_1;
a12 = term1_2*(b12-tau1)-(b12^2)*term2_2;
a13 = term1_3*(b13-tau1)-(b13^2)*term2_3;

a21 = term1_1*(b21-tau2)-(b21^2)*term2_1;
a22 = term1_2*(b22-tau2)-(b22^2)*term2_2;
a23 = term1_3*(b23-tau2)-(b23^2)*term2_3;

a31 = term1_1*(b31-tau3)-(b31^2)*term2_1;
a32 = term1_2*(b32-tau3)-(b32^2)*term2_2;
a33 = term1_3*(b33-tau3)-(b33^2)*term2_3;

a41 = term1_1*(b41-tau4)-(b41^2)*term2_1;
a42 = term1_2*(b42-tau4)-(b42^2)*term2_2;
a43 = term1_3*(b43-tau4)-(b43^2)*term2_3;

a51 = term1_1*(b51-tau5)-(b51^2)*term2_1;
a52 = term1_2*(b52-tau5)-(b52^2)*term2_2;
a53 = term1_3*(b53-tau5)-(b53^2)*term2_3;

a61 = term1_1*(b61-tau6)-(b61^2)*term2_1;
a62 = term1_2*(b62-tau6)-(b62^2)*term2_2;
a63 = term1_3*(b63-tau6)-(b63^2)*term2_3;

a71 = term1_1*(b71-tau7)-(b71^2)*term2_1;
a72 = term1_2*(b72-tau7)-(b72^2)*term2_2;
a73 = term1_3*(b73-tau7)-(b73^2)*term2_3;


A = [1-kappa1*deltat,0,0,0;0,1-kappa2*deltat,0,0;0,0,1-kappa3*deltat,0;0,0,0,1]; %% 4by4
B = [sqrt(deltat*sigma1^2),0,0,0;0,sqrt(deltat*sigma2^2),0,0;0,0,sqrt(deltat*sigma3^2),0;0,0,0,0]; %4by4
C = [b11/tau1,b12/tau1,b13/tau1, (-a11-a12-a13)/tau1+mu; ...
    b21/tau2,b22/tau2,b23/tau2, (-a21-a22-a23)/tau2+mu; ...
    b31/tau3,b32/tau3,b33/tau3,(-a31-a32-a33)/tau3+mu; ...
    b41/tau4,b42/tau4,b43/tau4,(-a41-a42-a43)/tau4+mu; ...
    b51/tau5,b52/tau5,b53/tau5,(-a51-a52-a53)/tau5+mu; ...
    b61/tau6,b62/tau6,b63/tau6,(-a61-a62-a63)/tau6+mu; ...
    b71/tau7,b72/tau7,b73/tau7,(-a71-a72-a73)/tau7+mu];% 7by4
D = [h1,0,0,0,0,0,0; ...
    0,h2,0,0,0,0,0; ...
    0,0,h3,0,0,0,0; ...
    0,0,0,h4,0,0,0; ...
    0,0,0,0,h5,0,0; ...
    0,0,0,0,0,h6,0; ...
    0,0,0,0,0,0,h7];% 7by7


Mean0 = [0;0;0;1];
Cov0 = eye(4);
Cov0(4,4) = 0;
StateType = [0;0;0;1];
end

