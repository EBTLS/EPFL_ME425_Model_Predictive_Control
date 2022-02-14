% ME-425 : Model Predictive Control
% Exercise sheet 5
% 
% Exercise Global Parameter Setting
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

A=[0.7115,-0.4345;
    0.4345,0.8853];

B=[0.2173;
    0.0573];

C=[0,1];

Bd=[0;
    0];

Cd=[1];

F_u=[1,0;
    0,-1];

f_u=[3;
    3];

umin=-3;
umax=3;

x0=[1;2];
d=0.2;
x_hat=[3;0];
d_hat=0;

nx = size(A, 1);
nu = size(B, 2);
nd = size(Bd, 2);

y=C*x0+d;