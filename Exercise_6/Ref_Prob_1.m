%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% EPFL | ME-425: Model Predictive Control | Exercise 6 Prob 1 %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Preliminaries

clear all
close all
clc

%% Parametric Solution

x1 = linspace(-5, -1.5, 20);
x2 = linspace(-1.5, 0.5, 20);
x3 = linspace(0.5, 5, 20);

f1 = 5*x1.^2+4*x1+5;
f2 = 3*x2.^2-2*x2+0.5;
f3 = 5*x3.^2-4*x3+1;

u1 = 2+0*x1;
u2 = 0.5-x2;
u3 = 0*x3;

x_test = [x1 x2 x3];
f_para = [f1 f2 f3];
u_para = [u1 u2 u3];

%% Optimizer Solution

x  = sdpvar();
xp = sdpvar();
u  = sdpvar();

obj = x^2 + u^2 + xp^2;
con = [xp == 2*x + u - 1, 0 <= u, u <= 2];
opti = optimizer(con, obj, sdpsettings('solver','sedumi'), x, {obj, u});

f_opti = zeros(size(x_test));
u_opti = zeros(size(x_test));

for i = 1:length(x_test)
    sol = opti{x_test(i)};
    f_opti(i) = sol{1};
    u_opti(i) = sol{2};
end

%% Plot results

figure;
hold on; grid on;

subplot(2,1,1);
plot(x_test, f_para,'r','linewidth',2);
hold on;
plot(x_test, f_opti,'--b','linewidth',2);
ylabel('f*');

subplot(2,1,2);
plot(x_test, u_para,'r','linewidth',2);
hold on;
plot(x_test, u_opti,'--b','linewidth',2);
ylabel('u*');