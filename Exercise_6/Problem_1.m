%% ME425 Model Predictive Control
% Exercise 6 Problem 2
%%%%%%%%%%%%%
%% Clear
clear all
close all

%% setting
x_set=[-2,0.3,1.0]
u_set=zeros(3,1);

%% Optimization Problem Modelling

for i=1:1:3

    x=x_set(i);
    H=[4];
    f=(4*x-2);
    
    lb=0;
    ub=2;
    
    u=quadprog(H,f,[],[],[],[],lb,ub);
    u_set(i)=u;

end

u_set


