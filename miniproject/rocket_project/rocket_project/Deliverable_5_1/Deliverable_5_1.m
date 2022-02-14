addpath(fullfile('..', 'src'));

%% TODO: This file should produce all the plots for the deliverable

clc
clear all

%% START
Ts = 1/20; % Sample time
rocket = Rocket(Ts);

[xs, us] = rocket.trim();
sys = rocket.linearize(xs, us);
[sys_x, sys_y, sys_z, sys_roll] = rocket.decompose(sys, xs, us); % Design MPC controller

H = 5; % Horizon length in seconds
% % subsystem x
mpc_x = MPC_Control_x(sys_x, Ts, H);

% subsystem y
mpc_y = MPC_Control_y(sys_y, Ts, H);

% subsystem z
mpc_z = MPC_Control_z(sys_z, Ts, H);

% subsystem roll
mpc_roll = MPC_Control_roll(sys_roll, Ts, H);

% merge 4 sub systems
mpc = rocket.merge_lin_controllers(xs, us, mpc_x, mpc_y, mpc_z, mpc_roll);
% Setup reference function
Tf = 30;
ref = @(t_, x_) rocket.MPC_ref(t_, Tf);
% ref = @(t_, x_) [0 0  1 0];


x0 = zeros(12, 1);
rocket.mass = 1.783;

%without estimator
% [T, X, U, Ref] = rocket.simulate_f(x0, Tf, mpc, ref);

% with estimator
[T, X, U, Ref, Z_hat] = rocket.simulate_f_est_z(x0, Tf, mpc, ref, mpc_z, sys_z);

% Plot pose
rocket.anim_rate = 5; % Increase this to make the animation faster
ph = rocket.plotvis(T, X, U, Ref);

