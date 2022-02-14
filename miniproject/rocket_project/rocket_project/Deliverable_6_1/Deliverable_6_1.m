addpath(fullfile('..', 'src'));
% close all;
%% TODO: This file should produce all the plots for the deliverable

Ts=1/10;
rocket=Rocket(Ts);

H=1;
nmpc=NMPC_Control(rocket,H);

% % MPC reference with default maximum roll = 15 deg
Tf=30;
ref=@(t_,x_) rocket.MPC_ref(t_,Tf);

% % MPC reference with specified maximum roll = 50 deg
% roll_max=deg2rad(50);
% ref=@(t_,x_) rocket.MPC_ref(t_,Tf,roll_max);

x0 = zeros(12, 1);
[T,X,U,Ref]=rocket.simulate_f(x0,Tf,nmpc,ref);

% Plot pose
rocket.anim_rate = 5; % Increase this to make the animation faster
ph = rocket.plotvis(T, X, U, Ref);

ph.fig.Name = 'Merged lin. MPC in nonlinear simulation'; % Set a figure title
%% linear 50 degree

Ts = 1/20; % Sample time
rocket = Rocket(Ts);

[xs, us] = rocket.trim();
sys = rocket.linearize(xs, us);
[sys_x, sys_y, sys_z, sys_roll] = rocket.decompose(sys, xs, us); % Design MPC controller

H=5;
% ---- subsystem x ----
mpc_x = MPC_Control_x(sys_x, Ts, H); 

% ---- subsystem y ----
mpc_y = MPC_Control_y(sys_y, Ts, H); 

% ---- subsystem z ----
mpc_z = MPC_Control_z(sys_z, Ts, H); 

% ---- subsystem roll ----
mpc_roll = MPC_Control_roll(sys_roll, Ts, H); 

% ---- Merge 4 controllers ----
% Merge four sub−system controllers into one full−system controller 
mpc = rocket.merge_lin_controllers(xs, us, mpc_x, mpc_y, mpc_z, mpc_roll);


% Setup reference function 
Tf = 30; 
roll_max=deg2rad(50);
ref = @(t_, x_) rocket.MPC_ref(t_, Tf, roll_max);

x0 = zeros(12,1); 
[T, X, U, Ref] = rocket.simulate_f(x0, Tf, mpc, ref);

% Plot pose 
rocket.anim_rate = 5; % Increase this to make the animation faster 
ph = rocket.plotvis(T, X, U, Ref); 

ph.fig.Name = 'Merged lin. MPC in linear simulation'; % Set a figure title