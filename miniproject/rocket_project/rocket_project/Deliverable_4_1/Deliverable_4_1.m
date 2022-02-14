addpath(fullfile('..', 'src'));

%% TODO: This file should produce all the plots for the deliverable


%% START
Ts = 1/20; % Sample time
rocket = Rocket(Ts);

[xs, us] = rocket.trim();
sys = rocket.linearize(xs, us);
[sys_x, sys_y, sys_z, sys_roll] = rocket.decompose(sys, xs, us); % Design MPC controller

%% Initialize 4 Sub COntrollers

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
ref = @(t_, x_) rocket.MPC_ref(t_, Tf);
% ref = @(t_, x_) [0 0 1 0];
% ref = [0 0 0, 0 0 0, 0 0 0, 0 0 -1];

x0 = zeros(12,1); 
% rocket.mass = 1.783;
[T, X, U, Ref] = rocket.simulate_f(x0, Tf, mpc, ref);

% Plot pose 
rocket.anim_rate = 5; % Increase this to make the animation faster 
ph = rocket.plotvis(T, X, U, Ref); 

ph.fig.Name = 'Merged lin. MPC in nonlinear simulation (wihout estimator)'; % Set a figure title
