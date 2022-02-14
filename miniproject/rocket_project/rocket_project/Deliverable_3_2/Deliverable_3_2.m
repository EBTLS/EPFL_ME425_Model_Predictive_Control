addpath(fullfile('..', 'src'));

%% TODO: This file should produce all the plots for the deliverable


Ts= 1/20; % Sample time 
rocket = Rocket(Ts);

[xs, us] = rocket.trim(); 
sys= rocket.linearize(xs, us); 
[sys_x, sys_y, sys_z, sys_roll] = rocket.decompose(sys, xs, us);% Design MPC controller 

H = 5; % Horizon length in seconds 

% % subsystem x
mpc_x = MPC_Control_x(sys_x, Ts, H); 
ux = mpc_x.get_u([0;0;0;0],-5); % Get control input 

[T,X_sub,U_sub] = rocket.simulate(sys_x,[0,0,0,0],8,@mpc_x.get_u,-5);
ph=rocket.plotvis_sub(T,X_sub,U_sub,sys_x,xs,us,-5);

% % subsystem y
mpc_y = MPC_Control_y(sys_y, Ts, H); 
uy = mpc_y.get_u([0;0;0;0],-5); % Get control input 

[T,Y_sub,U_sub] = rocket.simulate(sys_y,[0,0,0,0],8,@mpc_y.get_u,-5);
ph=rocket.plotvis_sub(T,Y_sub,U_sub,sys_y,xs,us,-5);

% % subsystem z
mpc_z = MPC_Control_z(sys_z, Ts, H); 
uz = mpc_z.get_u([0;0],-5); % Get control input 

[T,Z_sub,U_sub] = rocket.simulate(sys_z,[0,0],8,@mpc_z.get_u,-5);
ph=rocket.plotvis_sub(T,Z_sub,U_sub,sys_z,xs,us,-5);

% % subsystem roll
mpc_roll = MPC_Control_roll(sys_roll, Ts, H); 
uroll = mpc_roll.get_u([0;0],deg2rad(45)); % Get control input 

[T,Roll_sub,U_sub] = rocket.simulate(sys_roll,[0,0],8,@mpc_roll.get_u,deg2rad(45));
ph=rocket.plotvis_sub(T,Roll_sub,U_sub,sys_roll,xs,us,deg2rad(45));
