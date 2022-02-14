%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% EPFL | ME-425: Model Predictive Control | Exercise 6 Prob 2 & 3 %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Preliminaries

clear all
close all
clc

% System dynamics
A = [0.9752 1.4544; ... 
    -0.0327 0.9315];
B = [0.0248; 0.0327];

% Initial condition
x0 = [3; 0];

% Horizon and cost matrices
N = 10;
Q = 10 * eye(2);
R = 1;

% Constraints
umax = 1.75;
xmax = [5; 0.2];

% u in U = { u | Mu <= m }
M = [1;-1]; m = [umax; umax];
% x in X = { x | Fx <= f }
F = [1 0; 0 1; -1 0; 0 -1]; f = [xmax; xmax];

%% Defining the explicit MPC controller using MPT

% define system dynamics
sys = LTISystem('A', A, 'B', B);

% cost function
sys.x.penalty = QuadFunction(Q); 
sys.u.penalty = QuadFunction(R);

% constraints
sys.u.min = -umax;
sys.u.max = umax;
sys.x.min = -xmax;
sys.x.max = xmax;

% extract LQR gain and terminal set
Qf = sys.LQRPenalty.weight;
Xf = sys.LQRSet;

% terminal set and cost
sys.x.with('terminalPenalty');
sys.x.terminalPenalty = QuadFunction(Qf);
sys.x.with('terminalSet');
sys.x.terminalSet = Xf;

controller = MPCController(sys, N);
empc = controller.toExplicit();

figure;
empc.feedback.fplot();

%% Defining the MPC controller using YALMIP

x = sdpvar(2,N,'full');
u = sdpvar(1,N-1,'full');

con = (x(:,2) == A*x(:,1) + B*u(:,1)) + (M*u(:,1) <= m);
obj = u(:,1)'*R*u(:,1);
for i = 2:N-1
    con = con + (x(:,i+1) == A*x(:,i) + B*u(:,i));
    con = con + (F*x(:,i) <= f) + (M*u(:,i) <= m);
    obj = obj + x(:,i)'*Q*x(:,i) + u(:,i)'*R*u(:,i);
end
con = con + (Xf.A*x(:,N) <= Xf.b);
obj = obj + x(:,N)'*Qf*x(:,N);

% Compile the matrices
ctrl = optimizer(con, obj, sdpsettings('solver','sedumi'), x(:,1), u(:,1));
% ctrl = optimizer(con, obj, sdpsettings('solver','gurobi'), x(:,1), u(:,1));

%% Simulating the closed-loop system

solYALMIP.x(:,1) = x0;
solEMPC.x(:,1) = x0;

i = 1;
while norm(solYALMIP.x(:,end)) > 1e-3 % Simulate until convergence
    % YALMIP
    tic;
    [uopt,infeasible] = ctrl{solYALMIP.x(:,i)};
    solYALMIP.t(i) = toc;
    if infeasible == 1, error('Error in optimizer - could not solve the problem'); end
    solYALMIP.u(:,i) = uopt;
    solYALMIP.x(:,i+1) = A*solYALMIP.x(:,i) + B*solYALMIP.u(:,i);

    % eMPC
    tic;
    solEMPC.u(:,i) = empc.evaluate(solEMPC.x(:,i));
    solEMPC.t(i) = toc;
    solEMPC.x(:,i+1) = A*solEMPC.x(:,i) + B*solEMPC.u(:,i);

    if max(abs(solYALMIP.u(:,i) - solEMPC.u(:,i))) > 1e-4
        error('YALMIP and eMPC solution are different');
    end

    i = i + 1;
end

%% Plotting the results

figure;
hold on; grid on;

o = ones(1,size(solYALMIP.x,2));

subplot(3,1,1)
hold on; grid on;
plot(solYALMIP.x(1,:),'-k','markersize',20,'linewidth',2);
plot(1:size(solYALMIP.x,2),f(1)*o,'r','linewidth',2)
plot(1:size(solYALMIP.x,2),-f(3)*o,'r','linewidth',2)
ylabel('Position')

subplot(3,1,2)
hold on; grid on;
plot(solYALMIP.x(2,:),'-k','markersize',20,'linewidth',2);
plot(1:size(solYALMIP.x,2),f(2)*o,'r','linewidth',2)
plot(1:size(solYALMIP.x,2),-f(4)*o,'r','linewidth',2)
ylabel('Velocity')

subplot(3,1,3)
o = ones(1,size(solYALMIP.u,2));
hold on; grid on;
plot(solYALMIP.u,'k','markersize',20,'linewidth',2);
plot(1:size(solYALMIP.u,2),m(1)*o,'r','linewidth',2)
plot(1:size(solYALMIP.u,2),-m(2)*o,'r','linewidth',2)
ylabel('Input')

%% Plot solve times

figure;
semilogy(solYALMIP.t, 'r','linewidth',2);
hold on; grid on;
semilogy(solEMPC.t, 'g','linewidth',2);
legend('YALMIP','MPT3 eMPC')