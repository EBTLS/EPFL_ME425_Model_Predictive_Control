% ME-425 : Model Predictive Control
% Exercise sheet 5
%
% Exercise 2

%% Parameter Setting

yalmip('clear')
% clear all
% close all
% clc

global_setting;

%% optimizer settings

Q = eye(2);
R = 1;

N_horizon = 5;
N = 50;

%% Design the Oberserver
%% State Augment
A_bar = [A, Bd; ...
    zeros(1, nx), 1];
B_bar = [B; zeros(1, nd)];
C_bar = [C, Cd];

L = -place(A_bar', C_bar', [0.5, 0.6, 0.7])';
P = dlyap(A, Q);

eigs(A_bar + L * C_bar)

%% prepare to solve the problem
solution.x_hist = zeros(nx, N);
solution.x_hist(:, 1) = x0;
solution.x_hat_hist = zeros(nx, N);
solution.x_hat_hist(:, 1) = x_hat;
solution.d_hat_hist = zeros(1, N);
solution.d_hat_hist(:, 1) = d_hat;
solution.u_hist = zeros(nu, N - 1);

for i = 1:N - 1
    fprintf('step %i \n', i);
    [xs, us] = compute_sp(A, B, C, 1, solution.d_hat_hist(:, i), umin, umax);

    %% Set up the MPC cost and constraints using the computed set-point
    constraints = [];
    objective = 0;
    
    x = sdpvar(nx,N_horizon);
    u = sdpvar(nu,N_horizon);

    for k = 1:N_horizon - 1
        objective = objective + (x(:,k) - xs)' * Q * (x(:,k) - xs) + (u(:,k) - us)' * R * (u(:,k) - us);
        constraints = [constraints, umin <= u(:,k) <= umax, x(:,k+1) == A * x(:,k) + B * u(:,k)];
    end

    objective = objective + (x(:,N_horizon) - xs)' * P * (x(:,N_horizon) - xs);


%     x = sdpvar(repmat(nx,1,N_horizon),repmat(1,1,N_horizon));
%     u = sdpvar(repmat(nu,1,N_horizon),repmat(1,1,N_horizon));
% 
%     for i = 1:N_horizon - 1
%         objective = objective + (x{i} - xs)' * Q * (x{i} - xs) + (u{i} - us)' * R * (u{i} - us);
%         constraints = [constraints, umin <= u{i} <= umax, x{i + 1} == A * x{i} + B * u {i}];
%     end
% 
%     objective = objective + (x{N_horizon} - xs)' * P * (x{N_horizon} - xs);

    %% Optimize
    ops = sdpsettings('verbose', 0);
    ctrl = optimizer(constraints, objective, ops, x(:,1), u(:,1));

    [uopt,diagnostics] = ctrl{solution.x_hat_hist(:, i)};

    if diagnostics == 0
        fprintf("yalmip find available solution \n");
    else
        fprintf("no solution for yalmip \n");
    end

    solution.u_hist(:, i) = uopt;
    
    %% Apply to the plant (True x and y calculation)
    solution.x_hist(:, i + 1) = A * solution.x_hist(:, i) + B * solution.u_hist(:, i);

    %% Off-set observer (x, d estimate in the system)
    x_bar = [solution.x_hat_hist(:, i); ...
        solution.d_hat_hist(:, i)];

    x_bar_next = A_bar * x_bar + B_bar * solution.u_hist(:, i) + L * (C_bar * x_bar - y);

    solution.x_hat_hist(:, i + 1) = x_bar_next(1:nx);
    solution.d_hat_hist(:, i + 1) = x_bar_next(nx + 1);
    
    y = C * solution.x_hist(:, i+1) + d;


end

figure
plot(solution.u_hist); hold on
plot(umax * ones(size(solution.u_hist)), '--');
plot(umin * ones(size(solution.u_hist)), '--');
legend('u', 'u_max', 'u_min');

figure
plot(solution.x_hist(1, :), 'r'); hold on
plot(solution.x_hist(2, :), 'b');
plot(solution.x_hat_hist(1, :), 'r--'); hold on
plot(solution.x_hat_hist(2, :), 'b--');
legend('x_1', 'x_2', '\hat{x}_1', '\hat{x}_2')

figure
plot(d * ones(size(solution.d_hat_hist)), '--'); hold on
plot(solution.d_hat_hist);
legend('d', '\hat{d}')

figure
plot(C * solution.x_hist + d * ones(size(C * solution.x_hist))); hold on
plot(1 * ones(size(C * solution.x_hist)));
legend('y', 'r');

%% COMPUTE_SP compute a feasible set-point.
function [xs, us] = compute_sp(A, B, C, r, d, umin, umax)

    nx = size(A, 1);
    nu = size(B, 2);

    u = sdpvar(nu, 1);
    x = sdpvar(nx, 1);

    constraints = [umin <= u <= umax, ...
                    x == A * x + B * u, ...
                    r == C * x + d];

    objective = u^2;
    diagnostics = solvesdp(constraints, objective, sdpsettings('verbose', 0));

    if diagnostics.problem == 0
        % Good!
    elseif diagnostics.problem == 1
        throw(MException('', 'Infeasible'));
    else
        throw(MException('', 'Something else happened'));
    end

    xs = double(x);
    us = double(u);
end
