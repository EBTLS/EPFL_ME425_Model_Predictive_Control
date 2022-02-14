function opti_eval = NMPC_Control(rocket, H)

    import casadi.*
    opti = casadi.Opti(); % Optimization problem

    N = ceil(H / rocket.Ts); % MPC horizon
    nx = 12; % Number of states
    nu = 4; % Number of inputs

    % Decision variables (symbolic)
    X_sym = opti.variable(nx, N); % state trajectory
    U_sym = opti.variable(nu, N - 1); % control trajectory)

    % Parameters (symbolic)
    x0_sym = opti.parameter(nx, 1); % initial state
    ref_sym = opti.parameter(4, 1); % target position

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % YOUR CODE HERE YOUR CODE HERE YOUR CODE HERE YOUR CODE HERE

    % Define Terminal Cost
    [xs, us] = rocket.trim();
    sys_s = rocket.linearize(xs, us); % continuous
    sys_sd = c2d(sys_s, rocket.Ts); % discrete

%     compute terminal cost
    Q_f = 1 * ones(nx, nx);
%     Q_f(3,3) = 10; %wz
    Q_f(10, 10) = 5; % x
    Q_f(11, 11) = 5; % y
    Q_f(12, 12) = 15; % z
    Q_f(6, 6) = 20; % gamma
    R_f = ones(nu, nu);
    [K_f, S, ~] = dlqr(sys_sd.a, sys_sd.b, Q_f, R_f);

    % calculate steady state
    X_s = opti.variable(nx, 1);
    U_s = opti.variable(nu, 1);

    % State Decomposition
    angular_speed = X_sym([1:3], :); % angular velocites about the body axes
    angular = X_sym([4:6], :); % attitude of the body frame
    speed = X_sym([7:9], :); % velocity
    pos = X_sym([10:12], :); % position

    % Reference Decomposition
    ref_pos = ref_sym([1:3], :); % position reference
    ref_roll = ref_sym(4, :); % roll reference

    % Set weight for each state
    % 15 degree
    weight_angular_speed = 20*eye(3, 3);
    weight_angular_speed(2, 2) = 10;
    weight_angular_speed(3, 3) = 100;
    
    weight_angular = 40*eye(3, 3);
    weight_angular(3, 3) = 400;
    
    weight_speed = 200*eye(3, 3);
    weight_speed(2, 2) = 100;
    weight_speed(3, 3) = 200;
    
    weight_pos = 1000*eye(3, 3);
    weight_pos(3, 3) = 800;
    
    weight = blkdiag(weight_angular_speed, weight_angular, weight_speed, weight_pos)

    weight_input = eye(4, 4);
    weight_terminal = blkdiag(weight_angular_speed, weight_angular, weight_speed, weight_pos);

    
    % 50 degree
%     weight_angular_speed(1, 1) = 20; %wx
%     weight_angular_speed(2, 2) = 10; %wy
%     weight_angular_speed(3, 3) = 100; %wz
%     weight_angular = 40*eye(3, 3);
%     weight_angular(3, 3) = 400; %gamma
% 
%     weight_speed(1, 1) = 200;%vx
%     weight_speed(2, 2) = 100; %vy
%     weight_speed(3, 3) = 200; %vz
% 
%     weight_pos(1, 1) = 1000;
%     weight_pos(2, 2) = 1000;
%     weight_pos(3, 3) = 800;
% 
%     weight = blkdiag(weight_angular_speed, weight_angular, weight_speed, weight_pos)
% 
%     weight_input = eye(4, 4);
%     weight_terminal = blkdiag(weight_angular_speed, weight_angular, weight_speed, weight_pos);
%     
    % ---- objective ----
    cost = 0;
    % cost for state 1
    cost = cost + (U_sym(:, 1) - U_s)' * weight_input * (U_sym(:, 1) - U_s);
    % time 2 to N-1
    for i = 2:1:N - 1
        cost = cost + (X_sym(:, i) - X_s)' * weight * (X_sym(:, i) - X_s);
        cost = cost + (U_sym(:, i) - U_s)' * weight_input * (U_sym(:, i) - U_s);
    end

    % cost for terminal state
    cost = cost + (X_sym(:, N) - X_s)' * weight * (X_sym(:, N) - X_s);
%     cost = cost + (X_sym(:, N) - X_s)' * S * (X_sym(:, N) - X_s);

    opti.minimize(cost);

    % ---- multiple shooting ----
    %% with RK4 discretization
    f_discrete = @(x, u) RK4(x, u, rocket.Ts, @rocket.f);

    for k = 1:N - 1
        opti.subject_to(X_sym(:, k + 1) == f_discrete(X_sym(:, k), U_sym(:, k)));
    end

    % ---- input constraints ----
    U_lowerbound = [-0.26, -0.26, 50, -20]' * ones(1, N - 1);
    U_upperbound = [0.26, 0.26, 80, 20]' * ones(1, N - 1);
    opti.subject_to(U_lowerbound <= U_sym <= U_upperbound);

    % ---- state constraints ----
    angular_speed_lowerbound = -Inf * ones(3, N);
    angular_speed_upperbound = Inf * ones(3, N);
    angular_lowerbound = -Inf * ones(3, N);
    angular_upperbound = Inf * ones(3, N);
    angular_lowerbound(2, :) = deg2rad(-85) * ones(1, N);
    angular_upperbound(2, :) = deg2rad(85) * ones(1, N);
    speed_lowerbound = -Inf * ones(3, N);
    speed_upperbound = Inf * ones(3, N);
    pos_lowerbound = -Inf * ones(3, N);
    pos_upperbound = Inf * ones(3, N);

    opti.subject_to(angular_lowerbound <= angular <= angular_upperbound);

    % ---- steady state constraints ----
    opti.subject_to(angular_lowerbound <= X_s([4:6]) <= angular_upperbound);
    opti.subject_to(U_lowerbound <= U_s <= U_upperbound);
    opti.subject_to(X_s([10:12]) == ref_pos);
    opti.subject_to(X_s(6) == ref_roll);
    opti.subject_to(X_s == f_discrete(X_s, U_s)) % set steady state

    % ---- boundary conditions ----
    opti.subject_to(X_sym(:, 1) == x0_sym); % use initial state

    % YOUR CODE HERE YOUR CODE HERE YOUR CODE HERE YOUR CODE HERE
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    % ---- Setup solver ------
    ops = struct('ipopt', struct('print_level', 0, 'tol', 1e-3), 'print_time', false);
    opti.solver('ipopt', ops);

    % Create function to solve and evaluate opti
    opti_eval = @(x0_, ref_) solve(x0_, ref_, opti, x0_sym, ref_sym, U_sym);
end

function u = solve(x0, ref, opti, x0_sym, ref_sym, U_sym)

    % ---- Set the initial state and reference ----
    opti.set_value(x0_sym, x0);
    opti.set_value(ref_sym, ref);

    % ---- Solve the optimization problem ----
    sol = opti.solve();
    assert(sol.stats.success == 1, 'Error computing optimal input');

    u = opti.value(U_sym(:, 1));

    % Use the current solution to speed up the next optimization
    opti.set_initial(sol.value_variables());
    opti.set_initial(opti.lam_g, sol.value(opti.lam_g));
end
