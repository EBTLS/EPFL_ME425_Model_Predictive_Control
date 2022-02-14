classdef MPC_Control_z < MPC_Control

    properties
        A_bar, B_bar, C_bar % Augmented system for disturbance rejection
        L % Estimator gain for disturbance rejection
    end

    methods

        function mpc = MPC_Control_z(sys, Ts, H)
            mpc = mpc@MPC_Control(sys, Ts, H);

            [mpc.A_bar, mpc.B_bar, mpc.C_bar, mpc.L] = mpc.setup_estimator();
        end

        % Design a YALMIP optimizer object that takes a steady-state state
        % and input (xs, us) and returns a control input
        function ctrl_opti = setup_controller(mpc, Ts, H)

            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % INPUTS
            %   X(:,1)       - initial state (estimate)
            %   d_est        - disturbance estimate
            %   x_ref, u_ref - reference state/input
            % OUTPUTS
            %   U(:,1)       - input to apply to the system
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

            N = ceil(H / Ts); % Horizon steps

            [nx, nu] = size(mpc.B);

            % Targets (Ignore this before Todo 3.3)
            x_ref = sdpvar(nx, 1);
            u_ref = sdpvar(nu, 1);

            % Disturbance estimate (Ignore this before Part 5)
            d_est = sdpvar(1);

            % Predicted state and input trajectories
            X = sdpvar(nx, N);
            U = sdpvar(nu, N - 1);

            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % YOUR CODE HERE YOUR CODE HERE YOUR CODE HERE YOUR CODE HERE

            % NOTE: The matrices mpc.A, mpc.B, mpc.C and mpc.D are
            %       the DISCRETE-TIME MODEL of your system

            % SET THE PROBLEM CONSTRAINTS con AND THE OBJECTIVE obj HERE
            obj = 0;
            con = [];

            Q = 3*eye(nx);
            Q(2,2)=150;
            R = eye(nu);
            
            % ---- compute terminal cost ----
            % terminal weight
            Q_f = Q;
            R_f = R;

            % terminal controller
            [K_f, S, e] = dlqr(mpc.A, mpc.B, Q_f, R_f);

            % input constraints for time 1
            con = [con, (50 - 56.667) * ones(nu, 1) <= U(:, 1) <= (80 - 56.667) * ones(nu, 1)];

            % cost function for time 1
            obj = obj + (U(:, 1) - u_ref)' * R * (U(:, 1) - u_ref);

            % time 2 to N-1
            for i = 2:1:N - 1
                % input constraints
                con = [con, (50 - 56.667) * ones(nu, 1) <= U(:, i) <= (80 - 56.667) * ones(nu, 1)];
                % state constraints
                con = [con, [-inf; -inf] <= X(:, i) <= [inf; inf]];
                % state-space model constraints
                con = [con, X(:, i) == mpc.A * X(:, i - 1) + mpc.B * U(:, i - 1) + mpc.B * d_est];
                % object
                obj = obj + (X(:, i) - x_ref)' * Q * (X(:, i) - x_ref) + (U(:, i) - u_ref)' * R * (U(:, i) - u_ref);
            end

            % final state constraints
            con = [con, X(:, N) == mpc.A * X(:, N - 1) + mpc.B * U(:, N - 1) + mpc.B * d_est];

            % final state objective
            obj = obj + (X(:, N) - x_ref)' * S * (X(:, N) - x_ref);

            % YOUR CODE HERE YOUR CODE HERE YOUR CODE HERE YOUR CODE HERE
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

            % Return YALMIP optimizer object
            ctrl_opti = optimizer(con, obj, sdpsettings('solver', 'gurobi'), ...
            {X(:, 1), x_ref, u_ref, d_est}, U(:, 1));
        end

        % Design a YALMIP optimizer object that takes a position reference
        % and returns a feasible steady-state state and input (xs, us)
        function target_opti = setup_steady_state_target(mpc)

            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % INPUTS
            %   ref    - reference to track
            %   d_est  - disturbance estimate
            % OUTPUTS
            %   xs, us - steady-state target
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

            nx = size(mpc.A, 1);

            % Steady-state targets
            xs = sdpvar(nx, 1);
            us = sdpvar;

            % Reference position (Ignore this before Todo 3.3)
            ref = sdpvar;

            % Disturbance estimate (Ignore this before Part 5)
            d_est = sdpvar;

            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % YOUR CODE HERE YOUR CODE HERE YOUR CODE HERE YOUR CODE HERE
            % You can use the matrices mpc.A, mpc.B, mpc.C and mpc.D
            obj = 0;
            con = [50 - 56.6667 <= us <= 80 - 56.6667,
                xs == mpc.A * xs + mpc.B * us + mpc.B * d_est,
                ref == mpc.C * xs];
            obj = obj + us' * us;

            % YOUR CODE HERE YOUR CODE HERE YOUR CODE HERE YOUR CODE HERE
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

            % Compute the steady-state target
            target_opti = optimizer(con, obj, sdpsettings('solver', 'gurobi'), {ref, d_est}, {xs, us});
        end

        % Compute augmented system and estimator gain for input disturbance rejection
        function [A_bar, B_bar, C_bar, L] = setup_estimator(mpc)

            %%% Design the matrices A_bar, B_bar, L, and C_bar
            %%% so that the estimate x_bar_next [ x_hat; disturbance_hat ]
            %%% converges to the correct state and constant input disturbance
            %%%   x_bar_next = A_bar * x_bar + B_bar * u + L * (C_bar * x_bar - y);

            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % YOUR CODE HERE YOUR CODE HERE YOUR CODE HERE YOUR CODE HERE
            % You can use the matrices mpc.A, mpc.B, mpc.C and mpc.D
            nx = size(mpc.A, 1);
            nu = size(mpc.B, 2);
            nd = size(mpc.B, 2);
            ny = size(mpc.C, 1);
            A_bar = [mpc.A, mpc.B;
                zeros(nd, nx), ones(nd, nd)];
            B_bar = [mpc.B; zeros(nu, nd)];
            C_bar = [mpc.C, zeros(ny, nd)];
            poles_expected = [0.2, 0.23, 0.25];
            L = -place(A_bar', C_bar', poles_expected)';

            % YOUR CODE HERE YOUR CODE HERE YOUR CODE HERE YOUR CODE HERE
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        end

    end

end
