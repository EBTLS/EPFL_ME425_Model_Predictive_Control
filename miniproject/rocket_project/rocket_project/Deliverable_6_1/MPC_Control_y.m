classdef MPC_Control_y < MPC_Control

    methods
        % Design a YALMIP optimizer object that takes a steady-state state
        % and input (xs, us) and returns a control input
        function ctrl_opti = setup_controller(mpc, Ts, H)

            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % INPUTS
            %   X(:,1)       - initial state (estimate)
            %   x_ref, u_ref - reference state/input
            % OUTPUTS
            %   U(:,1)       - input to apply to the system
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

            N = ceil(H / Ts); % Horizon steps

            [nx, nu] = size(mpc.B);

            % Targets (Ignore this before Todo 3.2)
            x_ref = sdpvar(nx, 1);
            u_ref = sdpvar(nu, 1);

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
            
            Q = 3 * eye(nx);
            Q(1,1)=70;
            Q(2, 2) = 150;
            Q(4, 4) = 10;
            R = eye(nu);
%             Q = 3 * eye(nx);
%             Q(4, 4) = 6;
%             R = eye(nu);


            % ---- compute terminal set ----
            % terminal weight
            Q_f = Q;
            R_f = R;

            % terminal controller
            [K_f, S, e] = dlqr(mpc.A, mpc.B, Q_f, R_f);
            K_f = -K_f;

            % state and input constraints
            F_x = [1, 0, 0, 0;
                -1, 0, 0, 0;
                0, 1, 0, 0;
                0, -1, 0, 0;
                0, 0, 1, 0;
                0, 0, -1, 0;
                0, 0, 0, 1;
                0, 0, 0, -1];
            F_u = [1; -1];
            F = [F_x; F_u * K_f];

            f_x = [inf; inf; 0.0873; 0.0873; inf; inf; inf; inf];
            f_u = [0.26; 0.26];
            f = [f_x; f_u];

            % compute terminal invariant set
            X_inv = Polyhedron(F, f);
            current_set = X_inv;
            A_closed_loop = mpc.A + mpc.B * K_f;

            iteration = 0;
            while (1)
                fprintf("iteration %d \n", iteration)
                previous_set = current_set;
                F = current_set.A;
                f = current_set.b;
                preset = Polyhedron([F * A_closed_loop], f);
                current_set = intersect(previous_set, preset).minHRep;

                if (current_set == previous_set)
                    break;
                end

                iteration = iteration + 1;
            end

            X_inv = current_set;

            % input constraints for time 1
            con = [con, -0.26 * ones(nu, 1) <= U(:, 1) <= 0.26 * ones(nu, 1)];

            % cost function for time 1
            obj = obj + (U(:, 1) - u_ref)' * R * (U(:, 1) - u_ref);

            % time 2 to N-1
            for i = 2:1:N - 1
                % input constraints
                con = [con, -0.26 * ones(nu, 1) <= U(:, i) <= 0.26 * ones(nu, 1)];
                % state constraints
                con = [con, [-inf; -0.0873; -inf; -inf] <= X(:, i) <= [inf; 0.0873; inf; inf]];
                % state-space model constraints
                con = [con, X(:, i) == mpc.A * X(:, i - 1) + mpc.B * U(:, i - 1)];
                % object
                obj = obj + (X(:, i) - x_ref)' * Q * (X(:, i) - x_ref) + (U(:, i) - u_ref)' * R * (U(:, i) - u_ref);
            end

            % final state cost
            obj = obj + (X(:, N) - x_ref)' * S * (X(:, N) - x_ref);

            % final state constraints
            F_f = X_inv.A;
            f_f = X_inv.b;
            con = [con, F_f * X(:, N) <= f_f + F_f * x_ref];

            % YOUR CODE HERE YOUR CODE HERE YOUR CODE HERE YOUR CODE HERE
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

            % Return YALMIP optimizer object
            ctrl_opti = optimizer(con, obj, sdpsettings('solver', 'gurobi'), ...
            {X(:, 1), x_ref, u_ref}, U(:, 1));
        end

        % Design a YALMIP optimizer object that takes a position reference
        % and returns a feasible steady-state state and input (xs, us)
        function target_opti = setup_steady_state_target(mpc)

            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % INPUTS
            %   ref    - reference to track
            % OUTPUTS
            %   xs, us - steady-state target
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

            nx = size(mpc.A, 1);

            % Steady-state targets
            xs = sdpvar(nx, 1);
            us = sdpvar;

            % Reference position (Ignore this before Todo 3.2)
            ref = sdpvar;

            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % YOUR CODE HERE YOUR CODE HERE YOUR CODE HERE YOUR CODE HERE
            % You can use the matrices mpc.A, mpc.B, mpc.C and mpc.D
            obj = 0;
            con = [-0.26 <= us <= 0.26,
                xs == mpc.A * xs + mpc.B * us,
                ref == mpc.C * xs + mpc.D * us];
            obj = obj + us' * us;

            % YOUR CODE HERE YOUR CODE HERE YOUR CODE HERE YOUR CODE HERE
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

            % Compute the steady-state target
            target_opti = optimizer(con, obj, sdpsettings('solver', 'gurobi'), ref, {xs, us});
        end

    end

end
