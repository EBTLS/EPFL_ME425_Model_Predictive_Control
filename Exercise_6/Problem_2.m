%% ME425 Model Predictive Control
% Exercise 6 Problem 2
%%%%%%%%%%%%%
%% Clear
clear all
close all
%% System Parameters
A = [0.9752, 1.4544;
    -0.0327, 0.9315];

B = [0.0248;
    0.0327];

xmin = [-5; -0.2];

xmax = [5; 0.2];

umin = [-1.75];
umax = [1.75];

N = 10;
n_steps=50;

Q = 10 * eye(2);
R = [1];

%% Define MPC Problem

sys = LTISystem('A', A, 'B', B);
sys.x.max = xmax;
sys.x.min = xmin;
sys.u.max = umax;
sys.u.min = umin;

sys.x.penalty = QuadFunction(Q);
sys.u.penalty = QuadFunction(R);

% Extract Terminal state and terminal set
Qf = sys.LQRPenalty.weight;
Xf = sys.LQRSet;

% Define the controller
controller = MPCController(sys, N);

%% Generate explicit MPC
empc = controller.toExplicit();
empc.exportToC('output','directory')

empc.feedback.fplot();

%% simulate the result
x0 = [3; 0];

solution.x_hist=zeros(1,n_steps-1);
solution.x_hist=x0;
solution.u_hist=zeros(1,n_steps-1);

for i=1:1:n_steps-1
    
    [uopt,isfeasible]=empc.evaluate(solution.x_hist(:,end));
    
    if (isfeasible)
        
        solution.u_hist(:,i)=uopt;
        solution.x_hist(:,i+1)=A*solution.x_hist(:,i)+B*uopt;
       
    else
        fprintf("step %d solution not feasible", isfeasible);
    end
    
end

figure
subplot(3,1,1)
plot([1:1:n_steps],solution.x_hist(1,:));
hold on;
plot([1:1:n_steps],-5*ones(1,n_steps),'Color','r');
hold on;
plot([1:1:n_steps],5*ones(1,n_steps),'Color','r');
grid on;

subplot(3,1,2)
plot([1:1:n_steps],solution.x_hist(2,:));
hold on;
plot([1:1:n_steps],-0.2*ones(1,n_steps),'Color','r');
hold on;
plot([1:1:n_steps],0.2*ones(1,n_steps),'Color','r');
grid on;

subplot(3,1,3)
plot([0:1:n_steps-2],solution.u_hist);
hold on;
plot([0:1:n_steps-2],-1.75*ones(1,n_steps-1));
hold on;
plot([0:1:n_steps-2],1.75*ones(1, n_steps-1));
grid on;

