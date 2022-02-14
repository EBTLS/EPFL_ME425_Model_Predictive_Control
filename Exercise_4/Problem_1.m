% ME-425 : Model Predictive Control
% Exercise sheet 4
% 
% Exercise 1, Solution
%
%% Global Settings

A=[0.9752, 1.4544;
    -0.0327, 0.9315];

B=[0.0248;
    0.0327];

Q=10*eye(2);

R=eye(1);

H=[1,0;
   -1,0;
   0,1;
   0,-1];

h=[5;
   5;
   0.2;
   0.2];

G=[1;
   -1];

g=[1.75;
   1.75];

N=10;

x_0=[3;0];

%% Terminal Controller
% That is an infinite LQR Controller

Q_f=Q;
R_f=R;

[K_f,S,e]=dlqr(A,B,Q_f,R_f);
K_f=-K_f;


%% Define Terminal Set

figure('Name','Terminal Set')

G_f=G*K_f;
F=[H;G_f]
f=[h;g]

X=Polyhedron(F,f);
X.plot('Color','g','alpha',0.1);
X.volume

iteration=0;

current_set=X;

A_closed_loop=A+B*K_f

while (1)
   
    hold on;
    
    previous_set=current_set;
    F=current_set.A;
    f=current_set.b;
    preset=Polyhedron([F*A_closed_loop],f);
    current_set=intersect(previous_set,preset);
    current_set.plot('Color','b','alpha',0.2);
    
    if (current_set==previous_set)
    
        break;
        
    end
    
end

%% Validate Terminal Set

figure('Name','MPT result')

sys = LTISystem('A',A,'B',B);
sys.x.max=[5;0.2];
sys.x.min=[-5;-0.2];
sys.u.max=[1.75];
sys.u.min=[-1.75];
sys.x.penalty = QuadFunction(Q);
sys.u.penalty = QuadFunction(R);

T_K_mpt=sys.LQRGain
T_penalty_mpt=sys.LQRPenalty.weight
T_set_mpt=sys.LQRSet
T_set_mpt.plot;


%% simulate the system
sol.x=[x_0];
sol.u=[];
i=1;
while norm(sol.x(:,i))>1e-3
    
    [H_MPC,h_MPC,G_MPC,g_MPC,T_MPC,t_MPC]=ComputeMatrices(A,B,Q,R,S,H,h,G,g,N,sol.x(:,i));
    
    [z,fval,flag]=quadprog(H_MPC,h_MPC,G_MPC,g_MPC,T_MPC,t_MPC);
    
    if (flag)
        
        sol.u(:,i)=z(2*N+1);
        sol.x(:,i+1)=A*sol.x(:,i)+B*z(2*N+1);

    else
        
        disp("not optimal solution!!")
        
    end
    
    i=i+1;

end

subplot(3,1,1)
plot([1:1:i],sol.x(1,:));
hold on;
plot([1:1:i],-5*ones(1,i),'Color','r');
hold on;
plot([1:1:i],5*ones(1,i),'Color','r');
grid on;

subplot(3,1,2)
plot([1:1:i],sol.x(2,:));
hold on;
plot([1:1:i],-0.2*ones(1,i),'Color','r');
hold on;
plot([1:1:i],0.2*ones(1,i),'Color','r');
grid on;

subplot(3,1,3)
plot([0:1:i-2],sol.u);
hold on;
plot([0:1:i-2],-1.75*ones(1,i-1));
hold on;
plot([0:1:i-2],1.75*ones(1,i-1));
grid on;

%% function for computing Matrices for MPC

function [H_MPC,h_MPC,G_MPC,g_MPC,T_MPC,t_MPC]=ComputeMatrices(A,B,Q,R,S,H,h,G,g,N,x_0)

    H_MPC_1=kron(eye(N-1),Q);
    H_MPC_2=kron(eye(N),R);
    H_MPC=blkdiag(H_MPC_1,S,H_MPC_2);

    h_MPC=zeros(size(H_MPC,1),1);

    G_MPC_1=kron(eye(N),H);
    G_MPC_2=kron(eye(N),G);
    G_MPC=blkdiag(G_MPC_1,G_MPC_2);

    g_MPC_1=kron(ones(N,1),h);
    g_MPC_2=kron(ones(N,1),g);
    g_MPC=[g_MPC_1;g_MPC_2];

    T_MPC_1=kron(diag(ones(1,N-1),-1),A);
    T_MPC_2=kron(eye(N),B);
    T_MPC_3=zeros(N,size(T_MPC_1,2));
    T_MPC_4=eye(N);

    T_MPC=[T_MPC_1,T_MPC_2;
           T_MPC_3,T_MPC_4];

    T_MPC=eye(size(T_MPC))-T_MPC;

    t_MPC_1=kron([1;zeros(N-1,1)],A);
    t_MPC_2=zeros(N,size(A,2));
    t_MPC=[t_MPC_1;t_MPC_2]*x_0;

end



