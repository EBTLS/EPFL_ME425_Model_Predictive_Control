% ME-425 : Model Predictive Control
% Exercise sheet 4
% 
% Exercise 2, Solution
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

nx=2;
nu=1;

%% Terminal Controller
% That is an infinite LQR Controller

Q_f=Q;
R_f=R;

[K_f,S,e]=dlqr(A,B,Q_f,R_f);
K_f=-K_f;


%% Define Terminal Set

figure('Name','Terminal Set')

G_f=G*K_f;
F=[H;G_f];
f=[h;g];

X=Polyhedron(F,f);
X.plot('Color','g','alpha',0.1);
X.volume

iteration=0;

current_set=X;

A_closed_loop=A+B*K_f;

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

%% computing Matrices for MPC

x=sdpvar(repmat(nx,1,N),repmat(1,1,N));
u=sdpvar(repmat(nu,1,N-1),repmat(1,1,N-1));

con=[];
obj=0;

con=[con,x{2}==A*x{1}+B*u{1}];
% con=[con,H*x{1}<=h];
con=[con,G*u{1}<=g];

obj=obj+x{1}'*Q*x{1}+u{1}'*R*u{1};

for i=2:N-1

    con=[con,x{i+1}==A*x{i}+B*u{i}];
    con=[con,H*x{i}<=h];
    con=[con,G*u{i}<=g];
    obj=obj+x{i}'*Q*x{i}+u{i}'*R*u{i};

end

F_f=current_set.A;
f_f=current_set.b;

con=[con,F_f*x{N}<=f_f];
obj=obj+x{N}'*S*x{N};

ops = sdpsettings('solver','sedumi');
ctrl=optimizer(con,obj,ops,x{1},u{1});

sol.x=x_0;
sol.u=[];

i=1;

while norm(sol.x(:,end))>1e-3

    
    [uopt,isfeasible]=ctrl{sol.x(:,end)};
    
    if (isfeasible~=1)
       
        sol.u=[sol.u,uopt];
        temp_x=A*sol.x(:,end)+B*uopt;
        sol.x=[sol.x,temp_x];
        
    else    
        
        display("controller solution is not feasible");
        
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


