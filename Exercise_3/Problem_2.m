%%%
% ME-425 Exercise 3 Problem 2
%%%

%% Global Parameter Settings

alpha=pi/6;
beta=0.8;
A=[cos(alpha),sin(alpha);
    -sin(alpha),cos(alpha)]*beta;
B=[0.5;0.5];

H=[cos(pi/3),sin(pi/3);
    -cos(pi/3),-sin(pi/3);
    sin(pi/3),-cos(pi/3);
    -sin(pi/3),cos(pi/3)];
h=[2,1,2,5]';

u_lb=-0.5;
u_ub=0.5;

%% Input Constraints Transformation

H_u=[1;-1];
h_u=[0.5;
    0.5];

H_total=[H,zeros(size(H,1),size(H_u,2));
         zeros(size(H_u,1),size(H,2)),H_u];
h_total=[h;
        h_u];
    

%% Task 1

figure('Name','maximum controlled invariante set')

X=Polyhedron(H_total,h_total);
current_set=Polyhedron(H,h);
current_set.plot('Color','b','alpha',0.1);
hold on;

iteration=0;

while(1)

    previous_set=current_set;
    
    F=current_set.A;
    f=current_set.b;
    
    H_preset=[F*A,F*B;
              zeros(size(H_u,1),size(F*A,2)),H_u];
    h_preset=[f;h_u];
          
    preset=Polyhedron(H_preset,h_preset);
    current_set=intersect(previous_set,projection(preset,1:2));
%     current_set.plot('Color','g','alpha',0.3);
    hold on;
    
    iteration=iteration+1;
    if (previous_set==current_set)
        
        current_set.plot('Color','r','alpha',0.5)
    
        break;
    
    end
   
end

% figure("Name","Projection")
% X_projection=projection(X,[1:2]);
% X_projection.plot('Color','g','alpha',0.1);
% hold on;
% current_set.plot('Color','r','alpha',0.5)


%% Task 2: optimal LQR controller

figure('Name',"LQR Controller")

Q=eye(2)
R=1;

[K,S,e]=dlqr(A,B,Q,R);
A_lqr=A-B*K;

H_lqr=[H;H_u*K];
h_lqr=[h;h_u];

current_set=Polyhedron(H_lqr,h_lqr);
current_set.plot('Color','b','alpha',0.3);
hold on;

iteration=0;

while(1)
    
    previous_set=current_set;
    
    F=current_set.A;
    f=current_set.b;
    
    preset=Polyhedron(F*A_lqr,f);
    current_set=intersect(preset,previous_set);
    
    iteration=iteration+1;
    
    if (current_set==previous_set)
        current_set.plot('Color','r','alpha',0.5);
        break;
    end

end





