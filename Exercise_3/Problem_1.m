%%%
% MPC-425, Exercise 3 Probelm 1
%%%
%% Global Settings
alpha=pi/6;
beta=0.8;
A=[cos(alpha), sin(alpha);
    -sin(alpha), cos(alpha)]*beta;
H=[cos(pi/3),sin(pi/3);
    -cos(pi/3),-sin(pi/3);
    sin(pi/3),-cos(pi/3);
    -sin(pi/3),cos(pi/3)];
h=[2,1,2,5]';
%% Task 1: Compute the largest invariant set

figure('Name','task 1')

X=Polyhedron(H,h);
X.plot('Color','g','alpha',0.1)

current_set=Polyhedron(H,h);
previous_set=Polyhedron(H,h);

iteration=0;

while(1)
    
    hold on;
%     previous_set.plot('Color','b','alpha',0.2)
    previous_set=current_set;
    F=current_set.A;
    f=current_set.b;
    preset=Polyhedron([F*A],f);
    current_set=intersect(previous_set,preset);
    current_set.plot('Color','b','alpha',0.2)
    previous_set.plot('Color','r','alpha',0.2)
  
    if (current_set==previous_set)
        break;
    end
    
    iteration=iteration+1;

end

hold on;

current_set.plot('Colar','r','alpha',0.3)
title("Maximal Invariante Set")


%% Task 2: Plot state trajectories of the system

figure('Name', 'trajectories')

X.plot('Color','g','alpha',0.1)
hold on;
current_set.plot('Colar','r','alpha',0.3)

ss_model=ss(A,[0;0],eye(size(A)),[0;0],1);
step=50;
i=0;
while i<20

    x_0=rand(2,1)*8-4;
    
    if any(H*x_0>h)
        continue;
    else
        x_0

        i=i+1;
        trajectory=lsim(ss_model,zeros(1,step),[0:1:step-1],x_0)
        hdl = plot(trajectory(:,1),trajectory(:,2),'--k')
        hold on;
        set(hdl,'color','k')
        plot(x_0(1),x_0(2),'-o','markersize',10);
        grid on;
    end
        
end


%% Task 3: Plot the maximum invariant set