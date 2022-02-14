function [u,x,y,K,v]=DiscreteBellmanRecursion(A,B,C,Q,R,x_0,N)
    % Function to calculate optimal solution for a LQR with DiscrerteBellmanRecursion
    % A,B,C: system model
    % Q,R: QR matrices
    % x0: initial state
    % N: horizon


    x=x_0;
    K=[];
    u=[];
    y=[C*x_0];

    H_current=Q;

    for (i=N-1:-1:0)

        H_next=H_current;

        K_current=-inv(R+B'*H_next*B)*B'*H_next*A;
        H_current=Q+K_current'*R*K_current+(A+B*K_current)'*H_next*(A+B*K_current);

        K=[K_current,K];
    
    end

    for (i=0:1:N-1)
    
        x_current=x(:,i+1);
        K_current=K(:,[i*length(x_current)+1:(i+1)*length(x_current)]);
        u_current=K_current*x_current;
        u=[u,u_current];
        x_next=A*x_current+B*u_current;
        y_next=C*x_next;
        x=[x,x_next];
        y=[y,y_next];

    end

    v=x_0'*H_current*x_0;

end