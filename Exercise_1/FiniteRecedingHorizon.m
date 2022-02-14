function [u,x,K,v]=FiniteRecedingHorizon(A,B,C,Q,R,N,x_0)
    % Function to calculate optimal solution for a LQR with Finite Receding Horizon
    % A,B,C: system model
    % Q,R: QR matrices
    % x0: initial state
    % N: horizon


    B_math_base=eye(N);
    B_math=kron(B_math_base,B);

    A_math_base_1=diag(ones(1,N-1),-1);
    A_math_base_2=-diag(ones(1,N));
    A_math=kron(A_math_base_1,A)+kron(A_math_base_2,eye(size(A,1)));

    C_math_base=zeros(N,1);
    C_math_base(1)=1;
    C_math=kron(C_math_base,-A);

    Q_math_base=eye(N);
    R_math_base=eye(N);

    Q_math=kron(Q_math_base,Q);
    R_math=kron(R_math_base,R);

    F=-inv(A_math)*B_math;
    G=inv(A_math)*C_math;

    K=-inv(R_math+F'*Q_math*F)*F'*Q_math*G;
    u=K*x_0;
    x=inv(A_math)*(C_math*x_0-B_math*u);
    v=x'*Q_math*x+u'*R_math*u;

end