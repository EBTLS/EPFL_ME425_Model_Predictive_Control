
parameter_setting

Q=C'*C+0.001*eye(2);

R=0.001;

x_0=[10;10];

[u,x,y,K,v]=DiscreteBellmanRecursion(A,B,C,Q,R,x_0,10);

v

% [K,S,e]=dlqr(A, B, Q, R, 10)

