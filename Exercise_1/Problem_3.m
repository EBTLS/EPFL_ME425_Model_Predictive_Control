
parameter_setting;

[K,S,e]=dlqr(A, B, Q, R);

x_current=x_0;
v_infinite_controller=0;

for i=1:1:1000
    
    u=-K*x_current;
    x_next=A*x_current+B*u;
    temp_v=x_current'*Q*x_current+u'*R*u;
    v_infinite_controller=v_infinite_controller+temp_v;
    x_current=x_next;

end

v_infinite_controller

x_current=x_0;
v_finite_controller=0;

for i=1:1:1000

    [u,x,K,v]=FiniteRecedingHorizon(A,B,C,Q,R,10,x_current);
    
    u=reshape(u,[2,length(u)/2]);
    u_current=u(:,1);
    x=reshape(x,[2,length(x)/2]);
    x_next=x(:,1);
    temp_v=x_current'*Q*x_current+u_current'*R*u_current;
    
    v_finite_controller=v_finite_controller+temp_v; 
    x_current=x_next;
    
end

v_finite_controller