
parameter_setting;
% 
% x_trajectory=[];
% 
% scatter (x_0(1),x_0(2));
% 
% for (i=1:1:100)
%     
%     [u,x,K,v]=FiniteRecedingHorizon(A,B,C,Q,R,5,x_0);
% 
%     v=v+x_0'*Q*x_0
% 
%     x=reshape(x,[2,length(x)/2]);
%     x_previous=x_0;
%     x_0=x(:,1);
%    
%     scatter (x_0(1),x_0(2), 'b');
%     
%     plot(x_0,x_previous, 'r')
% 
%     grid on;
% 
%     hold on;
%     
%     input("wait")
%     
% end

[u,x,K,v]=FiniteRecedingHorizon(A,B,C,Q,R,10,x_0);
v=v+x_0'*Q*x_0
x=reshape(x,[2,length(x)/2]);
x=[x_0,x];
plot(x(1,:),x(2,:),'b')
grid on;
hold on;

% [u,x,K,v]=FiniteRecedingHorizon(A,B,C,Q,R,25,x_0);
% v=v+x_0'*Q*x_0
% x=reshape(x,[2,length(x)/2]);
% x=[x_0,x];
% plot(x(1,:),x(2,:),'r')
% grid on;
% hold on;
% 
% [u,x,K,v]=FiniteRecedingHorizon(A,B,C,Q,R,30,x_0);
% v=v+x_0'*Q*x_0
% x=reshape(x,[2,length(x)/2]);
% x=[x_0,x];
% plot(x(1,:),x(2,:),'g')
% grid on;