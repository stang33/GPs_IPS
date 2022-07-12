function [dx] =RHSfn(t,x,N,kernel)
% The RHS function of the ODE
% Input
%   x   - the state vector dN x 1, d=2; 
%        [x_1(1), ..., x_1(d), ... , x_N(1),..., x_N(d)]^T
%   t   - time 
% 
% Output
%   dx  -  right hand side of the ode  dN x 1  
%          [x_1'(1),..., x_1'(d), ... , x_N'(1),..., x_N'(d)]^T
%
% (c) Sui Tang, Fei Lu 

d = size(x,1)/N;
y = reshape(x,d,N); % [x_1   x_N]  d x N
aa= zeros(d,N);

for i=1:N
    temp      = y- repmat(y(:,i),1,N);  % d x N [x_1-x_i,\cdots, x_N-x_i]
    DD        = sqrt( sum(temp.^2,1) ); % 1 x N
    DD(DD==0) = 1;  %  this step is to avoid 0 x inf= NAN
    aa(:,i)   = temp* kernel(DD)';   % 1 x d
end
aa = aa/N;          % [(x_1'), \cdots, (x_N')]; d x N
dx = reshape(aa,[],1); % dN x 1  
end

