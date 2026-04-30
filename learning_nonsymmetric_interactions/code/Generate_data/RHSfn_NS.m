function [dx] =RHSfn_NS(t,x,N,kernel)
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
% Authored by Jinchao Feng and Sui Tang

d = size(x,1)/N;
y = reshape(x,d,N); % [x_1   x_N]  d x N
aa= zeros(d,N);

for i=1:N
    temp      = y(:,[1:i-1 i+1:N])- repmat(y(:,i),1,N-1);  % d x N-1 [x_1-x_i,\cdots,x_{i-1}-x_i,x_{i+1}-x_i,\cdot,x_N-x_i]
    DD        = sqrt( sum(temp.^2,1) ); % 1 x N-1
    DD(DD==0) = 1;  %  this step is to avoid 0 x inf= NAN

    if norm(kernel(DD))==0
        K =zeros(N-1,1);
    else
        K = kernel(DD)'/sum(kernel(DD));
    end
    
    aa(:,i)   = temp*K;   % d x 1
end
dx = reshape(aa,[],1); % dN x 1
end

