function aa = RHSfn_basis(kerFn,x,xdot)
% The right hand side of the self-organized dynamics
%      x_i'=sum_{j=1}^N a(|x_j-x_i|)(x_j-x_i 
% Output:
%      aa       - the RHS of equations Nd x 1
% Input:
%      x        - the state array   N x d
%      kerFn    - the interaction kernel
%      N        - number of agents

%  % test
% x   = randn(10,2);
% ker = @(r) sin(r).* (r<1) + exp(-r).*(r>=1);

[N,d] = size(x);

aa=zeros(N*d,1);
for i=1:N
    temp   = x- repmat(x(i,:),N,1);  % N*d
    DD     = sqrt( sum(temp.^2,2) ); % N*1
   % DD(DD==0)=1;
    aa((i-1)*d+1:i*d,1)= (kerFn(DD')* (xdot(i,:)-temp))';   % 1*N x Nxd 
end

return