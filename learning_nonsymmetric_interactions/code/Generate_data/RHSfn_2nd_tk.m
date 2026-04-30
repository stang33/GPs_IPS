function [dy] =RHSfn_2nd_tk(t,y,N,kernele, kernela, ncforce)
% The RHS function of the ODE
% Input
%   y=[x v]   - the state vector dN x 1, d=2; 
%        [x_1(1), ..., x_1(d), ..., x_N(1),..., x_N(d) v_1(1),...,v_1(d),...,v_N(1),...,v_N(d)]^T
%   t   - time 
% 
% Output
%   dy  -  right hand side of the ode  dN x 1  
%          [x_1'(1),..., x_1'(d), ... , x_N'(1),..., x_N'(d), v_1'(1),..., v_1'(d),...]^T
%
% (c) Sui Tang, Jinchao Feng

d = size(y,1)/(2*N);
y = reshape(y,d,[]); % [x_1,..., x_N,v_1 v_2..., v_N]  d x N
aa= zeros(d,2*N);

for i=1:N
    temp_position      = y(:,1:N)- repmat(y(:,i),1,N);  % d x N [x_1-x_i,\cdots, x_N-x_i]
    DD        = sqrt( sum(temp_position.^2,1) ); % 1 x N
    %DD(DD==0) = 1;  %  this step is to avoid 0 x inf= NAN
    temp_velocity      = y(:,N+1:2*N)- repmat(y(:,i+N),1,N);  % d x N [v_1-v_i,\cdots, v_N-v_i]
    %SS        = sqrt( sum(temp_velocity.^2,1) ); % 1 x N
    aa(:,i) = y(:,i+N);
    aa(:,i+N)   = ncforce(y(:,i+N),d) + temp_velocity* kernela(DD)'/N + temp_position* kernele(DD)'/N;   % 1 x d
end

dy = reshape(aa,[],1); % dN x 1  
end