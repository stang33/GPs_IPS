function [dy] =RHSfn_2nd_ncf_NS1(t,y,N,kernel,ncforce,type)
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
% Authored by Jinchao Feng and Sui Tang


d = size(y,1)/(2*N);
y = reshape(y,d,[]); % [x_1,..., x_N,v_1 v_2..., v_N]  d x 2N
aa= zeros(d,2*N);

switch type
    
    case 'E'

        for i=1:N
            temp_position      = y(:,1:N)- repmat(y(:,i),1,N);  % d x N-1 [x_1-x_i,\cdots,x_N-x_i]
            DD        = sqrt( sum(temp_position.^2,1) ); % 1 x N
            aa(:,i) = y(:,i+N);
            aa(:,i+N)   = ncforce(y(:,i+N),d) + temp_position * kernel(DD)'/sum(kernel(DD));   % 1 x d
        end
        dy = reshape(aa,[],1); % dN x 1  
        
        
     case 'A'
         
        for i=1:N
            temp_position      = y(:,1:N)- repmat(y(:,i),1,N);  % d x N [x_1-x_i,\cdot,x_N-x_i]
            DD        = sqrt( sum(temp_position.^2,1) ); % 1 x N
            temp_velocity      = y(:,N+1:2*N)- repmat(y(:,i+N),1,N);  % d x N [v_1-v_i,\cdot, v_N-v_i]
            aa(:,i) = y(:,i+N);
            aa(:,i+N)   = ncforce(y(:,i+N),d) + temp_velocity * kernel(DD)'/N;   % 1 x d
        end
        dy = reshape(aa,[],1); % dN x 1  
        
end
end