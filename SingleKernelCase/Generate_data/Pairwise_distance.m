function [c] = Pairwise_distance(x,N)
% to construct the measure rho_T, this function compute the pairwise
% distance of agents at different time instances
%
% Input  - x: position data Nd*L, where L is the number of time
%                 instances
% 
% Output     - pair-distance ((N-1)N/2)L*1  

% (c) XXXX

L = size(x,2); % total time instances
c = zeros((N-1)*N/2,L); 

for t=1:L
     temp  = reshape(x(:,t),[],N);% d x N 
     c(:,t)  = pdist(temp');% compute the pairwise distance between the x_i and x_j
end
c = reshape(c,[],1); % a column vector that constains nonzero pairwise distance at all time snapshots
end

