function [C] =BasisMatrix_construction(basis_funs,x,dotx,N)
% construct the basis tensor matrix C in the least square
% Input   - basis_funs: a cell array of basis functions handle 
%         - x: position data  Nd*L
%         - N: number of agents
%
% Output  - tensor matrix C dNL*n, n=# of basis
%
%idea C(:,:,i)=?


n=size(basis_funs,2); % number of basis functions
L=size(x,2);% number of time instances
d=size(x,1)/N;
C=zeros(d*N*L,n);


for i=1:n
    C(:,i)=Column_construction(basis_funs{i},x,dotx,N);
end
    
end
