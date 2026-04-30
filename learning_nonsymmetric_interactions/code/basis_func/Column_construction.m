function [temp] = Column_construction(ker,x,dotx,N)
%  Construct a column of the basis matrix C 
%  Input   - ker: function
%          - x: position data   Nd x L
%  Output - dNL x 1 matrix

L=size(x,2); % time steps
temp=[]; 
for j=1:L
    RHS  =  RHSfn_basis(ker,reshape(x(:,j),N,[]),reshape(dotx(:,j),N,[]));% Nd*1
    temp = [temp;RHS];
end
    
end
