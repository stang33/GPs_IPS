function [A] = Recfns(basis_funs,a,I)
% 

n=size(basis_funs,2); % number of basis functions 
A=zeros(size(I));

    for i=1:n
        temp=basis_funs{i};
        A=A+a(i)*temp(I);
    end

end
