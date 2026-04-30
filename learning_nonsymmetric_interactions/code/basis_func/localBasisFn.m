
function basis_funs = localBasisFn(edges,degree)
% piecewise polynomials with degree on the intervals specified by edges
% 
n = length(edges)-1;
if degree ==0
    basis_funs = cell(1, n); %  interval>Rmax is ingnored
    for i=1:n-1
        basis_funs{i}=@(x)(x<edges(i+1)).*(x>=edges(i));
    end
    basis_funs{n}=@(x)(x<=edges(n+1)).*(x>=edges(n));
elseif degree ==1
    basis_funs        = cell(1, 2*n); 
    for i=1:n
        basis_funs{2*i-1}=@(x)(x<edges(i+1)).*(x>=edges(i));  % on each interval, we have two basis:
        basis_funs{2*i}=@(x)x.*(x<edges(i+1)).*(x>=edges(i)); %  constant function+linear function
    end
end
end