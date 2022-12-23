function [R, Aplus] = NystNoK(K, rank)

%Initialize
Aplus = zeros(rank, rank);
A = K(1:rank, 1:rank);

[U,S,~] = svd(A,"vector");

for i = 1 : rank
    
    singval = S(i);
    col = U(:,i);
    Aplus = Aplus + (1.0 / singval) * (col * col.');

end

Rprime = K(1:rank,:);
R = Rprime';

end