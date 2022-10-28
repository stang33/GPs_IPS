function [M,P] = group_polarisation(data,sysInfo)
% compute the flocking score of the data 

N = sysInfo.N;
d = sysInfo.d;
L = size(data,2);
M=[];
P=[];
for l=1:L
    v=reshape(data(d*N+1:end,l),d,N);
    v= normc(v);% normalize the velocity vector to be 1
    vm = mean(v,2);

    M=[M vm];
    P= [P norm(vm)];

end

end

