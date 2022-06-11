%% compare the performance 
% 
% 

sysInfo.N =50;
sysInfo.d =2;

load('traj_real.mat');
xpath = xpath;
load('traj_hat_0.mat');
traj_hat_0 =traj_hat;
load('traj_hat_100.mat');
traj_hat_100 =traj_hat;

% the true flocking direction vectors and flocking socre 
[V,E]= flockingscore(xpath,sysInfo);

% the predicted one by the learned system with training step 100
[V1,E1]= flockingscore(traj_hat_100,sysInfo);
score1 = (diag(V'*V1));
% relative error in prediction of Y
err1 = norm(traj_hat_100-xpath)./norm(xpath);
% relative error in prediction of X
err_pos1 =  norm(traj_hat_100(1:sysInfo.d*sysInfo.N,:)-xpath(1:sysInfo.d*sysInfo.N,:))./norm(xpath(1:sysInfo.d*sysInfo.N,:));

fprintf('\n\n-----leanred system with training-----\n\n')

fprintf('Rel Position Err is %f\n', err_pos1);
fprintf('Rel State Err is %f\n', err1);



% the predicted one by the learned system without training step 
[V2,E2]= flockingscore(traj_hat_0,sysInfo);
score2 = (diag(V'*V1));

% relative error in prediction of Y
err2 = norm(traj_hat_0-xpath)./norm(xpath);
% relative error in prediction of X
err_pos2 =  norm(traj_hat_0(1:sysInfo.d*sysInfo.N,:)-xpath(1:sysInfo.d*sysInfo.N,:))./norm(xpath(1:sysInfo.d*sysInfo.N,:));


fprintf('\n\n-----leanred system without training-----\n\n')

fprintf('Rel Position Err is %f\n', err_pos2);
fprintf('Rel State Err is %f\n', err2);

function [V,E] = flockingscore(data,sysInfo)
% compute the flocking direction and score of the data 

N = sysInfo.N;
d = sysInfo.d;
L = size(data,2);
V=[];
E=[];
for l=1:L
v=reshape(data(d*N+1:end,l),d,N);
v= normc(v);% normalize the velocity vector to be 1
[w,D]=eigs(v*v',2);% find the top eigenvector
flag = sum(w(:,1)'*v(:,1));
if flag >0
V=[V w(:,1)];
E= [E; sum(w(:,1)'*v./N)];
else
    V=[V -w(:,1)];
    E= [E; sum(-w(:,1)'*v./N)]
end

end


end

