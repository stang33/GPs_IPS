function [traj_hat, traj_hat_std] = construct_traj_mcuq(learnInfo,sysInfo,test_time_vec,ICs)

% choose from existing sets of stable Initial Conditions
close all;

N = learnInfo.N;
d = learnInfo.d;
order = learnInfo.order;
dN = d*N*order;

% myODE1     = @(t,x)predict_ode(x,learnInfo);
myODE1     = @(t,x)predict_ode_mcuq(x,learnInfo,sysInfo);



L=length(test_time_vec);
M = size(ICs, 2);
traj_hat =zeros(dN,L,M);
traj_hat_std             = zeros (dN,L,size(ICs, 2)); % hat traj

solverInfo.time_span    = [0, test_time_vec(end)];                                                       % put it into the time_span vector, always starting from 0
solverInfo.option = odeset('RelTol',1e-5,'AbsTol',1e-6);

for m = 1 : M
 
 num_mc = 10;
 sol_hat = cell(1, num_mc);
 traj_hat_mc = cell(1, num_mc);

 for k = 1:num_mc
     
     sol_hat{k} = ode45(myODE1,solverInfo.time_span,ICs(:,m),solverInfo.option); % solu from adaptive solver
     traj_hat_mc{k} = deval(sol_hat{k},test_time_vec);
 end
 
 traj_hat_mc0            = zeros (dN,L,num_mc); % hat traj
 for k = 1:num_mc
     traj_hat_mc0(:,:,k) = traj_hat_mc{k};
 end
 
 traj_hat(:,:,m)   = mean(traj_hat_mc0,3);
 traj_hat_std(:,:,m)   = std(traj_hat_mc0,0,3);

end

end