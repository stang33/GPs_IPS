function traj_hat = construct_traj(learnInfo,test_time_vec,ICs)

% choose from existing sets of stable Initial Conditions
close all;

N = learnInfo.N;
d = learnInfo.d;
order = learnInfo.order;
dN = d*N*order;

% myODE1     = @(t,x)predict_ode(x,learnInfo);
myODE1     = @(t,x)predict_ode_interp(x,learnInfo);



L=length(test_time_vec);
M = size(ICs, 2);
traj_hat =zeros(dN,L,M);
solverInfo.time_span    = [0, test_time_vec(end)];                                                       % put it into the time_span vector, always starting from 0
solverInfo.option = odeset('RelTol',1e-5,'AbsTol',1e-6);
for m = 1 : M
 
 sol_hat = ode45(myODE1,solverInfo.time_span,ICs(:,m),solverInfo.option); % solu from adaptive solver
 traj_hat(:,:,m)   = deval(sol_hat,test_time_vec);

end

end
