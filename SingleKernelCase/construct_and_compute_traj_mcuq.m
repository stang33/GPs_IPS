function result = construct_and_compute_traj_mcuq(sysInfo, obsInfo,  solverInfo,learnInfo,ICs)
% function result = construct_and_compute_traj(, sys_info, syshat_info, ICs)

% (c) Sui Tang 

% choose from existing sets of stable Initial Conditions


%% basic setting of the system
N         = sysInfo.N;         % number of agents   
d         = sysInfo.d; 
order = sysInfo.ode_order;
dN = N*d*sysInfo.ode_order;
myODE =  @(t,x) RHSfn(t,x,N,sysInfo.phi{1});

if order ==1&&strcmp(learnInfo.name,'ODS')
    myODE =  @(t,x) RHSfn_c(t,x,N,sysInfo.phi{1},sysInfo.phi{2});
end


if order ==2&&strcmp(learnInfo.name,'CS')
    myODE =  @(t,x) RHSfn_2nd(t,x,N,sysInfo.phi{1});
end

if order ==2&&strcmp(learnInfo.name,'CSF')
    myODE =  @(t,y) RHSfn_2nd_ncf(t,y,N,sysInfo.phi{1},sysInfo.phi{2},sysInfo.phi_type);
end

if order ==2&&strcmp(learnInfo.name,'FM')
 myODE     = @(t,y) RHSfn_2nd_ncf(t,y,N,sysInfo.phi{1},sysInfo.phi{2},sysInfo.phi_type);
end


% myODE1     = @(t,x)predict_ode(x,learnInfo);
% myODE1     = @(t,x)predict_ode_interp(x,learnInfo,sysInfo);
myODE1     = @(t,x)predict_ode_mcuq(x,learnInfo,sysInfo);


train_time_vec=obsInfo.time_vec;
test_time_vec = [train_time_vec(1:end-1) linspace(train_time_vec(end),sysInfo.T_f,50)];% prediction on [0, T_f]                                                                     % final time the system will reach steady state
L = length(test_time_vec); %

traj_true            = zeros (dN,L,size(ICs, 2)); % true traj

traj_hat             = zeros (dN,L,size(ICs, 2)); % hat traj
%dtraj_hat             = zeros (dN,L,size(ICs, 2)); % hat traj
traj_hat_std             = zeros (dN,L,size(ICs, 2)); % hat traj


if learnInfo.dtraj 
    dtraj_true            = zeros (dN,L,size(ICs, 2)); % true traj
    dtraj_hat             = zeros (dN,L,size(ICs, 2));
end

traj_train_norm      = zeros(size(ICs, 2),1); % traj error over train time
traj_predict_norm    = zeros(size(ICs, 2),1); % traj error over prediction time
%C=zeros(d*N,length(test_time_vec),size(ICs,2));% convariance matrix for traj prediction

for m = 1 : size(ICs, 2)
  
 sol_true = ode45(myODE,solverInfo.time_span,ICs(:,m),solverInfo.option); % solu from adaptive solver
 traj_true(:,:,m)  =deval(sol_true,test_time_vec) ;% Nd x steps
 
 num_mc = 100;
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

 
 if learnInfo.dtraj
     dtraj_true(:,:,m) = Exact_derivative(traj_true(:,:,m),myODE);
     dtraj_hat(:,:,m)  = Exact_derivative(traj_hat(:,:,m),myODE1);
 end
 
 
 

%  for l=1:length(test_time_vec)
%         [~, C(:,l,m)] = predict_ode(traj_hat(:,l,m),learnInfo);
%  end
 
 %
 
  
 traj_train_norm(m)   = sqrt(max(ones(1,dN)*(traj_true(:,1:length(obsInfo.time_vec),m)-traj_hat(:,1:length(obsInfo.time_vec),m)).^2))./sqrt(max(ones(1,dN)*(traj_true(:,1:length(obsInfo.time_vec),m)).^2));
 traj_predict_norm(m) = sqrt(max(ones(1,dN)*(traj_true(:,length(obsInfo.time_vec)+1:end,m)-traj_hat(:,length(obsInfo.time_vec)+1:end,m)).^2))./sqrt(max(ones(1,dN)*(traj_true(:,length(obsInfo.time_vec)+1:end,m)).^2));

end




result.traj_true    = traj_true(1:d*N,:,:);
result.traj_hat     = traj_hat(1:d*N,:,:);

result.train_time_vec        = obsInfo.time_vec;
result.prediction_time_vec   = test_time_vec(length(obsInfo.time_vec)+1:end);

if learnInfo.dtraj
    result.dtraj_true    = dtraj_true(1:d*N,:,:);
    result.dtraj_hat    = dtraj_hat(1:d*N,:,:);
end

result.train_traj_error      = [mean(traj_train_norm) std(traj_train_norm)];
result.prediction_traj_error = [mean(traj_predict_norm) std(traj_predict_norm)];
%result.covmatrix = C;% Nd x L x M
end