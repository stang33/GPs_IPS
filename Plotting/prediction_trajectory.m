function [xpath_true,xpath_learn,xpath_true_new,xpath_learn_new,C,C_new] = prediction_trajectory(sysInfo,obsInfo,learnInfo,solverInfo)
% solve the ODE with a random inital condition 
% Input:
%   sysInfo  - parameters in the ODE and in its integrator
%            .N, d     : number of particles and dimension
%            .initdistr: initial distribution
%            .dt       : time step size
%            .t0,tEnd  : start and end time
%            .ODEoption: options for ODE solver
%   x0      - initial condition
% Output: 
%  xpath       - solution of the ODE    Nd x tN:   
%  dxpath   
% (c) XXXX
%
% ATTENTION:   
% % %  time   = dt:dt:tEnd;  ---- time instances of solution output NOT from t0 

%% basic setting of the system
N         = sysInfo.N;         % number of agents   
d         = sysInfo.d; 
myODE =  @(t,x) RHSfn(t,x,N,sysInfo.phi{1});
myODE1     = @(t,x)predict(x,learnInfo);

test_time_vec= [obsInfo.time_vec(1:end-1) linspace(obsInfo.time_vec(end),sysInfo.T_f,10)];                                                                     % final time the system will reach steady state



% trajectories with the same initial conditions as in training data
ICs = learnInfo.xpath_train(:,1,:);
xpath_true = zeros(d*N,length(test_time_vec),size(ICs,2));
xpath_learn =xpath_true;
for m=1:size(ICs,2)
sol = ode15s(myODE,solverInfo.time_span,ICs(:,m),solverInfo.option); % solu from adaptive solver
sol1 = ode15s(myODE1,solverInfo.time_span,ICs(:,m),solverInfo.option); % solu from adaptive solver
xpath_true(:,:,m) = deval(sol,test_time_vec);                   % interpolate solu for output
xpath_learn(:,:,m) =deval(sol1,test_time_vec);                   % interpolate solu for output
end

% # trajectories with random initial conditions for learning interaction kernel
xpath_true_new = xpath_true;
xpath_learn_new = xpath_learn;
for m=1:size(ICs,2)
x0 = sysInfo.mu0();
sol = ode15s(myODE,solverInfo.time_span,x0,solverInfo.option); % solu from adaptive solver
sol1= ode15s(myODE,solverInfo.time_span,x0,solverInfo.option);
xpath_true_new(:,:,m) = deval(sol,test_time_vec);                   % interpolate solu for output
xpath_learn_new(:,:,m) =  deval(sol1,test_time_vec);                    % interpolate solu for output
end

%%covariance estimation
C_new=zeros(d*N,length(test_time_vec),size(ICs,2));
C=zeros(d*N,length(test_time_vec),size(ICs,2));

for m=1:size(ICs,2)
    for l=1:length(test_time_vec)
        [~, covf1] = predict(xpath_learn_new(:,l,m),learnInfo);
        [~, covf2] = predict(xpath_learn_new(:,l,m),learnInfo);
        C(:,l,m) = diag(covf1);
        C_new(:,l,m) = diag(covf2);

    end
end



end    
    