function [dxpath_test,xpath_test, dxpath_train, xpath_train] = Generate_training_data(sysInfo,obsInfo,solverInfo)
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
%  xpath       - solution of the ODE    Nd x L:
%  dxpath      - the derivative function  Nd x L:
% (c) Sui Tang
%
% ATTENTION:
% % %  time   = dt:dt:tEnd;  ---- time instances of solution output NOT from t0

%% basic setting of the system
N         = sysInfo.N;         % number of agents
d         = sysInfo.d;         % dim of state vectors


switch sysInfo.ode_order
    case 1
        switch length(sysInfo.phi)
            case 1
               
                myODE     = @(t,x) RHSfn(t,x,N,sysInfo.phi{1});
                
            case 2
                
                myODE     = @(t,x) RHSfn_c(t,x,N,sysInfo.phi{1},sysInfo.phi{2});       
        
         end
        xpath_train = zeros(d*N,length(obsInfo.time_vec),obsInfo.M);
        dxpath_train = zeros(d*N,length(obsInfo.time_vec),obsInfo.M);

        test_time_vec= [obsInfo.time_vec(1:end-1) linspace(obsInfo.time_vec(end),sysInfo. T_f,200)];                                                                     % final time the system will reach steady state
        xpath_test = zeros(d*N,length(test_time_vec),obsInfo.M);
        dxpath_test = zeros(d*N,length(test_time_vec),obsInfo.M);
        
        if obsInfo.use_derivative
            
            for i = 1:obsInfo.M                                                                         % # trajectories with random initial conditions for learning interaction kernel
                x0 = sysInfo.mu0();
                sol = ode15s(myODE,solverInfo.time_span,x0,solverInfo.option); % solu from adaptive solver
                xpath_test(:,:,i) = deval(sol,test_time_vec);                   % interpolate solu for output
                dxpath_test(:,:,i) = Exact_derivative(xpath_test(:,:,i),myODE);
                dxpath_train(:,:,i) = dxpath_test(:,1:length(obsInfo.time_vec),i);
                xpath_train(:,:,i) =xpath_test(:,1:length(obsInfo.time_vec),i);                   % interpolate solu for output
            end
            
        else
            
            for i = 1:obsInfo.M                                                                         % # trajectories with random initial conditions for learning interaction kernel
                x0 = sysInfo.mu0();
                sol = ode15s(myODE,solverInfo.time_span,x0,solverInfo.option); % solu from adaptive solver
                xpath_test(:,:,i) = deval(sol,test_time_vec);                   % interpolate solu for output
                dxpath_test(:,:,i) = Exact_derivative(xpath_test,myODE);
                
                xpath_train(:,:,i) =xpath_test(:,1:length(obsInfo.time_vec),i);                   % interpolate solu for output
                dxpath_train(:,:,i) =  (deval(sol,obsInfo.time_vec+0.01)- deval(sol,obsInfo.time_vec))./0.01;%finite difference method
            end
            
        end
        
    case 2
        
        switch length(sysInfo.phi)
            case 1
               
                myODE     = @(t,y) RHSfn_2nd(t,y,N,sysInfo.phi{1});
                
            case 2
                
                myODE     = @(t,y) RHSfn_2nd_ncf(t,y,N,sysInfo.phi{1},sysInfo.phi{2},sysInfo.phi_type);
            
            case 3
                
                myODE     = @(t,y) RHSfn_2nd_tk(t,y,N,sysInfo.phi{1},sysInfo.phi{2},sysInfo.phi{3});

            case 4
                
                myODE     = @(t,y) RHSfn_2nd_ncf_xi(t,y,N,sysInfo.phi{1},sysInfo.phi{2},sysInfo.phi{3},sysInfo.phi{4});

        end
        xpath_train = zeros(d*2*N,length(obsInfo.time_vec),obsInfo.M);
        dxpath_train = zeros(d*2*N,length(obsInfo.time_vec),obsInfo.M);
        test_time_vec= [obsInfo.time_vec(1:end-1) linspace(obsInfo.time_vec(end),sysInfo. T_f,10)];                                                                     % final time the system will reach steady state
        
        xpath_test = zeros(d*2*N,length(test_time_vec),obsInfo.M);
        dxpath_test = zeros(d*2*N,length(test_time_vec),obsInfo.M);
        
        if obsInfo.use_derivative
            
            for i = 1:obsInfo.M                                                                         % # trajectories with random initial conditions for learning interaction kernel
                x0 = sysInfo.mu0();
                sol = ode15s(myODE,solverInfo.time_span,x0,solverInfo.option); % solu from adaptive solver
                xpath_test(:,:,i) = deval(sol,test_time_vec);                   % interpolate solu for output
                dxpath_test(:,:,i) = Exact_derivative(xpath_test(:,:,i),myODE);
                dxpath_train(:,:,i)= trajUnifNoiseAdditive(dxpath_test(:,1:length(obsInfo.time_vec),i), 0);
                xpath_train(:,:,i) =xpath_test(:,1:length(obsInfo.time_vec),i);                   % interpolate solu for output
            end
            
        else
            
            for i = 1:obsInfo.M                                                                         % # trajectories with random initial conditions for learning interaction kernel
                x0 = sysInfo.mu0();
                sol = ode15s(myODE,solverInfo.time_span,x0,solverInfo.option); % solu from adaptive solver
                xpath_test(:,:,i) = deval(sol,test_time_vec);                   % interpolate solu for output
                dxpath_test(:,:,i) = Exact_derivative(xpath_test,myODE);
                
                xpath_train(:,:,i) =xpath_test(:,1:length(obsInfo.time_vec),i);                   % interpolate solu for output
                dxpath_train(:,:,i) =  (deval(sol,obsInfo.time_vec+0.01)- deval(sol,obsInfo.time_vec))./0.01;%finite difference method
            end
            
        end
        
        
        
end


        
        
        
        
        
        
        
        
end













