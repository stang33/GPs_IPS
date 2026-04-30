function rhoLT = Generate_rhoT(sysInfo,obsInfo,solverInfo,saveON,plotON)
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
% Authored by Jinchao Feng and Sui Tang
%
% ATTENTION:
% % %  time   = dt:dt:tEnd;  ---- time instances of solution output NOT from t0

%% basic setting of the system


N         = sysInfo.N;         % number of agents
d         = sysInfo.d;         % dim of state vectors


switch sysInfo.ode_order
    case 1

        %myODE     = @(t,x) RHSfn_c(t,x,N,sysInfo.phi{1},sysInfo.phi{2});
       myODE     = @(t,x) RHSfn_NS(t,x,N,sysInfo.phi{1});
        %xpath_train = zeros(d*N,length(obsInfo.time_vec),obsInfo.MrhoT);
       xpath_train = zeros(d*N,length(obsInfo.time_vec),1);
        parfor i = 1:obsInfo.MrhoT                                                                         % # trajectories with random initial conditions for learning interaction kernel
            x0 = sysInfo.mu0();
            sol = ode15s(myODE,solverInfo.time_span,x0,solverInfo.option); % solu from adaptive solver
            xpath_train(:,:,i) = deval(sol,obsInfo.time_vec);                   % interpolate solu for output
            % interpolate solu for output
        end
        
        
        
    case 2
        
        switch length(sysInfo.phi)
            
            case 1
                if sysInfo.Noption == 1
                    myODE     = @(t,y) RHSfn_2nd_NS1(t,y,N,sysInfo.phi{1},sysInfo.phi_type);
                else
                    myODE     = @(t,y) RHSfn_2nd_NS(t,y,N,sysInfo.phi{1},sysInfo.phi_type);
                end
        end
        
        xpath_train = zeros(d*2*N,length(obsInfo.time_vec),obsInfo.M);

        parfor i = 1:obsInfo.MrhoT                                                                        % # trajectories with random initial conditions for learning interaction kernel
            x0 = sysInfo.mu0();
            sol = ode15s(myODE,solverInfo.time_span,x0,solverInfo.option); % solu from adaptive solver
            xpath_train(:,:,i) = deval(sol,obsInfo.time_vec);                   % interpolate solu for output
        end
        
end

rhoLT= rho_empirical(xpath_train,sysInfo,obsInfo,saveON,plotON);



end






















